"""
This script is to be used to define the Wheel model to calculate the Frictional forces acting on the
contact patches in x and y direction

Aman Tiwary
MSME - University of Washington
UW Formula Motorsports
Driverless
"""

import numpy as np
from vehicle_parameter import Vehicle_parameter
import utils

class Wheel_model:
    def __init__(self, del_w, beta, alpha, psi_dot, v_cg, F_z, F_y, model_type):

        """
        Initialize the variables for calculation of wheel velocity at contact patch

        method            : We have two methods for calculating contact patch velocities
                              - Method = 1 (Transformation of CoG Velocities)(Default)(use this method with 2 wheel model)
                              - Method = 2 (Calculating Individual curve radii)

        del_w             : tire turn angle is a vector of size (4 x 1)
                            - zero for rear tires if no rear wheel steering
                            - for bicycle mode, this is a scalar
        
        alpha             : Tire side slip angle
                          : is a vector of size (4 x 1)
                          : is a vector of size (2 x 1) for bicycle model

        beta              : Vehicle Body side slip angle
                          : is a single value

        psi_dot           : Yaw Rate
                          : is a single value

        v_cg              : Vehicle's Velocity at COG

        F_z               : Vertical Force acting at the wheel ground contact point
                          : is a vector of size (4 x 1)
                          : is a vector of size (2 x 1) for bicycle model

        F_y               : Lateral Force acting at the wheel ground contact point
                          : is a vector of size (4 x 1)
                          : is a vector of size (2 x 1) for bicycle model

        model_type       : 4 wheel model
                            - method_type = 4
                          : bicycle model
                            - method_type = 2 

        the above vectors store value for [Front left, Front Right, Rear Left, Rear Right]
        in the same order
        """
        
        #self.method             = method
        self.model_type         = model_type
        
        self.del_w              = del_w
        
        self.beta               = beta
        self.alpha              = alpha

        self.psi_dot            = psi_dot
        
        self.v_cg               = v_cg
        
        self.Fz                 = F_z
        self.Fy                 = F_y
        
        # Parameters for calculating caster
        self.l0                 = -0.03                 # in m
        self.l1                 = 0.12                  # in m
        self.Fz0                = 5000                  # in N
        self.c_press            = 230000                # in N/m

        # Vehicle Parameters
        car                     = Vehicle_parameter()
        self.lf                 = car.lf
        self.lr                 = car.lr
        self.bf                 = car.bf
        self.br                 = car.br


    def caster_distance(self):
        """
        Here we calculat the longitudinal (nl) and lateral(side)(ns) caster
        
        """
        #self.alpha              = utils.Slip_Angle(self.model_type, self.v_cg, self.psi_dot, self.lf, self.lr, self.beta, self.del_w)

        if self.model_type == 4:
            # Dynamic Caster nl
            self.nl = 0.5 * (self.l0 * np.ones((4,1)) + ((self.l1/self.Fz0) * self.Fz))
        
            # ns takes into consideration lateral force influence on pressure distribution
            self.ns = 3 * np.multiply(self.nl, np.tan(self.alpha)) + self.Fy/self.c_press
        elif self.model_type == 2:
            # Dynamic Caster nl
            self.nl = 0.5 * (self.l0 * np.ones((2,1)) + ((self.l1/self.Fz0) * self.Fz))
        
            # ns takes into consideration lateral force influence on pressure distribution
            self.ns = 3 * np.multiply(self.nl, np.tan(self.alpha)) + self.Fy/self.c_press

    def contact_patch_cog_distance(self):
        """
        This function calculates the distance of the contact
        patches from the vehicle's COG
        cp_cg_distance    : distance between contact patch and COG in XY plane
                          : if model_type = 4
                            - is a vector of size (4 x 1) where [FL, FR, RL, RR](tire order)
                          : if model_type = 2
                            - is a vector of size (2 x 1) where [F, R](tire order)
        """
        self.caster_distance()
        if self.model_type == 4:
            # distance between contact patch and COG in XY plane 
            self.cp_cg_distance     = np.zeros((4, 1))


            self.cp_cg_distance[0] = ((self.lf - self.nl[0] * np.cos(self.del_w[0]) + 
                        self.ns[0] * np.sin(self.del_w[0]))**2 + 
                        (self.bf/2 - self.ns[0] * np.cos(self.del_w[0]) -
                        self.nl[0] * np.sin(self.del_w[0]))**2)**(1/2)
            
            self.cp_cg_distance[1] = ((self.lf - self.nl[1] * np.cos(self.del_w[1]) + 
                        self.ns[1] * np.sin(self.del_w[1]))**2 + 
                        (self.bf/2 + self.ns[1] * np.cos(self.del_w[1]) + 
                        self.nl[1] * np.sin(self.del_w[1]))**2)**(1/2)

            self.cp_cg_distance[2] = ((self.lr - self.nl[2] * np.cos(self.del_w[2]) + 
                        self.ns[2] * np.sin(self.del_w[2]))**2 + 
                        (self.br/2 - self.ns[2] * np.cos(self.del_w[2]) -
                        self.nl[2] * np.sin(self.del_w[2]))**2)**(1/2)
            
            self.cp_cg_distance[3] = ((self.lr - self.nl[3] * np.cos(self.del_w[3]) + 
                        self.ns[3] * np.sin(self.del_w[3]))**2 + 
                        (self.br/2 + self.ns[3] * np.cos(self.del_w[3]) + 
                        self.nl[3] * np.sin(self.del_w[3]))**2)**(1/2)
        elif self.model_type == 2:
            # distance between contact patch and COG in XY plane 
            self.cp_cg_distance     = np.zeros((2, 1))

            self.cp_cg_distance[0] = ((self.lf - self.nl[0] * np.cos(self.del_w[0]) + 
                        self.ns[0] * np.sin(self.del_w[0]))**2 + 
                        (self.ns[0] * np.cos(self.del_w[0]) +
                        self.nl[0] * np.sin(self.del_w[0]))**2)**(1/2)
            
            self.cp_cg_distance[1] = ((self.lr + self.nl[1])**2 + 
                        (self.ns[1])**2)**(1/2)

        return self.cp_cg_distance
    
    def contact_patch_cog_angle(self):
        """
        This functions calculates the angle of the contact
        patches from the vehicle's COG in XY plane
        """

        self.caster_distance()
        if self.model_type == 4:
            
            # angle between contact patch and COG in XY plane
            self.cp_cg_angle        = np.zeros((4, 1))        
            
            self.cp_cg_angle[0] = np.arctan((self.bf/2 - self.ns[0] * np.cos(self.del_w[0])
                        - self.nl[0] * np.sin(self.del_w[0]))/(self.lf - self.nl[0] * np.cos(self.del_w[0])
                        + self.ns[0] ( np.sin(self.del_w[0]))))

            self.cp_cg_angle[1] = np.arctan((self.lf - self.nl[1] * np.cos(self.del_w[1]) 
                        + self.ns[0] * np.sin(self.del_w[0]))/(self.bf/2 - self.ns[0] * np.cos(self.del_w[0])
                        - self.nl[0] * np.sin(self.del_w[0])))

            self.cp_cg_angle[2] = np.arctan((self.lr + self.n[2])/(self.br/2- self.ns[2]))

            self.cp_cg_angle[3] = np.arctan((self.br/2 + self.ns[3])/(self.lr + self.nl[3]))
        
        elif self.model_type == 2:
            
            # angle between contact patch and COG in XY plane
            self.cp_cg_angle        = np.zeros((2, 1))        
            
            self.cp_cg_angle[0] = np.arctan((self.ns[0] * np.cos(self.del_w[0])
                        + self.nl[0] * np.sin(self.del_w[0]))/(self.lf - self.nl[0] * np.cos(self.del_w[0])
                        + self.ns[0] * np.sin(self.del_w[0])))

            self.cp_cg_angle[1] = np.arctan((self.ns[1])/(self.lr + self.nl[1]))

        return self.cp_cg_angle


    """Not Required
    def contact_patch_velocities(self):
        
        if self.method == 1:
            return utils.transform_cog_vel(self.v_cg, self.psi_dot, self.cp_cg_distance, self.cp_cg_angle, self.model_type)
        else:
            return utils.ICM(self.v_cg, self.psi_dot, self.bf, self.lf, self.beta, self.model_type)    
        
    """

    def validity_check(self):
        
        try:
            assert self.method == 1 or self.method == 0
        except:
            print(f"the method for contact patch velocities can either be 1(transform_cog_ve) or 0(ICM) The value entered was:{n}")
            return True
        
        try:
            assert self.model_type == 2 or self.model_type == 4
        except:
            print(f"You can use only 4 wheel model or 2 wheel(Bicycle) mode. The value entered was:{self.model_type}")
            return True
    