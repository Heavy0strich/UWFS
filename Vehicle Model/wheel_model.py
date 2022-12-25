"""
This script is to be used to define the Wheel model to calculate the forces acting on the
contact patches

Aman Tiwary
MSME - University of Washington
UW Formula

"""
import numpy as np
import math
from vehicle_parameter import Vehicle_parameter

class Wheel_model:
    def __init__(self, del_w, alpha, v_cg, F_z, F_y):

        """
        Initialize the variables for calculation of wheel velocity at contact patch

        del_w   : tire turn angle is a vector of size (4 x 1)
        
        alpha   : Tire side slip angle
                : is a vector of size (4 x 1)

        v_cg    : Vehicle's Velocity at COG

        F_z     : Vertical Force acting at the wheel ground contact point
                : is a vector of size (4 x 1)

        F_y     : Lateral Force acting at the wheel ground contact point
                : is a vector of size (4 x 1)

        the above vectors store value for [Front left, Front Right, Rear Left, Rear Right]
        in the same order
        """
        
        self.del_w = del_w
        self.alpha = alpha
        self.v_cg = v_cg
        self.Fz = F_z
        self.Fy = F_y
        
        # Parameters for calculating caster
        self.l0 = -0.03                 # in m
        self.l1 = 0.12                  # in m
        self.Fz0 = 5000                 # in N
        self.c_press = 230000           # in N/m

        # Vehicle Parameters
        car = Vehicle_parameter()
        self.lf = car.lf
        self.lr = car.lr
        self.bf = car.bf
        self.br = car.br
        
        # distance between contact patch and COG in XY plane 
        self.cp_cg_distance = np.zeros((4, 1))

        # angle between contact patch and COG in XY plane
        self.cp_cg_angle = np.zeros((4, 1))

    def caster_distance(self):
        """
        Here we calculat the longitudinal (nl) and lateral(side)(ns) caster
        
        """

        # Dynamic Caster nl
        self.nl = 0.5 * (self.l0 * np.ones((4,1)) + ((self.l1/self.Fz0) * self.Fz))
        
        # ns takes into consideration lateral force influence on pressure distribution
        self.ns = 3 * np.multiply(self.nl, np.tan(self.alpha)) + self.Fy/self.c_press

    def contact_patch_cog_distance(self):
        """
        This function calculates the distance of the contact
        patches from the vehicle's COG
        """

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

        return self.cp_cg_distance
    
    def contact_patch_cog_angle(self):
        """
        This functions calculates the angle of the contact
        patches from the vehicle's COG in XY plane
        """

        self.cp_cg_angle[0] = np.arctan((self.bf/2 - self.ns[0] * np.cos(self.del_w[0])
         - self.nl[0] * np.sin(self.del_w[0]))/(self.lf - self.nl[0] * np.cos(self.del_w[0])
          + self.ns[0] ( np.sin(self.del_w[0]))))

        self.cp_cg_angle[1] = np.arctan((self.lf - self.nl[1] * np.cos(self.del_w[1]) 
        + self.ns[0] * np.sin(self.del_w[0]))/(self.bf/2 - self.ns[0] * np.cos(self.del_w[0])
         - self.nl[0] * np.sin(self.del_w[0])))

        self.cp_cg_angle[2] = np.arctan((self.lr + self.n[2])/(self.br/2- self.ns[2]))

        self.cp_cg_angle[3] = np.arctan((self.br/2 + self.ns[3])/(self.lr + self.nl[3]))

        return self.cp_cg_angle
    
    