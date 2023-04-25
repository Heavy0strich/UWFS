"""
Bicycle Model

Aman Tiwary
MSME - University of Washington
UW Formula Motorsports
Driverless
"""

import numpy as np
import utils
from vehicle_parameter import Vehicle_parameter
from Pacejka import Pacejka
from wheel_model import Wheel_model

def numerical_simulation(f,t,x,t0=0.,dt=1e-4,ut=None,ux=None,utx=None,return_u=False):
    """
    simulate x' = f(t,x,u) 

    input:
    f : R x X x U --> X - vector field
        X - state space (must be vector space)
        U - control input set
    t - scalar - final simulation time
    x - initial condition; element of X

    (optional:)
    t0 - scalar - initial simulation time
    dt - scalar - stepsize parameter
    return_u - bool - whether to return u_

    (only one of:)
    ut : R --> U
    ux : X --> U
    utx : R x X --> U

    output:
    t_ - N array - time trajectory
    x_ - N x X array - state trajectory
    (if return_u:)
    u_ - N x U array - state trajectory
    """
    t_,x_,u_ = [t0],[x],[]

    inputs = sum([1 if u is not None else 0 for u in [ut,ux,utx]])
    assert inputs <= 1, "more than one of ut,ux,utx defined"

    if inputs == 0:
        assert not return_u, "no input supplied"
    else:
        if ut is not None:
            u = lambda t,x : ut(t)
        elif ux is not None:
            u = lambda t,x : ux(x)
        elif utx is not None:
            u = lambda t,x : utx(t,x)

    while t_[-1]+dt < t:
        if inputs == 0:
            _t,_x = t_[-1],x_[-1]
            dx = f(t_[-1],x_[-1]) * dt
        else:
            _t,_x,_u = t_[-1],x_[-1],u(t_[-1],x_[-1])
            dx = f(_t,_x,_u) * dt
            u_.append( _u )

        x_.append( _x + dx )
        t_.append( _t + dt )

    if return_u:
        return np.asarray(t_),np.asarray(x_),np.asarray(u_)
    else:
        return np.asarray(t_),np.asarray(x_)
    

class Dynamic_Bicycle_Model():

    def __init__(self):
        """
        
        Initialising the parameters of the dynamic bicycle model at each time step
        
        self.mass = mass
        self.I_zz = I_zz
        self.F_Lf = F_Lf
        self.F_Sf = F_Sf
        self.F_Lr = F_Lr
        self.F_Sr = F_Sr
        self.cp_cg_distance = cp_cg_distance
        
        """
        vp = Vehicle_parameter()

        self.mass         = vp.total_mass
        
        self.Izz          = vp.Izz

        self.lf           = vp.lf

        self.lr           = vp.lr

        # Aero Parameters
        self.Cl           = vp.Cl
        self.Cd           = vp.Cd
        self.A            = vp.frontArea
        self.COP_distance = vp.COP_distance
        self.COP_height   = vp.COP_height

        # Importing Tire Data
        self.tire       = utils.import_tire_data('R5_Hoosier 6.0-18.0-10 P12 LC0 Rim7.TIR')
        self.static_camber_f = self.tire['staticCamberF']
        self.static_camber_r = self.tire['staticCamberR']




    def nolinear_function(self, t, x, u):
        # TODO: define the states and dynamics!!
        """
        Input:
            x[0]        : X         X coordinate of the vehicle in global frame
            x[1]        : Y         Y coordinate of the vehicle in global frame
            x[2]        : psi       heading of the chassis in global frame
            x[3]        : v_x       velocity of the center of gravity in body frame
            x[4]        : v_y       velocity of the center of gravity in body frame
            x[5]        : psi_dot   angular velocity of the chassis in body frame

            u[0]        : delta     steering angle of the front wheels
            u[1]        : acceleraion of the vehicle
        Output:
            dx_dt[0]        : X_dot         Velocity in global frame
            dx_dt[1]        : Y_dot
            dx_dt[2]        : psi_dot
            dx_dt[3]        : v_x_dot       acceleration(along longitudinal axis) in body frame
            dx_dt[4]        : v_y_dot       acceleration in body frame
            dx_dt[5]        : psi_ddot         angular acceleration in body frame


        """

        
        X               = x[0]                                                  # X coordinate of the vehicle in global frame   

        Y               = x[1]                                                  # Y coordinate of the vehicle in global frame

        psi             = x[2]                                                  # angular velocity Chassis

        v_cg_x          = x[3]                                                  # velocity of the center of gravity in body frame

        v_cg_y          = x[4]                                                  # velocity of the center of gravity in body frame

        psi_dot         = x[5]                                                  # angular acceleration Chassis 

        beta            = np.arctan2(v_cg_y, v_cg_x)                            # Chassis Side Slip Angle

        v_cg            = np.sqrt(v_cg_x**2 + v_cg_y**2)                        # velocity of the center of gravity in body frame

        F_z             = utils.Fz(self.mass, v_cg, self.lf, self.lr, self.Cl, self.Cd, self.A, self.COP_distance, self.COP_height, u[1])

        beta            = utils.Vehicle_Body_Side_Slip_Angle(v_cg_x, v_cg_y, psi)

        alpha           = utils.Slip_Angle(2, v_cg, psi_dot, self.lf, self.lr, beta, u[0])

        wm              = Wheel_model(u[0], beta, alpha, psi_dot, v_cg, F_z, 2)
        cp_cg_distance  = wm.contact_patch_cog_distance()
        cp_cg_angle     = wm.contact_patch_cog_angle()

        cp_velocity     = utils.transform_cog_vel(v_cg, psi_dot, cp_cg_distance, cp_cg_angle, beta, 2)

        if v_cg == 0:
            #TODO Better way to update LSR and adding the affect of Torque Request on LSR
            self.LSR_f      = utils.Longitudinal_Slip_ratio(2, 0, 0, 0)
            self.LSR_r      = utils.Longitudinal_Slip_ratio(2, 0, 0, 0)
        else:
            self.LSR_f      = utils.Longitudinal_Slip_ratio(2, v_cg_x, cp_velocity[0], alpha[0])
            self.LSR_r      = utils.Longitudinal_Slip_ratio(2, v_cg_x, cp_velocity[1], alpha[1])

        pj_f            = Pacejka(self.tire, F_z[0], alpha[0], self.LSR_f, self.static_camber_f, 0, v_cg_x, 0)
        F_yf            = pj_f.request('Fy')   
        F_xf            = pj_f.request('Fx')

        pf_r            = Pacejka(self.tire, F_z[1], alpha[1], self.LSR_r, self.static_camber_r, 0, v_cg_x, 0)
        F_yr            = pf_r.request('Fy')
        F_xr            = pf_r.request('Fx')

        Fy              = np.array([F_yf, F_yr])
        
        
        ##TODO: Add Longitudinal Force generated from drivetrain

        dx_dt                     = [0, 0, 0, 0, 0, 0]

        dx_dt[0], dx_dt[1]        = utils.coordinate_transform_Interial_to_Vehicle(v_cg_x, v_cg_y, psi)

        dx_dt[2]                  = psi_dot
    
        dx_dt[3]                  = (F_xr + F_xf * np.cos(u[0]) - F_yf * np.sin(u[0]) + self.mass * v_cg_y * psi_dot)/self.mass
    
        dx_dt[4]                  = (F_yr + F_xf * np.sin(u[0]) + F_yf * np.cos(u[0]) - self.mass * v_cg_x * psi_dot)/self.mass
    
        dx_dt[5]                  = ((F_xf * np.sin(u[0]) + F_yf * np.cos(u[0])) * cp_cg_distance[0]  - F_yr * cp_cg_distance[1])/self.Izz

        return np.asarray(dx_dt)

