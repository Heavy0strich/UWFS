"""
Bicycle Model

Aman Tiwary
MSME - University of Washington
UW Formula Motorsports
Driverless
"""
import os
import sys
sys.path.append('C:\\Users\\tiwar\\Desktop\\Project\\FSDS\\UWFS\\Vehicle_Model')

import numpy as np
from Vehicle_Model import Pacejka, wheel_model, utils
from Vehicle_Model.vehicle_parameter import Vehicle_parameter
import do_mpc

# Defining Vehicle parameters
vp = Vehicle_parameter()

mass                    = vp.total_mass
        
Izz                     = vp.Izz
lf                      = vp.lf
lr                      = vp.lr
# Aero Parameters
Cl                      = vp.Cl
Cd                      = vp.Cd
A                       = vp.frontArea
COP_distance            = vp.COP_distance
COP_height              = vp.COP_height
# Importing Tire Data
tire                    = utils.import_tire_data('R5_Hoosier 6.0-18.0-10 P12 LC0 Rim7.TIR')
static_camber_f         = tire['staticCamberF']
static_camber_r         = tire['staticCamberR']



model_type = 'continuous'
model = do_mpc.model.Model(model_type)

"""
Input:
        x[0]            : X         X coordinate of the vehicle in global frame
        x[1]            : Y         Y coordinate of the vehicle in global frame
        x[2]            : psi       heading of the chassis in global frame
        x[3]            : v_x       velocity of the center of gravity in body frame
        x[4]            : v_y       velocity of the center of gravity in body frame
        x[5]            : psi_dot   angular velocity of the chassis in body frame

        u[0]            : delta     steering angle of the front wheels
        u[1]            : D         acceleraion of the vehicle(+/- 1)
Output:
        dx_dt[0]        : X_dot         Velocity in global frame
        dx_dt[1]        : Y_dot
        dx_dt[2]        : psi_dot
        dx_dt[3]        : v_x_dot       acceleration(along longitudinal axis) in body frame
        dx_dt[4]        : v_y_dot       acceleration in body frame
        dx_dt[5]        : psi_ddot         angular acceleration in body frame


"""

x               = model.set_variable(var_type='_x', var_name='X', shape=(6,1))

u               = model.set_variable(var_type='_u', var_name='u', shape=(2,1))

# define time varying parameters
beta            = model.set_variable(var_type='_tvp', var_name='beta', shape=(1,1))

alpha           = model.set_variable(var_type='_tvp', var_name='alpha', shape=(2,1))

v_cg            = model.set_variable(var_type='_tvp', var_name='v_cg', shape=(1,1))

F_z             = model.set_variable(var_type='_tvp', var_name='F_z', shape=(2,1))

LSR_f           = model.set_variable(var_type='_tvp', var_name='LSR_f', shape=(1,1))

LSR_r           = model.set_variable(var_type='_tvp', var_name='LSR_r', shape=(1,1))

F_yf            = model.set_variable(var_type='_tvp', var_name='F_yf', shape=(1,1))

F_xf            = model.set_variable(var_type='_tvp', var_name='F_xf', shape=(1,1))

F_yr            = model.set_variable(var_type='_tvp', var_name='F_yr', shape=(1,1))

F_xr            = model.set_variable(var_type='_tvp', var_name='F_xr', shape=(1,1))

cp_cg_distance  = model.set_variable(var_type='_tvp', var_name='cp_cg_distance', shape=(2,1))

cp_cg_angle     = model.set_variable(var_type='_tvp', var_name='cp_cg_angle', shape=(2,1))

cp_velocity     = model.set_variable(var_type='_tvp', var_name='cp_velocity', shape=(2,1))

# TODO
# set_rhs for the non-linear model

# TODO
# Linearize the model wrt some given reference state and input

# TODO
# need to define get_tvp_template, that assigns value to the tvp and a set_tvp_fun which returns 
# appropiate tvp template






