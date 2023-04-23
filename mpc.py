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
from Vehicle_Model import utils
from Vehicle_Model.vehicle_parameter import Vehicle_parameter
from Vehicle_Model.wheel_model import Wheel_model
from Vehicle_Model.Pacejka import Pacejka
import do_mpc
from casadi import vertcat

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
radius                  = vp.wheel_dia
max_wheel_angle         = vp.max_wheel_angle
steering_ratio          = vp.steering_ratio

C_m                     =  mass * 9.81
# Importing Tire Data
tire                    = utils.import_tire_data('R5_Hoosier 6.0-18.0-10 P12 LC0 Rim7.TIR')
static_camber_f         = tire['staticCamberF']
static_camber_r         = tire['staticCamberR']
tire_pressure           = vp.tire_pressure * 0.069                                                      # Converting from psi to bar    



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

#x              = model.set_variable(var_type='_x', var_name='x', shape=(6,1))
X               = model.set_variable(var_type='_x', var_name='X', shape=(1,1))
Y               = model.set_variable(var_type='_x', var_name='Y', shape=(1,1))
psi             = model.set_variable(var_type='_x', var_name='psi', shape=(1,1))
v_x             = model.set_variable(var_type='_x', var_name='v_x', shape=(1,1))
v_y             = model.set_variable(var_type='_x', var_name='v_y', shape=(1,1))
psi_dot         = model.set_variable(var_type='_x', var_name='psi_dot', shape=(1,1))


#u               = model.set_variable(var_type='_u', var_name='u', shape=(2,1))
delta           = model.set_variable(var_type='_u', var_name='delta', shape=(1,1))
D               = model.set_variable(var_type='_u', var_name='D', shape=(1,1))


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

tau             = model.set_variable(var_type='_tvp', var_name='tau', shape=(1,1))

cp_cg_distance  = model.set_variable(var_type='_tvp', var_name='cp_cg_distance', shape=(2,1))

cp_cg_angle     = model.set_variable(var_type='_tvp', var_name='cp_cg_angle', shape=(2,1))

cp_velocity     = model.set_variable(var_type='_tvp', var_name='cp_velocity', shape=(2,1))

del_e           = model.set_variable(var_type='_z', var_name='del_e', shape=(1,1))



# Defining variables for reference trajectory
_X               = model.set_variable(var_type='_tvp', var_name='_X', shape=(1,1))
_Y               = model.set_variable(var_type='_tvp', var_name='_Y', shape=(1,1))
_psi             = model.set_variable(var_type='_tvp', var_name='_psi', shape=(1,1))
_v_x             = model.set_variable(var_type='_tvp', var_name='_v_x', shape=(1,1))
_v_y             = model.set_variable(var_type='_tvp', var_name='_v_y', shape=(1,1))
_psi_dot         = model.set_variable(var_type='_tvp', var_name='_psi_dot', shape=(1,1))





# set_rhs for the non-linear model

model.set_rhs('X', utils.coordinate_transform_Interial_to_Vehicle(model.x['v_x'], model.x['v_y'], model.x['psi'])[0])
model.set_rhs('Y', utils.coordinate_transform_Interial_to_Vehicle(model.x['v_x'], model.x['v_y'], model.x['psi'])[1])
model.set_rhs('psi', psi_dot)
model.set_rhs('v_x', (1/mass) * (F_xr + F_xf * np.cos(model.u['delta']) - F_yf * np.sin(model.u['delta']) + mass * model.x['v_y'] * model.x['psi_dot'] + tau / radius ))
model.set_rhs('v_y', (1/mass) * (F_yr + F_xf * np.sin(model.u['delta']) + F_yf * np.cos(model.u['delta']) - mass * model.x['v_y'] * model.x['psi_dot']))
model.set_rhs('psi_dot', (1/Izz) * ((F_xf * np.sin(delta) + F_yf * np.cos(delta)) * cp_cg_distance[0] - F_yr * cp_cg_distance[1]))

model.set_rhs('del_e', np.sqrt((X - _X)**2 + (Y - _Y)**2))

x               = vertcat(X, Y, psi, v_x, v_y, psi_dot)
x_              = vertcat(_X, _Y, _psi, _v_x, _v_y, _psi_dot)

model.set_expression(expr_name='del_x', expr = (x - x_))

model.setup()

# Configuring The MPC Controller
mpc = do_mpc.controller.MPC(model)

# optimizer parameters

setup_mpc = {
        'n_horizon': 20,
        't_step': 0.1,
        'n_robust': 1,
        'store_full_solution': True
}
mpc.set_param(**setup_mpc)

# Objective Function
Q = np.diag([5, 5, 5, 1, 1, 1])
m_term = model.aux['del_x'].T @ Q @ model.aux['del_x']

l_term = model.aux['del_x'].T @ Q @ model.aux['del_x']

mpc.set_objective(m_term=m_term, l_term=l_term)

mpc.set_rterm(u = np.array([1, 1]).reshape(-1, 1))

# Define Constraints
# Lower and upper bounds on states 

mpc.bounds['lower', '_u', 'D'] = -1

mpc.bounds['upper', '_u', 'D'] = 1

mpc.bounds['lower', '_u', 'delta'] = -max_wheel_angle * 0.0174533
mpc.bounds['upper', '_u', 'delta'] = max_wheel_angle * 0.0174533

mpc.bounds['lower', '_u', 'D'] = -1
mpc.bounds['upper', '_u', 'D'] = 1

mpc.bounds['upper', '_z', 'del_e'] =  2.8  # Maximuum distance from the reference path

mpc.bounds['lower', '_tvp', 'tau'] = - mass * 9.81 * radius                                     # Maximum braking torque
mpc.bounds['upper', '_tvp', 'tau'] = mass * 9.81 * radius                                       # Maximum acceleration torque
mpc.setup()

# Configuring the simulator
simulator = do_mpc.simulator.Simulator(model)
simulator.set_param(t_step = 0.1)

# TODO
# need to define get_tvp_template, that assigns value to the tvp and a set_tvp_fun which returns 
# appropiate tvp template
tvp_template = simulator.get_tvp_template()


def tv_fun(t_now):
    tvp_template['beta'] = np.arctan2(model.x['v_y'], model.x['v_x'])
    tvp_template['v_cg'] = np.sqrt(model.x['v_x']**2 + model.x['v_y']**2)
    tvp_template['alpha'] = utils.Slip_Angle(2, model.tvp['v_cg'], model.x['psi_dot'], lf, lr, model.tvp['beta'], model.u['delta'])
    tvp_template['F_z']   = utils.Fz(mass, model.tvp['v_cg'], lf, lr, Cl, Cd, A, COP_distance, COP_height, model.u['D'])
    wm = Wheel_model(model.u['delta'], model.tvp['beta'], model.tvp['alpha'], model.x['psi_dot'], model.tvp['v_cg'], model.tvp['F_z'], 2)
    tvp_template['cp_cg_distance'] = wm.contact_point_cog_distance()
    tvp_template['cp_cg_angle'] = wm.contact_point_cog_angle()
    tvp_template['cp_velocity'] = utils.transform_cog_vel(model.tvp['v_cg'], model.x['psi_dot'], model.tvp['cp_cg_distance'], model.tvp['cp_cg_angle'], model.tvp['beta'], 2)
    tvp_template['LSR_f']  = utils.Longitudinal_Slip_ratio(2, model.x['v_x'], model.tvp['cp_velocity', 0], model.tvp['alpha', 0])
    tvp_template['LSR_r']  = utils.Longitudinal_Slip_ratio(2, model.x['v_x'], model.tvp['cp_velocity', 1], model.tvp['alpha', 1])
    tvp_template['tau'] = (C_m * model.u['D'] - utils.Rolling_resistance(model.tvp['v_cg'], tire_pressure) - utils.Aerodynamic_Force(model.tvp['v_cg'], A, Cd)) * radius
    pj_f  = Pacejka(tire, model.tvp['F_z', 0], model.tvp['alpha', 0], model.tvp['LSR_f'], static_camber_f, 0, model.x['v_x'], 0)
    tvp_template['F_yf'] = pj_f.request('Fy')
    tvp_template['F_xf'] = pj_f.request('Fx')

    pj_r = Pacejka(tire, model.tvp['F_z', 1], model.tvp['alpha', 1], model.tvp['LSR_r'], static_camber_r, 0, model.x['v_x'], 0)
    tvp_template['F_yr'] = pj_r.request('Fy')
    tvp_template['F_xr'] = pj_r.request('Fx')

    tvp_template['del_e']

    return tvp_template

simulator.set_tvp_fun(tv_fun)
simulator.setup()

# Assume direct state feedback, select an initial state
x0 = np.array([0, 0, 0, 1, 1, 0.1]).reshape(-1, 1)

simulator.x0 = x0
mpc.x0 = x0

mpc.set_initial_guess()


for i in range(200):
    u0 = mpc.make_step(x0)
    x0 = simulator.make_step(u0)

# TODO
# Linearize the model wrt some given reference state and input






