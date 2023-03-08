"""
Aman Tiwary
MSME - University of Washington
UW Formula Motorsports
Driverless
"""
import numpy as np

def coordinate_transform_Vehicle_to_Inertial(x, y, psi):
    """
    This function transforms the coordinates from vehicle frame to inertial frame
    """
    # Inertial frame coordinates
    X_inertial = np.array([[np.cos(psi), -np.sin(psi)], 
                           [np.sin(psi), np.cos(psi)]]) @ np.array([[x], [y]])

    return X_inertial

def coordinate_transform_Interial_to_Vehicle(x, y, psi):
    """
    This function transforms the coordinates from inertial frame to vehicle frame
    """
    # Inertial frame coordinates
    X_vehicle = np.array([[np.cos(psi), np.sin(psi)], 
                          [-np.sin(psi), np.cos(psi)]]) @ np.array([[x], [y]])

    return X_vehicle

def transform_cog_vel(v_cg, psi_dot, cp_cg_distance, cp_cg_angle, beta, model_type):
        """
        The Wheel ground contact patch velocities consists of 2 components:
            - The component due to CoG velocity
            - The component due to yaw(motion about vertical axis)
        
        The contact patch velocities are calculated in the vehicle frame
        """
        
        if model_type == 4:
            # contact patch velocities 
            #   cp_velocities[:, 0] : Longitudinal CoG direction
            #   cp_velocities[:, 1] : Lateral CoG direction
            cp_velocities       = np.zeros((4, 2))

            cp_velocities[0, 0] = v_cg * np.cos(beta) - psi_dot * cp_cg_distance[0] * np.sin(cp_cg_angle[0])
            cp_velocities[0, 1] = v_cg * np.sin(beta) + psi_dot * cp_cg_distance[0] * np.cos(cp_cg_angle[0])

            cp_velocities[1, 0] = v_cg * np.cos(beta) + psi_dot * cp_cg_distance[1] * np.cos(cp_cg_angle[1])
            cp_velocities[1, 1] = v_cg * np.sin(beta) + psi_dot * cp_cg_distance[1] * np.cos(cp_cg_angle[1])

            cp_velocities[2, 0] = v_cg* np.cos(beta) - psi_dot * cp_cg_distance[2] * np.cos(cp_cg_angle[2])
            cp_velocities[2, 1] = v_cg* np.sin(beta) - psi_dot * cp_cg_distance[2] * np.sin(cp_cg_angle[2])

            cp_velocities[3, 0] = v_cg * np.cos(beta) + psi_dot * cp_cg_distance[3] * np.cos(cp_cg_angle[3])
            cp_velocities[3, 1] = v_cg * np.sin(beta) - psi_dot * cp_cg_distance[3] * np.cos(cp_cg_angle[3])

        elif model_type == 2:
            cp_velocities = np.zeros((2, 2))

            cp_velocities[0, 0] = v_cg * np.cos(beta) - psi_dot * cp_cg_distance[0] * np.sin(cp_cg_angle[0])
            cp_velocities[0, 1] = v_cg * np.sin(beta) + psi_dot * cp_cg_distance[0] * np.cos(cp_cg_angle[0])

            cp_velocities[1, 0] = v_cg * np.cos(beta) + psi_dot * cp_cg_distance[1] * np.cos(cp_cg_angle[1])
            cp_velocities[1, 1] = v_cg * np.sin(beta) - psi_dot * cp_cg_distance[1] * np.cos(cp_cg_angle[1])

        return cp_velocities

def ICM(v_cg, psi_dot, bf, br, lf, lr, beta):
    """
    When a vehicle is seen from a birds eye view during a turn, each wheel
    follows an individual curve. The velocities of the CoG and the 
    instantaneous center of Motion(ICM) are perpendicular

    Assumption: 
        - R   : Distance from the CoG to ICM is much larger than
        individual curve radii R_ij
        - Disregard the caster effect

    We Calculate the four differential radii
    #TODO: Derive for bicycle model

    """

    cp_velocities_ICM = np.zeros((4, 1))

    cp_velocities_ICM[0] = v_cg - psi_dot * (bf/2 - lf * beta)
    
    cp_velocities_ICM[1] = v_cg + psi_dot * (bf/2 + lf * beta)
    
    cp_velocities_ICM[0] = v_cg - psi_dot * (br/2 + lr * beta)
    
    cp_velocities_ICM[0] = v_cg + psi_dot * (br/2 - lr * beta)

    return cp_velocities_ICM


def import_tire_data(tire_file):
    """
    Input:
        tire_file = 'R5_Hoosier 6.0-18.0-10 P12 LC0 Rim7.TIR'   <Tire data File>
    Output:
        tire dictionary with variable name as dictionary keys and values as 
        dctionary value
    """
    
    ident = open(tire_file)
    tire = {}
    i = 0
    while True:
        line = ident.readline()

        if not line:
            break

        if line[0] != '[' and line[0] != '$' and line[0] != '{':
            words = line.split('=', 2)
            if len(words) == 2:
                pos = words[1].find('$')
                try:
                    tire[words[0].rstrip()] = float(words[1][0:pos].rstrip())
                except:
                    tire[words[0].rstrip()] = words[1][0:pos].rstrip()

    return tire 

def Vehicle_Body_Side_Slip_Angle(v_cg_x, v_cg_y, psi):
    """
    This function calculates the body side slip angle beta
    v_cg : CoG velocity in the inertial frame([X, Y]) Directions
    """
    beta = 0

    try:
        beta = np.arctan(v_cg_y/v_cg_x) - psi
    except:
        print("The vehicle is not moving. The body side slip angle is undefined")

    return beta

def Slip_Angle(model_type, v_cg, psi_dot, lf, lr, beta, del_w):
    if model_type == 2:
        """
        This function calculates the slip angle of the front and rear wheels
        it is a vector of size 2. The first element is the slip angle of the front wheel and
        the second element is the slip angle of the rear wheel.
        """
        alpha = np.zeros((2, 1))

        try:
            alpha[0] = - np.arctan((v_cg * np.sin(beta) + lf * psi_dot)/v_cg * np.cos(beta)) + del_w
            alpha[1] = np.arctan((- v_cg * np.sin(beta) + lr * psi_dot)/v_cg * np.cos(beta))
        except:
            print("The vehicle is not moving. The slip angle is undefined")

    return alpha

def Longitudinal_Slip_ratio(model_type, v_r, v_cp, alpha):
    if model_type == 2:
        """
        This function calculates the Longitudinal slip ratio of the Wheels
        Calculates one wheel at a time
        """
        k = 0

        try:
            if v_r * np.cos(alpha) - v_cp <= 0:
                k = (v_r * np.cos(alpha) - v_cp) / v_cp
            elif v_r * np.cos(alpha) - v_cp > 0:
                k = (v_r * np.cos(alpha) - v_cp)/v_r * np.cos(alpha)
        except:
            print("The vehicle is not moving. The slip ratio is undefined")

    return k

def Aerodynamic_Force(v_cg, A, Cd):
    """
    This function calculates the aerodynamic Lift/Drag force on the vehicle depending on the 
    coefficient used
    """
    TempC           = 25                # Temperature in Celsius
    Baro            = 101400            # Barometric pressure in Pa
    Rair            = 287.04            # Gas constant for air in J/(kg*K) 
    TempK           = TempC + 273.15    # Temperature in Kelvin
    
    
    rho = Baro / (Rair * TempK) # Air density in Kg/m^3
    Fd = 0.5 * rho * A * Cd * v_cg**2

    return Fd


def Fz(total_mass, v_cg, lf, lr, Cl, Cd, A, COP_distance, COP_height, CG_height, v_cg_dot):
    """
    This function calculates the vertical load on each wheel
    """
    # Aerodynamic Forces
    F_aero_lift     = Aerodynamic_Force(v_cg, A, Cl)

    F_aero_drag     = Aerodynamic_Force(v_cg, A, Cd)

    Fz_aero_front   = (-(F_aero_lift * COP_distance) - F_aero_drag * COP_height/(lf + lr))

    Fz_aero_rear    = (-(F_aero_lift * (1 - COP_distance)) + F_aero_drag * COP_height/(lf + lr))

    # Load Transfer
    Fz_load_transfer_front = lr * total_mass * 9.81 / (lf + lr) - CG_height / (lf + lr) (total_mass/9.81 * v_cg_dot * CG_height)

    Fz_load_transfer_rear = lf * total_mass * 9.81 / (lf + lr) + CG_height / (lf + lr) (total_mass/9.81 * v_cg_dot * CG_height)

    # Total Veritical Load
    Fz_front = Fz_aero_front + Fz_load_transfer_front

    Fz_rear = Fz_aero_rear + Fz_load_transfer_rear

    return np.array([Fz_front, Fz_rear])

def Rotational_equivalent_Wheel_Velocity(v_cg, beta):
    """
    This function calculates the rotational equivalent of Wheel velocity. Here we assume the chassis and the wheel mount as rigid bodies.
    Thus, we have v_cog_x at the wheel center.
    ## TODO: Insted, use the torque request and calculate the wheel's angular velocity
    """

    v_r = v_cg * np.cos(beta)

    return v_r

    
