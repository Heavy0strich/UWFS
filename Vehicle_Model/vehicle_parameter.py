
"""
Initialising the vehicle's parameters

Aman Tiwary
MSME - University of Washington
UW Formula Motorsports
Driverless
"""
import numpy as np
import pandas as pd

class Vehicle_parameter():
    def __init__(self):
        """
        initialising the vehicle's parameters
        """
        data = pd.read_excel('C:\\Users\\tiwar\\Desktop\\Project\\FSDS\\UWFS\\Vehicle_Model\\mass2.xlsx', engine='openpyxl')
        mass                    = data['Mass'][:-2].to_numpy()
        cgx                     = data['CGx'][:-2].to_numpy()
        cgy                     = data['CGy'][:-2].to_numpy()
        cgz                     = data['CGz'][:-2].to_numpy()

        self.total_mass         = mass.sum()                                                        # Total Mass of the vehicle(kg)
        Ix                      = np.multiply(mass, (np.square(cgy) + np.square(cgz)))              # Moment of Inertia about X-Axis(kg-m^2) for each component
        Iy                      = np.multiply(mass, (np.square(cgx) + np.square(cgz)))              # Moment of Inertia about Y-Axis(kg-m^2) for each component
        Iz                      = np.multiply(mass, (np.square(cgx) + np.square(cgy)))              # Moment of Inertia about Z-Axis(kg-m^2) for each component

        self.Izz                = Iz.sum()                                                          # Moment of Inertia about Z-Axis(kg-m^2) for the whole vehicle
        self.Ixx                = Ix.sum()                                                          # Moment of Inertia about X-Axis(kg-m^2) for the whole vehicle
        self.Iyy                = Iy.sum()                                                          # Moment of Inertia about Y-Axis(kg-m^2) for the whole vehicle
        # finding the (x, y, z) coordinates of the center of gravity
        A                       = np.array([[0, 1, 1],
                                            [1, 0, 1],
                                            [1, 1, 0]])

        square_XYZ              = np.linalg.inv(A) @ np.array([Ix.sum(), Iy.sum(), Iz.sum()])/self.total_mass

        cd_coordinates          = np.zeros(3)                                                        # Center of Gravity Coordinates

        for i in range(3):
            cd_coordinates[i]   = np.sqrt(square_XYZ[i])                                             # (x, y, z) coordinates of the center of gravity

        self.cg_height          = cd_coordinates[2]                                                  # Height of CG from ground(m)

        # Suspension Parameters
        self.wheelbase          = 1.562                                                              # Vehicle Wheelbase(m)
        self.Ftrack             = 1.2446                                                             # Front Track(m)
        self.Rtrack             = 1.1938                                                             # Rear Track(m)
        self.AVGtrack           = (self.Ftrack + self.Rtrack)/2                                      # Average Track(m)
        self.wheel_dia          = 0.457                                                              # Wheel Diameter(m)
        self.tire_pressure      = 34                                                                 # Tire Pressure(psi)                               

        # distance of COG to front and rear axle
        self.lf                 = cd_coordinates[0]
        self.lr                 = self.wheelbase - cd_coordinates[0]

        # distance between the wheels of front and rear axle
        self.bf                 = self.Ftrack
        self.br                 = self.Rtrack

        # Steering Parameters
        self.steering_ratio     = 1.5                                                               # Steering Ratio
        self.max_wheel_angle    = 35                                                                # Tire Angle(deg)

        # Aero Parameters
        self.COP_distance       = 0.25                                                              # Distance of COP from CG(m)
        self.COP_height         = 0.508                                                             # Height of COP from CG(m)
        self.Cl                 = -3.05                                                             # Lift Coefficient
        self.Cd                 = 1.44                                                              # Drag Coefficient
        self.frontArea          = 1                                                                 # Frontal Area(m^2)  % Aero frontal area for BrkG.drag and lift calculation. Usually normalized to 1 by wind tunnel(ft^2)

        # Tire Parameters
        self.static_camber_f    = -1.9                                                                 # Static Camber of Front Tires(deg)
        self.static_camber_r    = -1.5                                                                 # Static Camber of Rear Tires(deg)






        