"""
Calculates resultant tire forces using Pacejka coefficients imp

Aman Tiwary
MSME - University of Washington
UW Formula Motorsports
Driverless

"""

import numpy as np

class Pacejka:
    def __init__(self, tire, Fz, SA, SR, IA, mu, V, tF):
        """
            tire:    dictionary with variable name as dictionary keys and values as 
                     dctionary value
            Fz  :   normal tire force (- is downward)
            SA  :   Slip angle in degrees
            SR  :   Slip Ratio in %
            IA  :   Inclination Angle
            mu  :   Friction Coefficient 
            V   :   Rotational Equivalent of wheel velocity ()
        """

        # Print Range warnings to screen
        self.warn_on         =   1

        # Check for range of requested inputs vs test range, warn if user exceeds the range
        self.check_range     =   1

        # Set Limit on peak finding and reverse lookup iteration counts
        self.max_itr         =   50

        # For safety and sanity, limit the output SA and SR peak finders to some
        # multiple of the max tested value
        self.peak_limit      =   2
        
        # Initialisation
        self.tire            =   tire
        #self.type            =   type
        self.Fz              =   Fz
        self.SA              =   SA
        self.SR              =   SR
        self.IA              =   IA
        self.mu              =   mu
        self.V               =   V
        self.tF              =   tF


        self.alpha           =   SA
        self.gamma           =   IA

        # if mu is non-zero, use given mu, otherwise use mu defined in the file
        if mu == 0:
            self.LMUY        =   tire['LMUY']        # Peak friction Coef 
            self.LMUX        =   tire['LMUX']        # Peak friction Coef
        else:
            self.LMUY        =   mu                  # Peak friction Coef
            self.LMUX        =   mu                  # Peak friction Coef

        self.LMUV            =   tire['LMUV']        # Speed with slip decaying function, 0 if not used
        self.AMu             =   10
        self.K               =   SR
        tire['PDX3']         =   0

        # GLobal Variables
        self.Eps             =   0.001               # TO prevent for small denom

    def validity_Check(self):
        """
        Check whether the values are within the range
        """
        warn_on = 0
        # Long Slip Range
        if np.max(self.SR) > self.tire['KPUMAX'] or np.min(self.SR) < self.tire['KPUMIN']:
            warn_on = 1
            print("WARNING: Slip ratio out of the test data range, verify validity")

        # Slip angle range
        if np.max(self.alpha) > self.tire['ALPMAX'] or np.min(self.alpha) < self.tire['ALPMIN']:
            warn_on = 1
            print("WARNING: Slip angle out of the test data range, verify validity")
        
        # Inclination angle range
        if self.gamma > self.tire['CAMMAX'] or self.gamma < self.tire['CAMMIN']:
            warn_on = 1
            print("WARNING: Camber out of test data range, verify validity")

        # Vertical force range
        if abs(np.max(self.Fz)) > self.tire['FZMAX'] or abs(np.min(self.Fz)) < self.tire['FZMIN']:
            warn_on = 1
            print("WARNING: Fz out of test data range, verify validity")

    
    def F_x(self):
        """
        Calculates F_x
        """
        #(4.E1)
        Fzp                 =   self.tire['LFZ0'] * self.Fz

        #(4.E2a)
        dfz                 =   (Fzp - self.tire['FNOMIN'] ) / self.tire['FNOMIN']

        #(4.E7)
        LStarMuX            =   self.LMUX / (1 + (self.LMUV * self.V)/ self.tire['LONGVL']) 

        #(4.E8)
        LPMuX               =   self.AMu * LStarMuX/(1 + (self.AMu - 1)* LStarMuX)

        #(4.E17)
        Shx                 =   (self.tire['PHX1'] + self.tire['PHX2'] * dfz)

        #(4.E10)
        Kx                  =   (self.K + Shx)

        #(4.E11)
        Cx                  =   self.tire['PCX1'] * self.tire['LCX']

        #(4.E13)
        muX                 =   (self.tire['PDX1'] + self.tire['PDX2'] * dfz)

        #(4.E12)
        Dx                  =   (muX * self.Fz)

        #(4.E14)
        Ex                  =  ((self.tire['PEX1'] + self.tire['PEX2'] * dfz + self.tire['PEX3'] * dfz**2) * (1 - self.tire['PEX4'])) * self.tire['LEX']

        #(4.E15)
        KxK                 =  (self.Fz * (self.tire['PKX1'] + self.tire['PKX2'] * dfz) * np.exp(self.tire['PKX3'] * dfz)) 

        #(4.E16)
        Bx                  =  (KxK / (Cx * Dx + self.Eps))

        #(4.E18)
        Svx                 =  (self.Fz * (self.tire['PVX1'] + self.tire['PVX2'] * dfz))

        #(4.E9)
        Fx                  =  (Dx * np.sin(Cx * np.arctan((Bx * Kx - Ex * (Bx * Kx - np.arctan(Bx * Kx))))) + Svx)

        return Fx 
    
    def F_y(self):
        """
        Calculates F_y
        """
        #(4.E1)
        Fzp                =  (self.tire['LFZ0'] * self.Fz)

        #(4.E2a)
        dfz                =  (Fzp - self.tire['FNOMIN'] / self.tire['FNOMIN'])

        #(4.E2b)
        #

        #(4.E3)
        alphastar          =  np.tan(self.alpha)

        #(4.E4)
        gammaStar          =  np.sin(self.gamma)

        #(4.E7)
        LStarMuY           =  self.LMUY / (1 + (self.LMUV * self.V)/ self.tire['LONGVL'])

        #(4.E8)
        LPMuY              =  self.AMu * LStarMuY / (1 + (self.AMu - 1) * LStarMuY)

        #(4.E23)
        muY                =  (self.tire['PDY1'] + self.tire['PDY2'] * dfz) * (1 - self.tire['PDY3'] * gammaStar**2) * LStarMuY

        #(4.E22)
        Dy                 =  (muY * self.Fz)

        #(4.E21)
        Cy                 =  (self.tire['PCY1'] * self.tire['LCY'])

        #(4.E28)
        Svygamma           =  self.Fz * (self.tire['PVY3'] + self.tire['PVY4'] * dfz) * gammaStar * self.tire['LKYC'] * LPMuY

        #(4.E29)
        Svy                =  self.Fz * (self.tire['PVY1'] + self.tire['PVY4'] * dfz) * self.tire['LVY'] * LPMuY + Svygamma

        #(4.E30)
        Kygamma0           =  self.Fz * (self.tire['PKY6'] + self.tire['PKY7'] * dfz) * self.tire['LKYC']

        #(4.E25)
        Kyalpha            =  self.tire['PKY1'] * Fzp * (1 - self.tire['PKY3'] * np.abs(gammaStar)) * np.sin(
                              self.tire['PKY4'] * np.arctan((self.Fz / Fzp) / ((self.tire['PKY2'] + self.tire['PKY5'] * gammaStar**2)))
                              ) * self.tire['LKY']
        
        #(4.E26)
        By                 =  Kyalpha / (Cy * Dy + self.Eps)

        #(4.E27)
        Shy                =  (self.tire['PHY1'] + self.tire['PHY2'] * dfz) * self.tire['LHY'] + ((Kygamma0 * gammaStar - Svygamma) / (Kyalpha + self.Eps))

        #(4.E20)
        alphay             =  alphastar

        #(4.E24)
        Ey                 =  (self.tire['PEY1'] + self.tire['PEY2'] * dfz) * (1 + self.tire['PEY5'] * gammaStar**2 - (self.tire['PEY3'] + self.tire['PEY4'] * gammaStar)) * self.tire['LEY']

        #(4.E19)
        Fy                 = (Dy * np.sin(Cy * np.arctan((By * alphay - Ey * (By * alphay - np.arctan(By * alphay))))) + Svy)

        return Fy

    def F_c(self):
        # Combined longitudinal and later (future scope)
        pass

    def request(self, str):
        # Settings
        self.Forces = {'Fx': self.F_x(), 'Fy': self.F_y(), 'Fc': self.F_c()}
        return self.Forces[str]