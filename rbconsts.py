from math import pi
from physconsts import *

## Rb87 constants
mRb = 1.4192261e-25 # [kg]
I = 3/2 # nuclear spin
nu_hf = 6.83468261090429 # [GHz]
gamma_D2 =2*pi*6.0659e6 ; # [rad/s]
lambda_D2 = 7.8e-7 # [m]
omega_D2 = 2*pi*c/lambda_D2
gS = 2.00023
gL = 1 # for ground state?
gI = -0.000995
Isat = (5/7)*(hbar*gamma_D2*omega_D2**3
             /(12*pi*c**2)) 
    # saturation intensity for D2 cooling W/m^2

## reduced matrix elements from Steck
D2_MatElem = 3.584e-29 # <J=1/2||er||J'=3/2> [C*m]
D1_MatElem = 2.537e-29 # <J=1/2||er||J'=1/2> [C*m]