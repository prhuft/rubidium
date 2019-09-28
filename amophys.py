"""
	AMO/QC functions!
	
	Preston Huft, 2019
"""

#### libraries
from sympy.physics.wigner import wigner_6j,wigner_3j,clebsch_gordan
from sympy import symbols,N,sympify,lambdify
from numpy import *

#### local files
from physconsts import *
from rbconsts import *

#### Electric and Magnetic Dipole Interactions 

def gF_fn(F,J,I,gJ,gI):
    """ Returns the F lande g factor """
    
    return (gJ*(F*(F+1)+J*(J+1)-I*(I+1))
            +gI*(F*(F+1)-J*(J+1)+I*(I+1)))/(2*F*(F+1))

def gJ_fn(J,L,S,gL,gS):
    """ Returns the J lande g factor """

    return (gL*(J*(J+1)+L*(L+1)-S*(S+1)) 
            +gS*(J*(J+1)-L*(L+1)+S*(S+1)))/(2*J*(J+1))

def hf_zeeman(states,gJ,gI,Bz=None,units=None):
    """ From Mark's notes for the general hf Zeeman matrix elements. Units
        are determined by UB. Could implement decorator function to change
        units."""
    global muB # Bohr magneton
    UB = symbols('U_B') # assume symbolic B field P.E. for now

    if Bz is not None:
        UB = muB*Bz# magnetic field potential energy

    if units == 'UB':
        UB = 1
        
    I,J,F,mF,FF,mFF = states
    q = 0 # assume B = Bz for now
    
    elem = 0
    if mF == mFF:
        elem += N(clebsch_gordan(F,1,FF,mF,q,mFF) \
                *sqrt(2*F+1)*(-1)**(1+J+I) \
                *(gJ*(-1)**F*sqrt(J*(J+1)*(2*J+1)) \
                *wigner_6j(J,I,F,FF,1,J) \
                +gI*(-1)**FF*sqrt(I*(I+1)*(2*I+1)) \
                *wigner_6j(I,J,F,FF,1,I))) 
        # N() is used to ensure matrix diagnolization doesn't get tripped up
    
    return UB*elem
	
def hf_coupling(F,mF,J,q,FF,mFF,JJ,I,RME=None):
    """ Returns the matrix element <F,mF,J|T_q|F',mF',J'>. 
        'RME': the reduced matrix element, e.g. the D2 line matrix
        element. If RME=None, the 
        matrix element is in units of [RME].
        
        I is the nuclear spin of the atom.
    """

    rme = 1
    if RME!=None:
        rme = RME

    ## From Mark's notes, eqs. A-50,51
    mat_elem = rme*pow(-1,F+JJ+1+I)*sqrt((2*F+1)*(2*JJ+1)) \
                *wigner_6j(J,I,F,FF,1,JJ) \
                *clebsch_gordan(1,F,FF,q,mF,mFF)
    
    return mat_elem
	
def scalar_ls(w_ab,w,I,RME=None):
    """ the scalar AC Stark shift, or light shift, seen by a two 
        level atom in an oscillating electric field.
        w_ab is the freq difference between a,b
        w is the applied light frequency
        I is the intensity of the light 
        RME is the reduced matrix element <J||er||J'>.
    """
    
    if RME is None:
        RME = 1 # the light shift is now in units of the RME
    
    return -(ee**2)*w_ab*cc(RME)*RME*I/(2*hbar*(w_ab**2-w**2))
	
#### Conversions

def radToTHz(w):
    return w/(2*pi*1e12)

def radToGHz(w):
    return w/(2*pi*1e9)

def radToMhz(w):
    return w/(2*pi*1e6) 

def radTokHz(w):
    return w/(2*pi*1e3)

def JToeV(u):
    global ee
    return u/ee

def eVToJ(u):
    global ee
    return u*ee

def GHzToeV(nu):
    global hbar
    return JToeV(2*pi*hbar*nu*1e9)

def eVToGHz(u):
    global hbar
    return eVToJ(u)/(2*pi*hbar*1e9)

	
#### Miscellaneous

def jmbasis(jlist):
    """ returns a numpy array of basis vectors {|J,mJ>} given a list of 
        J vals"""
    basis = empty(sum([2*j+1 for j in jlist]),list)
    i = 0
    for j in jlist:
        for m in range(-j,j+1):
            basis[i] = [j,m]
            i+=1 
    return flip(basis)

def cc(z):
    return conj(z)



