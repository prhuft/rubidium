"""
    AMO/QM functions!
    
    Preston Huft, 2019 - 
    
    Eventually bundle stuff into classes? idk, i'm usually fine importing everything
"""

#### libraries
from sympy.physics.wigner import wigner_6j,wigner_3j,clebsch_gordan
from sympy import symbols,N,sympify,lambdify
from sympy import MatrixSymbol,MatAdd,MatMul,Identity as eye,Matrix,zeros
from sympy.utilities.iterables import flatten
import numpy as np
# np.seterr(all='raise') # raise errors, don't just print warnings
from numpy import *
from numpy.linalg import eig
from numpy.random import normal
from scipy.optimize import curve_fit
from random import random as rand
import matplotlib.pyplot as plt


#### local files
from physconsts import *
from rbconsts import *
from rbensemble import RbEnsemble as ensemble


#### Electric and Magnetic Dipole Interactions 

def gF_fn(F,J,I,gJ,gI):
    """ Returns the F lande g factor """

    if F != 0:
        return (gJ*(F*(F+1)+J*(J+1)-I*(I+1))
                +gI*(F*(F+1)-J*(J+1)+I*(I+1)))/(2*F*(F+1))
    else:
        return np.nan

def gJ_fn(J,L,S,gL,gS):
    """ Returns the J lande g factor """

    return (gL*(J*(J+1)+L*(L+1)-S*(S+1)) 
            +gS*(J*(J+1)-L*(L+1)+S*(S+1)))/(2*J*(J+1))

def hf_zeeman(states,gJ,gI,Bz=None,units='Joules'):
    """ Return Zeeman Hamiltonian matrix element in hyperfine basis |L J I F mF>. Assumes the 
        field is along the z axis, i.e. q = 0.

        'states': list-like, pair of quantum states given by [I,J,F,mF,FF,mFF]
        'Bz': the magnetic field in T. If Bz is none initially, 
        the sympy free_symbol is the magnetic dipole energy, uB*B, not just B. 
        
        'units': The energy units. 'Joules' (default), 'eV', 'UB' (units of the magnetic dipole 
        energy). Note that Bz should still always be entered in T. This merely specifies an
        optional unit conversion that will be done before returning the Hamiltonian.
        
        From Mark's notes for the general hf Zeeman matrix elements are given by
        
        UB*N(clebsch_gordan(F,1,FF,mF,q,mFF) \
                *sqrt(2*F+1)*(-1)**(1+J+I) \
                *(gJ*(-1)**F*sqrt(J*(J+1)*(2*J+1)) \
                *wigner_6j(J,I,F,FF,1,J) \
                +gI*(-1)**FF*sqrt(I*(I+1)*(2*I+1)) \
                *wigner_6j(I,J,F,FF,1,I)))
                
        Units are determined by UB, which the the magnetic field times the Bohr magneton. 
    """
        
    ## TODO: build in better unit functionality or remove option to choose units
        
    I,J,F,mF,FF,mFF = states
    q = 0 # assume B = Bz for now

    elem = 0
    if mF == mFF: # q=0 => selection rule mF = mF'
        elem += N(clebsch_gordan(F,1,FF,mF,q,mFF) \
                *sqrt(2*F+1)*(-1)**(1+J+I) \
                *(gJ*(-1)**F*sqrt(J*(J+1)*(2*J+1)) \
                *wigner_6j(J,I,F,FF,1,J) \
                +gI*(-1)**FF*sqrt(I*(I+1)*(2*I+1)) \
                *wigner_6j(I,J,F,FF,1,I)))
        # N() is used to ensure diagnolization doesn't get tripped up

    if Bz is not None:
        elem *= uB*Bz # hmm check the sign
        if units == 'Joules':
            return elem
        elif units == 'eV':
            return JToeV(elem)
        elif units == 'GHz':
            return eVToGHz(JToeV(elem)) 
        else:
            print(f"invalid unit [{units}] for non-zero Bz. Result in 'J'.")
            return elem
    else: 
        if units == 'UB':
            return elem
        else:
            UB = symbols('U_B') # symbolic B field P.E. for now
            elem*= UB
            return elem
            
def hf_coupling(F,mF,J,q,FF,mFF,JJ,I):
    """
    Returns the constant relating matrix element <J,mJ|T_q|J',mJ'> to the 
    reduced fine structure matrix element. 
    
    Args:
        I is the nuclear spin of the atom.
    Returns:
        'coupling': 
        ## From Mark's notes, eqs. A-50,51. Also see Steck Rb datasheet.
        mat_elem = rme*pow(-1,F+JJ+1+I)*sqrt((2*F+1)*(2*JJ+1)) \
                   *wigner_6j(J,I,F,FF,1,JJ) \
                   *clebsch_gordan(1,F,FF,q,mF,mFF)
    """

    ## From Mark's notes, eqs. A-50,51
    coupling = pow(-1,F+JJ+1+I)*sqrt((2*F+1)*(2*JJ+1)) \
                *wigner_6j(J,I,F,FF,1,JJ) \
                *clebsch_gordan(1,F,FF,q,mF,mFF)

    return coupling
    
def hf_reduced_f(F,J,FF,JJ,I):
    """
    Returns the coupling constant relating the reduced matrix element 
    in the hyperfine structure basis to that in the fine structure basis
    """

    ## Mark's notes, eq. A-51 
    return (-1)**(1+I+JJ+F)*sqrt((2*F+1)*(2*FF+1))*wigner_6j(J,I,F,FF,1,JJ)
    
def f_reduced_ls(J,L,JJ,LL):
    """
    Returns the coupling constant relating the reduced matrix element 
    in the fine structure basis to that in the L,S basis
    """

    ## Mark's notes, eq. A-51
    S = 1/2
    return (-1)**(1+S+LL+J)*sqrt((2*J+1)*(2*JJ+1))*wigner_6j(L,I,J,JJ,1,LL)

def f_coupling(L,J,mJ,q,LL,JJ,mJJ,I):
    """ 
    Returns the constant relating matrix element <J,mJ|T_q|J',mJ'> to the 
    reduced fine structure matrix element. 
    
    Args:
        'RME': the reduced matrix element <alpha;J||r||alpha'J'> with alpha 
        including quantum numbers not relevant to the coupling, e.g. n. 
        If RME=None, the matrix element is in units of [RME].
        
        I is the nuclear spin of the atom.
    Returns:
        'mat_elem': 
        ## From Mark's notes, eqs. A-50,51. Also see Steck Rb datasheet.
        mat_elem = rme*pow(-1,F+JJ+1+I)*sqrt((2*F+1)*(2*JJ+1)) \
                   *wigner_6j(J,I,F,FF,1,JJ) \
                   *clebsch_gordan(1,F,FF,q,mF,mFF)
    """

    ## From Mark's notes, eqs. A-50,51
    mat_elem = clebsch_gordan(1,J,JJ,q,mJ,mJJ)

    return mat_elem
    
def coupling_f_to_j(jr,mjr,jp,fg,mg,q1,q2,deltap,I=I):
    """
    Coupling constant between state (jr,mjr) and (fg,mg) via a ladder transition
    through (jp, mjp) with two photons with projections q1,q2. Minho's thesis
    section 2.4.
    """
    coupling = 0
    
    # little_c lambda I,j,f,jj,ff = (-1)**(1+I+f+jj)*sqrt(2*f+1)*wigner_6j(j,1,f,ff,1,jj)
    pass
    
# def alpha0_hf(state, omega, nlist, atom, I, hf_states=hf_states, printterms=False):
#     """
#     Returns the scalar polarizability for a hyperfine level |n,j,l,f>
#     for field with angular frequency omega. nlist is a list or iterable
#     of the n values to use in the perturbative sum.
    
#     This function uses the hyperfine level frequencies for levels up to 
#     5p3/2 from hf_states in rbconsts.py, which are taken from Steck.
    
#     Args:
#         'state': list of quantum numbers n,l,j,f
#         'atom': an Atom object from the Alkali Rydberg Calculator module
#     Returns:
#         'alpha0': scalar hyperfine polarizability in S.I. units.
#     """
#     alpha0 = 0
#     terms = 0
    
#     n_a, l_a, j_a, f_a = state

#     # program:
#     # - loop over levels from nmin to nmax
#     #  - loop over l
#     #   - loop over j
#     #    - loop over f

#     for n_b in nlist:
#         for l_b in range(n_b): # runs through n_b - 1, inclusive
#             for j_b in j3_from_j1j2(s, l_b): 

#                 # triangle rule and -allowed conditions
#                 if abs(j_b - j_a) <= 1 and abs(l_b - l_a) == 1:

#                     try: 
#                         _ = hf_states[n_b][l_b][j_b]
#                         freq_from_dict = True
#                         # calculate from dictionary for each f later. 
#                         # assume all possible f for a given j are included in the dictionary

#                     except KeyError: # ignore the hf shift
#                         freq_from_dict = False
#                         omega_ba = 2*pi*((eVToGHz(atom.getEnergy(n_b, l_b, j_b) 
#                                      - atom.getEnergy(n_a, l_a, j_a))*1e9))

#                     matelemJ = atom.getReducedMatrixElementJ(n_a, l_a, j_a, n_b, l_b, j_b)*ee*a0
#     #                 print(f"< n={n_a}, l={l_a}, j={j_a} | x | n'={n_b}, l'={l_b}, j'={j_b} >")

#                     for f_b in j3_from_j1j2(I, j_b):

#                         if freq_from_dict:
#                             omega_ba = 2*pi*(hf_states[n_b][l_b][j_b][f_b] - hf_states[n_a][l_a][j_a][f_a])*1e9

#                         matelemF = hf_reduced_f(f_a,j_a,f_b,j_b,I)*matelemJ
#                         #*2*(f_a+1)
#                         summand = (2/(3*hbar*2*(f_a+1)))*omega_ba*abs(matelemF)**2/(omega_ba**2 - omega**2) 
#                         alpha0 += summand

#                         terms += 1
#                         # if printterms:
#                         #     print(f"alpha0 ~= {alpha0/(4*pi*e0*1e-30)} [Ang.^3], {terms} terms in sum")

#     return float(alpha0)
    
    # # this fudge factor seems to be necessary to get my values to agree with the wiki and 
    # # with the Delaware database, but previously I don't think I needed this. It's likely
    # # that my hf_reduced_f has changed.
    # return float(np.sqrt(2*j_a+1)*alpha0)

def alphaK_nJ(state, omega, nlist, atom, K):
    """
    Returns the reduced dynamic scalar polarizability for a fine structure state |n,J>
    for field with angular frequency omega, i.e., eq. (11) of the Rauschenbeutel paper below. nlist is a list or iterable of the n values to use in the perturbative sum. 

    References:
    ...Rauschenbeutel, "Dynamical polarizability of atoms in arbitrary light fields:general theory and application to cesium"
    
    Args:
        'state': list of quantum numbers n,l,j,f
        'atom': an Atom object from the Alkali Rydberg Calculator module
    Returns:
        'alpha1': vector hyperfine polarizability in S.I. units.
    """
    
    alphaK = 0.0
    terms = 0
    
    n_a, l_a, j_a, f_a = state

    # program: compute vector alphaJ
    # - loop over levels from nmin to nmax
    #  - loop over l
    #   - loop over j
    # multiply by F dependent factor to get alphaF
    
    for n_b in nlist:
        for l_b in range(n_b): # runs through n_b - 1, inclusive
            for j_b in j3_from_j1j2(s, l_b): 
                
                # triangle rule and dipole-allowed conditions
                if abs(j_b - j_a) <= 1 and abs(l_b - l_a) == 1:

                    omega_ba = 2*pi*abs((eVToGHz(atom.getEnergy(n_b, l_b, j_b) 
                                 - atom.getEnergy(n_a, l_a, j_a))*1e9))

                    Gamma_ba = 2*pi/atom.getStateLifetime(n_b, l_b, j_b) # as an approximation, ignore the branching into a particular state

                    matelemJ = atom.getReducedMatrixElementJ(n_a, l_a, j_a, n_b, l_b, j_b)*ee*a0

                    # print(f"< n={n_a}, l={l_a}, j={j_a} | x | n'={n_b}, l'={l_b}, j'={j_b} > = ",matelemJ)

                    frequency_part = np.real(1/(omega_ba - omega - 1j*Gamma_ba/2) - 
                                               (-1)**K/(omega_ba + omega + 1j*Gamma_ba/2))/hbar

                    # complex is used to reduce sympy expressions
                    alphaK += complex((-1)**j_b*wigner_6j(1,K,1,j_a,j_b,j_a)*abs(matelemJ)**2)*frequency_part
    
    # this is the complete alphaK
    alphaK *= (-1)**(K+j_a+1)*np.sqrt(2*K+1)
    
    return np.real(alphaK)


def alpha0_hf(state, omega, nlist, atom, I):
    """
    Returns the alpha_s_nJF, dynamic scalar polarizability for a hyperfine state |n,j,l,f>
    for field with angular frequency omega. nlist is a list or iterable
    of the n values to use in the perturbative sum. This definition 
    does not depend on mf. 

    References:
    ...Rauschenbeutel, "Dynamical polarizability of atoms in arbitrary light fields:general theory and application to cesium"
    
    Args:
        'state': list of quantum numbers n,l,j,f
        'atom': an Atom object from the Alkali Rydberg Calculator module
    Returns:
        'alpha_s_nJF': vector hyperfine polarizability in S.I. units.
    """

    n_a, l_a, j_a, f_a = state
    
    alpha_s_nJ = alphaK_nJ(state, omega, nlist, atom, K=0)/np.sqrt(3*(2*j_a+1))

    alpha_s_nJF = alpha_s_nJ
    return np.real(alpha_s_nJF)


def alpha1_hf(state, omega, nlist, atom, I):
    """
    Returns the dynamic vector polarizability for a hyperfine level |n,j,l,f>
    for field with angular frequency omega. nlist is a list or iterable
    of the n values to use in the perturbative sum. This definition 
    does not depend on mf. 

    References:
    ...Anderson, "Measurement and extinction of vector light shifts using interferometry of spinor condensates"
    ...Rauschenbeutel, "Dynamical polarizability of atoms in arbitrary light fields:general theory and application to cesium"
    
    Args:
        'state': list of quantum numbers n,l,j,f
        'atom': an Atom object from the Alkali Rydberg Calculator module
    Returns:
        'alpha_v_nJF': vector hyperfine polarizability in S.I. units.
    """

    n_a, l_a, j_a, f_a = state
    
    alpha_1_nJ = alphaK_nJ(state, omega, nlist, atom, K=1)

    alpha_v_nJF = alpha_1_nJ*complex((-1)**(j_a+I+f_a)*np.sqrt(2*f_a*(2*f_a+1)/(f_a+1))*wigner_6j(f_a, 1, f_a, j_a, I, j_a))
    
    return np.real(alpha_v_nJF)

def alpha2_hf(state, omega, nlist, atom, I):
    """
    Returns the dynamic tensor polarizability for a hyperfine level |n,j,l,f>
    for field with angular frequency omega. nlist is a list or iterable
    of the n values to use in the perturbative sum. This definition 
    does not depend on mf.

    References:
    ...Rauschenbeutel, "Dynamical polarizability of atoms in arbitrary light fields:general theory and application to cesium"
    
    Args:
        'state': list of quantum numbers n,l,j,f
        'atom': an Atom object from the Alkali Rydberg Calculator module
    Returns:
        'alpha_T_nJF': vector hyperfine polarizability in S.I. units.
    """

    n_a, l_a, j_a, f_a = state
    
    alpha_2_nJ = alphaK_nJ(state, omega, nlist, atom, K=2)

    alpha_T_nJF = (alpha_2_nJ*complex(-1*(-1)**(j_a+I+f_a)*np.sqrt(2*f_a*(2*f_a-1)*(2*f_a+1)/(2*(f_a+1)*(2*f_a+3)))*wigner_6j(f_a, 2, f_a, j_a, I, j_a)))
    
    return np.real(alpha_T_nJF)
    
    
def hamiltonian_z1(basis,gI,gJ,Bz=None,I=1.5,J=.5,units='Joules'):
    """ returns the hf zeeman hamiltonian Z1 for a provided basis.
        'Bz': the field strength [T]. None by default, then hamiltonian
        can be lambdified for the field energy UB = - muB*Bz. 
        
        'units': 'Joules','eV', more depending on matrix elem call 
        
        For now, assumes 87Rb (I = 3/2) ground states (J=1/2)
    """
    #TODO: make general. could specify a basis format by string, e.g. 
    # ['I', 'J', 'F', 'mF']

    dim = len(basis)
    H_Zz = np.empty((dim,dim),object)
        
    for i,state_i in enumerate(basis):
        F,mF = state_i
        for j,state_j in enumerate(basis):
            FF,mFF = state_j
            states = [I,J,F,mF,FF,mFF]
            try:
                H_Zz[i,j] = hf_zeeman(states,gJ,gI,Bz=Bz,units=units)
            except:
                print("Failed: %s" % states)
                print(gJ,gI,Bz)

    return H_Zz
    
def hamiltonian(basis,mat_elem):
    """ 
        returns a hamiltonian in the given basis, whose elements are to 
        be specified by a decorator function. 
    """
    #TODO: make general. could specify a basis format by string, e.g. 
    # ['I', 'J', 'F', 'mF']

    dim = len(basis)
    H = np.empty((dim,dim),object)
        
    for i,state_i in enumerate(basis):
        F,mF = state_i
        for j,state_j in enumerate(basis):
            FF,mFF = state_j
            states = [I,J,F,mF,FF,mFF]
            try:
                H[i,j] = mat_elem(states)
            except:
                print("Failed: %s" % states)
                print(gJ,gI,Bz)

    return H
    
def acstark_twolevel(O,D):
    """ The two level atom AC Stark shift from diagnolization of the full Hamiltonian,
        assuming the detuning D << E2-E1, and making the RWA. 
        
        Uac = -O**2/(4*D)
    """
    
    return -O**2/(4*D)

    if RME is None:
        RME = 1 # the light shift is now in units of the RME

    return -(ee**2)*w_ab*cc(RME)*RME*I/(2*hbar*(w_ab**2-w**2))
    
def rabi(I0,matelem):
    """ the rabi frequency for a beam of intensity I0 coupling states |e>,|g>
        such that matelem = <e|q x|g>. May need to include a coupling const. 
        as prefactor on matelem. 
    """
    return sqrt(2*I0/(c*e0))*matelem/hbar
    
#### Atomic State Evolution

def obe_derivs(y0,t,D,O,phi=0,t1=np.inf,t2=np.inf):
# def derivs(t,y0,params):
    """ Returns RHS of optical bloch eqs for current values at time t,
    
        drgg = ree/t1 - 1j/2*(O*cc(reg)-cc(O)*reg) 
        dree = -ree/t1 + 1j/2*(O*cc(reg)-cc(O)*reg)
        dreg = (1j*D-1/(2*t1))*reg+1j*O/2*(rgg-ree) # = cc(drge)
    
        'y0': [rgg,ree,reg]
        't': time
        'D': detuning
        'O': Rabi frequency
        't1': state lifetime
        't2': coherence
    """
    
    rgg,ree,reg = y0
    
    # time derivatives of density op elements
    curl = 1j/2*(O*cc(reg)-cc(O)*reg) 
    drgg = ree/t1 - curl 
    dree = -ree/t1 + curl
    dreg = (1j*D-1/(2*t1))*reg+1j*O/2*(rgg-ree) # actually reg tilda
    
    return array([drgg,dree,dreg])

#### Various classes

# TODO: make optical trap an ABC class from which dipole_trap and projected_dark_trap inherit

class dipole_trap:
    
    def __init__(self,lmbda,wx,Tdepth,Tatom,wy=None,m=mRb,units=1):
        """ A dipole trap object with the beams potential and distribution of 
            atoms specified by Tatom. 
            'wx': x beam waist in focal plane (z=0)
            'wy': Assumed equal to wx by default
            'Tdepth'
            'Tatom'
        """
        self.wx = wx
        self.Tdepth = Tdepth
        self.T = Tatom
        self.m = m
        self.units = units
        if wy is None:
            self.wy = wx
        else:
            self.wy = wy
        
        # FORT and atom parameter stuff
        self.umax = kB*self.Tdepth # the maximum FORT depth
        self.lmbda = lmbda # the trap wavelength [m]

        self.zR = pi*wx**2/self.lmbda
        self.omega_r = (2/sqrt((self.wx**2+self.wy**2)/2))*sqrt(kB*self.Tdepth/m) # radial trap frequency 
        self.omega_z = (1/self.zR)*sqrt(2*kB*self.Tdepth/m) # axial trap frequency
        # print(f"omega_r = {self.omega_r*1e-3:.3f} kHz, omega_z = {self.omega_z*1e-3:.3f} kHz")

    def U(self,x,y,z):
        """ the potential energy as a function of space in the dipole trap.
            ASTIGMATISM not correctly represented!!
        """
        zR = self.zR
        wx = self.wx
        wy = self.wy
        ww = (1+z**2/zR**2) 
        umax = self.umax
        return -umax*exp(-2*x**2/(wx**2*ww)-2*y**2/(wy**2*ww))/ww 
    
    def Uforce(self):
        """
        returns list of force components lambda x,y,z [fx,fy,fz]
        """
        zR = self.zR
        wx = self.wx
        wy = self.wy
        w0 = sqrt((wx**2+wy**2)/2)
        w = lambda z: w0*sqrt(1+(z/zR)**2)
        umax = self.umax
                
        fx = lambda x,y,z: -umax*(4*exp(-2*(x**2+y**2)/w(z)**2)*x
                             /(w(z)**2))
        fy = lambda x,y,z: -umax*(4*exp(-2*(x**2+y**2)/w(z)**2)*y
                             /(w(z)**2))
        fz = lambda x,y,z: umax*(4*exp(-2*(x**2+y**2)/w(z)**2)*(x**2+y**2)*z
                             /(w(z)**2*(1+(z/zR)**2)*zR**2)
                           - 2*exp(-2*(x**2+y**2)/w(z)**2)*z
                             /(w(z)**2*(zR/w0)**2))
        
        return fx,fy,fz
        
    
    def xdist(self,events=None,plane=None):
        """ velocity component distributions """
        # Grainger group method
        omega_r = self.omega_r
        omega_z = self.omega_z
        T = self.T
        m = self.m
        dx = dy = sqrt(kB*T/(m*omega_r**2))
        dz = sqrt(kB*T/(m*omega_z**2))
        zlist = normal(0,dz,size=events)
        xlist = normal(0,dx,size=events)
        ylist = normal(0,dy,size=events)
        
        if plane == 'xz':
            return xlist,zlist
        else:
            if events is None:
                return xlist[0],ylist[0],zlist[0]
            return xlist,ylist,zlist
    
    def vdist(self,events):
        """ maxwell boltzmann speeds """
        umax = self.umax
        atoms = ensemble(self.T,m=self.m,units=self.units)
        
        vlist = atoms.sampling_maxboltzv(events,[0,2/self.units]) # speeds

        vxlist = empty(events)
        vylist = empty(events)
        vzlist = empty(events)
        
        for i in range(events):
            ex = 2*rand()-1
            ey = 2*rand()-1
            ez = 2*rand()-1
            v = vlist[i]
            A = sqrt(ex**2+ey**2+ez**2)
            vxlist[i] = ex*v/A
            vylist[i] = ey*v/A
            vzlist[i] = ez*v/A
            
#         vlist = array([sqrt(vx**2+vy**2+vx**2) for vx,vy,vz in zip(vxlist,vylist,vzlist)])
#         plt.hist(vlist, 50, density=True) # show maxwell boltzmann speed dist
#         plt.xlabel('v')
#         plt.ylabel('occurences')
#         plt.show()

        return vxlist,vylist,vzlist
    
    def distplot(self,events, atoms=True, sx=1.5, sz=1.5, color="bone"):
        """ show atoms in FORT in z = 0 plane before drop and recapture """
        mu = self.units
        wx = self.wx
        zR = self.zR
        print(f"zr={zR/mu:.2f} [um], wx={wx/mu:.2f} [um]")
        
        xlist,ylist = self.xdist(events,plane='xz') # positions in [m]
        
        xpts = linspace(-sx*wx,sx*wx,100)
        zpts = linspace(-sz*zR,sz*zR,100)
        xx,zz = meshgrid(xpts,zpts)
        fpts = -self.U(xx,0,zz) # the fort intensity eval'd on the meshgrid
        
        cmap = plt.cm.get_cmap(color)
        
        plt.contourf(xpts/mu,zpts/mu,fpts, cmap=cmap)
        if atoms is True: # otherwise, just a dipole trap plot  
            plt.scatter(xlist/mu,ylist/mu,color='red')
        plt.xlabel("x")
        plt.ylabel("z")
#         plt.axes().set_aspect('equal')
        plt.show() 

    def drop_recap(self,tlist,T=None,events=None,base_retention=None,
            progress=False):
        """ Procedure for simulating a release ("drop") and recapture experiment
            to deduce the temperature of actual atoms in such an experiment. 
        
            Based on code by Mark, with some corrections
            'wx': waist
            'Tdepth': FORT temperature depth
            'T': atom temp
            'tmax': max time in units us
            'steps': number of FORT drop outs
            'events': number of release-recapture events per data pt
            'wy': optional waist for eliptical FORT 
        """
      
        Tdepth = self.Tdepth
        wx = self.wx
        umax = self.umax
        zR = self.zR
        tlist = 1e-6*tlist

        if T is None:
            T = self.T
        if events is None:
            events = 2000
        if base_retention is None:
            base_retention = 1 # the retention baseline with no fort drop

        retention = empty(len(tlist))
   
        xlist,ylist,zlist = self.xdist(events)
        vzlist,vxlist,vylist = self.vdist(events)

        for j,t in enumerate(tlist):

            escape = 0 
            nhot = 0 # this is an untrapped atom

            for i in range(events):     
                hot = 0
                KE = .5*mRb*((vxlist[i]-g*t)**2+vylist[i]**2
                              +vzlist[i]**2)
                PE0 = self.U(xlist[i],ylist[i],zlist[i])
                PE = self.U(xlist[i]+t*vxlist[i]+.5*g*(t)**2,
                       ylist[i]+t*vylist[i],
                       zlist[i]+t*vzlist[i])

                if KE + PE0 > 0:
                    hot = 1
                nhot += hot
                if KE + PE > 0:
                    escape += 1-hot
            retention[j] = base_retention*(1 - escape/events)
            
            if progress is not False:
                if j % 10 == 0:
                    print(f"timestep {j}: t = {t*1e6:.0f} [us], ret = {retention[j]:.2f}")     
        
        print(f"finished. T={T*1e6} [uK], r = {base_retention}")
        return tlist*1e6,retention
    
    def curvefit(self,tdata,rdata):
        """ For using the release_recap procedure to fit real data.
            
            tdata: time pts from data
            rdata: retention pts data
        """
        def f(tlist,T,r):
            """
                tlist: times [s]
                T: atom temp [K]
                r: baseline retention in (0,1]
            """
            
            t,y = self.drop_recap(tlist,T=T,base_retention=r)
            return y
        
        p0=[self.T,1] # scale up the temp 
        popt,pcov = curve_fit(f,tdata,rdata,p0,
                              absolute_sigma=True)
        
        Topt,ropt = popt # the optimum parameters
        
        print(f"T_fit = {Topt*1e6:.0f}")
        tlist,ret=self.drop_recap(tdata,T=Topt,base_retention=ropt)
        
        plt.plot(tlist,ret,label=f"T_fit = {Topt*1e6:.0f}")
        plt.scatter(tdata,rdata,label='real live data',color='r')
        plt.show()
        
        return tlist,ret
        
        
class projected_dark_trap:
    """
    a trap with profile U ~ |1 - E_gaussian|^2
    """
    
    def __init__(self,lmbda,wx,Tdepth,Tatom,m,units=1):
        """ A dipole trap object with the beams potential and distribution of 
            atoms specified by Tatom. 
            'wx': x beam waist in focal plane (z=0); wy = wx
            'Tdepth': [K]
            'Tatom': [K]
            'm': atomic mass
            units: length units. 1 for meters, 1e-3 for mm, etc. to be passed to the ensemble object
        """
        self.wx = wx
        self.wy = wx
        self.Tdepth = Tdepth
        self.T = Tatom
        self.m = m
        self.units = units
        
        # FORT and atom parameter stuff
        self.umax = kB*self.Tdepth # the maximum FORT depth
        self.lmbda = lmbda # the trap wavelength [m]

        self.zR = pi*wx**2/self.lmbda
#         self.omega_r = (1/sqrt((self.wx**2+self.wy**2)/2))*sqrt(2*kB*self.Tdepth/m) # radial trap frequency 
        self.omega_z = (1/self.zR)*sqrt(2*kB*self.Tdepth/m) # axial trap frequency
        # print(f"omega_r = {self.omega_r*1e-3:.3f} kHz, omega_z = {self.omega_z*1e-3:.3f} kHz")

    def U(self,x,y,z):
        """ 
        The optical potential evaluated at x,y,z
        """
        zR = self.zR
        w0 = self.wx
        w = w0*sqrt(1+(z/zR)**2) 
        rr = x**2+y**2 # rho squared
        umax = self.umax
        lmbda = self.lmbda
        
        eta = arctan2(z,zR)
        R = z*(1+(zR/z)**2)
        
        
        # this occasionally gives a RuntimeWarning of an underflow in the exp()
        # the reason for not simply having the if elif else is that the conditions
        # don't evaluate properly for arrays, but the only time arrays are passed
        # is when we're plotting U over a mesh, and we shouldn't ever get overflows
        # when doing that.
        try:
            Irho = exp(-2*rr/w**2)
        except OverflowError:
            if log10(2*rr) - 2*log10(w) == inf:
                arg = inf
            elif log10(2*rr) - 2*log10(w) == -inf:
                arg = -inf
            else:
                logval = log10(2*rr) - 2*log10(w)
                arg = exp(logval)
            Irho = exp(-arg)
        
        # the minus 1 at the end just adds an arbitrary offset,
        # such that the minimum is -umax, and umax(x=y=z=inf) = 0
        uxyz = umax*((1 - 2*(w0/w)*cos(eta-pi*rr/(R*lmbda))
                     *exp(-rr/w**2)
                     +(w0/w)**2*Irho) - 1) 
                   
        return uxyz
    
    def Uforce(self):
        """
        returns list of force components lambda x,y,z [fx,fy,fz]
        """
        zR = self.zR
        lmbda = self.lmbda
        w0 = self.wx
        w = lambda z: w0*sqrt(1+(z/zR)**2)
        umax = self.umax
        
        try:
            fx = lambda x,y,z: -umax*(4/(w0**2*(z**2+zR**2)**2)
                *exp(-2*(x**2+y**2)*zR**2/(w0**2*(z**2+zR**2)))
                *x*(-zR**4+(1/lmbda)*exp((x**2+y**2)*zR**2/(w0**2*(z**2+zR**2)))
                    *sqrt(1+(z/zR)**2)*zR**2*(
                        zR**2*lmbda*cos(pi*(x**2+y**2)*z/(lmbda*(z**2+zR**2))-arctan2(z,zR))
                        +pi*w0**2*z*sin(pi*(x**2+y**2)*z/(lmbda*(z**2+zR**2))-arctan2(z,zR))
                        )
                    )
                )/self.units
                    
            fy = lambda x,y,z: -umax*(4/(w0**2*(z**2+zR**2)**2)
                *exp(-2*(x**2+y**2)*zR**2/(w0**2*(z**2+zR**2)))
                *y*(-zR**4+(1/lmbda)*exp((x**2+y**2)*zR**2/(w0**2*(z**2+zR**2)))
                    *sqrt(1+(z/zR)**2)*zR**2*(
                        zR**2*lmbda*cos(pi*(x**2+y**2)*z/(lmbda*(z**2+zR**2))-arctan2(z,zR))
                        +pi*w0**2*z*sin(pi*(x**2+y**2)*z/(lmbda*(z**2+zR**2))-arctan2(z,zR))
                        )
                    )
                )/self.units
                
            fz = lambda x,y,z: umax*(2/(w0**2*(z**2+zR**2)**3)
                *exp(-2*(x**2+y**2)*zR**2/(w0**2*(z**2+zR**2)))*zR**2*(
                -2*(x**2+y**2)*z*zR**2+w0**2*z*(z**2+zR**2)
                +(1/lmbda)*exp((x**2+y**2)*zR**2/(w0**2*(z**2+zR**2)))
                *sqrt(1+(z/zR)**2)*(-z*(-2*(x**2+y**2)*zR**2+w0**2
                                             *(z**2+zR**2))*lmbda
                                          *cos(pi*(x**2+y**2)*z/(lmbda*(z**2+zR**2))-arctan2(z,zR))
                                          +w0**2*(pi*(x**2+y**2)*(z-zR)*(z+zR)
                                          +zR*(z**2+zR**2)*lmbda)
                        *sin(pi*(x**2+y**2)*z/(lmbda*(z**2+zR**2))-arctan2(z,zR))
                    )
                )
            )/self.units
            
            
        except RuntimeWarning as e:
            print(f"Warning: {e}")
            print(f"x={x}, y={y}, z={z}")
    
        
        return fx,fy,fz
        
    
    def xdist(self,events=None,plane=None):
        """ 
        position and velocity component distributions 
        
        the trap is harmonic along z, rho=0 but quartic along rho,z=0.
        the quartic radial distribution is derived in Isenhower et al
        "Atom trapping in an interferometrically generated bottle beam trap"
        """
        
        omega_z = self.omega_z
        T = self.T # atom temp
        Tdepth = self.Tdepth
        m = self.m
        w0 = self.wx
        
        dx = dy = (2/sqrt(3))*w0*(T/Tdepth)**(1/4)
        dz = sqrt(kB*T/(m*omega_z**2))
        zlist = normal(0,dz,size=events)
        xlist = normal(0,dx,size=events)
        ylist = normal(0,dy,size=events)
        
        if plane == 'xz':
            return xlist,zlist
        else:
            if events is None:
                return xlist[0],ylist[0],zlist[0]
            return xlist,ylist,zlist
    
    def vdist(self,events):
        """ maxwell boltzmann speeds """
        umax = self.umax
        atoms = ensemble(self.T,m=self.m,units=self.units)
        
        vlist = atoms.sampling_maxboltzv(events,[0,2/self.units]) # speeds

        vxlist = empty(events)
        vylist = empty(events)
        vzlist = empty(events)
        
        for i in range(events):
            ex = 2*rand()-1
            ey = 2*rand()-1
            ez = 2*rand()-1
            v = vlist[i]
            A = sqrt(ex**2+ey**2+ez**2)
            vxlist[i] = ex*v/A
            vylist[i] = ey*v/A
            vzlist[i] = ez*v/A

        return vxlist,vylist,vzlist
    
    def distplot(self,events, atoms=True, sx=1.5, sz=1.5, color="bone",normalize=True):
        """ show atoms in FORT in z = 0 plane before drop and recapture """
        mu = 1e-6/self.units
        wx = self.wx
        zR = self.zR
        print(f"zr={zR/mu:.2f} [um], wx={wx/mu:.2f} [um]")
        
        xlist,ylist = self.xdist(events,plane='xz') # positions in [m]
        
        xpts = linspace(-sx*wx,sx*wx,100)
        zpts = linspace(-sz*zR,sz*zR,100)
        xx,zz = meshgrid(xpts,zpts)
        if normalize:
            norm = abs(self.umax)
        else:
            norm = 1
        fpts = self.U(xx,0,zz)/norm # the fort intensity eval'd on the meshgrid
        
        cmap = plt.cm.get_cmap(color)
        
        plt.contourf(xpts/mu,zpts/mu,fpts, cmap=cmap)
        plt.colorbar()
        if atoms is True: # otherwise, just a dipole trap plot  
            plt.scatter(xlist/mu,ylist/mu,color='red')
        plt.xlabel("x")
        plt.ylabel("z")
#         plt.axes().set_aspect('equal')
        plt.show()

#### Optics

class gaussian_beam:
    
    ## TODO: figure out how to have methods not tied to object ref, 
    ## which can be called individually, and also ref each other
    
    
    def intensity(x,y,z,lmbda,wx,I0,wy=None,z_offset=0):
        if wy is None:
            wy = wx
            
        zRx = z_rayleigh(lmbda,wx)
        zRy = z_rayleigh(lmbda,wy)

        wzx = (1+(z/zRx)**2)
        wzy = (1+(z/zRy)**2)
        return I0*exp(-2*x**2/(wzx*(wx**2*wzx))-2*y**2/(wzy*(wy**2*wzy)))	
    
    
def z_rayleigh(lmbda,w0):
        return pi*w0**2/lmbda
        
#### Quantum Physics

def jmbasis(jlist, flip=False):
    """ returns a numpy array of basis vectors {|J,mJ>} given a list of 
        J vals. Output is flipped w.r.t. order of j's passed in.		
        
        flip: order by j descending if True. False by default.
    """
    
    # TODO: make flip optional
    numstates = sum([2*j+1 for j in jlist])
    assert numstates == int(numstates), "j values must be half integral or integral"
    basis = np.empty(int(numstates),dtype=list)
    i = 0
    for j in jlist:
        for m in np.arange(-j,j+1):
            basis[i] = [j,m]
            i+=1 
    if flip:
        return np.flip(basis)
    else:
        return basis

def j3_from_j1j2(j1,j2):
    """
    returns list of possible angular momentum values j3 from combining j1,j2
    """
    
    return [abs(j1-j2)+i for i in range(int(2*min(j1,j2))+1)]

def comm(A,B):
    """ Returns the commutator of A,B: [A,B]=A.B-B.A. Assumes 'A','B' are sympy 
        matrices."""
        
        # TODO: extend to check type, and work for numpy matrices/ndarrays
    return (MatMul(A,B)-MatMul(B,A))

#### Syntactic sugar and convenience functions

def czeros(num): # numpy can already do this with zeros(num,complex)
    """ return array of complex zeros """
    c = empty(num,complex)
    for i in range(num):
        c[i] = 0+0j
    return c

def cc(z):
    """ return numpy.conj(z)"""
    return conj(z)
    
def diagonal(mat):
    """ 
        diagonalize mat, but with eigenvalues ordered from largest (0,0) to
        smallest (N,N) for mat NxN. Warning: this casts the output to a np 
        array for now. 
        
        'mat' assumed to be sympy matrix. 
    """
    # TODO: extend to support numpy square ndarray diagnolization as well
    assert shape(mat)[0]==shape(mat)[1], "the matrix is not square"
    dim = shape(mat)[0]
    P,D = mat.diagonalize()
    D_copy = copy(D) # guarantees the return type the same as input type
    
    try:
        eigs = flip(sort([D[i,i] for i in range(dim)]))
        for i in range(dim):
            D_copy[i,i] = eigs[i]
        return Matrix(D_copy)
    except:
        print("make sure that the matrix is numeric")
        return D
        
#### Conversions
    
def cgsToSI(alpha):
    return 4*pi*e0*1e-6*alpha


def auToSI(alpha):
    """ for polarizability by default """
    return alpha/1.64877727436e-41


def radToTHz(w):
    return w/(2*pi*1e12)


def radToGHz(w):
    return w/(2*pi*1e9)


def radToMhz(w):
    return w/(2*pi*1e6) 


def radTokHz(w):
    return w/(2*pi*1e3)


def JToeV(u):
    return u/ee


def eVToJ(u):
    return u*ee


def GHzToeV(nu):
    return JToeV(2*pi*hbar*nu*1e9)
    

def eVToGHz(u):
    return eVToJ(u)/(2*pi*hbar*1e9)
    
    
def JToGHz(u):
    return eVToGHz(JToeV(u))


def GHzToJ(u):
    return eVToJ(GHzToeV(u))

