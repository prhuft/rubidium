"""
	AMO/QM functions!
	
	Preston Huft, 2019
"""

#### libraries
from sympy.physics.wigner import wigner_6j,wigner_3j,clebsch_gordan
from sympy import symbols,N,sympify,lambdify
from sympy import MatrixSymbol,MatAdd,MatMul,Identity as eye,Matrix,zeros
from sympy.utilities.iterables import flatten
import numpy as np
from numpy import *
from numpy import transpose,inf
from numpy.linalg import eig

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

def hf_zeeman(states,gJ,gI,Bz=None,units='Joules'):
	""" Return Zeeman Hamiltonian in hyperfine basis |L J I F mF>. Assumes the 
		field is along the z axis, i.e. q = 0.

		'states': list-like, pair of quantum states given by [I,J,F,mF,FF,mFF]
		If Bz is none initially, the sympy free_symbol is the magnetic dipole energy,
		uB*B, not just B. 
		
		'units': 'Joules' (default), 'eV', 'UB' (units of the magnetic dipole 
		energy).
		
		From Mark's notes for the general hf Zeeman matrix elements. Units
		are determined by UB. Could implement decorator function to change
		units.
	"""
		
	## TODO: build in better unit functionality
	
	# if Bz is not None:
		# UB = uB*Bz# magnetic field potential energy; [J] by default

	# if units == 'Joules':
		# pass
	# elif units == 'eV':
	# elif units == 'UB':
		# UB = 1
		
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
		# N() is used to ensure diagnolization doesn't get tripped up

	if Bz is not None:
		elem *= uB*Bz # hmm check the sign
		if units == 'Joules':
			return elem
		elif units == 'eV':
			return JToeV(elem)
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
	
def hamiltonian(basis):
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
				H_Zz[i,j] = matrix_elem(states)
			except:
				print("Failed: %s" % states)
				print(gJ,gI,Bz)

	return H
	
def acstark_scalar(w_ab,w,I,RME=None):
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
	
#### State Evolution

def obe_derivs(y0,t,D,O,t1=np.inf,t2=np.inf):
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
	dreg = (1j*D-1/(2*t1))*reg+1j*O/2*(rgg-ree) 

	return array([drgg,dree,dreg])

#### Quantum Physics

def jmbasis(jlist):
	""" returns a numpy array of basis vectors {|J,mJ>} given a list of 
		J vals. Output is flipped w.r.t. order of j's passed in.		
	"""
	
	# TODO: make flip optional
	basis = np.empty(sum([2*j+1 for j in jlist]),list)
	i = 0
	for j in jlist:
		for m in range(-j,j+1):
			basis[i] = [j,m]
			i+=1 
	return np.flip(basis)

def comm(A,B):
	""" Returns the commutator of A,B: [A,B]=A.B-B.A. Assumes 'A','B' are sympy 
		matrices."""
		
		# TODO: extend to check type, and work for numpy matrices/ndarrays
	return (MatMul(A,B)-MatMul(B,A))

#### Syntactic sugar and convenience functions

def czeros(num):
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
	assert shape(mat)[0]==shape(mat)[1]
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


