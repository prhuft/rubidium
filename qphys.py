"""
	Various functions for my quantum/amo simulations, etc.
	Preston Huft, Fall 2019
"""
# libraries

from numpy import *
from sympy import MatrixSymbol,MatAdd,MatMul,Identity,I,Matrix,symbols,Function
from sympy.physics.wigner import wigner_6j,wigner_3j,clebsch_gordan

# quantum physics functions

def comm(A,B):
    """ Returns the commutator of A,B: [A,B]=A.B-B.A. Assumes 'A','B' are sympy 
		matrices."""
    return (MatMul(A,B)-MatMul(B,A))
	
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
	

# syntactic sugar and convenience functions

def czeros(num):
    """ return array of complex zeros """
    c = empty(num,complex)
    for i in range(num):
        c[i] = 0+0j
    return c

def cc(z):
    """ return numpy.conj(z)"""
    return conj(z)

