from physconsts import *
from math import pi


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