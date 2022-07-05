""" 
Generic funcitonal tools for FEM 
Not object oriented

"""


import numpy as np


def insertFixedBCinModes(Qr, Tr):
    """
    Qr : (nr x nr) reduced modes
    Tr : (n x nr) reduction matrix that removed "0" fixed DOF
          such that  Mr = Tr' MM Tr
    """
    return Tr.dot(Qr)
