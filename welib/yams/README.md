# YAMS

YAMS (Yet Another Multibody Solver) has several submodules to handle the equation of motions of multibodies:

- yams symbolic kane: handles the equations of motion using symbolic calculation and Kane's method
- yams numeric recursive: handle the dynamics using a recursive velocity transformation method
- yams symbolic recursive: save as above, but uses symbolic calculations

The library also contains miscellaneous tools to handle the dynamics of flexible bodies:
- flexibility.py: tools for flexible bodies (mass matrix, stiffness matrix etc)
- sid.py: tools for the standard input data (SID) format by Walrapp
- windturbines: wrapper to compute main quantities of a wind turbine, mostly used for now to compute the numerical values necessary for "yams symbolic kane"
