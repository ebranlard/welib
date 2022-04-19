""" 
Example to extract rigid body mass and inertia of an entire turbine 
or of the rotor-nacelle-assembly (RNA) using OpenFAST linearization (or using YAMS).

See welib.fast.extract.py for more options


"""

import os
from welib.fast.extract import extractRigidBodyInertia, extractRNAInertia

if __name__ == '__main__':

    #MyDir=os.path.dirname(__file__)

    # --- Parameters
    fstFile='../../../data/NREL5MW/Main_Onshore.fst' # Main OpenFAST file describing model
    fastExe='../../../data/NREL5MW/openfast_x64d_2022-02-21_hdclean.exe' # TODO adapt. Path to OpenFAST executable
    workDir='./_temp/' # Optional, Path where files will be written (for debugging). Set to None otherwsie.
    useLin =True  # Use OpenFAST linearizations to compute properties (recommended)
    useYAMS=False # Use YAMS python library  (in beta)

    if useLin:
        # --- Extract rigid body inertia of full structure
        print('------------------ FULL WIND TURBINE ------------------')
        out, _ = extractRigidBodyInertia(fstFile, fastExe=fastExe, workDir=workDir, cleanUp=False, zeroInit=True)
        print('Mass:                       ',                 out['mass'])
        print('Center of mass position (from platform ref):', out['CM'])
        print('Inertia at center of mass:\n',                 out['J_G'])
        print('Inertia at body origin (platform ref):\n',     out['J_Ref'])


        # --- Extract rigid body inertia of the rotor nacelle assembly
        print('------------------------- RNA -------------------------')
        out,_ = extractRNAInertia(fstFile, fastExe, workDir=workDir, cleanUp=False, includeYawBearingMass=True)
        print('RNA Mass:                       ',              out['mass'])
        print('RNA Center of mass position (from tower top):', out['CM'])
        print('RNA Inertia at center of mass:\n',              out['J_G'])
        print('RNA Inertia at tower top:\n',                   out['J_TT'])



    if useYAMS:
        # --- Alternative Method using python YAMS library (work in progress)
        from welib.yams.windturbine import FASTWindTurbine
        print('------------------ FULL WIND TURBINE ------------------')
        WT = FASTWindTurbine(fstFile, algo='OpenFAST')
        body = WT.WT_rigid
        # print(body)
        print('Mass:                   ',                         body.mass) 
        print('Center of mass position (from platform ref): ',    body.masscenter)
        print('Center of mass position (from global origin):',    body.masscenter_pos_global)
        print('Inertia at center of mass:\n',                     body.masscenter_inertia)
        print('Inertia at body origin (platform ref):\n',         body.inertia)
        print('Inertia at global origin: \n',                     body.inertia_at((0,0,0)))
        print('------------------------- RNA -------------------------')
        RNA = WT.RNA
        # print(RNA)
        print('Mass:                   '    , RNA.mass) 
        print('Center of mass position :'   , RNA.masscenter,           '(wrt tower top)' )          
        print('Inertia at center of mass:\n', RNA.masscenter_inertia)
        print('Inertia at at tower top\n',    RNA.inertia)
        

