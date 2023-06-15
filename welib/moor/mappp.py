"""
Wrapper for MAP DLL


Adapted from:
  https://bitbucket.org/mmasciola/map-plus-plus/src/master/

  Copyright (C) 2014 mdm                                      
  marco[dot]masciola[at]gmail                                 
"""                                                              
# Licensed to the Apache Software Foundation (ASF) under one   
# or more contributor license agreements.  See the NOTICE file 
# distributed with this work for additional information        
# regarding copyright ownership.  The ASF licenses this file   
# to you under the Apache License, Version 2.0 (the            
# "License"); you may not use this file except in compliance   
# with the License.  You may obtain a copy of the License at   
#                                                              
#   http://www.apache.org/licenses/LICENSE-2.0                 
#                                                              
# Unless required by applicable law or agreed to in writing,   
# software distributed under the License is distributed on an  
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       
# KIND, either express or implied.  See the License for the    
# specific language governing permissions and limitations            
# under the License.   



import numpy as np
import sys
from ctypes import *
import os


from distutils.sysconfig import get_config_var

from sys import platform


def setupLib(libpath):
    """
    Load the MAP++ library and setup the data types for the C-routines
    """
    print('Loading MAP library:', libpath)
    lib = cdll.LoadLibrary(libpath)

    class ModelData_Type(Structure):
        _fields_ = []


    # void * object ;
    # double gravity ;
    # double seaDensity ;
    # double depth ;
    # char fileName[255] ;
    # char summaryFileName[255] ;
    # char libraryInputLine[255] ;
    # char nodeInputLine[255] ;
    # char elementInputLine[255] ;
    # char optionInputLine[255] ;
    class InitializationData_Type(Structure):
        _fields_= [("object",c_void_p),
                   ("gravity",c_double),
                   ("seaDensity",c_double),
                   ("depth",c_double),
                   ("fileName",c_char*255),
                   ("summaryFileName",c_char*255),
                   ("libraryInputLine",c_char*255),
                   ("nodeInputLine",c_char*255),
                   ("elementInputLine",c_char*255),
                   ("optionInputLine",c_char*255)]

        

    # void * object ;
    # char progName[99] ;
    # char version[99] ;
    # char compilingData[24] ;
    # char * writeOutputHdr ;     int writeOutputHdr_Len ;
    # char * writeOutputUnt ;     int writeOutputUnt_Len ;
    class InitializationOutputData_Type(Structure):
        _fields_ = [("object",c_void_p),
                    ("progName",c_char*99),
                    ("version",c_char*99),
                    ("CompilingData",c_char*99),
                    ("writeOutputHdr",c_char_p),
                    ("writeOutputHdr_Len",c_int),
                    ("writeOutputUnt",c_char_p),
                    ("writeOutputUnt_Len",c_int)]

    class InputData_Type(Structure):
        _fields_ = []


    class OutputData_Type(Structure):
        _fields_ = [("object", c_void_p),
                    ("Fx",POINTER(c_double)),
                    ("Fx_Len", c_int),
                    ("Fy",POINTER(c_double)),
                    ("Fy_Len", c_int),
                    ("Fz",POINTER(c_double)),
                    ("Fz_Len", c_int),
                    ("WriteOuput",POINTER(c_float)),
                    ("WriteOutput_Len", c_int),
                    ("wrtOuput",POINTER(c_double)),
                    ("wrtOutput_Len", c_int)]

        
    # void * object ;
    # double g ;
    # double depth ;
    # double rhoSea ;

    class ParameterData_Type(Structure):
        _fields_ = [("object",c_void_p),
                    ("g",c_double),
                    ("depth",c_double), 
                    ("rhoSea", c_double)]

    class ConstraintData_Type(Structure):
        _fields_ = []

    class ContinuousData_Type(Structure):
        _fields_ = []

    # fields for the fortran types
    # 
    # MAP_EXTERNCALL MAP_InitInputType_t* map_create_init_type( char* msg, MAP_ERROR_CODE* status );
    # MAP_EXTERNCALL MAP_InitOutputType_t* map_create_initout_type( char* msg, MAP_ERROR_CODE* status );
    # MAP_EXTERNCALL MAP_InputType_t* map_create_input_type( char* msg, MAP_ERROR_CODE* status );
    # MAP_EXTERNCALL MAP_ParameterType_t* map_create_parameter_type( char* msg, MAP_ERROR_CODE* status );
    # MAP_EXTERNCALL MAP_ConstraintStateType_t* map_create_constraint_type( char* msg, MAP_ERROR_CODE* status );
    # MAP_EXTERNCALL MAP_OtherStateType_t* map_create_other_type( char* msg, MAP_ERROR_CODE* status );
    # MAP_EXTERNCALL MAP_OutputType_t* map_create_output_type( char* msg, MAP_ERROR_CODE* status );
    # MAP_EXTERNCALL MAP_ContinuousStateType_t* map_create_continuous_type( char* msg, MAP_ERROR_CODE* status );

    MapData_Type       = POINTER(ModelData_Type)
    MapInit_Type       = POINTER(InitializationData_Type)
    MapInitOut_Type    = POINTER(InitializationOutputData_Type)
    MapInput_Type      = POINTER(InputData_Type)
    MapOutput_Type     = POINTER(OutputData_Type)
    MapParameter_Type  = POINTER(ParameterData_Type)
    MapConstraint_Type = POINTER(ConstraintData_Type)
    MapContinuous_Type = POINTER(ContinuousData_Type)

    # read file stuff
    lib.set_init_to_null.argtype=[MapInit_Type, c_char_p, POINTER(c_int) ]
    lib.map_set_summary_file_name.argtype=[MapInit_Type, c_char_p, POINTER(c_int) ]
    lib.map_add_cable_library_input_text.argtype=[MapInit_Type]
    lib.map_add_node_input_text.argtype=[MapInit_Type]
    lib.map_add_line_input_text.argtype=[MapInit_Type]
    lib.map_add_options_input_text.argtype=[MapInit_Type]

    lib.map_create_init_type.argtype       = [ c_char_p, POINTER(c_int) ]
    lib.map_create_initout_type.argtype    = [ c_char_p, POINTER(c_int) ]
    lib.map_create_input_type.argtype      = [ c_char_p, POINTER(c_int) ]
    lib.map_create_parameter_type.argtype  = [ c_char_p, POINTER(c_int) ]
    lib.map_create_constraint_type.argtype = [ c_char_p, POINTER(c_int) ]
    lib.map_create_other_type.argtype      = [ c_char_p, POINTER(c_int) ]
    lib.map_create_output_type.argtype     = [ c_char_p, POINTER(c_int) ]
    lib.map_create_continuous_type.argtype = [ c_char_p, POINTER(c_int) ]
    lib.map_create_continuous_type.argtype = [ MapData_Type ]
    
    lib.map_create_init_type.restype       = MapInit_Type
    lib.map_create_initout_type.restype    = MapInitOut_Type
    lib.map_create_input_type.restype      = MapInput_Type
    lib.map_create_parameter_type.restype  = MapParameter_Type
    lib.map_create_constraint_type.restype = MapConstraint_Type
    lib.map_create_other_type.restype      = MapData_Type
    lib.map_create_output_type.restype     = MapOutput_Type
    lib.map_create_continuous_type.restype = MapContinuous_Type

    lib.map_set_sea_depth.argtypes   = [ MapParameter_Type, c_double ]
    lib.map_set_gravity.argtypes     = [ MapParameter_Type, c_double ]
    lib.map_set_sea_density.argtypes = [ MapParameter_Type, c_double ]
    
    lib.map_size_lines.restype = c_int

    # numeric routines
    lib.map_residual_function_length.restype = c_double
    lib.map_residual_function_height.restype = c_double
    lib.map_jacobian_dxdh.restype            = c_double
    lib.map_jacobian_dxdv.restype            = c_double
    lib.map_jacobian_dzdh.restype            = c_double
    lib.map_jacobian_dzdv.restype            = c_double
 
    lib.map_residual_function_length.argtypes = [ MapData_Type, c_int, c_char_p, POINTER(c_int) ]
    lib.map_residual_function_height.argtypes = [ MapData_Type, c_int, c_char_p, POINTER(c_int) ]
    lib.map_jacobian_dxdh.argtypes            = [ MapData_Type, c_int, c_char_p, POINTER(c_int) ]
    lib.map_jacobian_dxdv.argtypes            = [ MapData_Type, c_int, c_char_p, POINTER(c_int) ]
    lib.map_jacobian_dzdh.argtypes            = [ MapData_Type, c_int, c_char_p, POINTER(c_int) ]
    lib.map_jacobian_dzdv.argtypes            = [ MapData_Type, c_int, c_char_p, POINTER(c_int) ]
      
    lib.map_get_fairlead_force_2d.argtypes = [POINTER(c_double), POINTER(c_double), MapData_Type, c_int, c_char_p, POINTER(c_int)]
    lib.map_get_fairlead_force_3d.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), MapData_Type, c_int, c_char_p, POINTER(c_int)]      
    try:
        lib.map_get_anchor_force_2d.argtypes = [POINTER(c_double), POINTER(c_double), MapData_Type, c_int, c_char_p, POINTER(c_int)]
        lib.map_get_anchor_force_3d.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), MapData_Type, c_int, c_char_p, POINTER(c_int)]
    except:
        print('[WARN] map_get_anchor_force_2d not available in this version of MAP')

    # plot routines
    lib.map_plot_x_array.argtypes = [ MapData_Type, c_int, c_int, c_char_p, POINTER(c_int) ]
    lib.map_plot_x_array.restype  = POINTER(c_double)
    lib.map_plot_y_array.argtypes = [ MapData_Type, c_int, c_int, c_char_p, POINTER(c_int) ]
    lib.map_plot_y_array.restype  = POINTER(c_double)
    lib.map_plot_z_array.argtypes = [ MapData_Type, c_int, c_int, c_char_p, POINTER(c_int) ]
    lib.map_plot_z_array.restype  = POINTER(c_double)
    lib.map_plot_array_free.argtypes = [ POINTER(c_double) ]
 
    # modifyers
    lib.map_offset_vessel.argtypes = [MapData_Type, MapInput_Type, c_double, c_double, c_double, c_double, c_double, c_double, c_char_p, POINTER(c_int)]        
    lib.map_linearize_matrix.argtypes = [MapInput_Type, MapParameter_Type, MapData_Type, MapOutput_Type, MapConstraint_Type, c_double, POINTER(c_int), c_char_p]        
    lib.map_linearize_matrix.restype  = POINTER(POINTER(c_double))
    lib.map_free_linearize_matrix.argtypes = [POINTER(POINTER(c_double))]

    try:
        lib.map_f_op.argtypes = [MapInput_Type, MapParameter_Type, MapData_Type, MapOutput_Type, MapConstraint_Type, POINTER(c_int), c_char_p]        
        lib.map_f_op.restype  = POINTER(c_double)
        lib.map_free_f_op.argtypes = [POINTER(c_double)]
    except:
        print('[WARN] map_f_op not available in this version of MAP')

    lib.map_init.argtypes = [ MapInit_Type,
                              MapInput_Type,
                              MapParameter_Type,
                              MapContinuous_Type,
                              c_void_p,
                              MapConstraint_Type,
                              MapData_Type,
                              MapOutput_Type,
                              MapInitOut_Type,
                              POINTER(c_int),
                              c_char_p]


    lib.map_update_states.argtypes = [ c_float,
                                       c_int,
                                       MapInput_Type,
                                       MapParameter_Type,
                                       MapContinuous_Type,
                                       c_void_p,
                                       MapConstraint_Type,
                                       MapData_Type,
                                       POINTER(c_int),
                                       c_char_p]

    lib.map_calc_output.argtypes = [c_float,
                                    MapInput_Type,
                                    MapParameter_Type,
                                    MapContinuous_Type,
                                    c_void_p,
                                    MapConstraint_Type,
                                    MapData_Type,
                                    MapOutput_Type,
                                    POINTER(c_int),
                                    c_char_p]            
    
    lib.map_end.argtypes = [ MapInput_Type,
                             MapParameter_Type,
                             MapContinuous_Type,
                             c_void_p,
                             MapConstraint_Type,
                             MapData_Type,
                             MapOutput_Type,
                             POINTER(c_int),
                             c_char_p]
    
    lib.map_initialize_msqs_base.argtypes = [MapInput_Type,
                                             MapParameter_Type,
                                             MapContinuous_Type,
                                             MapConstraint_Type,
                                             MapData_Type,
                                             MapOutput_Type,
                                             MapInitOut_Type]

    lib.map_size_lines.argtypes = [ MapData_Type,
                                    POINTER(c_int),
                                    c_char_p]

    lib.map_get_header_string.argtypes = [c_void_p, POINTER(c_char_p),   MapData_Type]
    lib.map_get_unit_string.argtypes = [c_void_p, POINTER(c_char_p),   MapData_Type]
    try:
        lib.map_offset_fairlead.argtypes = [MapInput_Type, c_int, c_double, c_double, c_double, c_char_p, POINTER(c_int)]                    
    except:
        print('[WARN] map_offset_fairlead not available in this version of MAP')
    return lib

class pyMAP(object):

    lib = None # Will hold the dll library handle

    # read file stuff
    @classmethod
    def initLib(cls, dllFileName=None):
        if dllFileName is None:
            MyDir=os.path.dirname(__file__)
            if platform == "win32":
                dllFileName=os.path.join(MyDir,'libmap.dll')
            else:
                dllFileName=os.path.join(MyDir,'libmap.so')
#             libext = get_config_var('EXT_SUFFIX')
#             if libext is None or libext == '':
#                 if platform == "linux" or platform == "linux2":
#                     libext = '.so'
#                 elif platform == "darwin":
#                     #libext = '.dyld'
#                     libext = '.so'
#                 elif platform == "win32":
#                     #libext = '.dll'
#                     libext = '.pyd'
#                 elif platform == "cygwin":
#                     libext = '.dll'
# 
#             maplib = '_libmap' + libext
#             dllFileName = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + os.path.sep + maplib

        # Load the library and setup C-interfaces
        lib = setupLib(dllFileName)

        # Store in object
        cls.lib = lib

        cls.ierr = c_int(0)
        cls.status = create_string_buffer(1024)
        cls.val = c_double


    def __init__( self, filename=None, lines=None, WtrDepth=None, gravity=None, WtrDens=None, sumFile=None, dllFileName=None ) :
        """

        """
        # Call Class method
        if pyMAP.lib is None:
            pyMAP.initLib(dllFileName)

        # Misc Data
        self._gravity=None
        self._WtrDepth=None
        self._WtrDens=None
        self._filename=None
        self._K_lin       = None # Linear stiffness matrix
        self._K_lin_point = None # Point where linear stiffness matrix was computed
        self.Nodes=[]
        self._initialized = False

        # Wrapper data
        self.f_type_d       = self.CreateDataState()
        self.f_type_u       = self.CreateInputState( )
        self.f_type_x       = self.CreateContinuousState( )
        self.f_type_p       = self.CreateParameterState( )
        self.f_type_y       = self.CreateOutputState( )
        self.f_type_z       = self.CreateConstraintState( )
        self.f_type_init    = self.CreateInitState( )
        self.f_type_initout = self.CreateInitoutState( )
        pyMAP.lib.set_init_to_null(self.f_type_init, self.status, pointer(self.ierr) )
        pyMAP.lib.map_initialize_msqs_base(self.f_type_u, self.f_type_p, self.f_type_x, self.f_type_z, self.f_type_d, self.f_type_y, self.f_type_initout)


        callInit = False
        # Read input file (either OpenFAST or MAP)
        if filename is not None:
            ext = os.path.splitext(filename)[1].lower()
            if ext=='.fst':
                from welib.weio import FASTInputFile
                fst = FASTInputFile(filename)
                if WtrDepth is None:
                    WtrDepth=fst['WtrDpth'] # m
                if gravity is None:
                    gravity=fst['gravity'] # m/s^2
                if WtrDens is None:
                    WtrDens=fst['WtrDens'] # kg/m^3
                filename = os.path.join(os.path.dirname(filename), fst['MooringFile'].replace('"',''))
            else:
                # Assume that it's a MAP input file
                pass

            self.read_file(filename)
            callInit = True
            sumFile = os.path.splitext(filename)[0]+'.map.sum'

        # Read input lines
        if lines is not None:
            self.read_list_input(lines)
            callInit = True

        # Set summary file
        if sumFile is None:
            sumFile = 'outlist.map.sum'
        self.summary_file(sumFile)

        # Set env conditions if provided or read from file
        if WtrDepth is not None:
            self.map_set_sea_depth(WtrDepth)  # m
        if gravity is not None:
            self.map_set_gravity(gravity)     # m/s^2
        if WtrDens is not None:
            self.map_set_sea_density(WtrDens) # kg/m^3

        # If all inputs have been provided, initialize
        if self._WtrDens is not None and self._WtrDepth is not None and self._gravity is not None:
            if callInit:
                self.init()

    def __repr__(self):
        s ='<{} object> with attributes:\n'.format(type(self).__name__)
        s+='|- Nodes: {}\n'.format(len(self.Nodes))
        for n in self.Nodes:
            s+='|   {}\n'.format(n)
        s+='|- _WtrDepth: {}\n'.format(self._WtrDepth)
        s+='|- _WtrDens : {}\n'.format(self._WtrDens)
        s+='|- _gravity : {}\n'.format(self._gravity)
        return s

    def init( self ):
        if self._initialized:
            raise Exception('Cannot call map_init twice')
        pyMAP.lib.map_init( self.f_type_init, self.f_type_u, self.f_type_p, self.f_type_x, None, self.f_type_z, self.f_type_d, self.f_type_y, self.f_type_initout, pointer(self.ierr), self.status )
        if self.ierr.value != 0 : print(self.status.value)
        self._initialized = True


    def size_lines(self):
        size = pyMAP.lib.map_size_lines(self.f_type_d, pointer(self.ierr), self.status )
        if self.ierr.value != 0 : print(self.status.value)
        return size


    def update_states(self, t, dt):
        pyMAP.lib.map_update_states(c_float(t), c_int(dt), self.f_type_u, self.f_type_p, self.f_type_x, None, self.f_type_z, self.f_type_d, pointer(self.ierr), self.status )
        if self.ierr.value != 0 : print(self.status.value)


    def calc_output(self, t):
        pyMAP.lib.map_calc_output(c_float(t), self.f_type_u, self.f_type_p, self.f_type_x, None, self.f_type_z, self.f_type_d, self.f_type_y, pointer(self.ierr), self.status )
        if self.ierr.value != 0 :
            print(self.status.value)


    def get_output(self):
        size_x, size_y, size_z = self.f_type_y.contents.Fx_Len, self.f_type_y.contents.Fy_Len, self.f_type_y.contents.Fz_Len
        arr_x, arr_y, arr_z = [None]*size_x, [None]*size_y, [None]*size_z
        fx = [self.f_type_y.contents.Fx[j] for j in range(size_x)]
        fy = [self.f_type_y.contents.Fy[j] for j in range(size_y)]
        fz = [self.f_type_y.contents.Fz[j] for j in range(size_z)]        
        return fx, fy, fz


    def get_output_labels(self):
        size = self.f_type_y.contents.WriteOutput_Len + self.f_type_y.contents.wrtOutput_Len
        string_buffers = [create_string_buffer(16) for i in range(size)]
        pointers = (c_char_p*size)(*map(addressof, string_buffers))
        pyMAP.lib.map_get_header_string(None, pointers, self.f_type_d)
        return [s.value for s in string_buffers]

    
    def get_output_units(self):
        size = self.f_type_y.contents.WriteOutput_Len + self.f_type_y.contents.wrtOutput_Len
        string_buffers = [create_string_buffer(16) for i in range(size)]
        pointers = (c_char_p*size)(*map(addressof, string_buffers))
        pyMAP.lib.map_get_unit_string(None, pointers, self.f_type_d)
        return [s.value for s in string_buffers]
    
    
    def get_output_buffer(self):
        size_float, size_double = self.f_type_y.contents.WriteOutput_Len, self.f_type_y.contents.wrtOutput_Len
        arr_float, arr_double = [None]*size_float, [None]*size_double
        arr_float = [self.f_type_y.contents.WriteOutput[j] for j in range(size_float)]
        arr_double = [self.f_type_y.contents.wrtOuput[j] for j in range(size_double)]
        return arr_float + arr_double
    
        
    # Calls function in main.c and fordatamanager.c to delete insteads of c structs. First, the malloc'ed arrays need to vanish
    # gracefully; we accomplish this by calling MAP_End(...) routine. Then, the structs themself are deleted. Order is important.
    # 
    # MAP_EXTERNCALL int MAP_End ( InputData *u, ParameterData *p, ContinuousData *x, ConstraintData *z, ModelData *data, OutputData *y, char *map_msg, MAP_ERROR_CODE *ierr )
    # MAP_EXTERNCALL void MAP_Input_Delete( InputData* u )
    # MAP_EXTERNCALL void MAP_Param_Delete( ParameterData* p )
    # MAP_EXTERNCALL void MAP_ContState_Delete( InputData* x )
    # MAP_EXTERNCALL void MAP_ConstrState_Delete( InputData* z )
    # MAP_EXTERNCALL void MAP_Output_Delete( InputData* y )
    # MAP_EXTERNCALL void MAP_OtherState_Delete( ModelData* data )
    def end(self):
        pyMAP.lib.map_end(self.f_type_u, self.f_type_p, self.f_type_x, None, self.f_type_z, self.f_type_d, self.f_type_y, pointer(self.ierr), self.status)


    # Set a name for the MAP summary file. Does not need to be called. If not called, the default name is 'outlist.map.sum'
    def summary_file(self, echo_file):
        self.f_type_init.contents.summaryFileName = echo_file.encode('utf-8')
        pyMAP.lib.map_set_summary_file_name(self.f_type_init, self.status, pointer(self.ierr) )


    # Calls function in fortdatamanager.c to create instance of c structs
    # MAP_EXTERNCALL InitializationData* MAP_InitInput_Create( char* map_msg, MAP_ERROR_CODE* ierr )
    def CreateInitState( self ) :
        obj = pyMAP.lib.map_create_init_type( self.status, pointer(self.ierr) )
        if self.ierr.value != 0 : print(self.status.value)
        return obj

    # Calls function in fortdatamanager.c to create instance of c structs
    # MAP_EXTERNCALL void MAP_InitOutput_Delete( InputData* io )
    def CreateInitoutState( self ) :
        obj = pyMAP.lib.map_create_initout_type( self.status, pointer(self.ierr) )
        if self.ierr.value != 0 : print(self.status.value)
        return obj


    # Calls function in fortdatamanager.c to create instance of c structs
    # MAP_EXTERNCALL ModelData *MAP_OtherState_Create( char *map_msg, MAP_ERROR_CODE *ierr )
    def CreateDataState( self ) :
        obj = pyMAP.lib.map_create_other_type( self.status, pointer(self.ierr) )
        if self.ierr.value != 0 : print(self.status.value)
        return obj


    # Calls function in fortdatamanager.c to create instance of c structs
    # MAP_EXTERNCALL InputData* MAP_Input_Create( char* map_msg, MAP_ERROR_CODE *ierr )
    def CreateInputState( self ) :
        obj = pyMAP.lib.map_create_input_type( self.status, pointer(self.ierr) )
        if self.ierr.value != 0 : print(self.status.value)
        return obj


    # Calls function in fortdatamanager.c to create instance of c structs
    # MAP_EXTERNCALL ContinuousData* MAP_ContState_Create( char* map_msg, MAP_ERROR_CODE *ierr )
    def CreateContinuousState( self ) :
        obj = pyMAP.lib.map_create_continuous_type( self.status, pointer(self.ierr) )
        if self.ierr.value != 0 : print(self.status.value)
        return obj

    # Calls function in fortdatamanager.c to create instance of c structs
    # MAP_EXTERNCALL OutputData *MAP_Output_Create( char *map_msg, MAP_ERROR_CODE *ierr )
    def CreateOutputState( self ) :
        obj = pyMAP.lib.map_create_output_type( self.status, pointer(self.ierr) )
        if self.ierr.value != 0 : print(self.status.value)
        return obj


    # Calls function in fortdatamanager.c to create instance of c structs
    # MAP_EXTERNCALL ConstraintData* MAP_ConstrState_Create( char* map_msg, MAP_ERROR_CODE *ierr )
    def CreateConstraintState( self ) :
        obj = pyMAP.lib.map_create_constraint_type( self.status, pointer(self.ierr) )
        if self.ierr.value != 0 : print(self.status.value)
        return obj


    # Calls function in fortdatamanager.c to create instance of c structs
    # MAP_EXTERNCALL ParameterData* MAP_Param_Create( char* map_msg, MAP_ERROR_CODE *ierr )
    def CreateParameterState( self ) :
        obj = pyMAP.lib.map_create_parameter_type( self.status, pointer(self.ierr) )
        if self.ierr.value != 0 : print(self.status.value)
        return obj


    def map_set_sea_depth( self, depth ):
         self._WtrDepth=depth
         pyMAP.lib.map_set_sea_depth( self.f_type_p, depth )

         
    def map_set_gravity( self, g ):
        self._gravity=g
        pyMAP.lib.map_set_gravity( self.f_type_p, g )

        
    def map_set_sea_density( self, rho ):
        self._WtrDens=rho
        pyMAP.lib.map_set_sea_density( self.f_type_p, rho )

    def plot_x( self, lineNum, length ) :
        arr = [None]*length
        array = POINTER(c_double)
        array = pyMAP.lib.map_plot_x_array( self.f_type_d, lineNum, length, self.status, pointer(self.ierr) )        
        if self.ierr.value != 0 :
            print(self.status.value)
            self.end( )
            pyMAP.lib.map_plot_array_free( array )        
            sys.exit('MAP terminated premature.')
        arr = [array[j] for j in range(length)]        
        pyMAP.lib.map_plot_array_free( array )        
        return arr 

    
    def plot_y( self, lineNum, length ) :
        arr = [None]*length
        array = POINTER(c_double)
        array = pyMAP.lib.map_plot_y_array( self.f_type_d, lineNum, length, self.status, pointer(self.ierr) )        
        if self.ierr.value != 0 :
            print(self.status.value)
            self.end( )
            pyMAP.lib.map_plot_array_free( array )        
            sys.exit('MAP terminated premature.')
        arr = [array[j] for j in range(length)]        
        pyMAP.lib.map_plot_array_free( array )        
        return arr 


    def plot_z( self, lineNum, length ) :
        arr = [None]*length
        array = POINTER(c_double)
        array = pyMAP.lib.map_plot_z_array( self.f_type_d, lineNum, length, self.status, pointer(self.ierr) )        
        if self.ierr.value != 0 :
            print(self.status.value)
            self.end( )
            pyMAP.lib.map_plot_array_free( array )        
            sys.exit('MAP terminated premature.')
        arr = [array[j] for j in range(length)]        
        pyMAP.lib.map_plot_array_free( array )        
        return arr
    

    def get_fairlead_force_2d(self, index):
        """Gets the horizontal and vertical fairlead force in a 2D plane along the 
        straight-line line. Must ensure update_states() is called before accessing 
        this function. The function will not solve the forces for a new vessel position
        if it updated. , otherwise the fairlead forces are not updated with the new 
        vessel position. Called C function:
        
        MAP_EXTERNCALL void map_get_fairlead_force_2d(double* H, double* V, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr);
    
        :param index: The line number the fairlead forces are being requested for. Zero indexed
        :returns: horizontal and vertical fairlead force [N]
    
        >>> H,V = print get_fairlead_force_2d(1)        
        """
        H_ref = c_double(-999.9)
        V_ref = c_double(-999.9)
        pyMAP.lib.map_get_fairlead_force_2d( pointer(H_ref), pointer(V_ref),self.f_type_d, index, self.status, pointer(self.ierr))
        return H_ref.value, V_ref.value
    
    
    def get_fairlead_force_3d(self, index):
        """Gets the horizontal and vertical fairlead force in a 3D frame along relative 
        referene global axis. Must ensure update_states() is called before accessing 
        this function. The function will not solve the forces for a new vessel position
        if it updated. , otherwise the fairlead forces are not updated with the new 
        vessel position. Called C function:
        
        MAP_EXTERNCALL void map_get_fairlead_force_3d(double* fx, double* fy, double* fz, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr);
    
        :param index: The line number the fairlead forces are being requested for. Zero indexed
        :returns: horizontal and vertical fairlead force [N]
    
        >>> fx,fy,fz = get_fairlead_force_3d(1)        
        """
        fx = c_double(-999.9)
        fy = c_double(-999.9)
        fz = c_double(-999.9)
        pyMAP.lib.map_get_fairlead_force_3d( pointer(fx), pointer(fy), pointer(fz), self.f_type_d, index, self.status, pointer(self.ierr))
        return fx.value, fy.value, fz.value
        

    def get_anchor_force_2d(self, index):
        """Gets the horizontal and vertical anchor force in a 2D plane along the 
        straight-line line. Must ensure update_states() is called before accessing 
        this function. The function will not solve the forces for a new vessel position
        if it updated. , otherwise the anchor forces are not updated with the new 
        vessel position. Called C function:
        
        MAP_EXTERNCALL void map_get_anchor_force_2d(double* H, double* V, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr);
    
        :param index: The line number the anchor forces are being requested for. Zero indexed
        :returns: horizontal and vertical anchor force [N]
    
        >>> H,V = print get_anchor_force_2d(1)        
        """
        Ha_ref = c_double(-999.9)
        Va_ref = c_double(-999.9)
        pyMAP.lib.map_get_anchor_force_2d( pointer(Ha_ref), pointer(Va_ref),self.f_type_d, index, self.status, pointer(self.ierr))
        return Ha_ref.value, Va_ref.value
    
    
    def get_anchor_force_3d(self, index):
        """Gets the horizontal and vertical anchor force in a 3D frame along relative 
        referene global axis. Must ensure update_states() is called before accessing 
        this function. The function will not solve the forces for a new vessel position
        if it updated. , otherwise the anchor forces are not updated with the new 
        vessel position. Called C function:
        
        MAP_EXTERNCALL void map_get_anchor_force_3d(double* fx, double* fy, double* fz, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr);
    
        :param index: The line number the anchor forces are being requested for. Zero indexed
        :returns: horizontal and vertical anchor force [N]
    
        >>> fx,fy,fz = get_anchor_force_3d(1)        
        """
        fxa = c_double(-999.9)
        fya = c_double(-999.9)
        fza = c_double(-999.9)
        pyMAP.lib.map_get_anchor_force_3d( pointer(fxa), pointer(fya), pointer(fza), self.f_type_d, index, self.status, pointer(self.ierr))
        return fxa.value, fya.value, fza.value
    
    def funcl( self, i ) :
        self.val = pyMAP.lib.map_residual_function_length(self.f_type_d, i, self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print(self.status.value)
            self.end( )
            sys.exit('MAP terminated premature.')
        return self.val


    def funch( self, i ) :
        self.val = pyMAP.lib.map_residual_function_height(self.f_type_d, i, self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print(self.status.value)
            self.end( )
            sys.exit('MAP terminated premature.')
        return self.val


    def dxdh( self, i ) :
        self.val = pyMAP.lib.map_jacobian_dxdh( self.f_type_d, i, self.status, pointer(self.ierr) )
        if self.ierr.value != 0 :
            print(self.status.value)
            self.end( )
            sys.exit('MAP terminated premature.')
        return self.val


    def dxdv( self, i ) :
        self.val = pyMAP.lib.map_jacobian_dxdv( self.f_type_d, i, self.status, pointer(self.ierr) )
        if self.ierr.value != 0 :
            print(self.status.value)
            self.end( )
            sys.exit('MAP terminated premature.')
        return self.val


    def dzdh( self, i ) :
        self.val = pyMAP.lib.map_jacobian_dzdh( self.f_type_d, i, self.status, pointer(self.ierr) )
        if self.ierr.value != 0 :
            print(self.status.value)
            self.end( )
            sys.exit('MAP terminated premature.')
        return self.val


    def dzdv( self, i ) :
        self.val = pyMAP.lib.map_jacobian_dzdv( self.f_type_d, i, self.status, pointer(self.ierr) )
        if self.ierr.value != 0 :
            print(self.status.value)
            self.end( )
            sys.exit('MAP terminated premature.')
        return self.val


    def linear( self, epsilon=1.e-3) :
        """
        Return linear matrix, transpose of stiffness matrix at (0,0,0)
        """
        array = POINTER(POINTER(c_double))
        array = pyMAP.lib.map_linearize_matrix( self.f_type_u, self.f_type_p, self.f_type_d, self.f_type_y, self.f_type_z, epsilon, pointer(self.ierr), self.status)        
        if self.ierr.value != 0 :
           print(self.status.value)
           self.end( )
           sys.exit('MAP terminated premature.')
        arr = [[array[j][i] for i in range(6)] for j in range(6)] 
        pyMAP.lib.map_free_linearize_matrix(array)        
        # Stiffness matrix is defined at (0,0,0)
        K_0   = np.array(arr) # NOTE: this is a transposed of a stiffness matrix !!!
        return K_0


    def stiffness_matrix(self, epsilon=1.e-3, point=None) :
        # Stiffness matrix at (0,0,0)
        K_0 = self.linear(epsilon).T # NOTE: transpose needed

        # Operating point loads
        F_op_0 = self.f_op()

        if point is None:
            K_D    = K_0
            F_op_D = F_op_0
        else:
            # Transfer to requested reference point
            r_S  = np.array((0,0,0)) # Source # TODO is it OK when vessel has been displaced?
            r_D = np.array(point)    # Destination
            r0  = r_S - r_D  # 
            K_D = - translateLoadsJacobian(-K_0, r0, F_op_0[:3]) # Jacobians are -K
            F_op_D = F_op_0.copy()
            F_op_D[3:6] += np.cross(r0,F_op_0[:3])
        self._K_lin       = K_D
        self._K_lin_point = point
        return K_D, F_op_D

    def f_op(self) :
        try:
            array = POINTER(c_double)
            array = pyMAP.lib.map_f_op( self.f_type_u, self.f_type_p, self.f_type_d, self.f_type_y, self.f_type_z, pointer(self.ierr), self.status)        
            if self.ierr.value != 0 :
               print(self.status.value)
               self.end( )
               sys.exit('MAP terminated premature.')
            arr = [array[i] for i in range(6)]
            pyMAP.lib.map_free_f_op(array)        
            Fop   = np.array(arr)
        except:
            print('[WARN] MAP f_op not available in this version of the library')
            Fop=np.zeros(6)
        return Fop
    

    def displace_vessel(self,x,y,z,phi,the,psi) :
        """
        NOTE:  
         - phi is rotation around x  
         - the is rotation around y  
         - psi is rotation around z  
        """
        pyMAP.lib.map_offset_vessel(self.f_type_d, self.f_type_u, x,y,z,phi,the,psi, self.status, pointer(self.ierr) )
        if self.ierr.value != 0 :
            print(self.status.value)
            self.end( )
            sys.exit('MAP terminated premature.')    

            
    #def offset_fairlead(self,x,y,z, idx):
    #    pyMAP.lib.map_offset_fairlead(self.f_type_u, c_int(idx), x, y, z, self.status, pointer(self.ierr) )
    #    if self.ierr.value != 0 :
    #        print(self.status.value)
    #        self.end( )
    #        sys.exit('MAP terminated premature.')   
    def read_file( self, fileName ):
        with open(fileName, 'r', errors="surrogateescape") as f:
            lines = f.read().splitlines()
        self.read_list_input( lines )
# 
                    
    def read_list_input(self, listIn):
        assert isinstance(listIn, list), 'Must input a python list of strings'
        assert len(listIn) >= 8, 'Must have at least 4 sections, 3 lines per section'
        assert type(listIn[0]) == type(''), 'List elements must be strings'
        charptr = POINTER(c_char)
        option_breaks = ("LINE DICTIONARY", "NODE PROPERTIES", "LINE PROPERTIES", "SOLVER OPTIONS")

        # --- Read file lines and store in different lists
        sCabLib = []
        sNodes = []
        sProps = []
        sOpts  = []
        listIter = iter(listIn)
        for line in listIter:
            line = line.strip()
            if "LINE DICTIONARY" in line.upper():
                for _ in range(3): line = next(listIter) # Process Header
                while not any(opt in line for opt in option_breaks):
                    sCabLib.append(line)
                    line = next(listIter)

            if "NODE PROPERTIES" in line.upper():
                for _ in range(3): line = next(listIter) # Process Header
                while not any(opt in line for opt in option_breaks):
                    sNodes.append(line)
                    line = next(listIter)

            if "LINE PROPERTIES" in line.upper():
                for _ in range(3): line = next(listIter) # Process Header
                while not any(opt in line for opt in option_breaks):
                    sProps.append(line)
                    line = next(listIter)

            if "SOLVER OPTIONS" in line.upper():
                for _ in range(2): line = next(listIter) # Process Header
                try:
                    line=next(listIter)
                    while not any(opt in line for opt in option_breaks):
                        sOpts.append(line)
                        line = next(listIter,"SOLVER OPTIONS")
                except StopIteration:
                    pass

        # --- Setup a list of Nodes. 
        self.Nodes=[]
        for line in sNodes:
            sp=line.split()
            n = {}
            n['ID']   = int(sp[0])
            n['type'] = sp[1]
            z=sp[4]
            if z.strip()=='depth':
                z=self._WtrDepth
            elif z.startswith('#'):
                print('>>> TODO mappp.py: z starting with #: ', z)
                z=float(z[1:])
            else:
                z=float(z)
            x=sp[2]
            y=sp[3]
            if x.startswith('#'):
                print('>>> TODO mappp.py: x starting with #: ', z)
                x=x[1:]
            if y.startswith('#'):
                print('>>> TODO mappp.py: y starting with #: ', z)
                y=y[1:]
            x = float(x)
            y = float(y)
            n['position'] = np.array((x, y, z))
            self.Nodes.append(n)

        # --- Load input lines into library
        for line in sCabLib:
            self.f_type_init.contents.libraryInputLine =  (line+'\n\0').encode('utf-8')
            pyMAP.lib.map_add_cable_library_input_text(self.f_type_init)                    
        for line in sNodes:
            self.f_type_init.contents.nodeInputLine = (line+'\n\0').encode('utf-8')
            pyMAP.lib.map_add_node_input_text(self.f_type_init)
        for line in sProps:
            self.f_type_init.contents.elementInputLine =(line+'\n\0').encode('utf-8')
            pyMAP.lib.map_add_line_input_text(self.f_type_init)
        for line in sOpts:
            self.f_type_init.contents.optionInputLine = (line+'\n\0').encode('utf-8')
            pyMAP.lib.map_add_options_input_text(self.f_type_init)            

    # --------------------------------------------------------------------------------}
    # --- Utils 
    # --------------------------------------------------------------------------------{
    # user function to plot the mooring profile and footprint
    def plot(self, numPoints, fig=None, ax=None, colors=None, ls='-'):
        """ plot the mooring profile """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        if fig is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
        if colors is None:
            colors = ['b']
        for i in range(self.size_lines()):
            x = self.plot_x( i, numPoints ) # i is the the line number, and 20 is the number of points plotted on the line 
            y = self.plot_y( i, numPoints)
            z = self.plot_z( i, numPoints)        
            ax.plot(x,y,z, ls=ls, color=colors[np.mod(i,len(colors))])     
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_zlabel('Z [m]')        
        return fig, ax

class Vessel:
    time = []
    x = []
    y = []
    z = []
    phi = []
    the = []
    psi = []
    out = []


def get_vessel_column_index(name,out_chanel):
    index = [0,0,0,0,0,0,0]
    fp = open(name) 
    for i,line in enumerate(fp):
        if i==6:
            words = line.split()
            for j in range(0,len(words)):
                if words[j]=='PtfmTDxi' or words[j]=='PtfmSurge':
                    index[0] = j;
                if words[j]=='PtfmTDyi' or words[j]=='PtfmSway':
                    index[1] = j;
                if words[j]=='PtfmTDzi' or words[j]=='PtfmHeave':
                    index[2] = j;
                if words[j]=='PtfmRDxi' or words[j]=='PtfmRoll':
                    index[3] = j;
                if words[j]=='PtfmRDyi' or words[j]=='PtfmPitch':
                    index[4] = j;
                if words[j]=='PtfmRDzi' or words[j]=='PtfmYaw':
                    index[5] = j;
                if words[j]==out_chanel:
                    index[6] = j;
    fp.close()
    return index


def set_vessel_prescribed_motion(table,index):
    vessel = Vessel()
    N = 1000#len(table)
    vessel.time = [float(table[i][0]) for i in range(8,N)]
    vessel.x = [float(table[i][index[0]]) for i in range(8,N)]
    vessel.y = [float(table[i][index[1]]) for i in range(8,N)]
    vessel.z = [float(table[i][index[2]]) for i in range(8,N)]
    vessel.phi = [float(table[i][index[3]]) for i in range(8,N)]
    vessel.the = [float(table[i][index[4]]) for i in range(8,N)]
    vessel.psi = [float(table[i][index[5]]) for i in range(8,N)]    
    vessel.out = [float(table[i][index[6]]) for i in range(8,N)]    
    return vessel



def translateLoadsJacobian(JS, r0, FS0):
    """ 
    Transfer Jacobians of loads at a source point "S" to a destination point "D"
    assuming a rigid body motion between the source and destination point

    INPUTS:
    - JS: Jacobian of loads at the source point, 6x6 array dF_S/dx_S
          where F_S  = (f_x, f_y, f_z, m_x, m_y, m_z)_S
          where dx_S = (dx, dy, dz, dt_x, dt_y, dt_z)_S (displacement and small rotations)
    - r0 = rS-rD : 3-vector from the destination point to the source point
    - FS0: 3-vector of the operating point forces at the source node

    OUTPUTS:
    - JD: Jacobian of loads at the source point, 6x6 array dF_S/dx_S

    Reference:
       See Branlard, Jonkman, Brown to be published 2022

    """
    def skew(x):
        """ Returns the skew symmetric matrix M, such that: cross(x,v) = M v 
        [ 0, -z , y]
        [ z,  0 ,-x]
        [-y,  x , 0]
        """
        x=np.asarray(x).ravel()
        return np.array([[0, -x[2], x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]])
    r0til  = skew(r0)
    FS0til = skew(FS0)
    I3 = np.eye(3)
    Z3 = np.zeros((3,3))
    T1 = np.block([ [ I3   ,  Z3 ],
                    [ r0til,  I3 ]])
    T2 = np.block([ [ I3   ,-r0til ],
                    [ Z3   ,  I3 ]])
    T3 = np.block([ [ Z3   ,  Z3 ],
                    [ Z3   ,  FS0til.dot(r0til) ] ])
    JD = T1.dot(JS.dot(T2)) + T3
    return JD


Map=pyMAP # backward compatibility
