from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from io import open
from builtins import range
from builtins import str
from future import standard_library
standard_library.install_aliases()
import os
import numpy as np
import re
import pandas as pd

from .fast_input_file import FASTInputFile

__all__  = ['FASTInputDeck']
# --------------------------------------------------------------------------------}
# --- Full FAST input deck
# --------------------------------------------------------------------------------{
class FASTInputDeck(object):
    """Container for input files that make up a FAST input deck"""

    def __init__(self,fstfile,readlist=['ED','AD'],silent=True):
        """Read FAST master file and read inputs for FAST modules that
        are used
        """
        # Backward compatibility:
        readlist= [rd.replace('Aero', 'AD') for rd in readlist]

        self.filename = fstfile
        self.modeldir = os.path.split(fstfile)[0]
        self.inputfiles = {}
        # read master file
        self.fst = FASTInputFile(fstfile)
        #print('Read',fstfile)
        self.Attributes=['fst']

        def sub_read(file_obj,store_obj,expected_keys,names):
            """ read sub files from object """
            for key,name in zip(expected_keys,names):
                # NOTE: fast input files are relative to their module file
                modeldir = os.path.split(file_obj.filename)[0]
                try:
                    keyval=file_obj[key]
                except:
                    print('Key not found',key)
                    continue
                if type(keyval) is list:
                    setattr(store_obj, name, [])
                    bFileList = True
                else:
                    bFileList = False
                    keyval=[keyval]
                for i,keyv in enumerate(keyval):
                    fpath = os.path.join(modeldir, os.path.normpath(keyv.strip('"').replace('\\','/')))
                    if os.path.isfile(fpath):
                        self.inputfiles[name] = fpath
                        modinput = FASTInputFile(fpath)
                        try:
                            modinput = FASTInputFile(fpath)
                        except:
                            print('Problem reading',fpath)
                        else:
                            if bFileList:
                                getattr(store_obj,name).append(modinput)
                                if not silent:
                                    print('Read',fpath,'as{}[{}]'.format(name,i))
                            else:
                                setattr(store_obj, name, modinput)
                                if not silent:
                                    print('Read',fpath,'as',name)
                                self.Attributes.append(name)
                    else:
                        if not silent:
                            print('Not a file:',fpath)

        FST_Keys= self.fst.keys()
    
        filekeys = [ key for key in FST_Keys if key.endswith('File') ]
        names    = [key[:-4] for key in filekeys]
        sub_read(self.fst,self,filekeys,names)
        # BD
        if 'BD' in readlist:
            try:
                if self.fst['CompElast']==2:
                    filekeys = ['BDBldFile(1)']
                    names    = ['BD']
                    sub_read(self.fst,self,filekeys,names)
            except:
                pass

        if 'InterpOrder' in FST_Keys:
            self.version='OF2'
            if self.fst['CompAero'] == 1:
                self.ADversion='AD14'
            elif self.fst['CompAero'] == 2:
                self.ADversion='AD15'
            else:
                self.ADversion='Unknown'
        elif 'TipRad' in FST_Keys:
            self.version='F7'
        else:
            self.version='Unknown'

        if hasattr(self,'Aero'):
            self.AD=self.Aero
            delattr(self,'Aero')
            self.Attributes = [at.replace('Aero', 'AD') for at in self.Attributes]

        if hasattr(self,'ED') and 'ED' in readlist:
            filekeys = ['BldFile(1)' , 'BldFile(2)' , 'BldFile(3)' , 'TwrFile']
            names    = ['Bld1'       , 'Bld2'       , 'Bld3'       , 'Twr']
            sub_read(self.ED,self.ED,filekeys,names)

        if hasattr(self,'AD') and 'AD' in readlist:
            if 'WakeMod' in self.AD.keys():
                self.ADversion='AD15'
                # NOTE airfoils are not read
                filekeys = ['ADBlFile(1)' , 'ADBlFile(2)' , 'ADBlFile(3)', 'AFNames']
                names    = ['Bld1'     , 'Bld2'     , 'Bld3'    , 'AF']
                sub_read(self.AD,self.AD,filekeys,names)
            else:
                self.ADversion='AD14'

        if self.version=='F7':
            filekeysAD14 = [ key for key in ['BldFile(1)','BldFile(2)','BldFile(3)'] if key in FST_Keys ]
            names        = ['Bld1'     , 'Bld2'     , 'Bld3'   ]
            sub_read(self.fst,self,filekeysAD14,names)


    def __repr__(self):
        s='<weio.FastInputDeck object>'+'\n'
        s+='filename   : '+self.filename+'\n'
        s+='version    : '+self.version+'\n'
        s+='AD version : '+self.ADversion+'\n'
        s+='available attributes: '+','.join(self.Attributes)
        s+='\n'
        return s

if __name__ == "__main__":
    fst=FASTInputDeck('NREL5MW.fst')
    print(fst)
