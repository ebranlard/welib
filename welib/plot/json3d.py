""" 

Relies on 
    from welib.FEM.graph import ConnectedObject
"""

import numpy as np
from collections import OrderedDict
import itertools

from welib.FEM.graph import ConnectedObject, displ2mat4


# --------------------------------------------------------------------------------}
# --- JSON3D 
# --------------------------------------------------------------------------------{
class JSON3DFile(dict):
    """ 
    JSON 3D file contains:
      - Nodes
      - Connectivity
      - ElemProp: element properties
      - Time Series
      - Modes

    """
    def __init__(self):
        self._Objects = OrderedDict()

    # --------------------------------------------------------------------------------}
    # --- Method at the Object level
    # --------------------------------------------------------------------------------{
    def addObject(self, obj):
        # TODO Check that object inherits from ConnectedObject
        # --- Sanity checks
        if obj.name in self._Objects.keys():
            raise Exception('Cannot add object `{}`, this name already in list of objects `{}` '.format(obj.name, self._Objects.keys()))
        # Add
        self._Objects[obj.name] = obj

    def addObjectMode(self, objname, displ, modename, freq=1, group='default', **kwargs):
        self._Objects[objname].addMode(displ, modename, freq=freq, group=group, **kwargs)

    def addObjectTimeSeries(self, objname, time, **kwargs):
        self._Objects[objname].addTimeSeries(time, **kwargs)

    # --------------------------------------------------------------------------------}
    # --- Time Series
    # --------------------------------------------------------------------------------{
    def initTimeSeries(self, time):
        pass
    # --------------------------------------------------------------------------------}
    # --- Main properties, cnocatenantion of all objects
    # --------------------------------------------------------------------------------{
    @property
    def Nodes(self):
        """ 
        Combine all Nodes from all objects
        """
        for i,k in enumerate(self._Objects.keys()):
            if i==0:
                Nodes = self._Objects[k].Nodes
            else:
                Nodes =np.vstack( (Nodes, self._Objects[k].Nodes) )
        return Nodes

    @property
    def Connectivity(self):
        """ 
        Combine all connectivity from all objects
        """
        iNode=0
        for i,k in enumerate(self._Objects.keys()):
            if i==0:
                Connectivity = self._Objects[k].Connectivity
                iNode=iNode+self._Objects[k].nNodes
            else:
                C =  self._Objects[k].Connectivity
                for l in C:
                    Connectivity.append([e+iNode for e in l])

        # --- Sanity
        allcon = list(itertools.chain.from_iterable(Connectivity))
        iMin = np.min(allcon)
        if iMin!=0:
            raise Exception('Connectivity table needs to have 0 as a minimum  node index')
        # Max number of nodes
        iMax = np.max(allcon)
        nNodes = len(self.Nodes)
        if iMax!=nNodes-1:
            raise Exception('Connectivity table needs to have {} as a maximum, corresponding to the maximum number of nodes'.format(nNodes-1))
        # Check that all nodes indices are present
        IC = np.unique(allcon)
        IA = np.arange(0,nNodes)
        ID = np.setdiff1d(IA,IC)
        if len(ID)>0:
            raise Exception('Connectivity table is missing the following indices : {}.'.format(ID))
        return Connectivity

    @property
    def ElemProps(self):
        """ 
        Combine all element properties from all objects
        """
        ElemProps=[p for key in self._Objects.keys() for p in self._Objects[key].ElemProps]
        return ElemProps

    def Modes(self, digits=2):
        """ 
        Combine all modes from all objects
        """
        allModeNames = np.unique([m['name'] for k,o in self._Objects.items()  for km, m in o.Modes.items()])
        allGroups    = np.unique([m['group'] for k,o in self._Objects.items() for km, m in o.Modes.items()])
        Modes=dict()
        for iMode,modeName in enumerate(allModeNames):
            displ, freq, group =  self.Mode(modeName, digits=digits)
            if group not in Modes.keys():
                Modes[group] = []
            Modes[group].append({'name':modeName, 'Displ':displ.tolist(), 'freq':freq})
        return Modes
    
    
    def Mode(self, modename, digits=2):
        """ Combine all nodal displacements for a given mode"""
        # Concatenate displacements
        for i,(k,o) in enumerate(self._Objects.items()):
            if modename not in o.Modes.keys():
                raise Exception('Object {} does not have mode {}'.format(o.name, modename))
            mode = o.Modes[modename]
            freq  = mode['freq']  # TODO which objet holds this..
            if i==0:
                group0 = mode['group']
                displ = np.asarray(mode['displ'])
            else:
                displ =np.vstack( (displ, np.asarray(mode['displ']) ) )
                group = mode['group']
                if group!=group0:
                    raise Exception('Object {}, mode {} has group {} instead of {}'.format(o.name, modename, group, group0))
        displ = np.around(displ, digits)
        return displ, freq, group0


    def TimeSeries(self, digits=2, flat=False, mat4=False):
        """ Combine all timeseries"""
        allTSNames = np.unique([ts['name']  for k,o in self._Objects.items() for kts, ts in o.TimeSeries.items()])
        allGroups  = np.unique([ts['group'] for k,o in self._Objects.items() for kts, ts in o.TimeSeries.items()])
        TS=dict()
        for _,tsName in enumerate(allTSNames):
            displ, m4, timeInfo, absolute, perElement, group =  self.TS(tsName, digits=digits)
            if group not in TS.keys():
                TS[group] = []
            ts= {'name':tsName, 'timeInfo':timeInfo, 'absolute':absolute, 'perElement':perElement}#, 'flat':flat}
            if mat4:
#                 if flat:
#                     raise NotImplementedError()
#                     m4=m4.flatten()
                ts['mat4']=m4.tolist()
            else:
#                 if flat:
#                     raise NotImplementedError()
#                     displ=displ.flatten()
                ts['Displ']= displ.tolist()
            TS[group].append(ts)
        return TS


    def TS(self, tsName, digits=2):
        """ Combine all nodal displacements for a given timeseries"""
        # Concatenate displacements
        for i,(k,o) in enumerate(self._Objects.items()):
            if tsName not in o.TimeSeries.keys():
                raise Exception('Object {} does not have ts {}'.format(o.name, tsName))
            ts = o.TimeSeries[tsName]
            absolute   = ts['absolute']   # TODO which objet holds this..
            perElement = ts['perElement'] # TODO which objet holds this..
            timeInfo   = ts.timeInfo      # TODO...
            if i==0:
                group0 = ts['group']
                displ = np.asarray(ts['displ'])
            else:
                displ0 = np.asarray(ts['displ'])
                displ  = np.concatenate( (displ,displ0), axis=1 )
                group = ts['group']
                if group!=group0:
                    raise Exception('Object {}, ts {} has group {} instead of {}'.format(o.name, tsName, group, group0))
        # Convert to mat4 only at the end for now. 
        # Concatenation of mat4 todo
        mat4  = np.around(displ2mat4(displ, flat = True), digits)
        displ = np.around(displ, digits)

        return displ, mat4, timeInfo, absolute, perElement, group0


    # --------------------------------------------------------------------------------}
    # --- Main Getters/Actions
    # --------------------------------------------------------------------------------{
    def toDict(self, digits=2, flat=False, mat4=False):
        d=dict()
        d['Nodes']        = self.Nodes.tolist()
        d['Connectivity'] = self.Connectivity
        d['ElemProps']    = self.ElemProps
        d['Modes']        = self.Modes(digits=digits)
        d['TimeSeries']   = self.TimeSeries(digits=digits, flat=flat, mat4=mat4)
        return d

    def write(self, jsonFile, digits=2, flat=False, mat4=False):
        import json
        from io import open
        d = self.toDict(digits=digits, flat=flat, mat4=mat4)
        # --- Sanity
        for ep in d['ElemProps']:
            if not all(key in ep for key in ['shape','type']):
                raise Exception('Element property must have `shape` and `type` as keys: {} '.format(ep))
            if ep['shape'] not in ['cylinder']:
                raise NotImplementedError()

#         allGroups= np.unique([ts['group'] for ts in self.TimeSeries]) 
#         d['TimeSeries']=dict()
#         for g in allGroups:
#             d['TimeSeries'][g] = [{ 
#                 'name':       TS['name'],
#                 'timeInfo':   TS.timeInfo,
#                 'absolute':   TS['absolute'],
#                 'perElement': TS['perElement'],
#                 'mat4':     TS.mat4(flat=True),
#                 }  for iTS,TS in enumerate(self.TimeSeries) if TS['group']==g]
# 
#         d['groundLevel']=np.min(Points[:,2]) # TODO

        # Writing 
        with open(jsonFile, 'w', encoding='utf-8') as f:
            try:
                #f.write(unicode(json.dumps(d, ensure_ascii=False))) #, indent=2)
                #f.write(json.dumps(d, ensure_ascii=False)) #, indent=2)
                f.write(json.dumps(d, ensure_ascii=False))
            except:
                print('>>> FAILED')
                json.dump(d, f, indent=0) 

    def __repr__(self):
        s='<{} object> with attributes:\n'.format(type(self).__name__)
        s+='- _Objects,        odict with keys: ({}):\n'.format(self._Objects.keys())
        s+='Objects:\n'
        for k in self._Objects.keys():
            s+=str(self._Objects[k])
        return s


