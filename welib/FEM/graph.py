""" 
Basics Classes for a "geometrical" graph model: 
  - nodes have a position x,y (,z) and some properties
  - links connect nodes

"""

import numpy as np
import pandas as pd


print('WARNING, using graph from welib.FEM, right now weio has priority')


# --------------------------------------------------------------------------------}
# --- Node
# --------------------------------------------------------------------------------{
class Node(object):
    def __init__(self, ID, x, y, z=0, **kwargs):
        self.ID = int(ID)
        self.x  = x
        self.y  = y
        self.z  = z
        self.data  = kwargs

    def setData(self, data_dict):
        for k,v in data_dict.items():
            #if k in self.data.keys():
            #    print('Warning overriding key {} for node {}'.format(k,self.ID))
            self.data[k]=v

    def __repr__(self):
        s='<Node{:4d}> x:{:7.2f} y:{:7.2f} z:{:7.2f} {:}'.format(self.ID, self.x, self.y, self.z, self.data)
        return s

# --------------------------------------------------------------------------------}
# --- Properties  
# --------------------------------------------------------------------------------{
class Property(dict):
    def __init__(self, ID, data=None, **kwargs):
        """ 
        data is a dictionary
        """
        dict.__init__(self)
        self['ID']= int(ID)
        self.update(kwargs)
        if data is not None:
            self.upadte(data)

    @property
    def data(self):
        return {k:v for k,v in self.items() if k is not 'ID'}

    def __repr__(self):
        s='<Prop{:4d}> {:}'.format(self['ID'], self.data)
        return s

class NodeProperty(Property):
    def __init__(self, ID, data=None, **kwargs):
        Property.__init__(self, ID, data, **kwargs)
    def __repr__(self):
        s='<NPrp{:4d}> {:}'.format(self['ID'], self.data)
        return s
    
class ElemProperty(Property):
    def __init__(self, ID, data=None, **kwargs):
        Property.__init__(self, ID, data, **kwargs)
    def __repr__(self):
        s='<EPrp{:4d}> {:}'.format(self['ID'], self.data)
        return s


# --------------------------------------------------------------------------------}
# --- Elements 
# --------------------------------------------------------------------------------{
class Element(dict):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        self.ID      = int(ID)
        self.nodeIDs = nodeIDs
        self.nodes   = nodes
        self.propset = propset
        self.propIDs = propIDs
        self.prop    = properties # Nodal properties
        self.data    = kwargs     # Nodal data
        if (self.propIDs is not None) and (self.propset is None):
            raise Exception('propset should be provided if propIDs are provided')
        
    @property
    def length(self):
        n1=self.nodes[0]
        n2=self.nodes[1]
        return np.sqrt((n1.x-n2.x)**2+(n1.y-n2.y)**2+(n1.z-n2.z)**2)

    def __repr__(self):
        s='<Elem{:4d}> NodeIDs: {} {}'.format(self.ID, self.nodeIDs, self.data)
        if self.nodes is not None:
            s+=' l={:.2f}'.format(self.length)
        return s


# --------------------------------------------------------------------------------}
# --- Mode 
# --------------------------------------------------------------------------------{
class Mode(dict):
    def __init__(self, data, name, freq=1, **kwargs):
        dict.__init__(self)

        self['name']=name
        self['freq']=freq
        self['data']=data # displacements nNodes x 3 assuming a given sorting of nodes

    def __repr__(self):
        s='<Mode> name:{:4s} freq:{:} '.format(self['name'], self['freq'])
        return s

    def reSort(self,I):
        self['data']=self['data'][I,:]

# --------------------------------------------------------------------------------}
# --- Graph
# --------------------------------------------------------------------------------{
class GraphModel(object):
    def __init__(self):
        self.Elements    = []
        self.Nodes       = []
        self.NodePropertySets= dict()
        self.ElemPropertySets= dict()
        self.MiscPropertySets= dict()
        # Dynamics
        self.Modes   = []
        self.Motions = []

    def addNode(self,node):
        self.Nodes.append(node)

    def addElement(self,elem):
        if elem.nodes is None:
            # Giving nodes to element if these were not provided
            elem.nodes=[self.getNode(i) for i in elem.nodeIDs]
#         if elem.prop is None:
#             # Giving props to element if these were not provided
#             if elem.propIDs is not None:
#                 elem.prop=[self.getElemProperty(elem.propset, i) for i in elem.propIDs]
        self.Elements.append(elem)

    # --- Getters
    def getNode(self,nodeID):
        for n in self.Nodes:
            if n.ID==nodeID:
                return n
        raise KeyError('NodeID {} not found in Nodes'.format(nodeID))

    def getElement(self, elemID):
        for e in self.Elements:
            if e.ID==elemID:
                return e
        raise KeyError('ElemID {} not found in Elements'.format(elemID))

    def getNodeProperty(self, setname, propID):
        for p in self.NodePropertySets[setname]:
            if p['ID']==propID:
                return p
        raise KeyError('PropID {} not found for Node propset {}'.format(propID,setname))

    def getElementProperty(self, setname, propID):
        for p in self.ElemPropertySets[setname]:
            if p['ID']==propID:
                return p
        raise KeyError('PropID {} not found for Element propset {}'.format(propID,setname))

    def getMiscProperty(self, setname, propID):
        for p in self.MiscPropertySets[setname]:
            if p['ID']==propID:
                return p
        raise KeyError('PropID {} not found for Misc propset {}'.format(propID,setname))

    # --- Properties
    def addElementPropertySet(self, setname):
        self.ElemPropertySets[setname]= []

    def addNodePropertySet(self, setname):
        self.NodePropertySets[setname]= []

    def addMiscPropertySet(self, setname):
        self.MiscPropertySets[setname]= []

    def addNodeProperty(self, setname, prop):
        if not isinstance(prop, NodeProperty):
            print(type(prop))
            raise Exception('Property needs to inherit from NodeProperty')
        self.PropertySets[setname].append(prop)

    def addNodeProperty(self, setname, prop):
        if not isinstance(prop, NodeProperty):
            print(type(prop))
            raise Exception('Property needs to inherit from NodeProperty')
        self.NodePropertySets[setname].append(prop)

    def addElementProperty(self, setname, prop):
        if not isinstance(prop, ElemProperty):
            print(type(prop))
            raise Exception('Property needs to inherit from ElementProperty')
        self.ElemPropertySets[setname].append(prop)

    def addMiscProperty(self, setname, prop):
        if not isinstance(prop, ElemProperty):
            print(type(prop))
            raise Exception('Property needs to inherit from Property')
        self.MiscPropertySets[setname].append(prop)

    # --- Data and node and element prop setters
    def setElementNodalProp(self, elem, propset, propIDs):
        """ 
        Set Nodal Properties to each node of an element
        """
        for node, pID in zip(elem.nodes, propIDs):
            node.setData(self.getNodeProperty(propset, pID).data)

    def setNodeNodalProp(self, node, propset, propID):
        """ 
        Set Nodal Properties to a node
        """
        node.setData(self.getNodeProperty(propset, propID).data)

    def setNodalData(self, nodeID, **data_dict):
        self.getNode(nodeID).setData(data_dict)

    def __repr__(self):
        s='<{} object> with keys:\n'.format(type(self).__name__)
        s+='- Nodes ({}):\n'.format(len(self.Nodes))
        s+='\n'.join(str(n) for n in self.Nodes)
        s+='\n- Elements ({}):\n'.format(len(self.Elements))
        s+='\n'.join(str(n) for n in self.Elements)
        s+='\n- NodePropertySets ({}):'.format(len(self.NodePropertySets))
        for k,v in self.NodePropertySets.items():
            s+='\n> {} ({}):\n'.format(k, len(v))
            s+='\n'.join(str(p) for p in v)
        s+='\n- ElementPropertySets ({}):'.format(len(self.ElemPropertySets))
        for k,v in self.ElemPropertySets.items():
            s+='\n> {} ({}):\n'.format(k, len(v))
            s+='\n'.join(str(p) for p in v)
        s+='\n- MiscPropertySets ({}):'.format(len(self.MiscPropertySets))
        for k,v in self.MiscPropertySets.items():
            s+='\n> {} ({}):\n'.format(k, len(v))
            s+='\n'.join(str(p) for p in v)
        s+='\n- Modes ({}):\n'.format(len(self.Modes))
        s+='\n'.join(str(m) for m in self.Modes)
        s+='\n- Motions ({}):'.format(len(self.Motions))
        for m in self.Motions:
            s+='\n> {}\n'.format({k:v for k,v in m.items() if not isintance(v,np.ndarray)})
        return s

    # --------------------------------------------------------------------------------}
    # --- Geometrical properties 
    # --------------------------------------------------------------------------------{
    @property
    def extent(self):
        xmax=np.max([node.x for node in self.Nodes])
        ymax=np.max([node.y for node in self.Nodes])
        zmax=np.max([node.z for node in self.Nodes])
        xmin=np.min([node.x for node in self.Nodes])
        ymin=np.min([node.y for node in self.Nodes])
        zmin=np.min([node.z for node in self.Nodes])
        return [xmin,ymin,zmin],[xmax,ymax,zmax],[xmax-xmin,ymax-ymin,zmax-zmin]

    @property
    def maxDimension(self):
        _,_,D=self.extent
        return np.max(D)

    @property
    def points(self):
        nNodes = len(self.Nodes)
        Points = np.zeros((nNodes,3))
        for i,n in enumerate(self.Nodes):
            Points[i,:]=(n.x, n.y, n.z)
        return Points

    @property
    def connectivity(self):
        """ returns connectivity, assuming points are indexed starting at 0 """
        conn=[]
        for ie, e in enumerate(self.Elements):
            elemConn=[]
            for n in e.nodes:
                elemConn.append(self.Nodes.index(n))
            conn.append(elemConn)
        return conn

    def toLines(self, output='coord'):
        if output=='coord':
            lines = np.zeros((len(self.Elements), 2, 3)) # 
            for ie, e in enumerate(self.Elements):
                n1=e.nodes[0]
                n2=e.nodes[-1]
                lines[ie, 0, : ] = (n1.x, n1.y, n1.z)
                lines[ie, 1, : ] = (n2.x, n2.y, n2.z)
        elif output=='lines3d':
            import mpl_toolkits.mplot3d as plt3d
            lines=[]
            for ie, e in enumerate(self.Elements):
                n1=e.nodes[0]
                n2=e.nodes[-1]
                line = plt3d.art3d.Line3D((n1.x,n2.x), (n1.y,n2.y), (n1.z,n2.z))
                lines.append(line)
        else:
            raise NotImplementedError()

        return lines

    # --------------------------------------------------------------------------------}
    # ---  
    # --------------------------------------------------------------------------------{
    def divideElement(self, elemID, nPerElement):
        """ divide a given element by nPerElement (add nodes and elements to graph) """ 
        if len(self.Modes)>0:
            raise Exception('Cannot divide graph when mode data is present')
        if len(self.Motions)>0:
            raise Exception('Cannot divide graph when motion data is present')


        maxNodeId=np.max([n.ID for n in self.Nodes])
        maxElemId=np.max([e.ID for e in self.Nodes])
        e = self.getElement(elemID)
        if len(e.nodes)==2:
            n1=e.nodes[0]
            n2=e.nodes[1]
            subNodes=[n1]
            for iSub in range(1,nPerElement):
                maxNodeId += 1
                #data_dict  = n1.data.copy()
                data_dict  = dict()
                fact       = float(iSub)/nPerElement
                x          = n1.x*(1-fact)+n2.x*fact
                y          = n1.y*(1-fact)+n2.y*fact
                z          = n1.z*(1-fact)+n2.z*fact
                for k,v in n1.data.items():
                    try:
                        data_dict[k] = n1.data[k]*(1-fact) + n2.data[k]*fact
                    except:
                        pass
                ni = Node(maxNodeId, x, y, z, **data_dict)
                subNodes.append(ni)
                self.addNode(ni)
            subNodes+=[n2]
            e.nodes  =subNodes[0:2]
            e.nodeIDs=[e.ID for e in e.nodes]
            for i in range(1,nPerElement):
                maxElemId+=1
                elem_dict = e.data.copy()
                elem= Element(maxElemId, [subNodes[i].ID, subNodes[i+1].ID], **elem_dict )
                self.addElement(elem)

    def sortNodesBy(self,key):
        """ 
        z
        """

        # TODO, that's quite doable
        if len(self.Modes)>0:
            raise Exception('Cannot sort nodes when mode data is present')
        if len(self.Motions)>0:
            raise Exception('Cannot sort nodes when motion data is present')

        nNodes = len(self.Nodes)
        if key=='x':
            values=[n.x for n in self.Nodes]
        elif key=='y':
            values=[n.y for n in self.Nodes]
        elif key=='z':
            values=[n.z for n in self.Nodes]
        elif key=='ID':
            values=[n.ID for n in self.Nodes]
        else:
            values=[n[key] for n in self.Nodes]
        I= np.argsort(values)
        self.Nodes=[self.Nodes[i] for i in I]

    def divideElements(self, nPerElement):
        """ divide all elements by nPerElement (add nodes and elements to graph) """ 
        maxNodeId=np.max([n.ID for n in self.Nodes])
        maxElemId=np.max([e.ID for e in self.Nodes])

        if nPerElement<=0:
            raise Exception('nPerElement should be more than 0')

        for ie in np.arange(len(self.Elements)): # cannot enumerate since length increases
            elemID = self.Elements[ie].ID
            self.divideElement(elemID, nPerElement)
                    
    # --------------------------------------------------------------------------------}
    # --- Dynamics
    # --------------------------------------------------------------------------------{
    def addMode(self,displ,name=None,freq=1):
        if name is None:
            name='Mode '+str(len(self.Modes))
        mode = Mode(data=displ, name=name, freq=freq)
        self.Modes.append(mode)


    # --------------------------------------------------------------------------------}
    # --- Ouputs / converters
    # --------------------------------------------------------------------------------{
    def nodalDataFrame(self, sortBy=None):
        """ return a DataFrame of all the nodal data """
        data=dict()
        nNodes=len(self.Nodes)
        for i,n in enumerate(self.Nodes):
            if i==0:
                data['ID'] = np.zeros(nNodes).astype(int)
                data['x']  = np.zeros(nNodes)
                data['y']  = np.zeros(nNodes)
                data['z']  = np.zeros(nNodes)

            data['ID'][i] = n.ID
            data['x'][i]  = n.x
            data['y'][i]  = n.y
            data['z'][i]  = n.z
            for k,v in n.data.items():
                if k not in data:
                    data[k] = np.zeros(nNodes)
                try:
                    data[k][i]=v
                except:
                    pass
        df = pd.DataFrame(data)
        # Sorting 
        if sortBy is not None:
            df.sort_values([sortBy],inplace=True,ascending=True)
            df.reset_index(drop=True,inplace=True) 
        return df


    def toJSON(self,outfile=None):
        d=dict();
        Points=self.points
        d['Connectivity'] = self.connectivity
        d['Nodes']        = Points.tolist()
        
        d['ElemProps']=list()
        for iElem,elem in enumerate(self.Elements):
            Shape = elem.data['shape'] if 'shape' in elem.data.keys() else 'cylinder'
            Type  = elem.data['Type'] if 'Type' in elem.data.keys() else 1
            Diam  = elem.data['D'] if 'D' in elem.data.keys() else 1
            if Shape=='cylinder':
                d['ElemProps'].append({'shape':'cylinder','type':Type, 'Diam':Diam})
            else:
                raise NotImplementedError()


        d['Modes']=[
                {
                    'name': self.Modes[iMode]['name'],
                    'omega':self.Modes[iMode]['freq']*2*np.pi, #in [rad/s]
                    'Displ':self.Modes[iMode]['data'].tolist()
                }  for iMode,mode in enumerate(self.Modes)]
        d['groundLevel']=np.min(Points[:,2]) # TODO

        if outfile is not None:
            import json
            from io import open
            jsonFile=outfile
            with open(jsonFile, 'w', encoding='utf-8') as f:
                try:
                    #f.write(unicode(json.dumps(d, ensure_ascii=False))) #, indent=2)
                    #f.write(json.dumps(d, ensure_ascii=False)) #, indent=2)
                    f.write(json.dumps(d, ensure_ascii=False))
                except:
                    print('>>> FAILED')
                    json.dump(d, f, indent=0) 
        return d

# 
