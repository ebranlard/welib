from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from io import open
from builtins import map
from builtins import range
from builtins import chr
from future import standard_library
standard_library.install_aliases()
import os
import numpy as np
import struct
from itertools import takewhile

def load_output(filename):
    """Load a FAST binary or ascii output file

    Parameters
    ----------
    filename : str
        filename

    Returns
    -------
    data : ndarray
        data values
    info : dict
        info containing:
            - name: filename
            - description: description of dataset
            - attribute_names: list of attribute names
            - attribute_units: list of attribute units
    """

    assert os.path.isfile(filename), "File, %s, does not exists" % filename
    with open(filename, 'r') as f:
        try:
            f.readline()
        except UnicodeDecodeError:
            return load_binary_output(filename)
    return load_ascii_output(filename)

def load_ascii_output(filename):
    with open(filename) as f:
        info = {}
        info['name'] = os.path.splitext(os.path.basename(filename))[0]
        # Header is whatever is before the keyword `time`
        in_header = True
        header = []
        while in_header:
            l = f.readline()
            if not l:
                raise Exception('Error finding the end of FAST out file header. Keyword Time missing.')
            in_header= (l+' dummy').lower().split()[0] != 'time'
            if in_header:
                header.append(l)
            else:
                info['description'] = header
                info['attribute_names'] = l.split()
                info['attribute_units'] = [unit[1:-1] for unit in f.readline().split()]

        # Data, up to end of file or empty line (potential comment line at the end)
        data = np.array([l.strip().split() for l in takewhile(lambda x: len(x.strip())>0, f.readlines())]).astype(np.float)
        return data, info


def load_binary_output(filename, use_buffer=True):
    """
    03/09/15: Ported from ReadFASTbinary.m by Mads M Pedersen, DTU Wind
    24/10/18: Low memory/buffered version by E. Branlard, NREL

    Info about ReadFASTbinary.m:
    % Author: Bonnie Jonkman, National Renewable Energy Laboratory
    % (c) 2012, National Renewable Energy Laboratory
    %
    %  Edited for FAST v7.02.00b-bjj  22-Oct-2012
    """
    def fread(fid, n, type):
        fmt, nbytes = {'uint8': ('B', 1), 'int16':('h', 2), 'int32':('i', 4), 'float32':('f', 4), 'float64':('d', 8)}[type]
        return struct.unpack(fmt * n, fid.read(nbytes * n))

    def freadRowOrderTableBuffered(fid, n, type_in, nCols, nOff=0, type_out='float64'):
        """ 
        Reads of row-ordered table from a binary file.

        Read `n` data of type `type_in`, assumed to be a row ordered table of `nCols` columns.
        Memory usage is optimized by allocating the data only once.
        Buffered reading is done for improved performances (in particular for 32bit python)

        `nOff` allows for additional column space at the begining of the storage table.
        Typically, `nOff=1`, provides a column at the beginning to store the time vector.

        @author E.Branlard, NREL

        """
        fmt, nbytes = {'uint8': ('B', 1), 'int16':('h', 2), 'int32':('i', 4), 'float32':('f', 4), 'float64':('d', 8)}[type_in]
        nLines          = int(n/nCols)
        GoodBufferSize  = 4096*40
        nLinesPerBuffer = int(GoodBufferSize/nCols)
        BufferSize      = nCols * nLinesPerBuffer
        nBuffer         = int(n/BufferSize)
        # Allocation of data
        data = np.zeros((nLines,nCols+nOff), dtype = type_out)
        # Reading
        try:
            nIntRead   = 0
            nLinesRead = 0
            while nIntRead<n:
                nIntToRead = min(n-nIntRead, BufferSize)
                nLinesToRead = int(nIntToRead/nCols)
                Buffer = np.array(struct.unpack(fmt * nIntToRead, fid.read(nbytes * nIntToRead)))
                Buffer = Buffer.reshape(-1,nCols)
                data[ nLinesRead:(nLinesRead+nLinesToRead),  nOff:(nOff+nCols)  ] = Buffer
                nLinesRead = nLinesRead + nLinesToRead
                nIntRead   = nIntRead   + nIntToRead
        except:
            raise Exception('Read only %d of %d values in file:' % (nIntRead, n, filename))
        return data


    FileFmtID_WithTime = 1  #% File identifiers used in FAST
    FileFmtID_WithoutTime = 2
    LenName = 10  #;  % number of characters per channel name
    LenUnit = 10  #;  % number of characters per unit name

    with open(filename, 'rb') as fid:
        FileID = fread(fid, 1, 'int16')[0]  #;             % FAST output file format, INT(2)
        if FileID not in [FileFmtID_WithTime, FileFmtID_WithoutTime]:
            raise Exception('FileID not supported {}. Is it a FAST binary file?'.format(FileID))

        NumOutChans = fread(fid, 1, 'int32')[0]  #;             % The number of output channels, INT(4)
        NT = fread(fid, 1, 'int32')[0]  #;             % The number of time steps, INT(4)


        if FileID == FileFmtID_WithTime:
            TimeScl = fread(fid, 1, 'float64')  #;           % The time slopes for scaling, REAL(8)
            TimeOff = fread(fid, 1, 'float64')  #;           % The time offsets for scaling, REAL(8)
        else:
            TimeOut1 = fread(fid, 1, 'float64')  #;           % The first time in the time series, REAL(8)
            TimeIncr = fread(fid, 1, 'float64')  #;           % The time increment, REAL(8)




        ColScl = fread(fid, NumOutChans, 'float32')  #; % The channel slopes for scaling, REAL(4)
        ColOff = fread(fid, NumOutChans, 'float32')  #; % The channel offsets for scaling, REAL(4)

        LenDesc = fread(fid, 1, 'int32')[0]  #;  % The number of characters in the description string, INT(4)
        DescStrASCII = fread(fid, LenDesc, 'uint8')  #;  % DescStr converted to ASCII
        DescStr = "".join(map(chr, DescStrASCII)).strip()



        ChanName = []  # initialize the ChanName cell array
        for iChan in range(NumOutChans + 1):
            ChanNameASCII = fread(fid, LenName, 'uint8')  #; % ChanName converted to numeric ASCII
            ChanName.append("".join(map(chr, ChanNameASCII)).strip())


        ChanUnit = []  # initialize the ChanUnit cell array
        for iChan in range(NumOutChans + 1):
            ChanUnitASCII = fread(fid, LenUnit, 'uint8')  #; % ChanUnit converted to numeric ASCII
            ChanUnit.append("".join(map(chr, ChanUnitASCII)).strip()[1:-1])


        #    %-------------------------
        #    % get the channel time series
        #    %-------------------------

        nPts = NT * NumOutChans  #;           % number of data points in the file
        #print('NT',NT)
        #print('NumOutChans',NumOutChans)


        if FileID == FileFmtID_WithTime:
            PackedTime = fread(fid, NT, 'int32')  #; % read the time data
            cnt = len(PackedTime)
            if cnt < NT:
                raise Exception('Could not read entire %s file: read %d of %d time values' % (filename, cnt, NT))

        if use_buffer:
            # Reading data using buffers, and allowing an offset for time column (nOff=1)
            data = freadRowOrderTableBuffered(fid, nPts, 'int16', NumOutChans, nOff=1, type_out='float64')
        else:
            # NOTE: unpacking huge data not possible on 32bit machines
            PackedData = fread(fid, nPts, 'int16')  #; % read the channel data
            cnt = len(PackedData)
            if cnt < nPts:
                raise Exception('Could not read entire %s file: read %d of %d values' % (filename, cnt, nPts))
            data = np.array(PackedData).reshape(NT, NumOutChans)
            del PackedData

    if FileID == FileFmtID_WithTime:
        time = (np.array(PackedTime) - TimeOff) / TimeScl;
    else:
        time = TimeOut1 + TimeIncr * np.arange(NT)

    #import pdb
    #pdb.set_trace()

    #    %-------------------------
    #    % Scale the packed binary to real data
    #    %-------------------------
    if use_buffer:
        # Scaling Data
        for iCol in range(NumOutChans):
            if np.isnan(ColScl[iCol]) and np.isnan(ColOff[iCol]):
                data[:,iCol+1] = 0 # probably due to a division by zero in Fortran
            else:
                data[:,iCol+1] = (data[:,iCol+1] - ColOff[iCol]) / ColScl[iCol]
        # Adding time column
        data[:,0] = time
    else:
        # NOTE: memory expensive due to time conversion, and concatenation
        data = (data - ColOff) / ColScl
        data = np.concatenate([time.reshape(NT, 1), data], 1)

    info = {'name': os.path.splitext(os.path.basename(filename))[0],
            'description': DescStr,
            'attribute_names': ChanName,
            'attribute_units': ChanUnit}
    return data, info

