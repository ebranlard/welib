from .File import File, WrongFormatError, BrokenFormatError, FileNotFoundError, EmptyFileError
from .FileFormats import FileFormat
# User defined formats
from .FASTInFile   import FASTInFile, FASTInputDeck
from .FASTOutFile  import FASTOutFile
from .FASTWndFile  import FASTWndFile
from .FASTLinFile  import FASTLinFile
from .FASTSummaryFile  import FASTSummaryFile
from .CSVFile      import CSVFile
from .HAWC2PCFile  import HAWC2PCFile
from .HAWC2AEFile  import HAWC2AEFile
from .HAWC2DatFile import HAWC2DatFile
from .HAWCStab2PwrFile import HAWCStab2PwrFile
from .HAWCStab2IndFile import HAWCStab2IndFile
from .FLEXBladeFile  import FLEXBladeFile
from .FLEXProfileFile  import FLEXProfileFile
from .FLEXOutFile  import FLEXOutFile
from .ExcelFile  import ExcelFile
from .TurbSimTSFile  import TurbSimTSFile
from .TurbSimFile  import TurbSimFile
from .NetCDFFile   import NetCDFFile
from .TDMSFile   import TDMSFile

class FormatNotDetectedError(Exception):
    pass

def fileFormats():
    formats = []
    formats.append(FileFormat(CSVFile))
    formats.append(FileFormat(ExcelFile))
    formats.append(FileFormat(FASTInFile))
    formats.append(FileFormat(FASTOutFile))
    formats.append(FileFormat(FASTWndFile))
    formats.append(FileFormat(FASTLinFile))
    formats.append(FileFormat(FASTSummaryFile))
    formats.append(FileFormat(TurbSimTSFile))
    formats.append(FileFormat(TurbSimFile))
    formats.append(FileFormat(HAWC2DatFile))
    formats.append(FileFormat(HAWC2PCFile))
    formats.append(FileFormat(HAWC2AEFile))
    formats.append(FileFormat(HAWCStab2PwrFile))
    formats.append(FileFormat(HAWCStab2IndFile))
    formats.append(FileFormat(FLEXBladeFile))
    formats.append(FileFormat(FLEXProfileFile))
    formats.append(FileFormat(FLEXOutFile))
    formats.append(FileFormat(NetCDFFile))
    formats.append(FileFormat(TDMSFile))
    return formats


def detectFormat(filename):
    """ Detect the file formats by looping through the known list. 
        The method may simply try to open the file, if that's the case
        the read file is returned. """
    import os
    formats=fileFormats()
    ext = os.path.splitext(filename.lower())[1]
    detected = False
    i = 0 
    while not detected and i<len(formats):
        myformat = formats[i]
        if ext in myformat.extensions:
            valid, F = myformat.isValid(filename)
            if valid:
                #print('File detected as :',myformat)
                detected=True
                return myformat,F

        i += 1

    if not detected:
        raise FormatNotDetectedError('The file format could not be detected for the file: '+filename)

def read(filename,fileformat=None):
    F = None
    # Detecting format if necessary
    if fileformat is None:
        fileformat,F = detectFormat(filename)
    # Reading the file with the appropriate class if necessary
    if not isinstance(F,fileformat.constructor):
        F=fileformat.constructor(filename=filename)
    return F




