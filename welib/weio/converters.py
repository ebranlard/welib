import os
from pathlib import Path
import argparse



# --------------------------------------------------------------------------------
# --- Writing pandas DataFrame to different formats
# --------------------------------------------------------------------------------
def writeDataFrameToFormat(df, filename, fformat):
    """  
    Write a dataframe to disk based on user-specified fileformat
    - df: pandas dataframe
    - filename: filename 
    - fformat: fileformat in: ['csv', 'outb', 'parquet']
    """

    if fformat=='outb':
        dataFrameToOUTB(df, filename)
    elif fformat=='parquet':
        dataFrameToParquet(df, filename)
    elif fformat=='csv':
        dataFrameToCSV(df, filename, sep=',', index=False)
    else:
        raise Exception('File format not supported for dataframe export `{}`'.format(fformat))

def writeDataFrameAutoFormat(df, filename, fformat=None):
    """ 
    Write a dataframe to disk based on extension
    - df: pandas dataframe
    - filename: filename 
    """
    if fformat is not None:
        raise Exception()
    base, ext = os.path.splitext(filename)
    ext = ext.lower()
    if ext in ['.outb']:
        fformat = 'outb'
    elif ext in ['.parquet', '.pq']:
        fformat = 'parquet'
    elif ext in ['.csv']:
        fformat = 'csv'
    else:
        print('[WARN] defaulting to csv, extension unknown: `{}`'.format(ext))
        fformat = 'csv'

    writeDataFrameToFormat(df, filename, fformat)

def writeFileDataFrames(fileObject, writer, extension='.conv', filename=None, **kwargs):
    """ 
    From a fileObejct, extract dataframes and write them to disk.

    - fileObject: object inheriting from weio.File with at least
                   - the attributes .filename
                   - the method     .toDataFrame()
    - writer: function with the interface:   writer ( dataframe, filename, **kwargs )
    """ 
    if filename is None:
        base, _ = os.path.splitext(fileObject.filename)
        filename = base + extension
    else:
        base, ext = os.path.splitext(filename)
        if len(ext)!=0:
            extension = ext
    if filename == fileObject.filename:
        raise Exception('Not overwritting {}. Specify a filename or an extension.'.format(filename))
        
    dfs = fileObject.toDataFrame()
    if isinstance(dfs, dict):
        for name,df in dfs.items():
            filename = base + name + extension
            if filename == fileObject.filename:
                raise Exception('Not overwritting {}. Specify a filename or an extension.'.format(filename))
            writeDataFrame(df=df, writer=writer, filename=filename, **kwargs)
    else:
        writeDataFrame(df=dfs, writer=writer, filename=filename, **kwargs)

def writeDataFrame(df, writer, filename, **kwargs):
    """ 
    Write a dataframe to disk based on a "writer" function. 
    - df: pandas dataframe
    - writer: function with the interface:   writer ( dataframe, filename, **kwargs )
    - filename: filename 
    """
    writer(df, filename, **kwargs)

# --- Low level writers
def dataFrameToCSV(df, filename, sep=',', index=False, **kwargs):
    df.to_csv(filename, sep=sep, index=index, **kwargs)

def dataFrameToOUTB(df, filename, **kwargs):
    from .fast_output_file import writeDataFrame as writeDataFrameToOUTB
    writeDataFrameToOUTB(df, filename, binary=True)

def dataFrameToParquet(df, filename, **kwargs):
    df.to_parquet(path=filename, **kwargs)





# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def fast_input_file_standardize(paths='./', target_extensions='', keys_for_exclusions='', overwrite=False, output_base_dir=None, verbose=True):
    from .fast_input_file import FASTInputFile
    from welib.tools.strings import OK, FAIL, INFO, WARN
    if isinstance(target_extensions, str):
        target_extensions = [ext.strip() for ext in target_extensions.split(',') if ext.strip()]
    if isinstance(keys_for_exclusions, str):
        keys_for_exclusions = [key.strip() for key in keys_for_exclusions.split(',') if key.strip()]
        
    print('paths:', paths)
    print('Target_extensions:', target_extensions)
    print('Keys_for_exclusions:', keys_for_exclusions)
    print('')
    
    if isinstance(paths, (str, Path)):
        path_str = str(paths)
        if ',' in path_str and not Path(path_str).exists():
            path_items = [Path(p.strip()) for p in path_str.split(',') if p.strip()]
        else:
            path_items = [Path(path_str)]
    else:
        path_items = [Path(p) for p in paths]
        
    if not overwrite and output_base_dir:
        out_base = Path(output_base_dir)
        out_base.mkdir(parents=True, exist_ok=True)

    collected_files = []
    for item in path_items:
        if not item.exists():
            FAIL(f"Path does not exist: {item}")
            continue
        if item.is_file():
            collected_files.append((item, item.parent))
        elif item.is_dir():
            for file_path in item.rglob('*'):
                if file_path.is_file():
                    collected_files.append((file_path, item))

    nFiles = len(collected_files)
    for file_path, base_root in collected_files:
        if file_path.suffix.lower() in [ext.lower() for ext in target_extensions]:
            try:
                if verbose:
                    INFO('Reading: ', file_path)
                fid = FASTInputFile(str(file_path))
                has_excluded_key = False
                if hasattr(fid, 'keys'):
                    has_excluded_key = any(key in fid.keys() for key in keys_for_exclusions)
                
                if not has_excluded_key:
                    if overwrite:
                        filename_out = str(file_path)
                    else:
                        if output_base_dir:
                            rel_path = file_path.relative_to(base_root)
                            filename_out = str(out_base / rel_path)
                            Path(filename_out).parent.mkdir(parents=True, exist_ok=True)
                        else:
                            filename_out = str(file_path.with_name(f"{file_path.stem}_out{file_path.suffix}"))
                    
                    fid.write(filename_out)
                    if verbose:
                        OK('Writing: ', filename_out)
                else:
                    if verbose:
                        WARN('Skipping:', file_path)
            except Exception as error:
                FAIL(f"Error processing {file_path}: {error}")


def fast_input_file_standardize_cli():
    """ Command line interface (CLI) for function above"""
    parser = argparse.ArgumentParser(description="Recursively standardize FAST input files.")
    parser.add_argument('paths', nargs='+', help='Files or directories to process.')
    parser.add_argument('--target_extensions', type=str, default='.fst,.dat', help='Comma-separated target extensions.')
    parser.add_argument('--keys_for_exclusions', type=str, default='NumCoords', help='Comma-separated exclusion keys.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite source files.')
    parser.add_argument('--output_base_dir', type=str, default=None, help='Base directory for output files.')
    
    args = parser.parse_args()
    fast_input_file_standardize(
        paths=args.paths,
        target_extensions=args.target_extensions,
        keys_for_exclusions=args.keys_for_exclusions,
        overwrite=args.overwrite,
        output_base_dir=args.output_base_dir
    )

