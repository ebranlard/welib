import os
import subprocess
import multiprocessing


def run_cmds(inputfiles, exe, parallel=True, ShowOutputs=True, nCores=None):
    """ Run a set of simple commands of the form `exe input_file`
        By default, the commands are run in parallel
    """
    ps=[]
    iProcess=0
    if nCores is None:
        nCores=multiprocessing.cpu_count()
    if nCores<0:
        nCores=len(inputfiles)+1
    for i,f in enumerate(inputfiles):
        #print('Process {}/{}: {}'.format(i+1,len(inputfiles),f))
        ps.append(run_cmd(f, exe, wait=(not parallel), ShowOutputs=ShowOutputs))
        iProcess += 1
        # waiting once we've filled the number of cores
        if parallel:
            if iProcess==nCores:
                for p in ps:
                    p.wait()
                ps=[]
                iProcess=0
    for p in ps:
        p.wait()
#

def run_cmd(input_file, exe, wait=True, ShowOutputs=False, debian=None):
    """ Run a simple command of the form `exe input_file`  """
    if not os.path.isabs(input_file):
        input_file=os.path.abspath(input_file)
    if not os.path.exists(exe):
        raise Exception('FAST Executable not found: {}'.format(exe))
    args= [exe,input_file]
    #args = 'cd '+workdir+' && '+ exe +' '+basename
    shell=False
    if debian is not None:
        input_file=input_file.replace('C:','/mnt/c')
        input_file=input_file.replace('Work','work')
        input_file=input_file.replace('\\','/')
        #print(input_file)
        workdir  = os.path.dirname(input_file)
        basename = os.path.basename(input_file)
        exe = './'+exe.strip()
        args = ['debian', 'run', 'cd '+workdir+' && '+exe+' '+basename]
        shell=True
    if ShowOutputs:
        STDOut= None
    else:
        STDOut= open(os.devnull, 'w') 
    #print(args)
    if wait:
        p=subprocess.call(args , stdout=STDOut, stderr=subprocess.STDOUT, shell=shell)
    else:
        p=subprocess.Popen(args, stdout=STDOut, stderr=subprocess.STDOUT, shell=shell)
    return p
