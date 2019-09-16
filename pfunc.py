import os, re, sys, shutil
import subprocess as sp
import random, string
import numpy as np
from .utils import *

DEBUG=False

# load package locations from yaml file, watch! global dict
package_locs = load_package_locations_from_yaml('user_default.yaml')

def pfunc(seq, package='vienna_2', T=37,
    constraint=None, motif=None,
    dangles=True, noncanonical=False,
    bpps=False, mfe=False, mea=False, param_file=None, coaxial=True):
    ''' Compute partition function for RNA sequence.

        Args:
        seq (str): nucleic acid sequence
        T (float): temperature (Celsius)
        constraint (str): structure constraints
        motif (str): argument to vienna motif 
        dangles (bool): dangles or not, specifiable for vienna, nupack
        coaxial (bool): coaxial stacking or not, specifiable for rnastructure, vfold
        noncanonical(bool): include noncanonical pairs or not (for contrafold, RNAstructure (Cyclefold))

        Possible packages: 
        'vienna_2', 'vienna_1','contrafold_1','contrafold_2','nupack_95','nupack_99','rnasoft_2007','rnasoft_1999','rnastructure','vfold_0','vfold_1'
        
    Returns
        str, float: secondary structure representation and free energy
    '''

    try:
        pkg, version = package.lower().split('_')
    except:
        pkg, version = package, None

    if not bpps: # if bpps, already printed these warnings
        if not dangles and pkg not in ['vienna','nupack']:
            print('Warning: %s does not support dangles options' % pkg)
        if not coaxial and pkg not in ['rnastructure','vfold']:
            print('Warning: %s does not support coaxial options' % pkg)

    if pkg=='vienna':
        Z, tmp_file = pfunc_vienna_(seq, version=version, T=T, dangles=dangles, constraint=constraint, motif=motif, bpps=bpps)

    elif pkg=='contrafold':
        Z, tmp_file = pfunc_contrafold_(seq, version=version, T=T, constraint=constraint, bpps=bpps, param_file=param_file)

    elif pkg=='rnastructure':
        Z, tmp_file = pfunc_rnastructure_(seq, version=version, T=T, coaxial=coaxial, constraint=constraint, bpps=bpps)

    elif pkg=='rnasoft':
        Z, tmp_file = pfunc_rnasoft_(seq, version=version, T=T, constraint=constraint, bpps=bpps)

    elif pkg=='nupack':
        Z, tmp_file = pfunc_nupack_(seq, version=version, dangles=dangles, T=T)

    elif pkg=='vfold':
        Z, tmp_file = pfunc_vfold_(seq, version=version, T=T, coaxial=coaxial)

    else:
        raise ValueError('package %s not understood.' % package)

    if bpps:
        return Z, tmp_file
    else:
        if tmp_file:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
        return Z

def pfunc_vienna_(seq, T=37, version='2', constraint=None, motif=None,
                                    dangles=True, bpps=False):
    """get partition function structure representation and Z

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif  
    Returns
        str, float: secondary structure representation and Z
    """

    if not version:
        version='2'

    if version.startswith('2'):
        LOC=package_locs['vienna_2']
    elif version.startswith('1'):
        LOC=package_locs['vienna_1']
    else:
        raise RuntimeError('Error, vienna version %s not present' % version)

    command = ['%s/RNAfold' % LOC, '-T', str(T), '-p']
    if motif is not None:
        command.append('--motif="%s"' % motif)

    if constraint is not None:
        fname = write([seq, constraint])
        command.append('-C')
    else:
        fname = write([seq])

    if not dangles:
        command.append('--dangles=0')

    with open(fname) as f:
        if DEBUG: print(fname)
        if DEBUG: print(' '.join(command))
        p = sp.Popen(command, stdin=f, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    if p.returncode:
        raise Exception('RNAfold failed: on %s\n%s' % (seq, stderr))
    os.remove(fname)
    m = re.search('([,|\(\.\)\]\[\{\}]+)\s+\[\s*(-*[0-9]+\.[0-9]+)', stdout.decode('utf-8'))

    free_energy = float(m.group(2))

    tmp_file= '%s.ps' % filename()
    shutil.move('dot.ps', tmp_file)

    return np.exp(-1*free_energy/(.0019*(273+T))), tmp_file

def pfunc_contrafold_(seq, T=37, version='se', constraint=None, bpps=False,param_file=None):
    """get partition function structure representation and free energy

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif  
    Returns
        float: partition function
    """
    if not version: version='se'


    fname = '%s.in' % filename()

    if version.startswith('2'):
        LOC=package_locs['contrafold_2']
    elif version.startswith('1'):
        LOC=package_locs['contrafold_1']
    elif version.startswith('se'):
        LOC=package_locs['contrafold_se']
    else:
        raise RuntimeError('Error, Contrafold version %s not present' % version)

    command = ['%s/contrafold' % LOC, 'predict', fname]

    if bpps:
        posterior_fname = '%s.posteriors' % filename()
        command = command + ['--posteriors', '0.001', posterior_fname]
    else:
        command.append('--partition')

    if param_file is not None:
        command = command + ['--params', param_file]


    if constraint is not None:
        convert_dbn_to_contrafold_input(seq, constraint, fname)
        command.append('--constraints')
    else:
        convert_dbn_to_contrafold_input(seq, ''.join(['.' for x in range(len(seq))]), fname)

    if DEBUG: print(' '.join(command))
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    if p.returncode:
        raise Exception('Contrafold failed: on %s\n%s' % (seq, stderr))

    os.remove(fname)

    if not bpps:
        logZ = float(stdout.decode('utf-8').rstrip().split()[-1])
        return np.exp(logZ), None
    else:
        return 0, posterior_fname

def pfunc_rnasoft_(seq, version='99', T=37, constraint=None, bpps=False, mfe=False):
    DIR = package_locs['rnasoft']

    if not version: version='99'

    if mfe: print('rnasoft mfe not implemented yet') #simfold instead of simfold pf

    # supported versions: 07, 99, 99-no-dangles, BL-no-dangles, BLstar, LAM-CG, NOM-CG
    param_locs = {'07': '%s/params/CG_best_parameters_ISMB2007.txt' % DIR,
    '99': '%s/params/turner_parameters_fm363_constrdangles.txt' % DIR,
    '99-no-dangles': '%s/params/turner_parameters_fm363_dangles0.txt' % DIR,
    'bl-no-dangles': '%s/params/BL-no-dangles.txt' % DIR,
    'blstar': '%s/params/BLstar.txt' % DIR,
    'lam-cg': '%s/params/LAM-CG.txt' % DIR,
    'nom-cg': '%s/params/NOM-CG.txt' % DIR}

    command = ['%s/simfold_pf' % DIR, '-s', seq, '-p', param_locs[version]]

    if DEBUG: print(' '.join(command))
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    bpps_fname = '%s.bpps' % filename()
    if bpps:
        with open(bpps_fname,'w') as f:
            for line in stdout.decode('utf-8').split('\n')[5:]:
                if not 'Glog' in line and len(line) > 1:
                    f.write(line+'\n')

    if p.returncode:
        raise Exception('RNAsoft partition failed: on %s\n%s' % (seq, stderr))

    return float(stdout.decode('utf-8').split('\n')[1].split()[-1]), bpps_fname

def pfunc_nupack_(seq, version='95', T=37, dangles=True):

    if not version: version='95'
    nupack_materials={'95': 'rna1995', '99': 'rna1999'}

    DIR = package_locs['nupack']

    if dangles:
        dangle_option='some'
    else:
        dangle_option='none'

    seqfile = write([seq])

    command=['%s/pfunc' % DIR, '%s' % seqfile.replace('.in',''),'-T', str(T), '-material', nupack_materials[version], '-dangles', dangle_option]

    if DEBUG: print(' '.join(command))
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    if p.returncode:
        raise Exception('Nupack pfunc failed: on %s\n%s' % (seq, stderr))

    Z=float(stdout.decode('utf-8').split('\n')[-2])

    os.remove(seqfile)

    return Z, None

def pfunc_rnastructure_(seq, version=None, T=37, constraint=None, coaxial=True,bpps=False):
    """get partition function structure representation and free energy

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif  
        coaxial (bool): Coaxial stacking or not (default True)
    Returns
        float: partition function
    """

    seqfile = write([seq])
    pfsfile = '%s.pfs' % filename()
    DIR = package_locs['rnastructure']
    command = ['%s/partition' % DIR, seqfile, pfsfile, '-T', str(T+273)]

    if not coaxial:
        command.extend(['--disablecoax'])

    if constraint is not None:
        fname = '%s.CON' % filename()
        #print(fname)
        convert_dbn_to_RNAstructure_input(seq, constraint, fname)
        command.extend(['--constraint', fname])

    if DEBUG: print(' '.join(command))
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    if p.returncode:
        raise Exception('RNAstructure partition failed: on %s\n%s' % (seq, stderr))

    os.remove(seqfile)

    if not bpps:
        command = ['%s/EnsembleEnergy' % DIR, pfsfile]

        if DEBUG: print(' '.join(command))
        p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

        stdout, stderr = p.communicate()

        if DEBUG:
            print('stdout')
            print(stdout)
            print('stderr')
            print(stderr)

        if p.returncode:
            raise Exception('RNAstructure EnsembleEnergy failed: on %s\n%s' % (seq, stderr))

        if DEBUG: print(stdout.decode('utf-8').split('\n')[3])
        free_energy = float(stdout.decode('utf-8').split('\n')[3].split(' ')[-2])

        return np.exp(-1*free_energy/(.0019*(273+T))), pfsfile
    else:
        return 0, pfsfile

def pfunc_vfold_(seq, version='0', T=37, coaxial=True, bpps=False):
    #available versions: 0 for Turner 04 params, 1 for Mfold 2.3 params
    #for bpps
        # command = ['%s/Vfold2d_npk_mac.o %d %d %s %s %d' % (DIR, int(coaxial),\
        #  T, infile, outfile, int(version))]

    DIR = package_locs["vfold"]

    cwd = os.getcwd()
    os.chdir(DIR) #vfold precompiled binaries don't work being called from elsewhere

    if DEBUG: print(os.getcwd())

    seqfile = write([seq])

    if sys.platform=="linux":
        platform='linux'
    elif sys.platform=="darwin":
        platform='mac'
    elif sys.platform=="win32":
        platform='win'
    else:
        raise RuntimeError('Vfold has binaries for linux, macOS, and win')


    command = ['./VfoldThermal_npk_%s.o %d %d %d %s tmp %d; cat tmp; rm tmp' % (platform, int(coaxial), T, T, seqfile, int(version))]

    if DEBUG: print(' '.join(command))

    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)

    stdout, stderr = p.communicate()
    os.chdir(cwd)

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)
    if p.returncode:
        raise Exception('VfoldThermal_npk failed: on %s\n%s' % (seq, stderr))

    Z=float(stdout.decode('utf-8').split('\n')[-2].split()[1])

    os.remove(seqfile)
    return Z, None
    #output: take second field of last line for Z 




