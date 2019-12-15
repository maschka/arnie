import os, re, sys
import subprocess as sp
import random, string
import numpy as np
from .utils import *
from .pfunc import pfunc

DEBUG=False

# load package locations from yaml file, watch! global dict
package_locs = load_package_locations()

def bpps(sequence, package='vienna', constraint=None, T=37, coaxial=True, dangles=True,param_file=None,reweight=None):
    ''' Compute base pairing probability matrix for RNA sequence.

    Args:
    sequence (str): nucleic acid sequence
    T (float): temperature (Celsius)
    constraint (str): structure constraint (functional in vienna, contrafold, rnastructure)
    motif (str): argument to vienna motif
    dangles (bool): dangles or not, specifiable for vienna, nupack
    coaxial (bool): coaxial stacking or not, specifiable for rnastructure, vfold
    noncanonical(bool): include noncanonical pairs or not (for contrafold, RNAstructure (Cyclefold))

    Possible packages: 'vienna_2', 'vienna_1','contrafold_1','contrafold_2','nupack_95','nupack_99','rnasoft_2007','rnasoft_1999','rnastructure','vfold_0','vfold_1'
 
    Returns
    array: NxN matrix of base pair probabilities
  '''
    try:
        pkg, version = package.lower().split('_')
    except:
        pkg, version = package, None

    if not dangles and pkg not in ['vienna','nupack']:
        print('Warning: %s does not support dangles options' % pkg)
    if not coaxial and pkg not in ['rnastructure','vfold']:
        print('Warning: %s does not support coaxial options' % pkg)

    if pkg=='nupack':
        return bpps_nupack_(sequence, version = version, dangles = dangles, T = T)
    elif pkg=='vfold':
    	return bpps_vfold_(sequence, version = version, T = T, coaxial = coaxial)
    else:
        _, tmp_file = pfunc(sequence, package=package, bpps=True, constraint=constraint, T=T, coaxial=coaxial, dangles=dangles, param_file=param_file,reweight=reweight)
        if 'contrafold' in package:
            return bpps_contrafold_(sequence, tmp_file)
        elif 'vienna' in package:
            return bpps_vienna_(sequence, tmp_file)
        elif 'rnasoft' in package:
            return bpps_rnasoft_(sequence, tmp_file)
        elif 'rnastructure' in package:
            return bpps_rnastructure_(sequence, tmp_file, coaxial=coaxial)
        else:
            raise RuntimeError('package not yet implemented')

def bpps_vienna_(sequence, tmp_file):

    dot_fname = tmp_file

    probs=np.zeros([len(sequence), len(sequence)])
    with open(dot_fname,'r') as f:
        for line in f.readlines():
            if 'ubox' in line:
                try:
                    i, j, p, _ = line.split()
                    i, j, p = int(i)-1, int(j)-1, float(p)**2
                    probs[i,j] = p
                    probs[j,i] = p
                except:
                    pass
    os.remove(dot_fname)
    return probs

def bpps_contrafold_(sequence, tmp_file):

    fname = tmp_file

    probs=np.zeros([len(sequence), len(sequence)])

    for line in open(fname).readlines():
        if len(line.split(':')) > 1:
            first_ind = int(line.split()[0])-1
            for x in line.split()[2:]:
                second_ind = int(x.split(':')[0])-1
                p = float(x.split(':')[1])
                probs[first_ind, second_ind] = p
                probs[second_ind, first_ind] = p

    os.remove(fname)

    return probs

def bpps_rnasoft_(sequence, tmp_file):
    fname = tmp_file

    probs=np.zeros([len(sequence), len(sequence)])
    for line in open(fname).readlines():
        i,j,p = int(line.split()[0]), int(line.split()[1]), float(line.split()[2])
        probs[i,j] = p
        probs[j,i] = p

    os.remove(fname)

    return probs

def bpps_nupack_(sequence, version='95', T=37, dangles=True):

    if not version: version='95'
    
    nupack_materials={'95': 'rna1995', '99': 'rna1999'}

    DIR = package_locs['nupack']

    if dangles:
        dangle_option='some'
    else:
        dangle_option='none'

    seqfile = write([sequence])

    command=['%s/pairs' % DIR, '%s' % seqfile.replace('.in',''),
      '-T', str(T), '-material', nupack_materials[version], '-dangles', dangle_option]

    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if p.returncode:
        raise Exception('Nupack pfunc failed: on %s\n%s' % (sequence, stderr))

    ppairs_file = '%s.ppairs' % seqfile.replace('.in','')
    os.remove(seqfile)

    probs=np.zeros([len(sequence), len(sequence)])

    with open(ppairs_file, 'r') as f:
        for line in f.readlines():
            if not line.startswith('%'):
                fields = line.split()
                if len(fields) > 1:
                    if int(fields[1]) <= len(sequence):
                        i, j, p = int(fields[0])-1, int(fields[1])-1, float(fields[2])
                        probs[i,j] = p
                        probs[j,i] = p


    return probs

def bpps_rnastructure_(sequence, tmp_file, coaxial=True):

    DIR = package_locs['rnastructure']

    pfsfile = tmp_file #'%s/rnastructtmp.pfs' % package_locs['TMP']
    outfile = '%s/rnastructtmp.probs' % package_locs['TMP']
    command = ['%s/ProbabilityPlot' % DIR, pfsfile, outfile, '-t']

    probs=np.zeros([len(sequence), len(sequence)])

    if DEBUG: print(' '.join(command))
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

    if DEBUG:
        print('stdout')
        print(stdout)
        print('stderr')
        print(stderr)

    if p.returncode:
        raise Exception('RNAstructure ProbabilityPlot failed: on %s\n%s' % (seq, stderr))

    with open(outfile, 'r') as f:
        for line in f.readlines()[2:]:
            fields = line.split()
            i, j, p = int(fields[0])-1, int(fields[1])-1, 10**(-1*float(fields[2]))
            probs[i,j] = p
            probs[j,i] = p
    
    os.remove(outfile)
    os.remove(pfsfile)
    return probs

def bpps_vfold_(sequence, version='0',T=37, coaxial=True):
    #available versions: 0 for Turner 04 params, 1 for Mfold 2.3 params

    DIR = package_locs["vfold"]

    cwd = os.getcwd()
    os.chdir(DIR) #vfold precompiled binaries don't work being called from elsewhere

    if DEBUG: print(os.getcwd())

    seqfile = write([sequence])

    outfile = filename()+'.pij'

    if sys.platform=="linux":
        platform='linux'
    elif sys.platform=="darwin":
        platform='mac'
    elif sys.platform=="win32":
        platform='win'
    else:
        raise RuntimeError('Vfold has binaries for linux, macOS, and win')

    command = ['./Vfold2d_npk_%s.o %d %d %s %s %d' % (platform, int(coaxial), T, seqfile, outfile, int(version))]

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
        raise Exception('Vfold2d_npk failed: on %s\n%s' % (sequence, stderr))

    os.remove(seqfile)
    probs = np.zeros([len(sequence),len(sequence)])
    p_ij_output = np.loadtxt(outfile,usecols=(0,2,3)) #col 0: set of inds 1, col 1: set of inds 2, col 2: bpp

    for i,j,p in p_ij_output:
    	probs[int(i-1),int(j-1)] = p
    	probs[int(j-1),int(i-1)] = p
    os.remove(outfile)

    return probs
    #output: take second field of last line for Z 



