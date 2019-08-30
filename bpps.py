import os, re, sys
import subprocess as sp
import random, string
import numpy as np
from package_toolkit.utils import *
import package_toolkit.pfunc as pfunc
import package_toolkit.settings as settings

DEBUG=False

def bpps(sequence, package='vienna', constraint=None, T=37, coaxial=False):
    '''Compute base pairing probability matrix for sequence.

    Args:
    sequence (str): nucleic acid sequence
    T (float): temperature
    constraint (str): structure constraint
    motif (str): argument to vienna motif
    coaxial (for vfold): True or false

  Returns
    array: NxN matrix of base pair probabilities
  '''
    try:
        pkg, version = package.lower().split('_')
    except:
        pkg, version = package, None

    if 'nupack' in package:
        return bpps_nupack_(sequence, version = version, T = T)
    elif 'vfold' in package:
    	return bpps_vfold_(sequence, version = version, T = T, coaxial = coaxial)
    else:
        _, tmp_file = pfunc.pfunc(sequence, package=package, bpps=True, constraint=constraint, T=T)
        if 'contrafold' in package:
            return bpps_contrafold_(sequence, tmp_file)
        elif 'vienna' in package:
            return bpps_vienna_(sequence, tmp_file)
        elif 'rnasoft' in package:
            return bpps_rnasoft_(sequence, tmp_file)
        elif 'rnastructure' in package:
            return bpps_rnastructure_(sequence, tmp_file)
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

    DIR = settings.LOC['nupack']

    if dangles:
        dangle_option='some'
    else:
        dangle_option='none'

    seqfile = utils.write([sequence])

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

def bpps_rnastructure_(sequence, tmp_file):

    DIR = settings.LOC['rnastructure']

    pfsfile = tmp_file #'%s/rnastructtmp.pfs' % settings.LOC['TMP']
    outfile = '%s/rnastructtmp.probs' % settings.LOC['TMP']
    command = ['%s/ProbabilityPlot' % DIR, pfsfile, outfile, '-t']
    probs=np.zeros([len(sequence), len(sequence)])

    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = p.communicate()

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

    DIR = settings.LOC["vfold_0"]

    cwd = os.getcwd()
    os.chdir(DIR) #vfold precompiled binaries don't work being called from elsewhere

    if DEBUG: print(os.getcwd())

    seqfile = write([sequence])

    outfile = filename()+'.pij'

    command = ['./Vfold2d_npk_mac.o %d %d %s %s %d' % (int(coaxial), T, seqfile, outfile, int(version))]

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

    os.remove(seqfile)
    probs = np.zeros([len(sequence),len(sequence)])
    p_ij_output = np.loadtxt(outfile,usecols=(0,2,3)) #col 0: set of inds 1, col 1: set of inds 2, col 2: bpp

    for i,j,p in p_ij_output:
    	probs[int(i-1),int(j-1)] = p
    	probs[int(j-1),int(i-1)] = p
    os.remove(outfile)

    return probs
    #output: take second field of last line for Z 



