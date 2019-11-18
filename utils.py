import os, re
import subprocess as sp
import random, string
import numpy as np
import arnie

def write_vector_to_file(vector, outfile):
  for x in vector:
    outfile.write('%.3f\n' % x)
  return

def write_matrix_to_file(matrix, outfile):
  for x in matrix:
    outfile.write('\t'.join(['%.3f' % y for y in x])+'\n')
  return

def convert_dotbracket_to_bp_list(s):
    m = {}
    bp1=[]
    bp2=[]
    for i, char in enumerate(s):
        if char=='(':
            bp1.append(i)
        if char==')':
            bp2.append(i)
    for i in list(reversed(bp1)):
        for j in bp2:
            if j > i:
                m[i]=j
                m[j]=i

                bp2.remove(j)
                break
    return m

def convert_dbn_to_RNAstructure_input(seq, constraints, filename):
  assert(len(seq) == len(constraints))

  bp_list = convert_dotbracket_to_bp_list(constraints)

  SS_list, pairs_list = [], []


  for i, (s, c) in enumerate(list(zip(seq, constraints))):
    if c=='x':
      SS_list.append(i+1)
    elif c=='.':
      pass
    elif c =='(':
      pairs_list.append([i+1, bp_list[i]+1])
    elif c ==')':
      pass
    else:
      print('Error reading constraint string', c)

  with open('%s' % filename, 'w') as out:
    out.write('SS:\n')
    for x in SS_list:
      out.write('%d\n' % x)
    out.write('-1\nPairs:\n')
    for x,y in pairs_list:
      out.write('%d %d\n' % (x,y))
    out.write('-1 -1')

def write_constraint_string(seq, constraint_dbn):
  '''write set of integers to represent constraints, i.e. for use in bpseq format.'''

  assert(len(seq) == len(constraint_dbn))

  bp_list = convert_dotbracket_to_bp_list(constraint_dbn)

  constraint_list = []

  for i, c in enumerate(constraint_dbn):
        if c=='x':
          constraint=0
        elif c=='.':
          constraint=-1 #or -1 if undefined
        elif c in ['(',')']:
          constraint=bp_list[i]+1
        else:
          print('Error reading constraint string', c)
        constraint_list.append(constraint)
  return constraint_list

def convert_dbn_to_contrafold_input(seq, constraints, filename):
  constraint_list = write_constraint_string(seq, constraints)
  with open('%s' % filename, 'w') as out:
    for i in range(len(seq)):
      out.write('%d\t%s\t%d\n'%(i+1, seq[i], constraint_list[i]))

def convert_multiple_dbns_to_eternafold_input(seq, list_of_constraint_strings, filename):
  '''hard-coded to have 3 constraints right now for use in eternafold training with kd-ligand data.'''
  constraint_lists=[]
  for constraint_string in list_of_constraint_strings:
    constraint_lists.append(write_constraint_string(constraint_string))
    
  with open('%s' % filename, 'w') as out:
    for i in range(len(seq)):
      out.write('%d\t%s\t%d\t%d\t%d\n' % (i+1, seq[i], constraint_list[0][i], constraint_list[1][i], constraint_list[2][i]))

def write_constraints(seq, MS2=False, LIG=False, lig1=('nAGGAUAU','(xxxxxx('), lig2=('AGAAGGn',')xxxxx)')):
  '''Inputs:
  seq: RNA sequence
  MS2: bool, whether to include MS2 constraint or not
  lig1: tuple (seq, struct) for 5' portion of ligand aptamer. Default is FMN.
  lig2: tuple (seq, struct) for 3' portion of ligand aptamer

  Outputs:
  dbn string, () for paired, x for unpaired, . for unspecified
  '''

  #when FMN aptamer and MS2 aptamer overlap, MS2 loses out on bp
  MS2_apt='ACAUGAGGAUCACCCAUGU'
  LIG_apt1=lig1[0].replace('n','')
  LIG_apt2=lig2[0].replace('n','')

  unpaired_list=[]
  bp_list={}

  dbn_string=['.']*len(seq)

  if LIG:
      if seq.find(LIG_apt1) == -1:
        raise RuntimeError("ligand 5' aptamer domain not found")

      else:
        if lig1[0].startswith('n'):
          start1 = seq.find(LIG_apt1) + len(LIG_apt1) - len(lig1[0])
          if start1 < 1:
            start1 = seq.find(LIG_apt1,start1+len(lig1[0])+1) + len(LIG_apt1) - len(lig1[0])
        else:
          start1 = seq.find(LIG_apt1)
          if start1 < 1:
            start1 = seq.find(LIG_apt1,start1+len(lig1[0])+1)

        finish1 = start1 + len(lig1[0])   

        if lig2[0].startswith('n'):       
          start2 = seq.find(LIG_apt2, finish1+1) + len(LIG_apt2) - len(lig2[0])
        else:
          start2 = seq.find(LIG_apt2, finish1+1)
        finish2 = start2 + len(lig2[0])         
        #print('start1, finish1, start2, finish2 FMN', start1, finish1, start2, finish2)
        dbn_string[start1:finish1] = list(lig1[1])
        dbn_string[start2:finish2] = list(lig2[1])

  if MS2:
      if seq.find(MS2_apt) == -1:
        raise RuntimeError("MS2 aptamer domain not found")
      else:
        start=seq.find(MS2_apt)
        finish=start+len(MS2_apt)
        #print('start, finish MS2', start, finish)

        if dbn_string[start] != ".":
          #print('warning, aptamer overlap')
          dbn_string[start+1:finish-1]=list('((((x((xxxx))))))')
        else:
          dbn_string[start:finish]=list('(((((x((xxxx)))))))')


  return ''.join(dbn_string)

def get_missing_motif_bases(seq):
  FMN_apt1='AGGAUAU'
  FMN_apt2='AGAAGG'
  a = seq.find(FMN_apt1) - 1
  #print(a)
  b = seq.find(FMN_apt2) + len(FMN_apt2)
  #print(seq.find(FMN_apt2), b)

  return seq[a], seq[b]

def filename(n=6):
  """generate random filename

  Args:
    n (int): number of characters
  """
  rand = ''.join([random.choice(string.ascii_lowercase) for _ in range(n)])
  tmpdir = load_package_locations()['TMP']
  return '%s/%s' % (tmpdir, rand)

def write(lines, fname=None):
  """write lines to file

  Args:
    lines (list): line(s) to write to file
    fname (str): filename to write to
  """
  if fname is None:
    fname = '%s.in' % filename()
  with open(fname, 'w') as f:
    for line in lines:
      f.write('%s\n' % line)
  return fname

def print_available_packages():
  package_dct = load_package_locations()
  for key,v in package_dct.items():
    if key != "TMP":
      print(key,v)

def package_list():
  pkg_list=[]
  package_dct = load_package_locations()
  for key,v in package_dct.items():
    if key != "TMP":
      if v != "None":
        pkg_list.append(key)
  return pkg_list

def load_package_locations():
    '''Read in user-supplied file to specify paths to RNA folding packages. Specify this in your ~/.bashrc as $ARNIEFILE'''
    return_dct={}
    package_path = os.path.dirname(arnie.__file__)
    with open("%s" % os.environ["ARNIEFILE"],'r') as f:
        for line in f.readlines():
            if line.strip():
              if not line.startswith('#'):
                key, string = line.split(':')
                string = string.strip()
                return_dct[key] = string
    return return_dct
