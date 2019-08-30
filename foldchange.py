import sys,os,subprocess, argparse
import numpy as np

from utils import *
from pfunc import *

def compute_AR(seq_file, z_file, ligand_conc=200e-6, kd=5e-7, ONswitch=True, package='vienna', motif_method=False, debug=True):

  seq = open(seq_file, 'r').read().rstrip()


  bonus = ligand_conc / kd

  MS2 = write_constraints(seq, MS2=True)
  FMN = write_constraints(seq, FMN=True)
  FMN_MS2 = write_constraints(seq, MS2=True, FMN=True)

  if debug:
    print(seq)
    print(FMN)
    print(MS2)
    print(FMN_MS2)

  _, Z = pfunc(seq, package = package)
  _, zMS2 = pfunc(seq, constraint = MS2, package = package)

  if motif_method:
    base_a, base_b = get_missing_motif_bases(seq)
    vienna_motif='%sAGGAUAU&AGAAGG%s,(xxxxxx(&)xxxxx),%.3f' % (base_a, base_b, -0.598*np.log(bonus))
    assert len(vienna_motif.split(',')[0]) == len(vienna_motif.split(',')[1])

    _, zFMN = pfunc(seq, motif = vienna_motif, package = package)
    _, zFMN_MS2 = pfunc(seq, constraint = MS2, motif = vienna_motif, package = package)

  else:
    _, zFMN = pfunc(seq, constraint = FMN, package = package)
    _, zFMN_MS2 = pfunc(seq, constraint = FMN_MS2, package = package)

  z_file.write("%.3f\t%.3f\t%.3f\t%.6f\n" % (Z, zFMN, zMS2, zFMN_MS2))
  if debug: print("%.3f\t%.3f\t%.3f\t%.6f\n" % (Z, zFMN, zMS2, zFMN_MS2))

  signal_off = zMS2 / Z

  if (package=='vienna' and motif_method):
    signal_on = zFMN_MS2 / zFMN
  else:
    signal_on = (zMS2 + zFMN_MS2*bonus) / (Z + zFMN*bonus)

  if ONswitch:
    return signal_on / signal_off

  else:
    return signal_off / signal_on

if __name__ == "__main__":

  '''Compute fold changes for R95, First command line arg is name of output file for partial partition functions.
  Second is number of sequences'''

  parser=argparse.ArgumentParser(
      description='''Compute foldchanges.''')
  parser.add_argument('-o', dest='o', help='output file name')
  parser.add_argument('-n', type=int, help='Number of constructs from R95')
  parser.add_argument('-p', dest='package', help='Package to use: vienna, contrafold, rnastruct, or rnasoft.')
  args=parser.parse_args()

  cd_hit_est_inds = np.loadtxt('sequences/R95_sequences/cd_hit_est_indices')
  sequence_files = ['sequences/R95_sequences/%d.seq' % x for x in cd_hit_est_inds[:args.n]]

  switch_inds = np.loadtxt('/scratch/users/hannahw1/das/data/MI_package/example_inputs/MS2_affinities/R95_onoff_ind')
  experimental_foldchanges = np.loadtxt('R95_foldchange')

  z_file = open('%s.z' % args.o,'w')

  f = open(args.o,'w')
  f.write('ind\texpt\t%s\tonoff\n' % args.package)

  for i, seq in enumerate(sequence_files):

    ind = int(cd_hit_est_inds[i])
    print(ind)
    if switch_inds[ind] == -1:
      ONswitch=False
    elif switch_inds[ind] == 1:
      ONswitch=True

    try:
      ARc = compute_AR(seq, z_file, ONswitch=ONswitch, package=args.package)
      f.write("%d\t%.2f\t%.2f\t%d\n" % (ind, experimental_foldchanges[ind], ARc, switch_inds[ind]))
      z_file.write('\n')
    except RuntimeError:
      pass

  z_file.close()
