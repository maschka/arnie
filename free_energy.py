import os, re, sys
import subprocess as sp
import random, string
import numpy as np
from .utils import *
from .pfunc import pfunc

DEBUG=False

# load package locations from yaml file, watch! global dict
package_locs = load_package_locations()

def free_energy(seq, structure=None, package='vienna_2', T=37, dangles=True, reweight=None):
	''' Compute free energy of RNA sequence. If structure is given, computes free energy of that structure. 
			Otherwise, returns MFE structure of sequence [NOT IMPLEMENTED YET].

		Args:
		seq (str): nucleic acid sequence
		structure (str): structure in dot bracket notation
		T (float): temperature (Celsius)
		constraint (str): structure constraints
		motif (str): argument to vienna motif 
		dangles (bool): dangles or not, specifiable for vienna, nupack
		coaxial (bool): coaxial stacking or not, specifiable for rnastructure, vfold
		noncanonical(bool): include noncanonical pairs or not (for contrafold, RNAstructure (Cyclefold))

		Possible packages: 
		'vienna_1', 'vienna_2'
		
	Returns
		free energy (float)
	'''


	if package.startswith('vienna'):
		return pfunc(seq, package=package, T=T, dangles=dangles, constraint=structure, reweight=reweight, return_free_energy=True)
	elif package.startswith('contrafold'):
		Z_constrained = pfunc(seq, package=package, T=T, dangles=dangles, constraint=structure)
		if DEBUG: print("Z_constr: %.3f" % Z_constrained)
		if DEBUG: print("dG constr: %.3f" % -0.0019899*(273+T) * np.log(Z_constrained))
		return -0.0019899*(273+T) * np.log(Z_constrained) # .00198 is k in kcal/mol
	else:
		raise RuntimeError("%s `free_energy` not implemented yet" % package)
