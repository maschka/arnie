# arnie
Python API to compute RNA energetics and do structure prediction from all available packages.

Das Lab, 2019

Hannah Wayment-Steele, with inspiration from MWu's `wuami_tools`


## Usage:
create a file that points to your builds of all the structure prediction packages you intend to make available.  An example file is provided in "user_default.yaml".  Create a variable in your .bashrc for this:

```
export ARNIEFILE="/path/to/arnie/<my_file.txt>"
```
NB: this file is technically yaml format, but isn't read in by yaml.

See `arnie_example.ipynb` for example syntax. In brief, comparing across packages is simple. For example, for computing base pairing probability matrices:

```
from arnie.bpps import bpps
%pylab inline
example_seq = 'GGGGAAAACCCC'

bpps={}

for pkg in ['vienna','contrafold','RNAsoft']:
    bpps[package] = bpps(example_seq, package=pkg)
    
imshow(bpps['vienna'])
```

Potentially helpful utilities in utils.py:

`write_constraints(seq, MS2=False, LIG=True)`: 
for riboswitches, write constraint string to then feed into contrafold or vienna to compute constrained partition functions.
I.E. for the three constrained states of a riboswitch, would do `MS2=False, LIG=True,`, `MS2=True, LIG=False`, `MS2=True, LIG=True`.

## Coming soon

Help for compiling packages
