import bpps, settings
import sys
sample_seq = 'GGGGAAAACCCC'

def test_bpps(pkg):
	p = bpps.bpps(sample_seq, package = pkg)
	print('test bpps %s' % pkg)
	print(p[0])
	return

if __name__=='__main__':
	print("Test: printing first row of bpp matrices")
	for pkg in sorted(settings.LOC.keys()):
		if pkg=='TMP':
			continue

		test_bpps(pkg.lower())
