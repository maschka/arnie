import pfunc
import settings

sample_seq = 'GGGGAAAACCCC'

def test_pkg(package):
	
	Z = pfunc.pfunc(sample_seq, package=package, bpps=False)
	print('test %s' % package, Z)
	return

def test_pkg_w_bpps(package):
	
	Z, tmp_file = pfunc.pfunc(sample_seq, package=package, bpps=True)
	print('test %s, tmp file for bpps %s' % (package, tmp_file), Z)
	return

if __name__=='__main__':
	for pkg in sorted(settings.LOC.keys()):
		if pkg=='TMP':
			continue
		test_pkg(pkg.lower())
		#test_pkg_w_bpps(pkg.lower())
