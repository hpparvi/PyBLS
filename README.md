PyBLS
=====

Python wrapper for a multithreaded F2003 version of the EEBLS routine by Kovacs et al. (Kovacs, Zucker & Mazeh 2002, A&A, Vol. 391, 369). The original F77 routine can be found from http://www.konkoly.hu/staff/kovacs/index.html.


Install
-------
With a fairly modern gcc

	setup.py config_fc --fcompiler=gnu95 --opt="-Ofast -ffast-math" --f90flags="-cpp -fopenmp -march=native" install --user

With a bit less modern gcc

	setup.py config_fc --fcompiler=gnu95 --opt="-O3 -ffast-math" --f90flags="-cpp -fopenmp -march=native" install --user

Usage
-----

Import the BLS class

	from pybls import BLS

Create an instance

	bls = BLS(time, flux, error, period_range=(300,400), nf=4500, nbin=2000, qmin=0.001, qmax=0.01)

or

	bls = BLS(time, flux, error, fmin=0.00222, nf=4500, df=2e-7, nbin=2000, qmin=0.001, qmax=0.01)
	
run the search	
	
	res = bls()

The abovementioned parameters search for periods between 320 and 450 days, and seem to work well for super-Earths.
