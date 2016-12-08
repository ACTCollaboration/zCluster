# zCluster

A code for measuring galaxy cluster photometric redshifts. For details of the algorithm, its
performance, and the output of the code, refer to the ACTPol D56 clusters paper.

zCluster has built-in support for querying large photometric surveys - currently:
    
* SDSS (DR7 -- DR12)
* SDSS Stripe 82
* CFHTLenS
* DECaLS (DR3 - experimental!)

zCluster can also run on user supplied .fits table photometric catalogs, provided that they have columns
named `ID`, `RADeg`, `decDeg`, and magnitude column names in the form `u_MAG_AUTO`, `u_MAGERR_AUTO` etc..

## Software needed

zCluster itself is written in pure python (2.7.x). It requires the following additional python modules 
(current versions used by the author are given in brackets, earlier and later versions also probably work):

* pyfits (3.3)
* numpy (1.11.1)
* scipy (0.17.1)
* matplotlib (1.5.2)
* astLib (git version probably: get it with `git clone http://git.code.sf.net/p/astlib/git astlib-git`)
* IPython (2.4.1)

In addition, for surveys other than SDSS, the Schlegel et al. `dust_getval` code and maps are needed, in order
to calculate and apply Galactic extinction corrections. The `dust_getval` needs to be in your $PATH, and
the appropriate environment variables set.

IPython isn't really required, but is used for debugging. Note that astropy could be used to replace some
of these dependencies in the future.

There is also an optional dependency, if you want to run the code in parallel:
    
* mpi4py (2.0.0)

Note that if you want to run the code on a cluster, the bottleneck will be fetching the photometric catalogs
over the internet. The MPI mode is still useful though on any machine with multiple cores.

## Installation

As root:
    
```
    sudo python setup.py install
```

Or, in your home directory:
    
```
   python setup.py install --prefix=$HOME/local
```

Then add `$HOME/local/bin` to $PATH, and e.g., `$HOME/local/lib/python2.5/site-packages` to $PYTHONPATH.

## Running zCluster

See the `examples/` dir. Here you will find two cluster catalog files (in .csv format rather than .fits),
`400SDAll.csv` and `400SDSmall.csv`. The latter contains just the first 20 rows of the former, and is enough
to check that the code is working.

This example runs zCluster using SDSS DR12 photometry on part of the [Burenin et al. (2007)](http://adsabs.harvard.edu/abs/2007ApJS..172..561B) 400SD cluster catalog.

`zCluster 400SDSmall.csv SDSSDR12`

This example will result in the creation of a catalog file called `zCluster_400SDSmall_SDSSDR12.fits`,
and a directory (where the results for individual clusters are cached) called `400SDSmall_SDSSDR12/`.

If you want to check how well zCluster is doing versus spectroscopic redshifts, you can use the
zClusterComparisonPlot code.

`zClusterComparisonPlot 400SDSmall.csv zCluster_400SDSmall_SDSSDR12.fits`

This will produce a plot named `comparison_400SDSmall_vs_zCluster_400SDSmall_SDSSDR12.pdf`, as well
outputting relevant statistics to the console.

You can find out about other options for both codes using the `-h` flag.

## Comments, bug reports, help, suggestions etc..

Please contact Matt Hilton <matt.hilton@mykolab.com>.

