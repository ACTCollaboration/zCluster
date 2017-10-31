# zCluster

A code for measuring galaxy cluster photometric redshifts. For details of the algorithm, its
performance, and the output of the code, refer to [Hilton et al. 2017](http://adsabs.harvard.edu/abs/2017arXiv170905600H).

zCluster has built-in support for querying large photometric surveys - currently:
    
* SDSS (DR7 - DR12)
* SDSS Stripe 82 (from SDSS DR7)
* CFHTLenS
* PS1 (DR1; experimental!)
* DECaLS (DR3; experimental!)

The code has only been tested on SDSS, S82, and CFHTLenS 
([see the paper](http://adsabs.harvard.edu/abs/2017arXiv170905600H)); PS1 and DECaLS are work in progress.

zCluster can also run on user supplied .fits table photometric catalogs, provided that they have columns
named `ID`, `RADeg`, `decDeg`, and magnitude column names in the form `u_MAG_AUTO`, `u_MAGERR_AUTO` etc..

## Software needed

zCluster itself is written in pure python (2.7.x). It requires the following additional python modules 
(current versions used by the author are given in brackets, earlier and later versions also probably work):

* pyfits (3.4)
* numpy (1.13.1)
* scipy (0.17.0)
* matplotlib (2.0.2)
* astLib (0.9.2+ or git version: get it with `git clone http://git.code.sf.net/p/astlib/git astlib-git`)
* astropy (1.1.1)
* IPython (2.4.1)

IPython isn't really required, but is used for debugging. Note that astropy could be used to replace some
of the other dependencies in future (e.g., pyfits).

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

Then add `$HOME/local/bin` to $PATH, and e.g., `$HOME/local/lib/python2.7/site-packages` to $PYTHONPATH.

```
export PATH=$HOME/local/bin:$PATH    
export PYTHONPATH=$HOME/local/lib/python2.7/site-packages:$PYTHONPATH
```

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

