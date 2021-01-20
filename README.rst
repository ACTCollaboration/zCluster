
zCluster is a code for measuring galaxy cluster photometric redshifts. For details of the algorithm, its
performance, and the output of the code, refer to `Hilton et al. (2018) <http://adsabs.harvard.edu/abs/2017arXiv170905600H>`_.
It has built-in support for querying large photometric surveys - currently:


* SDSS (DR7 - DR12)
* SDSS Stripe 82 (from SDSS DR7)
* CFHTLenS
* PS1 (DR2)
* DECaLS (DR7)
* DES (DR1 and Y3)
* KiDS (DR4)

The `published paper <http://adsabs.harvard.edu/abs/2017arXiv170905600H>`_ shows results based on SDSS, S82, 
and CFHTLenS only. The other surveys listed above are work in progress (so use with caution; PS1 in 
particular is problematic).

zCluster can also run on user-supplied .fits table photometric catalogs, provided that they have columns
named ``ID``\ , ``RADeg``\ , ``decDeg``\ , and magnitude column names in the form ``u_MAG_AUTO``\ , ``u_MAGERR_AUTO`` etc..

Software needed
===============

zCluster itself is written in pure Python (developed on Python 3.6; it may still run using Python 2.7 but 
this is not supported). It requires the following additional Python modules (current versions used by the 
author are given in brackets, earlier and later versions also probably work):


* numpy (1.13.3)
* scipy (0.19.1)
* matplotlib (2.0.2)
* astLib (0.10.2)
* astropy (3.0.5)
* IPython (7.2.0)
* mastcasjobs (for PS1; https://github.com/rlwastro/mastcasjobs)
* casjobs (for PS1; https://github.com/dfm/casjobs)

IPython isn't really required, but is used for debugging.

To run on DES photometry, there is an additional dependency:


* easyaccess (1.4.5)

If you want to run the code in parallel, you will also need:


* mpi4py (2.0.0)

Note that if you want to run the code on a cluster, the bottleneck will be fetching the photometric catalogs
over the internet. The MPI mode is still useful though on any machine with multiple cores.

All of the dependencies can be installed using ``pip``.

Installation
============

As root:

.. code-block::

   sudo python setup.py install

Or, in your home directory:

.. code-block::

   python setup.py install --prefix=$HOME/local

Then add ``$HOME/local/bin`` to $PATH, and e.g., ``$HOME/local/lib/python3.6/site-packages`` to $PYTHONPATH.

.. code-block::

   export PATH=$HOME/local/bin:$PATH    
   export PYTHONPATH=$HOME/local/lib/python3.6/site-packages:$PYTHONPATH

Alternatively, 

.. code-block::

   python setup.py install --user

will install ``zCluster`` under ``$HOME/.local`` (on Ubuntu), and in some other default location on Mac.

Running zCluster
================

See the ``examples/`` dir. Here you will find two cluster catalog files (in .csv format rather than .fits),
``400SDAll.csv`` and ``400SDSmall.csv``. The latter contains just the first 20 rows of the former, and is enough
to check that the code is working.

This example runs zCluster using SDSS DR12 photometry on part of the `Burenin et al. (2007) <http://adsabs.harvard.edu/abs/2007ApJS..172..561B>`_ 400SD cluster catalog.

``zCluster 400SDSmall.csv SDSSDR12``

This example will result in the creation of a catalog file called ``zCluster_400SDSmall_SDSSDR12.fits``\ ,
and a directory (where the results for individual clusters are cached) called ``400SDSmall_SDSSDR12/``.

If you want to check how well zCluster is doing versus spectroscopic redshifts, you can use the
zClusterComparisonPlot code.

``zClusterComparisonPlot 400SDSmall.csv zCluster_400SDSmall_SDSSDR12.fits``

This will produce a plot named ``comparison_400SDSmall_vs_zCluster_400SDSmall_SDSSDR12.pdf``\ , as well
outputting relevant statistics to the console.

You can find out about other options for both codes using the ``-h`` flag.

Comments, bug reports, help, suggestions etc..
==============================================

Please contact Matt Hilton matt.hilton@mykolab.com.
