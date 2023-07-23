.. image:: https://readthedocs.org/projects/zcluster/badge/?version=latest

**zCluster** is a package for measuring galaxy cluster photometric redshifts using
data from large public surveys. It can also produce photometric redshift estimates
and galaxy density maps for any point in the sky using the included `zField` tool.

* **Documentation:** https://zcluster.readthedocs.io
* **License:** `GPL v3 <COPYING>`_
* **Authors:** Matt Hilton, with contributions from Kabelo Kesebonye, Phumlani Phakathi,
  Denisha Pillay, and Damien Ragavan (not all reflected on GitHub).
* **Installation:** ``pip install zCluster``
* **Support:** Please use the `GitHub issues page <https://github.com/ACTCollaboration/zCluster/issues>`_, 
  and/or contact `Matt Hilton <mailto:matt.hilton@mykolab.com>`_.
  
**zCluster** has built-in support for querying large photometric surveys - currently:

* SDSS (DR7 - DR12)
* SDSS Stripe 82 (from SDSS DR7)
* CFHTLenS
* PS1 (DR2)
* DECaLS (DR8 - DR10)
* DES (DR1, DR2 and Y3 internal)
* KiDS (DR4)

For details of the algorithm, its performance, and the output of the code, refer to 
`Hilton et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018ApJS..235...20H/abstract>`_, which presents
results based on SDSS, S82, and CFHTLenS, and/or 
`Hilton et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJS..253....3H/abstract>`_, which presents
results based on DECaLS DR8. The other surveys listed above are work in progress (so use with caution; PS1 in 
particular is problematic). `Pillay et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021arXiv211104340P/abstract>`_
presents the first use of the package for producing projected galaxy density maps.

If you find **zCluster** useful in your work, please cite whichever one
of the above papers that you think is appropriate (together, of course, with the appropriate papers
for the optical/IR survey used).

**zCluster** can also run on user-supplied .fits table photometric catalogs, provided that they have columns
named ``ID``\ , ``RADeg``\ , ``decDeg``\ , and magnitude column names in the form ``u_MAG_AUTO``\ ,
``u_MAGERR_AUTO`` etc..

**zCluster** is under active development, and not all documentation is up to date. The package also
contains some experimental features that are not necessarily well tested.
