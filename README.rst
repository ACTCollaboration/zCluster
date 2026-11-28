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

=====================================================  =====================================================
Survey                                                 ``database`` name(s)
=====================================================  =====================================================
SDSS                                                   ``SDSSDR12``, ``SDSSDR16``
Pan-STARRS1 DR2 [#ps1]_                                ``PS1``
DESI Legacy Imaging Surveys (tractor catalogs, NERSC)  ``DECaLSDR8``, ``DECaLSDR9``, ``DECaLSDR10``
DESI Legacy Imaging Surveys (NOIRLab Data Lab)         ``DL_DECaLSDR10``, ``DL_DECaLSDR11``
DELVE DR2 (NOIRLab Data Lab)                           ``DELVEDR2``
DES public releases [#des]_                            ``DESDR1``, ``DESDR2``
KiDS DR4 (including VIKING)                            ``KiDSDR4``
VST ATLAS DR4                                          ``ATLASDR4``
Rubin Observatory LSST data previews [#rubin]_         ``RubinDP0``, ``RubinDP1``, ``RubinDP2``
=====================================================  =====================================================

.. [#ps1] Requires a `MAST CasJobs <https://mastweb.stsci.edu/ps1casjobs/>`_ login, given in the
   ``CASJOBS_WSID`` and ``CASJOBS_PW`` environment variables.
.. [#des] Requires `easyaccess <https://github.com/des-labs/easyaccess>`_ and a DES database login.
.. [#rubin] Requires Rubin data rights and a
   `Rubin Science Platform access token <https://rsp.lsst.io/guides/auth/creating-user-tokens.html>`_,
   saved in ``~/.rsp-tap.token``.

Run ``zCluster -h`` for the full list of options.
 
For details of the algorithm, its performance, and the output of the code, refer to 
`Hilton et al. (2018) <https://scixplorer.org/abs/2018ApJS..235...20H/abstract>`_, which presents
results based on SDSS, S82, and CFHTLenS, and/or 
`Hilton et al. (2021) <https://scixplorer.org/abs/2021ApJS..253....3H/abstract>`_, which presents
results based on DECaLS DR8. **zCluster** was also used in
`Aguena et al. (2026) <https://scixplorer.org/abs/2026OJAp....955863A/abstract>`_ with DECaLS DR10 data.
The other surveys listed above are work in progress (so use with caution; PS1 in
particular is problematic). `Pillay et al. (2021) <https://scixplorer.org/abs/2021arXiv211104340P/abstract>`_
presents the first use of the package for producing projected galaxy density maps.

If you find **zCluster** useful in your work, please cite whichever one
of the above papers that you think is appropriate (together, of course, with the appropriate papers
for the optical/IR survey used).

**zCluster** can also run on user-supplied .fits table photometric catalogs, provided that they have columns
named ``ID``\ , ``RADeg``\ , ``decDeg``\ , and magnitude column names in the form ``u_MAG_AUTO``\ ,
``u_MAGERR_AUTO`` etc..

**zCluster** is under active development, and not all documentation is up to date. The package also
contains some experimental features that are not necessarily well tested.
