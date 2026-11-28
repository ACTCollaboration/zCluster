.. _TestingPage:

=======
Testing
=======

**zCluster** uses the `Robot Framework <https://robotframework.org/>`_ for tests, which can be installed
with ``pip install robotframework``.

.. note::

   The files needed to run the tests are found in the
   `tests <https://github.com/ACTCollaboration/zCluster/tree/master/tests>`_ directory of the
   **zCluster** source distribution. The tests always run the code in the source tree that contains
   them (not whatever version of **zCluster** happens to be installed), so there is no need to
   re-install before testing a change.

There are two test suites:

* ``zCluster_regression.robot`` - runs **zCluster** on a single cluster from ``examples/400SDSmall.csv``
  for each of several databases, and checks the results against stored reference values. This is the
  suite to run routinely, e.g., after making changes to the code.
* ``retrievers.robot`` - checks that each of the supported photometric databases can still be queried.
  This is not part of the routine tests - see :ref:`RetrieverChecks` below.

All of the commands below are run from the root of the source distribution. Output (``log.html``,
``report.html``, and the files written by each test) is placed in the directory given by ``--outputdir``.


Regression tests
================

To run the regression tests:

.. code-block:: text

   robot --outputdir tests/results tests/zCluster_regression.robot

Each test checks that the photometric redshift of the cluster is close to its spectroscopic redshift,
and that the redshift and ``delta`` values match those stored in ``tests/reference.json``.

The galaxy photometry used by these tests is cached in ``tests/cache``, so only the first run needs to
query the remote databases. Subsequent runs therefore test the code, rather than the remote
services.

Tests are tagged with the name of the database they use, which can be used to run a subset of them, e.g.,

.. code-block:: text

   robot --outputdir tests/results --include SDSSDR12 tests/zCluster_regression.robot

To skip the tests that query the NOIRLab Data Lab, or the Rubin Science Platform, use
``--exclude datalab`` or ``--exclude rubin`` respectively. The Rubin tests are skipped automatically if
no access token is found in ``~/.rsp-tap.token`` (see
`here <https://rsp.lsst.io/guides/auth/creating-user-tokens.html>`__ for how to create one).

If a change to the code is expected to change the results (e.g., a change to the templates or the
redshift priors), the stored reference values can be updated with:

.. code-block:: text

   robot --outputdir tests/results --variable UPDATE_REFERENCE:True tests/zCluster_regression.robot

The updated ``tests/reference.json`` should then be committed along with the change. A reference
value is also written automatically the first time a test for a new database is run.


.. _RetrieverChecks:

Checking the database retrievers
================================

The services that host the photometric surveys supported by **zCluster** are occasionally retired, or
moved, or change their output format. The ``retrievers.robot`` suite checks for this, by querying each
database for galaxies in a small (4 arcmin) region at a position inside the survey footprint.

.. warning::

   These checks are meant to be run manually, on occasion - not as part of the routine tests. They
   always contact the remote services (nothing is cached between runs), and some take a long time to
   run (e.g., the DECaLS retrievers download whole tractor catalog files from NERSC). Because of this,
   always run the regression tests by giving the path to ``zCluster_regression.robot`` (as above),
   rather than running everything in the ``tests`` directory.

To run the retriever checks:

.. code-block:: text

   robot --outputdir tests/results tests/retrievers.robot

Each query times out after 90 seconds (600 seconds for the DECaLS tractor catalog downloads), which is
plenty of time for a working service. A test fails if the query times out, raises an error, or returns
fewer than 10 galaxies.

Some databases need a login or access token, and the corresponding tests are skipped if these are not
found:

* ``PS1`` - the ``CASJOBS_WSID`` and ``CASJOBS_PW`` environment variables
  (see `here <http://ps1images.stsci.edu/ps1_dr2_query.html>`__).
* ``DESDR1``, ``DESDR2`` - a DES database login in ``~/.desservices.ini`` (used by
  `easyaccess <https://github.com/des-labs/easyaccess>`_).
* ``RubinDP1``, ``RubinDP2`` - a Rubin Science Platform token in ``~/.rsp-tap.token``.

Tests are also skipped if a Python module that a retriever needs (e.g., ``easyaccess``) is not
installed.

As for the regression tests, tags can be used to run a subset of the checks. Each test is tagged with
its database name, and in addition:

* ``slow`` - the DECaLS tractor catalog downloads from NERSC.
* ``datalab`` - databases queried through the NOIRLab Data Lab.
* ``credentials`` - databases that need a login or token.
* ``rubin`` - the Rubin Observatory data previews.

For example, to run everything except the slow downloads:

.. code-block:: text

   robot --outputdir tests/results --exclude slow tests/retrievers.robot

A failure does not necessarily mean that a service has been retired - it may be down temporarily.
The output for each query is saved in ``tests/results/retrievers/<database>/retriever_stdout.log``,
which can help to tell the difference (e.g., an HTTP 502 error usually indicates a temporary outage,
while HTTP 404 suggests that the service has moved or been retired). If a check fails, it is worth
re-running it later before removing support for that database.
