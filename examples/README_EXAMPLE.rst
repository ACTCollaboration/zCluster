Here we provide a quick example of running **zCluster** on part of
the `Burenin et al. (2007) <https://ui.adsabs.harvard.edu/abs/2007ApJS..172..561B/abstract>`_ 400SD cluster catalog.

.. note::  The catalog files needed for this tutorial can be
           found in the `examples/ <https://github.com/ACTCollaboration/zCluster/tree/master/examples/>`_
           directory of the **zCluster** source distribution.

There are two cluster catalog files provided, ``400SDAll.csv`` and ``400SDSmall.csv``. The latter contains just the
first 20 rows of the former, and is enough to check that the code is working. 

To run :ref:`zClusterCommand` using galaxy photometry from SDSS DR12, use:

.. code-block::

   zCluster 400SDSmall.csv SDSSDR12

To plot a comparison of the results with the spectroscopic redshifts in the 400SD catalog, run:

.. code-block::

   zClusterComparisonPlot 400SDSmall.csv zCluster_400SDSmall_SDSSDR12.fits

These examples use the default options - to see others, run each command with the ``-h`` flag.
