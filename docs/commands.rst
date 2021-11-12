.. _Usage:

=================
zCluster Commands
=================

The **zCluster** package includes a number of command-line programs, each of which is described below.


.. _zClusterCommand:
    
zCluster
--------

.. argparse::
   :filename: ../bin/zCluster
   :func: makeParser
   :prog: zCluster
   
   :program:`zCluster` estimates galaxy cluster photometric redshifts for objects in a user-supplied
   cluster catalog.


.. _zFieldCommand:
    
zField
--------

.. argparse::
   :filename: ../bin/zField
   :func: makeParser
   :prog: zField
   
   :program:`zField` estimates galaxy photometric redshifts for objects in a user-supplied
   direction on the sky.


