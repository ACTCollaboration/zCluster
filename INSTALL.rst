**zCluster** is written in `Python <https://www.python.org/>`_ (3.6+), and requires the
following additional modules to be installed (currently used versions are given in
brackets, later versions also probably work):


* numpy (1.13.3)
* scipy (1.3.0)
* matplotlib (2.1.1)
* astLib (0.11.7)
* astropy (4.0)
* mastcasjobs (for PS1; https://github.com/rlwastro/mastcasjobs)
* casjobs (for PS1; https://github.com/dfm/casjobs)

To run on DES photometry, there is an additional dependency:

* easyaccess (1.4.5)

If you want to run the code in parallel, you will also need:

* mpi4py (3.0.0)

Note that if you want to run the code on a cluster, the bottleneck will be fetching the photometric catalogs
over the internet. The MPI mode is still useful though on any machine with multiple cores.

The latest tagged version of **zCluster** can be installed using ``pip``:
    
.. code-block::

   pip install zCluster

Other dependencies will be installed by ``pip``.

You may also install using the standard ``setup.py`` script, e.g., as root:

.. code-block::

   sudo python setup.py install

Alternatively, 

.. code-block::

   python setup.py install --user

will install ``zCluster`` under ``$HOME/.local`` (on Ubuntu), and in some other default location on Mac.

You can also use the ``--prefix`` option, e.g.,

.. code-block::

   python setup.py install --prefix=$HOME/local

and then add ``$HOME/local/bin`` to $PATH, and e.g., ``$HOME/local/lib/python3.6/site-packages`` to 
$PYTHONPATH (adjust the path according to your Python version number).

.. code-block::

   export PATH=$HOME/local/bin:$PATH    
   export PYTHONPATH=$HOME/local/lib/python3.6/site-packages:$PYTHONPATH

If **zCluster** has installed correctly, then you should find its command line tools are available, for
example,

.. code-block::
   
   zCluster -h
   
should display a helpful message about the command-line options for the main ``zCluster`` command.
