"""Robot Framework keyword library for zCluster regression tests.

Runs bin/zCluster from a given source tree (not whatever version happens to be installed) and compares
the output against stored reference results.

"""

import os
import sys
import json
import shutil
import subprocess
import numpy as np
import astropy.table as atpy
from robot.api import logger, SkipExecution
from robot.api.deco import keyword, library

@library(scope = 'GLOBAL')
class ZClusterLibrary:

    def __init__(self, zcluster_root = None):
        if zcluster_root is None:
            zcluster_root=os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
        self.root=os.path.abspath(zcluster_root)

    @keyword
    def run_zcluster(self, catalog, database, cluster_name, work_dir, cache_dir, timeout = 1800, *extra_args):
        """Runs zCluster on a single cluster, returning the path to the output .fits table.

        The work directory is wiped first, as zCluster otherwise skips clusters already in its .db files.
        stdin is /dev/null, so a stray IPython.embed() exits straight away instead of hanging the test.

        """
        label=database
        work_dir=os.path.abspath(work_dir)
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)
        os.makedirs(work_dir)
        os.makedirs(cache_dir, exist_ok = True)

        cmd=[sys.executable, os.path.join(self.root, 'bin', 'zCluster'), os.path.abspath(catalog), database,
             '-c', os.path.abspath(cache_dir), '-o', label, '-n', str(cluster_name)]+list(extra_args)
        env=dict(os.environ)
        env['PYTHONPATH']=self.root+os.pathsep+env.get('PYTHONPATH', '')
        env['MPLBACKEND']='Agg'
        logger.info("Running: %s" % (" ".join(cmd)))
        logger.info("Source tree: %s" % (self.root))
        try:
            proc=subprocess.run(cmd, cwd = work_dir, env = env, stdin = subprocess.DEVNULL,
                                stdout = subprocess.PIPE, stderr = subprocess.STDOUT, text = True,
                                timeout = float(timeout))
        except subprocess.TimeoutExpired as e:
            raise AssertionError("zCluster timed out after %s s" % (timeout))
        with open(os.path.join(work_dir, 'zCluster_stdout.log'), 'w') as outFile:
            outFile.write(proc.stdout)
        logger.info(proc.stdout[-5000:])

        if 'In [1]:' in proc.stdout:
            raise AssertionError("zCluster hit an IPython.embed() breakpoint - see zCluster_stdout.log")
        if proc.returncode != 0:
            raise AssertionError("zCluster exited with code %d - see zCluster_stdout.log" % (proc.returncode))
        outFileName=os.path.join(work_dir, 'zCluster_%s.fits' % (label))
        if not os.path.exists(outFileName):
            raise AssertionError("zCluster produced no output table (%s)" % (outFileName))

        return outFileName

    @keyword
    def retriever_should_return_galaxies(self, database, ra_deg, dec_deg, half_box_deg, work_dir,
                                         min_galaxies = 10, timeout = 600):
        """Queries the given database for galaxies in a small box around the given position, using a fresh
        (empty) cache, so that the remote service is always contacted. Fails if the query errors, or returns
        fewer than min_galaxies objects. Skips the test if a module needed by the retriever is not installed.

        Runs in a subprocess (against this source tree), so a hung query can be timed out.

        """
        work_dir=os.path.abspath(work_dir)
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)
        os.makedirs(work_dir)
        script="""
import sys, json
from zCluster import retrievers
retriever, options, passbandSet=retrievers.getRetriever(sys.argv[1])
if retriever is None:
    raise Exception("getRetriever() does not recognise database '%s'" % (sys.argv[1]))
options=dict(options) if options is not None else {}
options['altCacheDir']=sys.argv[5]
catalog=retriever(float(sys.argv[2]), float(sys.argv[3]), halfBoxSizeDeg = float(sys.argv[4]), optionsDict = options)
if isinstance(catalog, str):
    raise Exception("retriever returned '%s' instead of a catalog (bad response from the service?)" % (catalog))
keys=sorted(catalog[0].keys()) if catalog else []
print("RETRIEVER_RESULT "+json.dumps({'numGalaxies': None if catalog is None else len(catalog),
                                      'passbandSet': passbandSet, 'keys': keys}))
"""
        cmd=[sys.executable, '-c', script, database, str(ra_deg), str(dec_deg), str(half_box_deg),
             os.path.join(work_dir, 'cache')]
        env=dict(os.environ)
        env['PYTHONPATH']=self.root+os.pathsep+env.get('PYTHONPATH', '')
        env['MPLBACKEND']='Agg'
        env['PYTHONUNBUFFERED']='1'    # so that the traceback is the last thing in the output
        logger.info("Querying %s at (%s, %s), half box size %s deg" % (database, ra_deg, dec_deg, half_box_deg))
        try:
            proc=subprocess.run(cmd, cwd = work_dir, env = env, stdin = subprocess.DEVNULL,
                                stdout = subprocess.PIPE, stderr = subprocess.STDOUT, text = True,
                                timeout = float(timeout))
        except subprocess.TimeoutExpired as e:
            with open(os.path.join(work_dir, 'retriever_stdout.log'), 'w') as outFile:
                outFile.write(e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or ''))
            raise AssertionError("%s query timed out after %s s" % (database, timeout))
        with open(os.path.join(work_dir, 'retriever_stdout.log'), 'w') as outFile:
            outFile.write(proc.stdout)
        logger.info(proc.stdout[-5000:])

        if proc.returncode != 0:
            missing=[l for l in proc.stdout.splitlines() if l.startswith('ModuleNotFoundError')]
            if missing:
                raise SkipExecution("%s needs a module that is not installed: %s" % (database, missing[-1]))
            lastLine=proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else ''
            raise AssertionError("%s query failed: %s" % (database, lastLine))
        resultLines=[l for l in proc.stdout.splitlines() if l.startswith('RETRIEVER_RESULT ')]
        if not resultLines:
            raise AssertionError("%s query produced no result - see retriever_stdout.log" % (database))
        result=json.loads(resultLines[-1][len('RETRIEVER_RESULT '):])
        logger.info("Result: %s" % (result))
        if result['numGalaxies'] is None:
            raise AssertionError("%s retriever returned None (position outside the survey area?)" % (database))
        if result['numGalaxies'] < int(min_galaxies):
            raise AssertionError("%s returned %d galaxies (expected at least %s)"
                                 % (database, result['numGalaxies'], min_galaxies))
        for k in ['RADeg', 'decDeg']:
            if k not in result['keys']:
                raise AssertionError("%s catalog has no '%s' key" % (database, k))

        return result['numGalaxies']

    @keyword
    def get_cluster_result(self, fits_path, cluster_name):
        """Returns the row for the given cluster from a zCluster output table, as a dictionary.

        """
        tab=atpy.Table().read(fits_path)
        rows=tab[tab['name'].astype(str) == str(cluster_name)]
        if len(rows) == 0:
            raise AssertionError("No redshift for cluster '%s' in %s (zCluster drops clusters with z <= 0)"
                                 % (cluster_name, fits_path))
        result={k: float(rows[0][k]) for k in ['z', 'delta', 'errDelta', 'offsetArcmin', 'RADeg', 'decDeg']}
        logger.info("Result: %s" % (result))
        return result

    @keyword
    def get_catalog_redshift(self, catalog, cluster_name):
        """Returns the (spectroscopic) redshift listed for the cluster in the input .csv catalog.

        """
        with open(catalog) as inFile:
            for line in inFile:
                if line.startswith('#') or not line.strip():
                    continue
                bits=line.split()
                if bits[0] == str(cluster_name):
                    return float(bits[3])
        raise AssertionError("Cluster '%s' not found in %s" % (cluster_name, catalog))

    @keyword
    def photo_z_should_be_close_to(self, result, z_spec, max_dz_over_1pz = 0.1):
        """Sanity check: |z_phot - z_spec| / (1 + z_spec) must be below the given threshold.

        """
        dz=abs(result['z']-float(z_spec))/(1+float(z_spec))
        if dz > float(max_dz_over_1pz):
            raise AssertionError("Photo-z %.3f too far from z_spec %.3f (|dz|/(1+z) = %.3f > %s)"
                                 % (result['z'], float(z_spec), dz, max_dz_over_1pz))

    @keyword
    def result_should_match_reference(self, result, reference_file, database, z_tol = 0.02, delta_rel_tol = 0.1,
                                      update = False):
        """Compares result to the stored reference for this database.

        If update is True (or no reference exists yet for this database), the reference file is written
        instead, and the test passes with a warning.

        """
        if os.path.exists(reference_file):
            with open(reference_file) as inFile:
                refs=json.load(inFile)
        else:
            refs={}
        if str(update).lower() in ['true', '1', 'yes'] or database not in refs:
            refs[database]={'z': round(result['z'], 4), 'delta': round(result['delta'], 4),
                            'errDelta': round(result['errDelta'], 4)}
            with open(reference_file, 'w') as outFile:
                json.dump(refs, outFile, indent = 4, sort_keys = True)
                outFile.write("\n")
            logger.warn("Wrote reference result for %s to %s: %s" % (database, reference_file, refs[database]))
            return

        ref=refs[database]
        errors=[]
        if abs(result['z']-ref['z']) > float(z_tol):
            errors.append("z = %.4f, reference = %.4f (tol %s)" % (result['z'], ref['z'], z_tol))
        if not np.isclose(result['delta'], ref['delta'], rtol = float(delta_rel_tol), atol = 0):
            errors.append("delta = %.4f, reference = %.4f (rel. tol %s)" % (result['delta'], ref['delta'],
                                                                            delta_rel_tol))
        if errors:
            raise AssertionError("%s regression: %s" % (database, "; ".join(errors)))
