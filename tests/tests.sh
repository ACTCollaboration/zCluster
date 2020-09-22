# Before setting up a proper testing framework (e.g., using Robot), some quick things to run
# 1. Basic run on SDSS (ok... this first one is slow if want to re-fetch the whole cache)
zCluster ../examples/400SDAll.csv SDSSDR12 -c testCache -o test400SD_SDSS
# 2. MPI run
mpiexec zCluster ../examples/400SDAll.csv SDSSDR12 -c testCache -o test400SD_SDSS_MPI -M
# 3. Comparisons script check(s)
zClusterComparisonPlot ../examples/400SDAll.csv zCluster_test400SD_SDSS.fits
zClusterComparisonPlot ../examples/400SDAll.csv zCluster_test400SD_SDSS_MPI.fits
# 4. zField
zField 0.0 0.0 -r 0.2 -D DECaLS -d -z 0.3 -delta_z 0.05 -DS9 -c testCache
