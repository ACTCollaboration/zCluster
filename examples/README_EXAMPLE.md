To run zCluster on a small number of 400SD ([Burenin et al. 2007](http://adsabs.harvard.edu/abs/2007ApJS..172..561B)) 
cluster catalog positions:

```
zCluster 400SDSmall.csv SDSSDR12
```

To plot a comparison of the zCluster results with the redshifts in the 400SD catalog:

```
zClusterComparisonPlot 400SDSmall.csv zCluster_400SDSmall_SDSSDR12.fits
```

These examples use the default options - to see others, run each code with the `-h` flag.
