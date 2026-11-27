"""Convert Rubin/LSST baseline total throughput curves into zCluster passband files.

Download total_{u,g,r,i,z,y}.dat from a release of https://github.com/lsst/throughputs, e.g.

    for b in u g r i z y; do
        curl -sfLO https://raw.githubusercontent.com/lsst/throughputs/v1.9.2/baseline/total_$b.dat
    done

then run:

    python convert_lsst.py <dir containing total_*.dat> ../zCluster/passbands

Output files ({band}LSST.res) are in Angstroms, readable by astSED.Passband with default options.
If you use a different release, update the release tag in the header text written below.

"""
import sys, numpy as np
inDir, outDir = sys.argv[1], sys.argv[2]
for b in "ugrizy":
    hdr = [l for l in open(f"{inDir}/total_{b}.dat") if l.startswith("#")]
    d = np.loadtxt(f"{inDir}/total_{b}.dat")
    wl = np.round(d[:, 0] * 10.0, 3); tr = d[:, 1]
    assert (tr >= 0).all() and np.all(np.diff(wl) > 0)
    nz = np.where(tr > 0)[0]
    lo, hi = max(nz[0] - 1, 0), min(nz[-1] + 1, len(tr) - 1)
    with open(f"{outDir}/{b}LSST.res", "w") as f:
        f.write(f"# Rubin Observatory LSST {b}-band total system throughput (baseline)\n")
        f.write("# Source: https://github.com/lsst/throughputs release v1.9.2, baseline/total_%s.dat\n" % b)
        f.write("# Includes mirrors (3Ag), lenses, filter, detector, time-averaged losses, and X=1.2 standard atmosphere\n")
        f.write("# Converted from nm to Angstroms; leading/trailing zero-throughput padding removed\n")
        f.write("# Original header:\n")
        for l in hdr: f.write("#   " + l.lstrip("# ").rstrip() + "\n")
        f.write("# Wavelength(Angstrom)  Throughput(0-1)\n")
        for w, t in zip(wl[lo:hi+1], tr[lo:hi+1]):
            f.write(f"{w:.1f} {t:.8e}\n")
    print(b, wl[lo], wl[hi], hi - lo + 1, "rows; peak", tr.max().round(3))
