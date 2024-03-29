# xmascat 0.1.0 (Currently under construction)
*Xarray to Measurementset(v2) conversion package for ASte and other Common Astronomical radio Telescopes*

*xmascat refered to aste-xffts-merge (https://github.com/astropenguin/aste-xffts-merge) and b4rpipe (https://github.com/b4r-dev/b4rpipe).*

*May not work on MacOS (problem of casacore).*

## Requirement

* Python >=3.7
* pandas
* scipy
* xarray
* xarray_dataclasses
* astropy
* python-casacore
* netCDF4
* (joblib)

## Installation

```bash
pip install git+https://github.com/ShinjiFujita/xmascat.git
```

## Usage

```python
import xmascat

# create xarray from ASTE observation logs and XFFTS data. 
## OTF
xmascat.create_XFFTSxarray(path_startfile="./20211205031247/OKLDr.start.20211205031247", path_antlogfile="./20211205031247/antlog_20211205031247", path_XFFTSdata="./20211205031247/XFFTS.20211205031247.A02", path_messfiles="../astref2_2023/20*/mess*.xz")
xmascat.create_XFFTSxarray(path_startfile="./20211205031247/OKLDr.start.20211205031247", path_antlogfile="./20211205031247/antlog_20211205031247", path_XFFTSdata="./20211205031247/XFFTS.20211205031247.A04", path_messfiles="../astref2_2023/20*/mess*.xz")
## PS
xmascat.create_XFFTSxarray(path_startfile="./20211205031247/OKLDr.start.20211205031247", path_XFFTSdata="./20211205031247/XFFTS.20211205031247.A02", path_messfiles="../astref2_2023/20*/mess*.xz")
xmascat.create_XFFTSxarray(path_startfile="./20211205031247/OKLDr.start.20211205031247", path_XFFTSdata="./20211205031247/XFFTS.20211205031247.A04", path_messfiles="../astref2_2023/20*/mess*.xz")

# Plot xarray
xmascat.plotnc("./20211205031247/XFFTS.20211205031247.A02.nc", xmin=None, xmax=None, ymin=None, ymax=None)

# Flag
xmascat.delete_scans("./20211205031247/XFFTS.20211205031247.A02.nc", [5, 6]) # e.g., delete the scan 0005 and 0006 

# Slice
xmascat.slice_chs("./20211205031247/XFFTS.20211205031247.A02.nc", 345.79598990e9, 4096, "12CO32") # e.g., slice from 345.79598990e9 Hz - 4096ch to 345.79598990e9 Hz + 4096ch. "./20211205031247/XFFTS.20211205031247.A02.12CO32.nc" to be saved.

# Convert xarray to MS2 (MeasurementSet v2)
## OTF
xmascat.Xarray2MS2("./20211205031247/XFFTS.20211205031247.A02.nc")
xmascat.Xarray2MS2("./20211205031247/XFFTS.20211205031247.A04.nc")
```

## Notes

Examples: Convert MS2 to FITS using CASA

```bash
wget https://casa.nrao.edu/download/distro/casa/release/rhel/casa-6.5.5-21-py3.8.tar.xz
xz -dv casa-6.5.5-21-py3.8.tar.xz
tar xfv casa-6.5.5-21-py3.8.tar

./casa-6.5.5-21-py3.8/bin/casa
```

```python
### CASA
sdbaseline(infile="./20211205031247/XFFTS.20211205031247.A02.ms", outfile="./20211205031247/XFFTS.20211205031247.A02.bl.ms", datacolumn="float_data", spw="0:15000~16000;17000~18000", blfunc="poly", order=1, overwrite=True)

sdbaseline(infile="./20211205031247/XFFTS.20211205031247.A04.ms", outfile="./20211205031247/XFFTS.20211205031247.A04.bl.ms", datacolumn="float_data", spw="0:15000~16000;17000~18000", blfunc="poly", order=1, overwrite=True)

import math

gencal(vis="./20211205031247/XFFTS.20211205031247.A02.bl.ms", caltable="./20211205031247/XFFTS.20211205031247.A02.bl.mb.tbl", caltype="amp", parameter=[math.sqrt(0.45)])

gencal(vis="./20211205031247/XFFTS.20211205031247.A04.bl.ms", caltable="./20211205031247/XFFTS.20211205031247.A04.bl.mb.tbl", caltype="amp", parameter=[math.sqrt(0.45)])

applycal(vis="./20211205031247/XFFTS.20211205031247.A02.bl.ms", gaintable=["./20211205031247/XFFTS.20211205031247.A02.bl.mb.tbl"], calwt=[False])

applycal(vis="./20211205031247/XFFTS.20211205031247.A04.bl.ms", gaintable=["./20211205031247/XFFTS.20211205031247.A04.bl.mb.tbl"], calwt=[False])

sdimaging(infiles=["./20211205031247/XFFTS.20211205031247.A02.bl.ms", "./20211205031247/XFFTS.20211205031247.A04.bl.ms/"], outfile="./20211205031247/XFFTS.20211205031247.A02A04.bl.int", intent="*ON_SOURCE*", gridfunction="GAUSS", cell=["10arcsec", "10arcsec"], mode="velocity", nchan=201, start="-50.0km/s", width="0.5km/s", overwrite=True, imsize=[100, 100], phasecenter="J2000 5h35m14.16 -5d22m21.5", restfreq="345.795990GHz", outframe="lsrk")

exportfits(imagename="./20211205031247/XFFTS.20211205031247.A02A04.bl.int", fitsimage="./20211205031247/XFFTS.20211205031247.A02A04.bl.int.fits", velocity=True, dropstokes=True, overwrite=True)
```