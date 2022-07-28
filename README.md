# xmascat 0.1.0 (Currently under construction)
*Xarray to Measurementset(v2) conversion package for ASte and other Common Astronomical radio Telescopes*

*xmascat refered to aste-xffts-merge (https://github.com/astropenguin/aste-xffts-merge) and b4rpipe (https://github.com/b4r-dev/b4rpipe).*

*May not work on MacOS (problem of casacore).*

## Requirement

* Python >=3.7
* pandas
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
xmascat.create_XFFTSxarray("./20211205031247/OKLDr.start.20211205031247", "./20211205031247/antlog_20211205031247", "./20211205031247/XFFTS.20211205031247.A02")
xmascat.create_XFFTSxarray("./20211205031247/OKLDr.start.20211205031247", "./20211205031247/antlog_20211205031247", "./20211205031247/XFFTS.20211205031247.A04")

# Convert xarray to MS2 (MeasurementSet v2)
xmascat.Xarray2MS2("./20211205031247/XFFTS.20211205031247.A02.nc")
xmascat.Xarray2MS2("./20211205031247/XFFTS.20211205031247.A04.nc")
```

## Notes

Examples: Convert MS2 to FITS using CASA

```python
### CASA
sdbaseline(infile="./20211205031247/XFFTS.20211205031247.A02.ms", outfile="./20211205031247/XFFTS.20211205031247.A02.bl.ms", datacolumn="float_data", spw="0:15000~16000;17000~18000“, blfunc="poly", order=1, overwrite=True)

sdbaseline(infile="./20211205031247/XFFTS.20211205031247.A04.ms", outfile="./20211205031247/XFFTS.20211205031247.A04.bl.ms", datacolumn="float_data", spw="0:15000~16000;17000~18000“, blfunc="poly", order=1, overwrite=True)

import math

gencal(vis="./20211205031247/XFFTS.20211205031247.A02.bl.ms", caltable="./20211205031247/XFFTS.20211205031247.A02.bl.mb.tbl", caltype="amp", parameter=[math.sqrt(0.45)])

gencal(vis="./20211205031247/XFFTS.20211205031247.A04.bl.ms", caltable="./20211205031247/XFFTS.20211205031247.A04.bl.mb.tbl", caltype="amp", parameter=[math.sqrt(0.45)])

applycal(vis="./20211205031247/XFFTS.20211205031247.A02.bl.ms", gaintable=["./20211205031247/XFFTS.20211205031247.A02.bl.mb.tbl"], calwt=[False])

applycal(vis="./20211205031247/XFFTS.20211205031247.A04.bl.ms", gaintable=["./20211205031247/XFFTS.20211205031247.A04.bl.mb.tbl"], calwt=[False])

sdimaging(infiles=["./20211205031247/XFFTS.20211205031247.A02.bl.ms", "./20211205031247/XFFTS.20211205031247.A04.bl.ms/"], outfile="./20211205031247/XFFTS.20211205031247.A02A04.bl.int", intent="*ON_SOURCE*", gridfunction="GAUSS", cell=["10arcsec", "10arcsec"], mode="velocity", nchan=201, start="-50.0km/s", width="0.5km/s", overwrite=True, imsize=[100, 100], phasecenter="J2000 5h35m14.16 -5d22m21.5", restfreq="345.795990GHz")

exportfits(imagename="./20211205031247/XFFTS.20211205031247.A02A04.bl.int", fitsimage="./20211205031247/XFFTS.20211205031247.A02A04.bl.int.fits", velocity=True, dropstokes=True, overwrite=True)
```