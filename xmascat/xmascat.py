import numpy as np
import pandas as pd
import re

def del_dict_key(dict_inp, keys):
    import copy
    d = copy.deepcopy(dict_inp)
    for k in keys:
        try:
            del d[k]
        except:
            pass
    return d
    
def numbering_PTN_list(PTN_list):
    res = []
    num_R, num_SKY, num_ON, num_OFF, num_TRANS = 0, 0, 0, 0, 0
    for PTN in PTN_list:
        if PTN=="R":
            res.append("R_%s"%str(num_R).zfill(4))
            num_R += 1
        elif PTN=="SKY":
            res.append("SKY_%s"%str(num_SKY).zfill(4))
            num_SKY += 1
        elif PTN=="ON":
            res.append("ON_%s"%str(num_ON).zfill(4))
            num_ON += 1
        elif PTN=="OFF":
            res.append("OFF_%s"%str(num_OFF).zfill(4))
            num_OFF += 1
        elif PTN=="TRANS":
            res.append("TRANS_%s"%str(num_TRANS).zfill(4))
            num_TRANS += 1
        else:
            print("neko")
    return res

def read_startfile(filepath):
    list_all = []
    with open(filepath) as f:
        for line in f:
            list_all.append([_ for _ in re.split("[ \n]", line) if _!=""])
    list_XFFTS = [_[1:] for _ in list_all if "XFFTS" in _]
    
    SET_dict = {}
    for line in list_XFFTS:
        if line[0]=="OPEN":
            continue
        elif line[0]=="SET":
            SET_dict[line[2]] = line[3][1:-1]
            continue
        else:
            continue
    SET_dict = del_dict_key(SET_dict, ["DUMMY_MODE", "OTF_MODE", "INTEG_TIME", "CALB_INT"])
    
    PTN_list = []
    for line in list_XFFTS:
        if line[0]=="EXECUTE":
            if line[2]=="TYPE(ON)" or line[2]=="TYPE(OFF)" or line[2]=="TYPE(R)" or line[2]=="TYPE(SKY)":
                PTN_list.append(line[2][5:-1])
    PTN_list = numbering_PTN_list(PTN_list)
    return SET_dict, PTN_list








########################
#####  aste-xffts-merge

from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Tuple, Union
import os
import numpy as np
import datetime

import pandas as pd
import xarray as xr
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, FK5, SkyCoord
from xarray_dataclasses import AsDataset, Attr, Coordof, Data, Dataof

"""Submodule of common objects within the package."""


# standard library
from dataclasses import dataclass, field
from typing import Any, Literal, TypeVar


# dependencies
from xarray_dataclasses import Attr, Data


# constants
DIMS = "time", "chan"


# type hints
T = TypeVar("T")
time = Literal["time"]
chan = Literal["chan"]


# dataclasses
def const(default: T, **kwargs: Any) -> T:
    """Create a constant field for dataclasses."""
    return field(default=default, init=False, **kwargs)


@dataclass
class Time:
    """Time in UTC."""

    data: Data[time, Literal["M8[ns]"]]
    long_name: Attr[str] = const("Time in UTC")
    short_name: Attr[str] = const("Time")


@dataclass
class Chan:
    """Channel ID."""

    data: Data[chan, int]
    long_name: Attr[str] = const("Channel ID")
    short_name: Attr[str] = const("Channel")
    
"""Submodule of ASTE antenna logs."""


# standard library
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Tuple, Union


# dependencies
import pandas as pd
import xarray as xr
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, FK5, SkyCoord
from xarray_dataclasses import AsDataset, Attr, Coordof, Data, Dataof


# submodules
#from .common import Time, const, time


# constants
ASTE_SITE = EarthLocation.from_geodetic(
    lon="-67d42m11.89525s",
    lat="-22d58m17.69447s",
    height="4861.9m",
)
LOG_COLUMNS = (
    "time",
    "ra_prog",
    "dec_prog",
    "az_prog",
    "el_prog",
    "az_real",
    "el_real",
    "az_error",
    "el_error",
)
LOG_FRAME = "RADEC"
LOG_SEPARATOR = r"\s+"
LOG_TIMEFMT = "%y%m%d%H%M%S.%f"


# dataclasses
@dataclass
class Azimuth:
    """Antenna azimuth (degree)."""

    data: Data[time, float]
    long_name: Attr[str] = const("Antenna azimuth")
    short_name: Attr[str] = const("Azimuth")
    units: Attr[str] = const("degree")


@dataclass
class Elevation:
    """Antenna elevation (degree)."""

    data: Data[time, float]
    long_name: Attr[str] = const("Antenna elevation")
    short_name: Attr[str] = const("Elevation")
    units: Attr[str] = const("degree")


@dataclass
class Longitude:
    """Sky longitude (degree)."""

    data: Data[time, float]
    long_name: Attr[str] = const("Sky longitude")
    short_name: Attr[str] = const("Longitude")
    units: Attr[str] = const("degree")


@dataclass
class Latitude:
    """Sky latitude (degree)."""

    data: Data[time, float]
    long_name: Attr[str] = const("Sky latitude")
    short_name: Attr[str] = const("Latitude")
    units: Attr[str] = const("degree")


@dataclass
class Frame:
    """Sky coordinate frame."""

    data: Data[Tuple[()], str]
    long_name: Attr[str] = const("Sky coordinate frame")
    short_name: Attr[str] = const("Frame")


@dataclass
class Antenna(AsDataset):
    """ASTE antenna log."""

    time: Coordof[Time]
    """Time in UTC."""

    azimuth: Dataof[Azimuth]
    """Antenna azimuth (degree)."""

    elevation: Dataof[Elevation]
    """Antenna elevation (degree)."""

    longitude: Dataof[Longitude]
    """Sky longitude (degree)."""

    latitude: Dataof[Latitude]
    """Sky latitude (degree)."""

    frame: Attr[str]
    """Sky coordinate frame."""


# runtime functions
def read_antlogfile(path: Union[Path, str], save=False) -> xr.Dataset:
    """Read an antenna log and create a Dataset object.

    Args:
        path: Path of the antenna log.

    Returns:
        A Dataset object that follows ``Antenna``.

    """
    # check if the sky coordinate frame is supported
    with open(path) as f:
        frame = f.readline().split()[0]

    if not frame == LOG_FRAME:
        raise ValueError(f"RADEC is only supported. Got {frame}.")

    # read the antenna log
    date_parser = partial(pd.to_datetime, format=LOG_TIMEFMT)

    log = pd.read_csv(
        path,
        date_parser=date_parser,
        index_col=0,
        names=LOG_COLUMNS,
        sep=LOG_SEPARATOR,
        skiprows=1,
    )

    # calculate real-prog differences in sky
    sky_prog = SkyCoord(
        alt=log[LOG_COLUMNS[4]],
        az=log[LOG_COLUMNS[3]],
        frame=AltAz,
        location=ASTE_SITE,
        obstime=log.index,
        unit=u.deg,  # type: ignore
    ).transform_to(FK5)

    sky_real = SkyCoord(
        alt=log[LOG_COLUMNS[6]],
        az=log[LOG_COLUMNS[5]],
        frame=AltAz,
        location=ASTE_SITE,
        obstime=log.index,
        unit=u.deg,  # type: ignore
    ).transform_to(FK5)

    d_lon = (sky_real.ra - sky_prog.ra).deg  # type: ignore
    d_lat = (sky_real.dec - sky_prog.dec).deg  # type: ignore
    
    res = Antenna.new(
        time=log.index,
        azimuth=log[LOG_COLUMNS[5]],
        elevation=log[LOG_COLUMNS[6]],
        longitude=log[LOG_COLUMNS[1]] + d_lon,
        latitude=log[LOG_COLUMNS[2]] + d_lat,
        frame=FK5.name,  # type: ignore
    )
    
    if save==True:
        res.to_netcdf(path+".nc")
    return res
    
                                       #####
########################






import itertools
import os
import numpy as np
import datetime

def read_XFFTSdata(filename, PTN_list, nchan=32768):
    fp = open(filename, 'rb')
    filesize = os.path.getsize(filename)
    headersize = 20
    onedatasize = nchan*4
    num_data = filesize/(onedatasize+headersize)
    
    dict_list = []
    for i in range(int(num_data)):
        dict_tempo = {}
        header_tempo = fp.read(headersize)
        data_tempo = fp.read(onedatasize)
        timestamp_tempo = np.frombuffer(header_tempo[:header_tempo.find(b"\xca=")-2], dtype="float64")[0]
        integtime_tempo = np.frombuffer(header_tempo[header_tempo.find(b"\xca=")-2:header_tempo.find(b"\xca=")+2], dtype="float32")[0]
        scantype_start_ind = header_tempo.find(b"\xca=")+2
        scantype_end_ind = header_tempo.find(b"\x00\x00\x00\x00")
        if scantype_start_ind<scantype_end_ind:
            scantype_tempo = header_tempo[scantype_start_ind:scantype_end_ind]
        else:
            scantype_tempo = b"TRANS"
        dict_tempo["timestamp"] = datetime.datetime.utcfromtimestamp(timestamp_tempo).strftime("%Y-%m-%dT%H:%M:%S.%f%Z000")
        dict_tempo["integtime"] = integtime_tempo
        dict_tempo["scantype"] = scantype_tempo.decode()
        dict_tempo["data"] = np.frombuffer(data_tempo, dtype='f')
        dict_list.append(dict_tempo)
        
    dict_list_cw = [d for d in dict_list if (d["scantype"]!="ZERO")]
    
    timestamp_xffts_list = [d["timestamp"] for d in dict_list_cw if not d["scantype"]=="TRANS"]
    integtime_xffts_list = [d["integtime"] for d in dict_list_cw if not d["scantype"]=="TRANS"]
    scantype_xffts_list_pre = [d["scantype"] for d in dict_list_cw]
    data_xffts_list = [d["data"] for d in dict_list_cw if not d["scantype"]=="TRANS"]
    
    scantype_xffts_list = []
    count = 0
    for l in [(k, list(g)) for k, g in itertools.groupby(scantype_xffts_list_pre)]:
        if l[0]=="TRANS":
            continue
        else:
            for _ in range(len(l[1])):
                scantype_xffts_list.append(PTN_list[count])
            count += 1
    
    return timestamp_xffts_list, integtime_xffts_list, scantype_xffts_list, data_xffts_list










########
def create_XFFTSxarray(path_startfile, path_antlogtfile, path_XFFTSdata, Tamb=270.0):
	if (os.path.exists(path_startfile)==True and os.path.exists(path_antlogtfile) and os.path.exists(path_XFFTSdata)):
		SET_dict, PTN_list = read_startfile(filepath)
		antlog_xr = read_antlogfile(path_antlogtfile)
		timestamp_xffts_list, integtime_xffts_list, scantype_xffts_list, data_xffts_list = read_XFFTSdata(path_XFFTSdata, PTN_list)
		
		xr_cut = antlog_xr.sel(time=timestamp_xffts_list, method="nearest")
		xr_test_cut["integtime"] = (("time"), integtime_xffts_list)
		xr_test_cut["scantype"] = (("time"), scantype_xffts_list)
		xr_test_cut["ch"] = [i for i in range(nchan)]
		xr_test_cut["data"] = (("time", "ch"), data_xffts_list)
		
		A_num = int(path_XFFTSdata[-2:])
		for k in SET_dict.keys():
    		if k in ["REF_NUM", "RX_NAME", "REST_FREQ", "OBS_FREQ", "SIDBD_TYP", "FREQ_IF1"]:
        		xr_test_cut.attrs[k] = SET_dict[k].split(",")[A_num-1]
    		else:
        		xr_test_cut.attrs[k] = SET_dict[k]
        		
		if xr_test_cut.SIDBD_TYP=="USB":
    		xr_test_cut["freq"] = (("ch"), np.linspace(float(xr_test_cut.REST_FREQ) - 1.25e9, float(xr_test_cut.REST_FREQ) + 1.25e9, num=nchan))
		elif xr_test_cut.SIDBD_TYP=="LSB":
    		xr_test_cut["freq"] = (("ch"), np.linspace(float(xr_test_cut.REST_FREQ) + 1.25e9, float(xr_test_cut.REST_FREQ) - 1.25e9, num=nchan))
		else:
    		print("SIDBD_TYP is invalid. ")
    		
    	R_list = [_ for _ in PTN_list if _[:1]=="R"]
		SKY_list = [_ for _ in PTN_list if _[:3]=="SKY"]
		ON_list = [_ for _ in PTN_list if _[:2]=="ON"]
		OFF_list = [_ for _ in PTN_list if _[:3]=="OFF"]
		meantime_R_list = [np.mean((np.array(xr_test_cut["time"])[xr_test_cut.scantype==_]).astype("float64")).astype('datetime64[ns]') for _ in R_list]
		meantime_SKY_list = [np.mean((np.array(xr_test_cut["time"])[xr_test_cut.scantype==_]).astype("float64")).astype('datetime64[ns]') for _ in SKY_list]
		meantime_ON_list = [np.mean((np.array(xr_test_cut["time"])[xr_test_cut.scantype==_]).astype("float64")).astype('datetime64[ns]') for _ in ON_list]
		meantime_OFF_list = [np.mean((np.array(xr_test_cut["time"])[xr_test_cut.scantype==_]).astype("float64")).astype('datetime64[ns]') for _ in OFF_list]

		spe_array_list = []
		Tsys_median_list = []
		for i in range(len(ON_list)):
    		print("A%s"%str(A_num).zfill(2), ":", i+1, "/",  len(ON_list))
    		ON = ON_list[i]
    		time_ON = meantime_ON_list[i]
    		R = R_list[np.argmin(abs(meantime_R_list - meantime_ON_list[i]))]
    		SKY = SKY_list[np.argmin(abs(meantime_SKY_list - meantime_ON_list[i]))]
    		OFF = OFF_list[np.argmin(abs(meantime_OFF_list - meantime_ON_list[i]))]
    		ave_R = np.mean(np.array(xr_test_cut["data"])[xr_test_cut.scantype==R], axis=0)
    		ave_SKY = np.mean(np.array(xr_test_cut["data"])[xr_test_cut.scantype==SKY], axis=0)
    		ave_OFF = np.mean(np.array(xr_test_cut["data"])[xr_test_cut.scantype==OFF], axis=0)
    		Y = ave_R/ave_SKY
    		Tsys = Tamb/(Y-1.0)    
    		raw_ON_array = np.array(xr_test_cut["data"])[xr_test_cut.scantype==ON]
    		spe_array = np.array([Tamb*(raw_ON - ave_OFF)/(ave_R - ave_OFF) for raw_ON in raw_ON_array])
    		spe_array_list.append(spe_array)    
    		Tsys_median = [np.median(Tsys)]*len(raw_ON_array)
    		Tsys_median_list.append(Tsys_median)
    	ON_num = len([_ for _ in PTN_list_2 if _[:2]=="ON"])
		for i in range(ON_num):
    		xr_test_cut.data[np.array(xr_test_cut_cp.scantype=="ON_%s"%str(i).zfill(4))] = spe_array_list[i]
    	ON_mask = []
		for _ in np.array(xr_test_cut.scantype):
    		if "ON" in _:
        		ON_mask.append(True)
    		else:
        		ON_mask.append(False)
        xr_test_cut = xr_test_cut.isel(time=ON_mask)
        xr_test_cut["Tsys"] = (("time"), np.array([x for xs in Tsys_median_list for x in xs]))
        
        xr_test_cut.to_netcdf(path_XFFTSdata+".nc")
	else:
		print("Please check the paths. ")
