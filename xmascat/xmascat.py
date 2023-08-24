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
			str_tempo = str(line[3]).replace("'", "").replace("(", "").replace(")", "")
			SET_dict[line[2]] = str_tempo
			print(line[2], ": ", str_tempo)
			continue
		else:
			continue
	SET_dict = del_dict_key(SET_dict, ["DUMMY_MODE", "INTEG_TIME", "CALB_INT"])
	
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

	frame = "RADEC" #!!!!!!
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

def read_XFFTSdata(filename, PTN_list, nchan=32768, obsmode="OTF"):
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
		if obsmode=="OTF":
			dict_tempo["timestamp"] = datetime.datetime.utcfromtimestamp(timestamp_tempo).strftime("%Y-%m-%dT%H:%M:%S.%f%Z000")
		else:
			dict_tempo["timestamp"] = datetime.datetime.utcfromtimestamp(timestamp_tempo).strftime("%Y-%m-%dT%H:%M:%S.%f")
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
def create_XFFTSxarray(path_startfile=None, path_antlogfile=None, path_XFFTSdata=None, Tamb=280.0, nchan=32768, tBW=2.5e9):
	if path_startfile==None:
		print("Please specify the path_startfile. ")
		return
	elif os.path.exists(path_startfile)==False:
		print("Please check the path_startfile. ")
		return
	if path_XFFTSdata==None:
		print("Please specify the path_XFFTSdata. ")
		return
	elif os.path.exists(path_XFFTSdata)==False:
		print("Please check the path_XFFTSdata. ")
		return
	SET_dict, PTN_list = read_startfile(path_startfile)
	if SET_dict["OTF_MODE"] == "ON":
		obsmode = "OTF"
		if path_antlogfile==None:
			print("Please specify the path_antlogfile. ")
			return
		elif os.path.exists(path_antlogfile)==False:
		    print("Please check the path_antlogfile. ")
		    return
		xr_antlog = read_antlogfile(path_antlogfile)
		timestamp_xffts_list, integtime_xffts_list, scantype_xffts_list, data_xffts_list = read_XFFTSdata(path_XFFTSdata, PTN_list, nchan=nchan, obsmode=obsmode)
		#xr_data = xr_antlog.sel(time=timestamp_xffts_list, method="nearest") #### 今後分光データに座標を紐づけるようにする。
		#xr_data["ch"] = [i for i in range(nchan)]
		#xr_data["data"] = (("time", "ch"), data_xffts_list)
		xr_data = xr.Dataset(coords={"time":[datetime.datetime.fromisoformat(_[:-3]) for _ in timestamp_xffts_list], "ch":[i for i in range(nchan)]})
		xr_data["data"] = (("time", "ch"), data_xffts_list)
		#xr_data["integtime"] = (("time"), integtime_xffts_list)
		#xr_data["scantype"] = (("time"), scantype_xffts_list)
		#xr_antlog.to_netcdf(path_antlogfile+".nc")  ### !!!!!!!
		#xr_data.to_netcdf(path_XFFTSdata+".nc")  ### !!!!!!! 
		xr_antlog = xr_antlog.sel(time=~xr_antlog.get_index("time").duplicated())  ### !!!!!!! 
		xr_antlog_interp = xr_antlog.interp(time=xr_data.time)
		xr_data["azimuth"] = (("time"), np.array(xr_antlog_interp.azimuth))
		xr_data["elevation"] = (("time"), np.array(xr_antlog_interp.elevation))
		xr_data["longitude"] = (("time"), np.array(xr_antlog_interp.longitude))
		xr_data["latitude"] = (("time"), np.array(xr_antlog_interp.latitude))
		xr_data.attrs["frame"] = xr_antlog_interp.attrs["frame"]
	else:
		obsmode = "PS"
		PTN_list_ = [PTN_list[0]]
		num_dict = {}
		num_dict["R"], num_dict["SKY"], num_dict["ON"], num_dict["OFF"] = 1, 0, 0, 0
		for _ in PTN_list[1:]:
			PTN_before = PTN_list_[-1].split("_")[0]
			PTN_now = _.split("_")[0]
			num_now = num_dict[PTN_now]
			if PTN_before!=PTN_now:
				PTN_list_.append(PTN_now+"_"+str(num_now).zfill(4))
				num_dict[PTN_now] += 1
			else:
				pass
		PTN_list_ = PTN_list
		timestamp_xffts_list, integtime_xffts_list, scantype_xffts_list, data_xffts_list = read_XFFTSdata(path_XFFTSdata, PTN_list, nchan=nchan, obsmode=obsmode)
		xr_data = xr.Dataset(coords={"time":[datetime.datetime.fromisoformat(_) for _ in timestamp_xffts_list], "ch":[i for i in range(nchan)]})
		xr_data["data"] = (("time", "ch"), data_xffts_list)
	xr_data["integtime"] = (("time"), integtime_xffts_list)
	xr_data["scantype"] = (("time"), scantype_xffts_list)
		
	A_num = int(path_XFFTSdata[-2:])
	for k in SET_dict.keys():
		if k in ["REF_NUM", "RX_NAME", "REST_FREQ", "OBS_FREQ", "SIDBD_TYP", "FREQ_IF1"]:
			xr_data.attrs[k] = SET_dict[k].split(",")[A_num-1]
		else:
			xr_data.attrs[k] = SET_dict[k]
		
	# xr_data.VELO(m/s), xr_data.REST_FREQ(Hz)
	freq_offset = float(xr_data.VELO)/299792458.0*float(xr_data.REST_FREQ)
	print("freq_offset: ", freq_offset/1e9, " GHz")
	
	if xr_data.attrs["RX_NAME"] == "CAT8W": #######   !!!!!!!!!!!!!!!!!!!!!!!!
		if A_num==1 or A_num==3:
			xr_data["freq"] = (("ch"), np.linspace(float(xr_data.REST_FREQ) - tBW/2.0, float(xr_data.REST_FREQ) + tBW/2.0, num=nchan) - freq_offset)
		elif A_num==2 or A_num==4:
			xr_data["freq"] = (("ch"), np.linspace(float(xr_data.REST_FREQ) - tBW/2.0, float(xr_data.REST_FREQ) + tBW/2.0, num=nchan) - freq_offset)
	else:
		if xr_data.SIDBD_TYP=="USB":
			xr_data["freq"] = (("ch"), np.linspace(float(xr_data.REST_FREQ) - tBW/2.0, float(xr_data.REST_FREQ) + tBW/2.0, num=nchan) - freq_offset)
		elif xr_data.SIDBD_TYP=="LSB":
			xr_data["freq"] = (("ch"), np.linspace(float(xr_data.REST_FREQ) + tBW/2.0, float(xr_data.REST_FREQ) - tBW/2.0, num=nchan) - freq_offset)
		else:
			print("SIDBD_TYP is invalid. ")
	"""
	if xr_data.SIDBD_TYP=="USB":
		xr_data["freq"] = (("ch"), np.linspace(float(xr_data.REST_FREQ) - tBW/2.0, float(xr_data.REST_FREQ) + tBW/2.0, num=nchan) - freq_offset)
	elif xr_data.SIDBD_TYP=="LSB":
		xr_data["freq"] = (("ch"), np.linspace(float(xr_data.REST_FREQ) + tBW/2.0, float(xr_data.REST_FREQ) - tBW/2.0, num=nchan) - freq_offset)
	else:
		print("SIDBD_TYP is invalid. ")
	"""
	PTN_list = [_ for _ in PTN_list if sum(np.array(xr_data.scantype==_))!=0]
			
	R_list = [_ for _ in PTN_list if _[:1]=="R"]
	SKY_list = [_ for _ in PTN_list if _[:3]=="SKY"]
	ON_list = [_ for _ in PTN_list if _[:2]=="ON"]
	OFF_list = [_ for _ in PTN_list if _[:3]=="OFF"]
	meantime_R_list = [np.mean((np.array(xr_data["time"])[xr_data.scantype==_]).astype("float64")).astype('datetime64[ns]') for _ in R_list]
	meantime_SKY_list = [np.mean((np.array(xr_data["time"])[xr_data.scantype==_]).astype("float64")).astype('datetime64[ns]') for _ in SKY_list]
	meantime_ON_list = [np.mean((np.array(xr_data["time"])[xr_data.scantype==_]).astype("float64")).astype('datetime64[ns]') for _ in ON_list]
	meantime_OFF_list = [np.mean((np.array(xr_data["time"])[xr_data.scantype==_]).astype("float64")).astype('datetime64[ns]') for _ in OFF_list]

	spe_array_list = []
	Tsys_median_list = []
	for i in range(len(ON_list)):
		print("A%s"%str(A_num).zfill(2), ":", i+1, "/",  len(ON_list))
		ON = ON_list[i]
		time_ON = meantime_ON_list[i]
		R = R_list[np.argmin(abs(meantime_R_list - meantime_ON_list[i]))]
		SKY = SKY_list[np.argmin(abs(meantime_SKY_list - meantime_ON_list[i]))]
		OFF = OFF_list[np.argmin(abs(meantime_OFF_list - meantime_ON_list[i]))]
		ave_R = np.mean(np.array(xr_data["data"])[xr_data.scantype==R], axis=0)
		ave_SKY = np.mean(np.array(xr_data["data"])[xr_data.scantype==SKY], axis=0)
		ave_OFF = np.mean(np.array(xr_data["data"])[xr_data.scantype==OFF], axis=0)
		Y = ave_R/ave_SKY
		Tsys = Tamb/(Y-1.0)	
		raw_ON_array = np.array(xr_data["data"])[xr_data.scantype==ON]
		spe_array = np.array([Tamb*(raw_ON - ave_OFF)/(ave_R - ave_OFF) for raw_ON in raw_ON_array])
		spe_array_list.append(spe_array)	
		Tsys_median = [np.median(Tsys)]*len(raw_ON_array)
		Tsys_median_list.append(Tsys_median)
	ON_num = len([_ for _ in PTN_list if _[:2]=="ON"])
	for i in range(ON_num):
		xr_data.data[np.array(xr_data.scantype=="ON_%s"%str(i).zfill(4))] = spe_array_list[i]
	ON_mask = []
	for _ in np.array(xr_data.scantype):
		if "ON" in _:
			ON_mask.append(True)
		else:
	 		ON_mask.append(False)
	xr_data = xr_data.isel(time=ON_mask)
	xr_data["Tsys"] = (("time"), np.array([x for xs in Tsys_median_list for x in xs]))
	
	xr_data.to_netcdf(path_XFFTSdata+".nc")
	print("saved: ", path_XFFTSdata+".nc")
	return

def plotnc(path_xr, xmin=None, xmax=None, ymin=None, ymax=None):
	import matplotlib.pyplot as plt
	import os
	if path_xr[-2:]!="nc":
		print("Name of input Xarray file must be 'xxxxxx.nc'. ")
		return
	path_temp = path_xr[:-3]+"_plot"
	if not os.path.exists(path_temp):
		os.system('mkdir -p '+path_temp)
	print('mkdir -p '+path_temp)
	xr_data = xr.load_dataset(path_xr)
	scantype_array = np.unique(np.array(xr_data.scantype).astype("str"))
	x = np.array(xr_data.freq)
	np.save(os.path.join(path_temp, "freq.npy"), x)
	spe_list = []
	for scantype in scantype_array:
		y = np.nanmean(xr_data.data[xr_data.scantype==scantype], axis=0)
		plt.figure(figsize=(12, 6), facecolor="w")
		plt.plot(x/1e9, y)
		plt.xlabel("Frequency [GHz]")
		plt.ylabel("Ta* [K]")
		plt.grid()
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.title(scantype)
		plt.rcParams['figure.facecolor'] = "w"
		plt.rcParams['figure.edgecolor'] = "w"
		plt.savefig(os.path.join(path_temp, scantype+".png"))
		plt.clf()
		np.save(os.path.join(path_temp, scantype+".npy"), y)
		spe_list.append(y)
	spe_list = np.array(spe_list)
	np.save(os.path.join(path_temp, os.path.basename(path_xr)[:-3]+".spe_all_ave.npy"), np.nanmean(spe_list, axis=0))


def delete_scans(path_xr, delete_list):
	if path_xr[-2:]!="nc":
		print("Name of input Xarray file must be 'xxxxxx.nc'. ")
		return
	xr_data = xr.load_dataset(path_xr)
	for del_scan in delete_list:
		print("ON_%s"%str(del_scan).zfill(4))
		xr_data = xr_data.where(xr_data.scantype!="ON_%s"%str(del_scan).zfill(4), drop=True)
	xr_data["freq"] = (("ch"), np.array(xr_data.freq[:,0]))
	xr_data.to_netcdf(path_xr)
	return
	
def slice_chs(path_xr, freq_center, ch_pm, add_name):
	import numpy as np
	if path_xr[-2:]!="nc":
		print("Name of input Xarray file must be 'xxxxxx.nc'. ")
		return
	xr_data = xr.load_dataset(path_xr)
	ch_center = np.argmin(np.abs(np.array(xr_data.freq) - freq_center))
	print(ch_center-ch_pm, "ch : ", ch_center+ch_pm, "ch")
	xr_data = xr_data.sel(ch=slice(ch_center-ch_pm, ch_center+ch_pm))
	xr_data.to_netcdf(path_xr[:-3]+".%s.nc"%add_name)
	return
	



def Xarray2MS2(path_xr, removetemp=True):
	import casacore.tables as tb
	import xmascat.make_table
	if path_xr[-2:]!="nc":
		print("Name of input Xarray file must be 'xxxxxx.nc'. ")
		return
	xr_data = xr.load_dataset(path_xr)
	
	path_temp = os.path.dirname(path_xr)+"/temp/"
	MS2name = path_xr[:-3]+".ms"
	
	os.system('rm -rf '+path_temp)

	os.system('rm -rf '+MS2name)
	os.system('mkdir -p '+path_temp)
	
	A_num = int(path_xr[-5:-3])
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeMAIN_ASTEXFFTS")
	xmascat.make_table.makeMAIN_ASTEXFFTS(MS2name, xr_data, tBW=2.5e9)
	print("A%s"%str(A_num).zfill(2), ":", "END	:makeMAIN_ASTEXFFTS ")
	print()
	print("START: makeANTENNA")
	xmascat.make_table.makeANTENNA(MS2name+'/ANTENNA', path_temp+'ANTENNA')
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeANTENNA")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeDATA_DESCRIPTION")
	xmascat.make_table.makeDATA_DESCRIPTION(MS2name+'/DATA_DESCRIPTION', path_temp+'DATA_DESCRIPTION')
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeDATA_DESCRIPTION")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeDOPPLER")
	xmascat.make_table.makeDOPPLER(MS2name+'/DOPPLER', path_temp+'DOPPLER')
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeDOPPLER")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeFEED")
	xmascat.make_table.makeFEED(MS2name+'/FEED', path_temp+'FEED', xr_data)
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeFEED")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeFIELD")
	xmascat.make_table.makeFIELD(MS2name+'/FIELD', path_temp+'FIELD', xr_data)
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeFIELD")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeFLAG_CMD")
	xmascat.make_table.makeFLAG_CMD(MS2name+'/FLAG_CMD', path_temp+'FLAG_CMD')
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeFLAG_CMD")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeFREQ_OFFSET")
	xmascat.make_table.makeFREQ_OFFSET(MS2name+'/FREQ_OFFSET', path_temp+'FREQ_OFFSET', xr_data)
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeFREQ_OFFSET")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeHISTORY")
	xmascat.make_table.makeHISTORY(MS2name+'/HISTORY', path_temp+'HISTORY')
	print("END	: makeHISTORY")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeOBSERVATION")
	xmascat.make_table.makeOBSERVATION(MS2name+'/OBSERVATION', path_temp+'OBSERVATION', xr_data)
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeOBSERVATION")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makePOINTING")
	xmascat.make_table.makePOINTING(MS2name+'/POINTING', xr_data)
	print("A%s"%str(A_num).zfill(2), ":", "END	: makePOINTING")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makePOLARIZAION")
	xmascat.make_table.makePOLARIZAION(MS2name+'/POLARIZATION', path_temp+'POLARIZATION')
	print("A%s"%str(A_num).zfill(2), ":", "END	: makePOLARIZAION")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makePROCESSOR")
	xmascat.make_table.makePROCESSOR(MS2name+'/PROCESSOR', path_temp+'PROCESSOR')
	print("A%s"%str(A_num).zfill(2), ":", "END	: makePROCESSOR")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeSOURCE")
	xmascat.make_table.makeSOURCE(MS2name+'/SOURCE', path_temp+'SOURCE', xr_data)
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeSOURCE")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeSPECTRAL_WINDOW")
	xmascat.make_table.makeSPECTRAL_WINDOW(MS2name+'/SPECTRAL_WINDOW',path_temp+'SPECTRAL_WINDOW', xr_data)
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeSPECTRAL_WINDOW")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeSTATE")
	xmascat.make_table.makeSTATE(MS2name+'/STATE', path_temp+'STATE')
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeSTATE")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeSYSCAL")
	xmascat.make_table.makeSYSCAL(MS2name+'/SYSCAL', path_temp+'SYSCAL', xr_data)
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeSYSCAL")
	print()
	print("A%s"%str(A_num).zfill(2), ":", "START: makeWEATHER")
	xmascat.make_table.makeWEATHER(MS2name+'/WEATHER',path_temp+'WEATHER' )
	print("A%s"%str(A_num).zfill(2), ":", "END	: makeWEATHER")
	print()

	abs_path = os.path.abspath(MS2name)

	keywords = {'MS_VERSION': 2.0,
				'ANTENNA': 'Table: '+abs_path+'/ANTENNA',
				'DATA_DESCRIPTION': 'Table: '+abs_path+'/DATA_DESCRIPTION',
				'DOPPLER': 'Table: '+abs_path+'/DOPPLER',
				'FEED': 'Table: '+abs_path+'/FEED',
				'FIELD': 'Table: '+abs_path+'/FIELD',
				'FLAG_CMD': 'Table: '+abs_path+'/FLAG_CMD',
				'FREQ_OFFSET': 'Table: '+abs_path+'/FREQ_OFFSET',
				'HISTORY': 'Table: '+abs_path+'/HISTORY',
				'OBSERVATION': 'Table: '+abs_path+'/OBSERVATION',
				'POINTING': 'Table: '+abs_path+'/POINTING',
				'POLARIZATION': 'Table: '+abs_path+'/POLARIZATION',
				'PROCESSOR': 'Table: '+abs_path+'/PROCESSOR',
				'SOURCE': 'Table: '+abs_path+'/SOURCE',
				'SPECTRAL_WINDOW': 'Table: '+abs_path+'/SPECTRAL_WINDOW',
				'STATE': 'Table: '+abs_path+'/STATE',
				'SYSCAL': 'Table: '+abs_path+'/SYSCAL',
				'WEATHER': 'Table: '+abs_path+'/WEATHER',
				}
	returnedTable = tb.table(MS2name, readonly=False)
	returnedTable.putkeywords(keywords)

	returnedTable.flush(recursive=True)
	returnedTable.close()

	if removetemp==True:
		os.system('rm -rf '+path_temp)
	print("saved: ", abs_path)
	return