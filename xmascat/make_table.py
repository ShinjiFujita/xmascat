### common modules ###
import numpy as np
import os
import sys
from astropy.utils import iers
iers.conf.auto_download = False
iers.conf.auto_max_age = None

import xarray as xr
import casacore.tables as tb
from datetime import datetime

### functions ###

################################################################################
# addKeywords
def addKeywords(colname,keyword_list,type_list,value_list,f):

	if not isinstance(keyword_list, list):
		keyword_list = [keyword_list]

	if not isinstance(type_list, list):
		type_list = [type_list]

	if not isinstance(value_list, list):
		value_list = [value_list]

	f.write('.keywords '+colname+'\n')

	for keyword,_type,value in zip(keyword_list,type_list,value_list):
		f.write(keyword+' '+_type+' '+value+'\n')

	f.write('.endkeywords'+'\n')


#def makeMAIN_ASTEXFFTS(tablename,specdata,time,state_id,texp=0.1,tBW=2.5e9):
def makeMAIN_ASTEXFFTS(tablename, xr_data, tBW=2.5e9):

	# modules
	#import casacore.tables as tb

# 	# params
# 	nrow  = specdata.shape[0]
# 	nspw  = specdata.shape[1]
# 	npol = specdata.shape[2]
# 	nchan  = specdata.shape[3]
	
	# params
	nrow  = len(xr_data["time"])
# 	nspw  = 1
# 	npol = 1
	nchan  = len(xr_data["ch"])
	
	texp = float(xr_data["integtime"][0])

#	weight = tBW/float(nchan) * texp
#	sigma = (tBW/float(nchan) * texp)**-0.5
	weight = 1.0
	sigma = 1.0

	ind_spw = (np.linspace(0,2*nrow-1,2*nrow,dtype='int32') % 2)
	#ind_spw = np.zeros(nrow)

	# tables
	colnames = ['UVW',
				'FLAG',
				'FLAG_CATEGORY',
				'WEIGHT',
				'SIGMA',
				'ANTENNA1',
				'ANTENNA2',
				'ARRAY_ID',
				'DATA_DESC_ID',
				'EXPOSURE',
				'FEED1',
				'FEED2',
				'FIELD_ID',
				'FLAG_ROW',
				'INTERVAL',
				'OBSERVATION_ID',
				'PROCESSOR_ID',
				'SCAN_NUMBER',
				'STATE_ID',
				'TIME',
				'TIME_CENTROID',
				'FLOAT_DATA'
				]

	colkeywords = [
		{'MEASINFO': {'Ref': 'ITRF', 'type': 'uvw'},'QuantumUnits': np.array(['m', 'm', 'm'],dtype='|S2')},
		{},
		{'CATEGORY': np.array([],dtype='|S1')},
		{},{},{},{},{},{},
		{'QuantumUnits': np.array(['s'],dtype='|S2')},
		{},{},{},{},
		{'QuantumUnits': np.array(['s'],dtype='|S2')},
		{},{},{},{},
		{'MEASINFO': {'Ref': 'UTC', 'type': 'epoch'}, 'QuantumUnits': np.array(['s'],dtype='|S2')},
		{'MEASINFO': {'Ref': 'UTC', 'type': 'epoch'}, 'QuantumUnits': np.array(['s'],dtype='|S2')},
		{'UNIT': 'K'}
				  ]

	ndims = [1,2,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2]
	isarrays = [True,True,True,True,True,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,True]
	valuetypes = ['double','bool','bool','float','float','int','int','int','int','double','int','int','int','bool','double','int','int','int','int','double','double','float']
#	valuetypes = ['double','bool','bool','float','float','int','int','int','int','double','int','int','int','bool','double','int','int','int','int','string','string','float']

	descs = []
	for colname,colkeyword,ndim,valuetype,isarray in zip(colnames,colkeywords,ndims,valuetypes,isarrays):
		if colname=='UVW':
			descs.append(tb.makearrcoldesc(colname,0.0,datamanagertype='StandardStMan',datamanagergroup='StandardStMan',ndim=ndim,keywords=colkeyword,valuetype=valuetype,options=5,shape=np.array([3], dtype='int32')))
		elif isarray:
			descs.append(tb.makearrcoldesc(colname,0.0,datamanagertype='StandardStMan',datamanagergroup='StandardStMan',ndim=ndim,keywords=colkeyword,valuetype=valuetype))
		else:
			descs.append(tb.makescacoldesc(colname,0.0,datamanagertype='StandardStMan',datamanagergroup='StandardStMan',keywords=colkeyword,valuetype=valuetype))

	td = tb.maketabdesc(descs=descs)

#	returned_table = tb.table(tablename,tabledesc=td,nrow=2*nrow,readonly=False)
	returned_table = tb.table(tablename, tabledesc=td, nrow=nrow, readonly=False)

	# put values
#	value = np.zeros([2*nrow,3],dtype='float64')
	value = np.zeros([nrow, 3],dtype='float64')
	returned_table.putcol('UVW',value)

#	value = np.full([2*nrow,2],sigma)
	value = np.full([nrow, 1], sigma)
	returned_table.putcol('SIGMA', value)

#	value = np.full([2*nrow,2],weight)
	value = np.full([nrow, 1], weight)
	returned_table.putcol('WEIGHT', value)

#	value = np.zeros([2*nrow],dtype='int32')
	value = np.zeros([nrow], dtype='int32')
	returned_table.putcol('ANTENNA1',value)
	returned_table.putcol('ANTENNA2',value)
	returned_table.putcol('ARRAY_ID',value)
	returned_table.putcol('FEED1',value)
	returned_table.putcol('FEED2',value)
	returned_table.putcol('FIELD_ID',value)
	returned_table.putcol('OBSERVATION_ID',value)
	returned_table.putcol('PROCESSOR_ID',value)

#	value = np.zeros([2*nrow],dtype='bool')
	value = np.zeros([nrow], dtype='bool')
	returned_table.putcol('FLAG_ROW',value)

#	value = np.full([2*nrow], texp, dtype='float64')
	value = np.array(xr_data["integtime"])
	returned_table.putcol('EXPOSURE',value)
	returned_table.putcol('INTERVAL',value)

#	value = np.zeros_like(np.concatenate([specdata[:,0],specdata[:,1]],axis=0),dtype='bool')
#	value = value.transpose(0,2,1)
	value = np.zeros([nrow, nchan, 1], dtype='bool')
	returned_table.putcol('FLAG', value)
	
	value = np.zeros([nrow, 1, nchan, 1], dtype='bool')
	returned_table.putcol('FLAG_CATEGORY', value)

#	value  = np.zeros_like(np.concatenate([specdata[:,0],specdata[:,1]],axis=0),dtype='float64')
#	value  = np.zeros((nrow, nchan), dtype='float64')

# 	for i in range(2):
# 		value[ind_spw==i] = specdata[:,i].copy()
#	value = value.transpose(0,2,1)
	
	value = np.array(xr_data["data"])
	value = value.reshape(value.shape[0], value.shape[1], 1)######
	returned_table.putcol('FLOAT_DATA',value)

# 	value = np.zeros(2*nrow,dtype='int32')
# 	for i in range(2):
# 		value[ind_spw==i] = state_id

#	value = np.zeros(2*nrow,dtype='int32')
	value = np.zeros(nrow,dtype='int32')
#	for i in range(2):
#		value[ind_spw==i] = i
	returned_table.putcol('DATA_DESC_ID', value)

#	value = np.zeros(2*nrow,dtype='int32')
	value = np.zeros(nrow,dtype='int32')
#	for i in range(2):
#		value[ind_spw==i] = state_id
	returned_table.putcol('STATE_ID',value)

#	value = np.zeros(2*nrow,dtype='int32')
	value = np.linspace(0,nrow-1,nrow,dtype='int32')
#	for i in range(2):
#		value[ind_spw==i] = np.linspace(0,nrow-1,nrow,dtype='int32')
	returned_table.putcol('SCAN_NUMBER', value)

#	value = np.zeros(2*nrow,dtype='float64')
	value = np.array(xr_data["time"], dtype="float64")/1e9 + 3506750339.0
#	value = np.array([datetime.fromtimestamp(_).timestamp() for _ in np.array(xr_data["time"], dtype="float64")/1e9]).astype("float64")
#	value = np.zeros_like(np.array(xr_data["time"], dtype="float64"), dtype="float64")
#	for i in range(2):
#		value[ind_spw==i] = time.copy()
	returned_table.putcol('TIME', value)
	returned_table.putcol('TIME_CENTROID', value)

	returned_table.flush()
	returned_table.close()




# end makeMAIN_ASTEXFFTS
################################################################################
# makeANTENNA

def makeANTENNA(tablename,
				outputfilename,
				dish_diameter = 10.,
				nbeam = 1,
				beamname = 'ASTE',
				lon = "-67d42m11.89525s",
				lat = "-22d58m17.69447s",
				height = "4861.9m",
				type = 'GROUND-BASED',
				mounit='ALT-AZ'):

	# modules
	from astropy import coordinates
	#import casacore.tables as tb

	# params
	c = coordinates.EarthLocation.from_geodetic(lon=lon,lat=lat,height=height)
	x = c.value[0]
	y = c.value[1]
	z = c.value[2]
	beamname_list = []
	for i in range(nbeam):
		beamname_list.append(beamname+str(i))

	# header
	header1 = 'OFFSET;POSITION;TYPE;DISH_DIAMETER;FLAG_ROW;MOUNT;NAME;STATION'
	header2 = 'D3;D3;A;D;B;A;A;A'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	for i in range(nbeam):
		OFFSET = '0;0;0'
		POSITION = str(x)+';'+str(y)+';'+str(z)
		TYPE = '"'+type+'"'
		DISH_DIAMETER = str(dish_diameter)
		FLAG_ROW = '0'
		MOUNT = '"'+mounit+'"'
		NAME = '"'+';'.join(beamname_list)+'"'
		STATION = '""'

		f.write(OFFSET+';')
		f.write(POSITION+';')
		f.write(TYPE+';')
		f.write(DISH_DIAMETER+';')
		f.write(FLAG_ROW+';')
		f.write(MOUNT+';')
		f.write(NAME+';')
		f.write(STATION+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['m', 'm', 'm'],dtype='|S2'),
			 'MEASINFO': {'type': 'position','Ref': 'ITRF'},
			}
	returned_table.putcolkeywords('OFFSET',value)

	value = {'QuantumUnits': np.array(['m', 'm', 'm'],dtype='|S2'),
			 'MEASINFO': {'type': 'position','Ref': 'ITRF'},
			}
	returned_table.putcolkeywords('POSITION',value)

	value = {'QuantumUnits': np.array(['m'],dtype='|S2')}
	returned_table.putcolkeywords('DISH_DIAMETER',value)

# end makeANTENNA
################################################################################
# makeANTENNA

def makeDATA_DESCRIPTION(tablename,outputfilename,nspw=1):

	# modules
	#import casacore.tables as tb

	# params

	# header
	header1 = 'FLAG_ROW;POLARIZATION_ID;SPECTRAL_WINDOW_ID'
	header2 = 'B;I;I'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	for i in range(nspw):
		FLAG_ROW = '0'
		POLARIZATION_ID = '0'
		SPECTRAL_WINDOW_ID = str(i)
		f.write(FLAG_ROW+';')
		f.write(POLARIZATION_ID+';')
		f.write(SPECTRAL_WINDOW_ID+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

# end makeDATA_DESCRIPTION
################################################################################
# makeDOPPLER

def makeDOPPLER(tablename,outputfilename):

	# modules
	#import casacore.tables as tb

	# params

	# header
	header1 = 'DOPPLER_ID;SOURCE_ID;TRANSITION_ID;VELDEF'
	header2 = 'I;I;I;D'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	f.write('')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['m/s'],dtype='|S4'),
			 'MEASINFO': {'Ref': 'RADIO', 'type': 'doppler'},
			}
	returned_table.putcolkeywords('VELDEF',value)
	
	value = np.array([-1],dtype='int32')  ###!!!!!!!!
	returned_table.putcol('DOPPLER_ID',value)  ###!!!!!!!!

# end makeDOPPLER
################################################################################
# makeFEED

def makeFEED(tablename, outputfilename, xr_data, nspw=1):

	#modules
	#import casacore.tables as tb

	# params

	# header
	header1 = 'POSITION;BEAM_OFFSET;POLARIZATION_TYPE;POL_RESPONSE;RECEPTOR_ANGLE;ANTENNA_ID;BEAM_ID;FEED_ID;INTERVAL;NUM_RECEPTORS;SPECTRAL_WINDOW_ID;TIME'
	header2 = 'D3;D2,2;A2;X2,2;D2;I;I;I;D;I;I;D'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	for i in range(nspw):
		POSITION = '0;0;0'
		BEAM_OFFSET = '0;0;0;0'
		POLARIZATION_TYPE = '"X";"Y"'
		POL_RESPONSE = '0;0;0;0;0;0;0;0'
		RECEPTOR_ANGLE = '0;0'
		ANTENNA_ID = '0'
		BEAM_ID = '0'
		FEED_ID = '0'
		INTERVAL = '0'
		#NUM_RECEPTORS = '0'
		NUM_RECEPTORS = '2'
		SPECTRAL_WINDOW_ID = str(i)
		TIME = str(np.array(xr_data["time"].astype("float64")/1e9)[0])

		f.write(POSITION+';')
		f.write(BEAM_OFFSET+';')
		f.write(POLARIZATION_TYPE+';')
		f.write(POL_RESPONSE+';')
		f.write(RECEPTOR_ANGLE+';')
		f.write(ANTENNA_ID+';')
		f.write(BEAM_ID+';')
		f.write(FEED_ID+';')
		f.write(INTERVAL+';')
		f.write(NUM_RECEPTORS+';')
		f.write(SPECTRAL_WINDOW_ID+';')
		f.write(TIME+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['m', 'm', 'm'],dtype='|S2'),
			 'MEASINFO': {'type': 'position','Ref': 'ITRF'},
			}
	returned_table.putcolkeywords('POSITION',value)

	value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
			 'MEASINFO': {'type': 'direction','Ref': 'J2000'},
			}
	returned_table.putcolkeywords('BEAM_OFFSET',value)

	value = {'QuantumUnits': np.array(['rad'],dtype='|S4'),}
	returned_table.putcolkeywords('RECEPTOR_ANGLE',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
	returned_table.putcolkeywords('INTERVAL',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
			}
	returned_table.putcolkeywords('TIME',value)

# end makeFEED
################################################################################
# makeFIELD

#def makeFIELD(tablename,outputfilename,ra,dec,fieldname,time):
def makeFIELD(tablename, outputfilename, xr_data):

	#modules
	#import casacore.tables as tb

	# params
#	 ra_rad  = ra/180.*np.pi
#	 dec_rad = dec/180.*np.pi
#	 time_mean = np.mean(time)
	ra, dec = float(str(xr_data.SRC_POS).split(",")[0]), float(str(xr_data.SRC_POS).split(",")[1])
	ra_rad  = np.deg2rad(ra)
	dec_rad = np.deg2rad(dec)
	time_mean = np.mean(np.array(xr_data["time"], dtype="float64")/1e9) + 3506750339.0

	# header
	header1 = 'DELAY_DIR;PHASE_DIR;REFERENCE_DIR;CODE;FLAG_ROW;NAME;NUM_POLY;SOURCE_ID;TIME'
	header2 = 'D2,1;D2,1;D2,1;A;B;A;I;I;D'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	DELAY_DIR = str(ra_rad)+';'+str(dec_rad)
	PHASE_DIR = str(ra_rad)+';'+str(dec_rad)
	REFERENCE_DIR = str(ra_rad)+';'+str(dec_rad)
	CODE = '""'
	FLAG_ROW = '0'
	NAME = '"'+xr_data.SRC_NAME+'"'
	NUM_POLY = '0'
	SOURCE_ID = '0'
	TIME = str(time_mean)

	f.write(DELAY_DIR+';')
	f.write(PHASE_DIR+';')
	f.write(REFERENCE_DIR+';')
	f.write(CODE+';')
	f.write(FLAG_ROW+';')
	f.write(NAME+';')
	f.write(NUM_POLY+';')
	f.write(SOURCE_ID+';')
	f.write(TIME+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
			 'MEASINFO': {'type': 'direction','Ref': 'J2000'},
			}
	returned_table.putcolkeywords('DELAY_DIR',value)

	value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
			 'MEASINFO': {'type': 'direction','Ref': 'J2000'},
			}
	returned_table.putcolkeywords('PHASE_DIR',value)

	value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
			 'MEASINFO': {'type': 'direction','Ref': 'J2000'},
			}
	returned_table.putcolkeywords('REFERENCE_DIR',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
			}
	returned_table.putcolkeywords('TIME',value)

# end makeFIELD
################################################################################
# makeFLAG_CMD

def makeFLAG_CMD(tablename,outputfilename):

	# modules
	#import casacore.tables as tb

	# params

	# header
	header1 = 'APPLIED;COMMAND;INTERVAL;LEVEL;REASON;SEVERITY;TIME;TYPE'
	header2 = 'B;A;D;I;A;I;D;A'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	f.write('')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
			}
	returned_table.putcolkeywords('TIME',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
	returned_table.putcolkeywords('INTERVAL',value)

# end makeFLAG_CMD
################################################################################
# makeFREQ_OFFSET

def makeFREQ_OFFSET(tablename, outputfilename, xr_data):

	# modules
	#import casacore.tables as tb

	# params

	# header
	header1 = 'ANTENNA1;ANTENNA2;FEED_ID;INTERVAL;OFFSET;SPECTRAL_WINDOW_ID;TIME'
	header2 = 'I;I;I;D;D;I;D'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	#f.write('')
	
	APP_PARAMS = ""
	CLI_COMMAND = ""
	APPLICATION = ""
	MESSAGE = ""
	OBSERVER = ""
	OBJECT_ID = ""
	OBSERVATION_ID = ""
	ORIGIN = ""
	TIME = str(np.array(xr_data["time"].astype("float64")/1e9)[0])
	PRIORITY = ""
	f.write(APP_PARAMS+';')
	f.write(CLI_COMMAND+';')
	f.write(APPLICATION+';')
	f.write(MESSAGE+';')
	f.write(OBSERVER+';')
	f.write(OBJECT_ID+';')
	f.write(OBSERVATION_ID+';')
	f.write(ORIGIN+';')
	f.write(TIME+';')
	f.write(PRIORITY+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
			}
	returned_table.putcolkeywords('TIME',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
	returned_table.putcolkeywords('INTERVAL',value)

	value = {'QuantumUnits': np.array(['Hz'],dtype='|S3'),
			 'MEASINFO': {'type': 'frequency','Ref': 'TOPO'},
			}
	returned_table.putcolkeywords('OFFSET',value)

# end makeFREQ_OFFSET
################################################################################
# makeHISTORY

def makeHISTORY(tablename,outputfilename):

	# modules
	#import casacore.tables as tb

	# params

	# header
	header1 = 'APP_PARAMS;CLI_COMMAND;APPLICATION;MESSAGE;OBJECT_ID;OBSERVATION_ID;ORIGIN;PRIORITY;TIME'
	header2 = 'A1;A1;A;A;I;I;A;A;D'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	f.write('')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
			}
	returned_table.putcolkeywords('TIME',value)

# end makeHISTORY
################################################################################
# makeLMT_ARRAY

# def makeLMT_ARRAY(tablename,outputfilename,nbeam=1,npol=2,nspw=2):

#	 # modules
#	 #import casacore.tables as tb

#	 # params
#	 n_array = nbeam*npol*nspw

#	 # header
#	 header1 = 'ARRAY;BEAM;POLARIZATION;SPECTRAL_WINDOW'
#	 header2 = 'I;I;I;I'

#	 f = open(outputfilename+'.header','w')

#	 f.write(header1+'\n')
#	 f.write(header2+'\n')
#	 f.close()

#	 # table
#	 f = open(outputfilename+'.dat','w')

#	 for i in range(n_array):
#		 ARRAY = str(i)
#		 BEAM = str((i // npol) % nbeam)
#		 POLARIZATION = str((i % npol) * 3 + 9)
#		 SPECTRAL_WINDOW = str((i // (npol*nbeam)) % nspw)

#		 f.write(ARRAY+';')
#		 f.write(BEAM+';')
#		 f.write(POLARIZATION+';')
#		 f.write(SPECTRAL_WINDOW+'\n')

#	 f.close()
#	 # make table
#	 returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

# end makeLMT_ARRAY
################################################################################
# makeOBSERVATION

#def makeOBSERVATION(tablename,outputfilename,time,project,observer,telescope_name='NRO'):
def makeOBSERVATION(tablename, outputfilename, xr_data, telescope_name='ASTE'):

	# modules
	#import casacore.tables as tb

	# params
	#time = np.array([datetime.fromtimestamp(_).timestamp() for _ in np.array(xr_data["time"], dtype="float64")/1e9]).astype("float64")
	time = np.array(xr_data["time"], dtype="float64")/1e9 + 3506750339.0
	time_start = time.min()
	time_end   = time.max()

	# header
	header1 = 'TIME_RANGE;LOG;SCHEDULE;FLAG_ROW;OBSERVER;PROJECT;RELEASE_DATE;SCHEDULE_TYPE;TELESCOPE_NAME'
	header2 = 'D2;A1;A1;B;A;A;D;A;A'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	TIME_RANGE = str(time_start)+';'+str(time_end)
	LOG = ''
	SCHEDULE = ''
	FLAG_ROW = '0'
	OBSERVER = '"'+xr_data.OBS_USER+'"'
	PROJECT = '"'+xr_data.PROJECT+'"'
	RELEASE_DATE = '0'
	SCHEDULE_TYPE = ''
	TELESCOPE_NAME = '"'+telescope_name+'"'

	f.write(TIME_RANGE+';')
	f.write(LOG+';')
	f.write(SCHEDULE+';')
	f.write(FLAG_ROW+';')
	f.write(OBSERVER+';')
	f.write(PROJECT+';')
	f.write(RELEASE_DATE+';')
	f.write(SCHEDULE_TYPE+';')
	f.write(TELESCOPE_NAME+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
			}
	returned_table.putcolkeywords('TIME_RANGE',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
			}
	returned_table.putcolkeywords('RELEASE_DATE',value)

# end makeOBSERVATION
################################################################################
# makePOINTING

# def makePOINTING(tablename,outputfilename,direction,time,nbeam=1,texp=0.2):

# 	'''
# 	direction: (nscan*nbeam) array
# 	'''

# 	# modules
# 	#import casacore.tables as tb

# 	# params
# 	direction_rad = direction/180.*np.pi
# 	n_direction = direction.shape[0]
# 	time_end   = time.max()

# 	# header
# 	header1 = 'DIRECTION;ANTENNA_ID;INTERVAL;NAME;NUM_POLY;TARGET;TIME;TIME_ORIGIN;TRACKING'
# 	header2 = 'D2,1;I;D;A;I;D2,1;D;D;B'

# 	f = open(outputfilename+'.header','w')

# 	f.write(header1+'\n')
# 	f.write(header2+'\n')
# 	f.close()

# 	# table
# 	f = open(outputfilename+'.dat','w')

# 	for i in range(n_direction):
# 		#DIRECTION = str(direction_rad[i][0])+';'+str(direction_rad[i][1])
# 		DIRECTION = ';'
# 		ANTENNA_ID = str((n_direction // (n_direction/nbeam))-1)
# 		INTERVAL = str(texp)
# 		NAME = ''
# 		NUM_POLY= '0'
# 		TARGET = ';'
# 		TIME = time[i].astype('|U')
# 		TIME_ORIGIN = '0'
# 		TRACKING = '0'

# 		f.write(DIRECTION+';')
# 		f.write(ANTENNA_ID+';')
# 		f.write(INTERVAL+';')
# 		f.write(NAME+';')
# 		f.write(NUM_POLY+';')
# 		f.write(TARGET+';')
# 		f.write(TIME+';')
# 		f.write(TIME_ORIGIN+';')
# 		f.write(TRACKING+'\n')

# 	f.close()

# 	# make table
# 	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

# 	value  = np.zeros([n_direction,1,2],dtype='float64')

# 	value[:,0] = direction_rad.copy()

# 	returned_table.putcol('DIRECTION',value)

# 	value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
# 			 'MEASINFO': {'type': 'direction','Ref': 'J2000'},
# 			}
# 	returned_table.putcolkeywords('DIRECTION',value)

# 	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
# 	returned_table.putcolkeywords('INTERVAL',value)

# 	value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
# 			 'MEASINFO': {'type': 'direction','Ref': 'J2000'},
# 			}
# 	returned_table.putcolkeywords('TARGET',value)

# 	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
# 			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
# 			}
# 	returned_table.putcolkeywords('TIME',value)

# 	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
# 			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
# 			}
# 	returned_table.putcolkeywords('TIME_ORIGIN',value)

# 	returned_table.flush()
# 	returned_table.close()

# end makePOINTING
################################################################################
# makePOINTING2

#def makePOINTING2(tablename,direction,time,nbeam=1,texp=0.2):
def makePOINTING(tablename, xr_data, nbeam=1):

	# modules
	#import casacore.tables as tb

	# params
	time = np.array(xr_data["time"], dtype="float64")/1e9 + 3506750339.0
	#time = np.array([datetime.fromtimestamp(_).timestamp() for _ in np.array(xr_data["time"], dtype="float64")/1e9]).astype("float64")
	texp = float(xr_data["integtime"][0])
	n_direction = len(xr_data["azimuth"])
	#direction = np.array([xr_data["azimuth"], xr_data["elevation"]]).T
	#direction_rad = direction/180.*np.pi
	time_end   = time.max()
	NAME = str(xr_data.SRC_NAME)
	TARGET_ra = np.deg2rad(np.array(xr_data["longitude"]))
	TARGET_dec = np.deg2rad(np.array(xr_data["latitude"]))
	TARGET = np.array([[ra, dec] for ra, dec in zip(TARGET_ra, TARGET_dec)])

	# make table
	colnames = ['DIRECTION',
				'ANTENNA_ID',
				'INTERVAL',
				'NAME',
				'NUM_POLY',
				'TARGET',
				'TIME',
				'TIME_ORIGIN',
				'TRACKING']
	colkeywords = [{'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4')},
				   {},
				   {'QuantumUnits': np.array(['s'],dtype='|S2')},
				   {},{},
				   {'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4')},
				   {'MEASINFO': {'Ref': 'UTC', 'type': 'epoch'}, 'QuantumUnits': np.array(['s'],dtype='|S2')},
				   {'MEASINFO': {'Ref': 'UTC', 'type': 'epoch'}, 'QuantumUnits': np.array(['s'],dtype='|S2')},
				   {}
				  ]

	#ndims = [2,1,1,1,1,2,1,1,1]
	#isarrays = [True,False,False,False,False,True,False,False,False]
	valuetypes = ['double','int','double','string','int','double','double','double','bool']
#	valuetypes = ['double','int','double','string','int','double','Date','Date','bool']
	ndims = [2,1,1,1,1,-1,1,1,1]
	descs = []

	for colname,colkeyword,valuetype,ndim in zip(colnames,colkeywords,valuetypes,ndims):
		if (colname=='DIRECTION' or colname=='TARGET'):
			descs.append(tb.makearrcoldesc(colname,0.0,datamanagertype='StandardStMan',datamanagergroup='StandardStMan',ndim=ndim,keywords=colkeyword,valuetype=valuetype,options=0))
		else:
			descs.append(tb.makescacoldesc(colname,0.0,datamanagertype='StandardStMan',datamanagergroup='StandardStMan',keywords=colkeyword,valuetype=valuetype))

	td = tb.maketabdesc(descs=descs)

	returned_table = tb.table(tablename,tabledesc=td,nrow=n_direction,readonly=False)

#	value  = np.zeros([n_direction,1,2],dtype='float64')
#	value[:,0] = direction_rad.copy()
	value  = TARGET.reshape(n_direction, 1, 2)
	returned_table.putcol('DIRECTION',value)

	value = np.zeros([n_direction],dtype='int32')
	returned_table.putcol('ANTENNA_ID',value)
	returned_table.putcol('NUM_POLY',value)

#	value = np.zeros([n_direction],dtype='float64')
	value = time
	returned_table.putcol('TIME_ORIGIN',value)

	value = np.ones([n_direction], dtype='bool')
	returned_table.putcol('TRACKING',value)

	value = np.full([n_direction],texp,dtype='float64')
	returned_table.putcol('INTERVAL',value)
	
	value = np.full([n_direction], NAME)
	returned_table.putcol('NAME', value)

	value = time
	returned_table.putcol('TIME',value)
	
	value = TARGET.reshape(-1, 1, 2)
	returned_table.putcol('TARGET',value)

	returned_table.flush()
	returned_table.close()

# end makePOINTING2
################################################################################
# makePOINTING

def makePOLARIZAION(tablename,outputfilename,npol=1):

	# modules
	#import casacore.tables as tb

	# params

	# header
	header1 = 'CORR_TYPE;CORR_PRODUCT;FLAG_ROW;NUM_CORR'
	#header2 = 'I'+str(npol)+';I'+str(npol)*2+','+str(npol)+';B;I'
	header2 = 'I'+str(npol)+';I2'+','+str(npol)+';B;I'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	CORR_TYPE = ''
	for i in range(npol):
		CORR_TYPE = CORR_TYPE + str(3*i+9) + ';'
	CORR_TYPE = CORR_TYPE[0:-1]
	#CORR_PRODUCT = '0;'*(npol**2-1)+'0'
	CORR_PRODUCT = "0;0"
	FLAG_ROW = '0'
	NUM_CORR = str(npol)

	f.write(CORR_TYPE+';')
	f.write(CORR_PRODUCT+';')
	f.write(FLAG_ROW+';')
	f.write(NUM_CORR+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

# end makePOLARIZAION
################################################################################
# makePROCESSOR

def makePROCESSOR(tablename,outputfilename):

	# modules
	#import casacore.tables as tb

	# params

	# header
	header1 = 'FLAG_ROW;MODE_ID;TYPE;TYPE_ID;SUB_TYPE'
	header2 = 'B;I;A;I;A'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	#f.write('0;-1;"";-1;""')
	f.write('0;0;"";0;""')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

# end makePROCESSOR
################################################################################
# makeSOURCE

#def makeSOURCE(tablename,outputfilename,ra,dec,sysvel,sourcename,time,freq):
def makeSOURCE(tablename, outputfilename, xr_data):

	# modules
	#import casacore.tables as tb

	# params
#	 ra_rad  = ra/180.*np.pi
#	 dec_rad = dec/180.*np.pi
#	 time_interval = time.max() - time.min()
#	 time_mean = np.mean(time)
#	 nspw = freq.shape[0]

	ra, dec = xr_data.SRC_POS.split(",")
	ra, dec = float(ra), float(dec)
	ra_rad  = np.deg2rad(ra)
	dec_rad = np.deg2rad(dec)
	time = np.array(xr_data["time"], dtype="float64")/1e9 + 3506750339.0
	#time = np.array([datetime.fromtimestamp(_).timestamp() for _ in np.array(xr_data["time"], dtype="float64")/1e9]).astype("float64")
	time_interval = time.max() - time.min()
	time_interval = float(time_interval)/1e9
	time_mean = np.mean(time)
	nspw = 1
	sysvel = float(xr_data.VELO)
	sourcename = str(xr_data.SRC_NAME)
	freq = float(xr_data.REST_FREQ)

	# header
	header1 = 'DIRECTION;PROPER_MOTION;CALIBRATION_GROUP;CODE;INTERVAL;NAME;NUM_LINES;SOURCE_ID;SPECTRAL_WINDOW_ID;TIME;TRANSITION;REST_FREQUENCY;SYSVEL'
	header2 = 'D2;D2;I;A;D;A;I;I;I;D;A1;D1;D1'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	for i in range(nspw):
		DIRECTION = str(ra_rad)+';'+str(dec_rad)
		PROPER_MOTION = '0;0'
		#CALIBRATION_GROUP ='-1'
		CALIBRATION_GROUP ='0'
		CODE = '""'
		INTERVAL = str(time_interval)
		NAME = '"'+sourcename+'"'
		NUM_LINES = '1'
		SOURCE_ID = '0'
		SPECTRAL_WINDOW_ID = str(i)
		TIME = str(time_mean)
		TRANSITION = ''
		#REST_FREQUENCY = str(np.mean(freq[i]))
		REST_FREQUENCY = str(freq)
		SYSVEL = str(sysvel)

		f.write(DIRECTION+';')
		f.write(PROPER_MOTION+';')
		f.write(CALIBRATION_GROUP+';')
		f.write(CODE+';')
		f.write(INTERVAL+';')
		f.write(NAME+';')
		f.write(NUM_LINES+';')
		f.write(SOURCE_ID+';')
		f.write(SPECTRAL_WINDOW_ID+';')
		f.write(TIME+';')
		f.write(TRANSITION+';')
		f.write(REST_FREQUENCY+';')
		f.write(SYSVEL+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
			 'MEASINFO': {'type': 'direction','Ref': 'J2000'},
			}
	returned_table.putcolkeywords('DIRECTION',value)

	value = {'QuantumUnits': np.array(['rad/s'],dtype='|S6')}
	returned_table.putcolkeywords('PROPER_MOTION',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2')}
	returned_table.putcolkeywords('INTERVAL',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
			}
	returned_table.putcolkeywords('TIME',value)


	value = {'QuantumUnits': np.array(['Hz'],dtype='|S3'),
			 'MEASINFO': {'type': 'frequency','Ref': 'TOPO'},
			}
	returned_table.putcolkeywords('REST_FREQUENCY',value)

	value = {'QuantumUnits': np.array(['m/s'],dtype='|S4'),
			 'MEASINFO': {'Ref': 'LSRK', 'type': 'radialvelocity'},
			}
	returned_table.putcolkeywords('SYSVEL',value)


# end makeSOURCE
################################################################################
# makeSPECTRAL_WINDOW

#def makeSPECTRAL_WINDOW(tablename,outputfilename,freq,nspw=2):
def makeSPECTRAL_WINDOW(tablename, outputfilename, xr_data, nspw=1, tBW=2.5e9):

	# modules
	#import casacore.tables as tb

	# params
	nchan = len(xr_data["ch"])
	freq = np.array(xr_data["freq"]).reshape(nspw, nchan)
	tBW = abs(freq[0].max() - freq[0].min())

	# header
	header1 = 'MEAS_FREQ_REF;REF_FREQUENCY;FLAG_ROW;FREQ_GROUP;FREQ_GROUP_NAME;IF_CONV_CHAIN;NAME;NET_SIDEBAND;NUM_CHAN;TOTAL_BANDWIDTH;CHAN_FREQ;CHAN_WIDTH;EFFECTIVE_BW;RESOLUTION'
	header2 = 'I;D;B;I;A;I;A;I;I;D;D'+str(nchan)+';D'+str(nchan)+';D'+str(nchan)+';D'+str(nchan)

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	for i in range(nspw):
		MEAS_FREQ_REF = '5' #TOPO   ### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#MEAS_FREQ_REF = '1' 
		#CHAN_FREQ = ';'.join(map(str,freq[i]))
		#CHAN_FREQ = ''
		REF_FREQUENCY = str(freq[i][0])
		#CHAN_WIDTH = (str(np.diff(freq[i])[0])+';')*(nchan-1)+str(np.diff(freq[i])[0])
		#CHAN_WIDTH = ''
		#EFFECTIVE_BW = (str(np.abs(np.diff(freq[i]))[0])+';')*(nchan-1)+str(np.abs(np.diff(freq[i]))[0])
		#EFFECTIVE_BW = ''
		#RESOLUTION = (str(np.abs(np.diff(freq[i]))[0])+';')*(nchan-1)+str(np.abs(np.diff(freq[i]))[0])
		#RESOLUTION = ''
		FLAG_ROW = '0'
		FREQ_GROUP = '0'
		FREQ_GROUP_NAME = ''
		IF_CONV_CHAIN = '0'
		NAME = ''
		NET_SIDEBAND = str(i)
		NUM_CHAN = str(nchan)
		TOTAL_BANDWIDTH = str(tBW)

		f.write(MEAS_FREQ_REF+';')
		#f.write(CHAN_FREQ+';')
		f.write(REF_FREQUENCY+';')
		#f.write(CHAN_WIDTH+';')
		#f.write(EFFECTIVE_BW+';')
		#f.write(RESOLUTION+';')
		f.write(FLAG_ROW+';')
		f.write(FREQ_GROUP+';')
		f.write(FREQ_GROUP_NAME+';')
		f.write(IF_CONV_CHAIN+';')
		f.write(NAME+';')
		f.write(NET_SIDEBAND+';')
		f.write(NUM_CHAN+';')
		f.write(TOTAL_BANDWIDTH+';\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = freq.copy()
	returned_table.putcol('CHAN_FREQ',value)

	value = np.zeros_like(freq)
# 	for i in range(2):
# 		value[i] = np.full(nchan,np.diff(freq[i])[0])
	value[0] = np.full(nchan, np.diff(freq[0])[0])
	returned_table.putcol('CHAN_WIDTH', value)

	value = np.zeros_like(freq)
# 	for i in range(2):
# 		value[i] = np.full(nchan,abs(np.diff(freq[i])[0]))
	value[0] = np.full(nchan,abs(np.diff(freq[0])[0]))
	
	returned_table.putcol('EFFECTIVE_BW',value)
	returned_table.putcol('RESOLUTION',value)

	value = {'MEASINFO': {'TabRefCodes': np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8, 64], dtype='uint32'),
						  'TabRefTypes': np.array(['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP','CMB', 'Undefined'],dtype='|S10'),
						  'VarRefCol': 'MEAS_FREQ_REF',
						  'type': 'frequency'},
			 'QuantumUnits': np.array(['Hz'],dtype='|S3')}
	returned_table.putcolkeywords('CHAN_FREQ',value)

	value = {'MEASINFO': {'TabRefCodes': np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8, 64], dtype='uint32'),
						  'TabRefTypes': np.array(['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP','CMB', 'Undefined'],dtype='|S10'),
						  'VarRefCol': 'MEAS_FREQ_REF',
						  'type': 'frequency'},
			 'QuantumUnits': np.array(['Hz'],dtype='|S3')}
	returned_table.putcolkeywords('REF_FREQUENCY',value)

	value = {'QuantumUnits': np.array(['Hz'],dtype='S3')}
	returned_table.putcolkeywords('CHAN_WIDTH',value)

	value = {'QuantumUnits': np.array(['Hz'],dtype='S3')}
	returned_table.putcolkeywords('EFFECTIVE_BW',value)

	value = {'QuantumUnits': np.array(['Hz'],dtype='S3')}
	returned_table.putcolkeywords('RESOLUTION',value)

	value = {'QuantumUnits': np.array(['Hz'],dtype='S3')}
	returned_table.putcolkeywords('TOTAL_BANDWIDTH',value)

	returned_table.flush()
	returned_table.close()

# end makeSPECTRAL_WINDOW
################################################################################
# makeSTATE

def makeSTATE(tablename,outputfilename):

	# modules
	#import casacore.tables as tb

	# params

	# header
	header1 = 'CAL;FLAG_ROW;LOAD;OBS_MODE;REF;SIG;SUB_SCAN'
	header2 = 'D;B;D;A;B;B;I'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	#f.write('0;0;0;"OBSERVE_TARGET#OFF_SOURCE,POSITION_SWITCH";1;0;2'+'\n')##############
	#f.write('0;0;0;"OBSERVE_TARGET#ON_SOURCE,POSITION_SWITCH";1;0;1'+'\n')##############
	f.write('0;0;0;"OBSERVE_TARGET#ON_SOURCE,POSITION_SWITCH";0;1;0'+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['K'],dtype='S2')}
	returned_table.putcolkeywords('CAL',value)

	value = {'QuantumUnits': np.array(['K'],dtype='S2')}
	returned_table.putcolkeywords('LOAD',value)

# end makeSTATE
################################################################################
# makeSYSCAL

#def makeSYSCAL(tablename,outputfilename,Tsys,Tsys_time,time,nbeam=1,npol=2,nspw=2):
def makeSYSCAL(tablename, outputfilename, xr_data, nbeam=1,npol=1,nspw=1):

	'''
	Tsys: (nspw,npol,nchan) array
	'''

	# modules
	#import casacore.tables as tb

	# params
# 	nspw = Tsys.shape[0]
# 	npol = Tsys.shape[1]
# 	nchan = Tsys.shape[2]
	n_array = nbeam*nspw
	nchan = len(xr_data["ch"])
# 	time_interval = time.max() - Tsys_time
	time = np.array(xr_data["time"], dtype="float64")/1e9 + 3506750339.0
	time = np.array([datetime.fromtimestamp(_).timestamp() for _ in np.array(xr_data["time"], dtype="float64")/1e9]).astype("float64")
	time_interval = time.max() - time.min()
	Tsys = np.array(xr_data["Tsys"])
	Tsys_time = np.median(time).astype("float64")

	# header
	header1 = 'ANTENNA_ID;FEED_ID;INTERVAL;SPECTRAL_WINDOW_ID;TIME;TSYS' #;TSYS_SPECTRUM'
	header2 = 'I;I;D;I;D;R'+str(npol) #+';R'+str(npol)+','+str(nchan)
#	header1 = 'ANTENNA_ID;FEED_ID;INTERVAL;SPECTRAL_WINDOW_ID;TIME;TCAL;TSYS' #;TSYS_SPECTRUM'
#	header2 = 'I;I;D;I;D;I;R'+str(npol) #+';R'+str(npol)+','+str(nchan)

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	for i in range(n_array):
		ANTENNA_ID = str(i % nbeam)
		FEED_ID = '0'
		INTERVAL = str(time_interval)
		SPECTRAL_WINDOW_ID = str(i // nbeam)
		TIME = str(Tsys_time)
		#TSYS_SPECTRUM = ';'.join(map(str,Tsys[i // nbeam].T.reshape(npol*nchan)))
		TSYS_SPECTRUM = ''
#		TCAL = 1
		#TSYS = ';'.join(map(str,np.median(Tsys[(i // nbeam)],axis=(1,))))
		TSYS = str(np.median(Tsys))

		f.write(ANTENNA_ID+';')
		f.write(FEED_ID+';')
		f.write(INTERVAL+';')
		f.write(SPECTRAL_WINDOW_ID+';')
		f.write(TIME+';')
#		f.write(TCAL+';')
		f.write(TSYS+'\n')
		#f.write(TSYS_SPECTRUM+';\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	#value = Tsys.copy().transpose(0,2,1)
	#returned_table.putcol('TSYS_SPECTRUM',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
			}
	returned_table.putcolkeywords('TIME',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
	returned_table.putcolkeywords('INTERVAL',value)

	#value = {'QuantumUnits': np.array(['K'],dtype='S2')}
	#returned_table.putcolkeywords('TSYS_SPECTRUM',value)

	value = {'QuantumUnits': np.array(['K'],dtype='S2')}
	returned_table.putcolkeywords('TSYS',value)

# end makeSYSCAL
################################################################################
# makeWEATHER

def makeWEATHER(tablename, outputfilename):

	# modules
	#import casacore.tables as tb

	# params

	# header
	header1 = 'ANTENNA_ID;INTERVAL;TIME;TEMPERATURE;PRESSURE;REL_HUMIDITY;WIND_SPEED;WIND_DIRECTION'
	header2 = 'I;D;D;R;R;R;R;R'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	f.write('')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
			 'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
			}
	returned_table.putcolkeywords('TIME',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
	returned_table.putcolkeywords('INTERVAL',value)

	value = {'QuantumUnits': np.array(['K'],dtype='|S2'),}
	returned_table.putcolkeywords('TEMPERATURE',value)

	value = {'QuantumUnits': np.array(['hPa'],dtype='|S4'),}
	returned_table.putcolkeywords('PRESSURE',value)

	value = {'QuantumUnits': np.array(['%'],dtype='|S2'),}
	returned_table.putcolkeywords('REL_HUMIDITY',value)

	value = {'QuantumUnits': np.array(['m/s'],dtype='|S4'),}
	returned_table.putcolkeywords('WIND_SPEED',value)

	value = {'QuantumUnits': np.array(['rad'],dtype='|S4'),}
	returned_table.putcolkeywords('WIND_DIRECTION',value)

	returned_table.flush()
	returned_table.close()

# end makeWEATHER
################################################################################