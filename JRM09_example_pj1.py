

import matplotlib
import numpy as np
import sympy as sympy
import matplotlib.pyplot as plt
import sympy as sym
from scipy import special
from scipy import constants
from sympy.abc import x, y
from sympy import *
import spiceypy as spice 
import astropy
from astropy.time import Time, TimeISO, TimeDelta
from plotly.offline import download_plotlyjs, init_notebook_mode, plot
from plotly.graph_objs import *
import matplotlib.ticker as ticker
import chart_studio.plotly as py
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FixedLocator
import math
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
from matplotlib.dates import num2date as n2d
from matplotlib.dates import date2num as d2n
import matplotlib.patches as patches
from matplotlib import rc
import pandas as pd 
import os
from numpy.lib.recfunctions import stack_arrays
import datetime as dt
from matplotlib import gridspec
from math import factorial
import matplotlib.lines as mlines

#setting up order of time start with year and date, then time
class TimeYearDayTimeCustom(TimeISO):
	name = 'yday_custom'
	subfmts = (('date_hms',
		'%Y-%jT%H:%M:%S',
		'{year:d}{yday:03d}T{hour:02d}:{min:02d}:{sec:02d}'),
		('date_hm',
		'%Y-%jT%H:%M',
		'{year:d}{yday:03d}T{hour:02d}:{min:02d}'),
		('date',
		'%Y%j',
		'{year:d}{yday:03d}'),
		('time',
		'%Y-%H:%M:%S',
		'{year:d} {hour:02d}:{min:02d}:{sec:02d}'))

#setting up magnetic field
class BField(object):
	"""docstring for MagField"""

	def __init__(self,frame='PC', tres=1):
		self.frame = frame
		self.tres = tres
        
    
	def getFiles(self, tb):      
		"""docstring for getFiles"""
		fpath = './'
		fstr = {'RTN':'se', 'PC':'pc_'}[self.frame]
		rstr = '60s'.format(self.tres)
		phasedir = {'RTN':'CRUISE/SE_QL/', 'PC':''}[self.frame]
		#print(rstr)
	
# Get each doy between time bounds and make astropy Time array of doys
# Convert this to custom string of 'yyyydoy' to match file formats   
		if tb is not None:
			dt = tb[1] - tb[0]
			ndays = np.ceil(dt.jd)+1
			drange = tb[0]+TimeDelta(np.arange(ndays), format='jd')
			#print(drange)
			drange.out_subfmt = 'date'
			datestrs = drange.yday_custom
	    
    #gather the relevant files from the 1sec folder
			#self.files = ['/data/sol-ionosphere/juno/data/60sec/fgm_jno_l3_'+d+fstr+'r'+rstr+'_v01.sts' for d in datestrs]
			self.files = []
			for d in datestrs:
				if (os.path.isfile('/data/sol-ionosphere/juno/data/60sec/fgm_jno_l3_'+d+fstr+'pj01_'+'r'+rstr+'_v02.sts')) or (os.path.isfile('/data/sol-ionosphere/juno/data/60sec/fgm_jno_l3_'+d+fstr+'r'+rstr+'_v02.sts')):
					self.files.append('/data/sol-ionosphere/juno/data/60sec/fgm_jno_l3_'+d+fstr+'pj01_'+'r'+rstr+'_v02.sts')
					self.files.append('/data/sol-ionosphere/juno/data/60sec/fgm_jno_l3_'+d+fstr+'r'+rstr+'_v02.sts')
					
				#elif os.path.isfile('/data/sol-ionosphere/juno/data/60sec/fgm_jno_l3_'+d+fstr+'r'+rstr+'_v02.sts'):				
				#	self.files.append('/data/sol-ionosphere/juno/data/60sec/fgm_jno_l3_'+d+fstr+'r'+rstr+'_v02.sts')
				else:				
					self.files.append('/data/sol-ionosphere/juno/data/60sec/fgm_jno_l3_'+d+fstr+'r'+rstr+'_v01.sts')

	def readData(self):
		"""docstring for readData"""
		#Read each file and concatenate the recarrays
		nf = len(self.files)
		j = 0
		for i in np.arange(nf):
			dat, okflag = self.readFGMData(self.files[i])
			if okflag == 1:
				if j == 0:
					self.b = dat
					j += 1
				else:
					self.b = stack_arrays([self.b, dat],usemask=False, asrecarray=True)
		if j == 0:
			raise SystemExit ('No files were found for these dates. Returning...')

          
	def readFGMData(self, f):
		"""docstring for readFGMData"""
  
#Read one file
# First, check how long the header is
		try:
			with open(f, 'r') as fh:
				h = 0
				for line in fh:
					if line[:4] == '  20':
						break
					h+=1
    # Read the file and place into a recarray
			dat = np.recfromtxt(f, 
				names=['year','doy','h','m','s','ms','ddoy','x','y','z','rng','px','py','pz'],
				usemask=False, skip_header=h)
			okflag = 1
		except IOError:
			print("File doesn't exist: "+f)
			dat = 0
			okflag = 0
		return (dat, okflag)

	def getData(self, tbounds=None):
	#From the time bounds, read the appropriate files and trim to bounds
		if tbounds is not None:
			spice.furnsh('/data/sol-ionosphere/juno/naif0011.tls')
			tbet = [spice.str2et(tb) for tb in tbounds]
			tbiso = [spice.et2utc(tb, "ISOC", 1) for tb in tbet]
			spice.kclear()
			tb = Time(tbiso, format='isot', scale='utc')
		else:
			tb = tbounds   
		#print(tbiso)
		self.getFiles(tb)
		self.readData()
		nb = len(self.b)
		#nb = 10000
		# Convert decimal doy into datetime array and convert to astropy Time object 
		dtt = [dt.datetime(self.b.year[i], 1, 1) + dt.timedelta(self.b.ddoy[i] - 1) 
			for i in range(nb)]     
		t = Time(dtt, format='datetime')
		listed=np.argsort(t)
		t=t[listed]
		inx = np.where((t >= tb[0]) & (t < tb[1]))[0]        
		self.b=self.b[listed]
		self.b = self.b[inx]

		self.t = t
		self.t = self.t[inx]

		self.tjd=t.jd
		self.tjd=self.tjd[listed]
		self.tjd=self.tjd[inx]
		
		return self.b, self.t,self.tjd
def Cartmag_Sphmag(xrr, yr, zr, bxrr, byr, bzr):
	# ;first the position vectors
	xrr = xrr / 71492.
	yr = yr / 71492.
	zr = zr / 71492.
	r = (xrr**2 + yr**2 + zr**2)**0.5
    
	th = np.arccos(zr/r)
	ph = np.arctan2(yr,xrr)

	# ;then the magnetic vectors
	br = bxrr*np.sin(th)*np.cos(ph) + byr*np.sin(th)*np.sin(ph) + bzr*np.cos(th)
	btheta = bxrr*np.cos(th)*np.cos(ph) + byr*np.cos(th)*np.sin(ph) - bzr*np.sin(th)	
	bphi = byr*np.cos(ph) - bxrr*np.sin(ph)

	# ;trig functions give answers in radians, hence conversion back to degrees
	theta = np.rad2deg(th)
	phi = np.rad2deg(ph)
	
	spbmag = (br**2 + btheta**2 + bphi**2)**0.5
	#phi = phi % 360

	return (r, theta, phi, br, btheta, bphi, spbmag)



def legendre_poly(theta, max_degree):
	ang_theta = np.deg2rad(theta)
	num_of_data = len(theta)
	num_of_LenPol = 0.5 * (max_degree**2 + 3*max_degree) + 1
    
	LenPol = np.zeros((num_of_data, max_degree+1, max_degree+1))
	dLenPol = np.zeros((num_of_data, max_degree+1, max_degree+1))
	SP_pf = np.zeros((num_of_data, max_degree+1, max_degree+1))
	SP = np.zeros((num_of_data, max_degree+1, max_degree+1))
	dSP = np.zeros((num_of_data, max_degree+1, max_degree+1))
    
	#have to manually define the first 3 polynomials
	P00 = np.ones((len(ang_theta)))
	P10 = np.cos(ang_theta)
	P11 = -np.sin(ang_theta)
    
	#have to manually define the first 3 derivatives of the polynomials
	dP00 = np.zeros((len(ang_theta)))
	dP10 = -np.sin(ang_theta)
	dP11 = -np.cos(ang_theta)
    
	#initialise the matrix - put polynomials and their derivatives into their place in the arrays
	LenPol[:,0,0] = P00
	LenPol[:,1,0] = P10
	LenPol[:,1,1] = P11
    
	dLenPol[:,0,0] = dP00
	dLenPol[:,1,0] = dP10
	dLenPol[:,1,1] = dP11
    
	#Now need to calulate the normalisation components
	SP_pf[:,0,0] = 1.
	SP_pf[:,1,0] = 1.
	SP_pf[:,1,1] = ((-1)**1) * np.sqrt(2.0*factorial(0)/factorial(2))
    
	for i in range(0, max_degree+1):
		for j in range(0,i+1):
			if j == 0:
				SP_pf[:,i,j] = 1.
			else: 
				SP_pf[:,i,j] = ((-1)**j) * np.sqrt(2.0*factorial(i-j)/factorial(i+j))

    
	#calculate the polynomials and their first derivatives - supposedly like a Pascal's triangle method
	if max_degree >= 2:
		for i in range(2, max_degree+1):
			for j in range(0, i+1):
				if j <= i-2:
					LenPol[:, i, j] = (1.0/(i-j)) * (np.cos(ang_theta) * ((2*i)-1) * LenPol[:,(i-1),j] - (i+j-1) * LenPol[:,(i-2),j])
					dLenPol[:, i, j] = (1.0/(i-j)) * (((-np.sin(ang_theta) * ((2*i)-1) * LenPol[:, (i-1),j] + np.cos(ang_theta) * ((2*i)-1) * dLenPol[:,(i-1),j]) - ((i+j-1) * dLenPol[:,(i-2),j])))
				if j == i-1:
					LenPol[:,i,j] = np.cos(ang_theta) * ((2*i)-1) * LenPol[:,(i-1),(i-1)]
					dLenPol[:,i,j] = (-np.sin(ang_theta)) * ((2*i)-1) * LenPol[:, (i-1),(i-1)] + np.cos(ang_theta) * (2*i-1) * dLenPol[:,(i-1),(i-1)]
				if j == i:
					LenPol[:,i,j] = (-1.0) * ((2*i)-1) * np.sin(ang_theta) * LenPol[:,(i-1),(i-1)]
					dLenPol[:,i,j] = (-1.0) * ((2*i)-1) *(np.cos(ang_theta) * LenPol[:,(i-1),(i-1)] + np.sin(ang_theta) * dLenPol[:,(i-1),(i-1)])
    
	#end by combining the polynomials and normalisation
	SP = SP_pf * LenPol
	dSP = SP_pf * dLenPol
    
	return(LenPol, dLenPol, SP, dSP)
	
df = pd.read_table('jrm9.txt', header=None, names=['g,h','n','m','values'])
df.rename(columns=df.iloc[0])
df = df.dropna(axis='columns')
df.reset_index(inplace=True)
ghvals = df['values'].tolist()
gh = np.asarray(ghvals)

def int_field(r, theta, phi, gh, LenPol, dLenPol, SP, dSP):

	#defining array lengths and setting up angle values from degrees to radians
	len_gh = len(gh)
	l_max = np.sqrt(len_gh)-1
	max_degree = int(l_max)
	ang_theta = np.deg2rad(theta)
	ang_phi = np.deg2rad(phi)
	num_of_data = len(r)
	Br_Forward = np.zeros((num_of_data, len_gh))
	Bt_Forward = np.zeros((num_of_data, len_gh))
	Bp_Forward = np.zeros((num_of_data, len_gh))
    
	#create empty arrays for mag field values for each degree of model, r for rad, t for theta, p for phi
	#final arrays bring all values together
    
	gnmr = np.zeros((num_of_data, max_degree+1, max_degree+1))
	hnmr = np.zeros((num_of_data, max_degree+1, max_degree+1))
    
	gnmt = np.zeros((num_of_data, max_degree+1, max_degree+1))
	hnmt = np.zeros((num_of_data, max_degree+1, max_degree+1))
    
	gnmp = np.zeros((num_of_data, max_degree+1, max_degree+1))
	hnmp = np.zeros((num_of_data, max_degree+1, max_degree+1))
    
	#gnm = np.zeros((3*num_of_data, max_degree+1, max_degree+1))
	#hnm = np.zeros((3*num_of_data, max_degree+1, max_degree+1))

	#4.calculate the field for each coord dimension and degree of model
	#each of these dimensions are then put into the same array and separated later
    
	for i in range(0, max_degree+1):
    
		r_term = (1.0/r)**(i+2)
	 
		for j in range(0, i+1):

			gnmr[:,i,j] = (i+1) * r_term * SP[:,i,j] * np.cos(j*ang_phi)
			hnmr[:,i,j] = (i+1) * r_term * SP[:,i,j] * np.sin(j*ang_phi)
	    
			gnmt[:,i,j] = (-1) * r_term * dSP[:,i,j] * np.cos(j*ang_phi)
			hnmt[:,i,j] = (-1) * r_term * dSP[:,i,j] * np.sin(j*ang_phi)
	     
			gnmp[:,i,j] = ((-1) / np.sin(ang_theta)) * j * r_term * SP[:,i,j] * (-np.sin(j*ang_phi))
			hnmp[:,i,j] = ((-1) / np.sin(ang_theta)) * j * r_term * SP[:,i,j] * np.cos(j*ang_phi)
	     
    #5.sort coeffs for forward matrix
	for i in range(0, max_degree+1):
		for j in range(0, i+1):
			if j == 0:
				idx_coeff = (i-1)**2 +2*(i-1) + 1
				if i == 0:
					idx_coeff = 0
				Br_Forward[:,idx_coeff] = gnmr[:,i,j]
				Bt_Forward[:,idx_coeff] = gnmt[:,i,j]
				Bp_Forward[:,idx_coeff] = gnmp[:,i,j]
             
			else:
				idx_coeff_g = (i-1)**2 + 2*(i-1) + (2*j)
				idx_coeff_h = (i-1)**2 + 2*(i-1) + (2*j + 1)
		
				Br_Forward[:, idx_coeff_g] = gnmr[:,i,j]
				Br_Forward[:, idx_coeff_h] = hnmr[:,i,j]
				Bt_Forward[:, idx_coeff_g] = gnmt[:,i,j]
				Bt_Forward[:, idx_coeff_h] = hnmt[:,i,j]
				Bp_Forward[:, idx_coeff_g] = gnmp[:,i,j]
				Bp_Forward[:, idx_coeff_h] = hnmp[:,i,j]

	#6.Finish by putting all separate degree terms into model field
    
	brs = np.zeros((num_of_data, len_gh))
	bts = np.zeros((num_of_data, len_gh))
	bps = np.zeros((num_of_data, len_gh))

	for i in range(0, len_gh):
		brs[:,i] = gh[i] * Br_Forward[:, i]
		bts[:,i] = gh[i] * Bt_Forward[:, i]
		bps[:,i] = gh[i] * Bp_Forward[:, i]
    
	brtot = np.zeros((num_of_data))
	bttot = np.zeros((num_of_data))
	bptot = np.zeros((num_of_data))
    
    
	brtot = brs.sum(axis=1)
	bttot = bts.sum(axis=1)
	bptot = bps.sum(axis=1)
    
	bmag = (brtot**2 + bttot**2 + bptot**2)**0.5
	intbmag = bmag

	return(brtot, bttot, bptot, intbmag) 

# we are going to get positions between these two dates
utc1 = ['25 AUG 2016 12:00:00.000', '27 AUG 2016 12:51:20.000']
fgm = BField(frame='PC') #kwarg here is the frame name matching code in filename
b, bt, btjd = fgm.getData(utc1) 
max_degree = 10
#internal field values
inr, inth, inph, inbr, inbth, inbph, inbmag = Cartmag_Sphmag(b.px, b.py, b.pz, b.x, b.y, b.z)
LenPol, dLenPol, SP, dSP = legendre_poly(inth, max_degree)
intbr, intbt, intbp, intbmag = int_field(inr, inth, inph, gh, LenPol, dLenPol, SP, dSP)


