#!/usr/bin/env python
#
#
def main(**params):
	from numpy import linspace,sqrt,zeros,log10
	import matplotlib as mpl
	mpl.use('PS')
	import pylab as plt
	import sys

	include_crust = False
	include_mld = False
	
	HEIGHT=params["HEIGHT"]
	WIDTH=params["WIDTH"]

	
	SCAT_STRENGTH=params["SCAT_STRENGTH"]
	z0=HEIGHT-params["SCAT_DEPTH"]
	x0=WIDTH*0.50
	SCAT_RADIUS=params["SCAT_RADIUS"]	

	BUF=100.0
	
	ORIG_X=-BUF
	ORIG_Z=-BUF
#
	END_X=WIDTH+BUF
	END_Z=HEIGHT+BUF
#
	NX=1200
	NY=0
	NZ=600
#
	VP_MIN  =  5000.0
	VP_MAX  = 10000.0
	VS_MIN  = 2000.0
	VS_MAX  = 8000.0
	RHO_MIN = 3000.0
	RHO_MAX = 4000.0
#
	TMP=linspace(ORIG_X,END_X,num=NX)
	SPACING_X=TMP[1]-TMP[0]
	SPACING_Y=0.0
	TMP=linspace(ORIG_Z,END_Z,num=NZ)
	SPACING_Z=TMP[1]-TMP[0]
#
	fout=open('DATA/tomo_file.xyz','w')
	fout.write('%f  %f   %f   %f\n' % (ORIG_X,ORIG_Z,END_X,END_Z))
	fout.write('%f  %f\n'               % (SPACING_X,SPACING_Z))
	fout.write('%d     %d\n'           % (NX,NZ))
	fout.write('%f  %f   %f   %f   %f   %f\n' % (VP_MIN,VP_MAX,VS_MIN,VS_MAX,RHO_MIN,RHO_MAX))
#
	rho=3300.0

	Z_LAB=HEIGHT-params["DEPTH_TO_LAB"]
	Z_MOHO=HEIGHT-params["DEPTH_TO_MOHO"]
	
	Z_LAB_CRATON = HEIGHT - 200000.0
	Z_MOHO_CRATON= Z_MOHO #HEIGHT -  60000.0

	Z_MLD=Z_LAB

	VP_LITH=7920.0
	VP_ASTH=7280.0
	
	VS_LITH=4400.0
	VS_ASTH=4050.0


	VP_UPPER_LITH=VP_LITH
	VS_UPPER_LITH=VS_LITH
	VP_LOWER_LITH=(VP_LITH + VP_ASTH) / 2.0
	VS_LOWER_LITH=(VS_LITH + VS_ASTH) / 2.0

	VPMAT=zeros(NX*NZ).reshape(NX,NZ)
	VSMAT=zeros(NX*NZ).reshape(NX,NZ)
	RHOMAT=zeros(NX*NZ).reshape(NX,NZ) + rho
	
	XARR=linspace(ORIG_X,END_X,num=NX)
	ZARR=linspace(ORIG_Z,END_Z,num=NZ)
	
	DX=XARR[1]-XARR[0]
	DZ=ZARR[1]-ZARR[0]
	
	mode='hetero'
	
	ncol=2

	col = [dict() for x in range(ncol)] #list of dictionaries

	#Set Column 0 Properties
	extent_col_0 = WIDTH * 0.45
	edge={}
	col[0]["left"]=-BUF
	col[0]["right"]=extent_col_0
	col[0]["craton"]=False
	col[0]["ramp"]=False	

	#Set Column 1 Prop -- Non-Craton
	col[1]["left"]=col[0]["right"] #from Column 0
	col_width= (WIDTH-extent_col_0) + BUF 
	col[1]["right"]=col[1]["left"]+col_width
	col[1]["craton"]=False
	col[1]["ramp"]=False
	
	print '#### Column Properties'
	for icol in range(ncol):
		print icol,col[icol]["left"],col[icol]["right"]
	print '####'

	#Swap craton columns if idir = 2
	if params["ILLUMINATION_DIRECTION"] == 2:
		for icol in range(ncol):
			if col[icol]["craton"]:
				col[icol]["craton"]=False
			else:
				col[icol]["craton"]=True
			#be sure ramp stays non-cratonic
			if col[icol]["ramp"]:
				col[icol]["craton"]=False

	if mode=='homo':
		for iz,z in enumerate(ZARR):
			for ix,x in enumerate(XARR):
				vp=VP_ASTH
				vs=VS_ASTH
				fout.write('%10.2f %10.2f %7.2f %7.2f %7.2f\n' % (x,z,vp,vs,rho))
				
	else:
		for iz,z in enumerate(ZARR):
			for ix,x in enumerate(XARR):
				#Determine column
				thiscol = determine_column(x,col,ncol)
				if col[thiscol]["craton"]==False and col[thiscol]["ramp"]==False:
					if z > Z_LAB:
						vp = VP_LITH
						vs = VS_LITH
					else:
						vp = VP_ASTH
						vs = VS_ASTH
				elif col[thiscol]["craton"]==False and col[thiscol]["ramp"]==True:
					if params["ILLUMINATION_DIRECTION"] == 1:
						Z_RAMP=funx(x,x1=col[thiscol]["left"],x2=col[thiscol]["right"],y1=Z_LAB,y2=Z_LAB_CRATON)
					elif params["ILLUMINATION_DIRECTION"] == 2:
						Z_RAMP=funx(x,x1=col[thiscol]["left"],x2=col[thiscol]["right"],y2=Z_LAB,y1=Z_LAB_CRATON)
					if z > Z_RAMP:
						vp = VP_LOWER_LITH
						vs = VS_LOWER_LITH
					else:
						vp = VP_ASTH
						vs = VS_ASTH
				elif col[thiscol]["craton"]==True:
					if z > Z_LAB_CRATON:
						vp = VP_LOWER_LITH
						vs = VS_LOWER_LITH
					else:
						vp = VP_ASTH
						vs = VS_ASTH
				else:
					sys.exit("***ERROR: NOT IN A COLUMN, thiscol, x = %d   %f " % (thiscol, x))
					stop
					
				VPMAT[ix,iz]=vp
				VSMAT[ix,iz]=vs

				#overwrite as scatterer if in target range
				if point_in_range(x,z,x0,z0,SCAT_RADIUS):
					VPMAT[ix,iz]=vp*(1.+SCAT_STRENGTH)
					VSMAT[ix,iz]=vs*(1.+SCAT_STRENGTH)


		#Smooth columns of the craton
		icol_smo=[]
		for icol in range(ncol):
			if col[icol]["craton"]:
				icol_smo.append(icol)
			
		lambda_cutoff=params["VERTICAL_SMOOTHING_LAMBDA"]
		for icol in icol_smo:
			VPMAT=smooth_matrix(VPMAT,'vertical',XARR,ZARR,lambda_cutoff,xlimits=True,left_limit=col[icol]["left"],right_limit=col[icol]["right"])
			VSMAT=smooth_matrix(VSMAT,'vertical',XARR,ZARR,lambda_cutoff,xlimits=True,left_limit=col[icol]["left"],right_limit=col[icol]["right"])

		#Then apply horizontal smoothing
		lambda_cutoff=params["HORIZONTAL_SMOOTHING_LAMBDA"]
		VPMAT=smooth_matrix(VPMAT,'horizontal',XARR,ZARR,lambda_cutoff)
		VSMAT=smooth_matrix(VSMAT,'horizontal',XARR,ZARR,lambda_cutoff)	

		if include_crust:
			for iz,z in enumerate(ZARR):
				for ix,x in enumerate(XARR):
					#Determine column
					thiscol = determine_column(x,col,ncol)
					if (col[thiscol]["craton"]==False) and z>=Z_MOHO:
						VPMAT[ix,iz]=5660.0
						VSMAT[ix,iz]=3200.0
						RHOMAT[ix,iz]=2800.0
					elif (col[thiscol]["craton"]==True) and z>=Z_MOHO_CRATON:
						VPMAT[ix,iz]=5660.0
						VSMAT[ix,iz]=3200.0
						RHOMAT[ix,iz]=2800.0		

		if include_mld:
			for iz,z in enumerate(ZARR):
				for ix,x in enumerate(XARR):
					#Determine column
					thiscol = determine_column(x,col,ncol)
					#if ((col[thiscol]["craton"]==True) or (col[thiscol]["ramp"]==True)) and z<Z_MOHO_CRATON and z>Z_MLD:
					if z<Z_MOHO_CRATON and z>Z_MLD:
						VPMAT[ix,iz]=VP_UPPER_LITH
						VSMAT[ix,iz]=VS_UPPER_LITH
						RHOMAT[ix,iz]=2800.0		
		
	#write to file
		for iz,z in enumerate(ZARR):
			for ix,x in enumerate(XARR):
				vp=VPMAT[ix,iz]
				vs=VSMAT[ix,iz]
				rho=RHOMAT[ix,iz]
				fout.write('%10.2f %10.2f %7.2f %7.2f %7.2f\n' % (x,z,vp,vs,rho))

	extent=[END_X/1000.0,ORIG_X/1000.0,END_Z/1000.0,ORIG_Z/1000.0]
	plt.subplot(2,1,1)
	plt.imshow(VPMAT.T,origin='lower',cmap='rainbow_r',extent=extent,aspect='equal',interpolation='nearest')
	cbar1=plt.colorbar()
	plt.subplot(2,1,2)
	plt.imshow(VSMAT.T,origin='lower',cmap='rainbow_r',extent=extent,aspect='equal',interpolation='nearest')
	cbar2=plt.colorbar()
	
	cbar1.set_label(r'$\alpha$ (m/s) ')
	cbar2.set_label(r'$\beta$ (m/s) ')
	
	plt.xlabel('Distance (km)')
	plt.ylabel('Depth (km)')
	
	plt.savefig('tomo.eps')
	
def funx(x,x1,x2,y1,y2):
	"""
	"""
	y=y1+(y2-y1)*(x-x1)/(x2-x1)
	return y
	
def determine_column(x,col,ncol):
	#Determine column
	for icol in range(ncol):
		if x >= col[icol]["left"] and x <= col[icol]["right"]:
			thiscol=icol
			break
		else:
			thiscol=999
			
	return thiscol
	
def smooth_matrix(M,dir,XARR,ZARR,lambda_cutoff,xlimits=False,left_limit=0.0,right_limit=0.0):
	from numpy import shape
	
	nx,nz=shape(M)
	
	if dir == 'horizontal':
		dx=XARR[1]-XARR[0]
		for iz,z in enumerate(ZARR):
			tmp=M[:,iz]
			M[:,iz]=smooth(tmp,dx,lambda_cutoff)
	else:
		dz=ZARR[1]-ZARR[0]
		for ix,x in enumerate(XARR):
			if (xlimits and x>=left_limit and x<=right_limit) or (xlimits==False):
				tmp=M[ix,:]
				M[ix,:]=smooth(tmp,dz,lambda_cutoff)
		
	return M


def smooth(x,dx,lambda_cutoff,window='hamming'):
	from scipy.signal import firwin, convolve
	from numpy import mean, concatenate,zeros
	
	val1=x[0]
	val2=x[-1]
	
	NFIL=101
	
	a=zeros(NFIL)+val1
	c=zeros(NFIL)+val2
	
	cutoff= 2.0 * dx / lambda_cutoff

	if cutoff >= 1.0:
		return x
	
	#if (NFIL*DX < DX/cutoff):
	#	print "***WARNING: Smoothing away structure smaller than %f with a filter of length %f:" % (2.*DX/cutoff,NFIL*DX)
	
	h=firwin(NFIL,cutoff)
	
	x1=convolve(concatenate([a,x,c]),h,mode='same')
	
	x0=x1[NFIL:-NFIL]

	return x0

def point_in_range(x,z,x0,z0,SCAT_RADIUS):
	"""
	"""
	dd=(x-x0)**2 + (z-z0)**2
	if dd > SCAT_RADIUS**2:
		return False
	else:
		return True


