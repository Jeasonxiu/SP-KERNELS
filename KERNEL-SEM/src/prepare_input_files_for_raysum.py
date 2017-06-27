#!/usr/bin/env python
#
def write_file_geom(pos1,pos2,npos):
	fout=open('geom.in','w')
	fout.write('# Sample geometry file. Lines starting with "#" are ignored.\n')
	fout.write('# Columns: back-azimuth (deg), slowness (s/m), N-S shift (meters)\n')
	fout.write('# north), E-W shift (meters east).\n')
	fout.write('# Watch out with large shifts -- layers could pinch out, violating\n')
	fout.write('# modelling assumptions.\n')

	#input pos1,pos2,npos
	pos=pos2
	dpos=(pos1-pos2)/float(npos-1)
	for i in range(npos):
		fout.write('  0. 9.65E-5 %5.0f 0.\n' % (pos2+float(i)*dpos) )
	fout.close()

def write_file_mod(THICKNESS_LITHOS,DIP_ANGLE_DEGREES,THICKNESS_LAB):
	fout=open('mod.in','w')
	fout.write('# Sample model file. Lines starting with "#" are ignored.\n')
	fout.write('# Layers are listed from top to bottom. The bottom layer is\n')
	fout.write('# assumed to be a half-space. Interface strike and dip apply\n')
	fout.write('# to the upper interface of the layer.\n')
	fout.write('#\n')
	fout.write('# Format:\n')
	fout.write('#     Column   Contents\n')
	fout.write('#        1    Thickness (m)\n')
	fout.write('#        2    Density (kg/m^3)\n')
	fout.write('#        3    Average P-wave velocity (m/s)\n')
	fout.write('#        4    Average S-wave velocity (m/s)\n')
	fout.write('#        5    Isotropic-layer flag (1:isotropic, 0:anisotropic)\n')
	fout.write('#        6    %P anisotropy\n')
	fout.write('#        7    %S anisotropy (if 5 and 6 are zero, isotropic layer)\n')
	fout.write('#        8    Trend of fast axis (degrees)\n')
	fout.write('#        9    Plunge of fast axis (degrees)\n')
	fout.write('#       10    Interface strike (degrees)\n')
	fout.write('#       11    Interface dip (degrees)  \n')
	fout.write('#     Note that the percentages of anisotropy are peak-to-peak\n')
	fout.write('#     (the expressions used are from Farra et al. (1991))\n')
	fout.write('#\n')
	fout.write('# Layers: crust, anisotropic wedge, isotropic half-space.\n')
	fout.write('#thick rho  alph beta iso %P  %S  tr pl  st di\n')

	VP_TOP=7920.
	VS_TOP=4400.

	VP_BOT=7280.
	VS_BOT=4050.
	
	fout.write('   %d  3300 7920 4400  1  0   0    0  0   0  0\n' % (THICKNESS_LITHOS-THICKNESS_LAB/2.) )

	Z_TOP = 0.0
	Z_BOT = THICKNESS_LAB
	DELTA_Z = 1000.
	Z=0. - DELTA_Z/2.  # Interpolate to center of layer

	NLAY_TRANSITION = int(round(THICKNESS_LAB/DELTA_Z))
	for i in range(NLAY_TRANSITION):
		Z=Z+DELTA_Z
		VP=linear_interp(Z,Z_TOP,Z_BOT,VP_TOP,VP_BOT)
		VS=linear_interp(Z,Z_TOP,Z_BOT,VS_TOP,VS_BOT)
		fout.write(' %4d  3300 %4d %4d  1  0   0    0  0   90.0  %d\n' % (DELTA_Z,VP,VS,DIP_ANGLE_DEGREES) )

	fout.write('    0  3300 7280 4050  1  0   0    0  0   90.0  %d\n' % (DIP_ANGLE_DEGREES) )
	fout.close()

def write_raysum_params(SAMPLE_RATE):
	fout=open('raysum-params.in','w')
	fout.write(' # Multiples: 0 for none, 1 for Moho, 2 for all first-order, 3 to read file\n')
	fout.write('           0\n')
	fout.write(' # Number of samples per trace\n')
	fout.write('         2400\n')
	fout.write(' # Sample rate (seconds)\n')
	fout.write('   %E \n' % (SAMPLE_RATE))
	fout.write(' # Gaussian pulse width (seconds)\n')
	fout.write('   0.1    \n')  #pseudo delta function
	fout.write(' # Alignment: 0 is none, 1 aligns on P\n')
	fout.write('           1\n')
	fout.write(' # Shift of traces -- t=0 at this time (sec)\n')
	fout.write('   15.00000000    \n')
	fout.write(' # Rotation to output: 0 is NS/EW/Z, 1 is R/T/Z, 2 is P/SV/SH\n')
	fout.write('           2\n')
	fout.close()

def linear_interp(x,x0,x1,y0,y1):
	"""
	Wikipedia equation
	"""
	return y0+(y1-y0)*(x-x0)/(x1-x0)
