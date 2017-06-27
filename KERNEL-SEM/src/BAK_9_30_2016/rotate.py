#!/usr/bin/env python
#
def main():
	import os
	os.system('ls OUTPUT_FILES/AA.*.BXZ.semd > filelistZ')
	os.system('ls OUTPUT_FILES/AA.*.BXX.semd > filelistX')

	#os.system('ln -s /Users/mancinelli/PROG/ROTATE/xz2svp.py src/')

	fin=open('filelistX','r')
	filenames_X = fin.readlines()
	fin.close()

	fin=open('filelistZ','r')
	filenames_Z = fin.readlines()
	fin.close()

	for i in range(len(filenames_Z)):
		print filenames_X[i], filenames_Z[i]
		rotate(filenames_X[i].strip('\n'),filenames_Z[i].strip('\n'))

def rotate(filename_X, filename_Z, debug=False):
	"""
	"""
	import sys
	sys.path.append('/users/nmancine/PROG/ROTATE/')
	from xz2svp import xz2svp
	from numpy import arctan,tan,sin, pi
	from numpy import mean

	t,Z=readxy(filename_Z)
	t,R=readxy(filename_X)

	if (debug):
		print "DEBUGGING"

	vs_surf = 4400.
	vp_surf = 7920.
	
	vs_half = 4400. 
	vp_half = 7920.

	v_half=[vp_half,vs_half]

	#1 for P, 2 for S
	iphase=2

	incAngle = 23.0
	incAngle_rad = incAngle*pi/180.
	p = sin(incAngle_rad)/v_half[iphase-1]

	P,S = xz2svp(R,Z,vp_surf,vs_surf,p)

	tmp=filename_Z[16:21]
	writexy(t,S,'OUTPUT_FILES/AA.'+tmp+'.BXS.semd')
	writexy(t,P,'OUTPUT_FILES/AA.'+tmp+'.BXP.semd')

	if debug:
		import pylab as plt
		plt.figure(1,figsize=(12,6))
		plt.subplot(211)
		plt.plot(t,Z,color='black')
		plt.plot(t,R,color='red')
		plt.plot(t,P,'--',color='blue')
		plt.subplot(212)
		plt.plot(t,Z,color='black')
		plt.plot(t,R,color='red')
		plt.plot(t,S,'--',color='green')	

		plt.figure(2)
		plt.scatter(R,Z)
		linex=[-1,1]
		liney=[-1./tan(incAngle_rad),1./tan(incAngle_rad)]
		plt.plot(linex,liney)
		plt.ylim([-1,1])
		plt.xlim([-1,1])

		rise_over_runs=[]
		for i in range(len(R)):
			if R[i] >= 0.1:
				tmp=Z[i]/R[i]
				rise_over_runs.append(tmp)

		

		slope=mean(rise_over_runs)

		linex=[-1,1]
		liney=[-slope,slope]
		plt.plot(linex,liney,'--',color='black')

		print arctan(slope)*180./pi

		plt.show()

def get_slope(X,Y):
	sum=0.
	isum=0
	for i in range(len(X)):
		if X[i] != 0:
			isum=1+isum
			sum=Y[i]/X[i]+sum

	avg=sum/float(isum)

	return avg
		

def writexy(x,y,fname):
	
	fout=open(fname,'w')
	for i in range(len(x)):
		fout.write('%8.3f   %E\n' % (x[i],y[i]))
	fout.close()

def readxy(fname):
	from numpy import array
	fin=open(fname,'r')

	x,y=[],[]
	
	for line in fin.readlines():
		nfo=line.strip('\n').split()		
		tmp1=float(nfo[0])
		tmp2=float(nfo[1])

		x.append(tmp1)
		y.append(tmp2)

	fin.close()


	return array(x),array(y)

