def main():
	import matplotlib.gridspec as gridspec
	import pylab as plt
	import os
	
	fig = plt.figure(1,figsize=(17.5,15))
	
	gs = gridspec.GridSpec(3, 2)
	plt.subplot(gs[0, 0])
	plt.subplot(gs[1, 0])
	plt.subplot(gs[2, 0])
	plt.subplot(gs[0, 1])
	plt.subplot(gs[1, 1])
	plt.subplot(gs[2, 1])
	
	HORIZ=[40000,80000,160000]
	#VERT=[50000,100000,200000]
	DIR = [1,2]	

	list_of_axes=fig.axes
	
	count=0
	for D in DIR:
		for H in HORIZ:
			os.chdir('OUTPUT_FILES_'+str(H)+'-'+str(D))
			ax=list_of_axes[count]
			make_sub(ax,'BXP',title=str(D)+' direction',ylabel=str(H/1000.)+' km',scale_factor=90.0)
			os.chdir('..')
			count=count+1 #advance to next subplot
				
	#for ax in list_of_axes:
	#	ax.set_yticklabels([])
		
	#for ax in list_of_axes[3:9]:
	#	ax.set_ylabel('')
		
	#for ax in list_of_axes[12:]:
	#	ax.set_ylabel('')
		
	#for ax in list_of_axes[1::3]:
	#	ax.set_title('')
	#	ax.set_xlabel('')
		
	#for ax in list_of_axes[2::3]:
	#	ax.set_title('')
		
		
	#for ax in list_of_axes[::3]:
	#	ax.set_xlabel('')
		
	
	plt.savefig('mypost.eps')

def make_sub(ax,channel,title='',ylabel='',scale_factor=1.0):
	import pylab as plt
	from numpy import array
	from scipy.signal import correlate
	from numpy import argmax

	#from scipy.signal import firwin, convolve

	#dt=0.04
	#nyquist=1./2./dt

	#h=firwin(501,0.1,nyq=nyquist)

	fin1=open('_station_info.txt','r')
	
	x=0 #init
	seismos=[]

	for line in fin1.readlines():
		nfo=line.strip('\n').split()
	
		sta_name = nfo[0]
		tt_S = float(nfo[1])
		
		if sta_name == 'S0100':
			print 'Skipping '+ sta_name
			continue

		fname='OUTPUT_FILES/AA.'+sta_name+'.'+channel+'.semd'
	
		fin2=open(fname,'r')
		u=[]
		t=[]
		for line2 in fin2.readlines():
			nfo=line2.strip('\n').split()
			t.append(float(nfo[0]) - tt_S )
			u.append(float(nfo[1]))

		x=x+1
		u=array(u)
		#uf=convolve(u,h,mode='same')
		u=u*scale_factor + x		

		trace=[t,u]
		seismos.append(trace)
	
	trace_master = seismos[0]
	
	for trace in seismos:
		t = trace[0]
		u = trace[1]
	
		ax.plot(t,u,color='black')

	ax.set_title(title)
	ax.set_ylabel(ylabel)
	ax.set_xlabel('Time after $S$ (s)')

	ax.set_xlim([-45.,3.0])
	#plt.ylim([-1,32])

main()
