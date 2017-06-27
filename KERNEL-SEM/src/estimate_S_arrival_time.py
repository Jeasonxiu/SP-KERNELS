#!/usr/bin/env python
#
def main():
	import os,sys
	sys.path.append('src/')
	import rotate
	from numpy import argmax
	
	rotate.main() 

	os.system('ls OUTPUT_FILES/AA*BXS.semd > filelist')
	fin=open('filelist','r')
	fout=open('_station_info.txt','w')
	for line_raw in fin.readlines():
		line=line_raw.strip('\n')
		t,u = rotate.readxy(line)
		dt = t[1]-t[0]
		tS = t[argmax(u)]

		str = line[16:21] + '   %5.2f\n' %(tS)
		fout.write(str)

	fin.close()
	fout.close()


main()
