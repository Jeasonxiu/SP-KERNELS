#!/usr/bin/env python
#
def main():
	draw_row(1,4,'10000')
	draw_row(2,4,'20000')
	draw_row(3,4)
	draw_row(4,4)
	plt.savefig('mypost.eps')

def draw_row(irow,nrow):
	import pylab as plt
	make_sub(nrow,2,1,'BXS',tsfile='tsS.xy',title='S',scale_factor=1.0)
	make_sub(nrow,2,2,'BXP',tsfile='tsP.xy',title='Sp',scale_factor=30.0)

	for trace_number in range(1,31):
		tmp=str(trace_number).zfill(2)
		zoff=float(trace_number)

		add_raysum_seis_to_plot(nrow,2,1,'raysum_sv.'+tmp+'.xy',scale_factor=1.0,zoff=zoff)
		add_raysum_seis_to_plot(nrow,2,2,'raysum_p.'+tmp+'.xy',scale_factor=30.0,zoff=zoff)

