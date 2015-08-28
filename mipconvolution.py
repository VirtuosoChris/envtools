#!/usr/bin/env python

import argparse
import math
import os.path
from subprocess import call

def latlongResToPhongExponent (res, err = 1e-5):
	return (res * res) / (-8 * math.log(err))

if __name__=="__main__":
	parser=argparse.ArgumentParser(description="Generate mipmap convolution chains")

	parser.add_argument('--projection_format','-p',choices=['cube','dome','hemi','ball','latlong','llsquare'],help='The spherical projection used in the image input/output',default='latlong')
	parser.add_argument('--filter','-f',choices=['linear','nearest'],help='The interpolation filter for image operations',default='linear')
	parser.add_argument('--output_size','-s',default=256,help='The size of the output image',type=int)
	parser.add_argument('--sample_pattern','-sp',choices=['cent','rgss','box2','box3','box4'],default='rgss',help='The subsampling pattern for oversampling')
	parser.add_argument('--diffuse_size','-d',type=int,help='The resolution of the diffuse map (if specified)')
	parser.add_argument('--levels','-n',type=int,help='The number of mip levels to generate')
	parser.add_argument('--output_file','-o',default='result.hdr', help='The output file pattern (e.g. foo.hdr will become foo_0.hdr foo_1.hdr etc.)')
	parser.add_argument('--input_files', type = str, nargs=6, help='The input files, in the case you use 6 inputs for a cube')
	parser.add_argument('--input_file', type = str, help='The input file')

	args=parser.parse_args()
	if(not args.levels):
		args.levels=int(math.log(args.output_size,2))+1

	ofse=os.path.splitext(args.output_file)
	ofpattern=ofse[0]+'_%02d'+ofse[1]

	outargs=[]
	outargs.extend(['-i',args.projection_format])
	outargs.extend(['-p',args.sample_pattern])
	outargs.extend(['-f',args.filter])
	for level in range(args.levels):
		thisoarg=[]
		thisoarg.extend(['-o',ofpattern % (level)])

		dims = args.output_size >> level

		thisoarg.extend(['-n',str(dims)])

		if args.projection_format == 'cube':
			exponent = latlongResToPhongExponent(dims*2) #2 cube faces cover hemisphere from pole to pole, equivalent to height of latlong
		else:
			exponent = latlongResToPhongExponent(dims)

		print "Exponent for level " + str(level) + " : " + str(exponent)

		thisoarg.extend(['-s',str(exponent)])

		outargs.extend(thisoarg)

	if(args.diffuse_size):
		outargs.extend(['-o',ofse[0]+'_diffuse'+ofse[1]])
		outargs.extend(['-d','5'])
		outargs.extend(['-n',str(args.diffuse_size)])

	if args.input_files :
		outargs.extend(args.input_files)
	else:
		outargs.append(args.input_file)

	toolbin=os.path.join(os.path.dirname(os.path.realpath(__file__)),'envconv')
	call([toolbin]+outargs)
	#print([toolbin]+outargs)
	print("done")
