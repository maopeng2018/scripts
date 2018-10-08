#!/usr/bin/python
# Created by Miaomiao Zhou, 16-04-2012
# Updated by Mao Peng, 15-01-2016
# Updated by Victoria Aguilar, 18-01-2016
import sys, string, glob, os

def load_txt(file):
	input=open(file,'r')
	text=input.read()
	input.close()
	return text


# settings

#number = sys.argv[0]
# main prog

if not os.path.exists('Ttest'):os.makedirs('Ttest')

files = glob.glob('data/*.txt')
print files
template = load_txt('template.R')
for file in files:
	outfile =  file.replace('data','Ttest').replace('txt','csv')
	print outfile
	commandfile = 'rcommand.R'
	R_command = open(commandfile,'w')
	command =  template.replace("##IN##",file).replace("##OUT##",outfile)
	R_command.write(command)
	R_command.close()
	os.system('R --slave <rcommand.R')
