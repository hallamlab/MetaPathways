#!/bin/python


import sys


def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)

with open(sys.argv[1],'r') as infile:
   lines = infile.readlines()



synsamples={}
for _line in lines:
   synsample, sample = _line.strip().split('\t')

   if not synsample in synsamples:
      synsamples[synsample] = []

   synsamples[synsample].append(sample)

for synsample in synsamples:   
    printf(" -c  %s", synsample)
    for sample in synsamples[synsample]:
        printf(" -s  %s", sample)
    print 
   
   
   





