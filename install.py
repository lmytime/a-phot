#! /usr/bin/env python

from astropy. io import ascii
from subprocess import call
import sys,os,glob

### EDIT HERE, INSERTING THE PATH TO CFITSIO INSTALLATION FOLDER
cf="/your/path/to/cfitsio/"
###

lcf=cf+"lib"
icf=cf+"include"
nam="aphot_core"

cmd="chmod +x aphot.py"
stat=call(cmd,shell=True)
if stat:
    print "Could not chmod aphot.py"
    sys.exit()

srcs=""
srcsl=glob.glob("src/*.c")
for s in srcsl:
    srcs+=s+" "

cmd=" ".join(["gcc -o", nam, srcs, "-L", lcf,"-lcfitsio -lm -I", icf, "-I ./include"])
stat=call(cmd,shell=True)

if stat:
    print "Compilation aborted"
    sys.exit()
else:
    print "A-PHOT compiled"

if not os.path.isdir('bin'): 
    os.mkdir('bin')
cmd="cd bin ; mv ../aphot_core . ; ln -s -f ../aphot.py ./aphot ; cd .."
stat=call(cmd,shell=True)
if stat:
    print "Could not write executable in bin folder"
    sys.exit()
else:
    print "A-PHOT installed, have fun"
