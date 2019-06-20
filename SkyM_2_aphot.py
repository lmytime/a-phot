#! /usr/bin/env python
from astropy. io import ascii
from subprocess import call
import sys,os


if stat:
    print "Compilation aborted"
    sys.exit()
else:
    print "Compilation done"

cat=raw_input("Input SkyMaker list: ")
if not cat:
    sys.exit()
if not os.path.exists(cat):
    print cat,"does not exist, bye"
    sys.exit
c=ascii.read(cat)
o=open(cat+'_aphot','w')

ps=0.1
print "Pixel scale: ",ps
ido=0
for obj in c:
    ido+=1
    if obj['col9']>obj['col6']:
        r=10.0*obj['col9']/ps # ln(100), i.e. where SB is 1/100 the peak
        ell=1.0-obj['col10']
        pa=obj['col11']
    else:
        r=10.0*obj['col6']/ps
        ell=1.0-obj['col7']
        pa=obj['col8']
    o.write("%d %f %f %f %f %f\n"%(ido,obj['col2'],obj['col3'],r,1.0-ell,pa))
o.close()

print "All done"

