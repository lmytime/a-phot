#! /usr/bin/env python

from __future__ import print_function
from astropy. io import ascii
from subprocess import call
import sys,os
import numpy as np

print ("USAGE: ./SEx_2_aphot.py [INPUT_CAT [options]]")
print ("OPTIONS:") 
print (" -m : compute morphological parameters internally")
print (" -k <K> : change default value of multiplicative Kron factor to K (see manual)")
print (" -s <SM> : difference FWHM LRI-HRI, in pixels (will smooth lengths and radii accordingly)")


try:
    cat=sys.argv[1]
except:
    cat=raw_input("Input SEx catalog: ")
    if not cat:
        print ("Meh. Exiting")
        sys.exit()
    if not os.path.exists(cat):
        print (cat,"does not exist, bye")
        sys.exit

compute_morpho='n'
k=2.5
sm=0.0

for narg in range(np.size(sys.argv)):
    if sys.argv[narg][0]=="-":
        if sys.argv[narg][1]=="m":
            compute_morpho='y'
            print ("Computing morphological parameters internally")
        if sys.argv[narg][1]=="k":
            k=float(sys.argv[narg+1])
            print ("Changing Kron multiplicative factor")
        if sys.argv[narg][1]=="s":
            sm=float(sys.argv[narg+1])
            print ("Smoothing the lengths")
            
print ("k=",k)
print ("sm=",sm)

c=ascii.read(cat)
o=open(cat+'_aphot','w')
try:
    obj=c[0]
    a=obj['NUMBER']
    a=obj['X_IMAGE']
    a=obj['Y_IMAGE']
    a=obj['XMIN_IMAGE']
    a=obj['YMIN_IMAGE']
    a=obj['XMAX_IMAGE']
    a=obj['YMAX_IMAGE']     
    a=obj['KRON_RADIUS']
    a=obj['ELONGATION']
    a=obj['THETA_IMAGE']
except:
    print ("Something is wrong with the SExtractor catalogue / cannot find columns, please edit by hand")
    sys.exit()
    
for obj in c:
    
    try:
        ell=obj['ELLIPTICITY']
    except:
        ell=1.0-1.0/obj['ELONGATION']
        
    theta=obj['THETA_IMAGE']

    if compute_morpho=='y':        
        a=max(obj['X_IMAGE']-obj['XMIN_IMAGE'],obj['XMAX_IMAGE']-obj['X_IMAGE'])
        b=max(obj['Y_IMAGE']-obj['YMIN_IMAGE'],obj['YMAX_IMAGE']-obj['Y_IMAGE'])
        a=-(3.0*1.414213562*(max(a,b)+2))
        ell=0.0
        theta=501.0
    else:
        a=k*obj['KRON_RADIUS']*obj['A_IMAGE']/2.5
        
    if sm>0.0:
        # "smooth" to LRI FWHM
        if compute_morpho=='y':
            a=-(-a+sm)
        else:
            b=a*(1.0-ell)
            a=a+sm
            b=b+sm
            ell=1.0-b/a
        
    o.write("%d %f %f %f %f %f\n"%(obj['NUMBER'],obj['X_IMAGE'],obj['Y_IMAGE'],a,ell,theta))
    
o.close()

print ("All done")

