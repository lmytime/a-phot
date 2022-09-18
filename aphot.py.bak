#! /usr/bin/env python

#Copyright 2015 Emiliano Merlin
#
#This file is part of A-PHOT.
#
#A-PHOT is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#A-PHOT is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with A-PHOT.  If not, see <http://www.gnu.org/licenses/>.

import sys,os
import numpy as np
from subprocess import call
import time
from astropy.wcs import WCS
from astropy.io import ascii,fits

def wcs2pix(cat,img):

    nam=cat+'_xy.cat'
    o=open(nam,'w')
    h = fits.getheader(img)
    w = WCS(h)
    crv1,crv2=w.wcs.crval[0],w.wcs.crval[1]
    x0,y0=w.wcs_world2pix(crv1,crv2, 1)
    x1,y1=w.wcs_world2pix(crv1,crv2+0.0002777777778, 1) #this is 1 arcsec
    fac_wcs_conv=abs(np.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)))
    try:
        psi=h['CROTA2']
    except:
        psi1=np.arctan(-h['CD1_2']/h['CD2_2'])
        psi2=np.arctan(h['CD2_1']/h['CD1_1'])
        if (abs(psi1-psi2)/(psi1+1.e-8))<1.e-4:
            psi=psi1
        else:
            print "Sorry, can't properly compute the rotation angle. Aborting"
            sys.exit()
        psi*=57.2958
        
    c=ascii.read(cat)
                            
    for obj in c:

        ra,dec=obj[obj.colnames[1]],obj[obj.colnames[2]]
        x,y=w.wcs_world2pix(ra,dec, 1)
        a=fac_wcs_conv*obj[obj.colnames[3]]        
        theta=obj[obj.colnames[5]]
        if (theta<500.0):
            theta=theta-psi
        
        o.write("%d %f %f %f %f %f\n"%(obj[obj.colnames[0]],x,y,a,obj[obj.colnames[4]],theta))

    o.close()
    
    return(fac_wcs_conv,nam)

#---------------------------------------------------


def readparam(configfile,paramdict={},preservecase=False):

    f=open(configfile)
    for line in f:
        line=line.strip() #trim whitespace
        if not line.startswith('#') and len(line) > 0:
            #strip out trailing comments, if any
            try:
                ll,junk = line.split('#',1)
            except ValueError:
                ll=line
            try:
                keyword,vals=ll.split(None,1)
            except:
                print "Error in paramfile (bad entry):",ll
                print "Aborting"
                sys.exit()
            if not preservecase:
                keyword=keyword.lower()
            #Do we have multiple values?
            if len(vals.split()) == 1:
                paramdict[keyword] = truval(vals)
            else:
                paramdict[keyword] = []
                for x in vals.split():
                    paramdict[keyword].append(truval(x))

    return paramdict

#---------------------------------------------------

def truval(x):
  """ Return the properly typed value of x """

  boolvals={'true':True, 'false':False, 't': True, 'f': False,
            'yes': True, 'no': False, 'y':True, 'n': False}
  x=x.strip() #strip off leading & trailing whitespace first
  try: #is it an int?
    ans=int(x)
  except ValueError:
    try: #is it a float?
      ans = float(x)
    except ValueError: #is it a literal?
      if x.startswith("'") or x.startswith('"'):
        ans=x[1:-1]
      else:
        try: #is it a boolean?
          ans=boolvals[x.lower()]
        except: #must be a string
          ans=x
  return ans #  returns the new "properly typed" val

#-----------------------------------------#

def run(parfile,args):
    """
    Calls C core code
    """
    print "Start: ",time.time(),time.clock()
    
    d=readparam(parfile,{}) 

    n=0
    for arg in args:
        if arg.startswith('-'):
            try:
                kw=arg.strip('-')
                kw=kw.strip("'")
                kw=kw.strip('"')
                print "> Keyword",kw,"has now value",args[n+1]
                kw=kw.lower()
                d[kw]=args[n+1]
            except:
                print "Error in parsing command line arguments"
                sys.exit()
        n+=1

    apc=open('ap_circ.lst','w')
    ape=open('ap_ell.lst','w')

    if (d['wcscat']==True or d['wcscat']=='True'):
        try:
            fac_wcs_conv,incat=wcs2pix(d['input_cat'],d['image'])
        except:
            print "TERROR: Problems in converting from WCS to pixel space."
            print "Perhaps some header keyword is missing."
            print "Better to build the input catalogue in pixel space separately. Aborting"
            sys.exit()
            
        try:    
            while "," in d['aper_circ_list']:
                d1,d['aper_circ_list']=d['aper_circ_list'].split(",",1)
                apc.write("%f "%(fac_wcs_conv*float(d1)))
        except:
            pass
        d1=d['aper_circ_list']
        apc.write("%f\n"%(fac_wcs_conv*float(d1)))            
            
    else:
        incat=d['input_cat']
        try:
            apc.write("%s\n"%(d['aper_circ_list'].replace(',',' ')))
        except:
            apc.write("%s\n"%(d['aper_circ_list']))

    try:
        ape.write("%s\n"%(d['fac_ell_list'].replace(',',' ')))
    except:
        ape.write("%s\n"%(d['fac_ell_list']))
      
    apc.close()        
    ape.close()
    
    bkgdv=0
    if (d['bkgd']==True or d['bkgd']=='True'):
        bkgdv=1
    cl=0
    if (d['clip']==True or d['clip']=='True'):
        cl=1
    rpix=0
    if (d['replacepix']==True or d['replacepix']=='True'):
        rpix=1
    mJy=0
    if (d['microjy']==True or d['microjy']=='True'):
        mJy=1 
    force_kron_comput=0
    if (d['force_kron_comput']==True or d['force_kron_comput']=='True'):
        force_kron_comput=1
        
    cmd=' '.join(['aphot_core',d['image'],d['rms_image'],
                  d['seg_image'],incat,'ap_circ.lst','ap_ell.lst',
                  str(d['binsubpix']),str(d['gain']),str(bkgdv),str(cl),
                  d['output_cat'],d['output_cat']+"_morpho",
                  str(d['fmin_kron']),str(d['fmax_kron']),
                  str(d['fac_sigma']),str(d['fac_sigma2']),str(d['rbuf']),
                  str(d['fac_rbkgd']),str(d['zp']),str(d['pixscale']),
                  d['flag_image'],str(rpix),str(d['rms_tol']),str(mJy),str(force_kron_comput)])
    print cmd
    stat=call(cmd,shell=True)
    if stat != 0:
        print "A-PHOT failed, stat=",stat
        sys.exit()

    os.remove('ap_circ.lst')
    os.remove('ap_ell.lst')

    print "End: ",time.time()

#-----------------------------------------------#

def pf():
    """
    Prints template parameter file
    """
    pf=open('aphot.conf','w')
    
    pf.write("\n")
    pf.write("image             image.fits       # measurement image\n")
    pf.write("rms_image         image.rms.fits   # rms map image\n")
    pf.write("seg_image         None             # segmentation map image (optional)\n")
    pf.write("flag_image        None             # flag map image (optional)\n")
    pf.write("input_cat         image.cat_aphot  # input catalogue (aphot format)\n")
    pf.write("WCScat            False            # input catalogue in WCS [True/False]\n")
    pf.write("aper_circ_list    1.5,3.0,4.5      # list of circular apertures\n")
    pf.write("                                   # (diameters in pixels; input at least one)\n")
    pf.write("fac_ell_list      0.5,1.0          # list of Kron multiplicative factors\n")
    pf.write("                                   # for elliptical apertures; input at least one\n")
    pf.write("                                   # (input 0 for best S/N aperture)\n")
    pf.write("binsubpix         10               # subpixel binning\n")
    pf.write("gain              -1               # -1 to exclude\n")
    pf.write("bkgd              False            # background subtraction [True/False]\n")
    pf.write("clip              False            # flag pixels above fac_sigma2 [True/False]\n")
    pf.write("replacepix        False            # replace bad/contaminated pixels\n")
    pf.write("                                   #  with value of symmetric ones [True/False]\n")    
    pf.write("zp                0.0              # 0.0 to ignore\n")
    pf.write("pixscale          1.0              # arcsec; 1.0 to ignore\n")
    pf.write("force_kron_comput False          # compute only Kron radius from input morphology\n")
    pf.write("output_cat        image.aphot.out  # output catalogue\n")
    pf.write("microJy           False            # save output catalog in microJy (requires zp!=0)\n")
    pf.write("### DEFAULT PARAMETERS, do not change unless you know what you are doing :)\n")
    pf.write("fmin_kron         8.0\n")
    pf.write("fmax_kron         20.0\n")
    pf.write("fac_sigma         1.5\n")
    pf.write("fac_sigma2        3.0\n")
    pf.write("rbuf              10\n")
    pf.write("fac_rbkgd         1.2\n")
    pf.write("RMS_tol           1000.0\n")
    pf.close()


#----------------------------------------------#

if __name__ == '__main__':
    print " "
    print "*********** A-PHOT ************"
    print "**          v. 1.2           **"
    print "*******************************"

    if np.size(sys.argv)==1:
        print "Usage: aphot [paramfile] [options]"
        if os.path.exists('aphot.conf'):
            parfile='aphot.conf'
        else:
            print "No parameter file found. Aborting"
            sys.exit()          

    elif np.size(sys.argv)==2:

        parfile=sys.argv[1] # reads name of parameterfile
        if (parfile=="-h" or parfile=="h" or parfile=="H" or parfile=="-H" or ("help" in parfile)):
            print "Usage: aphot [paramfile] [options]"
            print " >> Options: to input a keyword add '-<keyword> <value>'"
            print " >> E.g.: aphot param.file -image img.fits -c_aper 1.5,6.0" 
            sys.exit()
        if (parfile=="-p"):
            print "Usage: aphot [paramfile] [options]"
            print "Print parameter file"
            pf()
            sys.exit()
        if not os.path.exists(parfile):
            print "Usage: aphot [paramfile] [options]"
            print "No parameter file "+parfile+" found. Aborting"
            sys.exit()     

    elif np.size(sys.argv)>=3:
        parfile=sys.argv[1]
        if not os.path.exists(parfile):
            print "Usage: aphot [-i] [paramfile] [options]"
            print "No parameter file "+parfile+" found. Aborting"
            sys.exit()          

    print "Parameter file:",parfile
    run(parfile,sys.argv)


        




