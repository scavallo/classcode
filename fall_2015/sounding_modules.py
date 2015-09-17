#!/usr/bin/python
''' sounding_modules

Computes thermodynamics diagnostics from soundings

Steven Cavallo
November 2014
'''

import weather_modules as wm
import numpy as np
from scipy import ndimage
from mstats import *
from numpy import trapz
from numpy import interp
import os, datetime





Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Rs_v=461.51           # Specific gas const for water vapour, J kg^{-1} K^{-1}
Cp_da=1004.6          # Specific heat at constant pressure for dry air
Cv_da=719.            # Specific heat at constant volume for dry air
Cp_v=1870.            # Specific heat at constant pressure for water vapour
Cv_v=1410.            # Specific heat at constant volume for water vapour
Cp_lw=4218	      # Specific heat at constant pressure for liquid water
Epsilon=0.622         # Epsilon=Rs_da/Rs_v; The ratio of the gas constants
degCtoK=273.15        # Temperature offset between K and C (deg C)
rho_w=1000.           # Liquid Water density kg m^{-3}
grav=9.81             # Gravity, m s^{-2}
Lv=2.5e6              # Latent Heat of vaporisation 
boltzmann=5.67e-8     # Stefan-Boltzmann constant
mv=18.0153            # Mean molar mass of water vapor(g/mol)


def get_cape(temp,pres,dewpt,hght,startp,startt,startdp,totalcape=False):
    """Wrapper for the numerics of calculating CAPE.

    INPUTS:                                                            
    temp,pres,dewpt,hght: Vertical profile of temperature (K), pressure (Pa), 
                     and Dewpoint temperature (K), heights (m)
    startp: Starting pressure level to lift parcel from (Pa)		     
    totalcape [=False]   : Flag defining method of identifying the so-
                           called "Equilibrium Level" (Reference).
                           If False  (default), use the first stable 
                           layer above the LFC, and ignore any CAPE in 
                           unstable layers above this. If True, use all
                           CAPE up to the highest equilibrium level.

    OUTPUTS:                                                           
    P_lcl                : The lifted condensation level (LCL)
    P_lfc                : The level of free convection (LFC). Can be
                           the same as the LCL, or can be NaN if there
                           are no unstable layers.
    P_el                 : The Equilibrium Level, used to determine the
                           CAPE. If totalcape=True, use the highest 
                           equilibrium level, otherwise use the first 
                           stable equilibrium level above the LFC.
    CAPE                 : CAPE calculated from virtual temperature
    CIN                  : CIN calculated from virtual temperature

    HINT:                     
    parcel=S.get_parcel('mu') 
    lcl,lfc,el,cape,cin=get_cape(*parcel)
    """    

    # Check units
    # Init temp is startt in C, Init dew point is stwrtdp,
    # pressure levels are in hPa   
    temp = temp - 273.15 # convert temperature to celsius
    dewpt = dewpt - 273.15 # convert dewpoint to celsius
    pres = pres/100 # convert pressure to hPa
    
    
    inds = np.where( (pres < startp) ) 
    tmp = pres[inds]
    del pres
    #pres = tmp[::-1]
    pres = tmp[:]
    del tmp    
    startp = startp/100
    
    tmp = temp[inds]
    del temp
    #temp = tmp[::-1]
    temp = tmp[:]
    del tmp    

    tmp = dewpt[inds]
    del dewpt
    #dewpt = tmp[::-1]
    dewpt = tmp[:]
    del tmp    

    tmp = hght[inds]
    del hght
    #hght = tmp[::-1]
    hght = tmp[:]
    del tmp    

    
    # Get Sub-LCL traces         
    presdry,tempdry,tempiso=dry_ascent(startp,startt-degCtoK,startdp-degCtoK)    
        

    # make lcl variables explicit
    P_lcl=presdry[-1]
    T_lcl=tempdry[-1]

    # Now lift a wet parcel from the intersection point
    # preswet=linspace(P_lcl,100,101)
    preswet,tempwet=moist_ascent(P_lcl,T_lcl)

    # tparcel is the concatenation of tempdry and 
    # tempwet, and so on.
    
    tparcel=np.concatenate((tempdry,tempwet[1:]))
    pparcel=np.concatenate((presdry,preswet[1:]))

    # Interpolating the environmental profile onto the 
    # parcel pressure coordinate
    # tempenv=interp(preswet,pres[::-1],temp[::-1])
    ## NEW, for total column:
    tempenv=interp(pparcel,pres[::-1],temp[::-1])


    # now solve for the equlibrium levels above LCL
    # (all of them, including unstable ones)
    # eqlev,stab=solve_eq(preswet[::-1],(tempwet-tempenv)[::-1])
    # NEW, for total column:
    # On second thought, we don't really want/need
    # any equilibrium levels below LCL
    # eqlev,stab=solve_eq(pparcel[::-1],(tparcel-tempenv)[::-1])
    # This is equivalent to the old statement :
    eqlev,stab=solve_eq(pparcel[pparcel<=P_lcl][::-1],\
            (tparcel-tempenv)[pparcel<=P_lcl][::-1])

    aa = tparcel-tempenv

    # Sorting index by decreasing pressure
    I=np.argsort(eqlev)[::-1]
    eqlev=eqlev[I]; stab=stab[I]

    # temperatures at the equilibrium level
    # tempeq=interp(eqlev,preswet[::-1],tempenv[::-1])
    ## NEW, for total column:
    tempeq=interp(eqlev,pparcel[::-1],tparcel[::-1])

    # This helps with debugging
    # for ii,eq in enumerate(eqlev):
        # print "%5.2f  %5.2f  %2d"%(eq,tempeq[ii],stab[ii])

    # need environmental temperature at LCL
    tenv_lcl=interp(P_lcl,pparcel[::-1],tempenv[::-1])

    isstab=np.where(stab==1.,True,False)
    unstab=np.where(stab==1.,False,True)        

    if eqlev.shape[0]==0:
        # no unstable layers in entire profile
        # because the parcel never crosses the tenv
        P_lfc=float('NaN')
        P_el=float('NaN')
    elif T_lcl>tenv_lcl:
        # check LCL to see if this is unstable
        P_lfc=P_lcl
        if totalcape:
            P_el=eqlev[isstab][-1]
        else:
            P_el=eqlev[isstab][0]
    elif eqlev.shape[0]>1:
        # Parcel is stable at LCL so LFC is the 
        # first unstable equilibrium level and 
        # "EQ" level is the first stable equilibrium 
        # level
        P_lfc=eqlev[unstab][0]
        if totalcape:
            P_el=eqlev[isstab][-1]
        else:
            P_el=eqlev[isstab][0]
    else:
        # catch a problem... if there is only
        # one eqlev and it's stable (this is 
        # unphysical), then it could be a vertical
        # resolution thing. This is a kind of 
        # "null" option
        try:
	    P_el=eqlev[isstab][0]
            P_lfc=eqlev[isstab][0]
        except:
	    P_el=eqlev[unstab][0]
            P_lfc=eqlev[unstab][0]	
	
    if np.isnan(P_lfc):
        return P_lcl,P_lfc,P_el,0,0

    # need to handle case where dwpt is not available 
    # above a certain level for any reason. Most simplest 
    # thing to do is set it to a reasonably low value; 
    # this should be a conservative approach!
    
    #dwpt=dewpt.copy().soften_mask()
    [inds] = np.where(np.isnan(dewpt))
    dwpt = dewpt
    dwpt[inds] = dwpt.min()
    
    # raise ValueError
    #if dwpt[(pres>=P_el).data*(pres<P_lfc).data].mask.any():
    #    print "WARNING: substituting dwpt.min() for masked values of DWPT in this sounding"
    #dwpt[dwpt.mask]=dwpt.min()
    # dwptenv=interp(preswet,pres[::-1],dwpt[::-1])
    # NEW:

    dwptenv=interp(pparcel,pres[::-1],dwpt[::-1])


    
    #if hght[(pres>=P_el).data].mask.any():
    #    raise NotImplementedError, "TODO: Implement standard atmosphere to substitute missing heights"
    # hghtenv=interp(preswet,pres[::-1],self.soundingdata['hght'][::-1])
    # NEW:
    hghtenv=interp(pparcel,pres[::-1],hght[::-1])
    

    # Areas of POSITIVE Bouyancy
    # cond1=(tempwet>=tempenv)*(preswet<=P_lfc)*(preswet>P_el)
    # NEW:
    cond1=(tparcel>=tempenv)*(pparcel<=P_lfc)*(pparcel>P_el)
    # Areas of NEGATIVE Bouyancy
    # cond2=(tempwet<tempenv)*(preswet<=P_lcl)*(preswet>P_el)
    # NEW:
    if totalcape:
        cond2=(tparcel<tempenv)*(pparcel>P_el)
    else:
        cond2=(tparcel<tempenv)*(pparcel>P_lfc)
    # Do CAPE calculation
    # 1. Virtual temperature of parcel... remember it's saturated above LCL.
    # e_parcel=SatVap(tempwet)
    # Tv_parcel=VirtualTemp(tempwet+degCtoK,preswet*100.,e_parcel)
    # e_env=SatVap(dwptenv)
    # Tv_env=VirtualTemp(tempenv+degCtoK,preswet*100.,e_env)
    # NEW:
    e_parcel=SatVap(tparcel)
    Tv_parcel=VirtualTemp(tparcel+degCtoK,pparcel*100.,e_parcel)
    e_env=SatVap(dwptenv)
    Tv_env=VirtualTemp(tempenv+degCtoK,pparcel*100.,e_env)

    CAPE=trapz(9.81*(Tv_parcel[cond1]-Tv_env[cond1])/Tv_env[cond1],hghtenv[cond1])
    CIN=trapz(9.81*(Tv_parcel[cond2]-Tv_env[cond2])/Tv_env[cond2],hghtenv[cond2])

    return P_lcl,P_lfc,P_el,CAPE,CIN

def dry_ascent(startp,startt,startdp):
    from numpy import interp
    #--------------------------------------------------------------------
    # Lift a parcel dry adiabatically from startp to LCL.
    # Init temp is startt in C, Init dew point is stwrtdp,
    # pressure levels are in hPa    
    #--------------------------------------------------------------------

    assert startdp<=startt

    if startdp==startt:
        return np.array([startp]),np.array([startt]),np.array([startdp]),

    Pres=np.linspace(startp,600)

    # Lift the dry parcel
    T_dry=(startt+degCtoK)*(Pres/startp)**(Rs_da/Cp_da)-degCtoK 

    # Mixing ratio isopleth
    starte=SatVap(startdp)
    startw=MixRatio(starte,startp*100)
    e=Pres*startw/(.622+startw)
    T_iso=243.5/(17.67/np.log(e/6.112)-1)

    # Solve for the intersection of these lines (LCL).
    # interp requires the x argument (argument 2)
    # to be ascending in order!
    P_lcl=interp(0,T_iso-T_dry,Pres)
    T_lcl=interp(P_lcl,Pres[::-1],T_dry[::-1])

    presdry=np.linspace(startp,P_lcl)
    tempdry=interp(presdry,Pres[::-1],T_dry[::-1])
    tempiso=interp(presdry,Pres[::-1],T_iso[::-1])


    return presdry,tempdry,tempiso

def moist_ascent(startp,startt,ptop=100):
    #--------------------------------------------------------------------
    # Lift a parcel moist adiabatically from startp to endp.
    # Init temp is startt in C, pressure levels are in hPa    
    #--------------------------------------------------------------------
    preswet=np.linspace(startp,ptop,101)
    temp=startt
    tempwet=np.zeros(preswet.shape);tempwet[0]=startt
    for ii in range(preswet.shape[0]-1):
        delp=preswet[ii]-preswet[ii+1]
        temp=temp-100*delp*GammaW(temp+degCtoK,(preswet[ii]-delp/2)*100)
        tempwet[ii+1]=temp

    return preswet,tempwet

def ctok(t):
    '''
    Convert temperature from Celsius to Kelvin
    Parameters
    ----------
    t : number, numpy array
    The temperature in Celsius
    Returns
    -------
    Temperature in Kelvin (number or numpy array)
    '''
    return t + 273.15

def SatVap(tempc,phase="liquid"):
    """Calculate saturation vapour pressure over liquid water and/or ice.

    INPUTS: 
    tempc: (C)
    phase: ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
    return saturation vapour pressure as follows:

    Tc>=0: es = es_liquid
    Tc <0: es = es_ice

   
    RETURNS: e_sat  (Pa)
    
    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)
    
    This formulation is chosen because of its appealing simplicity, 
    but it performs very well with respect to the reference forms
    at temperatures above -40 C. At some point I'll implement Goff-Gratch
    (from the same resource).
    """

    over_liquid=6.112*np.exp(17.67*tempc/(tempc+243.12))*100.
    over_ice=6.112*np.exp(22.46*tempc/(tempc+272.62))*100.
    # return where(tempc<0,over_ice,over_liquid)

    if phase=="liquid":
        # return 6.112*exp(17.67*tempc/(tempc+243.12))*100.
        return over_liquid
    elif phase=="ice":
        # return 6.112*exp(22.46*tempc/(tempc+272.62))*100.
        return np.where(tempc<0,over_ice,over_liquid)
    else:
        raise NotImplementedError

def MixRatio(e,p):
    """Mixing ratio of water vapour
    INPUTS
    e (Pa) Water vapor pressure
    p (Pa) Ambient pressure
          
    RETURNS
    qv (kg kg^-1) Water vapor mixing ratio`
    """

    return Epsilon*e/(p-e)

def MixR2VaporPress(qv,p):
    """Return Vapor Pressure given Mixing Ratio and Pressure
    INPUTS
    qv (kg kg^-1) Water vapor mixing ratio`
    p (Pa) Ambient pressure
          
    RETURNS
    e (Pa) Water vapor pressure
    """

    return qv*p/(Epsilon+qv)

def VaporPressure(dwpt):
    """Water vapor pressure
    INPUTS
    dwpt (C) Dew Point Temperature (for SATURATION vapor 
	     pressure use tempc)
          
    RETURNS
    e (Pa) Water Vapor Pressure

    SOURCE:
    Bolton, Monthly Weather Review, 1980, p 1047, eq. (10)
    """

    return 611.2*exp(17.67*dwpt/(243.5+dwpt))

def DewPoint(e):
    """ Use Bolton's (1980, MWR, p1047) formulae to find tdew.
    INPUTS:
    e (Pa) Water Vapor Pressure
    OUTPUTS:
    Td (C) 
      """

    ln_ratio=np.log(e/611.2)
    Td=((17.67-ln_ratio)*degCtoK+243.5*ln_ratio)/(17.67-ln_ratio)
    return Td-degCtoK

def GammaW(tempk,pres,e=None):
    """Function to calculate the moist adiabatic lapse rate (deg C/Pa) based
    on the temperature, pressure, and rh of the environment.

    INPUTS:
    tempk (K)
    pres (Pa)
    RH (%)

    RETURNS:
    GammaW: The moist adiabatic lapse rate (Dec C/Pa)
    """

    tempc=tempk-degCtoK
    es=SatVap(tempc)
    ws=MixRatio(es,pres)

    if e is None:
        # assume saturated
        e=es

    w=MixRatio(e,pres)

    tempv=VirtualTempFromMixR(tempk,w)
    latent=Latentc(tempc)

    A=1.0+latent*ws/(Rs_da*tempk)
    B=1.0+Epsilon*latent*latent*ws/(Cp_da*Rs_da*tempk*tempk)
    Rho=pres/(Rs_da*tempv)
    Gamma=(A/B)/(Cp_da*Rho)
    return Gamma
def VirtualTempFromMixR(tempk,mixr):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    mixr: Mixing Ratio (kg/kg)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia). This is an approximation
    based on a m
    """

    return tempk*(1.0+0.6*mixr)

def Latentc(tempc):
    """Latent heat of condensation (vapourisation)

    INPUTS:
    tempc (C)

    OUTPUTS:
    L_w (J/kg)

    SOURCE:
    http://en.wikipedia.org/wiki/Latent_heat#Latent_heat_for_condensation_of_water
    """
   
    return 1000*(2500.8 - 2.36*tempc + 0.0016*tempc**2 - 0.00006*tempc**3)
def solve_eq(preswet,func):
    """Solve the peicewise-linear stability of a parcel

    INPUTS: variables from the most ascent of a parcel
    preswet: pressure
    func   : piecewise linear function to solve (tw-te)

    OUTPUTS:
    solutions: zeros of the function (tw-te)
    stability: indication of the stability of this solution.

    NOTE ABOUT STABILITY
    Stability is the sign of (d(func)/dP). So if you have used tw-te
    like you were supposed to, d(tw-te)/dP>0 means this is a stbale 
    equilibrium level (flip the sign to envision d(tw-te)/dz).
    """

    from numpy import sign,diff

    # Sorry to be annoying but I'm going to force you to use
    # a monotonically increasing variable
    #assert (sign(diff(preswet))==1).all(), "Use a monotonically increasing abscissa"

    # Identify changes in sign of function
    dsign=sign(func)
    isdiff=np.zeros(dsign.shape,dtype=bool)
    isdiff[1:]=abs(diff(dsign)).astype(bool)

    # shift to get the value on the other side
    # of the x-axis
    shift=np.zeros(dsign.shape,dtype=bool)
    shift[:-1]=isdiff[1:]; shift[-1]=isdiff[0]

    # solve by linear interpolation between 
    # values points
    sols=np.zeros((isdiff.sum()))
    stab=np.zeros((isdiff.sum()))
    for ii in range(isdiff.sum()):
        f0=func[isdiff][ii]
        f1=func[shift][ii]
        p0=preswet[isdiff][ii]
        p1=preswet[shift][ii]
        slope=(f1-f0)/(p1-p0)
        sols[ii]=p0-f0/slope
        stab[ii]=sign(slope)

    ### Debug with plots
    #fig=plt.figure()
    #ax=fig.add_subplot(111)
    #ax.plot(preswet,func)
    #ax.plot(sols,zeros(sols.shape),ls='',marker='o')
    #ax.plot(preswet[isdiff],func[isdiff],ls='',marker='+',mew=2)
    #ax.plot(preswet[shift],func[shift],ls='',marker='x',mew=2)
    #ax.grid(True)
    #show()

    return sols,stab
def VirtualTemp(tempk,pres,e):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    e: vapour pressure (Pa)
    p: static pressure (Pa)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia)."""

    tempvk=tempk/(1-(e/pres)*(1-Epsilon))
    return tempvk
    
