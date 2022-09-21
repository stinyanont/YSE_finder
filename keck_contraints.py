import astropy.units as u
from astropy import coordinates
import sys, argparse
import numpy as np
import os
import time
import pandas as pd
from astropy.time import Time,TimeDelta
from datetime import datetime, timedelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun
from astroplan import Observer
import pdb
import math
# The figure below shows the telescope limits for the Keck telescopes. 
# In the regions where an elevation limit of ~35° is shown, 
# the telescope hits the Nasmyth deck. Note that this is a problem in the South and West for Keck II
# while for Keck I it is more of a problem in the North and East.

def keck1_rising_time(sky_coord, date_UT):
    date_UT = Time(date_UT)
    midnight = date_UT - utc_offset
    # start_time = date_UT - utc_offset - 6*u.hour #6pm
    # end_time =    date_UT - utc_offset + 6*u.hour #6am
    delta_midnight = np.linspace(-6, 6, 720*2)*u.hour #one per 0.5 minute
    observer_frame = AltAz(obstime=midnight+delta_midnight,
                              location=keck_location)

    obj_altaz = sky_coord.transform_to(observer_frame)
    rising = obj_altaz.az < 180*u.deg #select only times where the object is setting. 
    #now check for Az when Alt is either 

    nasmyth_platform = obj_altaz[rising].alt > 33.3*u.deg #time steps where object could be below the nasmyth platform, K1
    general_shutter  = obj_altaz[rising].alt > 18*u.deg

    #Object never rises
    if np.sum(general_shutter) == 0:
        rise_time = None
    else: #object do rise
        if obj_altaz[rising][nasmyth_platform][0].az >5.3*u.deg and obj_altaz[rising][nasmyth_platform][0].az < 146.2*u.deg:
            rise_time = obj_altaz[rising][nasmyth_platform][0].obstime
        else:
            rise_time = obj_altaz[rising][general_shutter][0].obstime
    return rise_time

def keck2_setting_time(sky_coord, date_UT):
    date_UT = Time(date_UT)
    midnight = date_UT - utc_offset
    # start_time = date_UT - utc_offset - 6*u.hour #6pm
    # end_time =    date_UT - utc_offset + 6*u.hour #6am
    delta_midnight = np.linspace(-6, 6, 720*2)*u.hour #one per 0.5 minute
    observer_frame = AltAz(obstime=midnight+delta_midnight,
                              location=keck_location)

    obj_altaz = sky_coord.transform_to(observer_frame)
    setting = obj_altaz.az > 180*u.deg #select only times where the object is setting. 
    #now check for Az when Alt is either 

    nasmyth_platform = obj_altaz[setting].alt < 36.8*u.deg #time steps where object could be below the nasmyth platform, K2
    general_shutter  = obj_altaz[setting].alt < 18*u.deg

    #object never sets
    if np.sum(general_shutter) == 0:
        set_time = None
    else: #object sets
        if obj_altaz[setting][nasmyth_platform][0].az >185.3*u.deg and obj_altaz[setting][nasmyth_platform][0].az < 332.8*u.deg:
            set_time = obj_altaz[setting][nasmyth_platform][0].obstime
        else:
            set_time = obj_altaz[setting][general_shutter][0].obstime
    return set_time


def keck2_wrap_time(sky_coord,nightstart,nightend):
    # Keck II Azimuth range  185.3° to 332.8° , Elevations accessible 36.8° to 89.5°°
    # Stable guiding cannot be guaranteed for targets transiting at elevations higher than 85°.

    t_altaz_start = keck.altaz(nightstart, sky_coord)
    t_altaz_end = keck.altaz(nightend, sky_coord)

    print('azimuth range %.0f-%.0f, dec %s'%(t_altaz_start.az.deg,t_altaz_end.az.deg,dec))
    
#keck2_wrap_time(source_coords,times['startime'],times['endtime'])

def keck1_wrap_time(sky_coord,nightstart,nightend):
    # Keck I   Azimuth range  5.3° to 146.2°   , Elevations accessible 33.3° to 88.9°
    # Stable guiding cannot be guaranteed for targets transiting at elevations higher than 85°.


    t_altaz_start = keck.altaz(nightstart, sky_coord)
    t_altaz_end = keck.altaz(nightend, sky_coord)
    print('azimuth range %.0f-%.0f, dec %s'%(t_altaz_start.az.deg,t_altaz_end.az.deg,dec))
    

def rtime(astroplan_datetime):
    # rounds astroplan datetime to nearest minute, and outputs time
    
    return()


def NIRES_telluric_exp_time(optical_mag):
    optical_mag=float(optical_mag)
    if optical_mag<7:
        exp_time_str='4x2s (bright, keep off slit centre)'
        exp_time_ind=2
    if 7.0 <= optical_mag< 7.5:
        exp_time_str='2 / 4x2'
        exp_time_ind=2
    if 7.5 <= optical_mag< 8:
        exp_time_str='2 / 4x3'
        exp_time_ind=3
    if 8 <= optical_mag< 8.5:
        exp_time_str='2 / 4x5'
        exp_time_ind=5
    if 8.5 <= optical_mag< 9:
        exp_time_str='2 / 4x10'
        exp_time_ind=10
    if 9 <= optical_mag< 9.5:
        exp_time_str='3 / 4x15'
        exp_time_ind=15
    if 9.5 <= optical_mag< 10:
        exp_time_str='3 / 4x20'
        exp_time_ind=20
    if 10 <= optical_mag< 11:
        exp_time_str='3 / 4x30'
        exp_time_ind=30
    if 11 <= optical_mag< 12:
        exp_time_str='3 / 4x40'
        exp_time_ind=40
    if optical_mag>= 12:
        exp_time_str='faint telluric, pick brighter one'
        exp_time_ind=60
    return exp_time_str,exp_time_ind

def NIRES_target_exp_time(optical_mag):
    optical_mag=float(optical_mag)
    if optical_mag<13:
        exp_time_str='20 / 4x60s'
        exp_time_ind=60
        dither="ABBA"
    if 13 <= optical_mag< 14:
        exp_time_str='20 / 4x80'
        exp_time_ind=80
        dither="ABBA"
    if 14 <= optical_mag< 15:
        exp_time_str='30 / 4x100'
        exp_time_ind=100
        dither="ABBA"
    if 15 <= optical_mag< 16:
        exp_time_str='30 / 4x120'
        exp_time_ind=120
        dither="ABBA"
    if 16 <= optical_mag< 16.5:
        exp_time_str='30 / 4x150'
        exp_time_ind=150
        dither="ABBA"
    if 16.5 <= optical_mag< 17:
        exp_time_str='30 / 4x200'
        exp_time_ind=200
        dither="ABBA"
    if 17 <= optical_mag< 17.5:
        exp_time_str='30 / 4x250'
        exp_time_ind=250
        dither="ABBA"
    if 17.5 <= optical_mag< 18:
        exp_time_str='30 / 4x300'
        exp_time_ind=300
        dither="ABBA"
    if 18.0 <= optical_mag< 18.5:
        exp_time_str='40 / 6x300'
        exp_time_ind=300
        dither="ABBAAB"
    if 18.5 <= optical_mag< 19.0:
        exp_time_str='40 / 8x300'
        exp_time_ind=300
        dither="ABBAABBA"
    if optical_mag>= 19:
        exp_time_str='faint, need ~1 hour on it, do you really want to do this?'
        exp_time_ind=300
        dither="ABBAABBA"

    return exp_time_str,exp_time_ind,dither

# Keck I   Azimuth range  5.3° to 146.2°   , Elevations accessible 33.3° to 88.9°
# Keck II Azimuth range  185.3° to 332.8° , Elevations accessible 36.8° to 89.5°°
# Stable guiding cannot be guaranteed for targets transiting at elevations higher than 85°.

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=\
        '''
        Calculates observing contraints for a given night
        
        Usage: observing_contraints.py target_list 
            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('filename',help = 'Final LRIS or NIRES target list.')
    parser.add_argument('instrument',help = 'Instrument. Options: LRIS, NIRES')
    parser.add_argument('--utdate', dest='utdate',type = Time,help = 'fmt: yyyy-mm-dd. Date of observations. Superseeds taking the time from the filename.')
    parser.add_argument("--fhalf", action="store_true",help="first half night")
    parser.add_argument("--shalf", action="store_true",help="second half night")

    args = parser.parse_args()

    lat = 19.8263
    lon = -155.47441
    height= 4145 #4159.581216 # metres
    
    cmt='#'

    #keck = Observer.at_site("Keck")
    keck_location = EarthLocation(lat=lat*u.deg,lon=lon*u.deg,height=height*u.m)
    keck = Observer(location=keck_location, name="keck", timezone="US/Hawaii")
    utc_offset = -10*u.hour #Hawaii standard time
    filename = args.filename

    if args.utdate:
        utdate=args.utdate
        date = Time(utdate)-TimeDelta(1)
        astropy_utdate = Time(utdate)
    else:
        if args.instrument=='NIRES':
         date = filename.split(".")[0].split("Keck_II_")[1].split("_")[0]
         utdate = Time(date)+TimeDelta(1)
        if args.instrument=='LRIS':
         date = filename.split(".")[0].split("Keck_I_")[1].split("_")[0]
         utdate = Time(date)+TimeDelta(1)
        print("Using UT date from the filename: ",utdate)


    # USNO, defintion of sunset time is when the solar disk center is at -0.8333 degrees altitude
    # to account for the solar radius and atmospheric refraction, I found -0.7 matches well with keck
    times = {'utdate':utdate,
            'caldate':date,
            'sunset': keck.sun_set_time(utdate, which='next',horizon=0*u.deg,n_grid_points=200).datetime,
            'etwi05': keck.sun_set_time(utdate, which='next',horizon=-0.5*u.deg,n_grid_points=200).datetime,
            'etwi06': keck.sun_set_time(utdate, which='next',horizon=-0.6*u.deg,n_grid_points=200).datetime,
            'etwi07': keck.sun_set_time(utdate, which='next',horizon=-0.7*u.deg,n_grid_points=200).datetime,
            'etwi4': keck.sun_set_time(utdate, which='next',horizon=-4.0*u.deg,n_grid_points=200).datetime,
            'etwi5': keck.sun_set_time(utdate, which='next',horizon=-5.0*u.deg,n_grid_points=200).datetime,
            'etwi6': keck.sun_set_time(utdate, which='next',horizon=-6.0*u.deg,n_grid_points=200).datetime,
            'etwi8': keck.sun_set_time(utdate, which='next',horizon=-8.0*u.deg,n_grid_points=200).datetime,
            'etwi10': keck.sun_set_time(utdate, which='next',horizon=-10.0*u.deg,n_grid_points=200).datetime,
            'etwi12': keck.sun_set_time(utdate, which='next',horizon=-12.0*u.deg,n_grid_points=200).datetime,
            'etwi18': keck.sun_set_time(utdate, which='next',horizon=-18.0*u.deg,n_grid_points=200).datetime,
            'midnight': keck.midnight(utdate, which='next').datetime,
            'mtwi18': keck.sun_rise_time(utdate, which='next',horizon=-18*u.deg,n_grid_points=200).datetime,
            'mtwi12': keck.sun_rise_time(utdate, which='next',horizon=-12*u.deg,n_grid_points=200).datetime,
            'mtwi10': keck.sun_rise_time(utdate, which='next',horizon=-10*u.deg,n_grid_points=200).datetime,
            'mtwi8': keck.sun_rise_time(utdate, which='next',horizon=-8*u.deg,n_grid_points=200).datetime,
            'mtwi6': keck.sun_rise_time(utdate, which='next',horizon=-6*u.deg,n_grid_points=200).datetime,
            'mtwi5': keck.sun_rise_time(utdate, which='next',horizon=-5*u.deg,n_grid_points=200).datetime,
            'mtwi4': keck.sun_rise_time(utdate, which='next',horizon=-4*u.deg,n_grid_points=200).datetime,
            'mtwi07': keck.sun_rise_time(utdate, which='next',horizon=-0.7*u.deg,n_grid_points=200).datetime,
            'mtwi06': keck.sun_rise_time(utdate, which='next',horizon=-0.6*u.deg,n_grid_points=200).datetime,
            'mtwi05': keck.sun_rise_time(utdate, which='next',horizon=-0.5*u.deg,n_grid_points=200).datetime,
            'sunrise': keck.sun_rise_time(utdate, which='next',horizon=0*u.deg,n_grid_points=200).datetime,

        }
    print('caldate  %s'%times['caldate'])
    print('utdate   %s'%times['utdate'])
    print('sunset   %s'%times['sunset'])
    print('etwi05   %s'%times['etwi05'])
    print('etwi06   %s'%times['etwi06'])
    print('etwi07   %s'%times['etwi07'])
    print('etwi4    %s'%times['etwi4'])
    print('etwi5    %s'%times['etwi5'])
    print('etwi6    %s'%times['etwi6'])
    print('etwi8    %s'%times['etwi8'])
    print('etwi10   %s'%times['etwi10'])
    print('etwi12   %s'%times['etwi12'])
    print('etwi18   %s'%times['etwi18'])
    print('midnight %s'%times['midnight'])
    print('mtwi18   %s'%times['mtwi18'])
    print('mtwi12   %s'%times['mtwi12'])
    print('mtwi10   %s'%times['mtwi10'])
    print('mtwi8    %s'%times['mtwi8'])
    print('mtwi6    %s'%times['mtwi6'])
    print('mtwi5    %s'%times['mtwi5'])
    print('mtwi4    %s'%times['mtwi4'])
    print('mtwi07   %s'%times['mtwi07'])
    print('mtwi06   %s'%times['mtwi06'])
    print('mtwi05   %s'%times['mtwi05'])
    print('sunrise  %s'%times['sunrise'])
    

    # NIRES should start/end ~20 min after/before sunset/sunrise, or 5 degree
    if args.instrument=='NIRES':
        if args.fhalf:
            times.update({'startime':times['etwi5']})
            times.update({'endtime':times['midnight']})
        if args.shalf:
            times.update({'startime':times['midnight']})
            times.update({'endtime':times['mtwi5']})
        if not (args.fhalf or args.shalf):
            print('NIRES FULL NIGHT')
            times.update({'startime':times['etwi5']})
            times.update({'endtime':times['mtwi5']})

    # LRIS can start at like 8 degree ish
    if args.instrument=='LRIS':
        if args.fhalf:
            times.update({'startime':times['etwi8']})
            times.update({'endtime':times['midnight']})
        if args.shalf:
            times.update({'startime':times['midnight']})
            times.update({'endtime':times['mtwi8']})
        if not (args.fhalf or args.shalf):
            print('NIRES FULL NIGHT')
            times.update({'startime':times['etwi8']})
            times.update({'endtime':times['mtwi8']})


    file = open(filename, 'r')
    lines = file.read().split('\n')
    out_file = open('spreadsheet_'+filename, 'w')


    if args.instrument=='NIRES':
     main_nires_header='\tName\tRa\tDec\tMag (B,V)\t Exposure seq (svc/spec)\t Dither \tTime(m)\tIndividual Exp\tTotal Exp(s)\tTotal Exp (min)\tSetting time(UT)\tAzimuth Range\t Wrap\n'
     shalf_nires_header='Start(UT)\tEnd(18,12,5,0 deg) (UT)'
     fhalf_nires_header='Start(0,5,12,18 deg)(UT)\tEnd (UT)'
     full_nires_header='Start(0,5,12,18 deg)(UT)\tEnd(18,12,5,0 deg) (UT)'
     
     if args.shalf:
      header=shalf_nires_header+main_nires_header
      out_file.write(header)
      out_file.write(str(times['midnight'].strftime("%H:%M:%S"))+'\t'+times['mtwi18'].strftime("%H:%M:%S")+','+times['mtwi12'].strftime("%H:%M:%S")+','+times['mtwi5'].strftime("%H:%M:%S")+','+times['mtwi07'].strftime("%H:%M:%S")+'\n')
      out_file.write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
 
     if args.fhalf:
      header=fhalf_nires_header+main_nires_header
      out_file.write(header)
      out_file.write(str(times['startime'])+'\t'+str(times['endtime']))
      out_file.write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
 
     if not (args.fhalf or args.shalf):
      header=full_nires_header+main_nires_header
      out_file.write(header)
      out_file.write(str(times['startime'])+'\t'+str(times['endtime']))
      out_file.write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')


    if args.instrument=='LRIS':
        # need to copy lris log...
        main_lris_header='\tName\tRa\tDec\tMag (B,V)\t Exposure seq (svc/spec)\t Dither \tTime(m)\tIndividual Exp\tTotal Exp(s)\tTotal Exp (min)\tSetting time(UT)\tAzimuth Range\t Wrap\n'
        #shalf_nires_header='Start(UT)\tEnd(18,12,5,0 deg) (UT)'
        #fhalf_nires_header='Start(0,5,12,18 deg)(UT)\tEnd (UT)'
        #full_nires_header='Start(0,5,12,18 deg)(UT)\tEnd(18,12,5,0 deg) (UT)'

       
        # To do, update this with LRIS specific



    for ind, i in enumerate(lines):
        if len(i.strip()) > 1:
            if i.strip()[0] != cmt : #The entire line is not commented
                if "_S" not in i.split()[0]: #this is a target, not an offset star

                    split=i.split()
                    #print(len(split))                    
                    name=split[0]
                    ra = split[1]+' '+split[2]+' '+split[3]
                    dec= split[4]+' '+split[5]+' '+split[6]
                    #print(len(split),name,ra,dec,split[7:])
                    ra_sheet = "'"+split[1]+':'+split[2]+':'+split[3]
                    dec_sheet= "'"+split[4]+':'+split[5]+':'+split[6]

                    source_coords = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg))
                    
                    #t_altaz_midnight = keck.altaz(times['midnight'], source_coords)
                    t_altaz_start = keck.altaz(times['startime'], source_coords)
                    t_altaz_end = keck.altaz(times['endtime'], source_coords)

                    print('azimuth range %.0f-%.0f, dec %s'%(t_altaz_start.az.deg,t_altaz_end.az.deg,dec))

                    if args.instrument=='LRIS':
                     k1_rise=keck1_rising_time(source_coords, times['utdate'])
                     #print(name+" K1 rising time is ", k1_rise)

                    if args.instrument=='NIRES':
                     k2_set=keck2_setting_time(source_coords, times['utdate'])
                     #print(name+" K2 setting time is ",k2_set)
    

                     #print(type(k2_set.datetime))
                     #print(k2_set.datetime.strftime("%H:%M:%S"))

                     if args.instrument=='NIRES':
                      if "HIP" in name:
                       vmag=split[17]
                       bmag=split[14]
                       dither='ABBA'
                       exp_time_str=NIRES_telluric_exp_time(vmag)[0]
                       exp_time=NIRES_telluric_exp_time(vmag)[1]
                       telluric_tot='00:10:00'
                       #total_exp_time_s=exp_time*len(dither)
                       #total_exp_time_m=total_exp_time_s/60.
                       out_file.write('\t\t %s \t %s \t %s \t %s,%s\t%s\t%s\t%s \t %s \t \t \t %s\t %.0f-%.0f \n'%\
                         (name,ra_sheet,dec_sheet,bmag,vmag,exp_time_str,dither,telluric_tot,exp_time,k2_set,t_altaz_start.az.deg,t_altaz_end.az.deg))
                      
                      if not "HIP" in name:
                       sn_mag=split[11].split('r=')[1]
                       exp_time_str=NIRES_target_exp_time(sn_mag)[0]
                       exp_time=NIRES_target_exp_time(sn_mag)[1]
                       dither=NIRES_target_exp_time(sn_mag)[2]

                       total_exp_time_s=exp_time*len(dither)
                       total_exp_time_m=math.ceil(total_exp_time_s/60.)
                       time_w_overheads=total_exp_time_m+5 
                       out_file.write('\t\t%s\t%s\t%s\t%s \t %s \t %s \t 00:%s:00 \t %s \t %s \t %.1f \t %s \t %.0f-%.0f \n'%\
                         (name,ra_sheet,dec_sheet,sn_mag,exp_time_str,dither,time_w_overheads,exp_time,total_exp_time_s,total_exp_time_m,k2_set,t_altaz_start.az.deg,t_altaz_end.az.deg))



