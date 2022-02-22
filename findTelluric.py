import numpy as np
from astropy.io import ascii as asci
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys

def getTelluric(coords, max_distance = 20, spec_type = 'A0V', Vmin = 4, Vmax = 11):
    """
    Given a source coordinates in SkyCoord format, query the Hipparcos catalog on Vizier to get the closest A0V star. 
    Return a result table with RA, Dec, Spectral type, Vmag, distance from source, and RA difference    
    """
    #Define what columns to query and what the conditions are
    hip_viz = Vizier(columns=['HIP', 'RAhms', 'DEdms', 'Vmag', 'SpType'], \
               column_filters={"Vmag":"%s..%s"%(str(Vmin), str(Vmax)), "SpType":spec_type})

    #query. Use twice the max distance and add flag if the closest one is too far
    result = hip_viz.query_region(coords, radius=2*max_distance*u.deg, catalog='I/239')
    # try:
    a0v_res = result[0] #[0] is Hipparcos; [1] is Tycho

    res_coords = SkyCoord(ra = a0v_res['RAhms'], dec = a0v_res['DEdms'], unit = (u.hourangle, u.deg))
    a0v_res['distance'] = coords.separation(res_coords)
    a0v_res['dRA'] = (coords.ra.deg - res_coords.ra.deg)
    a0v_res['abs_dRA'] = np.abs(coords.ra.deg - res_coords.ra.deg)

    #print(a0v_res)

    # if sortby == 'dra':
    #     a0v_res.sort('dRA_(source-cal)')
    # else:
    #     a0v_res.sort('distance')
    # except:
    #     print("No nearby A0V telluric")
    #     a0v_res = None

    return a0v_res

#Now the part we execute

if __name__ == "__main__":

    format = 'growth' #or 'iobserve'
    #growth is name hh mm ss dd mm ss !comments
    #iobserve is name hh:mm:ss dd:mm:ss #comments

    if len(sys.argv) == 1:
        file_name = input("Type file_name of coordinates ")
    elif len(sys.argv) == 2:
        file_name = sys.argv[1]
    elif len(sys.argv) == 3:
        file_name = sys.argv[1]
        format = sys.argv[2]

    if format == 'growth':
        cmt = '!'
        spc = ' '
        print("format is %s"%format)
    else:
        cmt = '#'
        spc = ':'
        print("format is %s"%format)
    #Open the file
    file = open(file_name, 'r')
    out_file = open('with_telluric_'+file_name, 'w')

    lines = file.read().split('\n')
    for i in lines:
    #     print(i.strip())
        if len(i.strip()) > 1:
            if i.strip()[0] != cmt :
                if "Telluric" not in (i.split(cmt)[1]): #Not a commented line and not a Telluric source
                    splitted = i.split()
                    if format == 'iobserve':
                        source_coords  = SkyCoord(ra =splitted[1], dec = splitted[2], unit = (u.hourangle, u.deg))
                    elif format == 'growth':
                        ra = splitted[1]+' '+splitted[2]+' '+splitted[3]
                        dec = splitted[4]+' '+splitted[5]+' '+splitted[6]
                        source_coords = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg))
                    tellurics = getTelluric(source_coords)
                    #Write source to output file
                    out_file.write(i+'\n')
                    #Write best telluric to output file
                    if tellurics is not None:
                        #Closest telluric in distance
                        #Closest telluric in RA
                        best_tel_dist = tellurics[tellurics['distance'] == np.min(tellurics['distance'])][0]
                        tellurics.sort('abs_dRA')
                        best_tel_ra = tellurics[0]

                        # print(best_tel_ra)
                        outstring1 = ('HIP'+str(best_tel_dist['HIP'])).ljust(20)+best_tel_dist['RAhms'].replace(' ',spc)+'  '+best_tel_dist['DEdms'].replace(' ',spc)+\
                                        '  2000.0  '+cmt+' Telluric V = %.2f distance = %.2f deg, dRA = %.2f deg \n'%(best_tel_dist['Vmag'], best_tel_dist['distance'], best_tel_dist['dRA'])
                        outstring2 = ('HIP'+str(best_tel_ra['HIP'])).ljust(20)+best_tel_ra['RAhms'].replace(' ',spc)+'  '+best_tel_ra['DEdms'].replace(' ',spc)+\
                                        '  2000.0  '+cmt+' Telluric V = %.2f distance = %.2f deg, dRA = %.2f deg \n'%(best_tel_ra['Vmag'], best_tel_ra['distance'], best_tel_ra['dRA'])   
                        out_file.write(outstring1)   
                        out_file.write(outstring2)                                 
                    else:
                        outstring = cmt+' No telluric found for the source above'
                        out_file.write(outstring+'\n')
            else: #for commented lines, just add it back
                out_file.write(i+'\n')
    out_file.close()
    file.close()