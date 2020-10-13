"""
This script takes the downloaded target list from YSE PZ and loops through the targets,
producing finder chart and finding offset stars for them. If host is indicated in the input
target list in a line below the target, the host and the PA are included in the finder chart. 

In the end, you get a finder chart for each target, and an updated target list with the PA 
(if host is provided) and the offset stars. 
"""
import astropy.io.ascii as asci
import astropy.units as u
from astropy import coordinates
import sys, argparse
import numpy as np

from finder_chart import get_finder




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=\
        '''
        Creates the finder charts for the whole night, provided a target list file.
        
        Usage: prepare_obs_run.py target_list_filename
            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    #Check if correct number of arguments are given
    if len(sys.argv) != 2:
    	print ("Usage: prepare_obs_run.py target_list_filename")
    	sys.exit()

    filename = sys.argv[1]

    targets = asci.read(filename, comment = '!', format = 'no_header')

    names = np.array(targets['col1'])

    RA = np.array([str(x['col2']) + ':' +str(x['col3']) +':' + str(x['col4']) for x in targets])
    Dec = np.array([str(x['col5']) + ':' +str(x['col6']) +':' + str(x['col7']) for x in targets])

    coords = coordinates.SkyCoord(RA, Dec, unit = (u.hourangle, u.deg))

    mags = np.array(targets['col11'])

    all_starlist = ''

    skip_host = False
    ######Now, call finder_chart.py
    for i in range(len(targets)):
        if skip_host == False: #this entry is not host, treat as target
            name = names[i]
            ra_deg = coords[i].ra.deg
            dec_deg = coords[i].dec.deg
            mag = mags[i]
            print("Preparing a finder chart for %s."%name)
            #Check if host is provided
            if i < len(targets)-1: #not the last item
                next_entry = names[i+1]
                if next_entry.split('_')[0] == name and next_entry.split('_')[1] == 'host':
                    host_ra = coords[i+1].ra.deg
                    host_dec = coords[i+1].dec.deg
                    print("Host galaxy coordinates provided. RA = %.5f Dec = %.5f"%(host_ra, host_dec))

                    skip_host = True #Because of this flag, the next item is skipped. 
                else:
                    host_ra = None
                    host_dec = None
            else: #For the last item, no host is provided. 
                host_ra = None
                host_dec = None
            
            # print(host_ra, host_dec)
            finder_size = 5/60 #5 arcmin, LRIS
            max_separation = 2*60 #2 arcmin, LRIS

            starlist_entry = get_finder( ra_deg, dec_deg, name,  finder_size, mag = mag, \
                            minmag=15, maxmag=18.5, num_offset_stars = 3, min_separation = 2, max_separation = max_separation,\
                            host_ra = host_ra, host_dec = host_dec, \
                            starlist=None, print_starlist = False,  return_starlist = True, debug = False)
            print(starlist_entry)
            all_starlist += starlist_entry
        elif skip_host == True: #In this case, the next item should be a target
            skip_host = False 
    print(all_starlist)
    ###Write to file
    out_file = open(filename.split('.')[0]+'_final.txt', "w")
    out_file.write(all_starlist)
    out_file.close()

