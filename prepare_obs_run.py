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
import os

from finder_chart import get_finder, get_host_PA_and_sep




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=\
        '''
        Creates the finder charts for the whole night, provided a target list file.
        
        Usage: prepare_obs_run.py target_list_filename telescope
            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    # parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str,
                        help="filename of the target list")

    parser.add_argument("telescope", type = str,
                        help="Keck or Lick")

    parser.add_argument("-r", "--rotate", action="store_true",
                        help="produce rotated finder chart")

    parser.add_argument("-d", "--debug", action="store_true",
                    help="debug mode")

    args = parser.parse_args()

    #Check if correct number of arguments are given
    # if len(sys.argv) != 2:
    # 	print ("Usage: prepare_obs_run.py target_list_filename")
    # 	sys.exit()

    # filename = sys.argv[1]
    filename = args.filename
    converters = {'col5': [asci.convert_numpy(np.str)]} #This is so that -00 gets preserved as -00

    targets = asci.read(filename, comment = '!', format = 'no_header', converters = converters)

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
            #Check if host is provided. This should be in the next line, with same name
            # Example
            # 2020pni        hh mm ss dd mm ss etc
            # 2020pni_host   hh mm ss dd mm ss etc
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
            if args.telescope == 'Keck':
                finder_size = 3/60 #3 arcmin, LRIS
                max_separation = 4*60 #2 arcmin, LRIS
                min_mag = 11
                max_mag = 21
            elif args.telescope == 'Lick':
                finder_size = 6/60 #4 arcmin, Kast
                max_separation = 6*60 #3 arcmin, Kast 
                min_mag = 5
                max_mag = 18   
            else:
                print("Telescope should be eihter Keck or Lick; default to Lick.")  
                finder_size = 4/60 #4 arcmin, Kast
                max_separation = 3*60 #3 arcmin, Kast 
                min_mag = 11
                max_mag = 17  

            #Minimum separation
            minsep = 60    

            #Obtain PA and separation from target ra/dec and host ra/dec
            #To do: make it not duplicate for the "get_finder" function. 
            if (host_ra is not None) and (host_dec is not None):
                host_pa, host_sep = get_host_PA_and_sep(ra_deg, dec_deg, host_ra, host_dec)
                pa_offset = 30 #30 degree offset in slit viewing camera PA

            starlist_entry = get_finder( ra_deg, dec_deg, name,  finder_size, mag = mag, \
                            minmag=min_mag, maxmag=max_mag, num_offset_stars = 3, min_separation = minsep, max_separation = None,\
                            host_ra = host_ra, host_dec = host_dec, \
                            starlist=None, print_starlist = False,  return_starlist = True, debug = args.debug, output_format='png')
            print(starlist_entry)
            all_starlist += starlist_entry
            ###If requesting rotated finder chart
            if args.rotate and ((host_ra is not None) and (host_dec is not None)):
                print("Make rotated finder chart using PA + 30 deg for LRIS.")
                try:
                    finder_name = name+'_finder.png'
                    rotated_finder = name+'_finder_rot.png'
                    os.system('convert %s -virtual-pixel white +distort SRT %.2f %s'%(finder_name, host_pa+pa_offset, rotated_finder))
                except:
                    print("Check if imagemagick is installed.")
        elif skip_host == True: #In this case, the next item should be a target
            skip_host = False 
    print(all_starlist)
    ###Write to file
    out_file = open(filename.split('.')[0]+'_final.txt', "w")
    out_file.write(all_starlist)
    out_file.close()

