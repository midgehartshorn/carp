'''
This is the Combined Astro Reduction Pipeline (CARP).
Contributors:
    Midge Hartshorn
    Sophie Booth

BASE FUNCTIONS:
    collect_files(datadir, cal_status)
    median_combine(filelist)
    bias_subtract(filename, path_to_master_bias, output_path)
    scale_master_dark(path_to_master_dark, desired_exptime)
    dark_subtract(filename, path_to_master_dark, outpath)
    norm_combine_flats(filelist)
    flatfield(filename, path_to_master_flat, outpath)
AUTOMATED FUNCTIONS:
    create_master_bias(biaspath)
    create_master_dark(darkpath, biaspath, exposuretime)
    create_master_flat(flatpath, darkpath, biaspath, targ, filter, exposuretime)
    reduce_images(datapath, biaspath, darkpath, flatpath, objectname, filtername, exposuretime)
    align_and_stack()
    
'''

import glob
import os
import numpy as np
from astropy.io import fits
import scipy.signal

def collect_files(datadir, cal_status):
    filelistB = glob.glob(datadir + cal_status + '*B.fit')
    filelistV = glob.glob(datadir + cal_status + '*V.fit')
    filelistR = glob.glob(datadir + cal_status + '*R.fit')
    return filelistB, filelistV, filelistR

def median_combine(filelist):
    '''
    a function that combines a list of files to create a median master 
    INPUT
    filelist: a list of files

    OUTPUT
    med_frame: a median image

    DEPENDENCIES
    astropy.io.fits.getdata
    numpy
    '''
    # number of input files
    n = len(filelist)
    
    # examining the first file data so we can find the dimensions
    first_frame_data = fits.getdata(filelist[0])
   
    # assigning the dimensions
    imsize_y, imsize_x = first_frame_data.shape
  
    # creating the empty 3d stack of the image dimensions x number of frames
    fits_stack = np.zeros((imsize_y, imsize_x , n)) 
                                                                
    # filling the empty stack with the pixel values from each data image
    for ii in range(0, n):
        im = fits.getdata(filelist[ii])
        fits_stack[:,:,ii] = im
    # taking the median of each pizel along the z axis (creating a median image of dimensions x,y)        
    med_frame = np.median(fits_stack, axis = 2)
    return med_frame

def create_master_bias(biaspath):
    '''
    A function to create a master bias from a list of bias files.
    
    INPUT
    biaspath: the LOCATION of the bias calibration files. These should be named 'cal*bias.fit'.
    
    OUTPUT
    a file named 'Master_Bias.fit'
    
    DEPENDENCIES
    astropy.io.fits
    glob
    '''
    biasfiles = glob.glob(biaspath+'cal*bias.fit')
    biasheader = fits.getheader(biasfiles[0])
    median_bias = median_combine(biasfiles)
    output_name = 'Master_Bias.fit'
    fits.writeto(biaspath + output_name, median_bias, biasheader, overwrite=True)
    print(f'Master Bias values are: {median_bias}')
    return
 
def bias_subtract(filename, path_to_masterbias, output_path):
    '''
    A function to subtract the master bias from another frame and write the subtracted file.

    INPUT
    filename: a single filename
    path_to_master_bias: path to the master bias calibration file
    output_path: path to write-to for subtracted file

    OUTPUT
    no return statement; writes calibrated file to output_path

    DEPENDENCIES
    astropy.io.fits
    '''
    raw_file = fits.getdata(filename)
    print('Processing file:', filename)
    masterbias = fits.getdata(path_to_masterbias)
    print('Subtracting bias...')
    cal_file = raw_file - masterbias
    # We'll just take the very last part of the filename, as described above
    output_filename = filename.split('/')[-1]
    
    # make sure to include header
    raw_header = fits.getheader(filename)

    # include a 'b' to indicate that this file has been bias subtracted
    fits.writeto(output_path + 'b_' + output_filename, cal_file, raw_header, overwrite=True)
    print('Wrote subtracted file to:', output_path,'b_',output_filename)
    return 
def create_master_dark(darkpath, biaspath, exposure_time):
    '''
    A function to create a master dark file, with a given exposure time
    
    INPUT
    darkpath: the LOCATION of the dark calibration files, filename format 'cal*{exposure_time}.fit'
    biaspath: the LOCATION of the master bias
    exposure_time: the unaltered exposure time of the darks (do NOT include '.' when defining exposure time)
    
    OUTPUT
    a master dark named 'Master_Dark_{exposure_time}s.fit'
    
    DEPENDENCIES
    astropy.io.fits
    glob
    '''
    darkfiles = glob.glob(darkpath + 'cal*' + exposure_time + '.fit')
    master_bias_path = biaspath + 'Master_Bias.fit'
    for i in range(len(darkfiles)):
        bias_subtract(darkfiles[i],master_bias_path, darkpath)
    darkheader = fits.getheader(darkfiles[0])
    b_darkframes = glob.glob(darkpath+'b*'+exposure_time+'.fit')
    median_darks = median_combine(b_darkframes)
    fits.writeto(darkpath + 'Master_Dark_'+exposure_time+'s.fit',median_darks,darkheader,overwrite=True)
    print(f'Master Dark of {exposure_time}s has values: {median_darks}')
    return

def scale_master_dark(path_to_master_dark, desired_exptime):
    '''
    A function to scale the master dark to the appropriate exposure time, and write to current directory.

    INPUT
    path_to_master_dark: absolute path to master dark frame
    desired_exptime: final exposure time (integer; as '.' will mess up output filename)

    OUTPUT
    no return statement; writes scaled dark to current directory
    
    DEPENDENCIES
    astropy.io.fits
    '''
    master_dark = fits.getdata(path_to_master_dark)
    d_head = fits.getheader(path_to_master_dark)
    # find the initial exposure time
    dark_exp = d_head['EXPTIME']
                            
    if dark_exp == desired_exptime:
        print('NO ADJUSTMENT NECESSARY')
    else:
        counts = master_dark / dark_exp
        scaled_dark = desired_exptime * counts
        # update the header
        d_head['EXPTIME'] = desired_exptime
        output_filename = 'Master_Dark_' 
        fits.writeto(output_filename + str(desired_exptime) +'s.fit', scaled_dark, d_head, overwrite=True)
    return


def dark_subtract(filename, path_to_master_dark, outpath):
    '''
    A function to subtract the master dark from a given frame
    
    INPUT
    filename: a single filename
    path_to_master_dark: absolute path to master dark frame
    outpath: path for write-to of calibrated file

    OUTPUT
    no return; writes calibrated file to specified outpath 

    DEPENDENCIES
    astropy.io.fits
    '''
    int_file = fits.getdata(filename)
    int_header = fits.getheader(filename)
    int_exptime = int_header['EXPTIME']
    masterdark = fits.getdata(path_to_master_dark)
    master_head = fits.getheader(path_to_master_dark)
    master_exptime = master_head['EXPTIME']
    print('Processing file: ', filename)
                                        
    # error handling:
    # use an if statement below to check to see whether the exposure time  
    # of the image matches the exposure time of the master dark.
    # if the exposure times don't match, throw an error!
    if master_exptime != int_exptime:
        masterdark = scale_master_dark(path_to_master_dark, master_exptime)
    print('Subtracting dark current...')
    cal_file = int_file - masterdark
    # We'll just take the very last part of the filename, as described above
    output_filename = filename.split('/')[-1]
   
    # Make sure it has the correct header! 'd' indicates file has been dark-subtracted
    fits.writeto(outpath + 'd' + output_filename, cal_file, int_header, overwrite=True) 
    print('Wrote subtracted file to: ', outpath, 'd', output_filename)
    return 

def norm_combine_flats(filelist):
    '''
    A function to normalize and combine a flat

    INPUT
    filelist: a list of files

    OUTPUT
    med_frame: a single, median combined, normalized frame; not written.

    DEPENDENCIES
    astropy.io.fits
    numpy
    '''
    # we need to know how many files, or how many times to iterate
    n = len(filelist)
    
    # get the data for the first frame
    first_frame_data = fits.getdata(filelist[0])
    
    # what are the dimensions of our frames?
    imsize_y, imsize_x = first_frame_data.shape
    
    # create an array of zeros with the dimensions of our images x the number of images
    fits_stack = np.zeros((imsize_y, imsize_x , n))
    
    # this step normalizes the flats and populates the earlier array
    for ii in range(0, n):
        im = fits.getdata(filelist[ii])
        norm_im =  im / np.median(im) # normalize flats
        fits_stack[:,:,ii] = norm_im
        
    # create our median!   
    med_frame = np.median(fits_stack, axis=2)
    return med_frame


def create_master_flat(flatpath,biaspath,darkpath,targname,filtername,exposuretime):
    '''
    A function to create a master flat in a given band
    
    INPUT
    flatpath: the LOCATION of the flats
    biaspath: the LOCATION of the master bias
    darkpath: the LOCATION of the master dark
    targname: the root/target of the filename, e.g. 'twilight', 'dome', etc.
    filtername: the band name, as a string
    exposuretime: the exposure time of the flats, as a string
    
    OUTPUT
    a master flat in the same location as the original flats, filename: 'Master_Flat_{filtername}band.fit'
    
    DEPENDENCIES
    glob
    astropy.io.fits
    '''
    flatfiles = glob.glob(flatpath+targname+'*_'+filtername+'.fit')
    master_bias_path = biaspath + 'Master_Bias.fit'
    master_dark_path = darkpath + 'Master_Dark_'+exposuretime+'s.fit'
    for i in range(len(flatfiles)):
        bias_subtract(flatfiles[i],master_bias_path, flatpath)
    b_flatfiles = glob.glob(flatpath+'b_'+targname+'*'+filtername+'*.fit')
    for i in range(len(b_flatfiles)):
        dark_subtract(b_flatfiles[i],master_dark_path, flatpath)
    db_flatfiles = glob.glob(flatpath+'db_'+targname+'*'+filtername+'*.fit')
    flats_combined = norm_combine_flats(db_flatfiles)
    flat_header = fits.getheader(db_flatfiles[0])
    fits.writeto(flatpath+'Master_Flat_'+filtername+'band.fit', flats_combined, flat_header, overwrite=True)
    print(f'The Master Flat of {exposuretime}s is: {flats_combined}')
    return

def flatfield(filename, path_to_master_flat, outpath):
    '''
    A function to divide a file by a master flatfield.

    INPUT
    filename: an individual file
    path_to_master_flat: absolute path to master flat
    outpath: where to write the calibrated file

    OUTPUT
    no return statement, but writes calibrated file to location specified by outpath

    DEPENDENCIES
    astropy.io.fits
    '''
    cal_file = fits.getdata(filename)
    print('Processing file: ', filename)
            
    masterflat = fits.getdata(path_to_master_flat)
    print('Flatfielding...')
    final_file = cal_file / masterflat
    raw_header = fits.getheader(filename)
    output_filename = filename.split('/')[-1]
    fits.writeto(outpath + 'f' + output_filename, final_file, raw_header, overwrite=True)
    print('Wrote file to: ', outpath,'f',output_filename)
    return

def reduce_images(datapath, biaspath, darkpath, flatpath, objectname, filtername, exposuretime):
    '''
    A function to perform arithmetic reduction on all raw files in a given directory
    
    INPUT 
    datapath: LOCATION of data files
    biaspath: LOCATION of master bias
    darkpath: LOCATION of master dark (match exp times)
    flatpath: LOCATION of master flat (match exp times)
    objectname: root name of object (only as much as needed to differentiate)
    filtername: name of filter, as string
    exposuretime: exposure time of raw images, avoid using '.'
    
    OUTPUT
    flatfielded, bias- and dark-subtracted files, written to same directory as raw files.
    DEPENDENCIES
    glob
    
    '''
    datafiles = glob.glob(datapath+objectname+'*_'+filtername+'.fit')
    master_bias_path = biaspath + 'Master_Bias.fit'
    master_dark_path = darkpath + 'Master_Dark_'+exposuretime+'s.fit'
    master_flat_path = flatpath + 'Master_Flat_'+filtername+'band.fit'
    
    for i in range(len(datafiles)):
        bias_subtract(datafiles[i],master_bias_path,datapath)
    b_data_files = glob.glob(datapath+'b_'+objectname+'*'+filtername+'.fit')
    for i in range(len(b_data_files)):
        dark_subtract(b_data_files[i],master_dark_path,datapath)
    db_data_files = glob.glob(datapath+'db_'+objectname+'*'+filtername+'.fit')
    for i in range(len(db_data_files)):
        flatfield(db_data_files[i],master_flat_path,datapath)
    fdb_data_files = glob.glob(datapath+'fdb_'+objectname+'*'+filtername+'.fit')
    print(f'List of fully reduced files of {objectname} with filter {filtername}: {fdb_data_files}')
    return 
 
def cross_image(im1, im2, **kwargs):
    """
    Takes two images, subtracts their averages. Creates a correlation image using FFT and locates the peak value
    and peak signal position. Returns the pixel shifts.

    DEPENDENCIES
    scipy.signal.fftconvolve
    """
    
    # The type cast into 'float' is to avoid overflows:
    im1_gray = im1.astype('float')
    im2_gray = im2.astype('float')

       # Enable a trimming capability using keyword argument option.
    if 'boxsize' in kwargs:
        im1_gray = im1_gray[0:kwargs['boxsize'],0:kwargs['boxsize']]
        im2_gray = im2_gray[0:kwargs['boxsize'],0:kwargs['boxsize']]

    # Subtract the averages of im1_gray and im2_gray from their respective arrays -- cross-correlation
    # works better that way.
    im1_gray -= np.mean(im1_gray)
    im2_gray -= np.mean(im2_gray)

    # Calculate the correlation image using fast Fourier transform (FFT)
    # Note the flipping of one of the images (the [::-1]) - this is how the convolution is done.
    corr_image = scipy.signal.fftconvolve(im1_gray, im2_gray[::-1,::-1], mode='same')
    # To determine the location of the peak value in the cross-correlated image, 
    # using np.argmax on the correlation image:
    peak_corr_index = np.argmax(corr_image)  

    # Find the peak signal position in the cross-correlation -- this gives the shift between the images.
    corr_tuple = np.unravel_index(peak_corr_index, corr_image.shape)
    
    # Calculate shifts (not cast to integer, but could be).
    xshift = corr_tuple[0] - corr_image.shape[0]/2.
    yshift = corr_tuple[1] - corr_image.shape[1]/2.

    return xshift,yshift

def shift_image(image,xshift,yshift):
    '''
    A function to shift the image by the previously determined x and y shifts
    
    INPUT
    image: the image to shift
    xshift: how many pixels to shift in the x-direction
    yshift: how many pixels to shift in the y-direction

    OUTPUT
    image, but rolled along both axes
    '''
    # Note that this will not do any trimming, 
    # so we'll want to  trim later the edges of the image using the maximum shift.
    return np.roll(np.roll(image,int(yshift),axis=1), int(xshift), axis=0)


def align_and_stack(targs, filters, datadir):
    '''
    A function to align and shift a series of images to a common reference frame.

    INPUT
    targs: list of target names
    filters: list of filters
    datadir: location of files to be aligned

    OUPUT

    DEPENDENCIES
    '''
    for targname in targs:
        print(' ')
        print('-----------------------------')      
        print('target: ', targname)
        print('-----------------------------')      

        # Using glob, make list of all reduced images of current target in all filters.
        # Complete the following line to create a list of the correct images to be shifted (use wildcards!):
        imlist = glob.glob(datadir + 'fdb*.fit')
        print(imlist)
        print(type(imlist))
        # sort the list to make sure it's in a reasonable file order:
        imlist.sort()    
        
        # Check to make sure that your new list has the right files:
        print("All files to be aligned: \n", imlist)
        print('\n') 
        
        # Open first image = master image; all other images of same target will be aligned to this one. 
        im1,hdr1 = fits.getdata(imlist[0],header=True)
        print("Aligning all images to:", imlist[0])
        
        print('\n') # adding some space to the print statements

        # working through the list of images, we get the data and header for each image, and calculate the x and y shifts
        # by running our cross-correlation function; we add these shifts to the respective shift-arrays, and print them
        # in a nice, neat manner.
          
        xshifts = {}
        yshifts = {}
        for index,filename in enumerate(imlist):
            im,hdr = fits.getdata(filename,header=True)
            xshifts[index], yshifts[index] = cross_image(im1, im)
            print("Shift for image", index, "is", xshifts[index], yshifts[index])

            # Calculate trim edges of new median stacked images so all stacked images of each target have same size 
        max_x_shift = int(np.max([xshifts[x] for x in xshifts.keys()]))
        max_y_shift = int(np.max([yshifts[x] for x in yshifts.keys()]))
        print('   Max x-shift={0}, max y-shift={1} (pixels)'.format(max_x_shift,max_y_shift))

            # Cycle through list of filters
        for filtername in filters:
        # for-loop + if-statement to create a list of FITS files matching *only* the selected filter:
            scilist = []
            for fitsfile in imlist:
                header = fits.getheader(fitsfile)
                if header['FILTER'] == filtername:
                    scilist.append(fitsfile)
            if len(scilist) < 1:
                print("Warning! No files in scilist. Your path is likely incorrect.")
                break
                        
                    # Complete the print statement below including the filename & ensuring each scilist entry has the right filter:
            for fitsfile in scilist:
                header = fits.getheader(fitsfile)
                print("filename:",fitsfile,'\t',"goal filter:",filtername,'\t',"image filter:",header['FILTER'])
                        
            nfiles = len(scilist)
            print('Stacking ', nfiles, filtername, ' science frames')

            # Define new array with same size as master image
            image_stack = np.zeros([im1.shape[0],im1.shape[1],len(scilist)])
                        
            # Populating the empty array with the shifted images.
            xshifts_filt = {}
            yshifts_filt = {}
            # Make a new directory in datadir for the new stacked fits files
            if os.path.isdir(datadir + 'Stacked') == False:
                os.mkdir(datadir + 'Stacked')
                print('\n Making new subdirectory for stacked images:', datadir + 'Stacked \n')
                           
            for index,filename in enumerate(scilist):
                im,hdr = fits.getdata(filename,header=True)
                xshifts_filt[index], yshifts_filt[index] = cross_image(im1, im)
                image_stack[:,:,index] = shift_image(im,xshifts_filt[index], yshifts_filt[index])
                # Save the final stacked images into the new folder:
                fits.writeto(datadir + 'Stacked/' + targname+ '_' + str(index) + filtername + '_stacked.fits', image_stack[:,:,index],hdr, overwrite=True)
                print('   Wrote FITS file ',targname+ '_' + str(index)+filtername+'_stacked.fits', 'in ',datadir + 'Stacked/','\n')
                                
                # take the median of the image stack (median combine the stacked images);
            median_image = np.nanmedian(image_stack,axis=2)

            # Sets the new image boundaries
            if (max_x_shift > 0) & (max_y_shift > 0): # don't apply cut if no shift!
                median_image = median_image[max_x_shift:-max_x_shift,max_y_shift:-max_y_shift]

               # Save the final combined image into your new folder:
            fits.writeto(datadir + 'Stacked/' + targname + '_' + filtername + 'stack.fits', median_image, hdr, overwrite=True)
            print('   Wrote FITS file ',targname+'_'+filtername+'stack.fits', 'in ',datadir + 'Stacked/','\n')
    print('\n Done stacking!')
    return