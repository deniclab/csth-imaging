#### batch_merge_channels.py ####
# a Jython script for use with Fiji to merge channels and generate maximum
# intensity projections from confocal microscopy images for viewing in
# Fiji/ImageJ.
#
# Author: Nicholas Weir, Denic Laboratory, Harvard University
# nweir@fas.harvard.edu
# v. 1.0.0 8/23/2016
# v. 1.0.1 8/24/2016 placing all outputs in one folder for easier movement
# 
# Please feel free to use, distribute, and modify this code however you see fit;
#however, if you publish anything generated using any of my code, please
# provide an acknowledgement!
#
####  IMPORTANT USAGE NOTES  ####
#
# This script is intended for use at the command line for headless operation of
# ImageJ/Fiji. To use this script, enter the following command:
#
# ./ImageJ-macosx --headless /path/to/script/batch_merge_channels.py img_directory first_letter_of_colors
#
# Replace ./ImageJ-macosx with your platform's ImageJ launcher.
# To indicate which channels are to be merged/projected, enter a string of the
# first character in each color name in any order you like; for example, to
# merge red, yellow, and cyan, you may enter ryc, cyr, etc. as an argument. use
# w for brightfield.
print "importing dependencies..."
from ij import IJ, ImagePlus, ImageStack, WindowManager
from ij.plugin import RGBStackMerge, ZProjector
import os
import sys

def mk_fname_ref(img_list, wlength_string, delimiter = '_'):
    '''Create a reference dict for matching stage position images from
    different channels.
    
    Keyword arguments:
    img_list: the list of image files being worked on. contains all
    wavelengths.
    wlength_string: string containing identifying information for the channel
    of interest, e.g. '488' for images with filenames containing
    'w3.488.laser.25'.
    delimiter (optional): the delimiter used to break up the filename for
    removing the wavelength information.

    Returns: a dictionary whose key:value pairs are:
        key: the filename with the wavelength identifying fragment removed
        value: the full filename of the image
        for all images in the wavelength defined.
        Example of a key:value pair:
            a_stage_pos_fname.tf:a_stage_pos_wavelengthinfo_fname.tif
    '''

    wlength_imlist = [i for i in img_list if wlength_string.lower() in i.lower()]
    wlength_rm = []
    for fname in wlength_imlist:
        split_fname = fname.split(delimiter)
        fname_no_wlength = '_'.join([x for x in split_fname if wlength_string.lower()
                                     not in x.lower()])
        wlength_rm.append(fname_no_wlength)
    return dict(zip(wlength_rm, wlength_imlist))

def maxZprojection(stackimp):
    '''copied from EMBL CMCI python-ImageJ cookbook.'''
    zp = ZProjector(stackimp)
    zp.setMethod(ZProjector.MAX_METHOD)
    zp.doProjection()
    zpimp = zp.getProjection()
    return zpimp

print "creating output directory..."
img_dir = sys.argv[1]
print "image directory: " + img_dir
colors = sys.argv[2]
print "colors: " + colors
output_dir = img_dir + '/fiji_processing'
merge_output_dir = output_dir+'/merges'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
if not os.path.exists(merge_output_dir):
    os.mkdir(merge_output_dir)
# retrieve the list of images
dir_contents = os.listdir(img_dir)
print "directory contents:"
for i in dir_contents:
    print i
im_list = [i for i in dir_contents if i.lower().endswith('.tif')]
blue = False
cyan = False
green = False
yellow = False
red = False
brightfield = False
# determine which colors are present. the colors input string provided at
# script call should be the first letter of each color to be included in the
# merged image. 
if 'b' in colors:
    blue = True
if 'c' in colors:
    cyan = True
if 'g' in colors:
    green = True
if 'y' in colors:
    yellow = True
if 'r' in colors:
    red = True
if 'w' in colors:
    brightfield = True
print "colors selected:"
if blue:
    print "blue"
if cyan:
    print "cyan"
if green:
    print "green"
if yellow:
    print "yellow"
if red:
    print "red"
if brightfield:
    print "brightfield"
if yellow and green:
    print 'This script does not support overlay of green and yellow wavelengths.'
    sys.exit()
if blue and cyan:
    print 'This script does not support overlay of blue and cyan wavelengths.'
    sys.exit()
# assign wavelength strings for colors
blue_wavelength = '405'
cyan_wavelength = '447'
green_wavelength = '488'
yellow_wavelength = '515'
red_wavelength = '594'
bf_delimiter = 'brightfield'
# create a dict of dicts. each sub-dictionary will contain key:value pairs
# whose keys are a shortened version of an image's filename with the wavelength
# information removed, and the value is the full filename. the top dictionary
# will contain one of these sub-dictionaries for each wavelength as the values,
# with the keys being the channel. see mk_fname_ref docstring for more
# clarification.  also define which color to use to search for matched images 
# of other channels (first_wavelength)
color_sublists = {}
first_wavelength = ''
if blue:
    if first_wavelength == '':
        first_wavelength = 'bfp'
    color_sublists['bfp'] = mk_fname_ref(im_list,blue_wavelength)
if cyan:
    if first_wavelength == '':
        first_wavelength = 'cfp'
    color_sublists['cfp'] = mk_fname_ref(im_list,cyan_wavelength)
if green:
    if first_wavelength == '':
        first_wavelength = 'gfp'
    color_sublists['gfp'] = mk_fname_ref(im_list,green_wavelength)
if yellow:
    if first_wavelength == '':
        first_wavelength = 'yfp'
    color_sublists['yfp'] = mk_fname_ref(im_list,yellow_wavelength)
if red:
    if first_wavelength == '':
        first_wavelength = 'rfp'
    color_sublists['rfp'] = mk_fname_ref(im_list,red_wavelength)
if brightfield:
    if first_wavelength == '':
        first_wavelength = 'bf'
    color_sublists['bf'] = mk_fname_ref(im_list, bf_delimiter)
if first_wavelength == '':
    raise ValueError('No wavelengths selected.')
print "color sublists:"
print color_sublists
# get the series of images identifier strings to be processed. see mk_fname_ref
# to see how these are made identical amongst images from different
# wavelengths.
im_series = color_sublists[first_wavelength].keys()
print "im_series: "
print im_series
z_proj_dir = output_dir + '/z_projections'
merge_z_dir = output_dir + '/z_projection_merges'
if not os.path.exists(z_proj_dir):
    os.mkdir(z_proj_dir)
if blue:
    if not os.path.exists(z_proj_dir+'/blue'):
        os.mkdir(z_proj_dir + '/blue')
if cyan:
    if not os.path.exists(z_proj_dir+'/cyan'):
        os.mkdir(z_proj_dir + '/cyan')
if green:
    if not os.path.exists(z_proj_dir+'/green'):
        os.mkdir(z_proj_dir+'/green')
if red:
    if not os.path.exists(z_proj_dir+'/red'):
        os.mkdir(z_proj_dir+'/red')
if yellow:
    if not os.path.exists(z_proj_dir+'/yellow'):
        os.mkdir(z_proj_dir+'/yellow')
if not os.path.exists(merge_z_dir):
    os.mkdir(merge_z_dir)
for stage_pos in im_series:
    print 'current stage position ID: ' + stage_pos
    imps_for_comp = [None, None, None, None, None, None, None]
    if blue:
        cyan_id = color_sublists['bfp'][stage_pos]
	print 'blue image file: ' + cyan_id
        c_imp = ImagePlus(img_dir+'/'+cyan_id)
        imps_for_comp[4] = c_imp
    if cyan:
        cyan_id = color_sublists['cfp'][stage_pos]
	print 'cyan image file: ' + cyan_id
        c_imp = ImagePlus(img_dir+'/'+cyan_id)
        imps_for_comp[4] = c_imp
    if green:
        green_id = color_sublists['gfp'][stage_pos]
 	print 'green image file: ' + green_id
        g_imp = ImagePlus(img_dir+'/'+green_id)
        imps_for_comp[1] = g_imp
    if yellow:
        green_id = color_sublists['yfp'][stage_pos]
 	print 'yellow image file: ' + green_id
        g_imp = ImagePlus(img_dir+'/'+green_id)
        imps_for_comp[1] = g_imp
    if red:
        red_id = color_sublists['rfp'][stage_pos]
  	print 'red image file: ' + red_id
        r_imp = ImagePlus(img_dir+'/'+red_id)
        imps_for_comp[5] = r_imp
    if brightfield:
        bf_id = color_sublists['bf'][stage_pos]
	print 'brightfield image file: ' + bf_id
        bf_imp = ImagePlus(img_dir+'/'+bf_id)
    print ''
    if 'w' in colors:
	if len(colors) > 2:
    	    composite = RGBStackMerge.mergeChannels(imps_for_comp, True)
    	    IJ.saveAsTiff(composite, merge_output_dir + '/' + stage_pos)
    	    composite.close()
    else:
	if len(colors) > 1:
    	    composite = RGBStackMerge.mergeChannels(imps_for_comp, True)
    	    IJ.saveAsTiff(composite, merge_output_dir + '/' + stage_pos)
    	    composite.close()
    imps_for_z_comp = [None, None, None, None, None, None, None]
    if blue:
        z_cyan = maxZprojection(c_imp)
        IJ.saveAsTiff(z_cyan, z_proj_dir+'/blue/'+stage_pos)
        imps_for_z_comp[4] = z_cyan
        c_imp.close()
    if cyan:
        z_cyan = maxZprojection(c_imp)
        IJ.saveAsTiff(z_cyan, z_proj_dir+'/cyan/'+stage_pos)
        imps_for_z_comp[4] = z_cyan
        c_imp.close()
    if green:
        z_green = maxZprojection(g_imp)
        IJ.saveAsTiff(z_green, z_proj_dir+'/green/'+stage_pos)
        imps_for_z_comp[1] = z_green
        g_imp.close()
    if yellow:
        z_green = maxZprojection(g_imp)
        IJ.saveAsTiff(z_green, z_proj_dir+'/yellow/'+stage_pos)
        imps_for_z_comp[1] = z_green
        g_imp.close()
    if red:
        z_red = maxZprojection(r_imp)
        IJ.saveAsTiff(z_red, z_proj_dir+'/red/'+stage_pos)
        imps_for_z_comp[5] = z_red
        r_imp.close()
    if brightfield:
        bf_NSlices = bf_imp.getNSlices()
        if bf_NSlices > 1:
            bf_stack = bf_imp.getStack()
            desired_slice = int(bf_NSlices/2)
            bf_imp = bf_stack.getProcessor(desired_slice)
        imps_for_z_comp[3] = bf_imp
    z_composite = RGBStackMerge.mergeChannels(imps_for_z_comp, True)
    IJ.saveAsTiff(z_composite, merge_z_dir+'/'+stage_pos)

