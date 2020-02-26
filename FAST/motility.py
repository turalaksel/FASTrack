#!/usr/bin/env python

#Invitro-motility (Actin sliding via myosins) analysis python module
#Tural Aksel

import matplotlib
matplotlib.use('TkAgg')

import os
import sys
import time
import math
import re

import numpy as np
import numpy
import matplotlib.pyplot as py
import matplotlib.cm as cm
import scipy.io

import plotparams as plotparams

from numpy import ma
from scipy.ndimage import label
from scipy.ndimage.morphology import binary_fill_holes,binary_closing,binary_opening
from scipy import stats

from skimage import img_as_uint
from skimage.filters import thresholding,rank,threshold_otsu,gaussian_filter
from skimage.morphology import disk, square, rectangle, skeletonize, dilation
from skimage.morphology.watershed import watershed
from skimage.io import imread_collection

from imageio import imwrite

from scipy.optimize import leastsq
from scipy.stats    import kde

import cv2

#Global variables/structures
sqr_1   = square(1)    #Square with a radius of 1 pixel
sqr_2   = square(2)    #Square with a radius of 2 pixels
sqr_3   = square(3)    #Square with a radius of 3 pixels
disk_1  = disk(1)      #Disk   with a radius of 1 pixel
disk_2  = disk(2)      #Disk   with a radius of 2 pixels
disk_3  = disk(3)      #Disk   with a radius of 3 pixels

#Very small numbers
ZERO    = 1E-100

#Utility functions
def make_N_colors(cmap_name, N): 
    cmap = cm.get_cmap(cmap_name, N) 
    return cmap(np.arange(N))

def stack_to_tiffs(fname,frame_rate=1.0):
    '''
    Read and convert tiff stack file to individual files
    '''
    #Find the directory the tiff stack file is located
    abs_path  = os.path.abspath(fname)
    head,tail = os.path.split(abs_path)
    base,ext  = os.path.splitext(tail)
    
    #Make the new directory
    new_dir   = head+'/'+('_'.join(base.split())).replace('#','')
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)
    
    #Read all the frames
    tiff_frames  = imread_collection(fname)
    num_frames   = len(tiff_frames)
    
    f = open(new_dir+'/metadata.txt','w')
    elapsed_time_ms = 0.0
    #Write out the individual image files
    for i in range(num_frames):
        fout = new_dir+'/img_000000%03d'%(i)+'__000.tif'
        imwrite(fout,tiff_frames[i])
        #Write elapsed times
        f.write('  "ElapsedTime-ms": %d,\n'%(elapsed_time_ms))
        elapsed_time_ms += 1000*1.0/frame_rate
    f.close()

#Functions for statistical analysis
def gaussian(X,amp,mu,stdev):
    '''
    Gaussian
    '''
    return amp*np.exp(-(X-mu)**2/(2*stdev**2))

def fit_gaussian(bin_centers,bin_amps):
    '''
    Fit gaussian to a histogram
    '''
    err = lambda params: params[0]*np.exp(-(bin_centers-params[1])**2/(2*params[2]**2)) - bin_amps
    params = [50,700,100]
    best_params,success = leastsq(err,params,maxfev=1000)
    
    return best_params[0],best_params[1],best_params[2]

def fit_length_velocity(length,velocity,fil_weights,weighted=False):
    '''
    Function for fitting a length velocity curve
    '''
    myosin_density = 1.0/36.0
    Neff           = length*myosin_density
    
    weights = np.ones(len(length))
    if weighted:
        weights = fil_weights
    err = lambda params: weights*(params[0]*(1.0-(1.0-params[1])**Neff) - velocity)
    
    params = [700,0.001]
    best_params,success = leastsq(err,params,maxfev=1000)
    residuals           = np.array(err(best_params)/weights)
    return best_params[0],best_params[1],residuals, success

def length_velocity(length,max_vel,f):
    '''
    Uyeda's simple length-velocity relationship
    '''
    myosin_density = 1.0/36.0
    Neff           = length*myosin_density
    return max_vel*(1.0-(1.0-f)**Neff)

def coupling_velocity(length,max_vel,amp,tau):
    '''
    Coupling relationship with single exponential decay
    '''
    return max_vel-amp*np.exp(-length/tau)

def fit_coupling_velocity(length,velocity,fil_weights,weighted=False):
    '''
    Function for fitting a length velocity curve
    '''
    weights = np.ones(len(length))
    if weighted:
        weights = fil_weights
    err = lambda params: weights*(params[0]-params[1]*np.exp(-length/params[2]) - velocity)
    
    params = [700,200,500]
    best_params,success = leastsq(err,params,maxfev=1000)
    residuals           = np.array(err(best_params)/weights)
    return best_params[0],best_params[1],best_params[2],residuals, success

def bin_length_velocity(length,velocity,dx=100):
    '''
    Bin legth vs.velocity profile
    '''
    max_len = np.max(length)
    bin_vel = []
    for i in np.arange(0,int(max_len/dx)):
        valid    = (length > i*dx)*(length <=(i+1)*dx)
        if np.sum(valid) > 0:
            mean_len = np.mean(length[valid])
            mean_vel = np.mean(velocity[valid])
            bin_vel.append([mean_len,mean_vel])
    
    return np.array(bin_vel)

def contour2contour(contour1,contour2,fil_direction):
    '''
    Find the distance between two counters
    '''
    
    short_contour = contour1
    long_contour  = contour2
    
    if len(short_contour) > len(long_contour):
        short_contour = contour2
        long_contour  = contour1
    
    #Length of the short and long vectors
    short_len     = len(short_contour)
    long_len      = len(long_contour)
    
    
    #Multiplicate measurements
    multiplicate_measures = long_len - short_len + 1
    distance_score = 0
    for i in range(multiplicate_measures):
        long_short_diff  = long_contour[i:i+short_len,:][::fil_direction] - short_contour
        distance_length  = np.mean(np.sqrt(np.sum(long_short_diff**2,axis=1)))
        distance_score  += distance_length
    
    distance_score /= multiplicate_measures
    
    return distance_score

#Global functions

def vec_length(vec):
    '''
    Returns the euclidian distance of a numpy array
    '''
    return np.sqrt(np.sum(vec**2,axis=1))

'''
Frame link class
'''
class Link:
    def __init__(self):
        self.frame1_no        = 0
        self.frame2_no        = 0
        self.filament1_label  = 0
        self.filament2_label  = 0
        self.filament1_length = 0
        self.filament2_length = 0
        self.filament1_contour= []
        self.filament2_contour= []
        self.filament1_cm     = []
        self.filament2_cm     = []
        
        self.average_length   = 0
        self.overlap_score    = 0
        self.area_score       = 0
        self.distance_score   = 0
        
        self.fil_direction    = 1
        self.mov_direction    = 1
        
        self.dt               = 0
        self.instant_velocity = 0
        
        self.forward_link     = None
        self.reverse_link     = None
        
        self.direct_link      = False

class Path:
    def __init__(self):
        self.links           = []
        self.first_frame_no  = 0
        self.path_length     = 0
        self.ave_fil_length  = 0
        self.ave_velocity    = 0
        self.std_velocity    = 0
        self.max_velocity    = 0
        self.min_velocity    = 0
        
        self.stuck           = False
        
class Motility:
    def __init__(self):
        self.elapsed_times    = []    #Elapsed times for the movie: read from metadata.txt file generated by micro manager 1.4
        self.dt               = 1.0   #Time difference between frames in seconds(s)
        self.dx               = 80.65 #Pixel length in nm (from calibration 12/11/13)
        self.frame            = None  #Single frame
        self.frame1           = None  #First  frame
        self.frame2           = None  #Second frame
        self.frame_links      = []    #Container for the frame links
        
        self.min_velocity     = 80    #Minimum average path velocity of the filament in nm/s (Default:80 nm/s)
        self.max_velocity     = 25    #Maximum velocity of the filament in pixels/frame (Default:25)
        self.min_fil_length   = 0     #Minimum filament length
        self.max_fil_length   = 125   #Maximum filament length
        self.max_fil_width    = 25    #Maximum filament width
        
        self.max_length_dif   = 5     #Maximum length difference for the same filament in two different frames (Default:5)
        self.max_velocity_dif = 5     #Maximum velocity difference between end-end, center-center, front-front (Default:5)
        
        self.directory        = ''    #Directory where all the frames are
        self.header           = ''    #Header name for the image files
        self.tail             = ''    #Tail   name for the image files
        self.num_frames       = 0     #Total number of frames starting from 0
        
        self.norm_len_vel     = []    #Corrected length-velocity after length dependent velocity fit
        self.full_len_vel     = []    #Uncorrected (raw) length-velocity profile
        self.max_len_vel      = []    #Keeps only the maximum instantaneous velocities for the paths
        
        self.corr_lens        = []    #Correlation length information of the filaments
        
        #Paths
        self.paths            = []    #Paths constructed from frame links
        self.path_img         = None  #Path image file
        
        #Default frame widths and heights
        self.width            = 1002
        self.height           = 1004
        
        #Cutoff values
        self.overlap_score_cutoff       = 0.4
        self.log_area_score_cutoff      = 1.0
        self.dif_log_area_score_cutoff  = 0.5
        
        #Force analysis parameter
        self.force_analysis  = False
        
    def min_length_filter(self, min_filament_length):
        '''
        Applies a minimum length filter (higher than or equal). Length is in nm.
        '''
        valid_length = np.nonzero(self.full_len_vel[:,0] >= min_filament_length)[0]
        self.full_len_vel = self.full_len_vel[valid_length,:]
        
    def read_metadata(self):
        '''
        Read metadata.txt file to retrieve elapsed time information
        '''
        fname = self.directory+'/metadata.txt'
        
        #Read metadata file
        if os.path.exists(fname):
            f     = open(fname,'r')
            lines = f.readlines()
            
            f.close()
            filtered_lines = filter(lambda x:x.find('"ElapsedTime-ms"') > 0,lines)
            
            #Elapsed times
            self.elapsed_times = []
            for line in filtered_lines:
                
                m = re.search('ElapsedTime-ms":\s+(\d+),', line)
                self.elapsed_times.append(float(m.group(1)))
            
            #Elapsed times in seconds
            self.elapsed_times = 0.001*np.array(self.elapsed_times)
            
            #Sort the time array
            self.elapsed_times = np.sort(self.elapsed_times)

    def calc_persistence_len(self):
        '''
        Calculate length correlation
        '''
        self.final_corr_len     = np.zeros(1000)
        self.final_corr_weight  = np.zeros(1000)
        
        for corr_len in self.corr_lens:
            if len(corr_len) == 0:
                continue
            new_corr_len    = np.zeros(1000)
            new_corr_weight = np.zeros(1000)
            
            new_corr_len[np.arange(corr_len.shape[0],dtype=int)]     =corr_len[:,1]
            new_corr_weight[np.arange(corr_len.shape[0],dtype=int)]  =corr_len[:,0]
            
            self.final_corr_len    += new_corr_len*new_corr_weight
            self.final_corr_weight += new_corr_weight
        
        #Finalize the arrays
        valid = self.final_corr_weight > 0
        self.final_corr_len[valid] = self.final_corr_len[valid]/self.final_corr_weight[valid]
        self.final_corr_len    = self.final_corr_len[valid]
        self.final_corr_weight = self.final_corr_weight[valid]
    
    def wire_frame_links(self,depth=5):
        '''
        Wire the frame links that are not connected to each other to create contigous paths
        '''
        num_frames = len(self.frame_links)
        for d in range(1,depth+1):
            for f1 in range(len(self.frame_links)):
                possible_links = filter(lambda link:link.forward_link == None,self.frame_links[f1])
                for l1 in range(len(possible_links)):
                    link1         = possible_links[l1]
                    new_filament1 = Filament()
                    new_filament1.contour    = link1.filament2_contour
                    new_filament1.fil_length = link1.filament2_length
                    new_filament1.cm         = link1.filament2_cm
                    new_filament1.time       = link1.filament2_time
                    
                    #Average distance score per time for link 1
                    avg_distance_score_1 = possible_links[l1].distance_score/possible_links[l1].dt
                    
                    if f1+d < num_frames:
                        forward_links = filter(lambda link:link.reverse_link == None and np.sqrt(np.sum((link.filament1_cm - new_filament1.cm)**2)) < d*self.max_velocity and np.fabs(link.filament1_length - new_filament1.fil_length) < self.max_length_dif ,self.frame_links[f1+d])
                        for l2 in range(len(forward_links)):
                            link2         = forward_links[l2]
                            new_filament2 = Filament()
                            new_filament2.contour    = link2.filament1_contour
                            new_filament2.fil_length = link2.filament1_length
                            new_filament2.cm         = link2.filament1_cm
                            new_filament2.time       = link2.filament1_time 
                            
                            #Average distance score per time for link 2
                            avg_distance_score_2 = forward_links[l2].distance_score/forward_links[l2].dt
                            
                            #Average distance score for the two links
                            avg_distance_score   = 0.5*(avg_distance_score_1+avg_distance_score_2)
                            
                            #Calculate time difference
                            dt = new_filament2.time - new_filament1.time
                            
                            #Calculate similarity scores
                            overlap_score,area_score,distance_score,fil_direction,mov_direction = new_filament1.sim_score(new_filament2)
                            if np.fabs(overlap_score) > self.overlap_score_cutoff and np.log10(area_score) < self.log_area_score_cutoff and distance_score/dt < 2*avg_distance_score and distance_score/dt > 0.5*avg_distance_score:
                                new_link = Link()
                                new_link.frame1_no          = link1.frame2_no
                                new_link.frame2_no          = link2.frame1_no
                                new_link.filament1_label    = link1.filament2_label
                                new_link.filament2_label    = link2.filament1_label
                                new_link.filament1_cm       = link1.filament2_cm
                                new_link.filament2_cm       = link2.filament1_cm
                                new_link.filament1_length   = link1.filament2_length
                                new_link.filament2_length   = link2.filament1_length
                                new_link.filament1_contour  = link1.filament2_contour
                                new_link.filament2_contour  = link2.filament1_contour
                                new_link.filament1_midpoint = link1.filament2_midpoint
                                new_link.filament2_midpoint = link2.filament1_midpoint
                                new_link.filament1_time     = link1.filament2_time
                                new_link.filament2_time     = link2.filament1_time
                                
                                new_link.fil_direction    = fil_direction
                                new_link.mov_direction    = mov_direction
                                
                                new_link.overlap_score    = overlap_score
                                new_link.area_score       = area_score
                                new_link.distance_score   = distance_score
                                
                                new_link.average_length   = 0.5*(new_link.filament1_length+new_link.filament2_length)
                                new_link.instant_velocity = new_link.distance_score/dt
                                new_link.dt               = dt
                                
                                #This is not a direct connection
                                new_link.direct_link      = False
                                
                                #Add new link to frame links
                                self.frame_links[new_link.frame1_no].append(new_link)
                                
                                #Make the link connections
                                link1.forward_link = new_link
                                link2.reverse_link = new_link
                                
                                new_link.reverse_link = link1
                                new_link.forward_link = link2
    
    def read_frame_links(self):
        '''
        Read frame links saved as npy
        '''
        if not self.force_analysis and os.path.exists(self.directory+'/links.npy'):
            try:
                self.frame_links      = np.load(self.directory+'/links.npy')
            #If links.npy is the output of an old version of motility
            except ImportError:
                print'Movie analysed previously with an old version of motility. Links will be regenerated.'
                return False
            return True
        else:
            return False
    
    def reconstruct_skeleton_images(self):
        '''
        Reconstruct images from reduced skeleton representations
        '''
        #Check if paths_2D.png figure exists - if not do not construct images

        if not os.path.isfile(self.directory+'/paths_2D.png'):
            return

        #filament size ratio
        ratio = self.width/1002.0

        for i in range(len(self.frame_links)):
            #Start with a frame
            new_frame = Frame()
            
            #Prepare the processed images
            new_frame.img_skeletons  = np.zeros((self.width,self.height),dtype=np.bool)
            for link in self.frame_links[i]:
                new_frame.img_skeletons[link.filament1_contour[:,0],link.filament1_contour[:,1]] = True
        
            #Dilate skeletons for better visualization
            new_frame.img_skeletons = dilation(new_frame.img_skeletons,selem=disk(ratio*6))
            
            #Make 0 values 2 for transparency
            new_frame.img_skeletons = ma.masked_where(new_frame.img_skeletons == 0, new_frame.img_skeletons)
            
            py.figure()
            py.imshow(new_frame.img_skeletons,cmap=cm.gray,interpolation='nearest',alpha=1.0)
            
            #Plot velocities and movement direction on skeleton images
            arrow_length = ratio*1.0
            for link in self.frame_links[i]:
                velocity = link.instant_velocity
                mp_1     = link.filament1_midpoint
                mp_2     = link.filament2_midpoint
                mp_diff  = mp_2 - mp_1
                
                py.arrow(mp_1[1],mp_1[0],arrow_length*mp_diff[1],arrow_length*mp_diff[0],color='r',head_width=ratio*20,head_length=ratio*30,alpha=1.0)
                py.text(mp_1[1],mp_1[0],"%.f"%(velocity),fontsize=10,color='k',alpha=1.0)
            
            #Change x,y-ticks to um scale
            ax = py.gca()
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            
            #Skeleton and paths image filename
            skeleton_fname = self.directory+'/skeletons_%03d.png'%(i)
            paths_fname    = self.directory+'/paths_2D.png'
            py.savefig(skeleton_fname,dpi=400,transparent=True)
            py.close()
            
            #Combine with the path image
            os.system('composite -compose src-over '+skeleton_fname+' -alpha on '+paths_fname+' '+skeleton_fname)
    
    def make_forward_links(self):
        '''
        Make the forward links for frame links
        '''
        #Traverse through all frame links starting from the last frame
        #Make the forward link assignments for the frame links
        
        for i in range(len(self.frame_links)-1,-1,-1):
            for link in self.frame_links[i]:
                prev_link    = link.reverse_link
                current_link = link
                if link.forward_link == None:
                    while not prev_link == None:
                        prev_link.forward_link = current_link
                        current_link = prev_link
                        prev_link    = prev_link.reverse_link
    
    def create_paths(self):
        '''
        Create paths from frame links
        '''
        self.paths = []
        for i in range(len(self.frame_links)-1,-1,-1):
            for link in self.frame_links[i]:
                new_path     = Path()
                prev_link    = link.reverse_link
                current_link = link
                if link.forward_link == None:
                    #Add the current filament2 cm to path 
                    new_path.links.append(current_link)
                    while not prev_link == None:
                        current_link = prev_link
                        #Add the current filament2 cm to path 
                        new_path.links.append(current_link)
                        prev_link    = prev_link.reverse_link
                    
                if len(new_path.links) > 0:
                    self.paths.append(new_path)
        
        #Correct velocity values - with respect to filament direction
        for path in self.paths:
            fil_direction = path.links[0].mov_direction
            for i in range(1,len(path.links)):
                fil_direction           *= path.links[i-1].fil_direction
                mov_direction            = fil_direction*path.links[i].mov_direction
                
                path.links[i].instant_velocity = path.links[i].instant_velocity*mov_direction
                
    def path_velocities(self,num_points=1):
        '''
        Calculate average velocities over a path with num_points points
        '''
        self.full_len_vel = []
        self.max_len_vel  = []
        
        for path in self.paths:
            
            #Minimum path length should higher than num points
            if len(path.links) < num_points:
                continue
            
            #Check if the filament has moved or not
            mp_diff   = path.links[-1].filament1_midpoint - path.links[0].filament2_midpoint
            time_diff = np.fabs(path.links[-1].filament1_time - path.links[0].filament2_time) 
            dist      = np.sqrt(np.sum(mp_diff**2))
            
            #Determine whether a filament is stuck
            if self.dx*dist/time_diff < self.min_velocity:
                path.stuck        = True
            
            #Determine filament length and velocities in nm and nm/s units
            for link in path.links:
                link.instant_velocity *= self.dx
                link.average_length   *= self.dx

            #Complete path analysis
            array_vel   = np.array([np.fabs(link.instant_velocity) for link in path.links])
            array_len   = np.array([np.fabs(link.average_length)   for link in path.links])
            path_length = len(array_vel)
            
            #Smooth velocity representation
            array_smooth  = []
            
            for i in range(len(path.links)-num_points+1):
                #Average length
                ave_len   = np.mean(array_len[i:i+num_points])
                if path.stuck == True:
                    ave_vel = 0
                    std_vel = 0
                else:
                    ave_vel   = np.mean(array_vel[i:i+num_points])
                    std_vel   = np.std(array_vel[i:i+num_points])
                    
                array_smooth.append([ave_len, ave_vel, std_vel, path_length])
                self.full_len_vel.append([ave_len, ave_vel, std_vel, path_length])
            
            #Pick the max velocity along the path
            array_smooth = np.array(array_smooth)
            max_i        = np.argmax(array_smooth[:,1])
            self.max_len_vel.append(array_smooth[max_i,:])
            
        self.full_len_vel = np.array(self.full_len_vel)
        self.max_len_vel  = np.array(self.max_len_vel)
        
        #If there is no path constructed-there won't ne any velocity data
        if len(self.full_len_vel) == 0:
            return
        
        #Sort the lists based on filament length
        sort_i  = np.argsort(self.full_len_vel[:,0])
        self.full_len_vel = self.full_len_vel[sort_i,:]
        
        sort_i  = np.argsort(self.max_len_vel[:,0])
        self.max_len_vel = self.max_len_vel[sort_i,:]
        
    def make_frame_links(self):
        '''
        Make links between two adjacent frames
        '''
        #Storage for link candidates
        link_candidates_major      = []
        
        #Storage for frame links between frame1 and frame2
        new_frame_links            = []
        
        #Pick filaments in frame1 with length higher than min_fil_length
        self.frame1.filaments = filter(lambda filament:filament.fil_length > self.min_fil_length,self.frame1.filaments)
        self.frame2.filaments = filter(lambda filament:filament.fil_length > self.min_fil_length,self.frame2.filaments)
        
        #Reset filament labels
        self.frame1.reset_filament_labels()
        self.frame2.reset_filament_labels()
        
        #Time difference between frame1 and frame2
        if len(self.elapsed_times) > 0:
            self.dt      = self.elapsed_times[self.frame2.frame_no] - self.elapsed_times[self.frame1.frame_no]
            frame1_time  = self.elapsed_times[self.frame1.frame_no]
            frame2_time  = self.elapsed_times[self.frame2.frame_no]
        else:
            frame1_time  = self.frame1.frame_no*self.dt
            frame2_time  = self.frame2.frame_no*self.dt
        
        for i in range(len(self.frame1.filaments)):
            filament1              = self.frame1.filaments[i]
            
            #First create two link candidate list
            link_candidates = []
            
            #Look only the filaments that are close in space with similar lengths (within 5 pixels)
            frame2_filaments = filter(lambda filament:np.sqrt(np.sum((filament.cm - filament1.cm)**2)) < self.max_velocity and np.fabs(filament.fil_length - filament1.fil_length) < self.max_length_dif ,self.frame2.filaments)
            
            for j in range(len(frame2_filaments)):
                filament2       = frame2_filaments[j]
                
                #Calculate similarity score filament1 vs. filament2
                overlap_score,area_score,distance_score,fil_direction,mov_direction = filament1.sim_score(filament2)
                link_candidates.append([filament2.label,area_score,overlap_score,distance_score,fil_direction,mov_direction])
            
            link_candidates = np.array(link_candidates)
            num_candidates  = len(link_candidates)
            
            if num_candidates > 0:
                #Sort candidate list based on area score
                sorted_i = np.argsort(link_candidates[:,2])
                link_candidates = link_candidates[sorted_i,:]
                
                #Take log scores
                area_score_list           = link_candidates[:,1]
                log_area_score_list       = np.log10(area_score_list)
                log_area_score_diff_list  = log_area_score_list[1:] - log_area_score_list[:-1]
                overlap_score_list        = link_candidates[:,2]
                distance_score_list       = link_candidates[:,3]
                fil_direction_list        = link_candidates[:,4]
                mov_direction_list        = link_candidates[:,5]
                
                #Acceptance criteria
                if  np.fabs(overlap_score_list[0]) > self.overlap_score_cutoff and log_area_score_list[0] < self.log_area_score_cutoff and ((num_candidates > 1 and log_area_score_diff_list[0] >= self.dif_log_area_score_cutoff) or (num_candidates == 1)):
                    new_link = Link()
                    new_link.frame1_no          = filament1.frame_no
                    new_link.frame2_no          = filament2.frame_no
                    new_link.filament1_label    = filament1.label
                    new_link.filament2_label    = filament2.label
                    new_link.filament1_cm       = filament1.cm
                    new_link.filament2_cm       = filament2.cm
                    new_link.filament1_length   = filament1.fil_length
                    new_link.filament2_length   = filament2.fil_length
                    new_link.filament1_contour  = filament1.contour
                    new_link.filament2_contour  = filament2.contour
                    new_link.filament1_midpoint = filament1.midpoint
                    new_link.filament2_midpoint = filament2.midpoint
                    new_link.filament1_time     = frame1_time
                    new_link.filament2_time     = frame2_time
                    
                    new_link.fil_direction    = fil_direction_list[0]
                    new_link.mov_direction    = mov_direction_list[0]
                    
                    new_link.overlap_score    = overlap_score_list[0]
                    new_link.area_score       = area_score_list[0]
                    new_link.distance_score   = distance_score_list[0]
                    
                    new_link.average_length   = 0.5*(filament1.fil_length+filament2.fil_length)
                    new_link.instant_velocity = new_link.distance_score/self.dt
                    new_link.dt               = self.dt
                    
                    #This is a direct connection
                    new_link.direct_link      = True
                    
                    new_link.reverse_link     = filament1.reverse_link
                    
                    #Pick only one confident link per filament
                    filament1.forward_link    = new_link
                    
                    if filament2.reverse_link == None:
                        filament2.reverse_link = new_link
                    elif new_link.overlap_score < filament2.reverse_link.overlap_score:
                        #Pick up the existing reverse link for filament2
                        prev_fil_label = int(filament2.reverse_link.filament1_label)
                        prev_filament1 = self.frame1.filaments[prev_fil_label]
                        prev_filament1.forward_link = None
                        
                        #Assign the reverse link for filament2
                        filament2.reverse_link      = new_link
        
        #Put all the new links in the frame link container
        for i in range(len(self.frame1.filaments)):
            filament1 = self.frame1.filaments[i]
            if not filament1.forward_link == None:
                new_frame_links.append(filament1.forward_link)
        
        #Finally add the links to the big list
        self.frame_links.append(new_frame_links)
    
    def plot_2D_path_data(self, num_points, extra_fname=None):
        '''
        Write path velocities and the image
        '''
        #Filament size ratio
        ratio = self.width/1002.0

        #Path data
        self.path_data  = []
        self.path_stats = []
        
        #Prepare the processed images
        self.path_img      = np.nan*np.ones((self.width,self.height),dtype=np.bool)
        
        #Filter paths - pick only paths longer-equal to num-points
        filtered_paths = filter(lambda x:len(x.links) >= num_points,self.paths)
        
        #If there are no filtered paths - continue
        if len(filtered_paths) == 0:
            return
            
        #Get colors from a colormap
        path_colors = make_N_colors('Accent',len(filtered_paths))
        
        py.figure(2000)
        py.imshow(self.path_img,cmap=cm.gray,alpha=1.0)
        
        py.figure(2001)
        py.imshow(self.path_img,cmap=cm.gray,alpha=1.0)
        
        #Go through each path and write velocities
        for i in range(len(filtered_paths)):
            path = filtered_paths[i]
            mp_mean     = np.mean(np.array([[link.filament1_midpoint[1],link.filament1_midpoint[0]] for link in path.links[::-1]]),axis=0)
            
            len_array   = np.array([np.fabs(link.average_length) for link in path.links[::-1]])
            vel_array   = np.array([np.fabs(link.instant_velocity) for link in path.links[::-1]])
            
            first_frame = path.links[-1].frame1_no
            path_length = len(path.links)
            
            stuck       = path.stuck
            
            #Keep data in the arrays
            self.path_data.append([first_frame ,stuck,vel_array])
            self.path_stats.append([first_frame,stuck,path_length,np.mean(len_array),np.mean(vel_array),np.std(vel_array)])
            
            #Mean velocity
            mean_velocity = np.fabs(np.mean(vel_array))
            if stuck:
                mean_velocity = 0
                
            for j in range(len(path.links)):
                mp_x1    = path.links[j].filament1_midpoint[1]
                mp_y1    = path.links[j].filament1_midpoint[0]
                mp_x2    = path.links[j].filament2_midpoint[1]
                mp_y2    = path.links[j].filament2_midpoint[0]
                
                #Plot the arrows connecting path
                py.figure(2000)
                py.arrow(mp_x2,mp_y2,mp_x1-mp_x2,mp_y1-mp_y2,color=path_colors[i],head_width=ratio*5,head_length=ratio*10,alpha=1.0)
                
                py.figure(2001)
                py.arrow(mp_x2,mp_y2,mp_x1-mp_x2,mp_y1-mp_y2,color=path_colors[i],head_width=ratio*5,head_length=ratio*10,alpha=1.0)
            
            #Write mean velocity for the path
            py.figure(2001)
            py.text(mp_mean[0],mp_mean[1],"%.f"%(mean_velocity),fontsize=10,color='k')
        
        #Convert path data to numpy array
        self.path_stats = np.array(self.path_stats)
        
        py.figure(2000)
        
        ax = py.gca()
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        py.savefig(self.directory+'/paths_2D.png',dpi=400,transparent=False)
        
        py.figure(2001)
        ax = py.gca()
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        
        if not extra_fname == None:
            py.figure(2001)
            py.savefig(extra_fname+'_2D.png',dpi=400,transparent=False)
        
        py.close('all')
        
        return self.path_data
    
    def write_path_data(self, extra_fname = None):
        #Write paths info
        f = open(self.directory+'/paths.txt','w')
        
        for data in self.path_data:
            f.write('%8d\t%8d'%(data[0],data[1]))
            for vel in data[2]:
                f.write('\t%8.f'%(vel))
            f.write('\n')
        f.close()
        
        if not extra_fname == None:
            #Write paths info in extra file
            f = open(extra_fname+'.txt','w')
            for data in self.path_data:
                f.write('%8d\t%8d'%(data[0],data[1]))
                for vel in data[2]:
                    f.write('\t%8.f'%(vel))
                f.write('\n')
            f.close()
    
    
    def plot_1D_path_data(self, min_path_length = 10, include_stuck = False, max_velocity = 1000, show_length = False, extra_fname=None):
        #Prepare two figures with and without velocity values
        py.figure(1000,figsize=(20,10))
        py.figure(1001,figsize=(20,10))
        
        #Number of bins
        nbins = 20
        
        py.figure(1000)
        valid = np.arange(len(self.path_stats[:,1]))
        hist_valid = np.nonzero((self.path_stats[:,1] == 0)*(self.path_stats[:,2] > min_path_length))[0]
        if not include_stuck:
            valid = np.nonzero(self.path_stats[:,1] == 0)[0]
            
        mean_vel = np.hstack(self.path_stats[hist_valid,4])
        std_vel  = np.hstack(self.path_stats[hist_valid,5])
        
        #Change maximum velocity if the preassigned value is less than the maximum of the data
        if max_velocity < np.max(mean_vel):
            max_velocity = np.max(mean_vel)
        
        #Histogram edges
        xedges = np.linspace(0,    max_velocity,nbins)
        yedges = np.linspace(0,0.5*max_velocity,nbins)
        
        data = np.vstack((mean_vel,std_vel))
        k = kde.gaussian_kde(data)
        X, Y = np.mgrid[0:max_velocity+1:max_velocity/nbins, 0:0.5*max_velocity+1:max_velocity/(2*nbins)]
        Z = k(np.vstack([X.flatten(), Y.flatten()]))
        
        x_bin_size = max_velocity/nbins
        y_bin_size = max_velocity/(2*nbins)
        py.pcolormesh(X, Y, Z.reshape(X.shape))
        CS = py.contour(X+0.5*x_bin_size, Y+0.5*y_bin_size, Z.reshape(X.shape))
        
        #Determine the peak mean values
        Z = Z.reshape(X.shape)
        high_point  = CS.levels[-1]
        high_valid  = np.nonzero(Z >= high_point)
        pvm,pvs,psm,pss = np.mean(X[high_valid[0],1]),np.std(X[high_valid[0],1]),np.mean(Y[1,high_valid[1]]),np.std(Y[1,high_valid[1]])
        
        #Calculate 1D projections - save them
        mean_vel_project = np.hstack((np.vstack(X[:,1]),np.vstack(np.sum(Z,axis=1))))
        std_vel_project  = np.hstack((np.vstack(Y[1,:]),np.vstack(np.sum(Z,axis=0))))
        
        #Show length information on 1D-path-plots
        if show_length:
            #Plot marker-legend
            for i in range(11):
                py.plot(max_velocity*i/11.0,0.5*max_velocity*10.0/11.0,markersize=i*3,mec=(1,1,1,0),mfc=(1, 1, 1, 0.7),marker='o')
            
            #Plot individual points
            for point in range(len(hist_valid)):
                valid_point = hist_valid[point]
                py.plot(self.path_stats[valid_point,4],self.path_stats[valid_point,5],markersize=30.0*self.path_stats[valid_point,3]/5000,mec=(1,1,1,0),mfc=(1, 1, 1, 0.7),marker='o')
        
        np.savetxt(self.directory+'/paths_mean_vel_projection.dat',mean_vel_project)
        np.savetxt(self.directory+'/paths_std_vel_projection.dat',std_vel_project)
        
        #Determine mean-values
        mvm = np.sum(Z*X)/np.sum(Z)
        mvs = np.sqrt(np.sum(Z*(X-mvm)**2)/np.sum(Z))
        
        msm = np.sum(Z*Y)/np.sum(Z)
        mss = np.sqrt(np.sum(Z*(Y-msm)**2)/np.sum(Z))
        
        #Peak value
        peak_stats = [pvm,pvs,psm,pss,mvm,mvs,msm,mss] 
        
        path_vel_mean = np.sum(self.path_stats[valid,4]*self.path_stats[valid,2])/np.sum(self.path_stats[valid,2])
        path_std_mean = np.sum(self.path_stats[valid,5]*self.path_stats[valid,2])/np.sum(self.path_stats[valid,2])
        
        path_vel_std  = np.sqrt(np.sum((self.path_stats[valid,4]- path_vel_mean)**2*self.path_stats[valid,2])/np.sum(self.path_stats[valid,2]))
        path_std_std  = np.sqrt(np.sum((self.path_stats[valid,5]- path_std_mean)**2*self.path_stats[valid,2])/np.sum(self.path_stats[valid,2]))
        
        py.xlim([0,max_velocity])
        py.ylim([0,max_velocity*0.5])
        py.grid(color='white',linewidth=5,linestyle='--',alpha=0.5)
        
        py.savefig(self.directory+'/paths_1D.png',dpi=200,transparent=True)
        
        #Set tick label sizes
        ax = py.gca()
        
        py.setp(ax.get_xticklabels() , fontsize=55, visible=True)
        py.setp(ax.get_yticklabels() , fontsize=55, visible=True)
        
        ax.set_xlabel("Mean path velocity (nm/s)",fontsize=55)
        ax.set_ylabel("Std Path velocity  (nm/s)",fontsize=55)
        
        #Save path stats
        np.savetxt(self.directory+'/paths_stats.dat',self.path_stats,header='first-frame stuck path-length mean-length mean-vel mean-std')
        
        if not extra_fname == None:
            py.figure(1001)
            valid      = np.nonzero(self.path_stats[:,2] > min_path_length)[0]
            hist_valid = np.nonzero((self.path_stats[:,1] == 0)*(self.path_stats[:,2] > min_path_length))[0]
            if not include_stuck:
                valid      = np.nonzero((self.path_stats[:,1] == 0)*(self.path_stats[:,2] > min_path_length))[0]
            
            mean_vel = np.hstack(self.path_stats[hist_valid,4])
            std_vel  = np.hstack(self.path_stats[hist_valid,5])
            data = np.vstack((mean_vel,std_vel))
            k = kde.gaussian_kde(data)
            X, Y = np.mgrid[0:max_velocity+1:max_velocity/nbins, 0:0.5*max_velocity+1:max_velocity/(2*nbins)]
            Z = k(np.vstack([X.flatten(), Y.flatten()]))
            
            x_bin_size = max_velocity/nbins
            y_bin_size = max_velocity/(2*nbins)
            py.pcolormesh(X, Y, Z.reshape(X.shape))
            CS = py.contour(X+0.5*x_bin_size, Y+0.5*y_bin_size, Z.reshape(X.shape))
            
            #Determine the peak mean values
            Z = Z.reshape(X.shape)
            high_point  = CS.levels[-1]
            high_valid  = np.nonzero(Z >= high_point)
            pvm,pvs,psm,pss = np.mean(X[high_valid[0],1]),np.std(X[high_valid[0],1]),np.mean(Y[1,high_valid[1]]),np.std(Y[1,high_valid[1]])
            
            #Calculate 1D projections - save them
            mean_vel_project = np.hstack((np.vstack(X[:,1]),np.vstack(np.sum(Z,axis=1))))
            std_vel_project  = np.hstack((np.vstack(Y[1,:]),np.vstack(np.sum(Z,axis=0))))
            
            #Show length information on 1D-path-plots
            if show_length:
                #Plot marker-legend
                for i in range(10):
                    py.plot(max_velocity*i/11.0,0.5*max_velocity*10.0/11.0,markersize=i*3,mec=(1,1,1,0),mfc=(1, 1, 1, 0.7),marker='o')
                
                #Plot individual points
                for point in range(len(hist_valid)):
                    valid_point = hist_valid[point]
                    py.plot(self.path_stats[valid_point,4],self.path_stats[valid_point,5],markersize=30.0*self.path_stats[valid_point,3]/5000,mec=(1,1,1,0),mfc=(1, 1, 1, 0.7),marker='o')
            
            np.savetxt(extra_fname+'_mean_vel_projection.dat',mean_vel_project)
            np.savetxt(extra_fname+'_std_vel_projection.dat',std_vel_project)
            
            #Save path stats
            np.savetxt(extra_fname+'_stats.dat',self.path_stats,header='first-frame stuck path-length mean-length mean-vel mean-std')
        
            #Determine mean-values
            mvm = np.sum(Z*X)/np.sum(Z)
            mvs = np.sqrt(np.sum(Z*(X-mvm)**2)/np.sum(Z))
            
            msm = np.sum(Z*Y)/np.sum(Z)
            mss = np.sqrt(np.sum(Z*(Y-msm)**2)/np.sum(Z))
            
            #Peak value
            peak_stats = [pvm,pvs,psm,pss,mvm,mvs,msm,mss] 
            
            #For histgram use only
            len_hist_data = len(self.path_stats[hist_valid,4])
            
            path_vel_mean = np.sum(self.path_stats[valid,4]*self.path_stats[valid,2])/np.sum(self.path_stats[valid,2])
            path_std_mean = np.sum(self.path_stats[valid,5]*self.path_stats[valid,2])/np.sum(self.path_stats[valid,2])
            
            path_vel_std  = np.sqrt(np.sum((self.path_stats[valid,4]- path_vel_mean)**2*self.path_stats[valid,2])/np.sum(self.path_stats[valid,2]))
            path_std_std  = np.sqrt(np.sum((self.path_stats[valid,5]- path_std_mean)**2*self.path_stats[valid,2])/np.sum(self.path_stats[valid,2]))
            
            py.xlim([0,max_velocity])
            py.ylim([0,max_velocity*0.5])
            py.grid(color='white',linewidth=5,linestyle='--',alpha=0.5)
            
            #Set tick label sizes
            ax = py.gca()
            
            py.setp(ax.get_xticklabels() , fontsize=55, visible=True)
            py.setp(ax.get_yticklabels() , fontsize=55, visible=True)
            
            #Set x and y-labels
            ax.set_xlabel("Mean path velocity (nm/s)",fontsize=55)
            ax.set_ylabel("Std  path velocity (nm/s)",fontsize=55)
            
            #Save figure
            py.savefig(extra_fname+'_1D.png',dpi=200,transparent=True)
            
        
        py.close('all')
        return self.path_stats,peak_stats
    
    def process_frame_links(self,num_points=5):
        '''
        Process frame links to extract length/velocity information
        '''
        #Make the forward links for frame links
        self.make_forward_links()
        
        #Connect disconnected paths
        self.wire_frame_links()
        
        #Create paths from frame links
        self.create_paths()
        
        #Calculate velocities from paths
        self.path_velocities(num_points)
        
    def make_movie(self,extra_fname=None):
        '''
        Make a movie file using ffmpeg
        '''
        
        #Check if paths_2D.png figure exists - if not do not make the movie
        if not os.path.isfile(self.directory+'/paths_2D.png'):
            return
        
        #Get current working directory
        cwd = os.getcwd()
        os.chdir(self.directory)
        
        avicommand = "avconv -y -r 1 -i skeletons_%03d.png -r 1 filament_tracks.avi"
        os.system(avicommand)
        
        #Copy the movie to output directory
        if not extra_fname == None:
            copy_command = "cp filament_tracks.avi "+extra_fname+"filament_tracks.avi"
            os.system(copy_command)
        
        #Delete all skeleton file
        os.system('rm -f skeletons_*.png')
        os.chdir(cwd)
    
    def read_frame(self,num_frame,force_read = False):
        '''
        Read single frame
        '''
        #Read the frame
        
        print 'Reading frame: %d'%(num_frame)
        self.frame           = Frame()
        self.frame.directory = self.directory
        self.frame.header    = self.header
        self.frame.tail      = self.tail
        self.frame.frame_no  = num_frame
        
        #If already exists load the saved array file
        filament_file = self.directory+'/filXYs%03d.npy'%num_frame
        if not force_read and os.path.isfile(filament_file):
            self.frame.read_filXYs()
            self.frame.filXY2filaments()
            return 1
        
        #If the file does not exist, exit the loop
        if not self.frame.read_frame(num_frame):
            sys.exit('File not found!')
        
        self.frame.low_pass_filter()

        self.frame.entropy_clusters()
   
        self.frame.filter_islands()
    
        self.frame.skeletonize_islands()
    
        self.frame.filaments2filXYs()
    
    def save_frame(self):
        '''
        Save the filament-xy's of frame
        '''
        self.frame.save_filXYs()
    
    def load_frame1(self, frame_no):
        '''
        Load frame 1
        '''
        self.frame1           = Frame()
        self.frame1.directory = self.directory
        self.frame1.header    = self.header
        self.frame1.tail      = self.tail
        self.frame1.frame_no  = frame_no
        
        #Read the filament XYs and reconstruct filaments
        self.frame1.read_filXYs()
        self.frame1.filXY2filaments()
        
    def load_frame2(self, frame_no):
        '''
        Load frame 2
        '''
        self.frame2           = Frame()
        self.frame2.directory = self.directory
        self.frame2.header    = self.header
        self.frame2.tail      = self.tail
        self.frame2.frame_no  = frame_no
        
        #Read the filament XYs and reconstruct filaments
        self.frame2.read_filXYs()
        self.frame2.filXY2filaments()
            
    def write_length_velocity(self, header='',extra_fname=None):
        '''
        Write length vs.velocity in a txt file
        '''
    
        np.savetxt(self.directory+'/'+header+'full_length_velocity.txt', self.full_len_vel)
        np.savetxt(self.directory+'/'+header+'max_length_velocity.txt', self.max_len_vel)
        
        if not extra_fname == None:
            np.savetxt(extra_fname+'full_length_velocity.txt',self.full_len_vel)
            np.savetxt(extra_fname+'max_length_velocity.txt',self.max_len_vel)
    
    def save_links(self):
        '''
        Save the results in working directory
        '''
        np.save(self.directory+'/links.npy' ,self.frame_links)
    
    def load_links(self):
        '''
        Load the results that have been saved earlier
        '''
        self.frame_links = np.load(self.directory+'/links.npy')
    
    def plot_length_velocity(self,header='',extra_fname=None,max_vel=2400, max_length= 10000,nbins=30, min_points=2, min_path_length = 5, weighted=True, percent_tolerance=500, print_plot=True, minimal_plot=False, maxvel_color='b', plot_xlabels = True, plot_ylabels = True, square_plot = True, plot_length_f=False, fit_f = 'exp',dpi_plot=200):
        '''
        Plot length-vs-velocity data
        '''
        
        #Pick only filaments shorter than maximum length
        valid = np.nonzero(self.full_len_vel[:,0] < max_length)[0]
        self.full_len_vel = self.full_len_vel[valid,:]
        
        #Calculate top 1%, top 5%, mean velocity excluding stuck
        tolerance_data = []
        tolerance_list = [2.5,5,10,20,40,80]
        valid_points   = np.nonzero(self.full_len_vel[:,1] >= 0)[0]
        for filter_value in tolerance_list[::-1]:
            filtered_data = self.full_len_vel[valid_points,:]
            if len(valid_points) > 10:
                stuck     = np.nonzero(filtered_data[:,1] == 0)[0]
                non_stuck = np.nonzero(filtered_data[:,1] != 0)[0]
                
                if len(non_stuck) > 0:
                    #Get filament velocities
                    fil_vel    = filtered_data[non_stuck,1]
                    
                    #Mean velocity of the moving filaments only
                    mean_vel_m = np.mean(fil_vel)
                    std_vel_m  = np.std(fil_vel)
                    
                    velocities_sorted = np.sort(fil_vel)[::-1]
                    top_1_num         = int(np.ceil(0.01*len(velocities_sorted)))
                    top_5_num         = int(np.ceil(0.05*len(velocities_sorted)))
                    
                    #Number of data points
                    num_filter_points = len(fil_vel)

                    top_1_velocity    = np.mean(velocities_sorted[:top_1_num])
                    top_5_velocity    = np.mean(velocities_sorted[:top_5_num])
                    
                    tolerance_data.append([filter_value*2, num_filter_points, top_1_velocity,top_5_velocity,mean_vel_m,std_vel_m])
                
                else:
                    tolerance_data.append([filter_value*2, 0.0, 0.0 ,0.0, 0.0, 0.0])
            
            else:
                tolerance_data.append([filter_value*2, 0.0, 0.0 ,0.0, 0.0, 0.0])
            #Valid points for the next iteration
            valid_points  = np.nonzero(self.full_len_vel[:,2] <= filter_value/100.0*self.full_len_vel[:,1])[0]
        
        #Convert tolerance data to array format
        tolerance_data = np.array(tolerance_data)
        
        #Calculate percent stuck
        percent_stuck = 100.0*np.sum(self.full_len_vel[:,1] == 0)/len(self.full_len_vel[:,0])
        
        #Text fontsize
        text_font_size = 30
        
        #Valid points - only moving, satisfying smooth movement tolerance
        valid_filtered      = np.nonzero((self.full_len_vel[:,1] > 0)*(self.full_len_vel[:,2] <= percent_tolerance/100.0*self.full_len_vel[:,1])*(self.full_len_vel[:,3] >= min_path_length))[0]
        num_points_filtered = len(valid_filtered)
        
        #Valid points - only moving
        valid_mobile        = np.nonzero((self.full_len_vel[:,1] > 0)*(self.full_len_vel[:,3] >= min_path_length))[0]
        num_points_mobile   = len(valid_mobile)
        
        #Valid points including stuck
        valid_all           = np.nonzero(self.full_len_vel[:,3] >= min_path_length)[0]
        num_points_all      = len(valid_all)
        
        #Valid points - only stuck filaments
        valid_stuck         = np.nonzero((self.full_len_vel[:,3] >= min_path_length)*(self.full_len_vel[:,1] == 0))[0]
        num_points_stuck    = len(valid_stuck) 
        
        if num_points_filtered < min_points:
            #There is no frame-link
            print 'Warning: There is not enough velocity data! - %d points'%(num_points_t)
            return -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
        
        #Statistics data
        MVEL_filtered      = np.mean(self.full_len_vel[valid_filtered,1])
        MVEL               = np.mean(self.full_len_vel[valid_mobile,1])
        MVIS               = np.mean(self.full_len_vel[valid_all,1])
        mean_len_filtered  = np.mean(self.full_len_vel[valid_filtered,0])
        mean_len_mobile    = np.mean(self.full_len_vel[valid_mobile,0])
        mean_len_stuck     = np.mean(self.full_len_vel[valid_stuck,0])
        mean_len_all       = np.mean(self.full_len_vel[valid_all,0])
        
        #Get filament length and velocities
        fil_len     = self.full_len_vel[valid_filtered,0]
        fil_vel     = self.full_len_vel[valid_filtered,1]
        
        #Length histogram bin edges-centers parameters
        l_bin_edges       = np.linspace(0,max_length*1E-3, nbins)
        l_bin_centers     = 0.5*(l_bin_edges[:-1]+l_bin_edges[1:])
        l_tick_locs       = l_bin_centers
        l_tick_labels     = ['%d'%round(x) for x in l_bin_centers]
        
        #Get length histogram values excluding stuck filaments
        l_bin_counts, l_bin_locs = np.histogram(self.full_len_vel[valid_filtered,0],bins=l_bin_edges,normed=False)
        
        #Velocity histogram bin edges-centers parameters
        v_bin_edges   = np.linspace(0,max_vel,nbins)
        v_bin_width   = v_bin_edges[1] - v_bin_edges[0]
        v_bin_centers = 0.5*(v_bin_edges[:-1]+v_bin_edges[1:])
        v_tick_locs   = [round(x) for x in v_bin_centers]
        v_tick_labels = ['%.0f'%round(x) for x in v_bin_centers]
        
        #Get velocity histogram values excluding stuck filaments
        v_bin_counts, v_bin_locs = np.histogram(self.full_len_vel[valid_filtered,1],bins=v_bin_edges,normed=True)
        
        #Determine top 1 and 5 percent velocities
        velocities_sorted = np.sort(fil_vel)[::-1]
        top_1_num         = np.ceil(0.01*len(velocities_sorted))
        top_5_num         = np.ceil(0.05*len(velocities_sorted))
        
        top_1_velocity    = np.mean(velocities_sorted[:top_1_num])
        top_5_velocity    = np.mean(velocities_sorted[:top_5_num])
        
        #Fit the data
        fil_len_digitized = np.digitize(1E-3*fil_len,l_bin_locs)
        fil_weights       = 1.0/l_bin_counts[fil_len_digitized-1]
        mean_vel_u,mean_vel_amp,mean_vel_tau,residuals,success = fit_coupling_velocity(fil_len,fil_vel,fil_weights,weighted=weighted)
        std_u                                                  = np.sqrt(np.mean(residuals**2))
        
        #Find the treshold filament length for plateou 
        bound_prob    = coupling_velocity(fil_len,0.0,-1.0,mean_vel_tau)
        plateu_valid  = np.nonzero(bound_prob > 0.95)[0]
        plateu_vel    = fil_vel[plateu_valid]
        mean_plateu   = np.mean(plateu_vel)
        std_plateu    = np.std(plateu_vel)
        roof_plateu   = mean_plateu+2.0*std_plateu
        
        #Find max probability
        max_index_t = np.argmax(v_bin_counts)
    
        #Peak velocity - moving
        peak_vel_t  = v_bin_centers[max_index_t]
        
        #Perform length velocity analysis on the maximum velocity data
        #Minimum path length is required for maximum lenght velocity analysis
        max_valid = np.nonzero((self.max_len_vel[:,1] > 0)*(self.max_len_vel[:,2] <= percent_tolerance/100.0*self.max_len_vel[:,1])*(self.max_len_vel[:,3] >= min_path_length))[0]
        
        #Function fitted can be either exponential or uyeda's equation
        if fit_f == 'exp':
            max_vel_u,max_vel_amp,max_vel_tau,residuals,success = fit_coupling_velocity(self.max_len_vel[max_valid,0],self.max_len_vel[max_valid,1],np.ones(len(self.max_len_vel[max_valid,0])),weighted=False)
            std_u                                               = np.sqrt(np.mean(residuals**2))
        elif fit_f == 'uyeda':
            max_vel_u,max_vel_r,residuals,success               = fit_length_velocity(self.max_len_vel[max_valid,0],self.max_len_vel[max_valid,1],np.ones(len(self.max_len_vel[max_valid,0])),weighted=False)
            std_u                                               = np.sqrt(np.mean(residuals**2))
        
        #Estimated curve from filament length data
        exp_len                               = np.linspace(np.min(fil_len),15000,1000)
        if fit_f == 'exp':
            exp_vel                           = coupling_velocity(exp_len,max_vel_u,max_vel_amp,max_vel_tau)
        elif fit_f == 'uyeda':
            exp_vel                           = length_velocity(exp_len,max_vel_u,max_vel_r)
        
        #Determine tolerance word
        if percent_tolerance == 500:
            tolerance_string = 'none'
        else:
            tolerance_string = str(percent_tolerance)        

        if print_plot:
            if minimal_plot:
                text_font_size = 55
                linewidth      = 10
                
                #Determine 0.99 probability length: myosin density=1/36
                
                if fit_f == 'exp':
                    length_f = max_vel_tau
                elif fit_f == 'uyeda':
                    length_f = np.log(0.01)/np.log(1-max_vel_r)*36.0
                
                #Plot square image
                x,y = plotparams.get_figsize(1080)
                if square_plot:
                    py.figure(0,figsize=(y,y))
                else:
                    py.figure(0,figsize=(x,y))
                
                py.plot(1E-3*fil_len,fil_vel,'.',markersize=10 ,color='gray')
                py.plot(1E-3*self.max_len_vel[max_valid,0],self.max_len_vel[max_valid,1],marker='^',markersize=10,mec=maxvel_color, mfc=maxvel_color,linestyle='None')
                
                #Do not plot plateaue line if user doesn't want data to be fitted
                if not fit_f == 'none': 
                    py.plot(1E-3*exp_len,exp_vel,'k-',linewidth=linewidth, alpha=0.7)
                
                py.plot(1E-3*exp_len,np.ones(len(exp_len))*top_5_velocity,'k-.'    ,linewidth=linewidth)
                
                if plot_length_f:
                    #Plot L-0.5
                    py.plot(1E-3*np.array([length_f,length_f]),[0,top_5_velocity],linestyle='dashed',color='k')
                    py.text(length_f*1E-3+0.1,10,'%.1f'%(length_f*1E-3),fontsize=text_font_size,color='k')
                
                #X/Y limits
                py.ylim([0,max_vel])
                py.xlim([0,max_length*1E-3+1.0])
                
                #Place the histogram ticks
                ax = py.gca()
                vel_ticks = ax.get_yticks()
                ax.set_yticks(vel_ticks[::2])
                ax.set_yticklabels(vel_ticks[::2]*1E-3)
                
                len_ticks = ax.get_xticks()
                ax.set_xticks(len_ticks[::2])
                ax.set_xticklabels([int(x) for x in len_ticks[::2]])
                
                #Change the padding between the axis and labels
                ax.tick_params(pad=10)
                
                #Set tick label size
                py.setp(ax.get_xticklabels(), fontsize=text_font_size, visible=plot_xlabels)
                py.setp(ax.get_yticklabels(), fontsize=text_font_size, visible=plot_ylabels)

            else:
                
                # definitions for the axes
                left,   width  = 0.1, 0.5
                bottom, height = 0.1, 0.5
                
                #Left boundary for velocity histograms
                left_h1        = left+width
                left_h2        = left_h1   +0.15
                
                #Lower boundary for the length histogram 
                bottom_v1      = bottom+0.27
                
                rect_scatter   = [left        , bottom_v1     , width  , height  ]
                rect_tolerance = [left_h1+0.01, bottom+0.02   , 0.29   , 0.24    ]
                rect_histy1    = [left_h1     , bottom_v1     , 0.15   , height  ]
                rect_histy2    = [left_h2     , bottom_v1     , 0.15   , height  ]
                rect_histx1    = [left        , bottom+0.02   , width  , 0.25    ]
                
                #Plot the data and fit
                py.figure(0,figsize=plotparams.get_figsize(1200))
                axScatter    = py.axes(rect_scatter)
                
                axHisty1     = py.axes(rect_histy1)
                axHisty2     = py.axes(rect_histy2)
                
                axHistx1     = py.axes(rect_histx1)
                
                axTolerance1 = py.axes(rect_tolerance)
                axTolerance2 = axTolerance1.twinx()
                
                #Plot tolerance data - top 5% and mean-velocity excluding stuck filaments
                max_tol_vel  = np.max(tolerance_data[:,3:5])
                min_tol_vel  = np.min(tolerance_data[:,3:5])
                
                #Plot the tolerance information
                #First plot lines
                axTolerance2.plot(tolerance_data[:,0],tolerance_data[:,3],color='k',linestyle='--',marker='.',linewidth=5,markersize=15)
                axTolerance2.plot(tolerance_data[:,0],tolerance_data[:,4],color='k',linestyle='-' ,marker='.',linewidth=5,markersize=15)
                axTolerance2.set_xscale('symlog')
                axTolerance1.set_xscale('symlog')
                
                tol_ymin = min_tol_vel-100
                tol_ymax = max_tol_vel+100
                tol_diff = max_tol_vel-min_tol_vel+200
                
                axTolerance2.set_ylim([tol_ymin, tol_ymax])
                axTolerance2.set_xlim([5,200])
                
                #Plot tolerance legend
                axTolerance2.plot([6,9],[tol_ymin+0.25*tol_diff,tol_ymin+0.25*tol_diff],color='k',linestyle='--',linewidth=5)
                axTolerance2.text(10,tol_ymin+0.20*tol_diff,r"%s"%('TOP5%') ,fontsize=text_font_size,color='k')
                
                axTolerance2.plot([6,9],[tol_ymin+0.10*tol_diff,tol_ymin+0.1*tol_diff],color='k',linestyle='-',linewidth=5)
                axTolerance2.text(10,tol_ymin+0.05*tol_diff,r"%s"%('Mean Velocity') ,fontsize=text_font_size,color='k')
                
                axTolerance2.set_xticks(tolerance_data[:,0])
                axTolerance2.set_xticklabels(['*']+[int(x) for x in tolerance_data[1:,0]])
                axTolerance1.set_xlabel('% Tolerance',fontsize=text_font_size,labelpad=20)
                
                #Set velocity labels
                vel_ticks = axTolerance2.get_yticks()
                axTolerance2.set_yticks(vel_ticks[1::2])
                axTolerance2.set_yticklabels(vel_ticks[1::2]*1E-3)
                
                ylim      = axTolerance2.get_ylim()
                tol_diff  = ylim[1]-ylim[0]
                
                axTolerance2.text(300,ylim[1]+0.1*tol_diff,r'$x10^3$',fontsize=25)
                
                py.setp(axTolerance2.get_yticklabels() , fontsize=text_font_size, visible=True)
                py.setp(axTolerance2.get_xticklabels() , fontsize=text_font_size, visible=True)
                
                py.setp(axTolerance1.get_yticklabels() , fontsize=text_font_size, visible=False)
                py.setp(axTolerance1.get_xticklabels() , fontsize=text_font_size, visible=True)
            
                #Plot Length histogram
                l_bin_counts, l_bin_locs, l_patches = axHistx1.hist(1E-3*self.full_len_vel[valid_filtered,0],bins=l_bin_edges,normed=False,orientation='vertical',color='gray')
                max_prob_l     = np.max(l_bin_counts)
            
                #Plot histogram with all mobile filaments
                v_bin_counts, v_bin_locs, v_patches = axHisty2.hist(self.full_len_vel[valid_mobile,1],bins=v_bin_edges,normed=True,orientation='horizontal',color='gray')
            
                #Get the y-limit for y-histogram-2
                max_prob_a = axHisty2.get_xlim()[1]
                
                #Plot histogram without stuck filaments
                v_bin_counts, v_bin_locs, v_patches = axHisty1.hist(self.full_len_vel[valid_filtered,1],bins=v_bin_edges,normed=True,orientation='horizontal',color='gray')
                
                #Get the y-limit for y-histogram-1
                max_prob_t = axHisty1.get_xlim()[1]
                
                #Plot the data and fit
                axScatter.plot(1E-3*fil_len,fil_vel,'.',markersize=5 ,color='gray')
                axScatter.plot(1E-3*self.max_len_vel[max_valid,0],self.max_len_vel[max_valid,1],marker='^',markersize=5,mec=maxvel_color, mfc=maxvel_color,linestyle='None')
                
                #Do not plot plateaue line if user doesn't want data to be fitted
                if not fit_f == 'none':
                    axScatter.plot(1E-3*exp_len,exp_vel,'k-',alpha=0.7)
                
                axScatter.plot(1E-3*exp_len,np.ones(len(exp_len))*top_5_velocity,'k--'    ,linewidth=5)
                
                axHisty1.plot([0,max_prob_t],np.ones(2)*MVEL_filtered  ,color='k',linestyle='-',linewidth=5)
                
                axHisty2.plot([0,max_prob_a],np.ones(2)*MVEL  ,'k-'    ,linewidth=5)
                
                axHistx1.plot([mean_len_filtered*1E-3,mean_len_filtered*1E-3],[0,max_prob_l],'k-'     ,linewidth=5)
                
                #X/Y limits
                axScatter.set_ylim([0,max_vel])
                axScatter.set_xlim([0,max_length*1E-3])
                
                #Place the histogram ticks
                vel_ticks = axScatter.get_yticks()
                vel_ticks = vel_ticks[::2]
                axScatter.set_yticks(vel_ticks)
                axScatter.set_yticklabels(vel_ticks*1E-3)
                
                len_ticks = axScatter.get_xticks()
                axScatter.set_xticks(len_ticks[:-1])
                axScatter.set_xticklabels([int(x) for x in len_ticks[:-1]])
                
                axHistx1.set_xticks(len_ticks[:-1])
                axHistx1.set_xticklabels([int(x) for x in len_ticks[:-1]])
                
                axHisty1.set_yticks(vel_ticks)
                axHisty1.set_yticklabels(vel_ticks*1E-3)
                
                axHisty2.set_yticks(vel_ticks)
                axHisty2.set_yticklabels(vel_ticks*1E-3)
                
                axScatter.text(0,max_vel,r'$x10^3$',fontsize=25)
                
                #X/Y limits
                axScatter.set_ylim([0,max_vel])
                axScatter.set_xlim([0,max_length*1E-3])
                
                #Format tick labels
                py.setp(axHisty1.get_yticklabels() , fontsize=text_font_size, visible=False)
                py.setp(axHistx1.get_xticklabels() , fontsize=text_font_size, visible=True)
                py.setp(axScatter.get_xticklabels(), fontsize=text_font_size, visible=False)
                py.setp(axScatter.get_yticklabels(), fontsize=text_font_size, visible=True)
                
                axHisty1.set_ylim([0,max_vel])
                axHisty2.set_ylim([0,max_vel])
                
                axHistx1.set_xlim([0,max_length*1E-3])
                
                #Remove the xtick labels for the histogram
                axHisty1.ticklabel_format(style='sci', axis='x', scilimits=(-5,5))
                axHisty2.ticklabel_format(style='sci', axis='x', scilimits=(-5,5))
                
                py.setp(axHisty1.get_xticklabels(), visible=False)
                py.setp(axHisty2.get_xticklabels(), visible=False)
                py.setp(axHisty2.get_yticklabels(), visible=False)
                py.setp(axHistx1.get_yticklabels(), visible=False)
    
                #Labels for the histogram plots
                axTolerance2.set_ylabel(r'Velocity (nm/s)',labelpad=20,fontsize=text_font_size)
                axScatter.set_ylabel(r'Velocity (nm/s)',labelpad=20,fontsize=text_font_size)
                axHistx1.set_xlabel(r'Actin filament length ($\mu m$)',labelpad=20,fontsize=text_font_size)
                
                axHisty1.text(0.1*max_prob_t,1.1*max_vel,'Filtered'  ,fontsize=text_font_size)
                axHisty2.text(0.1*max_prob_a,1.1*max_vel,'Unfiltered',fontsize=text_font_size)
                
                #Get maximum length for writing velocity on the figure
                axScatter.plot(max_length*1E-3*np.array([1,2])/15.0,[2150/2400.0*max_vel,2150/2400.0*max_vel],'k--'     ,linewidth=5)
                axScatter.text(max_length*1E-3*2.1/15.0            ,2100/2400.0*max_vel,r"%.f$^{TOP5\%%}$"%(top_5_velocity) ,fontsize=text_font_size,color='k')
                
                #Do not write plateaue velocity if user doesn't want data to be fitted
                if not fit_f == 'none':
                    axScatter.plot(max_length*1E-3*np.array([6,7])/15.0,[2150/2400.0*max_vel,2150/2400.0*max_vel],'k-'    ,linewidth=10)
                    axScatter.text(max_length*1E-3*7.1/15.0            ,2100/2400.0*max_vel ,r"%.f$^{PLATEAU}$"%(max_vel_u) ,fontsize=text_font_size,color='k')
                
                axHisty1.text(0.1*max_prob_t,1900/2400.0*max_vel,r"%.f$^{MVEL_{%s}}$"%(MVEL_filtered,tolerance_string)  ,fontsize=text_font_size,color='k')
                axHisty2.text(0.1*max_prob_a,1900/2400.0*max_vel,r"%.f$^{MVEL}$"%(MVEL)       ,fontsize=text_font_size,color='k')
                axHisty2.text(0.15*max_prob_a,1600/2400.0*max_vel,r"%.f$^{\%%STUCK}$"%(percent_stuck)   ,fontsize=text_font_size,color='k')
                axHistx1.text(mean_len_filtered*1E-3,max_prob_l*0.5,r"%.3f$^{<FIL-LENGTH>}$"%(mean_len_filtered*1E-3) ,fontsize=text_font_size,color='k')
            
            py.savefig(self.directory+'/'+header+'length_velocity.png',dpi=dpi_plot,transparent=False)
            if not extra_fname == None:
                py.savefig(extra_fname,dpi=dpi_plot,transparent=False)
            py.close()
        
        #List to be returned
        if not fit_f == 'none':
            return_list  = top_5_velocity, percent_stuck, MVEL, MVEL_filtered, max_vel_u, MVIS, mean_len_stuck, mean_len_filtered, mean_len_mobile, mean_len_all,num_points_filtered
        else:
            return_list  = top_5_velocity, percent_stuck, MVEL, MVEL_filtered, -1       , MVIS, mean_len_stuck, mean_len_filtered, mean_len_mobile, mean_len_all,num_points_filtered

        return return_list

    def plot_correlation_profile(self,extra_fname = None):
        py.figure(4,figsize=plotparams.get_figsize(1080))
        
        array_corr_len    = np.arange(len(self.final_corr_len))*self.dx
        array_corr_weight = np.arange(len(self.final_corr_weight))*self.dx 
        
        valid             = np.nonzero((self.final_corr_len > 0.7)*(array_corr_len<=1500))
        slope, intercept, r_value, p_value, std_err = stats.linregress(array_corr_len[valid],1.0*self.final_corr_len[valid])
        
        #The distance that leads to 0.7 correlation 
        length_0_7       = np.round((0.7 - intercept)/slope)
        mean_corr_1500   = np.mean(1.0*self.final_corr_len[valid])
        
        py.subplot(211)
        py.plot(array_corr_len   ,1.0*self.final_corr_len,'bo')
        py.plot(array_corr_len,array_corr_len*slope+intercept,'r-',linewidth=5)
        
        py.text(1500,0.9,r'l$_{0.7}$: %d'%(length_0_7),fontsize=50)
        py.xlim(0,3000)
        py.ylim(0.7,1.0)
        py.ylabel(r'c($\Delta$ nm)')
        
        py.subplot(212)
        py.plot(array_corr_weight,1.0*self.final_corr_weight)
        py.xlim(0,self.max_fil_length*self.dx)
        
        py.xlabel(r'$\Delta$ nm')
        py.ylabel(r'weight (#)')
        py.xlim(0,self.max_fil_length*self.dx)
        
        py.savefig(self.directory+'/correlation_length.png',dpi=200)
        if not extra_fname == None:
            py.savefig(extra_fname,dpi=200)
        py.close()
        
        return length_0_7, mean_corr_1500
    
class Frame:
    def __init__(self):
        #Global parameters - mainly cutoff-thresholding values
        self.frame_no                 = 0                           #Frame number
        self.window_island            = 15                          #The radius of the disk used for local entropy measurement (default:15)
        self.disk_win                 = disk(self.window_island)    #Disk used for local entropy calculation
        
        #The images to store
        self.img                      = None                        #Image in matrix form
        self.img_filaments            = None                        #Full  filaments' image
        self.img_skeletons            = None                        #Skeletonized filaments' image
        
        #Frame dimension properties
        self.width                    = 1002                        #Width of the image in pixels
        self.height                   = 1004                        #Height of the image in pixels
        
        #Conatainer for the islands and backward links
        self.backward_links           = []                          #Frame links in backward direction
        self.islands                  = []                          #Array that keeps the islands
        self.filaments                = []                          #Array that keeps the accepted filaments
        self.filXYs                   = []                          #Array that keeps X-Y positions of the filaments
        
        #File path properties
        self.directory                = ''                          #Directory for the frames
        self.header                   = ''                          #Header name for the tiff images
        self.tail                     = ''                          #Tail   name for the tiff images
        
        #Filament counter(label)
        self.filament_counter         = 0
    
    def reset_filament_labels(self):
        '''
        Resets filament labels
        '''
        for i in range(len(self.filaments)):
            self.filaments[i].label = i
        
    def read_frame(self, frame_no):
        '''
        Read the image file corresponding to a frame
        '''
        fname         = self.directory+'/'+self.header+'%03d'%frame_no+'_'+self.tail+'_000.tif' 
        self.frame_no = frame_no 
        
        #Use openCV to read 16bit images
        if not os.path.isfile(fname):
            return False
        
        self.img      = cv2.imread(fname,cv2.IMREAD_GRAYSCALE)
        
        #Convert to 8bits
        self.img      = img_as_uint(self.img)
        
        #Get the dimensions of the frame
        self.width,self.height = self.img.shape
        
        return True
    
    def filXY2filaments(self):
        '''
        Reconstruct filaments array from filXY
        '''
        self.filaments = []
        
        #Filament counter for filament label
        fil_counter    = 0
        
        for filXY,width,density,midpoint in self.filXYs:
            filament             = Filament()
            filament.frame_no    = self.frame_no
            filament.contour     = filXY
            filament.fil_width   = width
            filament.fil_density = density
            filament.label       = fil_counter
            filament.midpoint    = midpoint
            
            #Calculate filament properties
            filament.calc_props()
            
            #Add filament to list
            self.filaments.append(filament)
            
            #Pudate filament label
            fil_counter += 1
    
    
    def reconstruct_skeleton_images(self):
        #Prepare the processed images
        self.img_skeletons  = np.zeros((self.width,self.height),dtype=np.bool)
        
        for filament in self.filaments:
            self.img_skeletons[filament.contour[:,0],filament.contour[:,1]] = True
        

    def reconstruct_filament_images(self):
        '''
        Reconstruct images from reduced filament representations
        '''
        #Prepare the processed images
        self.img_filaments       = np.zeros((self.width,self.height),dtype=np.uint16)
    
        for filament in self.filaments:
            x_corrected = filament.xy[0]+filament.island.x_min
            y_corrected = filament.xy[1]+filament.island.y_min
            self.img_filaments[x_corrected,y_corrected] = filament.img_reduced[filament.xy_norm]
    
    def save_filament_img(self):
        '''
        Save skeleton image
        '''
        py.figure()
        py.imshow(self.img_filaments,cmap=cm.gray)
        
        #Change x,y-ticks to um scale
        ax = py.gca()
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        py.savefig(self.directory+'/filaments_%03d.png'%(self.frame_no),dpi=200)
        py.close()
    
    def calc_fil_corr_funcs(self):
        '''
        Calculate filament length-correlation profiles
        '''
        for filament in self.filaments:
            filament.correlation_function()
    
    def save_filXYs(self):
        '''
        Save contours of filaments
        '''
        #Save the filament-contours as npy file
        
        filament_file = self.directory+'/filXYs%03d'%self.frame_no
        np.save(filament_file,self.filXYs)
        
    def read_filXYs(self):
        '''
        Load the filXYs for filament reconstruction
        '''
        filament_file = self.directory+'/filXYs%03d.npy'%self.frame_no
        self.filXYs   = np.load(filament_file)
        
    def check_picture_quality(self):
        '''
        Check picture quality
        If the maximum value in difference image is less than 500, picture is of bad quality
        '''
        
        #Determine the 95-percentile intensity within a radius of 15 pixels
        img_1 = rank.percentile(self.img,self.disk_win,p0=0.95)
        
        #Determine the 5-percentile intensity within a radius of 15 pixels
        img_0 = rank.percentile(self.img,self.disk_win,p0=0.05)
        
        #Subtract the 5-percentile from 95-percentile map
        self.img_diff = img_1 - img_0
        
        #Maximum of the difference value
        max_diff = np.max(self.img_diff)
        relative_contrast = 1.0*np.max(self.img_diff)/np.max(img_1)
        
        if relative_contrast > 0.7 and max_diff > 1000:
            return 'good'
        else:
            return 'bad'
    
    def low_pass_filter(self,sigma=2):
        '''
        Low pass filter to remove high-frequency noise
        '''
        self.img = gaussian_filter(self.img,sigma=sigma)
        
    def entropy_clusters(self):
        '''
        Create entropy clusters
        '''
        
        #Determine the 95-percentile intensity within a radius of 15 pixels
        img_1 = rank.percentile(self.img,self.disk_win,p0=0.95)
        
        #Determine the 5-percentile intensity within a radius of 15 pixels
        img_0 = rank.percentile(self.img,self.disk_win,p0=0.05)
        
        #Subtract the 5-percentile from 95-percentile map
        self.img_diff                        = img_1 - img_0
        
        #Cutoff value is the mean intensity of the background subtracted image
        self.cutoff_diff                     = np.mean(self.img_diff)
        
        #Define a mask for the pixels with intensity higher than cut-off value
        self.mask_diff                       = self.img_diff > self.cutoff_diff
        
        #Label each separate area on the mask with a different number
        self.labels_island,self.num_island   = label(self.mask_diff)
        
        self.img_water                       = watershed(self.mask_diff,self.labels_island,mask=self.mask_diff)
        
        #Array to keep all the entropy islands
        self.islands                         = []
        
        #Prepare the entropy island objects
        for i in range(1,self.num_island+1):
            xy                     = np.nonzero(self.img_water == i)
            
            #If it is a single pixel island - ignore  
            if len(xy[0]) < 2:
                continue
            
            new_island             =  Island()
            new_island.reduce_image(xy,self.img)
            
            #Assign the frame for the island
            new_island.frame = self
            
            #Append island to the list
            self.islands.append(new_island)

    def filter_function(self,island):
        '''
        Filter function to filter islands
        '''
        area_constraint      = island.area      > self.window_island**2
        return area_constraint
    
    def filter_islands(self):
        '''
        Filter islands using the filter function
        '''
        self.islands = filter(self.filter_function,self.islands)
    
    def skeletonize_islands(self):
        '''
        Skeletonize the islands
        '''
        #Filament labels
        self.filament_counter = 0
        
        #First decompose to filaments
        [island.decompose_to_filaments()     for island in self.islands]
        
        #Second filter filaments in the islands
        [island.filter_filaments()           for island in self.islands]
        
        #Third skeletonize all the filaments in the islands
        [island.skeletonize_filaments()      for island in self.islands]
        
        #Forth remove crossing filaments
        [island.remove_crossing_filaments()  for island in self.islands]
        
        #Finally add these filaments to Frame's filament list
        self.filaments = [island.filaments for island in self.islands]
        self.filaments = [item for sub in self.filaments for item in sub]
    
    def filaments2filXYs(self):
        '''
        Convert filaments to filXY - that keeps the positions of the pixels
        '''
        self.filXYs = []
        for filament in self.filaments:
            self.filXYs.append([filament.contour,filament.fil_width,filament.fil_density,filament.midpoint])
            
class Island:
    def __init__(self):
        self.area                     = 0       #Area of the island
        self.xy                       = []      #Island x-y coordinates
        self.xy_norm                  = [[],[]] #Normalized x-y coordinates
        self.x_min                    = []      #Limits of x-axis
        self.y_min                    = []      #Limits of y-axis
        self.x_dim                    = 0       #Reduced x-dim
        self.y_dim                    = 0       #Reduced y-dim
        self.img_reduced              = None    #Reduced full image
        self.filaments                = []      #Container that keeps filament objects
        
        self.min_filament             = 4       #Minimum filament length accepted in the analysis (default:3)
        self.window_island            = 15      #The radius of the disk used for local entropy measurement (default:15)
        
        self.frame                    = None    #Pointer to the frame that the Island belongs to
    
    def reduce_image(self, xy, img):
        '''
        Reduce the image of a filament
        '''
        
        #Get the area
        self.area   =  len(xy[0])
        
        #Get the x,y coordinates
        self.xy     = xy
        self.x_min  =  np.min(self.xy[0])
        self.y_min  =  np.min(self.xy[1])
        self.x_dim  =  np.max(self.xy[0]) - np.min(self.xy[0])
        self.y_dim  =  np.max(self.xy[1]) - np.min(self.xy[1])
        
        #Find out the normalized x-y dimensions
        self.xy_norm    = [[],[]]
        self.xy_norm[0] = self.xy[0]-self.x_min
        self.xy_norm[1] = self.xy[1]-self.y_min
        
        #Get the reduced image
        self.img_reduced               = np.zeros((self.x_dim+1,self.y_dim+1))
        self.img_reduced[self.xy_norm] = img[self.xy]
    
    def decompose_to_filaments(self):
        '''
        Decompose entropy island to filaments
        '''
        
        #Pick only pixels with integer intensity values
        valid         = self.img_reduced > 0
        
        #Determine the Otsu threshold value
        cutoff             = threshold_otsu(self.img_reduced[valid])
        
        #Filament in coarse representation
        self.fil_reduced = self.img_reduced > cutoff
        self.img_fil     = self.fil_reduced*self.img_reduced
        
        #Label the filaments
        fil_labels, fil_features = label(self.fil_reduced)
        
        #In a cluster there may be more than a single filament - each cluster corresponds to a filament
        fine_clusters = watershed(self.fil_reduced,fil_labels,mask=self.fil_reduced)
        
        #Start with an empty list of filaments
        self.filaments = []
        
        for i in range(1,fil_features+1):
            
            xy_bool       = fine_clusters == i
            xy            = np.nonzero(xy_bool)
            
            #If it is a single pixel island - ignore  
            if len(xy[0]) < 2:
                continue
            
            new_filament        = Filament()
            
            #Assign the label
            new_filament.label  = self.frame.filament_counter 
            
            #Assign the current to island to the filament
            new_filament.island = self
            
            #Shrink the size of the filament image
            new_filament.reduce_image(xy)
            
            #Density of the filament in terms of intensity
            new_filament.fil_density = np.sum(1.0*self.img_reduced*xy_bool)/np.sum(xy_bool)
            
            #Fill in the holes in a filament
            new_filament.img_reduced = binary_fill_holes(new_filament.img_reduced)
            
            #Binary opening-closing needed for getting rid of extra branches
            new_filament.img_reduced = binary_closing(new_filament.img_reduced,structure=disk_1)
            
            #Add new filament to the list
            self.filaments.append(new_filament)
            
            #Increment filament counter
            self.frame.filament_counter += 1
    
    def filter_function(self,filament):
        '''
        Filter function to be used to filter bad filaments
        '''
        size_constraint    = np.sum(filament.img_reduced) > 10
        x_constraint       = (filament.x_min+self.x_min > 5) and (filament.x_min+self.x_min+filament.x_dim < self.frame.width )
        y_constraint       = (filament.y_min+self.y_min > 5) and (filament.y_min+self.y_min+filament.y_dim < self.frame.height)
        
        return size_constraint and x_constraint and y_constraint
    
    def filter_filaments(self):
        '''
        Filter the filaments using the filter function
        '''
        self.filaments = filter(self.filter_function,self.filaments)
    
    def skeletonize_filaments(self):
        '''
        Seperate routine for skeletonizing the filaments
        '''
        #Skeletonize
        [filament.make_skeleton() for filament in self.filaments]
        
        #Filter out filaments with length = 0
        self.filaments = filter(lambda fil:len(fil.contour) > 0,self.filaments)
        
        #Calculate filament properties
        [filament.calc_fil_stats() for filament in self.filaments]
        
    def remove_crossing_filaments(self):
        '''
        Remove crossing filaments
        '''
        self.filaments = filter(lambda filament:filament.num_tips == 2,self.filaments)
    
class Filament:
    def __init__(self):
        self.frame_no      = 0         #The number of the frame in which the filament is located in the current configuration
        self.label         = 0         #Filament label
        
        #Contour
        self.contour       = []        #The x-y positions of the filament pixels in ordered form (beginning to end/end to beginning)
        self.coarse        = []        #Coarse N point representation of the filament
        self.cm            = []        #Center of mass
        self.midpoint      = []        #Filament midpoint
        
        #Images
        self.edge          = 5         #Extra padding in the reduced img representation
        self.img_reduced   = []        #Reduced filament image
        self.img_skeleton  = []        #Skeletonized image
        
        #Tips
        self.tips          = []        #The coordinates of the tips
        self.num_tips      = 0         #Number of tips
        
        #Filament properties
        self.fil_length    = 0         #Filament length in pixel number
        self.fil_density   = 0         #Filament intensity density
        self.fil_area      = 0         #Filament area
        self.fil_width     = 0         #Filamnent width
        self.end2end       = 0         #Filament end-to-end distance
        
        #Filament links
        self.next_filament = None      #Next filament
        self.pre_filament  = None      #Previous filament
        
        #The Island that the filament belongs to
        self.island        = None      #Island that keeps the filament
    
        #Correlation function for the filament
        self.corr_len      = None      #Correlation length profile to calculate persistence length
        
        #Link to a filament in the next and previous frame
        self.forward_link  = None
        self.reverse_link  = None
        
        #Elapsed time
        self.time          = 0
        
    def reduce_image(self, xy):
        '''
        Reduce the image of a filament
        '''
        
        #Get the x,y coordinates
        self.xy     =  xy
        self.x_min  =  np.min(self.xy[0])
        self.y_min  =  np.min(self.xy[1])
        self.x_dim  =  np.max(self.xy[0]) - np.min(self.xy[0])
        self.y_dim  =  np.max(self.xy[1]) - np.min(self.xy[1])
        
        #Find out the normalized x-y dimensions
        self.xy_norm    = [[],[]]
        self.xy_norm[0] = self.xy[0]-self.x_min+self.edge
        self.xy_norm[1] = self.xy[1]-self.y_min+self.edge
        
        #Get the reduced image
        self.img_reduced               = np.zeros((self.x_dim+1+2*self.edge,self.y_dim+1+2*self.edge),dtype=np.uint16)
        self.img_reduced[self.xy_norm] = True
    
    def find_tips(self):
        '''
        Find the tips of a skeletonized filament
        '''
        neighbours    = rank.pop(self.img_skeleton,sqr_3,mask=self.img_skeleton)
        tips          = np.nonzero(neighbours*self.img_skeleton==2)
        if len(tips[0]) == 0:
            self.num_tips  = 0
            return np.array([])
        self.tips     = np.hstack((np.vstack(tips[0]),np.vstack(tips[1])))
        self.num_tips = len(self.tips)
        return self.tips
    
    def make_skeleton(self):
        #Skeletonize the image
        self.img_skeleton = skeletonize(self.img_reduced)
        
        #Find the tips
        self.find_tips()
        
        #Continue if only the number of tips is more than 2
        if self.num_tips ==2:
        
            #Remove bad pixels
            self.remove_bad_pixels()
            
            #Make the links in a filament
            self.make_links()
        
    def sim_score(self,fil_other):
        '''
        Similarity score of the first filament to another one
        '''
        short_current = False
        if len(self.contour) < len(fil_other.contour):
            short_current = True
        
        #Difference vectors
        contour1_diff = self.contour[1:] - self.contour[:-1]
        contour2_diff = fil_other.contour[1:]  - fil_other.contour[:-1]
        
        #Length of the short and long vectors
        contour1_len  = len(self.contour)
        contour2_len  = len(fil_other.contour)
        
        #Short contour length
        short_con_len = min((contour1_len,contour2_len))
        
        #Short filament length
        short_fil_len = min((self.fil_length,fil_other.fil_length))
        
        #Multiplicate measurements
        multiplicate_measures = abs(contour1_len - contour2_len) + 1
        
        #Keep the direction of the current filament with respect to each other and move direction
        fil_direction = 1
        mov_direction = 1
        
        #Calculate overlap score
        overlap_score  = 0
        
        for i in range(multiplicate_measures):
            if short_current:
                len1           = np.sum(np.sqrt(np.sum(contour2_diff[i:i+short_con_len-1,:]**2,axis=1)))
                len2           = np.sum(np.sqrt(np.sum(contour1_diff**2,axis=1)))
                if len1 > 0 and len2 > 0:
                    overlap_score += np.sum(contour2_diff[i:i+short_con_len-1,:]*contour1_diff)
            else:
                len1           = np.sum(np.sqrt(np.sum(contour1_diff[i:i+short_con_len-1,:]**2,axis=1)))
                len2           = np.sum(np.sqrt(np.sum(contour2_diff**2,axis=1)))
                if len1 > 0 and len2 > 0:
                    overlap_score += np.sum(contour1_diff[i:i+short_con_len-1,:]*contour2_diff)
        
        overlap_score /= 1.0*(multiplicate_measures*short_fil_len)
        
        #Find the directione
        if overlap_score > 0:
            fil_direction =  1
        else:
            fil_direction = -1
        
        #Calculate area and distance scores
        area_score     = 0
        distance_score = 0
        move_score     = 0

        for i in range(multiplicate_measures):
            if short_current:
                contour2_1_diff  = fil_other.contour[i:i+short_con_len,:][::fil_direction] - self.contour
                dot_prod         = contour2_1_diff[:-1,:]*contour1_diff
                cross_prod       = contour2_1_diff[:-1,:]*contour1_diff[:,[1,0]]
                cross_prod       = np.fabs(cross_prod[:,1]-cross_prod[:,0])
                contour1_diff_len= np.mean(np.sqrt(np.sum(contour1_diff**2,axis=1)))
            else:
                contour2_1_diff  = fil_other.contour[::fil_direction] - self.contour[i:i+short_con_len,:]
                dot_prod         = contour2_1_diff[:-1,:]*contour1_diff[i:i+short_con_len-1,:]
                cross_prod       = contour2_1_diff[:-1,:]*contour1_diff[i:i+short_con_len-1,[1,0]]
                cross_prod       = np.fabs(cross_prod[:,1]-cross_prod[:,0])
                contour1_diff_len= np.mean(np.sqrt(np.sum(contour1_diff[i:i+short_con_len-1,:]**2,axis=1)))
                
            area_score      += np.sum(cross_prod)
            distance_length  = np.mean(np.sqrt(np.sum(contour2_1_diff[:-1,:]**2,axis=1)))
            distance_score  += distance_length
            
            if distance_length > 0 and contour1_diff_len > 0:
                move_score      += np.mean(np.sum(dot_prod,axis=1))/(distance_length*contour1_diff_len)
        
        #ZERO added for taking logarithm
        area_score      = area_score/(short_fil_len*multiplicate_measures) + ZERO
        distance_score /= multiplicate_measures
        move_score     /= multiplicate_measures
        
        #Filament move direction
        if move_score > 0:
            mov_direction =  1
        else:
            mov_direction = -1

        return overlap_score, area_score, distance_score, fil_direction, mov_direction
    
    #Remove strongly connected paths - fortunately not very often used (suitable for tip systems only)
    def remove_bad_pixels(self):
        tip_s   = self.tips[0,:]
        tip_e   = self.tips[1,:]
        bad_fil = True
        while bad_fil == True:
            #Neighborhood information of all pixels
            nb_all = rank.pop(self.img_skeleton,disk_1,mask=self.img_skeleton) 
            
            #Find only pixels that are connected to 3 neighboring pixels on the faces
            nb_3   = np.nonzero(nb_all*self.img_skeleton==4)
            
            if len(nb_3[0]) > 0:
                bad_1 = [nb_3[0][0],nb_3[1][0]]
                self.img_skeleton[bad_1[0]][bad_1[1]] = 0
                new_tips = self.find_tips()
                for i in range(self.num_tips):
                    new_tip = new_tips[i,:]
                    if np.all(new_tip != tip_s) and np.all(new_tip != tip_e):
                        self.img_skeleton[new_tip[0]][new_tip[1]] = 0
            else:
                bad_fil = False
    
    #N point representation of the filament
    def N_point(self,N=3):
        points_floor = np.array([np.floor(x) for x in np.linspace(0,self.num_contour_pixels-1,N)],dtype=int)
        points_ceil  = np.array([np.ceil(x) for x in np.linspace(0,self.num_contour_pixels-1,N)],dtype=int)
        
        self.coarse  = 0.5*(1.0*self.contour[points_floor,:]+1.0*self.contour[points_ceil,:]) 
        
        return self.coarse
    
    #calcualte filament length
    def calc_fil_length(self):
        
        self.num_contour_pixels = len(self.contour[:,0])
        dist = self.contour[1:,:] - self.contour[:-1,:]
        self.fil_length = np.sum(np.sqrt(np.sum(dist**2,axis=1)))
    
    #Calculate filament stats
    def calc_fil_stats(self):
        #Calculate filament length
        self.calc_fil_length()
        
        #Calculate center of mass
        self.cm           = np.mean(1.0*self.contour,axis=0)
        
        #Filament midpoint
        self.midpoint     = self.contour[len(self.contour)/2-1]
        
        #Determine filament area
        self.fil_area   = np.sum(self.img_reduced)
        
        #Calculate the filament width
        self.fil_width  = 0.0
        if self.fil_length > 0:
            self.fil_width = 1.0*self.fil_area/self.fil_length
        
        #Prepare coarse representation
        self.N_point()
        
    #Calculate filament properties from contour
    def calc_props(self):
        '''
        Calculate filament properties from contour
        '''
        #Filament length
        self.calc_fil_length()
        
        #Calculate center of mass
        self.cm           = np.mean(self.contour,axis=0)
        
        #Filament midpoint
        self.midpoint     = self.contour[len(self.contour)/2-1]
        
        #Prepare coarse representation
        self.N_point()
    
    def correlation_function(self,P=3):
        '''
        Calculate the correlation function to determine persistence length
        '''
        #Tangent vectors
        tan_vecs      = self.contour[:-P,:] - self.contour[P:,:]
        
        #Length vectors
        len_vecs      = vec_length(tan_vecs)
        num_vecs      = len(tan_vecs)
        self.corr_len = [[num_vecs,1.0]]
        
        for i in range(1,len(tan_vecs)):
            if np.sum(len_vecs[:-i]*len_vecs[i:] == 0) == 0:
                self.corr_len.append([num_vecs-i,np.mean(np.sum(tan_vecs[:-i,:]*tan_vecs[i:,:],axis=1)/(len_vecs[:-i]*len_vecs[i:]))])
        
        #If there is invalid number in the array discard the correlation profile
        if np.sum(np.isnan(self.corr_len)) > 0:
            self.corr_len = []
        self.corr_len = np.array(self.corr_len)
        
    #Filament linked-list for a filament with two-tips 
    def make_links(self):
        self.contour            = []
        self.num_contour_pixels = int(np.sum(self.img_skeleton))
        img_c                   = self.img_skeleton.copy()
        tip_s                   = self.tips[0,:]
        tip_e                   = self.tips[1,:]
        self.contour.append(tip_s)
        
        for i in range(self.num_contour_pixels-1):
            img_c[tip_s[0],tip_s[1]] = 0
            
            if   img_c[tip_s[0]+1,tip_s[1]]   == 1:
                new_tip_s = [tip_s[0]+1,tip_s[1]]
            elif img_c[tip_s[0]  ,tip_s[1]+1] == 1:
                new_tip_s = [tip_s[0],tip_s[1]+1]
            elif img_c[tip_s[0]+1,tip_s[1]+1] == 1:
                new_tip_s = [tip_s[0]+1,tip_s[1]+1]
            elif img_c[tip_s[0]-1,tip_s[1]  ] == 1:
                new_tip_s = [tip_s[0]-1,tip_s[1]]
            elif img_c[tip_s[0]  ,tip_s[1]-1] == 1:
                new_tip_s = [tip_s[0] ,tip_s[1]-1]
            elif img_c[tip_s[0]-1,tip_s[1]-1] == 1:
                new_tip_s = [tip_s[0]-1,tip_s[1]-1]
            elif img_c[tip_s[0]+1,tip_s[1]-1] == 1:
                new_tip_s = [tip_s[0]+1,tip_s[1]-1]
            elif img_c[tip_s[0]-1,tip_s[1]+1] == 1:
                new_tip_s = [tip_s[0]-1,tip_s[1]+1]
            
            tip_s = new_tip_s
            self.contour.append(tip_s)
        self.contour = np.array(self.contour)
        
        #Correct the indices of contour in the context of full image
        
        #Determine the offsets
        offset_x     = self.x_min + self.island.x_min - self.edge
        offset_y     = self.y_min + self.island.y_min - self.edge
        
        self.contour = self.contour + np.array([offset_x,offset_y])
        
        return self.contour
