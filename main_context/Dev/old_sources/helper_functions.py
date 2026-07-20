#!/usr/bin/env python
# helper functions to be used for beam dynamics simulations
# Updated: 23.06.22 Haroon Rafique STFC ISIS Accelerator Division
# Requires cpymad tfs-pandas python libraries

import os
from sys import exit
import shutil
import subprocess

import math
import random
import numpy as np
import pandas as pd
from math import log10, floor
from itertools import combinations
from scipy.constants import c, m_p, e

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch, Rectangle

from cpymad.madx import Madx
from cpymad.madx import Sequence
from cpymad.madx import SequenceMap
from cpymad.types import Constraint
import tfs

########################################################################
# Matplotlib global plotting parameters
########################################################################
plt.rcParams['figure.figsize'] = [8.0, 5.0]
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 200

plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 14

plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

plt.rcParams['font.size'] = 10
plt.rcParams['legend.fontsize'] = 8

plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['lines.markersize'] = 5

########################################################################
# Used to calculate initial beam parameters for a given proton energy
########################################################################
class MADX_Proton_Beam_Parameters:                
    mass = 938.272E6 # in eV
    energy = -1. # in eV   
    beta = -1.
    gamma = -1.
    total_energy = -1.
    momentum = -1.

    def __init__(self, energy):
        self.energy = energy
        self.total_energy = self.get_total_energy()
        self.gamma = self.get_gamma()
        self.beta = self.get_beta()
        self.momentum = self.get_momentum()
        self.rigidity = self.get_rigidity()

    def get_total_energy(self): return (self.energy + self.mass)
    def get_gamma(self): return (self.total_energy / self.mass)
    def get_beta(self): return(np.sqrt( 1 - (1/self.gamma)**2 ))
    def get_momentum(self): return(self.gamma * self.mass * self.beta)
    def get_rigidity(self): return (self.momentum/299792458)
    
    def print_beam(self):
        print('M_proton = ', round_sig(self.mass/1E6) , 'MeV')
        print('Energy = ', round_sig(self.energy/1E9) , 'GeV')
        print('Total Energy = ', round_sig(self.total_energy/1.E9), 'GeV')
        print('Gamma = ', round_sig(self.gamma))
        print('Beta = ', round_sig(self.beta))
        print('Momentum = ', round_sig(self.momentum/1E9, 8), 'GeV/c')       
        print('Rigidity = ', round_sig(self.rigidity), 'Tm')


########################################################################
# Check if file exists
########################################################################
def check_if_file_exists(name):
    ret_val = False
    if os.path.isfile(name):
        print (name, ' exists')
        ret_val = True
    return ret_val
    
########################################################################
# Move file - check and overwrite
########################################################################
def move_file(file, to):    
    # Check initial file exists
    if os.path.isfile(file):
        print('move_file:: file ', file, ' exists.')
        
        # Check if final destination file already exists
        end_file = to + '/' + file
        if os.path.isfile(end_file):
            print('move_file:: file ', file, ' already exists in directory, overwriting.')
            os.replace(file, end_file)
        else:
            print('move_file:: file ', file, ' not in ',to,', moving.')
            shutil.move(file, to)           
    else: 
        print('move_file:: file ', file, ' does not exist. Doing nothing.')
    
########################################################################
# Make directory
########################################################################  
def make_directory(path, overwrite=False):
    if os.path.isdir(path):
        print ("Directory %s exists" % path)  
        
        if overwrite:
            os.rmdir(path)
            print ("Directory %s removed" % path)  
            try:
                os.mkdir(path)
            except OSError:
                print ("Creation of the directory %s failed" % path)
            else:
                print ("Successfully created the directory %s" % path)  
    else:
        try:
            os.mkdir(path)
        except OSError:
            print ("Creation of the directory %s failed" % path)
        else:
            print ("Successfully created the directory %s" % path)  
########################################################################
# Relativistic functions
########################################################################             
        
def LorentzGamma(E_tot, E_rest=938.27208816E6):
    return (E_tot / E_rest)
    
def LorentzGamma_from_beta(beta):
    return (1./np.sqrt(1.-beta**2))    

def LorentzBeta(gamma):
    return np.sqrt( 1. - (1./gamma**2) )

def RelativisticMomentum(gamma, E_rest=938.27208816E6):
    return (gamma * E_rest * LorentzBeta(gamma))

    
def E_from_gamma(gamma, E_rest=938.27208816E6):
    return (gamma*E_rest)

########################################################################
# Delta P over P from dE or vice versa
# dp_ov_p = dE_ov_E/beta^2
########################################################################
def dpp_from_dE(dE, E, beta):
    return (dE / (E * beta**2))
    
def dE_from_dpp(dpp, E, beta):
    return (dpp * E * beta**2)

def z_to_time(z, beta): 
    c = 299792458
    return z / (c * beta)
    
########################################################################
# Round number to n significant figures
########################################################################    
def round_sig(x, sig=3):
    if x == 0.0 : return 0.0
    else: return round(x, sig-int(floor(log10(abs(x))))-1)      

########################################################################
# lowercase all strings in a dataframe
########################################################################
def pandas_lowercase_all_strings(dataframe):    
    dataframe.columns = map(str.lower, dataframe.columns)

    for key in dataframe.keys():
        if dataframe[str(key)].dtype == 'O':
            #print(str(key+' column now lowercase'))
            dataframe[str(key)] = dataframe[str(key)].str.lower()
            
    return dataframe

########################################################################
# lowercase column names for a given dataframe
########################################################################
def pandas_save_to_file(dataframe, filename):
    with open(filename, 'w') as f:
        dfAsString = dataframe.to_string()
        f.write(dfAsString)
        
########################################################################
# lowercase column names for a given dataframe
########################################################################
def pandas_dataframe_lowercase_columns(dataframe):
    dataframe.columns = map(str.lower, dataframe.columns)

########################################################################
# Return ISIS synchrotron beam momentum at given time in ms - Billy Kyle
########################################################################    
def synchrotron_momentum(max_E, time):
    mpeV = m_p * c**2 / e           # Proton mass in eV
    R0 = 26                         # Mean machine radius
    n_dip = 10
    dip_l = 4.4                     # Dipole length
    
    dip_angle = 2 * np.pi / n_dip
    rho = dip_l / dip_angle      
    omega = 2 * np.pi * 50   
    
    Ek = np.array([70, max_E]) * 1e6
    E = Ek + mpeV
    p = np.sqrt(E**2 - mpeV**2)

    B = p / c / rho
    
    Bdip = lambda t: (B[1] + B[0] - (B[1] - B[0]) * np.cos(omega * t)) / 2  # Idealised B-field variation
    pdip = lambda t: Bdip(t) * rho * c                                      # Momentum from B-field in MeV
    
    return pdip(time*1E-3)

########################################################################
# Return normalised momentum pn
########################################################################  
def return_pn(time_array, pr=3.96):
    out_array = []
    # time in ms
    for t in time_array:
        out_array.append(0.5*(pr+1-(pr-1)*np.cos(0.1*np.pi*t)))
    return out_array

########################################################################
# Return ISIS synchrotron beam energy at given time in ms
########################################################################    
# Incorrect
# ~ def synchrotron_energy(max_E, time):
    # ~ b = (max_E-70)/2
    # ~ sol = (b*(math.sin(time*(np.pi/10)-(np.pi/2))))+b+70
    # ~ return sol
    
def synchrotron_energy(max_E, time):
    mpeV = m_p * c**2 / e           # Proton mass in eV    
    return (np.sqrt(synchrotron_momentum(max_E, time)**2 + mpeV**2) - mpeV)/1E6
    
def synchrotron_kinetic_energy(max_E, time):
    mpeV = m_p * c**2 / e           # Proton mass in eV    
    # Relativistic Kinetic Energy = Relativistic Energy - mass
    return (np.sqrt(synchrotron_momentum(max_E, time)**2 + mpeV**2) - mpeV) # Return array in eV
    #return (np.sqrt(synchrotron_momentum(max_E, time)**2 + mpeV**2) - mpeV)/1E6 # Return array in MeV
    
########################################################################
# Return ISIS synchrotron beam energy for an array of time in ms
########################################################################    
def synchrotron_energy_array(max_E, time_array):
    out_array = []
    for t in time_array:
        out_array.append(synchrotron_energy(max_E, t))
        
    return out_array
    
########################################################################
# Return ISIS synchrotron beam energy dataframe for a single time
# Includes additional information from MADX_Proton_Beam_Parameters class
########################################################################    
def synchrotron_energy_data(max_E, time):
    
    energy = synchrotron_energy(max_E, time)
    beam = MADX_Proton_Beam_Parameters(energy*1E6)

    gamma = []
    beta = []
    momentum = []
    rigidity = []
    e_gev = []    

    e_gev.append((beam.get_total_energy()-938.272E6)/1.E9)
    gamma.append(beam.get_gamma())
    beta.append(beam.get_beta())
    momentum.append(beam.get_momentum()/1E9)
    rigidity.append(beam.get_rigidity())

    df = pd.DataFrame({'Time [ms]':time, 'Energy [MeV]':energy, 'Energy [GeV]':e_gev, 'Momentum [GeV/c]': momentum, 'Gamma': gamma, 'Beta': beta, 'Rigidity [Tm]':rigidity})
        
    return df

########################################################################
# Return ISIS synchrotron beam energy dataframe for full cycle
# Includes additional information from MADX_Proton_Beam_Parameters class
########################################################################    
def synchrotron_energy_df(max_E, intervals=20):
    
    time_array = np.linspace(0., 10, intervals+1)
    energies = synchrotron_energy_array(max_E, time_array)   
    
    gamma = []
    beta = []
    momentum = []
    rigidity = []
    e_gev = []
    
    for E in energies:
        beam = MADX_Proton_Beam_Parameters(E*1E6) # convert MeV to eV
    
        e_gev.append((beam.get_total_energy()-938.272E6)/1.E9)
        gamma.append(beam.get_gamma())
        beta.append(beam.get_beta())
        momentum.append(beam.get_momentum()/1E9)
        rigidity.append(beam.get_rigidity())
    
    df = pd.DataFrame({'Time [ms]':time_array, 'Energy [MeV]':energies, 'Energy [GeV]':e_gev, 'Momentum [GeV/c]': momentum, 'Gamma': gamma, 'Beta': beta, 'Rigidity [Tm]':rigidity})
        
    return df

def synchrotron_kinetic_energy_df(max_E, time_array):
    
    energies = synchrotron_kinetic_energy(max_E, time_array)
    
    gamma = []
    beta = []
    momentum = []
    rigidity = []
    e_mev = []
    
    mpeV = m_p * c**2 / e           # Proton mass in eV
    
    for E in energies:
        beam = MADX_Proton_Beam_Parameters(E) 
    
        e_mev.append((beam.get_total_energy()-mpeV)/1.E6)
        gamma.append(beam.get_gamma())
        beta.append(beam.get_beta())
        momentum.append(beam.get_momentum()/1E9)
        rigidity.append(beam.get_rigidity())    
    
    df = pd.DataFrame({'Time [ms]':time_array, 'Energy [eV]':energies, 'Energy [MeV]':e_mev, 'Momentum [GeV/c]': momentum, 'Gamma': gamma, 'Beta': beta, 'Rigidity [Tm]':rigidity})    
    return df
########################################################################
# Return steering kick in mrad given the programmed kick in amperes, the
# measurement time, max energy, plane and super-period
########################################################################    
def calculate_steering_kick(amps, max_E, time, plane ='H', sp=0):
    
    sp_list = [0, 2, 3, 4, 5, 7, 9]
    if sp not in sp_list:
        print('calculate_steering_kick:: selected super-period has no steering magnet')
        exit(0)
    
    # Calibration provided by HVC 30.09.22
    calibration_data = {
        '0H' : 0.08350,
        '2H' : 0.09121,
        '3H' : 0.08,
        '4H' : 0.06600,
        '5H' : 0.07780,
        '7H' : 0.07580,
        '9H' : 0.07660,
        '0V' : 0.04620,
        '2V' : 0.04330,
        '3V' : 0.05210,
        '4V' : 0.04770,
        '5V' : 0.05400,
        '7V' : 0.05220,
        '9V' : 0.04510,    
    }
    
    df = synchrotron_energy_data(max_E, time)
        
    h_list = ['h', 'H', 'horizontal', 'Horizontal']
    if plane in h_list: key = str(sp) + 'H'   
    else: key = str(sp) + 'V'        
        
    return round_sig(float(amps*(calibration_data[key]/df['Rigidity [Tm]']))) # kick in milliradians
    
########################################################################
# Return steering current in amperes given the calculated kick in mrad,
# the measurement time, max energy, plane and super-period
########################################################################    
def calculate_steering_current(kick_mrad, max_E, time, plane ='H', sp=0):
    
    sp_list = [0, 2, 3, 4, 5, 7, 9]
    if sp not in sp_list:
        print('calculate_steering_kick:: selected super-period has no steering magnet')
        exit(0)
    
    # Calibration provided by HVC 30.09.22
    calibration_data = {
        '0H' : 0.08350,
        '2H' : 0.09121,
        '3H' : 0.08,
        '4H' : 0.06600,
        '5H' : 0.07780,
        '7H' : 0.07580,
        '9H' : 0.07660,
        '0V' : 0.04620,
        '2V' : 0.04330,
        '3V' : 0.05210,
        '4V' : 0.04770,
        '5V' : 0.05400,
        '7V' : 0.05220,
        '9V' : 0.04510,    
    }
    
    df = synchrotron_energy_data(max_E, time)
        
    h_list = ['h', 'H', 'horizontal', 'Horizontal']
    if plane in h_list: key = str(sp) + 'H'   
    else: key = str(sp) + 'V'        
        
    return round_sig(float((kick_mrad*1000)*(df['Rigidity [Tm]']/calibration_data[key]))) # current in amps

########################################################################
# Return steering kick in mrad given the programmed kick in amperes for 
# all times in the cycle (201 data points)
########################################################################  
def calculate_steering_kick_all_times(amps, max_E, plane='H', sp=0):
    
    time_array = np.linspace(0., 10, 201)
    kicks = []
    
    for t in time_array:
        kicks.append(calculate_steering_kick(amps, max_E, t, plane, sp))
    
    return (time_array, kicks)
    
########################################################################
# Return steering kick in mrad given the programmed kick in amperes for 
# all times in the cycle (201 data points)
########################################################################  
def is_idendity_matrix_6(m):
    
    identity = np.matrix([[1., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 0.], [0., 0., 1., 0., 0., 0.], [0., 0., 0., 1., 0., 0.], [0., 0., 0., 0., 1., 0.], [0., 0., 0., 0., 0., 1.]])
    comparison_matrix = (m == identity)
    
    #return comparison_matrix
    
    if comparison_matrix.all() == True:
        return True
    else:
        return False
        
########################################################################
# Decompose 6x6 matrix into horizontal and vertical 2x2 matrices
########################################################################  
def transport_matrix_transverse_decomposition(m6):
    
    mH = m6[np.ix_([0,1],[0,1])]
    mV = m6[np.ix_([2,3],[2,3])]
    
    return mH, mV        

########################################################################
# phase advance from 2x2 transport matrix
########################################################################  
def transverse_matrix_phase_advance(m2): 
    if m2.trace() >= 2.0:
        raise ValueError('transverse_matrix_phase_advance: Trace of matrix exceeds 2.0 at',m2.trace())
    try: phi = round_sig(np.arccos((m2[0,0] + m2[1,1])/2)/(2*np.pi),6)
    except ValueError:
        phi=0
    return phi
    
########################################################################
# Lattice functions (?) from one turn map - untested
########################################################################  
def transverse_matrix_lattice_functions(m2):
    output_params = {}
    phi = transverse_matrix_phase_advance(m2)
    output_params['beta'] = m2[0,1]/np.sin(phi)
    output_params['alpha'] = (m2[0,0]-m2[1,1,])/2*np.sin(phi)
    output_params['gamma'] = -m2[1,0]/np.sin(phi)
    output_params['gamma2'] = (1+output_params['alpha']**2)/output_params['beta']
    
    return output_params    

########################################################################
# multiply all matrices in a matrix dict
########################################################################  
def multiply_matrices(matrix_dict):
    first_element = True
    sum_matrix = np.empty(shape=(matrix_dict[list(matrix_dict.keys())[0]].shape))
    for m in matrix_dict.values():   
        if first_element:
            sum_matrix = m
            first_element = False
        else: sum_matrix = np.matmul(m, sum_matrix)
        
    return sum_matrix

########################################################################
# resonance_lines class from CERN
########################################################################  
class resonance_lines(object):

    def __init__(self, Qx_range, Qy_range, orders, periodicity, legend=False):

        if np.std(Qx_range):
            self.Qx_min = np.min(Qx_range)
            self.Qx_max = np.max(Qx_range)
        else:
            self.Qx_min = np.floor(Qx_range)-0.05
            self.Qx_max = np.floor(Qx_range)+1.05
        if np.std(Qy_range):
            self.Qy_min = np.min(Qy_range)
            self.Qy_max = np.max(Qy_range)
        else:
            self.Qy_min = np.floor(Qy_range)-0.05
            self.Qy_max = np.floor(Qy_range)+1.05

        self.periodicity = periodicity
        self.legend_flag = legend

        nx, ny = [], []

        for order in np.nditer(np.array(orders)):
            t = np.array(range(-order, order+1))
            nx.extend(order - np.abs(t))
            ny.extend(t)
        nx = np.array(nx)
        ny = np.array(ny)

        cextr = np.array([nx*np.floor(self.Qx_min)+ny*np.floor(self.Qy_min), \
                          nx*np.ceil(self.Qx_max)+ny*np.floor(self.Qy_min), \
                          nx*np.floor(self.Qx_min)+ny*np.ceil(self.Qy_max), \
                          nx*np.ceil(self.Qx_max)+ny*np.ceil(self.Qy_max)], dtype='int')
        cmin = np.min(cextr, axis=0)
        cmax = np.max(cextr, axis=0)
        res_sum = [range(cmin[i], cmax[i]+1) for i in range(cextr.shape[1])]
        self.resonance_list = zip(nx, ny, res_sum)

    def plot_resonance_ax(self, ax1):
        # Remove Borders
        #ax1.spines['top'].set_visible(False);
        #ax1.spines['bottom'].set_visible(False);
        #ax1.spines['left'].set_visible(False);
        #ax1.spines['right'].set_visible(False);
        #ax1.axes.get_yaxis().set_visible(False); 
        #ax1.axes.get_xaxis().set_visible(False);  
        
        Qx_min = self.Qx_min
        Qx_max = self.Qx_max
        Qy_min = self.Qy_min
        Qy_max = self.Qy_max 
        ax1.set_xlim(Qx_min, Qx_max);
        ax1.set_ylim(Qy_min, Qy_max);        
        ax1.set_xlabel(r'Horizontal Tune $Q_x$')
        ax1.set_ylabel(r'Vertical Tune $Q_y$')

        for resonance in self.resonance_list:
            nx = resonance[0]
            ny = resonance[1]
            for res_sum in resonance[2]:
                if ny:                    
                    
                    if ny%2:                                                  
                        if res_sum%self.periodicity:
                            ax1.plot([Qx_min, Qx_max], [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny], ls='--', color='b', lw=0.5, label='Non-Systematic Skew')
                        else:
                            ax1.plot([Qx_min, Qx_max], [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny], ls='--', color='r', lw=1, label='Systematic Skew')               
                        
                    else:
                        if res_sum%self.periodicity:
                            ax1.plot([Qx_min, Qx_max], [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny], color='b', lw=0.5, label='Non-Systematic Normal')
                        else:
                            ax1.plot([Qx_min, Qx_max], [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny], color='r', lw=1, label='Systematic Normal')
                else:                
                    if res_sum%self.periodicity:
                        ax1.plot([np.float(res_sum)/nx, np.float(res_sum)/nx],[Qy_min, Qy_max], color='b', lw=0.5, label='Non-Systematic Normal')
                    else:
                        ax1.plot([np.float(res_sum)/nx, np.float(res_sum)/nx],[Qy_min, Qy_max], color='r', lw=1, label='Systematic Normal')
    
        if self.legend_flag: 
            custom_lines = [Line2D([0], [0], color='r', lw=4),
                Line2D([0], [0], color='b', lw=4),
                Line2D([0], [2], color='r', lw=4, ls='--'),
                Line2D([0], [2], color='b', lw=4, ls='--')]
            ax1.legend(custom_lines, ['Systematic Normal', 'Non-Systematic Normal', 'Systematic Skew', 'Non-Systematic Skew'])
    
    
    def plot_resonance_fig(self, figure_object = None):
        plt.ion()
        if figure_object:
            fig = figure_object
            plt.figure(fig.number)
        else:
            fig = plt.figure()
        Qx_min = self.Qx_min
        Qx_max = self.Qx_max
        Qy_min = self.Qy_min
        Qy_max = self.Qy_max 
        plt.xlim(Qx_min, Qx_max)
        plt.ylim(Qy_min, Qy_max)
        plt.xlabel(r'Horizontal Tune $Q_x$')
        plt.ylabel(r'Vertical Tune $Q_y$')
        for resonance in self.resonance_list:
            nx = resonance[0]
            ny = resonance[1]
            for res_sum in resonance[2]:
                if ny:
                    line, = plt.plot([Qx_min, Qx_max], \
                        [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny])
                else:
                    line, = plt.plot([np.float(res_sum)/nx, np.float(res_sum)/nx],[Qy_min, Qy_max])
                    
                    
                    
                if ny%2:
                    plt.setp(line, linestyle='--') # for skew resonances
                if res_sum%self.periodicity:
                    plt.setp(line, color='b') # non-systematic resonances
                else:
                    plt.setp(line, color='r', linewidth=2.0) # systematic resonances
        plt.draw()
        return fig
    
    def print_resonances(self):
        for resonance in self.resonance_list:
            for res_sum in resonance[2]:
                print_string =  '%s %s%s = %s\t%s'%(str(resonance[0]).rjust(2), ("+", "-")[resonance[1]<0], \
                        str(abs(resonance[1])).rjust(2), str(res_sum).rjust(4), \
                        ("(non-systematic)", "(systematic)")[res_sum%self.periodicity==0])
                print(print_string)
                
########################################################################
# Read PTC Twiss and return dictionary of columns/values
# Obsolete due to tfs-pandas library
########################################################################
def Read_PTC_Twiss_Return_Dict(filename, verbose=True):
    # Dictionary for output
    d = dict()
    d['HEADER_FILENAME'] = filename
    keywords = ''
    
    # First we open and count header lines
    fin0=open(filename,'r').readlines()
    headerlines = 0
    for l in fin0:
        # Store each header line
        headerlines = headerlines + 1
        # Stop if we find the line starting '* NAME'
        if '* NAME' in l:
            keywords = l
            break
        # Store the headers as d['HEADER_<name>'] = <value>
        else:
            if '"' in l:
                d[str('HEADER_'+l.split()[1])]=[str(l.split('"')[1])]
            else:
                d[str('HEADER_'+l.split()[1])]=[float(l.split()[-1])]                 
    headerlines = headerlines + 1    
    
    if verbose: print ('\nRead_PTC_Twiss_Return_Dict found Keywords: \n',keywords)
    
    # Make a list of column keywords to return (as an aid to iterating)
    dict_keys = []
    for key in keywords.split():
        dict_keys.append(key)
    dict_keys.remove('*')
    
    if verbose: print ('\nRead_PTC_Twiss_Return_Dict Dict Keys: \n',dict_keys)
    
    # Initialise empty dictionary entries for column keywords 
    for key in dict_keys:
        d[key]=[]
        
    if verbose: print ('\nRead_PTC_Twiss_Return_Dict header only dictionary \n', d)
    
    # Strip header
    fin1=open(filename,'r').readlines()[headerlines:]   
    
    # Populate the dictionary line by line
    for l in fin1:
        i = -1        
        for value in l.split():
            i = i+1
            if 'KEYWORD' or 'NAME' in dict_keys[i]:
                d[dict_keys[i]].append(str(value))
            else:
                d[dict_keys[i]].append(float(value))    
                
    # Return list of column keywords 'dict_keys', and dictionary 'd'
    return dict_keys, d
   
# Adapted from Alex Huschauer's (CERN ABP) webtools found at:
# https://gitlab.cern.ch/acc-models/acc-models-ps/-/tree/2021/_scripts/web    
def plot_lattice_elements(ax, twiss_in, suppress_x=True):

    # Set plot limits
    ax.set_ylim(-1.5,1.5)
    
    # Suppress ticks on y-axis
    #ax.set_yticks([])
    #ax.set_yticklabels([])
    ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        
    if suppress_x:
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        #ax.set_xticks([])
        #ax.set_xticklabels([])

    # Extract Twiss Header
    twiss = tfs.read(twiss_in)
    twissHeader = dict(twiss.headers)
    print('plot_lattice_elements for sequence ', twissHeader['SEQUENCE'])
    
    # Positions and lengths of elements
    pos = twiss.S.values - twiss.L.values/2
    lengths = twiss.L.values
    total_length = (pos[-1]+lengths[-1])
    print('Full length of accelerator lattice = ', total_length, 'm')
    
    # modify lengths in order to plot zero-length elements
    lengths[np.where(lengths == 0)[0]] += 0.001
    
    # Plot line through centre
    ax.plot([0, total_length], [0., 0.], color='grey', linestyle='-', linewidth=0.5)
    
    # Markers - black 0.1m centred lines    
    idx = np.array([idx for idx, elem in enumerate(twiss.KEYWORD.values) if 'MARKER' in elem])
    for i in idx:
        ax.add_patch(Rectangle((pos[i], -0.5), width = 0.1, height = 1., angle=0.0, ec='k', fc='k', lw=0.0))
    
    # BENDS - blue centred rectangles   
    idx = np.array([idx for idx, elem in enumerate(twiss.KEYWORD.values) if 'BEND' in elem])
    for i in idx:
        ax.add_patch(Rectangle((pos[i]-1, -0.5), width = lengths[i], height = 1., angle=0.0, ec='k', fc='b', lw=0.0))
    
    # Kickers - cyan centred rectangles 
    idx = np.array([idx for idx, elem in enumerate(twiss.KEYWORD.values) if 'HKICKER' in elem])
    for i in idx:
        ax.add_patch(Rectangle((pos[i], -0.5), width = 0.1, height = 1., angle=0.0, ec='k', fc='c', lw=0.0))
    idx = np.array([idx for idx, elem in enumerate(twiss.KEYWORD.values) if 'VKICKER' in elem])       
    for i in idx:
        ax.add_patch(Rectangle((pos[i], -0.5), width = 0.1, height = 1., angle=0.0, ec='k', fc='c', lw=0.0))
              
    # QUADRUPOLES - red offset rectangles indicating Focussing or Defocussing
    idx = np.array([idx for idx, elem in enumerate(twiss.KEYWORD.values) if 'QUADRUPOLE' in elem])
    name = np.array(twiss.NAME.values)[idx]
    if (twissHeader['SEQUENCE'] == 'SYNCHROTRON'):        
        idx_1 = idx[np.array([i for i, n in enumerate(name) if 'QD' in n])]
        idx_2 = idx[np.array([i for i, n in enumerate(name) if 'QF' in n])]
        
    elif (twissHeader['SEQUENCE'] == 'RING'):    
        idx_1 = idx[np.array([i for i, n in enumerate(name) if (n.startswith('D') or 'QTD' in n)])]
        idx_2 = idx[np.array([i for i, n in enumerate(name) if (n.startswith('F') or 'QTF' in n)])]
        
    else: 
        idx_1 = idx[np.array([i for i, n in enumerate(name) if 'D' in n ])]
        idx_2 = idx[np.array([i for i, n in enumerate(name) if 'F' in n ])]        
        
    offset = [-0.5, 0.5]
    for i in idx_1:
        ax.add_patch(Rectangle((pos[i], (-0.5 + offset[0])), width = lengths[i], height = 1., angle=0.0, ec='k', fc='r', lw=0.0))
    for i in idx_2:
        ax.add_patch(Rectangle((pos[i], (-0.5 + offset[1])), width = lengths[i], height = 1., angle=0.0, ec='k', fc='r', lw=0.0))
   
