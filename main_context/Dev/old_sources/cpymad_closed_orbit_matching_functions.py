#!/usr/bin/env python
# cpymad functions to be used with modified ISIS lattice for closed orbit calculations
# Created: 17.06.22 Haroon Rafique STFC ISIS Accelerator Division
# Requires cpymad tfs-pandas python libraries

from helper_functions import *
from cpymad_helpers import *

########################################################################
########################################################################
#                         BPM Data Files
########################################################################
########################################################################


########################################################################
# Read BPM dat file and return pandas daataframe of raw data
########################################################################
def bpm_dat_to_df(filename):
    return pnd.read_table(filename, names=['Monitor', 0,1,2,3,4,5,6,7,8,9], delim_whitespace=True)

########################################################################
# Read BPM acquisition file + twiss file return prrocessed dataframe
########################################################################
def process_bpm_df(df_bpm, twiss_file, verbose=False):    
    
    # empty lists for new dataframe column data
    s_array = []
    monitor_array = []
    phase_array = []
    beta_array = []
    
    # Read Twiss file, lowercase all values and columns
    df_out = df_bpm
    df_myTwiss_in = tfs.read(twiss_file)
    df_myTwiss = pandas_lowercase_all_strings(df_myTwiss_in)
    
    # Iterate over BPM dataframe
    for index, row in df_out.iterrows(): 
        
        # Check if BPM name in Twiss file
        if df_myTwiss['name'].str.contains(str(row['Monitor']).lower()).any():

            monitor_str = str(row['Monitor']).lower()
            if verbose: print(row['Monitor']+ ' in lattice as '+ df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
            
            # Store monitor name and S from Twiss file
            monitor_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
            s_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].s[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
            
            if 'hm' in monitor_str:
                phase_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].mux[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
                beta_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].betx[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
            else:
                phase_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].muy[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
                beta_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].bety[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
        else:
            if verbose: print(row['Monitor']+ ' not in lattice')
            
    # Remove and replace old monitor column
    df_out = df_out.drop('Monitor', axis=1)
    
    # Add S column and reshuffle order of columns
    #if len(s_array) != len(df_out['S']):
    #    print(len(s_array), len(monitor_array), len(df_out['S']))
    #    print (df_out)
    #    print (s_array)
        
    df_out['S'] = s_array
    df_out['Monitor'] = monitor_array
    
    df_out['Mu'] = phase_array
    df_out['Beta'] = beta_array   
        
    mu_col = df_out.pop('Mu')
    df_out.insert(0, 'Mu', mu_col)
    bet_col = df_out.pop('Beta')
    df_out.insert(0, 'Beta', bet_col)        
    mon_col = df_out.pop('Monitor')
    df_out.insert(0, 'Monitor', mon_col)
    s_col = df_out.pop('S')
    df_out.insert(0, 'S', s_col)
    
    # ~ cols = df_out.columns.tolist()
    # ~ cols = cols[-1:] + cols[:-1]
    # ~ cols = cols[-1:] + cols[:-1]
    # ~ cols = cols[-1:] + cols[:-1]
    # ~ cols = cols[-1:] + cols[:-1]
    #df_out = df_out[cols]
    
    return df_out
    
def process_bpm_dat_df(dat_file, twiss_file, verbose=False):    
    df_bpm = bpm_dat_to_df(dat_file)    
    return process_bpm_df(df_bpm, twiss_file, verbose=False)

########################################################################
# Read BPM acquisition file + twiss file return prrocessed dataframe
########################################################################
def process_bpm_dat_df2(dat_file, twiss_file, verbose=False):
    
    df_bpm = bpm_dat_to_df(dat_file)
    horizontal = True
    
    # empty lists for new dataframe column data
    s_array = []
    monitor_array = []
    phase_array = []
    beta_array = []
    
    # Read Twiss file, lowercase all values and columns
    df_out = df_bpm
    df_myTwiss_in = tfs.read(twiss_file)
    df_myTwiss = pandas_lowercase_all_strings(df_myTwiss_in)
    
    # Iterate over BPM dataframe
    for index, row in df_out.iterrows(): 
        
        # Check if BPM name in Twiss file
        if df_myTwiss['name'].str.contains(str(row['Monitor']).lower()).any():

            monitor_str = str(row['Monitor']).lower()
            if verbose: print(row['Monitor']+ ' in lattice as '+ df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
            
            # Store monitor name and S from Twiss file
            monitor_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
            s_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].s[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
            
            if 'hm' in row['Monitor']:
                phase_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].mux[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
                beta_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].betx[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
            else:
                phase_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].muy[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
                beta_array.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].bety[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
        else:
            if verbose: print(row['Monitor']+ ' not in lattice')
            
    # Remove and replace old monitor column
    df_out = df_out.drop('Monitor', axis=1)
    
    # Add S column and reshuffle order of columns
    df_out['S'] = s_array
    df_out['Monitor'] = monitor_array
    df_out['Mu'] = phase_array
    df_out['Beta'] = beta_array   
        
    cols = df_out.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    cols = cols[-1:] + cols[:-1]
    cols = cols[-1:] + cols[:-1]
    cols = cols[-1:] + cols[:-1]
    df_out = df_out[cols]
    
    return df_out

########################################################################
# Read BPM acquisition file + twiss file return prrocessed dataframe H/V
########################################################################
def process_bpm_dat_df_h(dat_file, twiss_file, verbose=False):
    return process_bpm_dat_df_hv(dat_file, twiss_file, horizontal=True, verbose=verbose)

def process_bpm_dat_df_v(dat_file, twiss_file, verbose=False):
    return process_bpm_dat_df_hv(dat_file, twiss_file, horizontal=False, verbose=verbose)

def process_bpm_dat_df_hv(dat_file, twiss_file, horizontal=True, verbose=False):
    
    df_bpm = bpm_dat_to_df(dat_file)
    
    # empty lists for new dataframe column data
    s_list=[]
    monitor_list=[]
    
    # Read Twiss file, lowercase all values and columns
    df_out = df_bpm
    df_myTwiss_in = tfs.read(twiss_file)
    df_myTwiss = pandas_lowercase_all_strings(df_myTwiss_in)
    
    # Iterate over BPM dataframe
    for index, row in df_out.iterrows(): 
        
        # Check if BPM name in Twiss file
        if df_myTwiss['name'].str.contains(str(row['Monitor']).lower()).any():

            monitor_str = str(row['Monitor']).lower()
            if verbose: print(row['Monitor']+ ' in lattice as '+ df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
            
            # Store monitor name and S from Twiss file
            monitor_list.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
            s_list.append(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].s[int(df_myTwiss[df_myTwiss['name'].str.contains(monitor_str)].index[0])])
        else:
            if verbose: print(row['Monitor']+ ' not in lattice')
    
    # Remove and replace old monitor column
    df_out = df_out.drop('Monitor', axis=1)
    df_out['Monitor'] = monitor_list
    
    # Add S column and reshuffle order of columns
    df_out['S'] = s_list
    cols = df_out.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    cols = cols[-1:] + cols[:-1]
    df_out = df_out[cols]
    
    # Remove unwanted H or V data
    if horizontal: df_out.drop(df_out[df_out['Monitor'].str.contains('vm')].index, inplace=True)
    else: df_out.drop(df_out[df_out['Monitor'].str.contains('hm')].index, inplace=True)
    
    df_out.reset_index(inplace=True, drop=True)
    
    return df_out
    
########################################################################
# Read processed BPM dat file and return pandas daataframe of raw data
########################################################################
def processed_bpm_dat_to_df(filename):
    return pnd.read_table(filename, header=0, delim_whitespace=True)
    
########################################################################
# Difference between 2 processed bpm data files
########################################################################
def generate_processed_bpm_difference(processed_dat_file_1, processed_dat_file_2, output_filename, time='0'):
    
    df1 = processed_bpm_dat_to_df(processed_dat_file_1)
    df2 = processed_bpm_dat_to_df(processed_dat_file_2)
    
    return generate_processed_bpm_difference_df(df1, df2, output_filename, time)
    
########################################################################
# Difference between 2 processed bpm dataframes
########################################################################    
def generate_processed_bpm_difference_df(df1, df2, output_filename, time='0'):
    mon_list = []
    x_list = []
    s_list = []
    mu_list = []
    beta_list = []

    for index, row in df1.iterrows():   
        if df2['Monitor'].str.contains(str(row['Monitor'])).any():
            monitor_str = str(row['Monitor']).lower()
            mon_list.append(monitor_str)
            # Use try/excepts to handle cases where processed data file
            # stores bpm data under time '0' by default            
            try: diff = float(row[time] - df2[df2['Monitor'].str.contains(monitor_str)][time])
            except KeyError:
                try: diff = float(row[str(time)] - df2[df2['Monitor'].str.contains(monitor_str)][str(time)])
                except KeyError:
                    try: diff = float(row[float(time)] - df2[df2['Monitor'].str.contains(monitor_str)][float(time)])
                    except KeyError:
                        diff = float(row['0'] - df2[df2['Monitor'].str.contains(monitor_str)]['0'])
                    
            if not np.isnan(diff):
                x_list.append(diff)
                s_list.append(row['S'])
                mu_list.append(row['Mu'])
                beta_list.append(row['Beta'])

    out_df = pnd.DataFrame(list(zip(s_list, mon_list, mu_list, beta_list, x_list)), columns =['S', 'Monitor', 'Mu', 'Beta', time])
    pandas_save_to_file(out_df, output_filename)
    return out_df

########################################################################
# Sum of 2 processed bpm data files
########################################################################
def generate_processed_bpm_sum(processed_dat_file_1, processed_dat_file_2, output_filename, time1='0', time2='0'):
    
    df1 = processed_bpm_dat_to_df(processed_dat_file_1)
    df2 = processed_bpm_dat_to_df(processed_dat_file_2)
    
    mon_list = []
    x_list = []
    s_list = []
    mu_list = []
    beta_list = []

    for index, row in df1.iterrows():   
        if df2['Monitor'].str.contains(str(row['Monitor'])).any():
            monitor_str = str(row['Monitor']).lower()
            mon_list.append(monitor_str)
            # Use try/excepts to handle cases where processed data file
            # stores bpm data under time '0' by default
            summ = float(row[time1] + df2[df2['Monitor'].str.contains(monitor_str)][time2])
                    
            if not np.isnan(summ):
                x_list.append(summ)
                s_list.append(row['S'])
                mu_list.append(row['Mu'])
                beta_list.append(row['Beta'])

    out_df = pnd.DataFrame(list(zip(s_list, mon_list, mu_list, beta_list, x_list)), columns =['S', 'Monitor', 'Mu', 'Beta', time])
    pandas_save_to_file(out_df, output_filename)
    return out_df
########################################################################
# Generate n dummy kicks using kicker knob variables
########################################################################
def set_n_random_kicks_from_knobs(madx_instance, cpymad_logfile, num_kicks, vary_list, min_kick=-1E-3, max_kick=1E-3): 
    log_string = 'set_n_random_kicks_from_knobs: '
    random_kicks = np.random.uniform(size=num_kicks, low=min_kick, high=max_kick)
    knobs = []
    
    # Select n random knobs from vary list
    random_vary_list = random.sample(vary_list, num_kicks)
    
    # iterate over knobs, set to random value, add to log string
    i = 0
    for knob in random_vary_list:        
        exec_str = 'madx_instance.globals.'+str(knob)+' = '+str(random_kicks[i]) # note kick in mrad
        exec(exec_str) 
        knobs.append(knob)
        log_string += (exec_str+'\n')
        i += 1
        
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    zipped = list(zip(knobs, random_kicks*1E3))    
    df = pnd.DataFrame(zipped, columns=['knob', 'kick (mrad)'])

    return df

########################################################################
########################################################################
#                         Twiss Files
########################################################################
########################################################################


########################################################################
# Compare and plot twiss closed orbit with measurements from file
########################################################################
def compare_measurement_and_twiss(twiss_file, measurement_file, savename, horizontal=True, threshold=2.0, time='0'):
    
    # Read measurement file 
    df_meas = processed_bpm_dat_to_df(measurement_file)
    #s_list, meas_list = measurements_from_file(measurement_file, horizontal=horizontal)
    
    # Read Twiss file return dataframe #'Monitor', 'S', 'mu_x', 'x'
    df = read_twiss_generate_measured_co_df(twiss_file, horizontal, False, time)
    
    # Read Twiss file to check orbit
    df_myTwiss = tfs.read(twiss_file)
    df_myTwiss = pandas_lowercase_all_strings(df_myTwiss)
    if horizontal:
        if (df_myTwiss['x'] == 0).all() and df_myTwiss['x'].std()==0:
            print('WARNING: compare_measurement_and_twiss: twiss file '+twiss_file+' has perfect closed orbit in this plane')
            return 10
    else:
        if (df_myTwiss['y'] == 0).all() and df_myTwiss['y'].std()==0:
            print('WARNING: compare_measurement_and_twiss: twiss file '+twiss_file+' has perfect closed orbit in this plane')
            return 10  
        
    # Iterate over measurements, save difference, compute std
    diff_list = []
    for i in range(0,len(df_meas),1):
        diff = df_meas[time][i] - float(df[df['S']==df_meas['S'][i]][time])
        diff_list.append(diff)
    std_diff = np.std(diff_list)
    
    if std_diff >= threshold:
        return std_diff
    
    plot_str = r'$\sigma_{difference}$ = '+str(round_sig(std_diff,3))    
    title_str = 'Standard Deviation\n' + plot_str + ' mm'
    
    fig1 = plt.figure(facecolor='w', edgecolor='k', figsize=[9, 6])
    gs = fig1.add_gridspec(ncols=2, nrows=1, width_ratios=[1,1])
    gs.update(wspace=0.25, hspace=0.)

    ax1 = fig1.add_subplot(gs[0,0])
    ax2 = fig1.add_subplot(gs[0,1])    
    
    # plot closed orbit
    if horizontal: 
        ax1.plot(df_myTwiss.x*1E3, df_myTwiss.s, lw=0.5, color='g', label='cpymad Closed Orbit') 
        ax2.plot((df_myTwiss.x*1E3)/np.sqrt(df_myTwiss.betx), df_myTwiss.mux, lw=0.5, color='g', label='cpymad Closed Orbit')      
    else: 
        ax1.plot(df_myTwiss.y*1E3, df_myTwiss.s, lw=0.5, color='g', label='cpymad Closed Orbit')    
        ax2.plot((df_myTwiss.y*1E3)/np.sqrt(df_myTwiss.bety), df_myTwiss.muy, lw=0.5, color='g', label='cpymad Closed Orbit')  
    
    # plot all monitor positions from twiss    
    ax1.scatter(df[time], df.S, color='r', label='cpymad TWISS')  
    ax2.scatter(df[time]/np.sqrt(df.Beta), df.Mu, color='r', label='cpymad TWISS')    
    if horizontal: ax1.set_xlabel('x [mm]')
    else: ax1.set_xlabel('y [mm]')
        
    # plot measurements with error bars
    ax1.errorbar(df_meas[time], df_meas['S'], xerr = np.ones(len(df_meas['S'])), color='b', fmt=".", label='Measurements')
    ax2.errorbar(df_meas[time]/np.sqrt(df_meas['Beta']), df_meas['Mu'], xerr = np.ones(len(df_meas['Mu']))/np.sqrt(df_meas['Beta']), color='b', fmt=".", label='Measurements')
        
    #ax1.text(0, 0, plot_str)
    fig1.suptitle(title_str)
    
    ax1.set_ylabel('s [m]')    
    ax1.grid(which='both', ls=':', lw=0.5, color='grey')
    ax1.legend()        
    
    ax2.set_ylabel(r'$\mu_{x,y}$ [-]')    
    if horizontal: ax2.set_xlabel(r'$\frac{y}{\sqrt{\beta_y}}$ [$10^{-3}\sqrt{m}$]')
    else: ax2.set_xlabel(r'$\frac{x}{\sqrt{\beta_x}}$ [$10^{-3}\sqrt{m}$]')
    ax2.grid(which='both', ls=':', lw=0.5, color='grey')
    ax2.legend()      
    
    plt.savefig(savename, bbox_inches='tight')    
    
    return std_diff
    
def compare_measurement_and_twiss_times(twiss_file, measurement_file, savename, horizontal=True, threshold=2.0, time1='0', time2='0', title=None):
    
    # Read measurement file 
    df_meas = processed_bpm_dat_to_df(measurement_file)
    #s_list, meas_list = measurements_from_file(measurement_file, horizontal=horizontal)
    
    # Read Twiss file return dataframe #'Monitor', 'S', 'mu_x', 'x'
    df = read_twiss_generate_measured_co_df(twiss_file, horizontal, False, time1)
    
    # Read Twiss file to check orbit
    df_myTwiss = tfs.read(twiss_file)
    df_myTwiss = pandas_lowercase_all_strings(df_myTwiss)
    if horizontal:
        if (df_myTwiss['x'] == 0).all() and df_myTwiss['x'].std()==0:
            print('WARNING: compare_measurement_and_twiss: twiss file '+twiss_file+' has perfect closed orbit in this plane')
            return 10
    else:
        if (df_myTwiss['y'] == 0).all() and df_myTwiss['y'].std()==0:
            print('WARNING: compare_measurement_and_twiss: twiss file '+twiss_file+' has perfect closed orbit in this plane')
            return 10  
        
    # Iterate over measurements, save difference, compute std
    diff_list = []
    for i in range(0,len(df_meas),1):
        if np.isnan(df_meas[time2][i]): pass
        else:
            diff = df_meas[time2][i] - float(df[df['S']==df_meas['S'][i]][time1])        
            diff_list.append(diff)
    std_diff = np.std(diff_list)
        
    if std_diff >= threshold:
        return std_diff
    
    plot_str = r'$\sigma_{difference}$ = '+str(round_sig(std_diff,3))  
    if title: title_str = title + '\nStandard Deviation ' + plot_str + ' mm'
    else:    title_str = 'Standard Deviation\n' + plot_str + ' mm'
    
    fig1 = plt.figure(facecolor='w', edgecolor='k', figsize=[9, 6])
    gs = fig1.add_gridspec(ncols=2, nrows=1, width_ratios=[1,1])
    gs.update(wspace=0.25, hspace=0.)

    ax1 = fig1.add_subplot(gs[0,0])
    ax2 = fig1.add_subplot(gs[0,1])    
    
    # plot closed orbit
    if horizontal: 
        ax1.plot(df_myTwiss.x*1E3, df_myTwiss.s, lw=0.5, color='g', label='cpymad Closed Orbit') 
        ax2.plot((df_myTwiss.x*1E3)/np.sqrt(df_myTwiss.betx), df_myTwiss.mux, lw=0.5, color='g', label='cpymad Closed Orbit')      
    else: 
        ax1.plot(df_myTwiss.y*1E3, df_myTwiss.s, lw=0.5, color='g', label='cpymad Closed Orbit')    
        ax2.plot((df_myTwiss.y*1E3)/np.sqrt(df_myTwiss.bety), df_myTwiss.muy, lw=0.5, color='g', label='cpymad Closed Orbit')  
    
    # plot all monitor positions from twiss    
    ax1.scatter(df[time1], df.S, color='r', label='cpymad Expectation')  
    ax2.scatter(df[time1]/np.sqrt(df.Beta), df.Mu, color='r', label='cpymad TWISS')    
    if horizontal: ax1.set_xlabel('x [mm]')
    else: ax1.set_xlabel('y [mm]')
        
    # plot measurements with error bars
    ax1.errorbar(df_meas[time2], df_meas['S'], xerr = np.ones(len(df_meas['S'])), color='b', fmt=".", label='Measurements')
    ax2.errorbar(df_meas[time2]/np.sqrt(df_meas['Beta']), df_meas['Mu'], xerr = np.ones(len(df_meas['Mu']))/np.sqrt(df_meas['Beta']), color='b', fmt=".", label='Measurements')
        
    #ax1.text(0, 0, plot_str)
    fig1.suptitle(title_str)
    
    ax1.set_ylabel('s [m]')    
    ax1.grid(which='both', ls=':', lw=0.5, color='grey')
    ax1.legend()        
    
    ax2.set_ylabel(r'$\mu_{x,y}$ [-]')    
    if horizontal: ax2.set_xlabel(r'$\frac{y}{\sqrt{\beta_y}}$ [$10^{-3}\sqrt{m}$]')
    else: ax2.set_xlabel(r'$\frac{x}{\sqrt{\beta_x}}$ [$10^{-3}\sqrt{m}$]')
    ax2.grid(which='both', ls=':', lw=0.5, color='grey')
    ax2.legend()      
    
    plt.savefig(savename, bbox_inches='tight')    
    
    return std_diff
########################################################################
# Generate closed orbit dataframe (s, name, x, y) from twiss file
########################################################################
def read_twiss_generate_orbit_df(twissfile, horizontal=True, FMon=True):
        
    h_all = ['sp0_r0hm1','sp0_r0hm2','sp1_r1hm1','sp1_r1hm2','sp2_r2hm1','sp2_r2hm2','sp3_r3hm1','sp4_r4hm2','sp5_r5hm1','sp6_r6hm1','sp6_r6hm2','sp6_r6hm3','sp7_r7hm1','sp7_r7hm2','sp8_r8hm2','sp9_r9hm1']
    v_all = ['sp0_r0vm1','sp1_r1vm1','sp1_r1vm2','sp2_r2vm1','sp2_r2vm2','sp3_r3vm1','sp3_r3vm2','sp3_r3vm3','sp4_r4vm1','sp4_r4vm2','sp5_r5vm1','sp5_r5vm2','sp6_r6vm1','sp6_r6vm2','sp7_r7vm2','sp7_r7vm3','sp8_r8vm1','sp8_r8vm2','sp9_r9vm1','sp9_r9vm2','sp9_r9vm3']
    h_fmon = ['sp0_r0hm2','sp3_r3hm1','sp4_r4hm2','sp5_r5hm1','sp6_r6hm2','sp7_r7hm2','sp8_r8hm2','sp9_r9hm1']
    v_fmon = ['sp0_r0vm1','sp2_r2vm1','sp3_r3vm1','sp4_r4vm1','sp5_r5vm1','sp6_r6vm1','sp7_r7vm2','sp8_r8vm1','sp9_r9vm1']    
    
    t_s = []
    t_name = [] 
    t_x = []
    t_y = []
        
    df_myTwiss = tfs.read(twissfile)
    df_myTwiss = pandas_lowercase_all_strings(df_myTwiss)
    
    if horizontal:    
        if FMon:
            for element in h_fmon:
                t_s.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].s))
                t_name.append(df_myTwiss[df_myTwiss['name'].str.contains(element)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(element)].index[0])])
                t_x.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].x))
                t_y.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].y))
        else:
            for element in h_all:
                t_s.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].s))
                t_name.append(df_myTwiss[df_myTwiss['name'].str.contains(element)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(element)].index[0])])
                t_x.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].x))
                t_y.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].y))
                    
        zipped = list(zip(t_name, t_s, t_x, t_y))    
        df = pnd.DataFrame(zipped, columns=['Monitor', 'S', 'x', 'y'])
    else:    
        if FMon:
            for element in v_fmon:
                t_s.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].s))
                t_name.append(df_myTwiss[df_myTwiss['name'].str.contains(element)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(element)].index[0])])
                t_x.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].x))
                t_y.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].y))
        else:
            for element in v_all:
                t_s.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].s))
                t_name.append(df_myTwiss[df_myTwiss['name'].str.contains(element)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(element)].index[0])])
                t_x.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].x))
                t_y.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].y))
                    
        zipped = list(zip(t_name, t_s, t_x, t_y))    
        df = pnd.DataFrame(zipped, columns=['Monitor', 'S', 'x', 'y'])

    return df


########################################################################
# Generate mad-x readable closed orbit tfs table (s, name, x, y) from twiss file
########################################################################
def read_twiss_generate_orbit_tfs(twissfile, savefilename, horizontal=True, FMon=True):
    df = read_twiss_generate_orbit_df(twissfile, horizontal, FMon)
    f = open(savefilename, 'w')
    justification = 12
    
    header_string = '@ NAME             %08s "ORBIT_IN"\n@ TYPE             %05s "ORBIT"\n@ COMMENT          %05s "ORBIT"\n'
    f.write(header_string)
    
    f.write('*'+' '.join(map(lambda i: i.rjust(justification), ['NAME', 'X', 'Y'])) + '\n')
    f.write('$'+' '.join(map(lambda i: i.rjust(justification), ['%s', '%le', '%le'])) + '\n')
 
    for index, row in df.iterrows():
        mon_string = '\"'+str(row['Monitor'])+'\"'
        
        x_string = '{0:1.4e}'.format(row['x'])
        y_string = '{0:1.4e}'.format(row['y'])   
        tot_string = (' '.join(map(lambda i: i.rjust(justification), [mon_string, x_string, y_string])) + '\n')
    
        # old method
        #xy_string = (' '.rjust(justification).join(map(lambda i: ('%1.4e'%i).rjust(justification), [row['x'], row['y']])) + '\n')
        #tot_string = mon_string.rjust(justification)+xy_string
        f.write(tot_string)
    f.close()
    

########################################################################
# Generate mesurement (at monitors) closed orbit dataframe (s, mu, name, x/y) from twiss file
########################################################################
def read_twiss_generate_measured_co_df(twissfile, horizontal=True, FMon=True, time='0'):
    
    h_all = ['sp0_r0hm1','sp0_r0hm2','sp1_r1hm1','sp1_r1hm2','sp2_r2hm1','sp2_r2hm2','sp3_r3hm1','sp4_r4hm2','sp5_r5hm1','sp6_r6hm1','sp6_r6hm2','sp6_r6hm3','sp7_r7hm1','sp7_r7hm2','sp8_r8hm2','sp9_r9hm1']
    v_all = ['sp0_r0vm1','sp1_r1vm1','sp1_r1vm2','sp2_r2vm1','sp2_r2vm2','sp3_r3vm1','sp3_r3vm2','sp3_r3vm3','sp4_r4vm1','sp4_r4vm2','sp5_r5vm1','sp5_r5vm2','sp6_r6vm1','sp6_r6vm2','sp7_r7vm2','sp7_r7vm3','sp8_r8vm1','sp8_r8vm2','sp9_r9vm1','sp9_r9vm2','sp9_r9vm3']
    h_fmon = ['sp0_r0hm2','sp3_r3hm1','sp4_r4hm2','sp5_r5hm1','sp6_r6hm2','sp7_r7hm2','sp8_r8hm2','sp9_r9hm1']
    v_fmon = ['sp0_r0vm1','sp2_r2vm1','sp3_r3vm1','sp4_r4vm1','sp5_r5vm1','sp6_r6vm1','sp7_r7vm2','sp8_r8vm1','sp9_r9vm1']    
    
    t_s = []
    t_mu = []
    t_name = [] 
    t_x = []
    t_beta = []
    
    df_myTwiss = tfs.read(twissfile)
    df_myTwiss = pandas_lowercase_all_strings(df_myTwiss)
    
    if horizontal:    
        if FMon:
            for element in h_fmon:
                t_s.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].s))
                t_mu.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].mux))
                t_name.append(df_myTwiss[df_myTwiss['name'].str.contains(element)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(element)].index[0])])
                t_x.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].x)*1E3)
                t_beta.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].betx))
        else:
            for element in h_all:
                t_s.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].s))
                t_mu.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].mux))
                t_name.append(df_myTwiss[df_myTwiss['name'].str.contains(element)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(element)].index[0])])
                t_x.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].x)*1E3)
                t_beta.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].betx))
    else:    
        if FMon:
            for element in v_fmon:
                t_s.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].s))
                t_mu.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].muy))
                t_name.append(df_myTwiss[df_myTwiss['name'].str.contains(element)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(element)].index[0])])
                t_x.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].y)*1E3)
                t_beta.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].bety))
        else:
            for element in v_all:
                t_s.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].s))
                t_mu.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].muy))
                t_name.append(df_myTwiss[df_myTwiss['name'].str.contains(element)].name[int(df_myTwiss[df_myTwiss['name'].str.contains(element)].index[0])])
                t_x.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].y)*1E3)                
                t_beta.append(float(df_myTwiss[df_myTwiss['name'].str.contains(element)].bety))
                    
    zipped = list(zip(t_s, t_name, t_mu, t_beta, t_x))    
    df = pnd.DataFrame(zipped, columns=['S', 'Monitor', 'Mu', 'Beta', str(time)])

    return df

########################################################################
# Generate mesurement (at monitors) closed orbit file (s, name, x, y) from twiss file
########################################################################
def read_twiss_generate_measured_co(twissfile, savefilename, horizontal=True, FMon=True, time='0'):
    
    df = read_twiss_generate_measured_co_df(twissfile, horizontal, FMon, time)
    pandas_save_to_file(df, savefilename)
    return df
    
########################################################################
# Set all ISIS kickers to zero
########################################################################
def isis_reset_steering(madx_instance):   
    
    # Steering magnets
    madx_instance.globals.r0hd1_kick = 0.0
    madx_instance.globals.r2hd1_kick = 0.0
    madx_instance.globals.r3hd1_kick = 0.0
    madx_instance.globals.r4hd1_kick = 0.0
    madx_instance.globals.r5hd1_kick = 0.0
    madx_instance.globals.r7hd1_kick = 0.0
    madx_instance.globals.r9hd1_kick = 0.0      
    madx_instance.globals.r0vd1_kick = 0.0
    madx_instance.globals.r2vd1_kick = 0.0
    madx_instance.globals.r3vd1_kick = 0.0
    madx_instance.globals.r4vd1_kick = 0.0
    madx_instance.globals.r5vd1_kick = 0.0
    madx_instance.globals.r7vd1_kick = 0.0
    madx_instance.globals.r9vd1_kick = 0.0 
    
    # QTF Trim Quads
    madx_instance.globals.R0QTF_HKICK = 0.0
    madx_instance.globals.R1QTF_HKICK = 0.0
    madx_instance.globals.R2QTF_HKICK = 0.0
    madx_instance.globals.R3QTF_HKICK = 0.0
    madx_instance.globals.R4QTF_HKICK = 0.0
    madx_instance.globals.R5QTF_HKICK = 0.0
    madx_instance.globals.R6QTF_HKICK = 0.0      
    madx_instance.globals.R7QTF_HKICK = 0.0
    madx_instance.globals.R8QTF_HKICK = 0.0
    madx_instance.globals.R9QTF_HKICK = 0.0
    madx_instance.globals.R0QTF_VKICK = 0.0
    madx_instance.globals.R1QTF_VKICK = 0.0
    madx_instance.globals.R2QTF_VKICK = 0.0
    madx_instance.globals.R3QTF_VKICK = 0.0 
    madx_instance.globals.R4QTF_VKICK = 0.0
    madx_instance.globals.R5QTF_VKICK = 0.0
    madx_instance.globals.R6QTF_VKICK = 0.0
    madx_instance.globals.R7QTF_VKICK = 0.0
    madx_instance.globals.R8QTF_VKICK = 0.0
    madx_instance.globals.R9QTF_VKICK = 0.0
    
    # QTF Trim Quads
    madx_instance.globals.R0QTD_HKICK = 0.0
    madx_instance.globals.R1QTD_HKICK = 0.0
    madx_instance.globals.R2QTD_HKICK = 0.0
    madx_instance.globals.R3QTD_HKICK = 0.0
    madx_instance.globals.R4QTD_HKICK = 0.0
    madx_instance.globals.R5QTD_HKICK = 0.0
    madx_instance.globals.R6QTD_HKICK = 0.0      
    madx_instance.globals.R7QTD_HKICK = 0.0
    madx_instance.globals.R8QTD_HKICK = 0.0
    madx_instance.globals.R9QTD_HKICK = 0.0
    madx_instance.globals.R0QTD_VKICK = 0.0
    madx_instance.globals.R1QTD_VKICK = 0.0
    madx_instance.globals.R2QTD_VKICK = 0.0
    madx_instance.globals.R3QTD_VKICK = 0.0 
    madx_instance.globals.R4QTD_VKICK = 0.0
    madx_instance.globals.R5QTD_VKICK = 0.0
    madx_instance.globals.R6QTD_VKICK = 0.0
    madx_instance.globals.R7QTD_VKICK = 0.0
    madx_instance.globals.R8QTD_VKICK = 0.0
    madx_instance.globals.R9QTD_VKICK = 0.0
    
    # QD Doublet Quads
    madx_instance.globals.R0QF_HKICK = 0.0
    madx_instance.globals.R1QF_HKICK = 0.0
    madx_instance.globals.R2QF_HKICK = 0.0
    madx_instance.globals.R3QF_HKICK = 0.0
    madx_instance.globals.R4QF_HKICK = 0.0
    madx_instance.globals.R5QF_HKICK = 0.0
    madx_instance.globals.R6QF_HKICK = 0.0      
    madx_instance.globals.R7QF_HKICK = 0.0
    madx_instance.globals.R8QF_HKICK = 0.0
    madx_instance.globals.R9QF_HKICK = 0.0
    madx_instance.globals.R0QF_VKICK = 0.0
    madx_instance.globals.R1QF_VKICK = 0.0
    madx_instance.globals.R2QF_VKICK = 0.0
    madx_instance.globals.R3QF_VKICK = 0.0 
    madx_instance.globals.R4QF_VKICK = 0.0
    madx_instance.globals.R5QF_VKICK = 0.0
    madx_instance.globals.R6QF_VKICK = 0.0
    madx_instance.globals.R7QF_VKICK = 0.0
    madx_instance.globals.R8QF_VKICK = 0.0
    madx_instance.globals.R9QF_VKICK = 0.0
    
    # QF Doublet Quads
    madx_instance.globals.R0QD_HKICK = 0.0
    madx_instance.globals.R1QD_HKICK = 0.0
    madx_instance.globals.R2QD_HKICK = 0.0
    madx_instance.globals.R3QD_HKICK = 0.0
    madx_instance.globals.R4QD_HKICK = 0.0
    madx_instance.globals.R5QD_HKICK = 0.0
    madx_instance.globals.R6QD_HKICK = 0.0      
    madx_instance.globals.R7QD_HKICK = 0.0
    madx_instance.globals.R8QD_HKICK = 0.0
    madx_instance.globals.R9QD_HKICK = 0.0
    madx_instance.globals.R0QD_VKICK = 0.0
    madx_instance.globals.R1QD_VKICK = 0.0
    madx_instance.globals.R2QD_VKICK = 0.0
    madx_instance.globals.R3QD_VKICK = 0.0 
    madx_instance.globals.R4QD_VKICK = 0.0
    madx_instance.globals.R5QD_VKICK = 0.0
    madx_instance.globals.R6QD_VKICK = 0.0
    madx_instance.globals.R7QD_VKICK = 0.0
    madx_instance.globals.R8QD_VKICK = 0.0
    madx_instance.globals.R9QD_VKICK = 0.0
    
    # QDS Singlet Quads
    madx_instance.globals.R0QDS_HKICK = 0.0
    madx_instance.globals.R1QDS_HKICK = 0.0
    madx_instance.globals.R2QDS_HKICK = 0.0
    madx_instance.globals.R3QDS_HKICK = 0.0
    madx_instance.globals.R4QDS_HKICK = 0.0
    madx_instance.globals.R5QDS_HKICK = 0.0
    madx_instance.globals.R6QDS_HKICK = 0.0      
    madx_instance.globals.R7QDS_HKICK = 0.0
    madx_instance.globals.R8QDS_HKICK = 0.0
    madx_instance.globals.R9QDS_HKICK = 0.0
    madx_instance.globals.R0QDS_VKICK = 0.0
    madx_instance.globals.R1QDS_VKICK = 0.0
    madx_instance.globals.R2QDS_VKICK = 0.0
    madx_instance.globals.R3QDS_VKICK = 0.0 
    madx_instance.globals.R4QDS_VKICK = 0.0
    madx_instance.globals.R5QDS_VKICK = 0.0
    madx_instance.globals.R6QDS_VKICK = 0.0
    madx_instance.globals.R7QDS_VKICK = 0.0
    madx_instance.globals.R8QDS_VKICK = 0.0
    madx_instance.globals.R9QDS_VKICK = 0.0
    
    # MDS Start Main Dipole
    madx_instance.globals.R0MDS_HKICK = 0.0
    madx_instance.globals.R1MDS_HKICK = 0.0
    madx_instance.globals.R2MDS_HKICK = 0.0
    madx_instance.globals.R3MDS_HKICK = 0.0
    madx_instance.globals.R4MDS_HKICK = 0.0
    madx_instance.globals.R5MDS_HKICK = 0.0
    madx_instance.globals.R6MDS_HKICK = 0.0      
    madx_instance.globals.R7MDS_HKICK = 0.0
    madx_instance.globals.R8MDS_HKICK = 0.0
    madx_instance.globals.R9MDS_HKICK = 0.0
    madx_instance.globals.R0MDS_VKICK = 0.0
    madx_instance.globals.R1MDS_VKICK = 0.0
    madx_instance.globals.R2MDS_VKICK = 0.0
    madx_instance.globals.R3MDS_VKICK = 0.0 
    madx_instance.globals.R4MDS_VKICK = 0.0
    madx_instance.globals.R5MDS_VKICK = 0.0
    madx_instance.globals.R6MDS_VKICK = 0.0
    madx_instance.globals.R7MDS_VKICK = 0.0
    madx_instance.globals.R8MDS_VKICK = 0.0
    madx_instance.globals.R9MDS_VKICK = 0.0
    
    # MDE End Main Dipole
    madx_instance.globals.R0MDE_HKICK = 0.0
    madx_instance.globals.R1MDE_HKICK = 0.0
    madx_instance.globals.R2MDE_HKICK = 0.0
    madx_instance.globals.R3MDE_HKICK = 0.0
    madx_instance.globals.R4MDE_HKICK = 0.0
    madx_instance.globals.R5MDE_HKICK = 0.0
    madx_instance.globals.R6MDE_HKICK = 0.0      
    madx_instance.globals.R7MDE_HKICK = 0.0
    madx_instance.globals.R8MDE_HKICK = 0.0
    madx_instance.globals.R9MDE_HKICK = 0.0
    madx_instance.globals.R0MDE_VKICK = 0.0
    madx_instance.globals.R1MDE_VKICK = 0.0
    madx_instance.globals.R2MDE_VKICK = 0.0
    madx_instance.globals.R3MDE_VKICK = 0.0 
    madx_instance.globals.R4MDE_VKICK = 0.0
    madx_instance.globals.R5MDE_VKICK = 0.0
    madx_instance.globals.R6MDE_VKICK = 0.0
    madx_instance.globals.R7MDE_VKICK = 0.0
    madx_instance.globals.R8MDE_VKICK = 0.0
    madx_instance.globals.R9MDE_VKICK = 0.0

########################################################################
########################################################################
#                            MATCHING
########################################################################
########################################################################
def measurements_from_file(filename, time='0'):
    # ~ s_list = []
    # ~ meas_list = []    
    # ~ time_list = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    m_df = pnd.read_table(filename, index_col=0, header=0, delim_whitespace=True)
    # ~ time_label = time
    # ~ time_test = False
    
    # ~ if time not in m_df.keys():
        # ~ print('measurements_from_file selected time not in measurement file, finding available times')
        # ~ for t in time_list:
            # ~ if t in m_df.keys(): 
                # ~ time = t
                # ~ print('measurements_from_file: new time selected as ',t)
                # ~ break
            # ~ else:
                # ~ continue                   

    # ~ for index, row in m_df.iterrows():   
        # ~ if np.isnan(row[time]): print('Skipping constraint at ', row['Monitor'], ' as value is ', row[time])
        # ~ else:                 
            # ~ s_list.append(row['S'])
            # ~ meas_list.append(row[time])                        
    #return s_list, meas_list  
    return m_df
    
def measurements_from_file_old1(filename, horizontal=True):
    s_list = []
    meas_list = []    

    m_df = pnd.read_table(filename, index_col=0, header=0, delim_whitespace=True)

    for index, row in m_df.iterrows(): 
        if horizontal:                
            if 'h' in row['Monitor']:
                if np.isnan(row['x']): print('Skipping constraint at ', row['Monitor'], ' as value is ', row['x'])
                else:                 
                    s_list.append(row['S'])
                    meas_list.append(row['x'])
        else:
            if 'v' in row['Monitor']:
                if np.isnan(row['y']): print('Skipping constraint at ', row['Monitor'], ' as value is ', row['y'])
                else: 
                    s_list.append(row['S'])
                    meas_list.append(row['y'])
                        
    return s_list, meas_list     
    
def measurements_from_file_old2(filename, monitor_header='Monitor', horizontal=True, filetype='BPM', time='0'):
    s_list = []
    meas_list = []    
    
    if filetype == 'BPM':
        # Read file to dataframe
        bpm_df = pnd.read_table(filename, delim_whitespace=True)
        
        for index, row in bpm_df.iterrows(): 
            if horizontal:
                if 'h' in row[monitor_header]:
                    if np.isnan(row[time]): print('Skipping constraint at ', row[monitor_header], ' as value is ', row[time], 'at time ', time)
                    else: 
                        s_list.append(row['S'])
                        meas_list.append(row[time])
            else:
                if 'v' in row[monitor_header]:
                    if np.isnan(row[time]): print('Skipping constraint at ', row[monitor_header], ' as value is ', row[time], 'at time ', time)
                    else: 
                        s_list.append(row['S'])
                        meas_list.append(row[time])
    else:
        # Read file to dataframe
        m_df = pnd.read_table(filename, index_col=0, header=0, delim_whitespace=True)
        
        for index, row in m_df.iterrows(): 
            if horizontal:                
                if 'h' in row[monitor_header]:
                    if np.isnan(row['x']): print('Skipping constraint at ', row[monitor_header], ' as value is ', row['x'])
                    else:                 
                        s_list.append(row['S'])
                        meas_list.append(row['x'])
            else:
                if 'v' in row[monitor_header]:
                    if np.isnan(row['y']): print('Skipping constraint at ', row[monitor_header], ' as value is ', row['y'])
                    else: 
                        s_list.append(row['S'])
                        meas_list.append(row['y'])
                        
    return s_list, meas_list                                         


def constraints_from_file_range(filename, horizontal=True, time='0'):
    constraints_list = []
    
    # Read file to dataframe
    bpm_df = pnd.read_table(filename, delim_whitespace=True)
    
    for index, row in bpm_df.iterrows(): 
        if horizontal:
            if 'h' in row['Monitor']:
                if np.isnan(row[time]): print('Skipping constraint at ', row['Monitor'], ' as value is ', row[time], 'at time ', time)
                else: constraints_list.append(dict(range=row['Monitor'], x=Constraint(min=(row[time]-0.5)*1E-3, max=(row[time]+0.5)*1E-3)))
        else:
            if 'v' in row['Monitor']:
                if np.isnan(row[time]): print('Skipping constraint at ', row['Monitor'], ' as value is ', row[time], 'at time ', time)
                else: constraints_list.append(dict(range=row['Monitor'], y=Constraint(min=(row[time]-0.5)*1E-3, max=(row[time]+0.5)*1E-3)))

    return constraints_list
    
def constraints_from_file_range_old(filename, monitor_header='Monitor', horizontal=True, filetype='BPM', time='0'):
    constraints_list = []
    
    if filetype == 'BPM':
        # Read file to dataframe
        bpm_df = pnd.read_table(filename, delim_whitespace=True)
        
        for index, row in bpm_df.iterrows(): 
            if horizontal:
                if 'h' in row[monitor_header]:
                    if np.isnan(row[time]): print('Skipping constraint at ', row[monitor_header], ' as value is ', row[time], 'at time ', time)
                    else: constraints_list.append(dict(range=row[monitor_header], x=Constraint(min=(row[time]-0.5)*1E-3, max=(row[time]+0.5)*1E-3)))
            else:
                if 'v' in row[monitor_header]:
                    if np.isnan(row[time]): print('Skipping constraint at ', row[monitor_header], ' as value is ', row[time], 'at time ', time)
                    else: constraints_list.append(dict(range=row[monitor_header], y=Constraint(min=(row[time]-0.5)*1E-3, max=(row[time]+0.5)*1E-3)))

    else:
        # Read file to dataframe
        m_df = pnd.read_table(filename, index_col=0, header=0, delim_whitespace=True)
        
        for index, row in m_df.iterrows(): 

            if horizontal:                
                if 'h' in row[monitor_header]:
                    if np.isnan(row['x']): print('Skipping constraint at ', row[monitor_header], ' as value is ', row['x'])
                    else: constraints_list.append(dict(range=row[monitor_header], x=Constraint(min=(row['x']-0.5)*1E-3, max=(row['x']+0.5)*1E-3)))
            else:
                if 'v' in row[monitor_header]:
                    if np.isnan(row['y']): print('Skipping constraint at ', row[monitor_header], ' as value is ', row['y'])
                    else: constraints_list.append(dict(range=row[monitor_header], y=Constraint(min=(row['y']-0.5)*1E-3, max=(row['y']+0.5)*1E-3)))                     

    return constraints_list                                         
                    
def constraints_from_file(filename, horizontal=True):
    constraints_list = []

    # Read file to dataframe
    bpm_df = pnd.read_table(filename, delim_whitespace=True)
    
    for index, row in bpm_df.iterrows(): 
        if horizontal:
            if 'h' in row['Monitor']:
                if np.isnan(row[time]): print('Skipping constraint at ', row['Monitor'], ' as value is ', row[time], 'at time ', time)
                else: constraints_list.append(dict(range=row['Monitor'], x=row[time]*1E-3))
        else:
            if 'v' in row['Monitor']:
                if np.isnan(row[time]): print('Skipping constraint at ', row['Monitor'], ' as value is ', row[time], 'at time ', time)
                else: constraints_list.append(dict(range=row['Monitor'], y=row[time]*1E-3))

    return constraints_list                              

def constraints_from_file_old(filename, monitor_header='Monitor', horizontal=True, filetype='BPM', time='0'):
    constraints_list = []
    
    if filetype == 'BPM':
        # Read file to dataframe
        bpm_df = pnd.read_table(filename, delim_whitespace=True)
        
        for index, row in bpm_df.iterrows(): 
            if horizontal:
                if 'h' in row[monitor_header]:
                    if np.isnan(row[time]): print('Skipping constraint at ', row[monitor_header], ' as value is ', row[time], 'at time ', time)
                    else: constraints_list.append(dict(range=row[monitor_header], x=row[time]*1E-3))
            else:
                if 'v' in row[monitor_header]:
                    if np.isnan(row[time]): print('Skipping constraint at ', row[monitor_header], ' as value is ', row[time], 'at time ', time)
                    else: constraints_list.append(dict(range=row[monitor_header], y=row[time]*1E-3))

    else:
        # Read file to dataframe
        m_df = pnd.read_table(filename, index_col=0, header=0, delim_whitespace=True)
        
        for index, row in m_df.iterrows(): 

            if horizontal:                
                if 'h' in row[monitor_header]:
                    if np.isnan(row['x']): print('Skipping constraint at ', row[monitor_header], ' as value is ', row['x'])
                    else: constraints_list.append(dict(range=row[monitor_header], x=row['x']*1E-3))
            else:
                if 'v' in row[monitor_header]:                    
                    if np.isnan(row['y']): print('Skipping constraint at ', row[monitor_header], ' as value is ', row['y'])
                    else: constraints_list.append(dict(range=row[monitor_header], y=row['y']*1E-3))                   

    return constraints_list                                         
                                            

def match_co_from_file_range(madx_instance, filename, time='0', sp_list=['0','1','3','4','5','6','7','8','9'], horizontal=True, steering=True, dip_ends=False, doublet_quad=False, trim_quad=False, singlet_quad=False):
    
    madx_instance.command.match(chrom=False)    
    
    if horizontal:                
        if '0' in sp_list:
            if steering:     madx_instance.command.vary(name='r0hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R0MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R0MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R0QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R0QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R0QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R0QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R0QDS_HKICK', step=1E-4)
        if '1' in sp_list:
            if steering:     madx_instance.command.vary(name='r1hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R1MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R1MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R1QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R1QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R1QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R1QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R1QDS_HKICK', step=1E-4)         
        if '2' in sp_list:
            if steering:     madx_instance.command.vary(name='r2hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R2MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R2MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R2QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R2QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R2QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R2QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R2QDS_HKICK', step=1E-4)       
        if '3' in sp_list:
            if steering:     madx_instance.command.vary(name='r3hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R3MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R3MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R3QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R3QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R3QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R3QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R3QDS_HKICK', step=1E-4)         
        if '4' in sp_list:
            if steering:     madx_instance.command.vary(name='r4hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R4MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R4MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R4QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R4QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R4QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R4QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R4QDS_HKICK', step=1E-4)           
        if '5' in sp_list:
            if steering:     madx_instance.command.vary(name='r5hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R5MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R5MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R5QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R5QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R5QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R5QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R5QDS_HKICK', step=1E-4)          
        if '6' in sp_list:
            if steering:     madx_instance.command.vary(name='r6hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R6MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R6MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R6QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R6QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R6QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R6QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R6QDS_HKICK', step=1E-4)        
        if '7' in sp_list:
            if steering:     madx_instance.command.vary(name='r7hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R7MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R7MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R7QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R7QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R7QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R7QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R7QDS_HKICK', step=1E-4)                       
        if '8' in sp_list:
            if steering:     madx_instance.command.vary(name='r8hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R8MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R8MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R8QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R8QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R8QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R8QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R8QDS_HKICK', step=1E-4)                                     
        if '9' in sp_list:
            if steering:     madx_instance.command.vary(name='r9hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R9MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R9MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R9QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R9QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R9QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R9QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R9QDS_HKICK', step=1E-4)  
        
    else:
        if '0' in sp_list:
            if steering:     madx_instance.command.vary(name='r0vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R0MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R0MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R0QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R0QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R0QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R0QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R0QDS_VKICK', step=1E-4)
        if '1' in sp_list:
            if steering:     madx_instance.command.vary(name='r1vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R1MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R1MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R1QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R1QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R1QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R1QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R1QDS_VKICK', step=1E-4)           
        if '2' in sp_list:
            if steering:     madx_instance.command.vary(name='r2vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R2MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R2MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R2QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R2QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R2QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R2QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R2QDS_VKICK', step=1E-4)    
        if '3' in sp_list:
            if steering:     madx_instance.command.vary(name='r3vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R3MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R3MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R3QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R3QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R3QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R3QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R3QDS_VKICK', step=1E-4)            
        if '4' in sp_list:
            if steering:     madx_instance.command.vary(name='r4vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R4MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R4MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R4QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R4QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R4QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R4QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R4QDS_VKICK', step=1E-4)       
        if '5' in sp_list:
            if steering:     madx_instance.command.vary(name='r5vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R5MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R5MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R5QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R5QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R5QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R5QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R5QDS_VKICK', step=1E-4)          
        if '6' in sp_list:
            if steering:     madx_instance.command.vary(name='r6vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R6MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R6MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R6QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R6QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R6QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R6QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R6QDS_VKICK', step=1E-4)            
        if '7' in sp_list:
            if steering:     madx_instance.command.vary(name='r7vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R7MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R7MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R7QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R7QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R7QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R7QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R7QDS_VKICK', step=1E-4)                        
        if '8' in sp_list:
            if steering:     madx_instance.command.vary(name='r8vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R8MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R8MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R8QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R8QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R8QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R8QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R8QDS_VKICK', step=1E-4)                                    
        if '9' in sp_list:
            if steering:     madx_instance.command.vary(name='r9vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R9MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R9MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R9QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R9QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R9QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R9QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R9QDS_VKICK', step=1E-4)    

    constraints = constraints_from_file_range(filename=filename, time=time, horizontal=horizontal)
            
    # Apply constraints            
    for c in constraints: madx_instance.command.constraint(**c)
    
    # Jacobian Matching
    madx_instance.command.jacobian(calls=50000, tolerance=1e-6)
    madx_instance.command.endmatch()


def match_co_from_file(madx_instance, filename, sp_list=['0','1','3','4','5','6','7','8','9'], horizontal=True, steering=True, dip_ends=False, doublet_quad=False, trim_quad=False, singlet_quad=False):
    
    madx_instance.command.match(chrom=False)    
    
    if horizontal:                
        if '0' in sp_list:
            if steering:     madx_instance.command.vary(name='r0hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R0MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R0MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R0QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R0QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R0QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R0QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R0QDS_HKICK', step=1E-4)
        if '1' in sp_list:
            if steering:     madx_instance.command.vary(name='r1hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R1MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R1MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R1QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R1QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R1QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R1QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R1QDS_HKICK', step=1E-4)         
        if '2' in sp_list:
            if steering:     madx_instance.command.vary(name='r2hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R2MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R2MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R2QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R2QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R2QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R2QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R2QDS_HKICK', step=1E-4)       
        if '3' in sp_list:
            if steering:     madx_instance.command.vary(name='r3hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R3MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R3MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R3QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R3QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R3QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R3QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R3QDS_HKICK', step=1E-4)         
        if '4' in sp_list:
            if steering:     madx_instance.command.vary(name='r4hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R4MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R4MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R4QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R4QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R4QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R4QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R4QDS_HKICK', step=1E-4)           
        if '5' in sp_list:
            if steering:     madx_instance.command.vary(name='r5hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R5MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R5MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R5QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R5QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R5QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R5QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R5QDS_HKICK', step=1E-4)          
        if '6' in sp_list:
            if steering:     madx_instance.command.vary(name='r6hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R6MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R6MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R6QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R6QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R6QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R6QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R6QDS_HKICK', step=1E-4)        
        if '7' in sp_list:
            if steering:     madx_instance.command.vary(name='r7hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R7MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R7MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R7QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R7QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R7QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R7QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R7QDS_HKICK', step=1E-4)                       
        if '8' in sp_list:
            if steering:     madx_instance.command.vary(name='r8hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R8MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R8MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R8QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R8QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R8QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R8QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R8QDS_HKICK', step=1E-4)                                     
        if '9' in sp_list:
            if steering:     madx_instance.command.vary(name='r9hd1_kick', step=1E-4)
            if dip_ends:  
                madx_instance.command.vary(name='R9MDS_HKICK', step=1E-4)
                madx_instance.command.vary(name='R9MDE_HKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R9QTF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R9QTD_HKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R9QF_HKICK', step=1E-4)
                madx_instance.command.vary(name='R9QD_HKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R9QDS_HKICK', step=1E-4)  
        
    else:
        if '0' in sp_list:
            if steering:     madx_instance.command.vary(name='r0vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R0MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R0MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R0QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R0QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R0QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R0QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R0QDS_VKICK', step=1E-4)
        if '1' in sp_list:
            if steering:     madx_instance.command.vary(name='r1vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R1MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R1MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R1QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R1QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R1QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R1QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R1QDS_VKICK', step=1E-4)           
        if '2' in sp_list:
            if steering:     madx_instance.command.vary(name='r2vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R2MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R2MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R2QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R2QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R2QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R2QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R2QDS_VKICK', step=1E-4)    
        if '3' in sp_list:
            if steering:     madx_instance.command.vary(name='r3vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R3MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R3MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R3QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R3QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R3QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R3QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R3QDS_VKICK', step=1E-4)            
        if '4' in sp_list:
            if steering:     madx_instance.command.vary(name='r4vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R4MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R4MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R4QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R4QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R4QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R4QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R4QDS_VKICK', step=1E-4)       
        if '5' in sp_list:
            if steering:     madx_instance.command.vary(name='r5vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R5MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R5MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R5QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R5QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R5QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R5QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R5QDS_VKICK', step=1E-4)          
        if '6' in sp_list:
            if steering:     madx_instance.command.vary(name='r6vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R6MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R6MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R6QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R6QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R6QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R6QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R6QDS_VKICK', step=1E-4)            
        if '7' in sp_list:
            if steering:     madx_instance.command.vary(name='r7vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R7MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R7MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R7QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R7QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R7QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R7QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R7QDS_VKICK', step=1E-4)                        
        if '8' in sp_list:
            if steering:     madx_instance.command.vary(name='r8vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R8MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R8MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R8QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R8QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R8QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R8QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R8QDS_VKICK', step=1E-4)                                    
        if '9' in sp_list:
            if steering:     madx_instance.command.vary(name='r9vd1_kick', step=1E-4)
            if dip_ends:    
                madx_instance.command.vary(name='R9MDS_VKICK', step=1E-4)
                madx_instance.command.vary(name='R9MDE_VKICK', step=1E-4)
            if trim_quad:    
                madx_instance.command.vary(name='R9QTF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R9QTD_VKICK', step=1E-4)
            if doublet_quad:    
                madx_instance.command.vary(name='R9QF_VKICK', step=1E-4)
                madx_instance.command.vary(name='R9QD_VKICK', step=1E-4)
            if singlet_quad: madx_instance.command.vary(name='R9QDS_VKICK', step=1E-4)    

    constraints = constraints_from_file(filename, horizontal)
            
    # Apply constraints            
    for c in constraints: madx_instance.command.constraint(**c)
    
    # Jacobian Matching
    madx_instance.command.jacobian(calls=50000, tolerance=1e-6)
    madx_instance.command.endmatch()


def kick_name_list_from_twissfile(madx_instance, twissfile, horizontal=True):

    name_list = kicker_list_from_twissfile(twissfile, horizontal=horizontal)
    kicker_name_list = []
    
    for element in name_list:
        #print(element)
        input_str = 'kicker_name_list.append(str(madx_instance.elements.'+element+'.defs.kick))'
        try:
            exec(input_str)
        except AttributeError:
            pass
            if horizontal:
                input_str = 'kicker_name_list.append(str(madx_instance.elements.'+element+'.defs.hkick))'
                try:                    
                    exec(input_str)
                except AttributeError:
                    print('kick_name_list_from_twissfile::AttributeError: no kick or hkick')
            else:
                input_str = 'kicker_name_list.append(str(madx_instance.elements.'+element+'.defs.vkick))'
                try:
                    exec(input_str)
                except AttributeError:
                    print('kick_name_list_from_twissfile::AttributeError: no kick or vkick')
            
    return kicker_name_list


def kicker_list_from_twissfile(twissfile, horizontal=True):

    hname_list = []
    vname_list = []
    name_list = []

    df_myTwiss = tfs.read(twissfile)
    df_myTwiss = pandas_lowercase_all_strings(df_myTwiss)
    df_kickers = df_myTwiss[df_myTwiss['keyword'].str.contains('kicker')]
    df_vkickers = df_myTwiss[df_myTwiss['keyword'].str.contains('vkicker')]
    df_hkickers = df_myTwiss[df_myTwiss['keyword'].str.contains('hkicker')]

    # Make a list of hkickers and vkickers to use for exclusion
    for index, row in df_hkickers.iterrows(): 
        hname_list.append(row['name'])
    for index, row in df_vkickers.iterrows(): 
        vname_list.append(row['name'])      


    # Iterate over kicker dataframe, store kicker names where appropriate
    for index, row in df_kickers.iterrows(): 
        if horizontal:
            if row['name'] in vname_list:
                pass
            else:
                name_list.append(row['name'])  
        else:
            if row['name'] in hname_list:
                pass
            else:
                name_list.append(row['name'])    
            
    return name_list


def kick_name_list_from_smooth_twiss(madx_instance, smooth_twiss_df, min_kick_abs = 0.0001, horizontal=True):
    smooth_twiss_df=pandas_lowercase_all_strings(smooth_twiss_df)
    if horizontal:
        el_list = smooth_twiss_df[smooth_twiss_df['hkick'].abs()>min_kick_abs].index
    else:
        el_list = smooth_twiss_df[smooth_twiss_df['vkick'].abs()>min_kick_abs].index
    
    element_list = []
    for el in el_list:
        element_list.append(el)

    kicker_name_list = []
    for element in element_list:
        #print(element)
        input_str = 'kicker_name_list.append(str(madx_instance.elements.'+element+'.defs.kick))'
        try:
            exec(input_str)
        except AttributeError:
            pass
            if horizontal:
                input_str = 'kicker_name_list.append(str(madx_instance.elements.'+element+'.defs.hkick))'
                try:                    
                    exec(input_str)
                except AttributeError:
                    print('kick_name_list_from_twissfile::AttributeError: no kick or hkick')
            else:
                input_str = 'kicker_name_list.append(str(madx_instance.elements.'+element+'.defs.vkick))'
                try:
                    exec(input_str)
                except AttributeError:
                    print('kick_name_list_from_twissfile::AttributeError: no kick or vkick')
            
    return kicker_name_list


def match_single_kick_from_file(madx_instance, varyname, filename, horizontal=True):

    madx_instance.command.match(chrom=False)    
    
    madx_instance.command.vary(name=varyname, step=1E-4)
    
    constraints = constraints_from_file(filename, horizontal)
    #constraints = constraints_from_file_range(filename, horizontal)
    
    print(constraints)
    # Apply constraints            
    for c in constraints: madx_instance.command.constraint(**c)
    
    # Jacobian Matching
    madx_instance.command.jacobian(calls=50000, tolerance=1e-6)
    madx_instance.command.endmatch()

def match_single_kick_from_file_and_constraints(madx_instance, varyname, constraints):

    madx_instance.command.match(chrom=False)  
    madx_instance.command.vary(name=varyname, step=1E-4)

    # Apply constraints            
    for c in constraints: madx_instance.command.constraint(**c)
    
    # Jacobian Matching
    madx_instance.command.jacobian(calls=50000, tolerance=1e-6)
    madx_instance.command.endmatch()


def match_double_kick_from_file_and_constraints(madx_instance, varyname1, varyname2, constraints):

    madx_instance.command.match(chrom=False)  
    madx_instance.command.vary(name=varyname1, step=1E-4)
    madx_instance.command.vary(name=varyname2, step=1E-4)

    # Apply constraints            
    for c in constraints: madx_instance.command.constraint(**c)
    
    # Jacobian Matching
    madx_instance.command.jacobian(calls=50000, tolerance=1e-6)
    madx_instance.command.endmatch()


def find_single_kick_match(madx_instance, cpymad_logfile, sequence_name, vary_list, measurement_file, save_dir, horizontal=True, threshold=2.0):
    
    # For output dataframe
    odf_s = []
    odf_x = [] 
    odf_y = [] 
    odf_name=[]
    odf_knob=[]
    odf_kick=[]
    odf_std=[]
    
    # Make directory if it doesn't exist (function checks)
    make_directory(save_dir)

    # Make constraints from file
    constraints = constraints_from_file(measurement_file, horizontal)
    
    # Make kicker list
    twissfile = str(save_dir+'initial_twiss_.tfs')
    initial_twiss = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, twissfile)
    kicker_list = kicker_list_from_twissfile(twissfile, horizontal)
    
    # Iterate over knobs (variables used to vary kick strength)
    i = 0
    for knob in vary_list:
        print('\nMatching using element: '+knob.split('_')[0])
        
        # Reset all kickers
        isis_reset_steering(madx_instance)
        
        # Attempt match
        match_single_kick_from_file_and_constraints(madx_instance, knob, constraints)
        
        # Perform twiss
        twissfile = str(save_dir+'matched_twiss_'+str(i)+'_'+knob.split('_')[0]+'.tfs')
        matched_twiss = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, twissfile)
        matched_twiss = pandas_lowercase_all_strings(matched_twiss)
        
        # Check and plot difference
        match_std = compare_measurement_and_twiss(twissfile, measurement_file, str(save_dir+'cf_co_'+str(i)+'_'+knob.split('_')[0]+'.png'), horizontal, threshold);
        print('\t\tKnob '+knob.split('_')[0]+' standard deviation = '+str(match_std))
        
        if match_std >= threshold:
            os.remove(twissfile)
        else:
            # Store data for output
            odf_s.append(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].s[0])
            odf_name.append(kicker_list[i])
            odf_knob.append(knob)
            odf_std.append(match_std)
            if horizontal:
                odf_x.append(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].x[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])]*1E3)
                odf_kick.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].hkick[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])],3)*1E3)
            else:
                odf_y.append(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].y[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])]*1E3)
                odf_kick.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].vkick[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])],3)*1E3)
            # Plot CO with kick
            cpymad_plot_CO_kicks_from_file(madx, matched_twiss, ('Single Kick at '+knob.split('_')[0]), str(save_dir+'Kicks_'+str(i)+'_'+knob.split('_')[0]+'.png'), kicker_list, measurement_file, horizontal);  
            
        # increment label index
        i += 1
        
        print('End matching using element: '+knob.split('_')[0])
    
    if horizontal: 
        zipped = list(zip(odf_s, odf_name, odf_knob, odf_x, odf_kick, odf_std))    
        df = pnd.DataFrame(zipped, columns=['s', 'name', 'knob', 'x (mm)', 'kick (mrad)', 'std (mm)'])
    else:
        zipped = list(zip(odf_s, odf_name, odf_knob, odf_y, odf_kick, odf_std))    
        df = pnd.DataFrame(zipped, columns=['s', 'name', 'knob', 'y (mm)', 'kick (mrad)', 'std (mm)'])
    return df


def knob_to_name(knob):
    if 'hd' in knob:
        return('sp'+knob[1]+'_'+knob.split('_')[0])
    elif 'vd' in knob:
        return('sp'+knob[1]+'_'+knob.split('_')[0])
    else:
        return ('sp'+knob[1]+'_'+knob.split('_')[0][2:]+'k')


def find_double_kick_match(madx_instance, cpymad_logfile, sequence_name, vary_list, measurement_file, save_dir, horizontal=True, threshold=2.0):
    
    # For output dataframe
    odf_s1 = []
    odf_s2 = []
    odf_name1=[]
    odf_name2=[]
    odf_knob1=[]
    odf_knob2=[]
    odf_kick1=[]
    odf_kick2=[]
    odf_x1 = [] 
    odf_x2 = [] 
    odf_y1 = [] 
    odf_y2 = [] 
    odf_std=[]
    
    # Make directory if it doesn't exist (function checks)
    make_directory(save_dir)

    # Make constraints from file
    constraints = constraints_from_file(measurement_file, horizontal)
    
    # Make kicker list
    twissfile = str(save_dir+'initial_twiss_.tfs')
    initial_twiss = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, twissfile)
    #kicker_list = kicker_list_from_twissfile(twissfile, horizontal)
    
    i = 0
    itertool_combinations = combinations(vary_list, 2)
    for knobs in list(itertool_combinations):
        # first vary name is knobs[0], second is knobs[1]

        print('\nMatching using elements: '+knob_to_name(knobs[0])+' and ' +knob_to_name(knobs[1]))
        
        # Reset all kickers
        isis_reset_steering(madx_instance)
        
        # Attempt match
        match_double_kick_from_file_and_constraints(madx_instance, knobs[0], knobs[1], constraints)
        
        # Perform twiss
        twissfile = str(save_dir+'matched_twiss_'+str(i)+'_'+knobs[0].split('_')[0]+'_'+knobs[1].split('_')[0]+'.tfs')
        matched_twiss = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, twissfile)
        matched_twiss = pandas_lowercase_all_strings(matched_twiss)
        
        # Check and plot difference
        match_plot_name = str(save_dir+'cf_co_'+str(i)+'_'+knobs[0].split('_')[0]+'_'+knobs[1].split('_')[0]+'.png')
        match_std = compare_measurement_and_twiss(twissfile, measurement_file, match_plot_name, horizontal, threshold);
        print('\t\tKnobs '+knobs[0].split('_')[0]+' and ' +knobs[1].split('_')[0]+' standard deviation = '+str(match_std))
        
        if match_std >= threshold:
            os.remove(twissfile)
        else:
            # Store data for output
            odf_s1.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].s[0])
            odf_s2.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].s[0])
            odf_name1.append(knob_to_name(knobs[0]))
            odf_name2.append(knob_to_name(knobs[1]))
            odf_knob1.append(knobs[0])
            odf_knob2.append(knobs[1])
            odf_std.append(match_std)
            if horizontal:
                odf_x1.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].x[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].index[0])]*1E3)
                odf_x2.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].x[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].index[0])]*1E3)
                odf_kick1.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].hkick[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].index[0])],3)*1E3)
                odf_kick2.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].hkick[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].index[0])],3)*1E3)
            else:
                odf_y1.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].y[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].index[0])]*1E3)
                odf_y2.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].y[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].index[0])]*1E3)
                odf_kick1.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].vkick[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].index[0])],3)*1E3)
                odf_kick2.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].vkick[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].index[0])],3)*1E3)
            # Plot CO with kick
            kick_savename = str(save_dir+'Kicks_'+str(i)+'_'+knob_to_name(knobs[0])+'_' +knob_to_name(knobs[1])+'.png')
            kick_title = str('Double Kick at '+knob_to_name(knobs[0])+' and '+knob_to_name(knobs[1]))
            cpymad_plot_CO_kicks_from_file(madx, matched_twiss, kick_title, kick_savename, kicker_list, measurement_file, horizontal);  
            
        # increment label index
        i += 1
        
        print('End matching using element: '+knob_to_name(knobs[0])+' and ' +knob_to_name(knobs[1]))
    
    if horizontal: 
        zipped = list(zip(odf_s1, odf_s2, odf_name1, odf_name2, odf_knob1, odf_knob2, odf_x1, odf_x2, odf_kick1, odf_kick2, odf_std))    
        df = pnd.DataFrame(zipped, columns=['s1','s2', 'name1','name2', 'knob1','knob2', 'x1 (mm)','x2 (mm)', 'kick1 (mrad)','kick2 (mrad)', 'std (mm)'])
    else:
        zipped = list(zip(odf_s1, odf_s2, odf_name1, odf_name2, odf_knob1, odf_knob2, odf_y1, odf_y2, odf_kick1, odf_kick2, odf_std))    
        df = pnd.DataFrame(zipped, columns=['s1','s2', 'name1','name2', 'knob1','knob2', 'y1 (mm)','y2 (mm)', 'kick1 (mrad)','kick2 (mrad)', 'std (mm)'])
    return df


from itertools import combinations
def find_double_kick_match_tunes(madx_instance, cpymad_logfile, sequence_name, vary_list, measurement_file, save_dir, kqtf, kqtd, horizontal=True, threshold=2.0):
    
    # For output dataframe
    odf_s1 = []
    odf_s2 = []
    odf_name1=[]
    odf_name2=[]
    odf_knob1=[]
    odf_knob2=[]
    odf_kick1=[]
    odf_kick2=[]
    odf_x1 = [] 
    odf_x2 = [] 
    odf_y1 = [] 
    odf_y2 = [] 
    odf_std=[]
    
    # Make directory if it doesn't exist (function checks)
    make_directory(save_dir)

    # Make constraints from file
    constraints = constraints_from_file(measurement_file, horizontal)
    
    # Make kicker list
    isis_reset_steering(madx_instance)
    madx_instance.globals.kqtf = kqtf
    madx_instance.globals.kqtd = kqtd
    twissfile = str(save_dir+'initial_twiss_.tfs')
    initial_twiss = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, twissfile)
    #kicker_list = kicker_list_from_twissfile(twissfile, horizontal)
    
    i = 0
    itertool_combinations = combinations(vary_list, 2)
    for knobs in list(itertool_combinations):
        # first vary name is knobs[0], second is knobs[1]

        print('\nMatching using elements: '+knob_to_name(knobs[0])+' and ' +knob_to_name(knobs[1]))
        
        # Reset all kickers
        isis_reset_steering(madx_instance)
        madx_instance.globals.kqtf = kqtf
        madx_instance.globals.kqtd = kqtd
        
        # Attempt match
        match_double_kick_from_file_and_constraints(madx_instance, knobs[0], knobs[1], constraints)
        
        # Perform twiss
        twissfile = str(save_dir+'matched_twiss_'+str(i)+'_'+knobs[0].split('_')[0]+'_'+knobs[1].split('_')[0]+'.tfs')
        matched_twiss = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, twissfile)
        matched_twiss = pandas_lowercase_all_strings(matched_twiss)
        
        # Check and plot difference
        match_plot_name = str(save_dir+'cf_co_'+str(i)+'_'+knobs[0].split('_')[0]+'_'+knobs[1].split('_')[0]+'.png')
        match_std = compare_measurement_and_twiss(twissfile, measurement_file, match_plot_name, horizontal, threshold);
        print('\t\tKnobs '+knobs[0].split('_')[0]+' and ' +knobs[1].split('_')[0]+' standard deviation = '+str(match_std))
        
        if match_std >= threshold:
            os.remove(twissfile)
        else:
            # Store data for output
            odf_s1.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].s[0])
            odf_s2.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].s[0])
            odf_name1.append(knob_to_name(knobs[0]))
            odf_name2.append(knob_to_name(knobs[1]))
            odf_knob1.append(knobs[0])
            odf_knob2.append(knobs[1])
            odf_std.append(match_std)
            if horizontal:
                odf_x1.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].x[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].index[0])]*1E3)
                odf_x2.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].x[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].index[0])]*1E3)
                odf_kick1.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].hkick[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].index[0])],3)*1E3)
                odf_kick2.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].hkick[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].index[0])],3)*1E3)
            else:
                odf_y1.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].y[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].index[0])]*1E3)
                odf_y2.append(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].y[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].index[0])]*1E3)
                odf_kick1.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].vkick[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[0]))].index[0])],3)*1E3)
                odf_kick2.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].vkick[(matched_twiss[matched_twiss['name'].str.contains(knob_to_name(knobs[1]))].index[0])],3)*1E3)
            # Plot CO with kick
            kick_savename = str(save_dir+'Kicks_'+str(i)+'_'+knob_to_name(knobs[0])+'_' +knob_to_name(knobs[1])+'.png')
            kick_title = str('Double Kick at '+knob_to_name(knobs[0])+' and '+knob_to_name(knobs[1]))
            cpymad_plot_CO_kicks_from_file(madx, matched_twiss, kick_title, kick_savename, kicker_list, measurement_file, horizontal);  
            
        # increment label index
        i += 1
        
        print('End matching using element: '+knob_to_name(knobs[0])+' and ' +knob_to_name(knobs[1]))
    
    if horizontal: 
        zipped = list(zip(odf_s1, odf_s2, odf_name1, odf_name2, odf_knob1, odf_knob2, odf_x1, odf_x2, odf_kick1, odf_kick2, odf_std))    
        df = pnd.DataFrame(zipped, columns=['s1','s2', 'name1','name2', 'knob1','knob2', 'x1 (mm)','x2 (mm)', 'kick1 (mrad)','kick2 (mrad)', 'std (mm)'])
    else:
        zipped = list(zip(odf_s1, odf_s2, odf_name1, odf_name2, odf_knob1, odf_knob2, odf_y1, odf_y2, odf_kick1, odf_kick2, odf_std))    
        df = pnd.DataFrame(zipped, columns=['s1','s2', 'name1','name2', 'knob1','knob2', 'y1 (mm)','y2 (mm)', 'kick1 (mrad)','kick2 (mrad)', 'std (mm)'])
    return df


def return_best_matched_kick(matching_df):    
    index = matching_df['std (mm)'].idxmin()   
    print('return_best_matched_kick: Kicker '+str(matching_df['name'][index])+' STD = '+str(round_sig(matching_df['std (mm)'][index],3))+' Kick = '+str(matching_df['kick (mrad)'][index])+' mrad')
    return(str(matching_df['knob'][index]), float(matching_df['kick (mrad)'][index]))
  

def return_best_matched_kicks(matching_df):    
    index = matching_df['std (mm)'].idxmin()   
    print('return_best_matched_kick: Kicker STD = '+str(round_sig(matching_df['std (mm)'][index],3)))
    return(matching_df.iloc[index])


def find_second_kick_match(madx_instance, cpymad_logfile, sequence_name, vary_list, measurement_file, save_dir, first_knob, first_kick, horizontal=True, threshold=2.0):
    
    # For output dataframe
    odf_s = []
    odf_x = [] 
    odf_y = [] 
    odf_name=[]
    odf_knob=[]
    odf_kick=[]
    odf_std=[]
    
    # Make directory if it doesn't exist (function checks)
    make_directory(save_dir)

    # Make constraints from file
    constraints = constraints_from_file(measurement_file, 'Monitor', horizontal)
    
    # Make kicker list
    exec_str = 'madx_instance.globals.'+str(first_knob)+' = '+str(first_kick*1E-3) # note kick in mrad
    print(exec_str)
    exec(exec_str)
    twissfile = str(save_dir+'initial_twiss_.tfs')
    initial_twiss = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, twissfile)
    kicker_list = kicker_list_from_twissfile(twissfile, horizontal)
    kicker_list_x = kicker_list_from_twissfile(twissfile, True)
    kicker_list_y = kicker_list_from_twissfile(twissfile, False)
    
    # Iterate over knobs (variables used to vary kick strength)
    i = 0
    for knob in vary_list:
        print('\nMatching using element: '+knob.split('_')[0])
        
        # Reset all kickers
        isis_reset_steering(madx_instance)
        
        # Implement first kick
        exec_str = 'madx_instance.globals.'+str(first_knob)+' = '+str(first_kick*1E-3) # note kick in mrad
        exec(exec_str)
    
        # Attempt match
        match_single_kick_from_file_and_constraints(madx_instance, knob, constraints)
        
        # Perform twiss
        twissfile = str(save_dir+'matched_twiss_'+str(i)+'_'+knob.split('_')[0]+'.tfs')
        matched_twiss = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, twissfile)
        matched_twiss = pandas_lowercase_all_strings(matched_twiss)
        
        # Check and plot difference
        match_std = compare_measurement_and_twiss(twissfile, measurement_file, str(save_dir+'cf_co_'+str(i)+'_'+knob.split('_')[0]+'.png'), horizontal, threshold);
        print('\t\tKnob '+knob.split('_')[0]+' standard deviation = '+str(match_std))
        
        if match_std >= threshold:
            os.remove(twissfile)
        else:
            # Store data for output
            odf_s.append(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].s[0])
            odf_name.append(kicker_list[i])
            odf_knob.append(knob)
            odf_std.append(match_std)
            if horizontal:
                odf_x.append(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].x[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])]*1E3)
                odf_kick.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].hkick[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])],3)*1E3)
            else:
                odf_y.append(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].y[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])]*1E3)
                odf_kick.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].vkick[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])],3)*1E3)
            # Plot CO with kick
            cpymad_plot_CO_kicks_from_file_xy(madx, matched_twiss, ('Second Kick at '+knob.split('_')[0]), str(save_dir+'Kicks_'+str(i)+'_'+knob.split('_')[0]+'.png'), kicker_list_x, kicker_list_y, measurement_file, horizontal);  
            
        # increment label index
        i += 1
        
        print('End matching using elements: '+knob.split('_')[0])
    
    if horizontal: 
        zipped = list(zip(odf_s, odf_name, odf_knob, odf_x, odf_kick, odf_std))    
        df = pnd.DataFrame(zipped, columns=['s', 'name', 'knob', 'x (mm)', 'kick (mrad)', 'std (mm)'])
    else:
        zipped = list(zip(odf_s, odf_name, odf_knob, odf_y, odf_kick, odf_std))    
        df = pnd.DataFrame(zipped, columns=['s', 'name', 'knob', 'y (mm)', 'kick (mrad)', 'std (mm)'])
    return df


def find_single_kick_match_tunes(madx_instance, cpymad_logfile, sequence_name, vary_list, measurement_file, save_dir, kqtf, kqtd, horizontal=True, threshold=2.0):
    
    # For output dataframe
    odf_s = []
    odf_x = [] 
    odf_y = [] 
    odf_name=[]
    odf_knob=[]
    odf_kick=[]
    odf_std=[]
    
    # Make directory if it doesn't exist (function checks)
    make_directory(save_dir)

    # Make constraints from file
    constraints = constraints_from_file(measurement_file, horizontal)
    
    # Make kicker list
    isis_reset_steering(madx_instance)
    madx_instance.globals.kqtf = kqtf
    madx_instance.globals.kqtd = kqtd
    twissfile = str(save_dir+'initial_twiss_.tfs')
    initial_twiss = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, twissfile)
    kicker_list = kicker_list_from_twissfile(twissfile, horizontal)
    
    # Iterate over knobs (variables used to vary kick strength)
    i = 0
    for knob in vary_list:
        print('\nMatching using element: '+knob.split('_')[0])
        
        # Reset all kickers
        isis_reset_steering(madx_instance)
        madx_instance.globals.kqtf = kqtf
        madx_instance.globals.kqtd = kqtd
        
        # Attempt match
        match_single_kick_from_file_and_constraints(madx_instance, knob, constraints)
        
        # Perform twiss
        twissfile = str(save_dir+'matched_twiss_'+str(i)+'_'+knob.split('_')[0]+'.tfs')
        matched_twiss = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, twissfile)
        matched_twiss = pandas_lowercase_all_strings(matched_twiss)
        
        # Check and plot difference
        match_savename = str(save_dir+'cf_co_'+str(i)+'_'+knob.split('_')[0]+'.png')
        match_std = compare_measurement_and_twiss(twissfile, measurement_file, match_savename, horizontal, threshold);
        print('\t\tKnob '+knob.split('_')[0]+' standard deviation = '+str(match_std))
        
        if match_std >= threshold:
            os.remove(twissfile)
        else:
            # Store data for output
            odf_s.append(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].s[0])
            odf_name.append(kicker_list[i])
            odf_knob.append(knob)
            odf_std.append(match_std)
            if horizontal:
                odf_x.append(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].x[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])]*1E3)
                odf_kick.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].hkick[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])],3)*1E3)
            else:
                odf_y.append(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].y[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])]*1E3)
                odf_kick.append(round_sig(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].vkick[(matched_twiss[matched_twiss['name'].str.contains(kicker_list[i])].index[0])],3)*1E3)

            # Plot CO with kick
            kick_savename = str(save_dir+'Kicks_'+str(i)+'_'+knob.split('_')[0]+'.png')
            cpymad_plot_CO_kicks_from_file(madx, matched_twiss, ('Single Kick at '+knob.split('_')[0]), kick_savename, kicker_list, measurement_file, horizontal);  
            
        # increment label index
        i += 1
        
        print('End matching using element: '+knob.split('_')[0])
    
    if horizontal: 
        zipped = list(zip(odf_s, odf_name, odf_knob, odf_x, odf_kick, odf_std))    
        df = pnd.DataFrame(zipped, columns=['s', 'name', 'knob', 'x (mm)', 'kick (mrad)', 'std (mm)'])
    else:
        zipped = list(zip(odf_s, odf_name, odf_knob, odf_y, odf_kick, odf_std))    
        df = pnd.DataFrame(zipped, columns=['s', 'name', 'knob', 'y (mm)', 'kick (mrad)', 'std (mm)'])
    return df
    
def match_tune_ptc(madx_instance, sequence_name, requested_q1, requested_q2): #, beta_limits=None):

    ptc_command = '''  
    ptc_twiss_macro(order): macro = {
        select, flag=ptc_twiss, column=name,keyword,s,l,mu1,mu2,beta11,alfa11,beta22,alfa22,x,px,y,py,t,pt,disp1,disp2,disp3,disp4,energy,angle,K0L,K0SL,K1L,K1SL,K2L,K2SL,K3L,K3SL,K4L,K4SL,K5L,K5SL,vkick,hkick,tilt,slot_id,volt,lag,freq,harmon,gamma11,gamma22;
        ptc_create_universe;
        ptc_create_layout, time=false, model=2, exact=true, method=6, nst=3;
        ptc_twiss, closed_orbit, icase=56, no=order, table=ptc_twiss, summary_table=ptc_twiss_summary;  
        ptc_end;
    };

    !match_Tunes(QQx, QQy, sequence_name, betx_inj, bety_inj): macro={
    match_Tunes(QQx, QQy, sequence_name): macro={
        match, use_macro;
            vary, name = kqtd, step=1.0E-4;
            vary, name = kqtf, step=1.0E-4;
            use_macro, name = ptc_twiss_macro(2);
            !constraint, sequence=sequence_name, range=minj1_0, betx<betx_inj, bety<bety_inj, alfx>0, alfy>0;
            !constraint, sequence=sequence_name, range=minj2_0, betx<betx_inj, bety<bety_inj, alfx<0, alfy<0;
            !constraint, sequence=sequence_name, range=m1_0/m3_0, dx<1E-3;             
            constraint, expr = table(ptc_twiss_summary,q1)  = QQx;
            constraint, expr = table(ptc_twiss_summary,q2)  = QQy;
            jacobian,calls=50000,bisec=3;
        endmatch;
        value, kqtd, kqtf;
    };
    '''
    madx_instance.input(ptc_command)
    #ptc_command = 'exec, match_Tunes('+str(requested_q1)+', '+str(requested_q2)+', '+str(sequence_name)+', '+str(beta_limits[0])+', '+str(beta_limits[1])+');'
    ptc_command = 'exec, match_Tunes('+str(requested_q1)+', '+str(requested_q2)+', '+str(sequence_name)+');'
    madx_instance.input(ptc_command)
    
    ptc_twiss = tfs.read(file_out)
    pandas_dataframe_lowercase_columns(ptc_twiss)
    
    return ptc_twiss

def match_tune(madx_instance, sequence, requested_q1, requested_q2, requested_dq1=None, requested_dq2=None):
    
    madx_instance.command.match(chrom=True)

    madx_instance.command.vary(name='kqtd', step=1E-4)
    madx_instance.command.vary(name='kqtf', step=1E-4)

    if requested_dq1 or requested_dq2 == None:
        madx_instance.command.global_(sequence=sequence, q1=requested_q1, q2=requested_q2)
    else:
        madx_instance.command.global_(sequence=sequence, q1=requested_q1, q2=requested_q2, dq1=requested_dq1, dq2=requested_dq2)

    # Example constraints using markers m and minj
    #constraints=[ dict(range='m1_0/m3_0', dx=Constraint(max=0.)),
    #    dict(range='minj1_0', betx=Constraint(min=0.,max=beta_limits[0]), bety=Constraint(min=0.,max=beta_limits[1])),
    #    dict(range='minj2_0', betx=Constraint(min=0.,max=beta_limits[0]), bety=Constraint(min=0.,max=beta_limits[1]))]

    #for c in constraints:
            #madx_instance.command.constraint(**c)
        
    madx_instance.command.jacobian(calls=50000, tolerance=1e-6)
    madx_instance.command.endmatch()
    
########################################################################
########################################################################
#                            PLOTTING
########################################################################
########################################################################

########################################################################
# Plot H/V kicker values to axis as vertical bar plot with closed orbit from tfs file
########################################################################
def kicker_co_plot(ax1, tfs_file, horizontal=True, start_element=None, limits=None, ptc_twiss=False):
  
    # Read TFS file to dataframe using tfs
    df_myTwiss = tfs.read(tfs_file)
    # Make columns lowercase
    df_myTwiss.columns = map(str.lower, df_myTwiss.columns)
    # Make values lowercase
    for key in df_myTwiss.keys():
        if df_myTwiss[str(key)].dtype == 'O':
            df_myTwiss[str(key)] = df_myTwiss[str(key)].str.lower()

    # dataframe of steering magnet kickers
    hd_df = df_myTwiss[df_myTwiss['name'].str.contains('hd')]
    vd_df = df_myTwiss[df_myTwiss['name'].str.contains('vd')]

    # dataframe of main dipole start kickers
    mds_df = df_myTwiss[df_myTwiss['name'].str.contains('mdsk')]
    # dataframe of main dipole end kickers
    mde_df = df_myTwiss[df_myTwiss['name'].str.contains('mdek')]

    # dataframe of trim quad kickers
    qtf_df = df_myTwiss[df_myTwiss['name'].str.contains('qtfk')]
    qtd_df = df_myTwiss[df_myTwiss['name'].str.contains('qtdk')]

    # dataframe of singlet quad kickers
    qds_df = df_myTwiss[df_myTwiss['name'].str.contains('qdsk')]
            
    if horizontal:

        all_kicks_df = pnd.concat([hd_df.hkick, vd_df.hkick, mds_df.hkick, mde_df.hkick, qtf_df.hkick, qtd_df.hkick, qds_df.hkick])

        max_kick = np.max(all_kicks_df)*1E3
        min_kick = np.min(all_kicks_df)*1E3
        
        # Plot steering kicks
        ax1.bar(hd_df.s, hd_df.hkick*1E3, 1, color='g', alpha=0.5)
        
        # Plot main dipole start/end kicks
        ax1.bar(mds_df.s, mds_df.hkick*1E3, 1, color='magenta', alpha=0.5)
        ax1.bar(mde_df.s, mde_df.hkick*1E3, 1, color='magenta', alpha=0.5)
        
        # Plot trim quadrupole kicks
        ax1.bar(qtf_df.s, qtf_df.hkick*1E3, 1, color='orange', alpha=0.5)
        ax1.bar(qtd_df.s, qtd_df.hkick*1E3, 1, color='orange', alpha=0.5)
        
        # Plot singlet quadrupole kicks
        ax1.bar(qds_df.s, qds_df.hkick*1E3, 1, color='b', alpha=0.5)

        # Add text to plot
        integrated_kick = round_sig((np.sum(np.abs(hd_df.hkick)) + np.sum(np.abs(mds_df.hkick)) + np.sum(np.abs(mde_df.hkick)) + np.sum(np.abs(qtf_df.hkick)) + np.sum(np.abs(qtd_df.hkick)) + np.sum(np.abs(qds_df.hkick))),3) *1E3
        
        # by superperiod:
        sp0_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp0')].hkick)) ),3)*1E3
        sp1_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp1')].hkick)) ),3)*1E3
        sp2_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp2')].hkick)) ),3)*1E3
        sp3_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp3')].hkick)) ),3)*1E3
        sp4_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp4')].hkick)) ),3)*1E3
        sp5_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp5')].hkick)) ),3)*1E3
        sp6_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp6')].hkick)) ),3)*1E3
        sp7_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp7')].hkick)) ),3)*1E3
        sp8_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp8')].hkick)) ),3)*1E3
        sp9_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp9')].hkick)) ),3)*1E3

        #data = pnd.DataFrame({'0': sp0_kick, '1': sp1_kick, '2': sp2_kick, '3': sp3_kick, '4': sp4_kick, '5': sp5_kick, '6': sp6_kick, '7': sp7_kick, '8': sp8_kick, '9': sp9_kick}, index=[0])
        #widths = np.ones(10)*0.08
        #ax1.table(cellText=data.values, colLabels=data.columns, loc=3, colWidths=widths, edges='open', fontsize=12)
        
        #text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str(sp0_kick) + ', 1: ' + str(sp1_kick) + ', 2: ' + str(sp2_kick) + ', 3: ' + str(sp3_kick) + ', 4: ' + str(sp4_kick) + ', 5: ' + str(sp5_kick) + ', 6: ' + str(sp6_kick) + ', 7: ' + str(sp7_kick) + ', 8: ' + str(sp8_kick) + ', 9: ' + str(round_sig(sp9_kick,3))
        text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str("{:.2f}".format(sp0_kick)) + ', 1: ' + str("{:.2f}".format(sp1_kick)) + ', 2: ' + str("{:.2f}".format(sp2_kick)) + ', 3: ' + str("{:.2f}".format(sp3_kick)) + ', 4: ' + str("{:.2f}".format(sp4_kick)) + ', 5: ' + str("{:.2f}".format(sp5_kick)) + ', 6: ' + str("{:.2f}".format(sp6_kick)) + ', 7: ' + str("{:.2f}".format(sp7_kick)) + ', 8: ' + str("{:.2f}".format(sp8_kick)) + ', 9: ' + str("{:.2f}".format(sp9_kick))
        #text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad'
        ax1.text(0,min_kick-0.5,text_string, fontsize=8);
        
        ax1_1 = ax1.twinx()
        ax1_1.plot(df_myTwiss.s, df_myTwiss.x, lw=0.5, color='k')
        
    # Vertical
    else:
        
        all_kicks_df = pnd.concat([hd_df.vkick, vd_df.vkick, mds_df.vkick, mde_df.vkick, qtf_df.vkick, qtd_df.vkick, qds_df.vkick])

        max_kick = np.max(all_kicks_df)*1E3
        min_kick = np.min(all_kicks_df)*1E3
        
        # Plot steering kicks
        ax1.bar(vd_df.s, vd_df.vkick*1E3, 1, color='g', alpha=0.5)
        
        # Plot main dipole start/end kicks
        ax1.bar(mds_df.s, mds_df.vkick*1E3, 1, color='magenta', alpha=0.5)
        ax1.bar(mde_df.s, mde_df.vkick*1E3, 1, color='magenta', alpha=0.5)
        
        # Plot trim quadrupole kicks
        ax1.bar(qtf_df.s, qtf_df.vkick*1E3, 1, color='orange', alpha=0.5)
        ax1.bar(qtd_df.s, qtd_df.vkick*1E3, 1, color='orange', alpha=0.5)
        
        # Plot singlet quadrupole kicks
        ax1.bar(qds_df.s, qds_df.vkick*1E3, 1, color='b', alpha=0.5)

        # Add text to plot
        integrated_kick = round_sig((np.sum(np.abs(hd_df.vkick)) + np.sum(np.abs(mds_df.vkick)) + np.sum(np.abs(mde_df.vkick)) + np.sum(np.abs(qtf_df.vkick)) + np.sum(np.abs(qtd_df.vkick)) + np.sum(np.abs(qds_df.vkick))),3) *1E3
        
        # by superperiod:
        sp0_kick = round_sig(( np.sum(np.abs(hd_df[vd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp0')].vkick)) ),3)*1E3
        sp1_kick = round_sig(( np.sum(np.abs(hd_df[vd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp1')].vkick)) ),3)*1E3
        sp2_kick = round_sig(( np.sum(np.abs(hd_df[vd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp2')].vkick)) ),3)*1E3
        sp3_kick = round_sig(( np.sum(np.abs(hd_df[vd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp3')].vkick)) ),3)*1E3
        sp4_kick = round_sig(( np.sum(np.abs(hd_df[vd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp4')].vkick)) ),3)*1E3
        sp5_kick = round_sig(( np.sum(np.abs(hd_df[vd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp5')].vkick)) ),3)*1E3
        sp6_kick = round_sig(( np.sum(np.abs(hd_df[vd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp6')].vkick)) ),3)*1E3
        sp7_kick = round_sig(( np.sum(np.abs(hd_df[vd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp7')].vkick)) ),3)*1E3
        sp8_kick = round_sig(( np.sum(np.abs(hd_df[vd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp8')].vkick)) ),3)*1E3
        sp9_kick = round_sig(( np.sum(np.abs(hd_df[vd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp9')].vkick)) ),3)*1E3

        #text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str(sp0_kick) + ', 1: ' + str(sp1_kick) + ', 2: ' + str(sp2_kick) + ', 3: ' + str(sp3_kick) + ', 4: ' + str(sp4_kick) + ', 5: ' + str(sp5_kick) + ', 6: ' + str(sp6_kick) + ', 7: ' + str(sp7_kick) + ', 8: ' + str(sp8_kick) + ', 9: ' + str(round_sig(sp9_kick,3))
        text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str("{:.2f}".format(sp0_kick)) + ', 1: ' + str("{:.2f}".format(sp1_kick)) + ', 2: ' + str("{:.2f}".format(sp2_kick)) + ', 3: ' + str("{:.2f}".format(sp3_kick)) + ', 4: ' + str("{:.2f}".format(sp4_kick)) + ', 5: ' + str("{:.2f}".format(sp5_kick)) + ', 6: ' + str("{:.2f}".format(sp6_kick)) + ', 7: ' + str("{:.2f}".format(sp7_kick)) + ', 8: ' + str("{:.2f}".format(sp8_kick)) + ', 9: ' + str("{:.2f}".format(sp9_kick))
        ax1.text(0,min_kick-0.5,text_string, fontsize=6);
    
        custom_lines = [Line2D([0], [0], color='g', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='b', lw=4), Line2D([0], [0], color='magenta', lw=4)]
        ax1.legend(custom_lines, ['Steering Magnets', 'Trim Quads', 'Singlet Quad', 'Main Dipole Ends'], loc=1, fontsize=6)        
        
        ax1.grid(which='both', ls=':', lw=0.5, color='grey')
        ax1.set_ylim(min_kick-1, max_kick+1);
        
        ax1.set_ylabel('Kick [mrad]')
        
        ax1_1 = ax1.twinx()
        ax1_1.plot(df_myTwiss.s, df_myTwiss.y, lw=0.5, color='k')
    
    custom_lines = [Line2D([0], [0], color='g', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='b', lw=4), Line2D([0], [0], color='magenta', lw=4)]
    ax1.legend(custom_lines, ['Steering Magnets', 'Trim Quads', 'Singlet Quads', 'Main Dipole Ends'], loc=1, fontsize=6)        

    ax1.grid(which='both', ls=':', lw=0.5, color='grey')
    ax1.set_ylim(min_kick-1, max_kick+1);

    ax1.set_ylabel('Kick [mrad]')

########################################################################
# Plot H/V kicker values to axis as vertical bar plot from tfs file
########################################################################
def kicker_plot(ax_1, tfs_file, horizontal=True, start_element=None, limits=None, ptc_twiss=False):
  
    # Read TFS file to dataframe using tfs
    df_myTwiss = tfs.read(tfs_file)
    # Make columns lowercase
    df_myTwiss.columns = map(str.lower, df_myTwiss.columns)
    # Make values lowercase
    for key in df_myTwiss.keys():
        if df_myTwiss[str(key)].dtype == 'O':
            df_myTwiss[str(key)] = df_myTwiss[str(key)].str.lower()

    # dataframe of steering magnet kickers
    hd_df = df_myTwiss[df_myTwiss['name'].str.contains('hd')]
    vd_df = df_myTwiss[df_myTwiss['name'].str.contains('vd')]

    # dataframe of main dipole start kickers
    mds_df = df_myTwiss[df_myTwiss['name'].str.contains('mdsk')]
    # dataframe of main dipole end kickers
    mde_df = df_myTwiss[df_myTwiss['name'].str.contains('mdek')]

    # dataframe of trim quad kickers
    qtf_df = df_myTwiss[df_myTwiss['name'].str.contains('qtfk')]
    qtd_df = df_myTwiss[df_myTwiss['name'].str.contains('qtdk')]
    
    # dataframe of trim quad kickers
    qf_df = df_myTwiss[df_myTwiss['name'].str.contains('qfk')]
    qd_df = df_myTwiss[df_myTwiss['name'].str.contains('qdk')]
    
    # dataframe of singlet quad kickers
    qds_df = df_myTwiss[df_myTwiss['name'].str.contains('qdsk')]
            
    if horizontal:

        all_kicks_df = pnd.concat([hd_df.hkick, vd_df.hkick, mds_df.hkick, mde_df.hkick, qtf_df.hkick, qtd_df.hkick, qf_df.hkick, qd_df.hkick, qds_df.hkick])

        max_kick = np.max(all_kicks_df)*1E3
        min_kick = np.min(all_kicks_df)*1E3
        
        # Plot steering kicks
        ax_1.bar(hd_df.s, hd_df.hkick*1E3, 1, color='k', alpha=0.5)
        
        # Plot main dipole start/end kicks
        ax_1.bar(mds_df.s, mds_df.hkick*1E3, 1, color='yellow', alpha=0.5)
        ax_1.bar(mde_df.s, mde_df.hkick*1E3, 1, color='yellow', alpha=0.5)
        
        # Plot trim quadrupole kicks
        ax_1.bar(qtf_df.s, qtf_df.hkick*1E3, 1, color='magenta', alpha=0.5)
        ax_1.bar(qtd_df.s, qtd_df.hkick*1E3, 1, color='magenta', alpha=0.5)
        
        # Plot main quadrupole kicks
        ax_1.bar(qf_df.s, qf_df.hkick*1E3, 1, color='orange', alpha=0.5)
        ax_1.bar(qd_df.s, qd_df.hkick*1E3, 1, color='orange', alpha=0.5)
        
        # Plot singlet quadrupole kicks
        ax_1.bar(qds_df.s, qds_df.hkick*1E3, 1, color='green', alpha=0.5)

        # Add text to plot
        integrated_kick = round_sig((np.sum(np.abs(hd_df.hkick)) + np.sum(np.abs(mds_df.hkick)) + np.sum(np.abs(mde_df.hkick)) + np.sum(np.abs(qtf_df.hkick)) + np.sum(np.abs(qtd_df.hkick)) + np.sum(np.abs(qf_df.hkick)) + np.sum(np.abs(qd_df.hkick)) + np.sum(np.abs(qds_df.hkick))),3) *1E3
        
        # by superperiod:
        sp0_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp0')].hkick)) ),3)*1E3
        sp1_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp1')].hkick)) ),3)*1E3
        sp2_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp2')].hkick)) ),3)*1E3
        sp3_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp3')].hkick)) ),3)*1E3
        sp4_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp4')].hkick)) ),3)*1E3
        sp5_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp5')].hkick)) ),3)*1E3
        sp6_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp6')].hkick)) ),3)*1E3
        sp7_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp7')].hkick)) ),3)*1E3
        sp8_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp8')].hkick)) ),3)*1E3
        sp9_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp9')].hkick)) ),3)*1E3

        #data = pnd.DataFrame({'0': sp0_kick, '1': sp1_kick, '2': sp2_kick, '3': sp3_kick, '4': sp4_kick, '5': sp5_kick, '6': sp6_kick, '7': sp7_kick, '8': sp8_kick, '9': sp9_kick}, index=[0])
        #widths = np.ones(10)*0.08
        #ax_1.table(cellText=data.values, colLabels=data.columns, loc=3, colWidths=widths, edges='open', fontsize=12)
        
        #text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str(sp0_kick) + ', 1: ' + str(sp1_kick) + ', 2: ' + str(sp2_kick) + ', 3: ' + str(sp3_kick) + ', 4: ' + str(sp4_kick) + ', 5: ' + str(sp5_kick) + ', 6: ' + str(sp6_kick) + ', 7: ' + str(sp7_kick) + ', 8: ' + str(sp8_kick) + ', 9: ' + str(round_sig(sp9_kick,3))
        text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str("{:.2f}".format(sp0_kick)) + ', 1: ' + str("{:.2f}".format(sp1_kick)) + ', 2: ' + str("{:.2f}".format(sp2_kick)) + ', 3: ' + str("{:.2f}".format(sp3_kick)) + ', 4: ' + str("{:.2f}".format(sp4_kick)) + ', 5: ' + str("{:.2f}".format(sp5_kick)) + ', 6: ' + str("{:.2f}".format(sp6_kick)) + ', 7: ' + str("{:.2f}".format(sp7_kick)) + ', 8: ' + str("{:.2f}".format(sp8_kick)) + ', 9: ' + str("{:.2f}".format(sp9_kick))
        #text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad'
        ax_1.text(0,-0.95,text_string, fontsize=8);
        
    # Vertical
    else:
        
        all_kicks_df = pnd.concat([hd_df.vkick, vd_df.vkick, mds_df.vkick, mde_df.vkick, qtf_df.vkick, qtd_df.vkick, qf_df.vkick, qd_df.vkick, qds_df.vkick])

        max_kick = np.max(all_kicks_df)*1E3
        min_kick = np.min(all_kicks_df)*1E3
        
        # Plot steering kicks
        ax_1.bar(vd_df.s, vd_df.vkick*1E3, 1, color='k', alpha=0.5)
        
        # Plot main dipole start/end kicks
        ax_1.bar(mds_df.s, mds_df.vkick*1E3, 1, color='yellow', alpha=0.5)
        ax_1.bar(mde_df.s, mde_df.vkick*1E3, 1, color='yellow', alpha=0.5)
        
        # Plot trim quadrupole kicks
        ax_1.bar(qtf_df.s, qtf_df.vkick*1E3, 1, color='magenta', alpha=0.5)
        ax_1.bar(qtd_df.s, qtd_df.vkick*1E3, 1, color='magenta', alpha=0.5)
        
        # Plot main quadrupole kicks
        ax_1.bar(qf_df.s, qf_df.vkick*1E3, 1, color='orange', alpha=0.5)
        ax_1.bar(qd_df.s, qd_df.vkick*1E3, 1, color='orange', alpha=0.5)
        
        # Plot singlet quadrupole kicks
        ax_1.bar(qds_df.s, qds_df.vkick*1E3, 1, color='green', alpha=0.5)

        # Add text to plot
        integrated_kick = round_sig((np.sum(np.abs(vd_df.vkick)) + np.sum(np.abs(mds_df.vkick)) + np.sum(np.abs(mde_df.vkick)) + np.sum(np.abs(qtf_df.vkick)) + np.sum(np.abs(qtd_df.vkick)) + np.sum(np.abs(qf_df.vkick)) + np.sum(np.abs(qd_df.vkick)) + np.sum(np.abs(qds_df.vkick))),3) *1E3
        
        # by superperiod:
        sp0_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp0')].vkick)) ),3)*1E3
        sp1_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp1')].vkick)) ),3)*1E3
        sp2_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp2')].vkick)) ),3)*1E3
        sp3_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp3')].vkick)) ),3)*1E3
        sp4_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp4')].vkick)) ),3)*1E3
        sp5_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp5')].vkick)) ),3)*1E3
        sp6_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp6')].vkick)) ),3)*1E3
        sp7_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp7')].vkick)) ),3)*1E3
        sp8_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp8')].vkick)) ),3)*1E3
        sp9_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp9')].vkick)) ),3)*1E3

        #text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str(sp0_kick) + ', 1: ' + str(sp1_kick) + ', 2: ' + str(sp2_kick) + ', 3: ' + str(sp3_kick) + ', 4: ' + str(sp4_kick) + ', 5: ' + str(sp5_kick) + ', 6: ' + str(sp6_kick) + ', 7: ' + str(sp7_kick) + ', 8: ' + str(sp8_kick) + ', 9: ' + str(round_sig(sp9_kick,3))
        text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str("{:.2f}".format(sp0_kick)) + ', 1: ' + str("{:.2f}".format(sp1_kick)) + ', 2: ' + str("{:.2f}".format(sp2_kick)) + ', 3: ' + str("{:.2f}".format(sp3_kick)) + ', 4: ' + str("{:.2f}".format(sp4_kick)) + ', 5: ' + str("{:.2f}".format(sp5_kick)) + ', 6: ' + str("{:.2f}".format(sp6_kick)) + ', 7: ' + str("{:.2f}".format(sp7_kick)) + ', 8: ' + str("{:.2f}".format(sp8_kick)) + ', 9: ' + str("{:.2f}".format(sp9_kick))
        ax_1.text(0,-0.95,text_string, fontsize=8);    

        ax_1.set_ylim(min_kick-1, max_kick+1);
        ax_1.set_ylabel('Kick [mrad]')
    
    custom_lines = [Line2D([0], [0], color='k', lw=4), Line2D([0], [0], color='yellow', lw=4), Line2D([0], [0], color='magenta', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='green', lw=4)]
    ax_1.legend(custom_lines, ['Steering Magnets', 'Main Dipole Ends', 'Trim Quads', 'Main Quads', 'Singlet Quad'], loc=1, fontsize=6)        

    #ax_1.grid(which='both', ls=':', lw=0.5, color='grey')
    ax_1.set_ylim(-1, 1);

    ax_1.set_ylabel('Kick [mrad]')

########################################################################
# Plot H/V kicker difference between two twiss files
########################################################################
def kicker_plot_diff(ax_1, tfs_file_0, tfs_file, horizontal=True, start_element=None, limits=None, ptc_twiss=False):
    
    # tfs_file_0 is the reference
    
    # Read TFS file to dataframe using tfs
    df_myTwiss_0 = tfs.read(tfs_file_0)
    # Make columns lowercase
    df_myTwiss_0.columns = map(str.lower, df_myTwiss_0.columns)
    # Make values lowercase
    for key in df_myTwiss_0.keys():
        if df_myTwiss_0[str(key)].dtype == 'O':
            df_myTwiss_0[str(key)] = df_myTwiss_0[str(key)].str.lower()
            
    # dataframe of steering magnet kickers
    hd_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('hd')]
    vd_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('vd')]

    # dataframe of main dipole start kickers
    mds_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('mdsk')]
    # dataframe of main dipole end kickers
    mde_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('mdek')]

    # dataframe of trim quad kickers
    qtf_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('qtfk')]
    qtd_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('qtdk')]
    
    # dataframe of main quad kickers
    qf_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('qfk')]
    qd_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('qdk')]

    # dataframe of singlet quad kickers
    qds_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('qdsk')]    
    
    
    # Read TFS file to dataframe using tfs
    df_myTwiss = tfs.read(tfs_file)
    # Make columns lowercase
    df_myTwiss.columns = map(str.lower, df_myTwiss.columns)
    # Make values lowercase
    for key in df_myTwiss.keys():
        if df_myTwiss[str(key)].dtype == 'O':
            df_myTwiss[str(key)] = df_myTwiss[str(key)].str.lower()
            
    # dataframe of steering magnet kickers
    hd_df = df_myTwiss[df_myTwiss['name'].str.contains('hd')]
    vd_df = df_myTwiss[df_myTwiss['name'].str.contains('vd')]

    # dataframe of main dipole start kickers
    mds_df = df_myTwiss[df_myTwiss['name'].str.contains('mdsk')]
    mde_df = df_myTwiss[df_myTwiss['name'].str.contains('mdek')]

    # dataframe of trim quad kickers
    qtf_df = df_myTwiss[df_myTwiss['name'].str.contains('qtfk')]
    qtd_df = df_myTwiss[df_myTwiss['name'].str.contains('qtdk')]

    # dataframe of trim quad kickers
    qf_df = df_myTwiss[df_myTwiss['name'].str.contains('qfk')]
    qd_df = df_myTwiss[df_myTwiss['name'].str.contains('qdk')]

    # dataframe of singlet quad kickers
    qds_df = df_myTwiss[df_myTwiss['name'].str.contains('qdsk')]
            
    if horizontal:

        all_kicks_df_0 = pnd.concat([hd_df_0.hkick, vd_df_0.hkick, mds_df_0.hkick, mde_df_0.hkick, qtf_df_0.hkick, qtd_df_0.hkick, qf_df_0.hkick, qd_df_0.hkick, qds_df_0.hkick])
        all_kicks_df = pnd.concat([hd_df.hkick, vd_df.hkick, mds_df.hkick, mde_df.hkick, qtf_df.hkick, qtd_df.hkick, qf_df.hkick, qd_df.hkick, qds_df.hkick])

        max_kick = np.max(all_kicks_df)*1E3
        min_kick = np.min(all_kicks_df)*1E3
        
        # Plot steering kicks
        ax_1.bar(hd_df.s, (hd_df.hkick-hd_df_0.hkick)*1E3, 1, color='k', alpha=0.5)
        
        # Plot main dipole start/end kicks
        ax_1.bar(mds_df.s, (mds_df.hkick-mds_df_0.hkick)*1E3, 1, color='yellow', alpha=0.5)
        ax_1.bar(mde_df.s, (mde_df.hkick-mde_df_0.hkick)*1E3, 1, color='yellow', alpha=0.5)
        
        # Plot trim quadrupole kicks
        ax_1.bar(qtf_df.s, (qtf_df.hkick-qtf_df_0.hkick)*1E3, 1, color='magenta', alpha=0.5)
        ax_1.bar(qtd_df.s, (qtd_df.hkick-qtd_df_0.hkick)*1E3, 1, color='magenta', alpha=0.5)
        
        # Plot main quadrupole kicks
        ax_1.bar(qf_df.s, (qf_df.hkick-qf_df_0.hkick)*1E3, 1, color='orange', alpha=0.5)
        ax_1.bar(qd_df.s, (qd_df.hkick-qd_df_0.hkick)*1E3, 1, color='orange', alpha=0.5)
                
        # Plot singlet quadrupole kicks
        ax_1.bar(qds_df.s, (qds_df.hkick-qds_df_0.hkick)*1E3, 1, color='green', alpha=0.5)

        # Add text to plot
        integrated_kick_0 = round_sig((np.sum(np.abs(hd_df_0.hkick)) + np.sum(np.abs(mds_df_0.hkick)) + np.sum(np.abs(mde_df_0.hkick)) + np.sum(np.abs(qtf_df_0.hkick)) + np.sum(np.abs(qtd_df_0.hkick)) + np.sum(np.abs(qf_df_0.hkick)) + np.sum(np.abs(qd_df_0.hkick)) + np.sum(np.abs(qds_df_0.hkick))),3) *1E3
        
        # by superperiod:
        sp0_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp0')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp0')].hkick)) ),3)*1E3
        sp1_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp1')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp1')].hkick)) ),3)*1E3
        sp2_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp2')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp2')].hkick)) ),3)*1E3
        sp3_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp3')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp3')].hkick)) ),3)*1E3
        sp4_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp4')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp4')].hkick)) ),3)*1E3
        sp5_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp5')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp5')].hkick)) ),3)*1E3
        sp6_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp6')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp6')].hkick)) ),3)*1E3
        sp7_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp7')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp7')].hkick)) ),3)*1E3
        sp8_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp8')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp8')].hkick)) ),3)*1E3
        sp9_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp9')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp9')].hkick)) ),3)*1E3

        integrated_kick = round_sig((np.sum(np.abs(hd_df.hkick)) + np.sum(np.abs(mds_df.hkick)) + np.sum(np.abs(mde_df.hkick)) + np.sum(np.abs(qtf_df.hkick)) + np.sum(np.abs(qtd_df.hkick)) + np.sum(np.abs(qf_df.hkick)) + np.sum(np.abs(qd_df.hkick)) + np.sum(np.abs(qds_df.hkick))),3) *1E3
        
        # by superperiod:
        sp0_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp0')].hkick)) ),3)*1E3
        sp1_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp1')].hkick)) ),3)*1E3
        sp2_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp2')].hkick)) ),3)*1E3
        sp3_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp3')].hkick)) ),3)*1E3
        sp4_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp4')].hkick)) ),3)*1E3
        sp5_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp5')].hkick)) ),3)*1E3
        sp6_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp6')].hkick)) ),3)*1E3
        sp7_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp7')].hkick)) ),3)*1E3
        sp8_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp8')].hkick)) ),3)*1E3
        sp9_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp9')].hkick)) ),3)*1E3

        
        #text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str("{:.2f}".format(sp0_kick)) + ', 1: ' + str("{:.2f}".format(sp1_kick)) + ', 2: ' + str("{:.2f}".format(sp2_kick)) + ', 3: ' + str("{:.2f}".format(sp3_kick)) + ', 4: ' + str("{:.2f}".format(sp4_kick)) + ', 5: ' + str("{:.2f}".format(sp5_kick)) + ', 6: ' + str("{:.2f}".format(sp6_kick)) + ', 7: ' + str("{:.2f}".format(sp7_kick)) + ', 8: ' + str("{:.2f}".format(sp8_kick)) + ', 9: ' + str("{:.2f}".format(sp9_kick))
        text_string = 'Difference Integrated Kick: ' + str(integrated_kick-integrated_kick_0)  + ' mrad\nBy SP: 0: ' + str("{:.2f}".format(sp0_kick-sp0_kick_0)) + ', 1: ' + str("{:.2f}".format(sp1_kick-sp1_kick_0)) + ', 2: ' + str("{:.2f}".format(sp2_kick-sp2_kick_0)) + ', 3: ' + str("{:.2f}".format(sp3_kick-sp3_kick_0)) + ', 4: ' + str("{:.2f}".format(sp4_kick-sp4_kick_0)) + ', 5: ' + str("{:.2f}".format(sp5_kick-sp5_kick_0)) + ', 6: ' + str("{:.2f}".format(sp6_kick-sp6_kick_0)) + ', 7: ' + str("{:.2f}".format(sp7_kick-sp7_kick_0)) + ', 8: ' + str("{:.2f}".format(sp8_kick-sp8_kick_0)) + ', 9: ' + str("{:.2f}".format(sp9_kick-sp9_kick_0))
        ax_1.text(0,min_kick-0.5,text_string, fontsize=8);
        
    # Vertical
    else:
        
        all_kicks_df = pnd.concat([hd_df.vkick, vd_df.vkick, mds_df.vkick, mde_df.vkick, qtf_df.vkick, qtd_df.vkick, qds_df.vkick])

        max_kick = np.max(all_kicks_df)*1E3
        min_kick = np.min(all_kicks_df)*1E3
        
        # Plot steering kicks
        ax_1.bar(vd_df.s, (vd_df.vkick-vd_df_0.vkick)*1E3, 1, color='k', alpha=0.5)
        
        # Plot main dipole start/end kicks
        ax_1.bar(mds_df.s, (mds_df.vkick-mds_df_0.vkick)*1E3, 1, color='yellow', alpha=0.5)
        ax_1.bar(mde_df.s, (mde_df.vkick-mde_df_0.vkick)*1E3, 1, color='yellow', alpha=0.5)
        
        # Plot trim quadrupole kicks
        ax_1.bar(qtf_df.s, (qtf_df.vkick-qtf_df_0.vkick)*1E3, 1, color='magenta', alpha=0.5)
        ax_1.bar(qtd_df.s, (qtd_df.vkick-qtd_df_0.vkick)*1E3, 1, color='magenta', alpha=0.5)
        
        # Plot main quadrupole kicks
        ax_1.bar(qf_df.s, (qf_df.vkick-qf_df_0.vkick)*1E3, 1, color='orange', alpha=0.5)
        ax_1.bar(qd_df.s, (qd_df.vkick-qd_df_0.vkick)*1E3, 1, color='orange', alpha=0.5)
                
        # Plot singlet quadrupole kicks
        ax_1.bar(qds_df.s, (qds_df.vkick-qds_df_0.vkick)*1E3, 1, color='green', alpha=0.5)
        
        # Add text to plot
        integrated_kick_0 = round_sig((np.sum(np.abs(vd_df_0.vkick)) + np.sum(np.abs(mds_df_0.vkick)) + np.sum(np.abs(mde_df_0.vkick)) + np.sum(np.abs(qtf_df_0.vkick)) + np.sum(np.abs(qtd_df_0.vkick)) + np.sum(np.abs(qf_df_0.vkick)) + np.sum(np.abs(qd_df_0.vkick)) + np.sum(np.abs(qds_df_0.vkick))),3) *1E3
        
        # by superperiod:
        sp0_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp0')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp0')].vkick)) ),3)*1E3
        sp1_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp1')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp1')].vkick)) ),3)*1E3
        sp2_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp2')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp2')].vkick)) ),3)*1E3
        sp3_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp3')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp3')].vkick)) ),3)*1E3
        sp4_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp4')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp4')].vkick)) ),3)*1E3
        sp5_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp5')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp5')].vkick)) ),3)*1E3
        sp6_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp6')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp6')].vkick)) ),3)*1E3
        sp7_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp7')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp7')].vkick)) ),3)*1E3
        sp8_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp8')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp8')].vkick)) ),3)*1E3
        sp9_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp9')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp9')].vkick)) ),3)*1E3

        # Add text to plot
        integrated_kick = round_sig((np.sum(np.abs(vd_df.vkick)) + np.sum(np.abs(mds_df.vkick)) + np.sum(np.abs(mde_df.vkick)) + np.sum(np.abs(qtf_df.vkick)) + np.sum(np.abs(qtd_df.vkick)) + np.sum(np.abs(qf_df.vkick)) + np.sum(np.abs(qd_df.vkick)) + np.sum(np.abs(qds_df.vkick))),3) *1E3
        
        # by superperiod:
        sp0_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp0')].vkick)) ),3)*1E3
        sp1_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp1')].vkick)) ),3)*1E3
        sp2_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp2')].vkick)) ),3)*1E3
        sp3_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp3')].vkick)) ),3)*1E3
        sp4_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp4')].vkick)) ),3)*1E3
        sp5_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp5')].vkick)) ),3)*1E3
        sp6_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp6')].vkick)) ),3)*1E3
        sp7_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp7')].vkick)) ),3)*1E3
        sp8_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp8')].vkick)) ),3)*1E3
        sp9_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp9')].vkick)) ),3)*1E3

        #text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str("{:.2f}".format(sp0_kick)) + ', 1: ' + str("{:.2f}".format(sp1_kick)) + ', 2: ' + str("{:.2f}".format(sp2_kick)) + ', 3: ' + str("{:.2f}".format(sp3_kick)) + ', 4: ' + str("{:.2f}".format(sp4_kick)) + ', 5: ' + str("{:.2f}".format(sp5_kick)) + ', 6: ' + str("{:.2f}".format(sp6_kick)) + ', 7: ' + str("{:.2f}".format(sp7_kick)) + ', 8: ' + str("{:.2f}".format(sp8_kick)) + ', 9: ' + str("{:.2f}".format(sp9_kick))
        text_string = 'Difference Integrated Kick: ' + str(integrated_kick-integrated_kick_0)  + ' mrad\nBy SP: 0: ' + str("{:.2f}".format(sp0_kick-sp0_kick_0)) + ', 1: ' + str("{:.2f}".format(sp1_kick-sp1_kick_0)) + ', 2: ' + str("{:.2f}".format(sp2_kick-sp2_kick_0)) + ', 3: ' + str("{:.2f}".format(sp3_kick-sp3_kick_0)) + ', 4: ' + str("{:.2f}".format(sp4_kick-sp4_kick_0)) + ', 5: ' + str("{:.2f}".format(sp5_kick-sp5_kick_0)) + ', 6: ' + str("{:.2f}".format(sp6_kick-sp6_kick_0)) + ', 7: ' + str("{:.2f}".format(sp7_kick-sp7_kick_0)) + ', 8: ' + str("{:.2f}".format(sp8_kick-sp8_kick_0)) + ', 9: ' + str("{:.2f}".format(sp9_kick-sp9_kick_0))
        ax_1.text(0,min_kick-0.5,text_string, fontsize=8);
         
        #ax_1.grid(which='both', ls=':', lw=0.5, color='grey')
        ax_1.set_ylim(min_kick-1, max_kick+1);
        
        ax_1.set_ylabel('Kick [mrad]')
    
    custom_lines = [Line2D([0], [0], color='k', lw=4), Line2D([0], [0], color='yellow', lw=4), Line2D([0], [0], color='magenta', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='green', lw=4)]
    ax_1.legend(custom_lines, ['Steering Magnets', 'Main Dipole Ends', 'Trim Quads', 'Main Quads', 'Singlet Quad'], loc=1, fontsize=6)        

    #ax_1.grid(which='both', ls=':', lw=0.5, color='grey')
    ax_1.set_ylim(-1, 1);

    ax_1.set_ylabel('Kick [mrad]')

########################################################################
# Plot H/V kicker difference between two twiss files (old)
########################################################################
def kicker_plot_cf(ax_1, tfs_file_0, tfs_file, horizontal=True, start_element=None, limits=None, ptc_twiss=False):
    
    # tfs_file_0 is the reference
    
    # Read TFS file to dataframe using tfs
    df_myTwiss_0 = tfs.read(tfs_file_0)
    # Make columns lowercase
    df_myTwiss_0.columns = map(str.lower, df_myTwiss_0.columns)
    # Make values lowercase
    for key in df_myTwiss_0.keys():
        if df_myTwiss_0[str(key)].dtype == 'O':
            df_myTwiss_0[str(key)] = df_myTwiss_0[str(key)].str.lower()
            
    # dataframe of steering magnet kickers
    hd_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('hd')]
    vd_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('vd')]

    # dataframe of main dipole start kickers
    mds_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('mdsk')]
    # dataframe of main dipole end kickers
    mde_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('mdek')]

    # dataframe of trim quad kickers
    qtf_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('qtfk')]
    qtd_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('qtdk')]
    
    # dataframe of main quad kickers
    qf_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('qfk')]
    qd_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('qdk')]

    # dataframe of singlet quad kickers
    qds_df_0 = df_myTwiss_0[df_myTwiss_0['name'].str.contains('qdsk')]    
    
    
    # Read TFS file to dataframe using tfs
    df_myTwiss = tfs.read(tfs_file)
    # Make columns lowercase
    df_myTwiss.columns = map(str.lower, df_myTwiss.columns)
    # Make values lowercase
    for key in df_myTwiss.keys():
        if df_myTwiss[str(key)].dtype == 'O':
            df_myTwiss[str(key)] = df_myTwiss[str(key)].str.lower()
            
    # dataframe of steering magnet kickers
    hd_df = df_myTwiss[df_myTwiss['name'].str.contains('hd')]
    vd_df = df_myTwiss[df_myTwiss['name'].str.contains('vd')]

    # dataframe of main dipole start kickers
    mds_df = df_myTwiss[df_myTwiss['name'].str.contains('mdsk')]
    mde_df = df_myTwiss[df_myTwiss['name'].str.contains('mdek')]

    # dataframe of trim quad kickers
    qtf_df = df_myTwiss[df_myTwiss['name'].str.contains('qtfk')]
    qtd_df = df_myTwiss[df_myTwiss['name'].str.contains('qtdk')]

    # dataframe of trim quad kickers
    qf_df = df_myTwiss[df_myTwiss['name'].str.contains('qfk')]
    qd_df = df_myTwiss[df_myTwiss['name'].str.contains('qdk')]

    # dataframe of singlet quad kickers
    qds_df = df_myTwiss[df_myTwiss['name'].str.contains('qdsk')]
            
    if horizontal:

        all_kicks_df = pnd.concat([hd_df.hkick, vd_df.hkick, mds_df.hkick, mde_df.hkick, qtf_df.hkick, qtd_df.hkick, qf_df.hkick, qd_df.hkick, qds_df.hkick])

        max_kick = np.max(all_kicks_df)*1E3
        min_kick = np.min(all_kicks_df)*1E3
        
        # Plot steering kicks
        ax_1.bar(hd_df.s, hd_df.hkick*1E3, 1, color='k', alpha=0.5)
        
        # Plot main dipole start/end kicks
        ax_1.bar(mds_df.s, mds_df.hkick*1E3, 1, color='yellow', alpha=0.5)
        ax_1.bar(mde_df.s, mde_df.hkick*1E3, 1, color='yellow', alpha=0.5)
        
        # Plot trim quadrupole kicks
        ax_1.bar(qtf_df.s, qtf_df.hkick*1E3, 1, color='magenta', alpha=0.5)
        ax_1.bar(qtd_df.s, qtd_df.hkick*1E3, 1, color='magenta', alpha=0.5)
        
        # Plot main quadrupole kicks
        ax_1.bar(qf_df.s, qf_df.hkick*1E3, 1, color='orange', alpha=0.5)
        ax_1.bar(qd_df.s, qd_df.hkick*1E3, 1, color='orange', alpha=0.5)
                
        # Plot singlet quadrupole kicks
        ax_1.bar(qds_df.s, qds_df.hkick*1E3, 1, color='green', alpha=0.5)

        # Add text to plot
        integrated_kick_0 = round_sig((np.sum(np.abs(hd_df_0.hkick)) + np.sum(np.abs(mds_df_0.hkick)) + np.sum(np.abs(mde_df_0.hkick)) + np.sum(np.abs(qtf_df_0.hkick)) + np.sum(np.abs(qtd_df_0.hkick)) + np.sum(np.abs(qf_df_0.hkick)) + np.sum(np.abs(qd_df_0.hkick)) + np.sum(np.abs(qds_df_0.hkick))),3) *1E3
        
        # by superperiod:
        sp0_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp0')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp0')].hkick)) ),3)*1E3
        sp1_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp1')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp1')].hkick)) ),3)*1E3
        sp2_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp2')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp2')].hkick)) ),3)*1E3
        sp3_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp3')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp3')].hkick)) ),3)*1E3
        sp4_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp4')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp4')].hkick)) ),3)*1E3
        sp5_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp5')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp5')].hkick)) ),3)*1E3
        sp6_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp6')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp6')].hkick)) ),3)*1E3
        sp7_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp7')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp7')].hkick)) ),3)*1E3
        sp8_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp8')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp8')].hkick)) ),3)*1E3
        sp9_kick_0 = round_sig(( np.sum(np.abs(hd_df_0[hd_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp9')].hkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp9')].hkick)) ),3)*1E3

        integrated_kick = round_sig((np.sum(np.abs(hd_df.hkick)) + np.sum(np.abs(mds_df.hkick)) + np.sum(np.abs(mde_df.hkick)) + np.sum(np.abs(qtf_df.hkick)) + np.sum(np.abs(qtd_df.hkick)) + np.sum(np.abs(qf_df.hkick)) + np.sum(np.abs(qd_df.hkick)) + np.sum(np.abs(qds_df.hkick))),3) *1E3
        
        # by superperiod:
        sp0_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp0')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp0')].hkick)) ),3)*1E3
        sp1_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp1')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp1')].hkick)) ),3)*1E3
        sp2_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp2')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp2')].hkick)) ),3)*1E3
        sp3_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp3')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp3')].hkick)) ),3)*1E3
        sp4_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp4')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp4')].hkick)) ),3)*1E3
        sp5_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp5')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp5')].hkick)) ),3)*1E3
        sp6_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp6')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp6')].hkick)) ),3)*1E3
        sp7_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp7')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp7')].hkick)) ),3)*1E3
        sp8_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp8')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp8')].hkick)) ),3)*1E3
        sp9_kick = round_sig(( np.sum(np.abs(hd_df[hd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].hkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp9')].hkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp9')].hkick)) ),3)*1E3

        
        text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str("{:.2f}".format(sp0_kick)) + ', 1: ' + str("{:.2f}".format(sp1_kick)) + ', 2: ' + str("{:.2f}".format(sp2_kick)) + ', 3: ' + str("{:.2f}".format(sp3_kick)) + ', 4: ' + str("{:.2f}".format(sp4_kick)) + ', 5: ' + str("{:.2f}".format(sp5_kick)) + ', 6: ' + str("{:.2f}".format(sp6_kick)) + ', 7: ' + str("{:.2f}".format(sp7_kick)) + ', 8: ' + str("{:.2f}".format(sp8_kick)) + ', 9: ' + str("{:.2f}".format(sp9_kick))
        text_string += '\n Difference Integrated Kick: ' + str(integrated_kick-integrated_kick_0)  + ' mrad\nBy SP: 0: ' + str("{:.2f}".format(sp0_kick-sp0_kick_0)) + ', 1: ' + str("{:.2f}".format(sp1_kick-sp1_kick_0)) + ', 2: ' + str("{:.2f}".format(sp2_kick-sp2_kick_0)) + ', 3: ' + str("{:.2f}".format(sp3_kick-sp3_kick_0)) + ', 4: ' + str("{:.2f}".format(sp4_kick-sp4_kick_0)) + ', 5: ' + str("{:.2f}".format(sp5_kick-sp5_kick_0)) + ', 6: ' + str("{:.2f}".format(sp6_kick-sp6_kick_0)) + ', 7: ' + str("{:.2f}".format(sp7_kick-sp7_kick_0)) + ', 8: ' + str("{:.2f}".format(sp8_kick-sp8_kick_0)) + ', 9: ' + str("{:.2f}".format(sp9_kick-sp9_kick_0))
        
        ax_1.text(0,-0.95,text_string, fontsize=6);
        
    # Vertical
    else:
        
        all_kicks_df = pnd.concat([hd_df.vkick, vd_df.vkick, mds_df.vkick, mde_df.vkick, qtf_df.vkick, qtd_df.vkick, qds_df.vkick])

        max_kick = np.max(all_kicks_df)*1E3
        min_kick = np.min(all_kicks_df)*1E3
        
        # Plot steering kicks
        ax_1.bar(vd_df.s, vd_df.vkick*1E3, 1, color='k', alpha=0.5)
        
        # Plot main dipole start/end kicks
        ax_1.bar(mds_df.s, mds_df.vkick*1E3, 1, color='yellow', alpha=0.5)
        ax_1.bar(mde_df.s, mde_df.vkick*1E3, 1, color='yellow', alpha=0.5)
        
        # Plot trim quadrupole kicks
        ax_1.bar(qtf_df.s, qtf_df.vkick*1E3, 1, color='magenta', alpha=0.5)
        ax_1.bar(qtd_df.s, qtd_df.vkick*1E3, 1, color='magenta', alpha=0.5)
        
        # Plot main quadrupole kicks
        ax_1.bar(qf_df.s, qf_df.vkick*1E3, 1, color='orange', alpha=0.5)
        ax_1.bar(qd_df.s, qd_df.vkick*1E3, 1, color='orange', alpha=0.5)
        
        # Plot singlet quadrupole kicks
        ax_1.bar(qds_df.s, qds_df.vkick*1E3, 1, color='green', alpha=0.5)
        
        # Add text to plot
        integrated_kick_0 = round_sig((np.sum(np.abs(vd_df_0.vkick)) + np.sum(np.abs(mds_df_0.vkick)) + np.sum(np.abs(mde_df_0.vkick)) + np.sum(np.abs(qtf_df_0.vkick)) + np.sum(np.abs(qtd_df_0.vkick)) + np.sum(np.abs(qf_df_0.vkick)) + np.sum(np.abs(qd_df_0.vkick)) + np.sum(np.abs(qds_df_0.vkick))),3) *1E3
        
        # by superperiod:
        sp0_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp0')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp0')].vkick)) ),3)*1E3
        sp1_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp1')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp1')].vkick)) ),3)*1E3
        sp2_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp2')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp2')].vkick)) ),3)*1E3
        sp3_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp3')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp3')].vkick)) ),3)*1E3
        sp4_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp4')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp4')].vkick)) ),3)*1E3
        sp5_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp5')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp5')].vkick)) ),3)*1E3
        sp6_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp6')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp6')].vkick)) ),3)*1E3
        sp7_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp7')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp7')].vkick)) ),3)*1E3
        sp8_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp8')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp8')].vkick)) ),3)*1E3
        sp9_kick_0 = round_sig(( np.sum(np.abs(vd_df_0[vd_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp9')].vkick)) +  np.sum(np.abs(mde_df_0[mde_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtf_df_0[qtf_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtd_df_0[qtd_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qf_df_0[qf_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qd_df_0[qd_df_0['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qds_df_0[qds_df_0['name'].str.contains('sp9')].vkick)) ),3)*1E3

        # Add text to plot
        integrated_kick = round_sig((np.sum(np.abs(vd_df.vkick)) + np.sum(np.abs(mds_df.vkick)) + np.sum(np.abs(mde_df.vkick)) + np.sum(np.abs(qtf_df.vkick)) + np.sum(np.abs(qtd_df.vkick)) + np.sum(np.abs(qf_df.vkick)) + np.sum(np.abs(qd_df.vkick)) + np.sum(np.abs(qds_df.vkick))),3) *1E3
        
        # by superperiod:
        sp0_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp0')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp0')].vkick)) ),3)*1E3
        sp1_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp1')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp1')].vkick)) ),3)*1E3
        sp2_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp2')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp2')].vkick)) ),3)*1E3
        sp3_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp3')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp3')].vkick)) ),3)*1E3
        sp4_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp4')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp4')].vkick)) ),3)*1E3
        sp5_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp5')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp5')].vkick)) ),3)*1E3
        sp6_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp6')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp6')].vkick)) ),3)*1E3
        sp7_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp7')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp7')].vkick)) ),3)*1E3
        sp8_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp8')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp8')].vkick)) ),3)*1E3
        sp9_kick = round_sig(( np.sum(np.abs(vd_df[vd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].vkick)) +  np.sum(np.abs(mde_df[mde_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtf_df[qtf_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qtd_df[qtd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qf_df[qf_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qd_df[qd_df['name'].str.contains('sp9')].vkick)) + np.sum(np.abs(qds_df[qds_df['name'].str.contains('sp9')].vkick)) ),3)*1E3

        text_string = 'Total Integrated Kick: ' + str(integrated_kick) + ' mrad\n' + 'By SP: 0: ' + str("{:.2f}".format(sp0_kick)) + ', 1: ' + str("{:.2f}".format(sp1_kick)) + ', 2: ' + str("{:.2f}".format(sp2_kick)) + ', 3: ' + str("{:.2f}".format(sp3_kick)) + ', 4: ' + str("{:.2f}".format(sp4_kick)) + ', 5: ' + str("{:.2f}".format(sp5_kick)) + ', 6: ' + str("{:.2f}".format(sp6_kick)) + ', 7: ' + str("{:.2f}".format(sp7_kick)) + ', 8: ' + str("{:.2f}".format(sp8_kick)) + ', 9: ' + str("{:.2f}".format(sp9_kick))
        text_string += '\n Difference Integrated Kick: ' + str(integrated_kick-integrated_kick_0)  + ' mrad\nBy SP: 0: ' + str("{:.2f}".format(sp0_kick-sp0_kick_0)) + ', 1: ' + str("{:.2f}".format(sp1_kick-sp1_kick_0)) + ', 2: ' + str("{:.2f}".format(sp2_kick-sp2_kick_0)) + ', 3: ' + str("{:.2f}".format(sp3_kick-sp3_kick_0)) + ', 4: ' + str("{:.2f}".format(sp4_kick-sp4_kick_0)) + ', 5: ' + str("{:.2f}".format(sp5_kick-sp5_kick_0)) + ', 6: ' + str("{:.2f}".format(sp6_kick-sp6_kick_0)) + ', 7: ' + str("{:.2f}".format(sp7_kick-sp7_kick_0)) + ', 8: ' + str("{:.2f}".format(sp8_kick-sp8_kick_0)) + ', 9: ' + str("{:.2f}".format(sp9_kick-sp9_kick_0))
        ax_1.text(0,-0.95,text_string, fontsize=6);
         
        ax_1.grid(which='both', ls=':', lw=0.5, color='grey')
        ax_1.set_ylim(min_kick-1, max_kick+1);
        
        ax_1.set_ylabel('Kick [mrad]')
    
    custom_lines = [Line2D([0], [0], color='k', lw=4), Line2D([0], [0], color='yellow', lw=4), Line2D([0], [0], color='magenta', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='green', lw=4)]
    ax_1.legend(custom_lines, ['Steering Magnets', 'Main Dipole Ends', 'Trim Quads', 'Main Quads', 'Singlet Quad'], loc=1, fontsize=6)        

    #ax_1.grid(which='both', ls=':', lw=0.5, color='grey')
    ax_1.set_ylim(-1, 1);

    ax_1.set_ylabel('Kick [mrad]')


def plot_kick_co_comparison_twissfiles(twissfile1, twissfile2, mfile1, mfile2, savename, horizontal=True, ftype1='BPM', ftype2='BPM', ylims=None, titles=None):
    
    fig1 = plt.figure(facecolor='w', edgecolor='k', figsize=[8, 10])
    gs = fig1.add_gridspec(ncols=1, nrows=3, height_ratios=[1,1,1])
    gs.update(wspace=0.025, hspace=0.)

    ax1 = fig1.add_subplot(gs[0,0])
    ax2 = fig1.add_subplot(gs[1,0], sharex=ax1)
    ax3 = fig1.add_subplot(gs[2,0], sharex=ax1)

    kicker_plot(ax1, twissfile1, horizontal)
    kicker_plot(ax2, twissfile2, horizontal)    
    kicker_plot_diff(ax3, twissfile1, twissfile2, horizontal)

    if ylims is None:
        ylims=[-20,20]
    else:
        ylims=ylims
        
    ax1_1 = ax1.twinx()
    co_meas_plot(ax1_1, twissfile1, mfile1, horizontal, ylims)

    ax2_1 = ax2.twinx()
    co_meas_plot(ax2_1, twissfile2, mfile2, horizontal, ylims)

    ax3_1 = ax3.twinx()
    if ylims is None:
        co_plot_diff(ax3_1, twissfile1, twissfile2, horizontal=horizontal)
    else:
        co_plot_diff(ax3_1, twissfile1, twissfile2, horizontal=horizontal, ylims=ylims)
    
    ax3.set_xlabel('s [m]');

    if titles != None:
        ax1.set_title(titles[0], y=1.0, pad=-14)
        ax2.set_title(titles[1], y=1.0, pad=-14)
        ax3.set_title(str(titles[0]+' -> '+titles[1]), y=1.0, pad=-14)
    
    plt.savefig(savename, bbox_inches='tight')    


def co_plot_diff(ax_1, tfs_file, tfs_file2, horizontal=True, ylims=None):
  
    # Read TFS file to dataframe using tfs
    df_myTwiss = tfs.read(tfs_file)
    # Make columns lowercase
    df_myTwiss.columns = map(str.lower, df_myTwiss.columns)
    # Make values lowercase
    for key in df_myTwiss.keys():
        if df_myTwiss[str(key)].dtype == 'O':
            df_myTwiss[str(key)] = df_myTwiss[str(key)].str.lower()
            
    # Read TFS file to dataframe using tfs
    df_myTwiss2 = tfs.read(tfs_file2)
    # Make columns lowercase
    df_myTwiss2.columns = map(str.lower, df_myTwiss2.columns)
    # Make values lowercase
    for key in df_myTwiss2.keys():
        if df_myTwiss2[str(key)].dtype == 'O':
            df_myTwiss2[str(key)] = df_myTwiss2[str(key)].str.lower()
            
    if horizontal: 
        ax_1.plot(df_myTwiss.s, (df_myTwiss.x-df_myTwiss2.x)*1E3, lw=0.5, color='k', label='Closed Orbit Difference')
        ax_1.legend(fontsize=6, loc=2)        
        ax_1.set_ylabel('x [mm]')
        
    # Vertical
    else: 
        ax_1.plot(df_myTwiss.s, (df_myTwiss.y-df_myTwiss2.y)*1E3, lw=0.5, color='k', label='Closed Orbit Difference')
        ax_1.legend(fontsize=6, loc=2)        
        ax_1.set_ylabel('y [mm]')

    if ylims:
        ax_1.set_ylim(ylims[0], ylims[1]);
        
    ax_1.grid(which='both', ls=':', lw=0.5, color='grey')
    #ax_1.set_ylim(min_kick-1, max_kick+1);


def co_meas_plot(ax_1, tfs_file, measurement_file, time='0', ylims=None):
  
    meas_df = processed_bpm_dat_to_df(measurement_file)
    #s, meas = measurements_from_file(measurement_file, time=time)
    
    # Read TFS file to dataframe using tfs
    df_myTwiss = tfs.read(tfs_file)
    # Make columns lowercase
    df_myTwiss.columns = map(str.lower, df_myTwiss.columns)
    # Make values lowercase
    for key in df_myTwiss.keys():
        if df_myTwiss[str(key)].dtype == 'O':
            df_myTwiss[str(key)] = df_myTwiss[str(key)].str.lower()
            
    if horizontal: 
        ax_1.plot(df_myTwiss.s, df_myTwiss.x*1E3, lw=0.5, color='k', label='Closed Orbit')
        ax_1.errorbar(meas_df['S'], meas_df[time], yerr = np.ones(len(meas)), color='b', fmt=".", label='Measurements')
        ax_1.legend(fontsize=6, loc=2)        
        ax_1.set_ylabel('x [mm]')
        
    # Vertical
    else: 
        ax_1.plot(df_myTwiss.s, df_myTwiss.y*1E3, lw=0.5, color='k', label='Closed Orbit')
        ax_1.errorbar(meas_df['S'], meas_df[time], yerr = np.ones(len(meas)), color='b', fmt=".", label='Measurements')
        ax_1.legend(fontsize=6, loc=2)        
        ax_1.set_ylabel('y [mm]')

    if ylims:
        ax_1.set_ylim(ylims[0], ylims[1]);
        
    ax_1.grid(which='both', ls=':', lw=0.5, color='grey')
    #ax_1.set_ylim(min_kick-1, max_kick+1);

def cpymad_plot_CO_kicks_from_file(madx_instance, df_myTwiss, title, save_file, kicker_list, filename, horizontal=True, time='0', limits = None, ptc_twiss=False):
    
    if ptc_twiss:
        gamma_key = 'GAMMA'; pc_key='PC';
        ptc_twiss_read_Header = dict(df_myTwiss.headers)
        gamma_rel = ptc_twiss_read_Header[gamma_key]
        beta_rel = np.sqrt( 1. - (1./gamma_rel**2) )
        p_mass_GeV = 0.93827208816 #Proton mass GeV
        tot_energy = gamma_rel * p_mass_GeV
        kin_energy = tot_energy - p_mass_GeV
        momentum = ptc_twiss_read_Header[pc_key]

        print('Relativistic Gamma = ', round(gamma_rel,3))
        print('Relativistic Beta = ', round(beta_rel,3))
        print('Total Energy = ', round(tot_energy,4), 'GeV')
        print('Kinetic Energy = ', round(kin_energy*1E3,3), 'MeV')
        print('momentum = ', round(momentum,3), 'GeV/c')
        
        qx = ptc_twiss_read_Header['Q1']
        qy = ptc_twiss_read_Header['Q2']
    
    else:
        # Plot title = title + tunes
        qx = madx_instance.table.summ.q1[0]
        qy = madx_instance.table.summ.q2[0]     
    
    plot_title = title +' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
        
    # Start Plot
    fig = plt.figure(figsize=(10,8),facecolor='w', edgecolor='k',constrained_layout=True);
        
    heights = [1, 3, 4]
    spec = gridspec.GridSpec(ncols=1, nrows=3, figure=fig, height_ratios=heights)
    
    # Block diagram
    f_ax1 = fig.add_subplot(spec[0]);
    f_ax1.set_title(plot_title);  
        
    if limits is not None:
        if len(limits) != 2:
            print('cpymad_plot_CO::ERROR, limits must be given as a 2 variable list such as [0., 1.]')
            exit()
        if ptc_twiss:
            block_diagram(f_ax1, df_myTwiss, limits, ptc_twiss=True);
        else: 
            block_diagram(f_ax1, df_myTwiss, limits, ptc_twiss=False);
    else:
        if ptc_twiss:
            block_diagram(f_ax1, df_myTwiss, ptc_twiss=True);
        else: 
            block_diagram(f_ax1, df_myTwiss, ptc_twiss=False);
    
    # Plot betas   
    #-----------------------
    f_ax2 = fig.add_subplot(spec[1], sharex=f_ax1); 
    f_ax2.plot(df_myTwiss['s'], df_myTwiss['betx'],'b', label='$\\beta_x$');
    f_ax2.plot(df_myTwiss['s'], df_myTwiss['bety'],'r', label='$\\beta_y$');    
    
    f_ax2.legend(loc=2);
    f_ax2.set_ylabel(r'$\beta_{x,y}$[m]');
    f_ax2.grid(which='both', ls=':', lw=0.5, color='k');
    #f_ax2.set_xlabel('s [m]')
    #f_ax2.set_xticklabels([])
    
    if np.min(df_myTwiss['bety']) < np.min(df_myTwiss['betx']): bet_min = round_down_n(np.min(df_myTwiss['bety']),5)
    else: bet_min = round_down_n(np.min(df_myTwiss['betx']),5)
    if np.max(df_myTwiss['bety']) > np.max(df_myTwiss['betx']): bet_max = round_up_n(np.max(df_myTwiss['bety']),10)
    else: bet_max = round_up_n(np.max(df_myTwiss['betx']),10)        
    f_ax2.set_ylim(bet_min,bet_max)
    
    # Plot dispersion   
    #-----------------------
    ax2 = f_ax2.twinx()   # instantiate a second axes that shares the same x-axis
    if ptc_twiss:     
        ax2.plot(df_myTwiss['s'], df_myTwiss['disp1']/beta_rel,'green', label='$D_x$');
        ax2.plot(df_myTwiss['s'], df_myTwiss['disp3']/beta_rel,'purple', label='$D_y$');
        key_dx = 'disp1';        key_dy = 'disp3';  
    
    else:
        ax2.plot(df_myTwiss['s'], df_myTwiss['dx'],'green', label='$D_x$');
        ax2.plot(df_myTwiss['s'], df_myTwiss['dy'],'purple', label='$D_y$');
        key_dx = 'dx';        key_dy = 'dy';  
        
    ax2.legend(loc=1);
    ax2.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor='green')
    ax2.grid(which='both', ls=':', lw=0.5, color='green')

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down_n(np.min(df_myTwiss[key_dy]),1)
    else: d_min = round_down_n(np.min(df_myTwiss[key_dx]),1)    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_n(np.max(df_myTwiss[key_dy]),10)
    else: d_max = round_up_n(np.max(df_myTwiss[key_dx]),10) 
    ax2.set_ylim(d_min,d_max)   
    
    if horizontal:
    # Plot closed orbit x   
    #-----------------------
        f_ax3 = fig.add_subplot(spec[2], sharex=f_ax1);
        f_ax3.plot(df_myTwiss['s'], df_myTwiss['x']*1E3,'k', lw=1.5, label='$x$');      

        df_meas = processed_bpm_dat_to_df(filename)
        #s_list, meas_list = measurements_from_file(filename, time)
        f_ax3.errorbar(df_meas['S'], df_meas[time], yerr = np.ones(len(df_meas['S'])), color='b', fmt=".", label='Measurements');

        f_ax3.legend(loc=0);
        f_ax3.set_ylabel('x [mm]')
        f_ax3.grid(which='both', ls=':', lw=0.5, color='k');
        f_ax3.set_ylim(-40,40)
        f_ax3.set_xlabel('s [m]')
            
    # Plot kicks x  
    #----------------------- 
        f_ax5 = f_ax3.twinx()
        
        for kname in kicker_list:     
            if 'qtfk' in kname:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='magenta', alpha=0.5);
            elif 'qtdk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='magenta', alpha=0.5);
            elif 'qfk' in kname:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='orange', alpha=0.5);
            elif 'qdk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='orange', alpha=0.5);
            elif 'mdsk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='yellow', alpha=0.5);
            elif 'mdek' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='yellow', alpha=0.5);
            elif 'qdsk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='green', alpha=0.5);
            else:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='k', alpha=0.5);

        custom_lines = [Line2D([0], [0], color='k', lw=4), Line2D([0], [0], color='yellow', lw=4), Line2D([0], [0], color='magenta', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='green', lw=4)]
        f_ax5.legend(custom_lines, ['Steering Magnets', 'Main Dipole Ends', 'Trim Quads', 'Main Quads', 'Singlet Quad'], loc=1);
   
        f_ax5.set_ylim(-5,5)
        f_ax5.set_ylabel('Kick [mrad]', color='g')  # we already handled the x-label with ax1
        f_ax5.tick_params(axis='y', labelcolor='g')

    else:
    # Plot closed orbit y   
    #-----------------------   
        f_ax3 = fig.add_subplot(spec[2], sharex=f_ax1);
        f_ax3.plot(df_myTwiss['s'], df_myTwiss['y']*1E3,'k', lw=1.5, label='$y$');    

        df_meas = processed_bpm_dat_to_df(filename)
        #s_list, meas_list = measurements_from_file(filename, time)
        f_ax3.errorbar(df_meas['S'], df_meas[time], yerr = np.ones(len(df_meas['S'])), color='b', fmt=".", label='Measurements');

        f_ax3.legend(loc=0)
        f_ax3.set_ylabel('y [mm]')
        f_ax3.grid(which='both', ls=':', lw=0.5, color='k');
        f_ax3.set_ylim(-40,40)
        f_ax3.set_xlabel('s [m]')
        
    # Plot kicks y  
    #-----------------------          
        f_ax5 = f_ax3.twinx()    
        
        for kname in kicker_list:     
            if 'qtfk' in kname:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='magenta', alpha=0.5);
            elif 'qtdk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='magenta', alpha=0.5);
            elif 'qfk' in kname:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='orange', alpha=0.5);
            elif 'qdk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='orange', alpha=0.5);
            elif 'mdsk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='yellow', alpha=0.5);
            elif 'mdek' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='yellow', alpha=0.5);
            elif 'qdsk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='green', alpha=0.5);
            else:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='k', alpha=0.5);
        
        custom_lines = [Line2D([0], [0], color='k', lw=4), Line2D([0], [0], color='yellow', lw=4), Line2D([0], [0], color='magenta', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='green', lw=4)]
        f_ax5.legend(custom_lines, ['Steering Magnets', 'Main Dipole Ends', 'Trim Quads', 'Main Quads', 'Singlet Quad'], loc=1);
    
        f_ax5.set_ylim(-5,5)
        f_ax5.set_ylabel('Kick [mrad]', color='g')  # we already handled the x-label with ax1
        f_ax5.tick_params(axis='y', labelcolor='g')

    # Auto-limits mess
    #if np.min(df_myTwiss['y']) < np.min(df_myTwiss['x']): co_min = round_down_n(np.min(df_myTwiss['y']*1E3),10)
    #else: co_min = round_down_n(np.min(df_myTwiss['x']*1E3),10)
    #if np.max(df_myTwiss['y']) > np.max(df_myTwiss['x']): co_max = round_up_n(np.max(df_myTwiss['y']*1E3),10)
    #else: co_max = round_up_n(np.max(df_myTwiss['x']*1E3),10)
    #f_ax4.set_ylim(co_min,co_max)
    #f_ax3.set_ylim(co_min,co_max)
    
    if save_file != None: plt.savefig(save_file)
    
def cpymad_plot_CO_kicks_from_file_phase(madx_instance, df_myTwiss, title, save_file, kicker_list, filename=None, horizontal=True, time='0', limits = None, ptc_twiss=False):
    
    if ptc_twiss:
        gamma_key = 'GAMMA'; pc_key='PC';
        ptc_twiss_read_Header = dict(df_myTwiss.headers)
        gamma_rel = ptc_twiss_read_Header[gamma_key]
        beta_rel = np.sqrt( 1. - (1./gamma_rel**2) )
        p_mass_GeV = 0.93827208816 #Proton mass GeV
        tot_energy = gamma_rel * p_mass_GeV
        kin_energy = tot_energy - p_mass_GeV
        momentum = ptc_twiss_read_Header[pc_key]

        print('Relativistic Gamma = ', round(gamma_rel,3))
        print('Relativistic Beta = ', round(beta_rel,3))
        print('Total Energy = ', round(tot_energy,4), 'GeV')
        print('Kinetic Energy = ', round(kin_energy*1E3,3), 'MeV')
        print('momentum = ', round(momentum,3), 'GeV/c')
        
        qx = ptc_twiss_read_Header['Q1']
        qy = ptc_twiss_read_Header['Q2']
    
    else:
        # Plot title = title + tunes
        qx = madx_instance.table.summ.q1[0]
        qy = madx_instance.table.summ.q2[0]   
        
    if qx <= 1.: qx += 4.
    if qy <= 1.: qy += 3.
    
    plot_title = title +' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
        
    # Start Plot
    fig = plt.figure(figsize=(10,8),facecolor='w', edgecolor='k',constrained_layout=True);
        
    heights = [1, 3, 4]
    spec = gridspec.GridSpec(ncols=1, nrows=3, figure=fig, height_ratios=heights)
    
    # Block diagram
    f_ax1 = fig.add_subplot(spec[0]);
    f_ax1.set_title(plot_title);  
        
    if limits is not None:
        if len(limits) != 2:
            print('cpymad_plot_CO::ERROR, limits must be given as a 2 variable list such as [0., 1.]')
            exit()
        if ptc_twiss:
            block_diagram(f_ax1, df_myTwiss, limits, ptc_twiss=True);
        else: 
            block_diagram(f_ax1, df_myTwiss, limits, ptc_twiss=False);
    else:
        if ptc_twiss:
            block_diagram(f_ax1, df_myTwiss, ptc_twiss=True);
        else: 
            block_diagram(f_ax1, df_myTwiss, ptc_twiss=False);
    
    # Plot betas   
    #-----------------------
    f_ax2 = fig.add_subplot(spec[1]);#, sharex=f_ax1); 
    f_ax2.plot(df_myTwiss['mux'], df_myTwiss['betx'],'b', label='$\\beta_x$');
    f_ax2.plot(df_myTwiss['muy'], df_myTwiss['bety'],'r', label='$\\beta_y$');    
    
    f_ax2.legend(loc=2);
    f_ax2.set_ylabel(r'$\beta_{x,y}$[m]');
    f_ax2.grid(which='both', ls=':', lw=0.5, color='k');
    #f_ax2.set_xlabel('s [m]')
    #f_ax2.set_xticklabels([])
    
    if np.min(df_myTwiss['bety']) < np.min(df_myTwiss['betx']): bet_min = round_down_n(np.min(df_myTwiss['bety']),5)
    else: bet_min = round_down_n(np.min(df_myTwiss['betx']),5)
    if np.max(df_myTwiss['bety']) > np.max(df_myTwiss['betx']): bet_max = round_up_n(np.max(df_myTwiss['bety']),10)
    else: bet_max = round_up_n(np.max(df_myTwiss['betx']),10)        
    f_ax2.set_ylim(bet_min,bet_max)
    
    # Plot dispersion   
    #-----------------------
    ax2 = f_ax2.twinx()   # instantiate a second axes that shares the same x-axis
    if ptc_twiss:     
        ax2.plot(df_myTwiss['mux'], df_myTwiss['disp1']/beta_rel,'green', label='$D_x$');
        ax2.plot(df_myTwiss['muy'], df_myTwiss['disp3']/beta_rel,'purple', label='$D_y$');
        key_dx = 'disp1';        key_dy = 'disp3';  
    
    else:
        ax2.plot(df_myTwiss['mux'], df_myTwiss['dx'],'green', label='$D_x$');
        ax2.plot(df_myTwiss['muy'], df_myTwiss['dy'],'purple', label='$D_y$');
        key_dx = 'dx';        key_dy = 'dy';  
        
    ax2.legend(loc=1);
    ax2.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor='green')
    ax2.grid(which='both', ls=':', lw=0.5, color='green')

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down_n(np.min(df_myTwiss[key_dy]),1)
    else: d_min = round_down_n(np.min(df_myTwiss[key_dx]),1)    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_n(np.max(df_myTwiss[key_dy]),10)
    else: d_max = round_up_n(np.max(df_myTwiss[key_dx]),10) 
    ax2.set_ylim(d_min,d_max)   
    
    if horizontal:
    # Plot closed orbit x   
    #-----------------------
        f_ax3 = fig.add_subplot(spec[2], sharex=f_ax2);
        f_ax3.plot(df_myTwiss['mux'], df_myTwiss['x']*1E3/np.sqrt(df_myTwiss['betx']),'k', lw=1.5, label='$x$');      

        if filename:
            df_meas = processed_bpm_dat_to_df(filename)
            #s_list, meas_list = measurements_from_file(filename, time)
            f_ax3.errorbar(df_meas['Mu'], df_meas[time]/np.sqrt(df_meas['Beta']), yerr = np.ones(len(df_meas['Mu']))/np.sqrt(df_meas['Beta']), color='b', fmt=".", label='Measurements');

        f_ax3.legend(loc=0);
        f_ax3.set_ylabel(r'$\frac{x}{\sqrt{\beta_x}}$ [$10^{-3}\sqrt{m}$]')
        f_ax3.grid(which='both', ls=':', lw=0.5, color='k');
        #f_ax3.set_ylim(-40,40)
        f_ax3.set_xlabel(r'$\mu_x$ [-]')
            
    # Plot kicks x  
    #----------------------- 
        f_ax5 = f_ax3.twinx()
        
        for kname in kicker_list:     
            if 'qtfk' in kname:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].mux, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, qx/163., color='magenta', alpha=0.5);
            elif 'qtdk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].mux, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, qx/163., color='magenta', alpha=0.5);
            elif 'qfk' in kname:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].mux, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, qx/163., color='orange', alpha=0.5);
            elif 'qdk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].mux, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, qx/163., color='orange', alpha=0.5);
            elif 'mdsk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].mux, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, qx/163., color='yellow', alpha=0.5);
            elif 'mdek' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].mux, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, qx/163., color='yellow', alpha=0.5);
            elif 'qdsk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].mux, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, qx/163., color='green', alpha=0.5);
            else:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].mux, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, qx/163., color='k', alpha=0.5);

        custom_lines = [Line2D([0], [0], color='k', lw=4), Line2D([0], [0], color='yellow', lw=4), Line2D([0], [0], color='magenta', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='green', lw=4)]
        f_ax5.legend(custom_lines, ['Steering Magnets', 'Main Dipole Ends', 'Trim Quads', 'Main Quads', 'Singlet Quad'], loc=1);
   
        f_ax5.set_ylim(-5,5)
        f_ax5.set_ylabel('Kick [mrad]', color='g')  # we already handled the x-label with ax1
        f_ax5.tick_params(axis='y', labelcolor='g')

    else:
    # Plot closed orbit y   
    #-----------------------   
        f_ax3 = fig.add_subplot(spec[2], sharex=f_ax2);
        f_ax3.plot(df_myTwiss['muy'], df_myTwiss['y']*1E3/np.sqrt(df_myTwiss['bety']),'k', lw=1.5, label='$y$');    

        df_meas = processed_bpm_dat_to_df(filename)
        #s_list, meas_list = measurements_from_file(filename, time)
        f_ax3.errorbar(df_meas['Mu'], df_meas[time]/np.sqrt(df_meas['Beta']), yerr = np.ones(len(df_meas['Mu']))/np.sqrt(df_meas['Beta']), color='b', fmt=".", label='Measurements');

        f_ax3.legend(loc=0)
        f_ax3.set_ylabel(r'$\frac{y}{\sqrt{\beta_y}}$ [$10^{-3}\sqrt{m}$]')
        f_ax3.grid(which='both', ls=':', lw=0.5, color='k');
        #f_ax3.set_ylim(-40,40)
        f_ax3.set_xlabel(r'$\mu_y$ [-]')
        
    # Plot kicks y  
    #-----------------------          
        f_ax5 = f_ax3.twinx()    
        
        for kname in kicker_list:     
            if 'qtfk' in kname:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].muy, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, qy/163., color='magenta', alpha=0.5);
            elif 'qtdk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].muy, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, qy/163., color='magenta', alpha=0.5);
            elif 'qfk' in kname:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].muy, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, qy/163., color='orange', alpha=0.5);
            elif 'qdk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].muy, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, qy/163., color='orange', alpha=0.5);
            elif 'mdsk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].muy, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, qy/163., color='yellow', alpha=0.5);
            elif 'mdek' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].muy, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, qy/163., color='yellow', alpha=0.5);
            elif 'qdsk' in kname:    
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].muy, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, qy/163., color='green', alpha=0.5);
            else:
                f_ax5.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].muy, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, qy/163., color='k', alpha=0.5);
        
        custom_lines = [Line2D([0], [0], color='k', lw=4), Line2D([0], [0], color='yellow', lw=4), Line2D([0], [0], color='magenta', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='green', lw=4)]
        f_ax5.legend(custom_lines, ['Steering Magnets', 'Main Dipole Ends', 'Trim Quads', 'Main Quads', 'Singlet Quad'], loc=1);
    
        f_ax5.set_ylim(-5,5)
        f_ax5.set_ylabel('Kick [mrad]', color='g')  # we already handled the x-label with ax1
        f_ax5.tick_params(axis='y', labelcolor='g')

    # Auto-limits mess
    #if np.min(df_myTwiss['y']) < np.min(df_myTwiss['x']): co_min = round_down_n(np.min(df_myTwiss['y']*1E3),10)
    #else: co_min = round_down_n(np.min(df_myTwiss['x']*1E3),10)
    #if np.max(df_myTwiss['y']) > np.max(df_myTwiss['x']): co_max = round_up_n(np.max(df_myTwiss['y']*1E3),10)
    #else: co_max = round_up_n(np.max(df_myTwiss['x']*1E3),10)
    #f_ax4.set_ylim(co_min,co_max)
    #f_ax3.set_ylim(co_min,co_max)
    
    if save_file != None: plt.savefig(save_file)

def cpymad_plot_CO_kicks_from_file_xy(madx_instance, df_myTwiss, title, save_file, kicker_list_x, kicker_list_y, filename, time='0', horizontal=True, limits = None, ptc_twiss=False):
    
    if ptc_twiss:
        gamma_key = 'GAMMA'; pc_key='PC';
        ptc_twiss_read_Header = dict(df_myTwiss.headers)
        gamma_rel = ptc_twiss_read_Header[gamma_key]
        beta_rel = np.sqrt( 1. - (1./gamma_rel**2) )
        p_mass_GeV = 0.93827208816 #Proton mass GeV
        tot_energy = gamma_rel * p_mass_GeV
        kin_energy = tot_energy - p_mass_GeV
        momentum = ptc_twiss_read_Header[pc_key]

        print('Relativistic Gamma = ', round(gamma_rel,3))
        print('Relativistic Beta = ', round(beta_rel,3))
        print('Total Energy = ', round(tot_energy,4), 'GeV')
        print('Kinetic Energy = ', round(kin_energy*1E3,3), 'MeV')
        print('momentum = ', round(momentum,3), 'GeV/c')
        
        qx = ptc_twiss_read_Header['Q1']
        qy = ptc_twiss_read_Header['Q2']
    
    else:
        # Plot title = title + tunes
        qx = madx_instance.table.summ.q1[0]
        qy = madx_instance.table.summ.q2[0]     
    
    plot_title = title +' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
        
    # Start Plot
    fig = plt.figure(figsize=(10,8),facecolor='w', edgecolor='k',constrained_layout=True);
        
    heights = [1, 1, 3, 3]
    spec = gridspec.GridSpec(ncols=1, nrows=4, figure=fig, height_ratios=heights)
    
    # Block diagram
    f_ax1 = fig.add_subplot(spec[0]);
    f_ax1.set_title(plot_title);  
        
    if limits is not None:
        if len(limits) != 2:
            print('cpymad_plot_CO::ERROR, limits must be given as a 2 variable list such as [0., 1.]')
            exit()
        if ptc_twiss:
            block_diagram(f_ax1, df_myTwiss, limits, ptc_twiss=True);
        else: 
            block_diagram(f_ax1, df_myTwiss, limits, ptc_twiss=False);
    else:
        if ptc_twiss:
            block_diagram(f_ax1, df_myTwiss, ptc_twiss=True);
        else: 
            block_diagram(f_ax1, df_myTwiss, ptc_twiss=False);
    
    # Plot betas   
    #-----------------------
    f_ax2 = fig.add_subplot(spec[1], sharex=f_ax1); 
    f_ax2.plot(df_myTwiss['s'], df_myTwiss['betx'],'b', label='$\\beta_x$');
    f_ax2.plot(df_myTwiss['s'], df_myTwiss['bety'],'r', label='$\\beta_y$');    
    
    f_ax2.legend(loc=2);
    f_ax2.set_ylabel(r'$\beta_{x,y}$[m]');
    f_ax2.grid(which='both', ls=':', lw=0.5, color='k');
    #f_ax2.set_xlabel('s [m]')
    #f_ax2.set_xticklabels([])
    
    if np.min(df_myTwiss['bety']) < np.min(df_myTwiss['betx']): bet_min = round_down_n(np.min(df_myTwiss['bety']),5)
    else: bet_min = round_down_n(np.min(df_myTwiss['betx']),5)
    if np.max(df_myTwiss['bety']) > np.max(df_myTwiss['betx']): bet_max = round_up_n(np.max(df_myTwiss['bety']),10)
    else: bet_max = round_up_n(np.max(df_myTwiss['betx']),10)        
    f_ax2.set_ylim(bet_min,bet_max)
    
    # Plot dispersion   
    #-----------------------
    ax2 = f_ax2.twinx()   # instantiate a second axes that shares the same x-axis
    if ptc_twiss:     
        ax2.plot(df_myTwiss['s'], df_myTwiss['disp1']/beta_rel,'green', label='$D_x$');
        ax2.plot(df_myTwiss['s'], df_myTwiss['disp3']/beta_rel,'purple', label='$D_y$');
        key_dx = 'disp1';        key_dy = 'disp3';  
    
    else:
        ax2.plot(df_myTwiss['s'], df_myTwiss['dx'],'green', label='$D_x$');
        ax2.plot(df_myTwiss['s'], df_myTwiss['dy'],'purple', label='$D_y$');
        key_dx = 'dx';        key_dy = 'dy';  
        
    ax2.legend(loc=1);
    ax2.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor='green')
    ax2.grid(which='both', ls=':', lw=0.5, color='green')

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down_n(np.min(df_myTwiss[key_dy]),1)
    else: d_min = round_down_n(np.min(df_myTwiss[key_dx]),1)    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_n(np.max(df_myTwiss[key_dy]),10)
    else: d_max = round_up_n(np.max(df_myTwiss[key_dx]),10) 
    ax2.set_ylim(d_min,d_max)   
    
    # Plot closed orbit x   
    #-----------------------
    f_ax3 = fig.add_subplot(spec[2], sharex=f_ax1);
    if (df_myTwiss['x'] == 0).all() and df_myTwiss['x'].std()==0:
            print('WARNING: cpymad_plot_CO_kicks_from_file_xy: twiss has perfect closed orbit in this plane')
            
    f_ax3.plot(df_myTwiss['s'], df_myTwiss['x']*1E3,'k', lw=1.5, label='$x$');      
    
    if horizontal:
        df_meas = processed_bpm_dat_to_df(filename)
        #s_list, meas_list = measurements_from_file(filename, time)
        f_ax3.errorbar(df_meas['S'], df_meas[time], yerr = np.ones(len(df_meas['S'])), color='b', fmt=".", label='Measurements');

    f_ax3.legend(loc=0);
    f_ax3.set_ylabel('x [mm]')
    f_ax3.grid(which='both', ls=':', lw=0.5, color='k');
    f_ax3.set_ylim(-40,40)
    f_ax3.set_xlabel('s [m]')

    # Plot kicks x  
    #----------------------- 
    f_ax4 = f_ax3.twinx()

    for kname in kicker_list_x:     
        if 'qtfk' in kname:
            f_ax4.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='magenta', alpha=0.5);
        elif 'qtdk' in kname:    
            f_ax4.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='magenta', alpha=0.5);
        elif 'qfk' in kname:
            f_ax4.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='orange', alpha=0.5);
        elif 'qdk' in kname:    
            f_ax4.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='orange', alpha=0.5);
        elif 'mdsk' in kname:    
            f_ax4.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='yellow', alpha=0.5);
        elif 'mdek' in kname:    
            f_ax4.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='yellow', alpha=0.5);
        elif 'qdsk' in kname:    
            f_ax4.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='green', alpha=0.5);
        else:
            f_ax4.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].hkick*1E3, 1, color='k', alpha=0.5);

    custom_lines = [Line2D([0], [0], color='k', lw=4), Line2D([0], [0], color='yellow', lw=4), Line2D([0], [0], color='magenta', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='green', lw=4)]
    f_ax4.legend(custom_lines, ['Steering Magnets', 'Main Dipole Ends', 'Trim Quads', 'Main Quads', 'Singlet Quad'], loc=1);

    f_ax4.set_ylim(-5,5)
    f_ax4.set_ylabel('Kick [mrad]', color='g')  # we already handled the x-label with ax1
    f_ax4.tick_params(axis='y', labelcolor='g')

    # Plot closed orbit y   
    #-----------------------   
    f_ax5 = fig.add_subplot(spec[3], sharex=f_ax1);
    f_ax5.plot(df_myTwiss['s'], df_myTwiss['y']*1E3,'k', lw=1.5, label='$y$');    

    if not horizontal:
        df_meas = processed_bpm_dat_to_df(filename)
        #s_list, meas_list = measurements_from_file(filename, 'Monitor', False)
        f_ax5.errorbar(df_meas['S'], df_meas[time], yerr = np.ones(len(df_meas['S'])), color='b', fmt=".", label='Measurements');

    f_ax5.legend(loc=0)
    f_ax5.set_ylabel('y [mm]')
    f_ax5.grid(which='both', ls=':', lw=0.5, color='k');
    f_ax5.set_ylim(-40,40)
    f_ax5.set_xlabel('s [m]')

    # Plot kicks y  
    #-----------------------          
    f_ax6 = f_ax5.twinx()    

    for kname in kicker_list_y:     
        if 'qtfk' in kname:
            f_ax6.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='magenta', alpha=0.5);
        elif 'qtdk' in kname:    
            f_ax6.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='magenta', alpha=0.5);
        elif 'qfk' in kname:
            f_ax6.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='orange', alpha=0.5);
        elif 'qdk' in kname:    
            f_ax6.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='orange', alpha=0.5);
        elif 'mdsk' in kname:    
            f_ax6.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='yellow', alpha=0.5);
        elif 'mdek' in kname:    
            f_ax6.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='yellow', alpha=0.5);
        elif 'qdsk' in kname:    
            f_ax6.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='green', alpha=0.5);
        else:
            f_ax6.bar(df_myTwiss[df_myTwiss.name.str.contains(kname)].s, df_myTwiss[df_myTwiss.name.str.contains(kname)].vkick*1E3, 1, color='k', alpha=0.5);

    custom_lines = [Line2D([0], [0], color='k', lw=4), Line2D([0], [0], color='yellow', lw=4), Line2D([0], [0], color='magenta', lw=4), Line2D([0], [0], color='orange', lw=4), Line2D([0], [0], color='green', lw=4)]
    f_ax6.legend(custom_lines, ['Steering Magnets', 'Main Dipole Ends', 'Trim Quads', 'Main Quads', 'Singlet Quad'], loc=1);

    f_ax6.set_ylim(-5,5)
    f_ax6.set_ylabel('Kick [mrad]', color='g')  # we already handled the x-label with ax1
    f_ax6.tick_params(axis='y', labelcolor='g')

    # Auto-limits mess
    #if np.min(df_myTwiss['y']) < np.min(df_myTwiss['x']): co_min = round_down_n(np.min(df_myTwiss['y']*1E3),10)
    #else: co_min = round_down_n(np.min(df_myTwiss['x']*1E3),10)
    #if np.max(df_myTwiss['y']) > np.max(df_myTwiss['x']): co_max = round_up_n(np.max(df_myTwiss['y']*1E3),10)
    #else: co_max = round_up_n(np.max(df_myTwiss['x']*1E3),10)
    #f_ax4.set_ylim(co_min,co_max)
    #f_ax3.set_ylim(co_min,co_max)
    
    if save_file != None: plt.savefig(save_file)
