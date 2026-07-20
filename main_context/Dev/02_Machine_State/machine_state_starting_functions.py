import numpy as np
from cpymad.types import Constraint


# Loading of ISIS lattice files to MAD-X

cpymad_logfile = 'cpymad_logfile.txt'
lattice_folder = '../../Lattice_Files/00_Simplified_Lattice/'
# lattice_folder = '../../Lattice_Files/01_Original_Lattice/'
#lattice_folder = '../../Lattice_Files/02_Aperture_Lattice/'
#lattice_folder = '../../Lattice_Files/03_CO_Kick_Lattice/'
#lattice_folder = '../../Lattice_Files/04_New_Harmonics/'
sequence_name = 'synchrotron'
madx = cpymad_start(cpymad_logfile)
madx.call(file=lattice_folder+'ISIS.injected_beam')
madx.call(file=lattice_folder+'ISIS.strength')
madx.call(file=lattice_folder+'2023.strength') # including this file essentially switches from RCS to SRM
madx.call(file=lattice_folder+'ISIS.elements')
madx.call(file=lattice_folder+'ISIS.sequence')


########################################################################
# check if sequence exists and use it (required for madx)
########################################################################    
def cpymad_check_and_use_sequence(madx_instance, cpymad_logfile, sequence_name):     
        if sequence_name in cpymad_get_active_sequence(madx_instance):
            madx_instance.use(sequence=sequence_name)
            print('Sequence ', str(sequence_name), ' is active.')
            return True
        else:         
            madx_instance.use(sequence=sequence_name)
            if 'warning' and sequence_name in cpymad_get_output(cpymad_logfile)[0][-1]:
                print(cpymad_get_output(cpymad_logfile)[0][-1])
                print('cpymad_check_and_use_sequence::Sequence not valid in this instance of MAD-X')           
                log_string = '! cpymad_check_and_use_sequence called for sequence ' + sequence_name
                cpymad_write_to_logfile(cpymad_logfile, log_string)  
                return False
            else: 
                print('Sequence',sequence_name,'exists in this instance of MAD-X. Active sequences:')
                print(cpymad_get_active_sequence(madx_instance))       
                log_string = '! cpymad_check_and_use_sequence called for sequence ' + sequence_name
                cpymad_write_to_logfile(cpymad_logfile, log_string)  
                return True
                
# plot the twiss and magnet strengths:
def cpymad_plot_madx_twiss_quads(madx_instance, df_myTwiss, title=None, savename=None, limits=None, use_caps=False):
        
    fig = plt.figure(figsize=(13,8),facecolor='w', edgecolor='k')
    
    ptc = False
    if ptc:        
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
        
        if 'ptc_twiss_summary' in list(madx_instance.table):
            qx = madx_instance.table.ptc_twiss_summary.q1[0]
            qy = madx_instance.table.ptc_twiss_summary.q2[0]
            dqx = madx_instance.table.ptc_twiss_summary.dq1[0]
            dqy = madx_instance.table.ptc_twiss_summary.dq2[0] 
              
            if title is None:        
                active_seq = str(cpymad_get_active_sequence(madx_instance)).split('\'')[1]
                plot_title = active_seq +r' Q$_x$='+format(qx,'2.3f')+r', Q$_y$='+ format(qy,'2.3f')+r', $\xi_x$='+ format(dqx,'2.3f') + r', $\xi_y$='+ format(dqy,'2.3f') 
            else: plot_title = title + r' Q$_x$='+format(qx,'2.3f')+r', Q$_y$='+ format(qy,'2.3f')+r', $\xi_x$='+ format(dqx,'2.3f') + r', $\xi_y$='+ format(dqy,'2.3f') 
        else:
            if title is None: plot_title = active_seq
            else: plot_title = title                 
    else:           
        if 'summ' in list(madx_instance.table):
            qx = madx_instance.table.summ.q1[0]
            qy = madx_instance.table.summ.q2[0]    
            dqx = madx_instance.table.summ.dq1[0]
            dqy = madx_instance.table.summ.dq2[0]    
        
            if title is None:        
                active_seq = str(cpymad_get_active_sequence(madx_instance)).split('\'')[1]
                plot_title = active_seq +r' Q$_x$='+format(qx,'2.3f')+r', Q$_y$='+ format(qy,'2.3f')+r', $\xi_x$='+ format(dqx,'2.3f') + r', $\xi_y$='+ format(dqy,'2.3f') 
            else: plot_title = title + r' Q$_x$='+format(qx,'2.3f')+r', Q$_y$='+ format(qy,'2.3f')+r', $\xi_x$='+ format(dqx,'2.3f') + r', $\xi_y$='+ format(dqy,'2.3f') 
        else:
            if title is None: plot_title = active_seq
            else: plot_title = title 
        

    ax1=plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=1)      
    plt.title(plot_title) 
    
    if use_caps is False:
        try:    plt.plot(df_myTwiss['s'], 0*df_myTwiss['s'],'k')
        except KeyError: 
            try: 
                plt.plot(df_myTwiss['S'], 0*df_myTwiss['S'],'k')
                use_caps = True
                print('cpymad_plotTwiss::use_caps = True')
            except: print('cpymad_plotTwiss::unkown bug')
    else:
        plt.plot(df_myTwiss['S'], 0*df_myTwiss['S'],'k')        

    if use_caps: 
        s_key = 'S'
        keyword = 'KEYWORD'        
    else: 
        s_key =  's'
        keyword = 'keyword'   
    
    quad_max = 0.
    dipole_max = 0.
    
    if use_caps: key = 'QUADRUPOLE'
    else: key =  'quadrupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            plotLatticeSeriesCaps(plt.gca(),aux, height=aux.K1L, v_offset=aux.K1L/2, color='r')
            if np.max(abs(aux.K1L)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.K1L)))
        else: 
            plotLatticeSeries(plt.gca(),aux, height=aux.k1l, v_offset=aux.k1l/2, color='r')
            if np.max(abs(aux.k1l)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.k1l)))
    
    if use_caps: key = 'MULTIPOLE' 
    else: key =  'multipole' 
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            plotLatticeSeriesCaps(plt.gca(),aux, height=aux.K1L, v_offset=aux.K1L/2, color='r')
            if np.max(abs(aux.K1L)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.K1L)))
        else: 
            plotLatticeSeries(plt.gca(),aux, height=aux.k1l, v_offset=aux.k1l/2, color='r')
            if np.max(abs(aux.k1l)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.k1l)))
 
    if use_caps: key = 'SBEND' 
    else: key =  'sbend' 
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            plotLatticeSeriesCaps(plt.gca(),aux, height=aux.K1L, v_offset=aux.K1L/2, color='r')
            if np.max(abs(aux.K1L)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.K1L)))
        else: 
            plotLatticeSeries(plt.gca(),aux, height=aux.k1l, v_offset=aux.k1l/2, color='r')
            if np.max(abs(aux.k1l)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.k1l)))
   

    #plt.ylim(-.065,0.065)
    color = 'red'
    ax1.set_ylabel('Main Quads 1/f=K1L [m$^{-1}$]', color=color, fontsize='small')  # we already handled the x-label with ax1
    ax1.tick_params(axis='y', labelcolor=color)
    plt.grid()
    
    plt.ylim(-quad_max,quad_max)
    
    if limits is not None:
        if len(limits) != 2:
            print('cpymad_plot_madx_twiss::ERROR, limits must be given as a 2 variable list such as [0., 1.]')
            exit()
        ax1.set_xlim(limits[0], limits[1]);
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis   
    
    quad_max = 0.    
    color = 'blue'
    ax2.set_ylabel('Trim Quads 1/f=K1L [m$^{-1}$]', color=color, fontsize='small')  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    
    if use_caps: key = 'QUADRUPOLE'
    else: key =  'quadrupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps:             
            if 'QT' in aux.NAME:
                plotLatticeSeriesCaps(plt.gca(),aux, height=aux.K1L, v_offset=aux.K1L/2, color='b')
                if np.max(abs(aux.K1L)) > quad_max: quad_max = round_up_p01(np.max(abs(aux.K1L)))
        else: 
            if 'qt' in aux.name:
                plotLatticeSeries(plt.gca(),aux, height=aux.k1l, v_offset=aux.k1l/2, color='b')
                if np.max(abs(aux.k1l)) > quad_max: quad_max = round_up_p01(np.max(abs(aux.k1l)))
    if quad_max == 0.0: plt.ylim(-0.01, 0.01)
    else: plt.ylim(-quad_max, quad_max)


    # large subplot
    plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2,sharex=ax1)
    if use_caps: key_betx = 'BETX';        key_bety = 'BETY';
    else:        key_betx = 'betx';        key_bety = 'bety';        
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_betx],'b', label='$\\beta_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_bety],'r', label='$\\beta_y$')
    plt.legend(loc=2)
    plt.ylabel(r'$\beta_{x,y}$[m]')
    plt.xlabel('s [m]')
    plt.grid(which='both', ls=':', lw=0.5, color='k')
    
    if np.min(df_myTwiss[key_bety]) < np.min(df_myTwiss[key_betx]): bet_min = round_down_10(np.min(df_myTwiss[key_bety]))
    else: bet_min = round_down_10(np.min(df_myTwiss[key_betx]))
    if np.max(df_myTwiss[key_bety]) > np.max(df_myTwiss[key_betx]): bet_max = round_up_10(np.max(df_myTwiss[key_bety]))
    else: bet_max = round_up_10(np.max(df_myTwiss[key_betx]))        
    plt.ylim(bet_min,bet_max)

    ax3 = plt.gca().twinx()   # instantiate a second axes that shares the same x-axis
    if use_caps: key_dx = 'DX';        key_dy = 'DY';
    else:        key_dx = 'dx';        key_dy = 'dy';  
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dx],'green', label='$D_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dy],'purple', label='$D_y$')
    ax3.legend(loc=1)
    ax3.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax3.tick_params(axis='y', labelcolor='green')
    plt.grid(which='both', ls=':', lw=0.5, color='green')

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down(np.min(df_myTwiss[key_dy]))
    else: d_min = round_down(np.min(df_myTwiss[key_dx]))    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_10(np.max(df_myTwiss[key_dy]))
    else: d_max = round_up_10(np.max(df_myTwiss[key_dx]))        
    plt.ylim(d_min,d_max)
    #plt.ylim(round_down(np.min(df_myTwiss[key_dx])), round_up_10(np.max(df_myTwiss[key_dx])))
     
    if savename is None: pass
    else: plt.savefig(savename)
    

# Orbit correction functions


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
        
    #return round_sig(float(amps*(calibration_data[key]/df['Rigidity [Tm]']))) # kick in milliradians
    return round_sig(float(amps * (calibration_data[key] / df['Rigidity [Tm]'].iloc[0])))  # kick in milliradians


def cpymad_set_v_correctors(madx_instance, cpymad_logfile, corrector_dict, max_E=800., time=0.0):
    """
    Applies the vertical corrector kick values (converted from Amperes to mrad) 
    to a cpymad MAD-X instance.

    Parameters:
    madx_instance (Madx): An instance of cpymad's MAD-X.
    cpymad_logfile (str): Path to the cpymad log file (not used in function).
    corrector_dict (dict): Dictionary with keys as MAD-X variable names 
                           and values as the programmed kicks in Amperes.
    """

    for key, amps in corrector_dict.items():
        # Extract plane ('V' from 'vd1') and super-period (from 'rX' where X is 0, 2, etc.)
        sp = int(key[1])  # Extract the second character as the super-period
        plane = 'V' if 'vd' in key else 'H'  # Determine plane from key name

        # Convert kick from Amperes to milliradians
        kick_mrad = calculate_steering_kick(amps, max_E, time, plane, sp)

        # Print key, amps, and converted kick
        print(f"{key}: {amps:.6f} A -> {kick_mrad:.6f} mrad")

        # Apply the converted kick to MAD-X
        kick_mrad *= 1E-3 # convert from millirad to radians
        madx_instance.input(f"{key} := {kick_mrad};")

def calculate_corrector_current(kick_mrad, max_E, time, plane='H', sp=0):
    """
    Returns the corrector current in amperes given the desired steering kick in milliradians.

    Parameters:
    kick_mrad (float): Desired steering kick in milliradians.
    max_E (float): Maximum energy.
    time (float): Measurement time.
    plane (str): 'H' for horizontal or 'V' for vertical.
    sp (int): Super-period number.

    Returns:
    float: Corrector current in amperes.
    """
    sp_list = [0, 2, 3, 4, 5, 7, 9]
    if sp not in sp_list:
        print('calculate_corrector_current:: selected super-period has no steering magnet')
        exit(0)

    # Calibration provided by HVC 30.09.22
    calibration_data = {
        '0H': 0.08350, '2H': 0.09121, '3H': 0.08, '4H': 0.06600,
        '5H': 0.07780, '7H': 0.07580, '9H': 0.07660, '0V': 0.04620,
        '2V': 0.04330, '3V': 0.05210, '4V': 0.04770, '5V': 0.05400,
        '7V': 0.05220, '9V': 0.04510
    }

    df = synchrotron_energy_data(max_E, time)

    h_list = ['h', 'H', 'horizontal', 'Horizontal']
    key = f"{sp}{'H' if plane in h_list else 'V'}"

    # Compute the current in amperes
    amps = kick_mrad * df['Rigidity [Tm]'].iloc[0] / calibration_data[key]

    return round_sig(amps)

# Note that the madx lattice uses kicks but the controls system uses currents in amperes so both methods are required
def convert_kicks_to_currents(v_corrected_dict, max_E=800., time=0.0):
    """
    Converts a dictionary of corrector kicks in milliradians to corrector settings in amperes.

    Parameters:
    v_corrected_dict (dict): Dictionary where keys are corrector names and values are kicks in milliradians.
    max_E (float): Maximum energy, default is 800.
    time (float): Measurement time, default is 0.0.

    Returns:
    dict: A new dictionary with the same keys but values converted to amperes.
    """
    current_dict = {}

    for key, kick_mrad in v_corrected_dict.items():
        # Extract plane ('V' from 'vd1') and super-period (from 'rX' where X is 0, 2, etc.)
        sp = int(key[1])  # Extract super-period from the second character
        plane = 'V' if 'vd' in key else 'H'  # Determine plane based on key name

        # Convert the kick from milliradians to amperes
        amps = calculate_corrector_current(kick_mrad, max_E, time, plane, sp)

        # Store the converted value in the new dictionary
        current_dict[key] = amps

    return current_dict


def cpymad_set_isis_cycle_time(madx_instance, max_E, time):
    # Ensure time is a float and in valid increments
    if not isinstance(time, float) or time < 0.0 or time > 10.0 or (time * 10) % 5 != 0:
        print(f"Error: time must be a float between 0.0 and 10.0 in 0.5 increments. Received: {time}")
        return

    # Generate dataframe of synchrotron energy and related info
    energy_df = synchrotron_energy_df(max_E, intervals=20)

    # store some values for this time point
    try:
        energy = energy_df[energy_df['Time [ms]'] == time]['Energy [MeV]'].iloc[0]
        pc = energy_df[energy_df['Time [ms]'] == time]['Momentum [GeV/c]'].iloc[0]
    except IndexError:
        print(f"Error: No matching time value found in energy dataframe for time = {time} ms")
        return

    # set the beam to this energy in cpymad
    madx_instance.input(f'beam, particle = proton, pc = {pc};')

    # print confirmation
    print(f'ISIS cpymad run, energy set to {energy} MeV, pc = {pc}')

# Misalignment functions:

def get_madx_table_df(madx, table_name="efield", nonzero=True):
    """
    Extracts a MAD-X table as a pandas DataFrame.

    Parameters:
        madx : cpymad.madx.Madx
            The active MAD-X instance.
        table_name : str
            Name of the MAD-X table to extract (default: "efield").
        nonzero : bool
            If True, return only rows where any of the key misalignment/rotation columns are non-zero.

    Returns:
        pd.DataFrame
            DataFrame containing the filtered or full table data.
    """
    if table_name not in list(madx.table):
        raise ValueError(f"MAD-X table '{table_name}' not found. Available tables: {list(madx.table)}")

    raw_table = getattr(madx.table, table_name)
    raw_data = raw_table.copy()
    df = pd.DataFrame(raw_data, columns=raw_data.keys())

    if nonzero:
        cols = ["ds", "dx", "dy", "dtheta", "dphi", "dpsi"]
        df = df[(df[cols] != 0).any(axis=1)][["name"] + cols]

    return df

def cpymad_apply_error_table(madx_instance, error_table_file):
    """
    Apply a MAD-X error table file using madx_instance.input().

    Parameters:
        madx_instance: The cpymad.Madx instance.
        error_table_file: Path to the .tfs error table file.
    """
    madx_instance.input(f'READMYTABLE, file="{error_table_file}", table=efield;')
    madx_instance.input('SETERR, TABLE=efield;')
    
def cpymad_apply_error_tfs(madx_instance, error_tfs, atol=1e-10, rtol=1e-12):
    """
    Apply a dumped MAD-X error table (.tfs) to the active MAD-X sequence,
    then read back the resulting MAD-X error table and verify that it matches
    the contents of the input .tfs file.

    Parameters
    ----------
    madx_instance : cpymad.madx.Madx
        Active MAD-X instance.

    error_tfs : str
        Path to the dumped MAD-X error table file produced by
        build_and_save_madx_error_tables().

    atol : float, optional
        Absolute tolerance for comparison.

    rtol : float, optional
        Relative tolerance for comparison.

    Returns
    -------
    bool
        True if the applied MAD-X error table matches the input .tfs file
        within tolerance, otherwise raises AssertionError.
    """

    def _read_error_tfs_as_df(tfs_file):
        """
        Read a MAD-X TFS error table into a pandas DataFrame.
        Returns only the standard error columns in MAD-X units (m/rad).
        """
        with open(tfs_file, "r") as f:
            lines = f.readlines()

        header_line = None
        data_start = None

        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith("*"):
                header_line = stripped
            elif stripped.startswith("$") and header_line is not None:
                data_start = i + 1
                break

        if header_line is None or data_start is None:
            raise ValueError(f"Could not parse MAD-X TFS file: {tfs_file}")

        columns = header_line.lstrip("*").split()

        data_lines = []
        for line in lines[data_start:]:
            stripped = line.strip()
            if not stripped or stripped.startswith("@") or stripped.startswith("*") or stripped.startswith("$"):
                continue
            data_lines.append(stripped.split())

        tfs_df = pd.DataFrame(data_lines, columns=columns)

        required_cols = ["NAME", "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI"]
        missing = [col for col in required_cols if col not in tfs_df.columns]
        if missing:
            raise ValueError(f"Missing required columns in TFS file: {missing}")

        tfs_df = tfs_df[required_cols].copy()
        tfs_df.columns = [c.lower() for c in tfs_df.columns]
        tfs_df["name"] = tfs_df["name"].astype(str).str.strip('"').str.lower()

        for col in ["dx", "dy", "ds", "dphi", "dtheta", "dpsi"]:
            tfs_df[col] = pd.to_numeric(tfs_df[col], errors="coerce")

        return tfs_df

    # Apply error table to MAD-X
    cpymad_apply_error_table(madx_instance, error_tfs)

    # Read back MAD-X table in m/rad
    madx_df = get_madx_table_df(madx_instance, nonzero=False)
    madx_df = madx_df[["name", "dx", "dy", "ds", "dphi", "dtheta", "dpsi"]].copy()
    madx_df["name"] = madx_df["name"].astype(str).str.lower()

    for col in ["dx", "dy", "ds", "dphi", "dtheta", "dpsi"]:
        madx_df[col] = pd.to_numeric(madx_df[col], errors="coerce")

    # Read original TFS table in m/rad
    tfs_df = _read_error_tfs_as_df(error_tfs)

    # Sort both DataFrames
    madx_sorted = madx_df.sort_values("name").reset_index(drop=True)
    tfs_sorted = tfs_df.sort_values("name").reset_index(drop=True)

    # Compare with tolerance
    pd.testing.assert_frame_equal(
        madx_sorted,
        tfs_sorted,
        check_dtype=False,
        rtol=rtol,
        atol=atol,
    )

    return True


# Function to write error table:

class cpymad_ErrorTableBuilder:
    def __init__(self, twiss_df):
        """
        Initialize the error table builder with a MAD-X twiss table DataFrame.
        """
        self.twiss_df = twiss_df.copy()
        self.error_df = pd.DataFrame(columns=["NAME", "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI", "S"])

    def _get_twiss_lookup(self):
        """
        Prepare Twiss lookup with cleaned names (without colon suffixes).
        """
        return (
            self.twiss_df.copy()
            .assign(name_clean=self.twiss_df["name"].astype(str).str.split(":").str[0])
            .drop_duplicates(subset="name_clean")
            .set_index("name_clean")["s"]
        )

    def _match_central_magnet_only(self, base_name):
        """
        Return only the central element matching base_name, excluding fringe pieces.

        Allowed:
        - exact match, e.g. sp0_qd
        - colon suffixed exact form, e.g. sp0_qd:1 -> stored as sp0_qd

        Excluded:
        - sp0_qdfr1
        - sp0_dipfr8
        - sp0_qds when base_name is sp0_qd
        """
        clean_names = self.twiss_df["name"].astype(str).str.split(":").str[0]
        matched = clean_names[clean_names == base_name]
        return matched.drop_duplicates().tolist()

    def match_main_magnet_parts(self, base_name):
        """
        Return only the central magnet element matching base_name.
        """
        return self._match_central_magnet_only(base_name)

    def add_dipole_misalignment(self, base_name, misalignment_type, value_mm):
        """
        Add a misalignment value (in mm or mrad) for the central dipole element only.

        Parameters:
            - base_name: base name of dipole (e.g. "sp6_dip")
            - misalignment_type: one of "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI"
            - value_mm: float, Magnitude in mm (for DX/DY/DS) or mrad (for DPHI/DTHETA/DPSI)
        """
        misalignment_type = misalignment_type.upper()
        assert misalignment_type in ["DX", "DY", "DS", "DPHI", "DTHETA", "DPSI"], f"Invalid type {misalignment_type}"

        parts = self._match_central_magnet_only(base_name)
        twiss_lookup = self._get_twiss_lookup()

        for name in parts:
            match = self.error_df["NAME"] == name
            if match.any():
                self.error_df.loc[match, misalignment_type] = float(value_mm)
            else:
                try:
                    s_val = float(twiss_lookup.get(name))
                except Exception:
                    s_val = float("inf")

                columns = ["NAME", "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI", "S"]
                row_full = {col: 0.0 for col in columns}
                row_full.update({
                    "NAME": name,
                    misalignment_type: float(value_mm),
                    "S": s_val
                })
                self.error_df = pd.concat([self.error_df, pd.DataFrame([row_full])], ignore_index=True)

    def add_quadrupole_misalignment(self, base_name, misalignment_type, value_mm):
        """
        Add a misalignment value (in mm or mrad) for the central magnet element only.

        Parameters:
            - base_name: base name of magnet (e.g. "sp6_qf")
            - misalignment_type: one of "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI"
            - value_mm: float, Magnitude in mm (for DX/DY/DS) or mrad (for DPHI/DTHETA/DPSI)
        """
        misalignment_type = misalignment_type.upper()
        assert misalignment_type in ["DX", "DY", "DS", "DPHI", "DTHETA", "DPSI"], f"Invalid type {misalignment_type}"

        parts = self.match_main_magnet_parts(base_name)
        twiss_lookup = self._get_twiss_lookup()

        for name in parts:
            match = self.error_df["NAME"] == name
            if match.any():
                self.error_df.loc[match, misalignment_type] = float(value_mm)
            else:
                try:
                    s_val = float(twiss_lookup.get(name))
                except Exception:
                    s_val = float("inf")

                columns = ["NAME", "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI", "S"]
                row_full = {col: 0.0 for col in columns}
                row_full.update({
                    "NAME": name,
                    misalignment_type: float(value_mm),
                    "S": s_val
                })
                self.error_df = pd.concat([self.error_df, pd.DataFrame([row_full])], ignore_index=True)

    def save_to_tfs(self, filename, origin=None):
        now = datetime.now()
        date_str = now.strftime("%d/%m/%y")
        time_str = now.strftime("%H.%M.%S")
        if origin is None:
            origin = "cpymad"

        header_lines = [
            '@ NAME             %06s "EFIELD"',
            '@ TYPE             %06s "EFIELD"',
            '@ TITLE            %08s "no-title"',
            f'@ ORIGIN           %16s "{origin}"',
            f'@ DATE             %08s "{date_str}"',
            f'@ TIME             %08s "{time_str}"',
        ]

        col_names = []
        for i in range(21):
            col_names.append(f'K{i}L')
            col_names.append(f'K{i}SL')
        col_names += [
            'DX', 'DY', 'DS', 'DPHI', 'DTHETA', 'DPSI',
            'MREX', 'MREY', 'MREDX', 'MREDY', 'AREX', 'AREY',
            'MSCALX', 'MSCALY', 'RFM_FREQ', 'RFM_HARMON', 'RFM_LAG'
        ]
        for i in range(21):
            col_names.append(f'P{i}L')
            col_names.append(f'P{i}SL')

        col_headers = "* NAME                        " + " ".join(f"{col:<12}" for col in col_names)
        col_types = "$ %s                          " + " ".join("%le".rjust(12) for _ in col_names)

        df = self.error_df.copy().sort_values("S", na_position="last")
        if "S" in df.columns:
            df = df.drop(columns=["S"])

        with open(filename, "w") as f:
            for line in header_lines:
                f.write(line + "\n")
            f.write(col_headers + "\n")
            f.write(col_types + "\n")

            for _, row in df.iterrows():
                name = f'"{row["NAME"].upper()}"'
                line = f' {name:<28}'

                for col in col_names:
                    val = row.get(col, 0.0)
                    try:
                        num = float(val)
                    except Exception:
                        num = 0.0

                    if col in ['DX', 'DY', 'DS']:
                        num *= 1e-3
                    elif col in ['DPHI', 'DTHETA', 'DPSI']:
                        num *= 1e-3

                    line += f"{num:.12f} "

                f.write(line.rstrip() + "\n")

    def add_vertical_misalignments_from_dataframe(self, df):
        """
        Process a DataFrame with survey data and apply DY and DPHI misalignments.

        Expected columns:
        - magnet, S_start, S_end, S_centre, angle, offset_centre

        Rules:
        - 'Dipole #' -> calls add_dipole_misalignment("sp#_dip", ...)
        - 'QD #'     -> calls add_quadrupole_misalignment("sp#_qd", ...)
        - 'QF #'     -> calls add_quadrupole_misalignment("sp#_qf", ...)
        - 'QC #'     -> calls add_quadrupole_misalignment("sp#_qds", ...)
        """
        def map_name_and_type(magnet):
            if magnet.startswith("Dipole "):
                return f"sp{magnet.split()[-1]}_dip", "dipole"
            elif magnet.startswith("QD "):
                return f"sp{magnet.split()[-1]}_qd", "quad"
            elif magnet.startswith("QF "):
                return f"sp{magnet.split()[-1]}_qf", "quad"
            elif magnet.startswith("QC "):
                return f"sp{magnet.split()[-1]}_qds", "quad"
            else:
                raise ValueError(f"Unrecognised magnet label: {magnet}")

        df = df.copy()
        df[["name", "type"]] = df["magnet"].apply(lambda m: pd.Series(map_name_and_type(m)))

        for _, row in df.iterrows():
            dy = row["offset_centre"]
            dphi = row["angle"]
            if row["type"] == "dipole":
                self.add_dipole_misalignment(row["name"], "DY", dy)
                self.add_dipole_misalignment(row["name"], "DPHI", -dphi)
            else:
                self.add_quadrupole_misalignment(row["name"], "DY", dy)
                self.add_quadrupole_misalignment(row["name"], "DPHI", -dphi)

    def add_horizontal_misalignments_from_dataframe(self, df):
        """
        Process a DataFrame with survey data and apply DX and DTHETA misalignments.

        Expected columns:
        - magnet, S_start, S_end, S_centre, angle, offset_corrected

        Rules:
        - 'Dipole #' -> calls add_dipole_misalignment("sp#_dip", ...)
        - 'QD #'     -> calls add_quadrupole_misalignment("sp#_qd", ...)
        - 'QF #'     -> calls add_quadrupole_misalignment("sp#_qf", ...)
        - 'QC #'     -> calls add_quadrupole_misalignment("sp#_qds", ...)
        """
        def map_name_and_type(magnet):
            if magnet.startswith("Dipole "):
                return f"sp{magnet.split()[-1]}_dip", "dipole"
            elif magnet.startswith("QD "):
                return f"sp{magnet.split()[-1]}_qd", "quad"
            elif magnet.startswith("QF "):
                return f"sp{magnet.split()[-1]}_qf", "quad"
            elif magnet.startswith("QC "):
                return f"sp{magnet.split()[-1]}_qds", "quad"
            else:
                raise ValueError(f"Unrecognised magnet label: {magnet}")

        df = df.copy()
        df[["name", "type"]] = df["magnet"].apply(lambda m: pd.Series(map_name_and_type(m)))

        for _, row in df.iterrows():
            dx = row["offset_corrected"]
            dtheta = row["angle"]
            if row["type"] == "dipole":
                self.add_dipole_misalignment(row["name"], "DX", dx)
                self.add_dipole_misalignment(row["name"], "DTHETA", dtheta)
            else:
                self.add_quadrupole_misalignment(row["name"], "DX", dx)
                self.add_quadrupole_misalignment(row["name"], "DTHETA", dtheta)


# Tune matching and harmonic tune functions below:

def get_isis_tunes(madx_instance, cpymad_logfile):
    initial_tunes = {}
    
    # Helper function for repeated tasks
    def get_tune_data(madx, cpymad_log, seq_name, is_ptc=False):
        if is_ptc:
            cpymad_ptc_twiss(madx, cpymad_log, seq_name)
            qx, qy = madx.table.ptc_twiss_summary.q1[0], madx.table.ptc_twiss_summary.q2[0]
        else:
            cpymad_madx_twiss(madx, cpymad_log, seq_name)
            qx, qy = madx.table.summ.q1[0], madx.table.summ.q2[0]
            
        return qx, qy
    
    # Fetch tunes from superperiod matrix phase advance
    qx, qy = superperiod_matrix_phase_advance(madx_instance, cpymad_logfile, madx_instance.sequence.superperiod)
    initial_tunes['SP Matrix Qx'], initial_tunes['SP Matrix Qy'] = round(qx * 10, 4), round(qy * 10, 4)
    
    # Define list of operation settings for looping
    operations = [
        ('superperiod', 'MAD-X sp Qx', 'MAD-X sp Qy', False),
        ('superperiod', 'PTC sp Qx', 'PTC sp Qy', True),
        ('synchrotron', 'MAD-X Qx', 'MAD-X Qy', False),
        ('synchrotron', 'PTC Qx', 'PTC Qy', True)
    ]

    # Get tune data for each operation
    for seq_name, qx_key, qy_key, is_ptc in operations:
        qx, qy = get_tune_data(madx_instance, cpymad_logfile, seq_name, is_ptc)
        
        if 'sp' in qx_key:
            qx *= 10
            qy *= 10
        elif qx_key == 'PTC Qx':
            qx += 4
            qy += 3
        
        # Rounding the values to 4 decimal places
        initial_tunes[qx_key], initial_tunes[qy_key] = round(qx, 4), round(qy, 4)
    
    return initial_tunes

def match_tune(madx_instance, sequence, requested_q1, requested_q2, requested_dq1=None, requested_dq2=None):
    
    madx_instance.command.match(chrom=True)

    madx_instance.command.vary(name='kqtd', step=1E-4)
    madx_instance.command.vary(name='kqtf', step=1E-4)

    if requested_dq1 is None or requested_dq2 is None:
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
    
from cpymad.types import Constraint
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
    
def tune_to_trim_quad_current_di(Qh=4.331, Qv=3.731,
                                 baseQh=4.331, baseQv=3.731, pn=1.0,
                                 z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03])):
    '''
    Calculates the trim quad currents required to obtain the given tune under
    Di's tune control system. The default values for the base tunes and
    normalised momentum are those at t = 0.0 ms.

    Parameters
    ----------
    Qh : Float, optional
        Horizontal tune required.
        Can be a float array of tunes. Must be of same length as Qv.
        The default is 4.331.
    Qv : Float, optional
        Vertical tune required.
        Can be a float array of tunes. Must be of same length as Qh
        The default is 3.731.
    baseQh : Float, optional
        Base horizontal tune.
        Can be a float array of tunes. Must be of same length as Qh, Qv.
        The default is 4.331.
    baseQv : Float, optional
        Base vertical tune.
        Can be a float array of tunes. Must be of same length as Qh, Qv.
        The default is 3.731.
    pn : Float, optional
        Normalised momentum.
        Can be a float array of momenta. Must be of same length as Qh, Qv.
        The default is 1.0.
    z : Float array of length 4, optional
        Coefficients of the tune control system.
        The default is np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03]).

    Returns
    -------
    I_F : Float, or array of floats with length same as Qh, Qv
        Current to apply to QTF to obtain Qh, Qv according to Di controls.
    I_D : Float, or array of floats with length same as Qh, Qv
        Current to apply to QTD to obtain Qh, Qv according to Di controls.

    '''

    # Get the control coefficients
    z1, z2, z3, z4 = z

    # Calculate the change in tune required
    dQh = Qh - baseQh
    dQv = Qv - baseQv

    # Calculate the currents required
    I_F = pn * (z1 * dQv - z3 * dQh) / (z1 * z4 - z2 * z3)
    I_D = pn * (z4 * dQh - z2 * dQv) / (z1 * z4 - z2 * z3)

    # Return the currents
    return -I_F, -I_D
    
def current_to_strength(I, Gcal=1.997e-3, Brho=1.23, pn=1.0):
    k = I * Gcal / Brho / pn
    return k
    
def tune_di_df(Qh=4.331, Qv=3.731, baseQh=4.331, baseQv=3.731, time_array=None, E_Max=800, z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03])):
    if time_array is None: 
        print('tune_to_trim_quad_current_di_df::error: time_array is None')
        return False    
    else:
        pn_array = return_pn(time_array)
        df_KE = synchrotron_kinetic_energy_df(E_Max, time_array)
        
    Iqtf_array = []
    Iqtd_array = []
    Kqtf_array = []
    Kqtd_array = []
    
    for pn, (_, ke_row) in zip(pn_array, df_KE.iterrows()):
        Iqtf, Iqtd = tune_to_trim_quad_current_di(Qh, Qv, baseQh, baseQv, pn, z)
        Iqtf_array.append(Iqtf)
        Iqtd_array.append(Iqtd)
        Kqtf_array.append(current_to_strength(Iqtf, Gcal=1.997e-3, Brho=ke_row['Rigidity [Tm]'], pn=pn))
        Kqtd_array.append(current_to_strength(Iqtd, Gcal=1.997e-3, Brho=ke_row['Rigidity [Tm]'], pn=pn))
        
    data = {
        'time': time_array,
        'pn': pn_array,
        'I_qtf': Iqtf_array,
        'I_qtd': Iqtd_array,
        'K_qtf': Kqtf_array,
        'K_qtd': Kqtd_array
    }
    
    df = pd.DataFrame(data)
    return df

import numpy as np
import matplotlib.pyplot as plt

def isis_plot_harmonic_tune_expectation(D7SIN=0.0, D7COS=0.0, D8SIN=0.0, D8COS=0.0, F8SIN=0.0, F8COS=0.0, F9SIN=0.0, F9COS=0.0):
    # Constants
    TQGCAL = 1.997E-3  # TQ gradient calibration per amp in T m-1 A-1
    HQCAL = 1.25  # Scaling factor for HQ for HQ=1 (multiply by value -10 -> +10 for D7cos etc)
    Brho = 1.23  # Injection value - change for later

    # Amplitudes and Harmonics for QTD
    da1 = D7COS * HQCAL * TQGCAL / Brho
    db1 = D7SIN * HQCAL * TQGCAL / Brho
    dhn1 = 7

    da2 = D8COS * HQCAL * TQGCAL / Brho
    db2 = D8SIN * HQCAL * TQGCAL / Brho
    dhn2 = 8

    # Amplitude and Harmonic for QTF
    fa1 = F8COS * HQCAL * TQGCAL / Brho
    fb1 = F8SIN * HQCAL * TQGCAL / Brho
    fhn1 = 8

    fa2 = F9COS * HQCAL * TQGCAL / Brho  # Doesn't exist in control system, but useful for future/models
    fb2 = F9SIN * HQCAL * TQGCAL / Brho
    fhn2 = 9

    # HER_qtd calculations
    HER_qtd = [
        (da1 * np.cos(dhn1 * (i / 10) * 2 * np.pi)) + (db1 * np.sin(dhn1 * (i / 10) * 2 * np.pi)) +
        (da2 * np.cos(dhn2 * (i / 10) * 2 * np.pi)) + (db2 * np.sin(dhn2 * (i / 10) * 2 * np.pi))
        for i in range(10)
    ]

    # HER_qtf calculations
    HER_qtf = [
        (fa1 * np.cos(fhn1 * (i / 10) * 2 * np.pi)) + (fb1 * np.sin(fhn1 * (i / 10) * 2 * np.pi)) +
        (fa2 * np.cos(fhn2 * (i / 10) * 2 * np.pi)) + (fb2 * np.sin(fhn2 * (i / 10) * 2 * np.pi))
        for i in range(10)
    ]

    # Generate dynamic title based on non-zero variables
    title_parts = []
    if D7SIN != 0.0:
        title_parts.append(f"D7SIN={D7SIN}")
    if D7COS != 0.0:
        title_parts.append(f"D7COS={D7COS}")
    if D8SIN != 0.0:
        title_parts.append(f"D8SIN={D8SIN}")
    if D8COS != 0.0:
        title_parts.append(f"D8COS={D8COS}")
    if F8SIN != 0.0:
        title_parts.append(f"F8SIN={F8SIN}")
    if F8COS != 0.0:
        title_parts.append(f"F8COS={F8COS}")
    if F9SIN != 0.0:
        title_parts.append(f"F9SIN={F9SIN}")
    if F9COS != 0.0:
        title_parts.append(f"F9COS={F9COS}")

    # Combine the title parts
    plot_title = "Harmonic Tunes: " + ", ".join(title_parts)

    # HER_qtd and HER_qtf values combined for the plot
    HER_values = HER_qtd + HER_qtf

    # Labels for each variable
    HER_labels = [
        f"R{i}_QTD" for i in range(10)
    ] + [
        f"R{i}_QTF" for i in range(10)
    ]

    # Indices for the bar plot
    indices = np.arange(len(HER_values))

    # Create the bar plot
    plt.figure(figsize=(12, 6))
    plt.bar(indices, HER_values, color=['blue' if 'qtd' in label else 'orange' for label in HER_labels])

    # Add labels, title, and legend
    plt.xlabel("Trim Quad")
    plt.ylabel(r"Quadrupole Strength $K$ [$m^{-2}$]")
    plt.title(plot_title)
    plt.xticks(indices, HER_labels, rotation=45, ha="right")
    plt.tight_layout()
    plt.grid(which='both', ls=':', lw=0.5, c='grey')

    # Display the plot
    plt.show()
    
    
def isis_print_harmonic_tunes(madx_instance):
    """
    Print the current harmonic tune values from the MAD-X instance.

    Parameters:
        madx_instance: Instance of MAD-X to retrieve harmonic tune values from.
    """
    # Retrieve and print values for all harmonics
    harmonics = ["D7SIN", "D7COS", "D8SIN", "D8COS", "F8SIN", "F8COS", "F9SIN", "F9COS"]

    for harmonic in harmonics:
        if harmonic in madx_instance.globals.defs:
            value = madx_instance.globals.defs[harmonic]
            print(f"{harmonic} = {value}")
        else:
            print(f"{harmonic} is not defined in the MAD-X instance.")
            

def isis_set_harmonic_tune(madx_instance, D7SIN=0.0, D7COS=0.0, D8SIN=0.0, D8COS=0.0, F8SIN=0.0, F8COS=0.0, F9SIN=0.0, F9COS=0.0):
    """
    Set harmonic tunes in a MAD-X instance with the provided values.
    
    Parameters:
        madx_instance: Instance of MAD-X to send commands to.
        D7SIN, D7COS, D8SIN, D8COS, F8SIN, F8COS, F9SIN, F9COS: Harmonic tune values.
    """
    # Set commands for all harmonics
    harmonics = {
        "D7SIN": D7SIN,
        "D7COS": D7COS,
        "D8SIN": D8SIN,
        "D8COS": D8COS,
        "F8SIN": F8SIN,
        "F8COS": F8COS,
        "F9SIN": F9SIN,
        "F9COS": F9COS,
    }

    for harmonic, value in harmonics.items():
        command = f"{harmonic} = {value};"
        madx_instance.input(command)
