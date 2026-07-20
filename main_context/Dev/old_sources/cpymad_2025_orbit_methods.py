from helper_functions import *
from cpymad_helpers import *

import os
import shutil
import numpy as np
import pandas as pd

from cpymad.madx import Madx        
from cpymad.types import Constraint

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
from matplotlib.lines import Line2D

#from cpymad.madx import Constraint

all_kicker_names = [
    # Steering correctors (horizontal)
    #'r0hd1_kick', 'r2hd1_kick', 'r3hd1_kick', 'r4hd1_kick', 'r5hd1_kick', 'r7hd1_kick', 'r9hd1_kick',
    
    # Steering correctors (vertical)
    #'r0vd1_kick', 'r2vd1_kick', 'r3vd1_kick', 'r4vd1_kick', 'r5vd1_kick', 'r7vd1_kick', 'r9vd1_kick',
    
    # Dipole kickers (horizontal)
    'R0DIP_HKICK', 'R1DIP_HKICK', 'R2DIP_HKICK', 'R3DIP_HKICK', 'R4DIP_HKICK',
    'R5DIP_HKICK', 'R6DIP_HKICK', 'R7DIP_HKICK', 'R8DIP_HKICK', 'R9DIP_HKICK',
    
    # Dipole kickers (vertical)
    'R0DIP_VKICK', 'R1DIP_VKICK', 'R2DIP_VKICK', 'R3DIP_VKICK', 'R4DIP_VKICK',
    'R5DIP_VKICK', 'R6DIP_VKICK', 'R7DIP_VKICK', 'R8DIP_VKICK', 'R9DIP_VKICK',
    
    # Quad Focus (QF) kickers (horizontal)
    'R0QF_HKICK', 'R1QF_HKICK', 'R2QF_HKICK', 'R3QF_HKICK', 'R4QF_HKICK',
    'R5QF_HKICK', 'R6QF_HKICK', 'R7QF_HKICK', 'R8QF_HKICK', 'R9QF_HKICK',
    
    # Quad Focus (QF) kickers (vertical)
    'R0QF_VKICK', 'R1QF_VKICK', 'R2QF_VKICK', 'R3QF_VKICK', 'R4QF_VKICK',
    'R5QF_VKICK', 'R6QF_VKICK', 'R7QF_VKICK', 'R8QF_VKICK', 'R9QF_VKICK',
    
    # Quad Defocus (QD) kickers (horizontal)
    'R0QD_HKICK', 'R1QD_HKICK', 'R2QD_HKICK', 'R3QD_HKICK', 'R4QD_HKICK',
    'R5QD_HKICK', 'R6QD_HKICK', 'R7QD_HKICK', 'R8QD_HKICK', 'R9QD_HKICK',
    
    # Quad Defocus (QD) kickers (vertical)
    'R0QD_VKICK', 'R1QD_VKICK', 'R2QD_VKICK', 'R3QD_VKICK', 'R4QD_VKICK',
    'R5QD_VKICK', 'R6QD_VKICK', 'R7QD_VKICK', 'R8QD_VKICK', 'R9QD_VKICK',
    
    # Trim Quads QTF (horizontal)
    'R0QTF_HKICK', 'R1QTF_HKICK', 'R2QTF_HKICK', 'R3QTF_HKICK', 'R4QTF_HKICK',
    'R5QTF_HKICK', 'R6QTF_HKICK', 'R7QTF_HKICK', 'R8QTF_HKICK', 'R9QTF_HKICK',
    
    # Trim Quads QTF (vertical)
    'R0QTF_VKICK', 'R1QTF_VKICK', 'R2QTF_VKICK', 'R3QTF_VKICK', 'R4QTF_VKICK',
    'R5QTF_VKICK', 'R6QTF_VKICK', 'R7QTF_VKICK', 'R8QTF_VKICK', 'R9QTF_VKICK',
    
    # Trim Quads QTD (horizontal)
    'R0QTD_HKICK', 'R1QTD_HKICK', 'R2QTD_HKICK', 'R3QTD_HKICK', 'R4QTD_HKICK',
    'R5QTD_HKICK', 'R6QTD_HKICK', 'R7QTD_HKICK', 'R8QTD_HKICK', 'R9QTD_HKICK',
    
    # Trim Quads QTD (vertical)
    'R0QTD_VKICK', 'R1QTD_VKICK', 'R2QTD_VKICK', 'R3QTD_VKICK', 'R4QTD_VKICK',
    'R5QTD_VKICK', 'R6QTD_VKICK', 'R7QTD_VKICK', 'R8QTD_VKICK', 'R9QTD_VKICK',
    
    # Singlet Quad (QDS) kickers (horizontal)
    'R0QDS_HKICK', 'R1QDS_HKICK', 'R2QDS_HKICK', 'R3QDS_HKICK', 'R4QDS_HKICK',
    'R5QDS_HKICK', 'R6QDS_HKICK', 'R7QDS_HKICK', 'R8QDS_HKICK', 'R9QDS_HKICK',
    
    # Singlet Quad (QDS) kickers (vertical)
    'R0QDS_VKICK', 'R1QDS_VKICK', 'R2QDS_VKICK', 'R3QDS_VKICK', 'R4QDS_VKICK',
    'R5QDS_VKICK', 'R6QDS_VKICK', 'R7QDS_VKICK', 'R8QDS_VKICK', 'R9QDS_VKICK'
]

def isis_reset_steering(madx_instance):
    kicker_groups = [
        # Steering correctors
        ['r0hd1_kick', 'r2hd1_kick', 'r3hd1_kick', 'r4hd1_kick', 'r5hd1_kick', 'r7hd1_kick', 'r9hd1_kick',
         'r0vd1_kick', 'r2vd1_kick', 'r3vd1_kick', 'r4vd1_kick', 'r5vd1_kick', 'r7vd1_kick', 'r9vd1_kick'],
        
        # QTF trim quad kickers
        [f'R{i}QTF_HKICK' for i in range(10)] + [f'R{i}QTF_VKICK' for i in range(10)],
        
        # QTD trim quad kickers
        [f'R{i}QTD_HKICK' for i in range(10)] + [f'R{i}QTD_VKICK' for i in range(10)],
        
        # QF doublet quad kickers
        [f'R{i}QF_HKICK' for i in range(10)] + [f'R{i}QF_VKICK' for i in range(10)],
        
        # QD doublet quad kickers
        [f'R{i}QD_HKICK' for i in range(10)] + [f'R{i}QD_VKICK' for i in range(10)],
        
        # QDS singlet quad kickers
        [f'R{i}QDS_HKICK' for i in range(10)] + [f'R{i}QDS_VKICK' for i in range(10)],
        
        # Main dipole kickers
        [f'R{i}DIP_HKICK' for i in range(10)] + [f'R{i}DIP_VKICK' for i in range(10)],
    ]

    # Set all to zero
    for group in kicker_groups:
        for kicker in group:
            madx_instance.globals[kicker] = 0.0

def match_co(madx_instance, constraints, sp_list=None,
             horizontal=True, steering=True, dip=False,
             doublet_quad=False, trim_quad=False, singlet_quad=False,
             manual_vary_list=None, step_size=1e-4):
    """
    Setup orbit correction matching in MAD-X using corrector, dipole and quadrupole kicks.

    If manual_vary_list is given, will vary only those variables.
    """

    if sp_list is None:
        sp_list = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

    madx_instance.command.match(chrom=False)

    if manual_vary_list is not None:
        for kicker in manual_vary_list:
            madx_instance.command.vary(name=kicker, step=step_size)
    else:
        # Regular auto selection by plane and type
        vary_elements = {
            'steering': {'H': 'r{sp}hd1_kick', 'V': 'r{sp}vd1_kick'},
            'dip': {'H': 'R{sp}DIP_HKICK', 'V': 'R{sp}DIP_VKICK'},
            'trim_quad': {'H': ['R{sp}QTF_HKICK', 'R{sp}QTD_HKICK'], 'V': ['R{sp}QTF_VKICK', 'R{sp}QTD_VKICK']},
            'doublet_quad': {'H': ['R{sp}QF_HKICK', 'R{sp}QD_HKICK'], 'V': ['R{sp}QF_VKICK', 'R{sp}QD_VKICK']},
            'singlet_quad': {'H': 'R{sp}QDS_HKICK', 'V': 'R{sp}QDS_VKICK'}
        }

        knob_flags = {
            'steering': steering,
            'dip': dip,
            'trim_quad': trim_quad,
            'doublet_quad': doublet_quad,
            'singlet_quad': singlet_quad
        }

        plane = 'H' if horizontal else 'V'

        for sp in sp_list:
            for knob_type, active in knob_flags.items():
                if active:
                    element_names = vary_elements[knob_type][plane]
                    if isinstance(element_names, list):
                        for elem in element_names:
                            madx_instance.command.vary(name=elem.format(sp=sp), step=step_size)
                    else:
                        madx_instance.command.vary(name=element_names.format(sp=sp), step=step_size)

    for c in constraints:
        madx_instance.command.constraint(**c)

    madx_instance.command.jacobian(calls=50000, tolerance=1e-6)
    madx_instance.command.endmatch()


def brute_force_find_best_kicks(madx_instance, constraints, all_kicker_names, N,
                                 horizontal=True, step_size=1e-4):
    """
    Brute-force search for the N most likely kickers matching orbit measurements.

    Parameters:
    - madx_instance: cpymad MAD-X instance
    - constraints: list of constraint dictionaries
    - all_kicker_names: list of all kicker knob names available (e.g., ['r0hd1_kick', 'R0QF_HKICK', ...])
    - N: number of kickers to select
    - horizontal: whether matching horizontal plane (affects knob templates)
    - step_size: vary step size (default 1e-4)

    Returns:
    - best_combo: list of best kicker names
    - best_residual: residual norm for best solution
    """

    best_residual = np.inf
    best_combo = None

    for combo in itertools.combinations(all_kicker_names, N):
        try:
            # Reset sequence if needed
            madx_instance.input("use, period=synchrotron;")

            # Match with only these kickers varied
            match_co(
                madx_instance=madx_instance,
                constraints=constraints,
                manual_vary_list=list(combo),
                step_size=step_size
            )

            # Measure residual
            residual = calculate_constraint_residual(madx_instance, constraints)

            if residual < best_residual:
                best_residual = residual
                best_combo = combo

        except Exception as e:
            print(f"Skipping combo {combo} due to error: {e}")
            continue

    return best_combo, best_residual


def calculate_constraint_residual(madx_instance, constraints):
    """
    Calculates the total residual from all constraint violations.

    Very simple: sum of squared deviations from targets.
    """
    res = 0.0
    for c in constraints:
        range_name = c['range']
        if 'x' in c:
            model_value = madx_instance.table.twiss.x[madx_instance.table.twiss.name.index(range_name)]
            res += (model_value - c['x'])**2
        if 'y' in c:
            model_value = madx_instance.table.twiss.y[madx_instance.table.twiss.name.index(range_name)]
            res += (model_value - c['y'])**2
    return np.sqrt(res)
    
class BPMFitResultsLoader:
    def __init__(self, filepath, reverse_co=False):
        self.filepath = filepath
        self.reverse_co = reverse_co
        self.df = self._load_file()
        
    def _load_file(self):
        return pd.read_csv(self.filepath, sep=r"\s+")
    
    def get_bpm_list(self):
        return self.df['bpm'].unique().tolist()

    def get_available_parameters(self):
        return self.df.columns

    def get_parameter(self, param_name, bpm=None):
        if param_name not in self.df.columns:
            raise ValueError(f"Parameter '{param_name}' not found in data.")
        if bpm:
            return self.df[self.df['bpm'] == bpm][param_name].values
        return self.df[param_name].values

    def get_plane_parameters(self, plane="H", parameters=None):
        if parameters is None:
            parameters = self.get_available_parameters()
        missing = [p for p in parameters if p not in self.df.columns]
        if missing:
            raise ValueError(f"Missing parameters in file: {missing}")
        df_plane = self.df[self.df['plane'] == plane].copy()
        return df_plane[["bpm"] + parameters]

    def get_co_data(self, plane="H", twiss_df=None):
        required_columns = ["closed_orbit_mm", "closed_orbit_mm_err"]
        missing = [col for col in required_columns if col not in self.df.columns]
        if missing:
            raise ValueError(f"Missing required CO data columns: {missing}")
        
        df_plane = self.df[self.df['plane'] == plane].copy()
    
        if self.reverse_co:
            df_plane["closed_orbit_mm"] *= -1
    
        # Optional Twiss 's' lookup by partial BPM name
        if twiss_df is not None:
            twiss_lookup = twiss_df.copy()
            twiss_lookup["name_stripped"] = twiss_lookup["name"].astype(str).str.lower()
            bpm_to_s = {}
    
            for bpm in df_plane["bpm"]:
                bpm_key = bpm.lower()
                matched = twiss_lookup[twiss_lookup["name_stripped"].str.contains(bpm_key)]
                if len(matched) == 1:
                    bpm_to_s[bpm] = matched["s"].values[0]
                else:
                    bpm_to_s[bpm] = float("nan")  # ambiguous or not found
    
            df_plane["s"] = df_plane["bpm"].map(bpm_to_s)
    
        return df_plane[["bpm", "closed_orbit_mm", "closed_orbit_mm_err"] + (["s"] if "s" in df_plane.columns else [])]

    def get_dataframe(self):
        return self.df

def constraints_from_co_data(df_h=None, df_v=None, err=None):
    constraints_list = []
    
    err = err or 0.1  # If None, default to 0.1

    if df_h is not None:
        for _, row in df_h.iterrows():
            bpm = row['bpm']
            value_m = row['closed_orbit_mm'] * 1e-3
            error_m = max(row['closed_orbit_mm_err'], err) * 1e-3
            constraints_list.append(dict(
                range=bpm,
                x=Constraint(min=value_m - error_m, max=value_m + error_m)
            ))

    if df_v is not None:
        for _, row in df_v.iterrows():
            bpm = row['bpm']
            value_m = row['closed_orbit_mm'] * 1e-3
            error_m = max(row['closed_orbit_mm_err'], err) * 1e-3
            constraints_list.append(dict(
                range=bpm,
                y=Constraint(min=value_m - error_m, max=value_m + error_m)
            ))

    return constraints_list
    
def check_constraints_against_twiss(twiss_df, constraints):
    """
    Fix constraint names by matching them to the real names in the twiss table,
    removing any ':1' or similar from names to be valid in MAD-X constraints.

    Args:
        twiss_df (DataFrame): MAD-X twiss() output as a pandas DataFrame.
        constraints (list): List of constraint dictionaries.

    Returns:
        List of corrected constraint dictionaries.
    """
    corrected_constraints = []
    available_names = twiss_df['name'].str.lower().values

    for constraint in constraints:
        target_range = constraint.get('range', '').lower()

        # Find matching element in twiss names
        matching_names = [name for name in available_names if target_range in name]

        if not matching_names:
            print(f"Warning: No match found for constraint '{target_range}'. Constraint will be skipped.")
            continue

        if len(matching_names) > 1:
            print(f"Warning: Multiple matches found for constraint '{target_range}': {matching_names}. Using first match.")

        matched_name = matching_names[0]

        # Remove colon part, e.g. sp2_r2vm2:1 -> sp2_r2vm2
        matched_name_clean = matched_name.split(':')[0]

        # Build new constraint with corrected range
        corrected_constraint = constraint.copy()
        corrected_constraint['range'] = matched_name_clean
        corrected_constraints.append(corrected_constraint)

    return corrected_constraints
    
def print_constraints_readable(constraints, mm=True):
    """
    Nicely print a list of constraints for human readability.

    Parameters:
    -----------
    constraints : list
        List of constraint dictionaries.
    mm : bool, optional
        If True, prints values in mm instead of m.
    """
    unit = 'mm' if mm else 'm'
    scale = 1e3 if mm else 1.0

    for c in constraints:
        var_type = 'x' if 'x' in c else 'y'
        constraint_obj = c[var_type]
        min_val = constraint_obj.min * scale
        max_val = constraint_obj.max * scale
        print(f"Range: {c['range']:20}  Type: {var_type}  Min: {min_val:.6g} {unit}  Max: {max_val:.6g} {unit}")


def cpymad_get_nonzero_kicks(madx):
    kick_globals = {}
    for name, value in madx.globals.defs.items():
        if 'kick' in name.lower() and abs(value) > 0.0:
            value_mrad = value * 1e3  # convert from rad to mrad
            formatted_value = float(f"{value_mrad:.3g}")  # 3 significant figures
            kick_globals[name] = formatted_value
    return kick_globals


class TwissClosedOrbitPlotter:
    def __init__(self, sequence_name=None, ptc_twiss=False):
        self.sequence_name = sequence_name
        self.ptc_twiss = ptc_twiss
        self.twiss_list = []  # (df_myTwiss, label)
        self.bpm_data_list = []  # (df_bpm_data, label, plane)

    def add_twiss(self, df_myTwiss, label=""):
        self.twiss_list.append((df_myTwiss, label))

    def add_data(self, df_bpm_data, label="", plane="H"):
        self.bpm_data_list.append((df_bpm_data.copy(), label, plane.upper()))

    def plot(self, save_file=None, xlimits=None, ylimits=None, heights=[1, 3, 2, 2], tunes=None):
        if not self.twiss_list:
            raise ValueError("No Twiss data provided. Use add_twiss() first.")
    
        # Validate all Twiss
        ref_len = len(self.twiss_list[0][0])
        ref_start = self.twiss_list[0][0]['name'].iloc[0]
        for df, lbl in self.twiss_list:
            if len(df) != ref_len or df['name'].iloc[0] != ref_start:
                raise ValueError(f"TWISS mismatch: {lbl} has incompatible length or start point")
    
        df_ref, _ = self.twiss_list[0]
    
        if self.ptc_twiss:
            gamma_rel = df_ref.headers['GAMMA']
            beta_rel = np.sqrt(1 - (1 / gamma_rel**2))
            qx = df_ref.headers['Q1']
            qy = df_ref.headers['Q2']
        else:
            beta_rel = 1.0
            if tunes is not None:
                qx, qy = tunes
            else:
                # Try find q1, q2 in df_ref
                if 'q1' in df_ref.columns and 'q2' in df_ref.columns:
                    qx = df_ref['q1'].iloc[0]
                    qy = df_ref['q2'].iloc[0]
                elif 'Q1' in df_ref.columns and 'Q2' in df_ref.columns:
                    qx = df_ref['Q1'].iloc[0]
                    qy = df_ref['Q2'].iloc[0]
                else:
                    # If nothing, default to zero
                    qx = 0.0
                    qy = 0.0

    
        plot_title = f"{self.sequence_name} Q1={qx:.3f}, Q2={qy:.3f}" if self.sequence_name else f"Q1={qx:.3f}, Q2={qy:.3f}"
    
        fig = plt.figure(figsize=(16, 8), constrained_layout=True)
        spec = gridspec.GridSpec(nrows=4, ncols=1, figure=fig, height_ratios=heights)
    
        ax1 = fig.add_subplot(spec[0])
        self._block_diagram(ax1, df_ref, limits=xlimits)
        ax1.set_title(plot_title)
    
        ax2 = fig.add_subplot(spec[1], sharex=ax1)
        ax2.plot(df_ref['s'], df_ref['betx'], 'b', label=r'$\beta_x$')
        ax2.plot(df_ref['s'], df_ref['bety'], 'r', label=r'$\beta_y$')
        ax2.set_ylabel(r'$\beta_{x,y}$ [m]')
        ax2.grid(which='both', ls=':', lw=0.5, color='k')
        ax2.legend(loc='center left', bbox_to_anchor=(1.05, 0.7))
    
        ax2r = ax2.twinx()
        if self.ptc_twiss:
            ax2r.plot(df_ref['s'], df_ref['disp1'] / beta_rel, 'g', label=r'$D_x$')
            ax2r.plot(df_ref['s'], df_ref['disp3'] / beta_rel, 'purple', label=r'$D_y$')
        else:
            ax2r.plot(df_ref['s'], df_ref['dx'], 'g', label=r'$D_x$')
            ax2r.plot(df_ref['s'], df_ref['dy'], 'purple', label=r'$D_y$')
        ax2r.set_ylabel(r'$D_{x,y}$ [m]')
        ax2r.legend(loc='center left', bbox_to_anchor=(1.15, 0.7))
        ax2r.grid(which='both', ls=':', lw=0.5, color='green')
    
        # Generate colormap
        import matplotlib.cm as cm
        cmap = plt.colormaps['jet'].resampled(len(self.twiss_list) + len(self.bpm_data_list))
        colors = [cmap(i) for i in range(len(self.twiss_list) + len(self.bpm_data_list))]
        twiss_colors = colors[:len(self.twiss_list)]
        bpm_colors = colors[len(self.twiss_list):]
    
        ax3 = fig.add_subplot(spec[2], sharex=ax1)
        x_vals_all = df_ref['x'].values * 1e3
        for i, (df_tw, label) in enumerate(self.twiss_list):
            ax3.plot(df_tw['s'], df_tw['x'] * 1e3, lw=1.5, label=label, color=twiss_colors[i])
            x_vals_all = np.concatenate([x_vals_all, df_tw['x'].values * 1e3])
    
        for i, (df, label, plane) in enumerate(self.bpm_data_list):
            if plane == 'H':
                s_vals = df['bpm'].apply(lambda bpm_name: self._find_s_for_bpm(df_ref, bpm_name))
                y_vals = df['closed_orbit_mm']
                y_errs = df['closed_orbit_mm_err']
                ax3.errorbar(s_vals, y_vals, yerr=y_errs, fmt='o', label=label, color=bpm_colors[i])
                x_vals_all = np.concatenate([x_vals_all, y_vals.values + y_errs, y_vals.values - y_errs])
        
        # Injection marker
        ax3.axvline(x=8.5, color='black', linestyle='--', lw=1.2)
        ax3.text(8.5, 0, 'Injection', ha='center', va='bottom', fontsize=8, color='black')

        ax3.set_ylabel('x [mm]')
        ax3.grid(which='both', ls=':', lw=0.5, color='k')
        ax3.legend(loc='center left', bbox_to_anchor=(1.01, 0.5))
    
        for i in range(10):
            x_pos = 16.336282 * i
            ax3.axvline(x=x_pos, color='gray', linestyle='--', lw=0.8)
            ax3.text(x_pos + 16.336282 / 2, 0, f'SP{i}',
                     ha='center', va='bottom', fontsize=8, color='gray')
    
        ax4 = fig.add_subplot(spec[3], sharex=ax1)
        y_vals_all = df_ref['y'].values * 1e3
        for i, (df_tw, label) in enumerate(self.twiss_list):
            ax4.plot(df_tw['s'], df_tw['y'] * 1e3, lw=1.5, label=label, color=twiss_colors[i])
            y_vals_all = np.concatenate([y_vals_all, df_tw['y'].values * 1e3])
    
        for i, (df, label, plane) in enumerate(self.bpm_data_list):
            if plane == 'V':
                s_vals = df['bpm'].apply(lambda bpm_name: self._find_s_for_bpm(df_ref, bpm_name))
                y_vals = df['closed_orbit_mm']
                y_errs = df['closed_orbit_mm_err']
                ax4.errorbar(s_vals, y_vals, yerr=y_errs, fmt='o', label=label, color=bpm_colors[i])
                y_vals_all = np.concatenate([y_vals_all, y_vals.values + y_errs, y_vals.values - y_errs])
        
        # Injection marker
        ax4.axvline(x=8.5, color='black', linestyle='--', lw=1.2)
        ax4.text(8.5, 0, 'Injection', ha='center', va='bottom', fontsize=8, color='black')

        ax4.set_ylabel('y [mm]')
        ax4.set_xlabel('s [m]')
        ax4.grid(which='both', ls=':', lw=0.5, color='k')
        ax4.legend(loc='center left', bbox_to_anchor=(1.01, 0.5))
    
        for i in range(10):
            x_pos = 16.336282 * i
            ax4.axvline(x=x_pos, color='gray', linestyle='--', lw=0.8)
            ax4.text(x_pos + 16.336282 / 2, 0, f'SP{i}',
                     ha='center', va='bottom', fontsize=8, color='gray')
    
        co_max = np.ceil(max(np.max(np.abs(x_vals_all)), np.max(np.abs(y_vals_all))) / 20) * 20
        co_min = -co_max
        ax3.set_ylim(co_min, co_max)
        ax4.set_ylim(co_min, co_max)
    
        if save_file:
            plt.savefig(save_file, bbox_inches='tight', dpi=200)
        plt.show()
       
    def _find_s_for_bpm(self, df_twiss, bpm_name):
        matches = df_twiss[df_twiss['name'].str.contains(bpm_name, case=False, na=False)]
        if not matches.empty:
            return matches.iloc[0]['s']
        return np.nan

    def _block_diagram(self, ax1, df_myTwiss, limits=None):
        ax1.axis('off')
        keyword = 'keyword'
        elem_defs = {
            'marker': ('k', 0.1),
            'kicker': ('c', 0.5), 'hkicker': ('c', 0.5), 'vkicker': ('c', 0.5),
            'sextupole': ('green', 'l'), 'quadrupole': ('r', 'l'),
            'sbend': ('b', 'l'), 'rbend': ('b', 'l')
        }
        for key, (color, length) in elem_defs.items():
            key_to_match = key.upper() if self.ptc_twiss else key
            if keyword not in df_myTwiss.columns:
                continue
            df_elem = df_myTwiss[df_myTwiss[keyword] == key_to_match]
            for _, elem in df_elem.iterrows():
                l = length if isinstance(length, (int, float)) else elem[length]
                x0 = elem['s'] - l if isinstance(length, str) else elem['s']
                ax1.add_patch(patches.Rectangle((x0, 0.), l, 1.0, color=color, alpha=0.5))

        if limits:
            ax1.set_xlim(limits[0], limits[1])
        else:
            L = df_myTwiss.headers['LENGTH'] if self.ptc_twiss else df_myTwiss.iloc[-1].s
            ax1.set_xlim(0, L)

        custom_lines = [
            Line2D([0], [0], color='b', lw=4, alpha=0.5),
            Line2D([0], [0], color='r', lw=4, alpha=0.5),
            Line2D([0], [0], color='green', lw=4, alpha=0.5),
            Line2D([0], [0], color='cyan', lw=4, alpha=0.5),
            Line2D([0], [0], color='k', lw=4, alpha=0.5)
        ]
        ax1.legend(custom_lines, ['Dipole', 'Quadrupole', 'Sextupole', 'Kicker', 'Marker'],
           loc='center left', bbox_to_anchor=(1.01, 0.5))

def cpymad_vertical_offset_error(madx_instance, cpymad_logfile, sequence_name, magnet_name, offset=0.0, clear=False):
    if clear:
        cpymad_clear_errors(madx_instance)
    madx_instance.input(f'''
        select, flag=error, pattern={magnet_name};
        ealign, dy := {offset} * 1e-3;
    ''')
    cpymad_helpers.cpymad_write_to_logfile(cpymad_logfile, f"Misalignment applied to {magnet_name}: dy = {offset} mm")


def cpymad_vertical_angle_error(madx_instance, cpymad_logfile, sequence_name, magnet_name, angle=0.0, clear=False):
    if clear:
        cpymad_clear_errors(madx_instance)
    madx_instance.input(f'''
        select, flag=error, pattern={magnet_name};
        ealign, dpsi := {angle} * 1e-3;
    ''')
    cpymad_helpers.cpymad_write_to_logfile(cpymad_logfile, f"Misalignment applied to {magnet_name}: dpsi = {angle} mrad")

def expand_vertical_errors_from_table(
    madx_instance,
    sequence_name,
    error_table,
    offset_col="offset_mm",
    angle_col="angle_mrad",
    magnet_col="magnet_name",
    s_distance_warning=1.0,
    check_s_distance=False
):
    """
    For each entry in the input error table, find all matching MAD-X elements,
    and return a new DataFrame with one row per matched element, including
    dy and dpsi values for each.

    Parameters:
    - madx_instance: cpymad MAD-X instance
    - sequence_name: name of the MAD-X sequence
    - error_table: input DataFrame with misalignment data
    - offset_col: name of column for vertical offsets (in mm)
    - angle_col: name of column for vertical angles (in mrad)
    - magnet_col: name of column for matching magnets
    - s_distance_warning: warn if matched elements are far apart
    - check_s_distance: if True, check and warn on large s spreads

    Returns:
    - A new DataFrame with columns: ['element_name', 's', 'dy_mm', 'dpsi_mrad', 'magnet_name']
    """
    madx_instance.use(sequence=sequence_name)
    tw_df = madx_instance.twiss().dframe()

    expanded_rows = []

    # Pre-match element names for each unique magnet
    name_map = {}
    for magnet_base in error_table[magnet_col].unique():
        # Strip any :1, :2, etc. when saving in expanded list
        matched = tw_df[tw_df['name'].str.contains(magnet_base, case=False, regex=False)]
        name_map[magnet_base] = matched

    for _, row in error_table.iterrows():
        magnet_base = row[magnet_col]
        offset = row[offset_col]
        angle = row[angle_col]

        matched = name_map.get(magnet_base, pd.DataFrame())
        if matched.empty:
            print(f"[Warning] No match found for {magnet_base}")
            continue

        if check_s_distance and len(matched) > 1:
            s_span = matched['s'].max() - matched['s'].min()
            if s_span > s_distance_warning:
                print(f"[Warning] Elements for {magnet_base} span {s_span:.3f} m in s.")

        for _, match_row in matched.iterrows():
            # Strip instance suffix like ":1" from element name
            clean_name = match_row["name"].split(":")[0]
            expanded_rows.append({
                "element_name": clean_name,
                "s": match_row["s"],
                "dy_mm": offset,
                "dpsi_mrad": angle,
                "magnet_name": magnet_base
            })

    return pd.DataFrame(expanded_rows)


def apply_expanded_vertical_errors(
    madx_instance,
    cpymad_logfile: str,
    sequence_name: str,
    expanded_df: pd.DataFrame,
    error_filename: str = "cpymad_error_file.tfs"
):
    """
    Apply vertical offset and angle errors from an expanded dataframe using existing helper functions.

    Parameters:
    - madx_instance: cpymad MAD-X instance
    - cpymad_logfile: Path to logfile
    - sequence_name: MAD-X sequence name
    - expanded_df: DataFrame with columns ['element_name', 'dy_mm', 'dpsi_mrad']
    - error_filename: Filename to save MAD-X error configuration
    """
    import cpymad_helpers  # assumed to be defined

    # Initialise error system to ADD mode
    cpymad_eoption(madx_instance, cpymad_logfile, add_errors=True)
    madx_instance.use(sequence=sequence_name)

    # Apply each error to each matching element
    for _, row in expanded_df.iterrows():
        name = row['element_name']
        if 'dy_mm' in row and pd.notna(row['dy_mm']):
            cpymad_vertical_offset_error(
                madx_instance, cpymad_logfile, sequence_name,
                magnet_name=name, offset=row['dy_mm'], clear=False
            )
        if 'dpsi_mrad' in row and pd.notna(row['dpsi_mrad']):
            cpymad_vertical_angle_error(
                madx_instance, cpymad_logfile, sequence_name,
                magnet_name=name, angle=row['dpsi_mrad'], clear=False
            )

    # Save error table
    cpymad_save_error_table(madx_instance, error_filename=error_filename)

def cpymad_eoption(madx_instance, cpymad_logfile, add_errors = True):
    seed = np.random.randint(0,1E6)
    if add_errors: 
        madx_instance.input(f"eoption, seed={seed}, add=True;")
        log_string = f"Errors set to add, i.e. not overwrite"
        cpymad_helpers.cpymad_write_to_logfile(cpymad_logfile, log_string)
    else:
        madx_instance.input(f"eoption, seed={seed}, add=False;")
        log_string = f"Errors set to overwrite, i.e. not sum"
        cpymad_helpers.cpymad_write_to_logfile(cpymad_logfile, log_string)    

def cpymad_save_error_table(madx_instance, error_filename=None):
    # Save the error configuration to a file (optional)
    if error_filename is None: error_filename = 'cpymad_error_file.tfs'
    madx_instance.input(f'esave, file="{error_filename}";')

def cpymad_set_correctors(madx_instance, cpymad_logfile, corrector_dict, max_E=800., time=0.0):
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
