from qol.tools.integrate import trapz_cond
import qol.mesa.const as const

import h5py
import os
import warnings

import numpy as np
from astropy.table import Table, Column
from scipy.interpolate import interp1d

# MESA to GYRE input version formats; colnames
version_number_colname_dict = \
        {1: ['k', 'r', 'M_r/(M_*-M_r)', 'L_r', 'P', 'T', 'rho', # Version 0.01
             'Grad', 'N2', 'c_V', 'c_P', 'chi_T', 'chi_rho',
             'kappa', 'kappa_T', 'kappa_rho',
             'eps', 'eps_n eps_n,T', 'eps_n eps_n,rho'],
        19: ['k', 'r', 'M_r/(M_*-M_r)', 'L_r', 'P', 'T', 'rho', # Version 0.19
             'Grad', 'N2', 'Gamma_1', 'Grad_ad', 'nu_T',
             'kappa', 'kappa_T', 'kappa_rho',
             'eps', 'eps_n eps_n,T', 'eps_n eps_n,rho',
             'Omega_rot'],
        100:['k', 'r', 'M_r','L_r', 'P', 'T', 'rho',
             'Grad', 'N2', 'Gamma_1', 'Grad_ad', 'nu_T', # Version 1.00
             'kappa', 'kappa kappa_T', 'kappa kappa_rho',
             'eps', 'eps_n eps_n,T', 'eps_n eps_n,rho',
             'Omega_rot'],
        101:['k', 'r', 'M_r', 'L_r', 'P', 'T', 'rho', # Version 1.01
             'Grad', 'N2', 'Gamma_1', 'Grad_ad', 'nu_T',
             'kappa', 'kappa kappa_T', 'kappa kappa_rho',
             'eps_n', 'eps_n eps_n,T', 'eps_n,rho',
             'Omega_rot'],
        120:['k', 'r', 'M_r', 'L_r', 'P', 'T', 'rho', # Version 1.20
             'Grad', 'N2', 'Gamma_1', 'Grad_ad', 'nu_T',
             'kappa', 'kappa kappa_T', 'kappa kappa_rho',
             'eps_n', 'eps_n eps_n,T', 'eps_n,rho', 'eps_grav',
             'Omega_rot'],
        }


class GyreRun:
    """
    Reads the contents of a GYRE directory and implements some helpful methods
    """
    def __init__(self, gyre_path, input_profile_rel_path=None, summary_file='summary.h5'):
        """
        for now, assumes all files are h5. TODO also accommodate txt

        input_profile_rel_path = relative path to input MESA profile in GYRE format
        """
        # store gyre_path
        self.gyre_path = gyre_path
        
        assert os.path.exists(os.path.join(gyre_path, summary_file))
        self.summary_file = summary_file

        # Read summary file
        self.summary = self.read_h5(summary_file, cache_key=None)

        # If given, read input file
        self.input_profile_rel_path = input_profile_rel_path
        if input_profile_rel_path is not None:
            self.read_input_profile(input_profile_rel_path)

        self.detail_cache_dict = {}
    
    def read_input_profile(self, rel_path):
        input_profile_abs_path = os.path.join(self.gyre_path, self.input_profile_rel_path)
        assert os.path.exists(input_profile_abs_path)

        with open(input_profile_abs_path, 'r') as f:
            lines = f.readlines()

        header_items = lines[0].split()
        self.header_N = header_N = int(header_items[0])
        self.header_Mstar = header_Mstar = float(header_items[1])
        self.header_Rstar = header_Rstar = float(header_items[2])
        self.header_Lstar = header_Lstar = float(header_items[3])

        # Determine version number and grab relevant column names
        version_number = int(header_items[4]) if len(header_items) == 5 else 1
        colnames = version_number_colname_dict[version_number]

        # Grab data and form a table
        coldatas = [np.array([]) for _ in colnames]

        for line in lines[1:-1]: # exclude first (header) and last (blank) lines
            items = line.split()
            
            assert len(colnames) == len(items)
            for jj, item in enumerate(items):
                if jj == 0: # if colname is the first one, 'k'
                    item = int(item)
                else:
                    item = float(item)

                coldatas[jj] = np.append(coldatas[jj], item)

        cols = [Column(coldata, name=colname) for coldata, colname in zip(coldatas, colnames)]
        data = Table(cols)

        self.input_profile = data

    def read_detail(self, l, n_pg):
        rel_path = f'detail.l{l}.n{n_pg:+d}.h5'
        cache_key = f'l{l}.n{n_pg:+d}'

        data = self.read_h5(rel_path=rel_path, cache_key=cache_key)

        return data

    def read_h5(self, rel_path, cache_key=None):
        """
        rel_path is relative to gyre_path -- read and store data as astropy.table.Table

        if cache_key is not None:
        - if key exists, just read from detail_cache_dict
        - if not, store it in detail_cache_dict with cache_as as key
        """
        # check if cached
        if (cache_key is not None) and (cache_key in self.detail_cache_dict.keys()):
            data = self.detail_cache_dict[cache_key]
            return data

        abs_path = os.path.join(self.gyre_path}, {rel_path})
        assert os.path.exists(abs_path)

        cols = []

        with h5py.File(abs_path, 'r') as f:
            for key in f.keys():
                col = Column(f[key][:], name=key)
                cols.append(col)
        
        data = Table(cols)

        if cache_key is not None:
            self.detail_cache_dict[cache_key] = data
        
        return data
    
    def get_ω(self, l, n_pg, freq_units='UHZ'):
        """
        Get an angular frequency for a given mode
        """
        assert ((self.summary['l'] == l) & (self.summary['n_pg'] == n_pg)).any()
        freq_tuple = self.summary['freq'][(self.summary['l'] == l) & (self.summary['n_pg'] == n_pg)][0]
        freq = freq_tuple[0] + 1j * freq_tuple[1]
        freq = np.real(freq) # assume real (ok for adiabatic)

        match freq_units: # convert to angular frequency
            case 'UHZ':
                ω = 2. * const.pi * 1e-6 * freq
            case _:
                raise NotImplementedError()
        
        return ω

    def get_magnetic_scriptI_nonasympt(self, l, n_pg, freq_units='UHZ'):
        """
        calculate scriptI using Equation B1 in Rui, Fuller, Hermes 2025 (ApJ)
        """
        # Need input file to get dimensionful radius
        assert self.input_profile_rel_path is not None

        data = self.read_detail(l=l, n_pg=n_pg)
        assert 'x' in data.colnames
        assert 'xi_h' in data.colnames
        assert 'prop_type' in data.colnames

        # Interpolate grid to get rho
        r_model = self.input_profile['r']
        rho_model = self.input_profile['rho'] # rho on model grid
        rho_interp = interp1d(r_model, rho_model, bounds_error=False, fill_value=(rho_model[0], rho_model[-1]))

        r = self.header_Rstar * data['x']
        rho = rho_interp(r) # rho on wavefunction grid

        # Grab frequency
        ω = self.get_ω(l=l, n_pg=n_pg, freq_units=freq_units)

        # Calculate integrals
        warnings.warn('magnetic_scriptI_nonasympt currently uses a kludge which assumes imaginary part of xi_h is zero. Fine for adiabatic. TODO: generalize this')
        xi_h = np.array([xi_h_tuple[0] + 1j * xi_h_tuple[1] for xi_h_tuple in data['xi_h']])
        xi_h = np.real(xi_h)
        is_g = (data['prop_type'] == -1) # propagating gravity wave

        numer_integral = trapz_cond(x=r, y=np.gradient(r * xi_h, r) ** 2, c=is_g)
        denom_integral = trapz_cond(x=r, y=xi_h ** 2 * rho * r ** 2, c=is_g)

        magnetic_scriptI_nonasympt = ω ** 2 * numer_integral / denom_integral / l / (l + 1)

        return magnetic_scriptI_nonasympt

    def get_magnetic_scriptI_nonasympt_grid(self, l, freq_units='UHZ'):
        """
        run magnetic_scriptI_nonasympt for all modes of a given l
        """
        good = (self.summary['l'] == l) & (self.summary['n_pg'] < 0) # assert correct l and that it needs to have some g-mode part
        l_arr = np.array(self.summary['l'][good])
        n_pg_arr = np.array(self.summary['n_pg'][good])

        ω_arr = np.array([self.get_ω(l, n_pg, freq_units=freq_units) for n_pg in n_pg_arr])
        magnetic_scriptI_nonasympt_arr = np.array([self.get_magnetic_scriptI_nonasympt(l, n_pg, freq_units=freq_units) for n_pg in n_pg_arr])

        return ω_arr, magnetic_scriptI_nonasympt_arr

    def get_Brshift_nonasympt_RFH25(self, νB_uHz, l, n_pg, freq_units='UHZ'):
        ω = self.get_ω(l=l, n_pg=n_pg, freq_units=freq_units)
        P = 2 * const.pi / ω

        νB = 1e-6 * νB_uHz
        magnetic_scriptI_nonasympt = self.get_magnetic_scriptI_nonasympt(l=l, n_pg=n_pg, freq_units=freq_units)

        Al = l * (l + 1) / 2
        Brshift_nonasympt = 8 * const.pi ** 2.5 * np.sqrt(νB / (Al * magnetic_scriptI_nonasympt * P ** 3))

        return Brshift_nonasympt

    def get_Brshift_nonasympt_RFH25_grid(self, νB_uHz, l, freq_units='UHZ'):
        """
        run Brshift_nonasympt for all modes of a given l
        """
        good = (self.summary['l'] == l) & (self.summary['n_pg'] < 0)
        l_arr = np.array(self.summary['l'][good])
        n_pg_arr = np.array(self.summary['n_pg'][good])

        ω_arr = np.array([self.get_ω(l, n_pg, freq_units=freq_units) for n_pg in n_pg_arr])
        Brshift_nonasympt_arr = np.array([self.get_Brshift_nonasympt_RFH25(νB_uHz=νB_uHz, l=l, n_pg=n_pg, freq_units=freq_units) for n_pg in n_pg_arr])

        return ω_arr, Brshift_nonasympt_arr
