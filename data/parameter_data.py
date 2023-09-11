# Copied from https://github.com/sringe/catmap-1/blob/master/catmap/data/parameter_data.py

#Define default ideal gas parameters
#first is symmetry number (how many possible rotations (not reflections) lead to 
#  indistinguishable conformation)
#second is linear or nonlinear
#third is spin
ideal_gas_params = {
        'H2_g':[2,'linear',0],
        'N2_g':[2,'linear',0],
        'O2_g':[2,'linear',1.0],
        'H2O_g':[2,'nonlinear',0],
        'CO_g':[1,'linear',0],
        'CH4_g':[12,'nonlinear',0],
        'NH3_g':[3,'nonlinear',0],
        'N2O_g':[2,'nonlinear',0],
        'NH2OH_g':[1,'nonlinear',0],
        'NO_g':[2,'linear',0],
        'N2_g':[2,'linear',0],
        'CH3OH_g':[1,'nonlinear',0],
        'CH4O2_g':[2,'nonlinear',0], #methanediol
        'CH3CH2OH_g':[1,'nonlinear',0],
        'CO2_g':[2,'linear',0],
        'CH2O_g':[2,'nonlinear',0],
        'HCOOH_g':[1,'nonlinear',0],
        'CH2CH2_g':[4,'nonlinear',0],
        'CH3CHCH2_g':[1,'nonlinear',0], #propene
        'CH3CH2CHCH2_g':[1,'nonlinear',0], #1-butene
        'CH3CHCHCH3_g':[2,'nonlinear',0], #2-butene, ok for both trans and cis
        'CH3CH3CCH2_g':[2,'nonlinear',0], #isobutene
        'pe_g':[2,'linear',0], # fictitious gas molecule corresponding to proton electron pair
        'C2H2_g':[2,'linear',0],
        'C2H4_g':[4,'nonlinear',0],
        'C2H6_g':[6,'nonlinear',0],
        'C3H6_g':[1,'nonlinear',0],
        'CH3COOH_g':[1,'nonlinear',0],
        'CH3CHO_g':[1,'nonlinear',0],
        'C5H4O2_g':[1,'nonlinear',0], #Added by N. Shan, KSU from NIST/CRC
        'C5H6O2_g':[1,'nonlinear',0], #Added by N. Shan, KSU from NIST/CRC
        'C5H6O_g':[1,'nonlinear',0], #Added by N. Shan, KSU from NIST/CRC
        'HCl_g':[1,'linear',0], #Added by M. Andersen, TUM from NIST
        'Cl2_g':[2,'linear',0], #Added by M. Andersen, TUM from NIST
}

def generate_hbond_dict():
    """Returns a dictionary of generic, surface-agnostic hydrogen bond stabilizations, in eV."""
    OH = -0.50  # OH directly on the surface
    ROH = -0.25  # a 'floppy' OH group
    CO = -0.1  # carbon monoxide
    d = {'COOH'   : ROH,
         'OCHO'   : 0.,
         'CO'     : CO,
         'CHO'    : CO,
         'CH2O'   : 0.,
         'CO2'    : 0.,
         'OCH3'   : 0.,
         'O'      : 0.,
         'OH'     : OH,
         'H'      : 0.,
         'COH'    : ROH,
         'C'      : 0.,
         'CH'     : 0.,
         'CH2'    : 0.,
         'CH3'    : 0.,
         'CHOH'   : ROH,
         'COHOH'  : ROH,
         'OCH2O'  : 0.,
         'CH2OH'  : ROH,
         'OCHCH2' : 0.,
         'OCHCHO' : 0.,
        }
    return d

hbond_dict = generate_hbond_dict()

#Atmos. Chem. Phys., 15, 4399–4981, 2015
henry_consts = {  # mol/m^3/Pa
    "CH4"       : 1.4e-5,
    "C2H6"      : 1.9e-5,
    "CH3OH"     : 2.0,
    "CH3CH2OH"  : 1.9,
    "CO"        : 9.7e-6,
    "CO2"       : 3.3e-4,
    "N2"        : 6.4e-6,
    "H2"        : 7.8e-6,
    "NH3"       : 5.9e-1,
    "O2"        : 1.2e-5,
    "CH2O"      : 3.2e1,
    "NO"        : 1.9e-5,
}