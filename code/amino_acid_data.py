"""
Source data and some conversion functions for amino acids. 
The content is adapted from the localCIDER package. 
http://pappulab.github.io/localCIDER/
https://github.com/Pappulab/localCIDER/blob/master/localcider/backend/data/aminoacids.py
"""
#.............................................................................#

# One letter code to full name

ONE_TO_FULL = {'A': 'Alanine',
               'C': 'Cysteine',
               'D': 'Aspartic acid',
               'E': 'Glutamic acid',
               'F': 'Phenylalanine',
               'G': 'Glycine',
               'H': 'Histidine',
               'I': 'Isoleucine',
               'K': 'Lysine',
               'L': 'Leucine',
               'M': 'Methionine',
               'N': 'Asparagine',
               'P': 'Proline',
               'Q': 'Glutamine',
               'R': 'Arginine',
               'S': 'Serine',
               'T': 'Threonine',
               'V': 'Valine',
               'W': 'Tryptophan',
               'Y': 'Tyrosine'}

#.............................................................................#

# One letter code to three letters code

ONE_TO_THREE = {'A': 'ALA',
                'C': 'CYS',
                'D': 'ASP',
                'E': 'GLU',
                'F': 'PHE',
                'G': 'GLY',
                'H': 'HIS',
                'I': 'ILE',
                'K': 'LYS',
                'L': 'LEU',
                'M': 'MET',
                'N': 'ASN',
                'P': 'PRO',
                'Q': 'GLN',
                'R': 'ARG',
                'S': 'SER',
                'T': 'THR',
                'V': 'VAL',
                'W': 'TRP',
                'Y': 'TYR'}

#.............................................................................#

# Three letters code to one letter code

THREE_TO_ONE = {'ALA': 'A',
                'CYS': 'C',
                'ASP': 'D',
                'GLU': 'E',
                'PHE': 'F',
                'GLY': 'G',
                'HIS': 'H',
                'ILE': 'I',
                'LYS': 'K',
                'LEU': 'L',
                'MET': 'M',
                'ASN': 'N',
                'PRO': 'P',
                'GLN': 'Q',
                'ARG': 'R',
                'SER': 'S',
                'THR': 'T',
                'VAL': 'V',
                'TRP': 'W',
                'TYR': 'Y'}

# Kyte–Doolittle hydrophobicity values of the amino acids
# https://web.expasy.org/protscale/pscale/Hydropath.Doolittle.html

KD_HYDROPHOBICITY = {'ILE': 4.5,
                     'VAL': 4.2,
                     'LEU': 3.8,
                     'PHE': 2.8,
                     'CYS': 2.5,
                     'MET': 1.9,
                     'ALA': 1.8,
                     'GLY': -0.4,
                     'THR': -0.7,
                     'SER': -0.8,
                     'TRP': -0.9,
                     'TYR': -1.3,
                     'PRO': -1.6,
                     'HIS': -3.2,
                     'GLU': -3.5,
                     'GLN': -3.5,
                     'ASP': -3.5,
                     'ASN': -3.5,
                     'LYS': -3.9,
                     'ARG': -4.5}
#.............................................................................#

# Charges of amino acid residues

RESIDUE_CHARGE = {'ILE': 0,
                  'VAL': 0,
                  'LEU': 0,
                  'PHE': 0,
                  'CYS': 0,
                  'MET': 0,
                  'ALA': 0,
                  'GLY': 0,
                  'THR': 0,
                  'SER': 0,
                  'TRP': 0,
                  'TYR': 0,
                  'PRO': 0,
                  'HIS': 1,
                  'GLU': -1,
                  'GLN': 0,
                  'ASP': -1,
                  'ASN': 0,
                  'LYS': 1,
                  'ARG': 1}
#.............................................................................#

# Twenty standard amino acids

TWENTY_AAs = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q',
              'C', 'G', 'P', 'A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']

#.............................................................................#


def one_to_three(one_let):
    return ONE_TO_THREE[one_let]


def one_to_full(one_let):
    return ONE_TO_FULL[one_let]


def three_to_one(three_let):
    return THREE_TO_ONE[three_let]


def get_residue_hydrophobicity(one_let):
    """
    Takes the one letter code of an amino acid and returns its hydrophobicity in
    the Kyte–Doolittle scale.
    A simple method for displaying the hydropathic character of a protein. 
    Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32.
    """
    three_let = one_to_three(one_let)
    return KD_HYDROPHOBICITY[three_let]


def get_residue_charge(one_let):
    """
    Takes the one letter code of an amino acid and returns its side-chain charge.
    """
    three_let = one_to_three(one_let)
    return RESIDUE_CHARGE[three_let]


def get_CGed_hyd_chrg(one_let, cut_off=0.0):
    """
    Takes a one letter code of an amino acid and returns its coarse-grained
    hydrophobicity-charge string. 
    Ref:
    Schilling, C., Mack, T., Lickfett, S., Sieste, S., Ruggeri, F. S., 
    Sneideris, T., Dutta, A., Bereau, T., Naraghi, R., Sinske, D., Knowles,
    T. P. J., Synatschke, C. V., Weil, T., Knöll, B., Sequence‐Optimized 
    Peptide Nanofibers as Growth Stimulators for Regeneration of Peripheral 
    Neurons. Adv. Funct. Mater. 2019, 29, 1809112. 
    """
    hyd_val = get_residue_hydrophobicity(one_let)
    chrg_val = get_residue_charge(one_let)

    if (hyd_val > cut_off):
        hyd_type = 'P'
    else:
        hyd_type = 'N'

    if(chrg_val > 0.0):
        chrg_type = 'P'
    elif(chrg_val == 0.0):
        chrg_type = 'Z'
    else:
        chrg_type = 'N'

    return hyd_type+chrg_type

# print(get_CGed_hyd_chrg('H'))