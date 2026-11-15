"""AUTO-GENERATED from reference/src/MGN2MECH/INCLDIR/SPC_NOCONVER.EXT.

Do not edit by hand; regenerate with tools/parse_spca.py. Holds MEGAN's
201-species speciation table: names, molar weights (g/mol), and the map from
each of the 201 speciated species to its parent MEGAN class (0-based index)
with the speciation fraction EFFS_ALL.
"""

N_SPCA = 201

SPCA_NAMES = [
    'isoprene', 'MBO_2m3e2ol', 'pinene_a', 'pinene_b', 'myrcene', 'ocimene_al', 'ocimene_c_b',
    'ocimene_t_b', 'bornene', 'camphene', 'carene_3', 'fenchene_a', 'fenchene_b', 'tricyclene',
    'limonene', 'phellandrene_b', 'phellandrene_a', 'sabinene', 'terpinene_a', 'terpinene_g',
    'terpinolene', 'thujene_a', 'verbenene', 'cymene_o', 'cymene_p', 'meta-cymenene',
    'p-cymenene', 'borneol', 'bornyl_ACT', 'camphor', 'estragole', 'fenchone', 'ionone_b',
    'ipsenol', 'myrtenal', 'octanol', 'octenol_1e3ol', 'piperitone', 'terpineol_4',
    'terpineol_a', 'thujone_a', 'thujone_b', 'acoradiene', 'bergamotene_a', 'bergamotene_b',
    'bisabolene_a', 'bisabolene_b', 'cadinene_d', 'cadinene_g', 'calarene', 'caryophyllene_b',
    'caryophyllene_c_i', 'elemene_b', 'farnescene_a', 'farnescene_b', 'farnesol',
    'germacrene_B', 'germacrene_D', 'humulene_a', 'humulene_g', 'muurolene_a', 'muurolene_g',
    'nerolidol_c', 'nerolidol_t', 'patchoulene_b', 'selinene_b', 'selinene_d',
    '2ethyl_hexyl_salate', 'aromadendrene', 'bourbonene_b', 'cadinol_a', 'cedrene_a', 'cedrol',
    'copaene_a', 'cubebene_a', 'cubebene_b', 'gurjunene_b', 'homosalate', 'isolongifolene',
    'kaur_16_ene', 'longifolene', 'longipinene', 'thujopsene', 'ylangene_a', 'zingiberene_a',
    'copaene_b', 'longicyclene', 'methanol', 'acetone', 'acetaldehyde', 'ethanol',
    'acetic_acid', 'formic_acid', 'pyruvic_acid', '2_methylprop_1_ene', 'butane', 'butene',
    'butene_c2', 'butene_t2', 'ethane', 'ethene', 'propane', 'propene', 'decanal',
    'geranyl_acetone', 'met_heptenone', 'neryl_acetone', 'nonanal', 'nonenal', 'octanal',
    'pentanal', 't-2-octenal', 'formaldehyde', 'oxopentanal', '1_pentanol', '1_penten_3_one',
    '2_met_butanenitrile', '2_pentene', '3_metfuran', '3_metthiophene', '3_pentanone',
    '3_penten_2_ol', '4_penten_1_ol', '8_heptadecene', 'acetophenone', 'anisole',
    'benzaldehyde', 'benzene', 'benzyl_benzoate', 'benzyl_acetate', 'benzyl_alcohol',
    'chavicol', 'cineole_1_8', 'cinnamaldehyde', 'cinnamic_acid', 'coniferyl_alcohol',
    '2met_nonatriene', 'ethyl_cinnamate', 'eugenol', 'guaiacol', 'hexanal', 'hexanol_1',
    'hexenal_c3', 'hexenal_t2', 'hexenol_c3', 'hexenyl_ACT_c3', 'hydrogen_cyanide',
    'isobutyronitrile', 'jasmone', 'linalool', 'linalool_OXD_c', 'linalool_OXD_t',
    'methacrolein', 'met_benzoate', 'met_jasmonate', 'met_salicylate', 'met_acetate',
    'met_vinylketone', 'p_coumaric_acid', 'phenol', 'phenylacetaldehyde', 'safrole',
    'salicylaldehyde', 'toluene', '3met_3DCTT', 'xylene', 'z2_penten_1ol', 'z2_penten_3ol',
    'z3_hexen_1yl_butyrate', 'butanone_2', 'carbon_2s', 'carbonyl_s', 'decane', 'diallyl_2s',
    '2met_2s', '2met_s', 'dodecene_1', 'heptanal', 'heptane', 'heptanone', 'hexane',
    'hydrogen_s', 'indole', 'methane', 'met_bromide', 'met_chloride', 'met_iodide',
    'met_mercaptan', 'met_propenyl_2s', 'MBO_3m2e1ol', 'MBO_3m3e1ol', 'naphthalene', 'pentane',
    'PPPP_2s', 'santene', 'skatole', 'terpinyl_ACT_a', 'tetradecene_1', 'undecane',
    'carbon_monoxide', 'nitric_OXD'
]

# molar weight (g/mol) per SPCA species
SPCA_MWT = [
    68.12, 84.12, 136.24, 136.24, 136.24, 136.24, 136.24, 136.24, 136.24, 136.24, 136.24,
    136.24, 136.24, 136.24, 136.24, 136.24, 136.24, 136.24, 136.24, 136.24, 136.24, 136.24,
    134.22, 134.22, 134.22, 132.21, 132.21, 154.25, 196.29, 152.24, 148.2, 152.24, 192.3,
    154.25, 150.22, 128.21, 128.21, 152.24, 154.25, 154.25, 152.24, 152.24, 204.36, 204.36,
    204.36, 204.36, 204.36, 204.36, 204.36, 204.36, 204.36, 204.36, 204.36, 204.36, 204.36,
    222.37, 204.36, 204.36, 204.36, 204.36, 204.36, 204.36, 222.37, 222.37, 204.36, 204.36,
    204.36, 250.34, 204.36, 204.36, 222.37, 204.36, 222.37, 204.36, 204.36, 204.36, 204.36,
    262.35, 204.36, 272.47, 204.36, 204.36, 204.36, 204.36, 204.36, 204.36, 204.36, 32.04,
    58.08, 44.05, 46.07, 60.05, 46.03, 88.06, 56.11, 58.12, 56.11, 56.11, 56.11, 30.07, 28.05,
    44.1, 42.08, 156.27, 194.32, 126.2, 194.32, 142.24, 140.23, 128.21, 86.13, 126.2, 30.03,
    100.12, 88.15, 84.12, 83.13, 70.13, 82.1, 98.17, 86.13, 86.13, 86.13, 238.46, 120.15,
    108.14, 106.12, 78.11, 212.25, 150.18, 108.14, 134.18, 154.25, 132.16, 148.16, 180.2,
    150.26, 176.22, 164.2, 124.14, 100.16, 102.18, 98.14, 98.14, 100.16, 142.2, 27.03, 69.11,
    164.25, 154.25, 170.25, 170.25, 70.09, 136.15, 224.3, 152.15, 74.08, 70.09, 164.16, 94.11,
    120.15, 162.19, 122.12, 92.14, 218.38, 106.17, 86.13, 86.13, 170.25, 72.11, 76.14, 60.08,
    142.28, 146.28, 94.2, 62.14, 168.32, 114.19, 100.2, 114.19, 86.18, 34.08, 117.15, 16.04,
    94.94, 50.49, 141.94, 48.11, 120.24, 86.13, 86.13, 128.17, 72.15, 148.29, 122.21, 131.18,
    196.29, 196.38, 156.31, 28.01, 30.01
]

# 0-based parent MEGAN class (0..18) for each of the 201 species
MAP_MG20 = [
    0, 1, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 18, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    10, 11, 12, 12, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 17, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
    17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 7
]

# 0-based SPCA species index each map entry writes to
MAP_SPCA = [
    0, 1, 2, 4, 5, 6, 7, 3, 9, 8, 11, 10, 12, 13, 21, 17, 39, 38, 196, 16, 19, 18, 14, 15, 20,
    22, 24, 23, 25, 26, 29, 28, 37, 30, 27, 41, 40, 149, 194, 135, 132, 34, 31, 32, 33, 199, 50,
    51, 47, 47, 46, 53, 64, 52, 62, 55, 58, 60, 43, 56, 65, 49, 54, 42, 59, 61, 57, 48, 63, 44,
    45, 72, 77, 67, 71, 82, 80, 70, 84, 78, 86, 73, 69, 81, 75, 83, 74, 85, 79, 76, 68, 87, 88,
    90, 89, 92, 91, 93, 99, 100, 101, 95, 96, 102, 94, 97, 98, 112, 110, 177, 105, 103, 109,
    107, 113, 108, 111, 104, 106, 140, 139, 138, 160, 125, 137, 133, 134, 148, 131, 162, 158,
    150, 151, 136, 164, 114, 127, 146, 147, 152, 157, 156, 153, 161, 120, 124, 130, 126, 163,
    159, 117, 141, 155, 128, 129, 143, 119, 122, 144, 118, 154, 165, 166, 115, 145, 121, 142,
    182, 123, 168, 116, 167, 183, 184, 185, 186, 187, 36, 35, 170, 175, 169, 195, 191, 192, 179,
    180, 176, 172, 178, 171, 189, 174, 190, 198, 197, 173, 193, 188, 181, 200
]

# speciation fraction of the parent class assigned to this species
MAP_EFFS = [
    1.0, 1.0, 1.0, 0.5, 0.01, 0.09, 0.4, 0.33, 0.18, 0.005, 0.005, 0.2, 0.005, 0.02, 0.04, 0.14,
    0.05, 0.02, 0.005, 0.06, 0.09, 0.07, 0.6, 0.1, 0.075, 0.005, 0.2, 0.07, 0.02, 0.05, 0.08,
    0.02, 0.01, 0.02, 0.05, 0.02, 0.03, 0.2, 0.005, 0.005, 0.2, 0.005, 0.005, 0.005, 0.005, 1.0,
    0.5, 0.005, 0.005, 0.005, 0.02, 0.05, 0.005, 0.01, 0.005, 0.005, 0.1, 0.005, 0.04, 0.005,
    0.08, 0.005, 0.05, 0.01, 0.01, 0.03, 0.02, 0.005, 0.005, 0.02, 0.005, 0.01, 0.15, 0.08,
    0.03, 0.01, 0.15, 0.01, 0.01, 0.01, 0.01, 0.2, 0.01, 0.03, 0.02, 0.01, 0.05, 0.03, 0.15,
    0.01, 0.02, 1.0, 1.0, 0.5, 0.5, 0.8, 0.19, 0.01, 0.18, 0.36, 0.12, 0.005, 0.03, 0.28, 0.01,
    0.01, 0.005, 0.9, 0.005, 0.005, 0.01, 0.005, 0.005, 0.005, 0.03, 0.005, 0.005, 0.02, 0.005,
    0.06, 0.005, 0.01, 0.005, 0.01, 0.005, 0.005, 0.005, 0.005, 0.02, 0.005, 0.005, 0.005,
    0.005, 0.05, 0.01, 0.005, 0.005, 0.05, 0.005, 0.05, 0.05, 0.01, 0.01, 0.005, 0.005, 0.01,
    0.01, 0.01, 0.01, 0.01, 0.005, 0.02, 0.09, 0.005, 0.005, 0.05, 0.005, 0.005, 0.1, 0.01,
    0.01, 0.005, 0.005, 0.005, 0.09, 0.005, 0.1, 0.01, 0.005, 0.005, 0.005, 0.005, 0.15, 0.03,
    0.15, 0.22, 0.01, 0.01, 0.01, 0.01, 0.1, 0.07, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
    0.01, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05, 1.0
]
