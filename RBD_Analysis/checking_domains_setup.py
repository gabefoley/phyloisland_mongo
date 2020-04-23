import os

# Specify the base directory
base_dir = "/Users/gabefoley/Dropbox/Code/Python_Workspace/phyloisland_2019/RBD_Analysis"

# Specify the A region file
filename = 'first20'

tag_dict = {}

# Specify the outgroup
outgroup = 'NZ_AMBZ01000025.1_information_Paenibacillus_alvei_region_TcdA1_expanded_1274632_1281550_backward'

# Specify the profile name and load the profile_dict
profile_name = 'leidreiter_profiles'

# Define the colours to use for each type - list of colours are here - http://etetoolkit.org/docs/latest/reference/reference_treeview.html
colour_dict = {'type_1' :'dodgerblue', 'type_2B': 'gold', 'type_2A':'green', 'type_3' :'purple', 'various' : 'red', 'unknown' : 'black'}