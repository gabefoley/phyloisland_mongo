from Bio.PDB import *

# # Instance variables
# def __init__(self, id = str(), path = '.') :
#     self.id = id
#     self.path = path
#     self.lines = list()

# methods : downloading and reading flat files (PDB file or CSV)
def download_pdb(self, pdb_id, output_dir = '.') :
    ''' Download a PDB file with Biopython PDBList class. Returns the donwloaded
    file path.
    /!\ the Biopython fonction assings the format name : 'pdb<pdb_id>.ent'
    '''
    pdb_file = PDBList()
    pdb_file.retrieve_pdb_file(pdb_id, pdir = output_dir, file_format = 'pdb')
    file_name = "pdb"+pdb_id.lower()+".ent"
    self.id = pdb_id
    self.path = output_dir + file_name

def read_file(self, path = '.') :
    ''' Read a flat file. Assigns a lines list to a lines attribute. This
    fonction is used by CSV and PDB files.
    '''
    if path != '.' :
        self.path = path
    f = open(self.path, "r")
    lines = f.readlines()
    f.close()
    self.lines = lines

def split_PDBfile_by_chains(self, output_dir = '.', chains = 'all', all_sections = True ) :
    ''' Split a pdb file in different pdb files by chains. data is a list of
    pdb file lines. chains must be a list of PDB ids (e.g. ['A', 'B'])
    '''
    pdblines = self.lines
    # file split :
    initial_sections = list()
    dict_chains = dict()
    final_sections = list()
    i = 0
    while i < len(pdblines) :
        line = pdblines[i]
        if line[0:4] != 'ATOM' and line[0:3] != 'TER' :
            initial_sections.append(line)
            i += 1
        else :
            break
    while i < len(pdblines) :
        line = pdblines[i]
        possible_sections = ['ATOM  ', 'ANISOU', 'TER   ', 'HETATM']
        if line[0:6]in possible_sections:
            chain_id = line[21]
            if not(chain_id in dict_chains) :
                dict_chains[chain_id] = [line]
            else :
                dict_chains[chain_id].append(line)
            i += 1
        else :
            break
    while i < len(pdblines) :
        line = pdblines[i]
        final_sections.append(line)
        i += 1

    # Chains selection :
    if chains == 'all' :
        chains_id_list = dict_chains.keys()
        print('esto va si all en split', dict_chains.keys())
    else :
        chains_id_list = sorted(chains)
    pdb_id = self.id
    self.id = list()
    self.path = list()

    # Write the different files
    for chain_id in chains_id_list :
        sub_file_id = pdb_id +  '_' + chain_id
        sub_file_name = 'pdb' + sub_file_id + '.ent'
        sub_file_path = output_dir + sub_file_name
        f = open(sub_file_path, 'w')
        if all_sections :
            f.writelines(initial_sections)
        f.writelines(dict_chains[chain_id])
        if all_sections :
            f.writelines(final_sections)
        f.close()
        self.id.append((pdb_id, chain_id))
        self.path.append(sub_file_path)

