import importlib.util
import pymatgen as pmg
import pandas as pd
from pymatgen.io.xyz import XYZ
import os
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor as pm
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.graphs import MoleculeGraph
import pymatgen.analysis.local_env as env
from pymatgen.analysis.dimensionality import get_structure_components, get_dimensionality_larsen
import time

import shutil

import openbabel
import pybel as pb
import numpy as np

#Molecule = importlib.util.spec_from_file_location("Molecule", "/home/rgur/tools/ob-2.4.1/lib/python3.7/site-packages/pybel.py")
from pathlib import Path
from pymatgen.io.babel import BabelMolAdaptor
import sys
from multiprocessing import Pool

args = sys.argv

parent_dir = args[1]
bond_df_path = args[2]
n_core = int(args[3])

try:
    past_df_path = args[4]
except:
    past_df_path = None

try:
    cifs_dir = args[5]
except:
    cifs_dir = parent_dir    
    
m_dict = {0: {"Zn": 4, "O": 1}, 1: {"Zn": 2}, 2: {"Cu": 2}, 3: {"V": 3, "O": 3}, 4: {"Zr": 6, "O": 6}}

md = {0: "Zn", 1: "Zn", 2: "Cu", 3: "V", 4: "Zr"}


def readCIF(path):
    f = open(path, "r")
    if f.mode == 'r':
        cif = f.read()
        return cif

def saveCIF(new_cif, path):
    '''
    This function saves new_cif to path
    '''
    f = open(path, "w")
    if f.mode == 'w':
        f.write(new_cif)

def rejoin(new_cif_l):
    '''
    This file turns a list of strings into a file 
    '''
    return "\n".join(new_cif_l)

def startsWithMetal(m_id, line):
    '''
    This function returns True if the line indicates a metal atom position
    '''
    if m_id == 0:
        if line.startswith('Zn'):
            return True
        else:
            return False
    else:
        m = md[m_id]
        if line.startswith(m):
            return True
        else:
            return False

def splitCif(cif):
    '''
    This function splits the cif file into a list of lines
    '''
    return cif.split('\n')

def saveMol(mol, write_path):
    '''
    This function save a molecule as an xyz to path
    '''
    XYZ(mol).write_file(write_path)

def getSmiles(mol):
    '''
    This function returns the SMILES string for a molecule
    '''
    a = BabelMolAdaptor(mol)
    pba = pb.Molecule(a.openbabel_mol)
    return pba.write("can").split()[0]

def subgraph_translate3(read_path, write_dir):
    '''
    This function takes the metal-deficient structure and saves the all connected components (i.e. linkers)
    '''
    ase_struct = read(read_path)
    struct = pm.get_structure(ase_struct)
    n_atoms = len(struct.sites)
    sg = StructureGraph.with_local_env_strategy(struct, env.JmolNN(tol=.3))
    count = 1
    for i in range(n_atoms):
        for j in range(n_atoms):
            if i != j:
                try:
                    sg.alter_edge(i, j, new_weight=count)
                except:
                    pass
                count += 1
    
    A = sg.get_subgraphs_as_molecules(use_weights=True)
    
    
    ls = []
    count = 0
    for mol in A:
        if len(mol.sites) > 4:
            smiles = getSmiles(mol)
            ls.append(smiles)
            path = write_dir + 'mol_' + str(count) + '.xyz'
            saveMol(mol, path)
            count += 1
    return count, ls

def removeMetal(m_id, cif, path, path_to_pristine_mof):
    '''
    This function filters out lines from cif which indicate the presence of a metal ion and saves it in path
    '''
    new_cif_l = []
    for line in splitCif(readCIF(path_to_pristine_mof)):
        if not (startsWithMetal(m_id, line)):
            new_cif_l.append(line)
    new_cif = rejoin(new_cif_l)
    saveCIF(new_cif, path)

def getMetalID(filename):
    '''
    This function will determine the Metal ID from the filename (does not work for full path)
    '''
    
    return int(filename.split('i')[2][1])

def getLinkerCols(df):
    '''
    This function returns the list of columns relevant to Linker SMILES
    '''
    
    return [col for col in df.columns if 'Smiles' in col]

def getSmilesSet(df):
    '''
    This function returns a set of all non-NA linkers in df
    '''
    cols = getLinkerCols(df)
    
    my_set = set()
    for col in cols:
        my_set.update(df[col])
    
    #my_set.remove('NA')
    try:
        my_set.remove('NA')
    except:
        try:
            my_set.remove(np.nan)
        except:
            pass
    return my_set    

def getSmilesMap(df):
    '''
    This function returns a map of the linker SMILES in the form of a dataframe to use in PG
    '''
    linker_set = getSmilesSet(df)
    n_linkers = len(linker_set)
    df = pd.DataFrame({"SMILES": list(linker_set)})
    a = df.rename_axis("ID")
    return a.dropna().reset_index().drop('ID', axis=1).rename_axis('ID')

def handle_file(tup):
    '''
    This function is carried out by each worker in the pool
    '''
    ind = tup[0]
    f = tup[1]
    n_mof = tup[2]
    parent_path = tup[3]
    cifs_dir = tup[4]

    if ind % 10 == 0:
            print("\n" + str(ind+1) + " out of " + str(n_mof) + " complete.")
    m_id = getMetalID(f)
    whole_path = parent_path+f
    f_name = whole_path.split('/')[-1][:-4]
    write_dir = ('/').join(whole_path.split('/')[:-1]) + '/' + f_name + '/'
    Path(write_dir).mkdir(parents=True, exist_ok=True)
    split_path = write_dir + 'split_' + f
    
    path_to_pristine_mof = whole_path = cifs_dir+f
    
    removeMetal(m_id, whole_path, split_path, path_to_pristine_mof)
    n_mol, smiles_l = subgraph_translate3(split_path, write_dir)
    return (n_mol, smiles_l, m_id, f)    

def getShortPath(path):
    '''
    This function returns just the filename from the absolute path
    '''
    
    return path.split('/')[-1]

def main4(parent_path, cifs_dir, bond_df_path, past_df_path):
    '''
    This is the main function for a multi-threaded process
    '''
    start = time.time()

    if not cifs_dir.endswith('/'):
        cifs_dir += '/'
    
    if not parent_path.endswith('/'):
        parent_path += '/'


    bond_df = pd.read_csv(bond_df_path, delim_whitespace=True)

    bonds = [bond.split('-') for bond in list(bond_df['Name'])]

    allowed_bonds = [set([i[0], i[1]]) for i in bonds]

    all_files = os.listdir(cifs_dir)
    #files = [f for f in all_files if ('cat' not in f and '.cif' in f)]
    files = [f for f in all_files if ('.cif' in f)]
    
    if past_df_path != None: #check to see if some files need to be ignored
        past_df = pd.read_csv(past_df_path, index_col=0)
        already_tested = past_df['filename'].tolist()
        all_short = [getShortPath(f) for f in files]
        files = [files[i] for i,f in enumerate(all_short) if f not in already_tested]
        
    n_mof = len(files)
    n_linker = []
    L_Smiles = [['NA']*n_mof for i in range(20)]
    m_ids = []

    data1 = list(enumerate(files))
    data2 = [(i[0], i[1], n_mof, parent_path, cifs_dir) for i in data1]

    pool = Pool(n_core)
    result = pool.map(handle_file, data2)
    pool.close()
    pool.join()

    for ind, i in enumerate(result):
        n_mol = i[0]
        smiles_l = i[1]
        m_id = i[2]
        f = i[3]
        for j in range(n_mol):
            if j < 20:
                try:
                    L_Smiles[j][ind] = smiles_l[j]
                except:
                    print("filename: ", f)
                    print("Number of Linkers: ", n_mol)
                    print("Current Linker index: ", i)
                    print("MOF number: ", ind)
                    print("Total number of MOFs", n_mof)
            else:
                n_mol = 20

        n_linker.append(n_mol)
        m_ids.append(m_id)

    end = time.time()
    print("Runtime: ", end-start)

    df = pd.DataFrame({"filename": files, "Metal ID": m_ids, "# of Linkers": n_linker,
                        "L0_Smiles": L_Smiles[0],
                        "L1_Smiles": L_Smiles[1],
                        "L2_Smiles": L_Smiles[2],
                        "L3_Smiles": L_Smiles[3],
                        "L4_Smiles": L_Smiles[4],
                        "L5_Smiles": L_Smiles[5],
                        "L6_Smiles": L_Smiles[6],
                        "L7_Smiles": L_Smiles[7],
                        "L8_Smiles": L_Smiles[8],
                        "L9_Smiles": L_Smiles[9],
                        "L10_Smiles": L_Smiles[10],
                        "L11_Smiles": L_Smiles[11],
                        "L12_Smiles": L_Smiles[12],
                        "L13_Smiles": L_Smiles[13],
                        "L14_Smiles": L_Smiles[14],
                        "L15_Smiles": L_Smiles[15],
                        "L16_Smiles": L_Smiles[16],
                        "L17_Smiles": L_Smiles[17],
                        "L18_Smiles": L_Smiles[18],
                        "L19_Smiles": L_Smiles[19]})
    
    if past_df_path != None:
        try:
            df = pd.concat([df, past_df], ignore_index=True)
        except:
            print("Data Frame concatentation failed. Only returning the new entries.")
            pass
    
    df.to_csv(parent_path + 'SMILESofMofs.csv')

    smiles_df = getSmilesMap(df)

    smiles_df.to_csv(parent_path + 'SmilesMap.csv')
    

main4(parent_dir, cifs_dir, bond_df_path, past_df_path)
