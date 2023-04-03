#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LATEST VERSION 27 JAN 2022


"""

import os
import sys
import PeptideBuilder.PeptideBuilder as PeptideBuilder
import Bio.PDB
from biopandas.pdb import PandasPdb

def xml_transform(geo_list, outTempPdb):
    if isinstance(geo_list, str):
        geo_list = list(geo_list)
        
    structure = PeptideBuilder.make_structure_from_geos(geo_list)

    # add terminal oxygen (OXT) to the final glycine
    PeptideBuilder.add_terminal_OXT(structure)
        
    # read the pdb file
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save(outTempPdb)
 
def centerPDB(inPDB,outPDB):
    pd = PandasPdb()
    data1 = pd.read_pdb(inPDB)  
    
    rows = data1.df['ATOM'].shape[0]
    columns = data1.df['ATOM'].shape[1]
    
    # find median
    x_median = data1.df['ATOM']['x_coord'].median()
    y_median = data1.df['ATOM']['y_coord'].median()
    z_median = data1.df['ATOM']['z_coord'].median()
    
    # move data frame
    data1.df['ATOM']['x_coord'] += -x_median
    data1.df['ATOM']['y_coord'] += -y_median
    data1.df['ATOM']['z_coord'] += -z_median
        
    #print(data1.df['ATOM'])
    
    data1.to_pdb(path = outPDB,
                records = None,
                gz = False,
                append_newline = True)
    
    return rows

def shiftPDB(inPDB,outPDB,x_move=5, y_move=5, z_move=5, chain_id="A", atom_number_move=0, residue_number_move = 0):
    pd = PandasPdb()
    data = pd.read_pdb(inPDB)
    
    rows = data.df['ATOM'].shape[0]
    columns = data.df['ATOM'].shape[1]
    
    # move data frame
    data.df['ATOM']['x_coord'] += x_move
    data.df['ATOM']['y_coord'] += y_move
    data.df['ATOM']['z_coord'] += z_move
    data.df['ATOM']['atom_number'] += atom_number_move
    data.df['ATOM']['chain_id'] = chain_id
    data.df['ATOM']['residue_number'] += residue_number_move
    
    data.to_pdb(path = outPDB,
                records = None,
                gz = False,
                append_newline = True)

    
    
def calBoxSize(concentration_mM=30, chains=0):
    
    concentration = concentration_mM * 10**(-3)
    chains = 2
    volume = chains * 10**27 / (concentration * 6.02*10**23)
    cube_root = volume ** (1/3) 
    
    cryst = open('cryst1.pdb', 'w')
    cryst.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f \n" %(cube_root,cube_root,cube_root,90,90,90))
    cryst.close()   
    
    return (cube_root)
  
 
    

def concatenate_files(flist):
    """
    Iterates over a list of files and yields each line sequentially.
    """

    for fhandle in flist:
        for line in fhandle:
            yield line
        fhandle.close()


def merger(headerFiles, PDBfiles, filename):
    fl = []
    new = open(filename, 'w')
    
    # add header
    if len(headerFiles) >= 1:
        for fn in headerFiles:
            fh = open(fn, 'r')
            
            for line in fh:
                new.write(line)
            
    # add pdb except last two lines
    if len(PDBfiles) >= 1:
        for fn in PDBfiles:
            fh = open(fn, 'r')
            
            for line in fh.readlines()[:-2]:
                new.write(line)
                
    new.write("END\n")
    new.close()


     
        

def amino_to_pdb_twoChains(geo_list1, geo_list2, filename):
    xml_transform("G"+geo_list1+"G",'temp1.pdb')
    xml_transform("G"+geo_list2+"G",'temp2.pdb')
    rows1 = centerPDB("temp1.pdb","center1.pdb")
    roww2 = centerPDB("temp2.pdb","center2.pdb")
    
    shiftPDB("center1.pdb", "modified1.pdb", x_move = 30, y_move = 30, z_move = 30)
    shiftPDB("center2.pdb", "modified2.pdb", x_move = 20, y_move = 20, z_move = 20, chain_id="B", atom_number_move = rows1, residue_number_move=6)
   
    calBoxSize(30)
    merger(["cryst1.pdb"], ["modified1.pdb", "modified2.pdb"], filename)
    
    if os.path.exists("cryst1.pdb"):
        os.remove("cryst1.pdb")
    if os.path.exists("modified1.pdb"):
        os.remove("modified1.pdb")
    if os.path.exists("temp1.pdb"):
        os.remove('temp1.pdb')
    if os.path.exists("modified2.pdb"):
        os.remove("modified2.pdb")
    if os.path.exists("temp2.pdb"):
        os.remove('temp2.pdb') 
    if os.path.exists("center1.pdb"):
        os.remove('center1.pdb')
    if os.path.exists("center2.pdb"):
        os.remove('center2.pdb')  



        
if __name__ == '__main__':
    geo_list1 = "HGAVIL"  # change to the amino acid sequence of interest
    geo_list2 = "HGAVIL" # change the aminio acid sequence of interest
    start_pdb = geo_list1 + "_" + geo_list2 + "_start.pdb"
    amino_to_pdb_twoChains(geo_list1, geo_list2, start_pdb)

