# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at the University of the Western Cape
# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV Integrase Subtype C and CRF02_AG" 
# Currently this project is funded by the Poliomyeletis Research Foundation
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties

#Author:	Keaghan Brown (3687524) - MSc Candidate at South African National Bioinformatics Institute (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

#Automated Mutation Introduction and Analysis - AMIA (Version1.8.2)
from tkinter import *
from tkinter.filedialog import askopenfilename,  askdirectory
from tkinter import messagebox
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import ntpath
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
import pymol
import os

class Mutations:

    def __init__(self):
        """Suplies dictionaries for global use"""
        self.drug_mutation= {}
        self.three_letter ={'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
    
    def get_mut_list(self, mutations):
        """ Retrieves the mutations and appends it to the dictionary for further processing """
        self.mutations = mutations
        if self.mutations.endswith(".txt"):
            my_file = open(self.mutations)
            file_data = my_file.read()
            my_file.close()
            mutation_data = file_data.split("\n")
            while "" in mutation_data:
                mutation_data.remove("")
            mutation_list_counter = 0
            for i in mutation_data:
                if i.startswith('***'):
                    mutation_data.remove(i)
            for j in mutation_data:
                mutation_list_counter+=1
                if j == '###':                                                                                                                                                             
                    drug_name = mutation_list_counter                                                                                                                   
                elif j == "####":                                                                                                                                                                   
                    first_mutation = mutation_list_counter                                                                                                                       
                elif j == "##":                                                                                                                                                                       
                    last_mutation = mutation_list_counter - 1                                                                                                                 
                    self.drug_mutation[mutation_data[drug_name]] = (str(first_mutation) + ','  + str(last_mutation))      
            for key in self.drug_mutation:                                                                                                                                     
                start_stop_pos = self.drug_mutation[key].split(',')
                self.drug_mutation[key] = mutation_data[int(start_stop_pos [0]) : int(start_stop_pos [1])] 
        return (self.drug_mutation)
    	
    def single_mutation_intro(self, structure):   
        """ Introduces single mutations from the dictionary, into individual structures and saves it to a specified directory """
        self.structure =   structure
        PDB_storage=  askdirectory(title='Select Mutant Storage Folder')
        for key in self.drug_mutation:
            for mutation in  self.drug_mutation[key]:
                pymol.cmd.reinitialize()
                mutated_residue = mutation[len(mutation)-1]
                residue_pos = int(mutation[1:len(mutation)-1])
                os.chdir(os.path.dirname(self.structure))
                parser = PDBParser(PERMISSIVE=1)
                structure_id = ntpath.basename(self.structure).split('.')[0]
                filename = ntpath.basename(self.structure)
                structure = parser.get_structure(structure_id, filename)
                ppb=PPBuilder()
                ppb.build_peptides(structure)
                start_pos = []
                for pp in ppb.build_peptides(structure):
                    start_pos = int((pp[0].get_id()[1]))
                pymol.cmd.load(self.structure)
                pymol.cmd.wizard('mutagenesis')
                pymol.cmd.do("refresh_wizard")
                pymol.cmd.get_wizard().set_mode(self.three_letter[mutated_residue])
                pymol.cmd.get_wizard().do_select(str((start_pos + residue_pos) -1)+ '/')
                pymol.cmd.get_wizard().apply()
                pymol.cmd.select('mutation', 'resi ' + str(residue_pos) + ' around 8')
                pymol.cmd.clean('mutation')
                os.chdir(PDB_storage)
                pymol.cmd.save(str(key) + mutation[0] + str((start_pos + residue_pos) -1) + mutated_residue + '.pdb')
            messagebox.showinfo("Information","Mutated " + str(key.split(':')[0]) + " Structures Generated")
    	
    def multi_mutation_intro(self, structure):   
        """ Introduces multiple mutations from the dictionary and into a single structure and saves it to a specified directory """
        self.structure =   structure
        PDB_storage=  askdirectory(title='Select Mutant Storage Folder')
        position_record = []
        pymol.cmd.load(self.structure)
        for key in self.drug_mutation:
            for mutation in self.drug_mutation[key]:
                mutated_residue = mutation[len(mutation)-1]
                residue_pos = int(mutation[1:len(mutation)-1])
                if residue_pos in position_record:
                    messagebox.showwarning("Warning","A mutation already exists at position " + str(residue_pos) + " and will be replaced!")
                position_record.append(residue_pos)
                os.chdir(os.path.dirname(self.structure))
                parser = PDBParser(PERMISSIVE=1)
                structure_id = ntpath.basename(self.structure).split('.')[0]
                filename = ntpath.basename(self.structure)
                structure = parser.get_structure(structure_id, filename)
                ppb=PPBuilder()
                ppb.build_peptides(structure)
                start_pos = []
                for pp in ppb.build_peptides(structure):
                    start_pos = int((pp[0].get_id()[1]))
                pymol.cmd.wizard('mutagenesis')
                pymol.cmd.do("refresh_wizard")
                pymol.cmd.get_wizard().set_mode(self.three_letter[mutated_residue])
                pymol.cmd.get_wizard().do_select(str((start_pos + residue_pos) -1)+ '/')
                pymol.cmd.get_wizard().apply()
                pymol.cmd.select('mutation', 'resi ' + str(residue_pos) + ' around 8')
                pymol.cmd.clean('mutation')
        os.chdir(PDB_storage)
        pymol.cmd.save("Multiple Mutation Intorudction.pdb")
        messagebox.showinfo("Information","Mutated Structure Generated")      
