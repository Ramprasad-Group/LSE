# LSE
LSE = Linker SMILES Extractor. Return the SMILES strings of several input MOFs' constituent linkers. This work was done by Rishi Gurnani.

# Installation
Step 1: Use Python3.7 <br />
Step 2: Install OpenBabel with python bindings <br />
Step 3: Install all packages listed in requirements.txt

# Usage
To use the package try the following command:
```
python main.py [path to folder to dump output into] [number of cpu cores] [path to folder containing CIF files]
```

A few things to note:

1) All .cif files which you want to extract from MUST be located in the same folder.
2) The outputs of the code are as follows: <br />
  -SMILESofMofs.csv: contains filename, # of Linkers, and SMILES string of each linker in MOF (up to 20 SMILES) <br />
  -SmilesMap.csv: contains all unique linkers in all MOFs input to the code. <br />
  -For each .cif file: <br />
  &ensp;A folder with name matching the prefix of the .cif file. The folder contains the .cif file with metallic atoms removed (split_*.cif) and the .xyz files of each linker
  
To test the package try the following command. Please allow up to five minutes for results. Most MOFs do not require this much amount of time.
```
python main.py ./ 2 test/inputs/
```

The files which are produced should match the files in ```test/outputs/```

# Disclaimer
All .cif files located in ```test/inputs/``` are courtesy of the following publication: C.E. Wilmer, M. Leaf, C.Y. Lee, O.K. Farha, B.G. Hauser, J.T. Hupp, and R.Q. Snurr, Nat. Chem. 4, 83 (2012).

# License
This work is protected under the license titled 'LICENSE' in this repository.
