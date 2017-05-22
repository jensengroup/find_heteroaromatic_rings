import sys, os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

rings_mol = []
rings_smiles = []

smiles_file_name = "/Users/jan/Dropbox/Lundbeck/big.smiles"

smiles_file = open(smiles_file_name, "r")

for line in smiles_file:
    words = line.split()
    name = words[0]
    smiles = words[1]
    
    mol =  Chem.MolFromSmiles(smiles)

    bis = mol.GetSubstructMatches(Chem.MolFromSmarts('[R]-[*]'))
    bs = [mol.GetBondBetweenAtoms(x,y).GetIdx() for x,y in bis]
    
    if len(bs) == 0: 
        if smiles not in rings_smiles: 
            rings_smiles.append(smiles) 
            rings_mol.append(Chem.MolFromSmiles(smiles)) 
        continue
    
    fragments_mol = Chem.FragmentOnBonds(mol,bs,addDummies=False)
    
    fragments = Chem.MolToSmiles(fragments_mol,True).split(".")
    
    ring = Chem.MolFromSmarts('[R]')
    
    for fragment in fragments:
        if Chem.MolFromSmiles(fragment).HasSubstructMatch(ring):
            if fragment not in rings_smiles and "n" in fragment: 
                rings_mol.append(Chem.MolFromSmiles(fragment))
                rings_smiles.append(fragment) 

img = Draw.MolsToGridImage(rings_mol,molsPerRow=4,subImgSize=(200,200),useSVG=True)

svg_file_name = "/Users/jan/Dropbox/Lundbeck/rings.svg"
svg_file = open(svg_file_name, 'w')
svg_file.write(img.data)
svg_file.close()
os.system('sed -i "" "s/xmlns:svg/xmlns/" '+svg_file_name)
