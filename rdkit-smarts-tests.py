from rdkit import Chem
from rdkit.Chem import Draw

import argparse
import pathlib

def mkdirp(location):
    sciezka = pathlib.Path(location)
    sciezka.mkdir(parents=True, exist_ok=True)


parser = argparse.ArgumentParser(description='out file prefix')
parser.add_argument('--prefix',
                    dest='prefix',
                    required=True,
                    help='input pefix')

args = parser.parse_args()
prefix = args.prefix


mkdirp("outputs")

##### FROM SMILES

# here provide the SMILES for the test molecule(s). If you don't have a better idea, leave as it is
testing_molecule = 'CC(C)(F)c1ccc(cc1)C[C@H](N)CC[NH3+].[H]O[C@@]1([H])[C@@]([H])(O[C@]([H])(C([H])([H])OP(O)([O-])=O)[C@@]1([H])O[H])n1c([H])nc2c1n([H])c(nc2=O)N([H])[H]'
smarts_pattern_to_test = '[!$([#1,#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]'

mother_mol = Chem.MolFromSmiles(testing_molecule)

mother_img = Chem.Draw.MolToFile(mother_mol,
                                 "outputs/%s--molecule-1.png" % (prefix),
                                 size=(400, 400),
                                 kekulize=False)

mol = Chem.MolFromSmiles(testing_molecule)

# add explicit hydrogens (optional)
# m = Chem.AddHs(m)

# ... and highlight desired atoms
substructure = Chem.MolFromSmarts(smarts_pattern_to_test)

triggeredatoms = mol.GetSubstructMatches(substructure)
triggeredatomsList = [x[0] for x in triggeredatoms]
print(triggeredatomsList)

print("Activated for these atoms: ", triggeredatoms)

# save m here again
img = Chem.Draw.MolToFile(mol,
                          "outputs/%s--molecule-2.png" % (prefix),
                          size=(400, 400),
                          kekulize=False,
                          highlightAtoms=triggeredatomsList)





#### FROM PDB

molPdb = Chem.rdmolfiles.MolFromPDBFile("3d2v.pdb")

# pdb_mother_img = Chem.Draw.MolToFile(molPdb,
#                                  "outputs/%s--molecule-pdb-1.png" % (prefix),
#                                  size=(400, 400),
#                                  kekulize=False)


substructure = Chem.MolFromSmarts(smarts_pattern_to_test)

triggeredatoms = molPdb.GetSubstructMatches(substructure)
triggeredatomsList = [x[0] for x in triggeredatoms]
print(triggeredatomsList)

# print("Activated for these atoms: ", triggeredatoms)


# save m here again
# img = Chem.Draw.MolToFile(molPdb,
#                           "outputs/%s--molecule-pdb-2.png" % (prefix),
#                           size=(400, 400),
#                           kekulize=False,
#                           highlightAtoms=triggeredatomsList)
