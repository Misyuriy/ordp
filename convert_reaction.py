import numpy as np
import pandas as pd

import joblib

from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def fingerprint(smiles, generator):
    mol = Chem.MolFromSmiles(smiles)
    fp = generator.GetFingerprint(mol)

    return np.unpackbits(np.frombuffer(DataStructs.BitVectToBinaryText(fp), dtype=np.uint8), bitorder='little')


def convert_equation(equation: str, max_reactants=8, max_products=6):
    r, p = equation.split('>>>')
    reactants = r.split('&')
    products = p.split('&')

    result = {}
    for i in range(max_reactants):
        result['reactant_' + str(i)] = -1
    for i in range(max_products):
        result['product_' + str(i)] = -1

    clustering = joblib.load('models/compound_clustering.joblib')
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    for index, reactant in enumerate(reactants):
        compound = reactant.strip(' ')
        result['reactant_' + str(index)] = clustering.predict([fingerprint(compound, generator)])[0]

    for index, product in enumerate(products):
        compound = product.strip(' ')
        result['product_' + str(index)] = clustering.predict([fingerprint(compound, generator)])[0]

    return pd.DataFrame([pd.Series(result)])
