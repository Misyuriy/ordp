import numpy as np
import pandas as pd

import joblib

from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def fingerprint(smiles: str, generator):
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


def get_known_equation(equation: str, max_reactants=8, max_products=6):
    r, p = equation.split('>>>')
    reactants = r.split('&')
    products = p.split('&')

    result = {}
    for i in range(max_reactants):
        result['reactant_' + str(i)] = -1
    for i in range(max_products):
        result['product_' + str(i)] = -1

    for index, reactant in enumerate(reactants):
        reactants[index] = Chem.CanonSmiles(reactant.strip(' '))
    for index, product in enumerate(products):
        products[index] = Chem.CanonSmiles(product.strip(' '))

    products = set([Chem.CanonSmiles(i) for i in products])
    reactants = set([Chem.CanonSmiles(i) for i in reactants])

    paths = ['known/known_more_products.csv',
             'known/known_2_product.csv',
             'known/known_1_product.csv']

    for path in paths:
        known = pd.read_csv(path)

        for row in known.iloc:
            known_reactants = [row['reactant_' + str(i)] for i in range(max_reactants)]
            known_reactants = set([(Chem.CanonSmiles(i) if type(i) == str else i) for i in known_reactants])

            if None in known_reactants:
                known_reactants.remove(None)
            if np.nan in known_reactants:
                known_reactants.remove(np.nan)

            if known_reactants == set(reactants):
                return row['temperature'], row['solvent_0'], row['catalyst_0']

    return None, None, None
