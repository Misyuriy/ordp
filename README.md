# ORDP condition prediction

Python utility for chemical reaction planning.

Predicts:
- optimal temperature
- solvent
- catalyst

## Installation:
    git clone https://github.com/Misyuriy/ordp.git

## Running the project:
    python main.py

## Usage from terminal:
The program accepts **Equation** syntax:

    reagent & reagent >>> product & product

Each compound (reagent or product) must be a valid SMILES string.
To get valid SMILES representation, find compound on PubChem: https://pubchem.ncbi.nlm.nih.gov/

Example reaction: benzaldehyde with acetic acid

    C(C)(=O)O & C(C1=CC=CC=C1)=O >>> C=C=O & C=C1CC(=O)O1

## Usage in code:
Prediction of catalyst for benzaldehyde reaction with acetic acid:

    from ordp.predict import predict_catalyst
    from ordp.convert_reaction import convert_equation
    
    data = convert_equation('C(C)(=O)O & C(C1=CC=CC=C1)=O >>> C=C=O & C=C1CC(=O)O1')
    
    print(predict_catalyst(data))

Alternative forms of reaction input are WIP


## Project status:
Work in progress:
- Alternative forms of reaction input (for example, dataframes)
- Predictions based on reaction clustering
- Improved compound clustering
- More features for usage in code

NOTE: temperature prediction accuracy is low.

NOTE 2: the models are almost useless for inorganic reactions. Predictions will be random
