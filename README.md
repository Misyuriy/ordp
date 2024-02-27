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

Example reaction: steam methane reforming

    C & O >>> [HH] & [C-]#[O+]

## Usage in code:
Prediction of catalyst for steam methane reforming:

    from ordp.predict import predict_catalyst
    from ordp.convert_reaction import convert_equation
    
    data = convert_equation('C & O >>> [HH] & [C-]#[O+]')
    
    print(predict_catalyst(data))

Alternative forms of reaction input are WIP


## Project status:
Work in progress:
- Alternative forms of reaction input (for example, dataframes)
- Predictions based on reaction clustering
- Improved compound clustering
- More features for usage in code

NOTE: temperature prediction accuracy is low.
