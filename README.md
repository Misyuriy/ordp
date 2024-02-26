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

## Usage:
The program accepts **Equation** syntax:

    reagent & reagent >>> product & product

Each compound (reagent or product) must be a valid SMILES string.
To get valid SMILES representation, find compound on PubChem: https://pubchem.ncbi.nlm.nih.gov/

Example reaction: steam methane reforming

    C & O >>> [HH] & [C-]#[O+]

## Project status:
Currently, the models are a work in progress. 
Temperature prediction accuracy is low.
