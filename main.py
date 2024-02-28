import predict
from convert_reaction import convert_equation, get_known_equation

import warnings
warnings.filterwarnings("ignore")


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def main():
    print('Equation syntax:')
    print(bcolors.OKCYAN + 'reagent & reagent & reagent >>> product & product' + bcolors.ENDC)
    print('reagents and products should be entered in SMILES')
    encoding = None

    equation = 'm'

    while encoding is None:
        print()
        print('Enter equation:')

        equation = input()

        try:
            encoding = convert_equation(equation)
        except ValueError:
            print(bcolors.FAIL + 'Equation syntax invalid.' + bcolors.ENDC)
            print('Valid equation must follow this syntax:')
            print(bcolors.OKCYAN + 'reagent & reagent & reagent >>> product & product' + bcolors.ENDC)
        except:
            print(bcolors.FAIL + 'Compounds not readable. Ensure compound SMILES strings are valid.' + bcolors.ENDC)
            print('To get valid SMILES representation, find compound on PubChem: https://pubchem.ncbi.nlm.nih.gov/')

    print('Equation valid. Predicting conditions...')

    #temperature, solvent, catalyst = None, None, None
    temperature, solvent, catalyst = get_known_equation(equation)

    print('Temperature prediction:', end=' ')
    p_temperature = predict.predict_temperature(encoding)
    if not temperature:
        temperature = p_temperature

    print(bcolors.OKCYAN + str(round(temperature, 2)) + ' K (' + str(
        round(temperature - 273.15, 2)) + '°C)' + bcolors.ENDC)

    print('Solvent prediction:', end=' ')
    p_solvent = predict.predict_solvent(encoding)
    if not solvent:
        solvent = p_solvent

    print(bcolors.OKCYAN + solvent + bcolors.ENDC)

    print('Catalyst prediction:', end=' ')
    p_catalyst = predict.predict_catalyst(encoding)
    if not catalyst:
        catalyst = p_catalyst

    print(bcolors.OKCYAN + catalyst + bcolors.ENDC)


if __name__ == '__main__':
    main()
