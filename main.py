import warnings

import predict
from convert_reaction import convert_equation


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


if __name__ == '__main__':
    print('Equation syntax:')
    print(bcolors.OKCYAN + 'reagent & reagent & reagent >>> product & product' + bcolors.ENDC)
    print('reagents and products should be entered in SMILES')
    encoding = None

    while encoding is None:
        print()
        print('Enter equation:')

        equation = input()

        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore')

            try:
                encoding = convert_equation(equation)
            except ValueError:
                print(bcolors.FAIL + 'Equation syntax invalid.' + bcolors.ENDC)
                print('Valid equation must follow this syntax:')
                print(bcolors.OKCYAN + 'reagent & reagent & reagent >>> product & product' + bcolors.ENDC)
            except:
                print(bcolors.FAIL + 'Compounds not readable. Ensure compound SMILES strings are valid.' + bcolors.ENDC)
                print('To get valid SMILES representation, find compound on PubChem: https://pubchem.ncbi.nlm.nih.gov/')

    print()
    print('Equation valid. Predicting conditions...')

    print('Temperature prediction:', end=' ')
    temperature = predict.predict_temperature(encoding)
    print(bcolors.OKCYAN + str(round(temperature, 2)) + ' K (' + str(round(temperature - 273.15, 2)) + '°C)' + bcolors.ENDC)

    print('Solvent prediction:', end=' ')
    solvent = predict.predict_solvent(encoding)
    print(bcolors.OKCYAN + solvent + bcolors.ENDC)

    print('Catalyst prediction:', end=' ')
    catalyst = predict.predict_catalyst(encoding)
    print(bcolors.OKCYAN + catalyst + bcolors.ENDC)