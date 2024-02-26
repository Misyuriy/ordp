import pandas as pd
import numpy as np

import joblib

from xgboost import XGBRegressor, XGBClassifier


def predict_solvent(data: pd.Series()) -> str:
    classifier = XGBClassifier()
    classifier.load('models/solvent_XGB')
    encoder = joblib.load('models/solvent_encoder.joblib')

    prediction_encoded = classifier.predict(data)

    return encoder.inverse_transform(prediction_encoded)


def predict_catalyst(data: pd.Series()) -> str:
    classifier = XGBClassifier()
    classifier.load('models/catalyst_XGB')
    encoder = joblib.load('models/catalyst_encoder.joblib')

    prediction_encoded = classifier.predict(data)

    return encoder.inverse_transform(prediction_encoded)


def predict_temperature(data: pd.Series()) -> float:
    regressor = XGBRegressor()
    regressor.load('models/temperature_XGB')

    return regressor.predict(data)
