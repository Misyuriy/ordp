import pandas as pd
import numpy as np

import joblib

from xgboost import XGBRegressor, XGBClassifier


def predict_solvent(data: pd.Series()) -> str:
    classifier = XGBClassifier()
    classifier.load_model('models/solvent_XGB.json')
    encoder = joblib.load('models/solvent_encoder.joblib')

    data = data[classifier.get_booster().feature_names]
    prediction_encoded = classifier.predict(data)

    return encoder.inverse_transform(prediction_encoded)[0]


def predict_catalyst(data: pd.Series()) -> str:
    classifier = XGBClassifier()
    classifier.load_model('models/catalyst_XGB.json')
    encoder = joblib.load('models/catalyst_encoder.joblib')

    data = data[classifier.get_booster().feature_names]
    prediction_encoded = classifier.predict(data)

    return encoder.inverse_transform(prediction_encoded)[0]


def predict_temperature(data: pd.Series()) -> float:
    regressor = XGBRegressor()
    regressor.load('models/temperature_XGB')

    data = data[regressor.get_booster().feature_names]

    return regressor.predict(data)
