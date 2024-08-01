import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from preprocess import preprocess_data


def predict(model: RandomForestRegressor, file_path: str, output_path: str):
    # Load and preprocess test data
    df = pd.read_csv(file_path)
    X_test = df.drop(columns=['target'])

    # Make predictions
    predictions = model.predict(X_test)

    # Save predictions
    results = pd.DataFrame({'id': df['id'], 'predictions': predictions})
    results.to_csv(output_path, index=False)
