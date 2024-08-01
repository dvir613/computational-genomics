import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from preprocess import preprocess_data


def train_model(file_path):
    # Load and preprocess data
    X, y = preprocess_data(file_path)

    # Split data into training and validation sets
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)

    # Train the model
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)

    # Evaluate the model
    score = model.score(X_val, y_val)
    print(f"Validation Score: {score:.2f}")

    return model
