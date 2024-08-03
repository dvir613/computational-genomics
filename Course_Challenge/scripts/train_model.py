import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.multioutput import MultiOutputRegressor
# from lazypredict.Supervised import LazyRegressor
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import mean_squared_error, r2_score

# from preprocess import preprocess_data


def remove_outliers(X: pd.DataFrame, Y: pd.DataFrame):
    q1, q3 = np.percentile(Y, [25, 75])
    iqr = q3 - q1
    lower_fence = q1 - (1.5 * iqr)
    higher_fence = q3 + (1.5 * iqr)
    mask = np.all((Y > lower_fence) & (Y < higher_fence), axis=1)
    X = X.loc[mask]
    Y = Y.loc[mask]
    # X = X[(y > lower_fence) & (y < higher_fence)]
    # y = y[(y > lower_fence) & (y < higher_fence)]
    return X, Y


def scale(X1, X2):
    scaler = StandardScaler()
    scaler.fit(X1)
    X1 = pd.DataFrame(scaler.transform(X1), columns=X1.columns)
    X2 = pd.DataFrame(scaler.transform(X2), columns=X2.columns)
    return X1, X2


def prepare_model_data(X: pd.DataFrame, Y: pd.DataFrame, outliers=False):
    # remove the first row (the control sequence)
    X.drop(index=0, inplace=True)
    X.reset_index(inplace=True)  # update the indexing to start from 0 and not from 1
    # and remove the first column (the variant number/name)
    X.drop(columns='Variant number', inplace=True)
    Y.drop(columns='name', inplace=True)

    if outliers:
        X, Y = remove_outliers(X, Y)

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)

    X_train, X_test = scale(X_train, X_test)

    return X_train, X_test, Y_train, Y_test


def train_model(X_path, Y_path):
    X = pd.read_csv(X_path)
    Y = pd.read_csv(Y_path)
    X_train, X_val, Y_train, Y_val = prepare_model_data(X, Y, outliers=True)

    # Evaluate models
    results = evaluate_models(X_train, X_val, Y_train, Y_val)

    # Print the results
    for name, metrics in results.items():
        print(f"Model: {name}")
        for i, target in enumerate(Y.columns):
            print(f"  Target: {target}")
            print(f"    MSE: {metrics['MSE'][i]:.4f}")
            print(f"    R2: {metrics['R2'][i]:.4f}")

    return results


def evaluate_models(X_train, X_val, Y_train, Y_val):
    models = {
        "RandomForestRegressor": MultiOutputRegressor(RandomForestRegressor(n_estimators=100, random_state=42)),
        "LinearRegression": MultiOutputRegressor(LinearRegression())
        # "KNeighborsRegressor": MultiOutputRegressor(KNeighborsRegressor()),
        # "DecisionTreeRegressor": MultiOutputRegressor(DecisionTreeRegressor())
    }

    results = {}
    for name, model in tqdm(models.items()):
        model.fit(X_train, Y_train)
        predictions = model.predict(X_val)
        mse = mean_squared_error(Y_val, predictions, multioutput='raw_values')
        r2 = r2_score(Y_val, predictions, multioutput='raw_values')
        results[name] = {"MSE": mse, "R2": r2}

    return results


train_model("../data/Train_data_with_features.csv", "../data/Target_data.csv")
