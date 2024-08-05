import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression, Lasso, Ridge
from xgboost import XGBRegressor
from sklearn.metrics import mean_squared_error, r2_score, make_scorer
from sklearn.feature_selection import SequentialFeatureSelector
from scipy.stats import spearmanr
from joblib import Parallel, delayed
from tqdm import tqdm

def remove_outliers(X: pd.DataFrame, Y: pd.DataFrame):
    q1, q3 = np.percentile(Y, [25, 75])
    iqr = q3 - q1
    lower_fence = q1 - (1.5 * iqr)
    higher_fence = q3 + (1.5 * iqr)
    mask = np.all((Y > lower_fence) & (Y < higher_fence), axis=1)
    X = X.loc[mask]
    Y = Y.loc[mask]
    return X, Y

def scale(X1, X2):
    scaler = StandardScaler()
    scaler.fit(X1)
    X1 = pd.DataFrame(scaler.transform(X1), columns=X1.columns)
    X2 = pd.DataFrame(scaler.transform(X2), columns=X2.columns)
    return X1, X2

def prepare_model_data(X: pd.DataFrame, Y: pd.DataFrame, outliers=False):
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
    X_train, X_val, Y_train, Y_val = prepare_model_data(X, Y, outliers=False)

    # Evaluate models
    results = evaluate_models(X_train, X_val, Y_train, Y_val)

    # Print the results
    for name, metrics in results.items():
        print(f"Model: {name}")
        for i, target in enumerate(Y.columns):
            print(f"  Target: {target}")
            print(f"    MSE: {metrics['MSE'][i]:.4f}")
            print(f"    R2: {metrics['R2'][i]:.4f}")
            print(f"    Spearman Correlation: {metrics['Spearman'][i]:.4f}")

    return results

def spearman_scorer(y_true, y_pred):
    return spearmanr(y_true, y_pred)[0]

def parallel_feature_selection(model, X_train, Y_train, n_features=40):
    sfs = SequentialFeatureSelector(
        model.estimator,
        n_features_to_select=n_features,
        direction='forward',
        scoring=make_scorer(mean_squared_error),
        n_jobs=-1
    )
    sfs.fit(X_train, Y_train)
    return sfs

def evaluate_models(X_train, X_val, Y_train, Y_val):
    models = {
        "RandomForestRegressor": MultiOutputRegressor(RandomForestRegressor(n_estimators=350, random_state=42)),
        "LinearRegression": MultiOutputRegressor(LinearRegression()),
        "Lasso": MultiOutputRegressor(Lasso()),
        "Ridge": MultiOutputRegressor(Ridge()),
        "XGBRegressor": MultiOutputRegressor(XGBRegressor(objective='reg:squarederror', n_estimators=100, random_state=42, tree_method='gpu_hist'))
    }

    def train_and_evaluate_model(name, model):
        sfs = parallel_feature_selection(model, X_train, Y_train, n_features=40)
        X_train_sfs = sfs.transform(X_train)
        X_val_sfs = sfs.transform(X_val)
        model.fit(X_train_sfs, Y_train)
        predictions = model.predict(X_val_sfs)
        mse = mean_squared_error(Y_val, predictions, multioutput='raw_values')
        r2 = r2_score(Y_val, predictions, multioutput='raw_values')
        spearman_corr = [spearmanr(Y_val.iloc[:, i], predictions[:, i])[0] for i in range(Y_val.shape[1])]
        return name, {"MSE": mse, "R2": r2, "Spearman": spearman_corr}

    results = Parallel(n_jobs=-1)(delayed(train_and_evaluate_model)(name, model) for name, model in tqdm(models.items()))
    return dict(results)


train_model("../data/Train_data_with_features.csv", "../data/Target_data.csv")
