import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from xgboost import XGBRegressor
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from joblib import Parallel, delayed, dump
from tqdm import tqdm

# Function to identify extreme values by percentiles
def identify_extreme_percentiles(Y):
    numeric_Y = Y.select_dtypes(include=[np.number])
    lower_bound = numeric_Y.quantile(0.10)
    upper_bound = numeric_Y.quantile(0.90)
    Y['is_extreme'] = numeric_Y.apply(lambda row: any(row < lower_bound) or any(row > upper_bound), axis=1)
    return Y

# Function to split data with stratification based on extreme values
def split_data_with_extremes(X, Y, test_size=0.2, random_state=42):
    Y_extreme = identify_extreme_percentiles(Y)
    X_train, X_test, Y_train, Y_test = train_test_split(
        X, Y_extreme, test_size=test_size, random_state=random_state, stratify=Y_extreme['is_extreme']
    )
    Y_train.drop(columns='is_extreme', inplace=True)
    Y_test.drop(columns='is_extreme', inplace=True)
    return X_train, X_test, Y_train, Y_test

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

def apply_umap(X_train, X_test, X_TEST, n_components=30):
    umap = UMAP(n_components=n_components, random_state=42)
    X_train_umap = umap.fit_transform(X_train)
    X_test_umap = umap.transform(X_test)
    X_TEST_umap = umap.transform(X_TEST)
    return X_train_umap, X_test_umap, X_TEST_umap

def apply_pca(X_train, X_test, n_components):
    pca = PCA(n_components=n_components)
    X_train_pca = pca.fit_transform(X_train)
    X_test_pca = pca.transform(X_test)
    return X_train_pca, X_test_pca

def apply_tsne(X_train, X_test, n_components):
    combined_data = np.vstack((X_train, X_test))
    tsne = TSNE(n_components=n_components)
    combined_tsne = tsne.fit_transform(combined_data)
    X_train_tsne = combined_tsne[:X_train.shape[0], :]
    X_test_tsne = combined_tsne[X_train.shape[0]:, :]
    return X_train_tsne, X_test_tsne

def prepare_model_data(X: pd.DataFrame, Y: pd.DataFrame, method, n_components, outliers=False):
    X_TEST = pd.read_csv("../data/Test_data_with_features.csv")
    # Save predictions with DT metrics
    MODEL_df = pd.DataFrame()
    MODEL_df['name'] = X_TEST['Variant number']
    X_TEST.drop(columns='Variant number', inplace=True)
    if 'Variant number' in X.columns:
        X.drop(columns='Variant number', inplace=True)
    if 'name' in Y.columns:
        Y.drop(columns='name', inplace=True)

    if outliers:
        X, Y = remove_outliers(X, Y)

    # X_train, X_test, Y_train, Y_test = split_data_with_extremes(X, Y, test_size=0.2, random_state=42)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.15, random_state=42)
    X_train, X_test = scale(X_train, X_test)

    if method == 'UMAP':
        print(X_TEST.shape)
        X_train, X_test, X_TEST = apply_umap(X_train, X_test, X_TEST, n_components)
    elif method == 'PCA':
        X_train, X_test = apply_pca(X_train, X_test, n_components)
    elif method == 'TSNE':
        X_train, X_test = apply_tsne(X_train, X_test, n_components)

    return X_train, X_test, X_TEST, MODEL_df, Y_train, Y_test

def train_and_evaluate_model(X_train, X_test, X_TEST, MODEL_df, Y_train, Y_test, model, save_test_predictions=False, model_name=None):
    model.fit(X_train, Y_train)
    predictions = model.predict(X_test)

    if save_test_predictions:
        predictions_test = model.predict(X_TEST)
        MODEL_df['dt_average'] = predictions_test[:, 0]
        MODEL_df['DT_max'] = predictions_test[:, 1]
        MODEL_df.to_csv(f'{model_name}_Predictions.csv', index=False)

    mse = mean_squared_error(Y_test, predictions, multioutput='raw_values')
    r2 = r2_score(Y_test, predictions, multioutput='raw_values')
    spearman_corr = [spearmanr(Y_test.iloc[:, i], predictions[:, i])[0] for i in range(Y_test.shape[1])]
    return {"MSE": mse, "R2": r2, "Spearman": spearman_corr, "Model": model}

def try_dimensionality_reduction_methods(X, Y, methods, n_components_range):
    best_model = None
    best_score = 0
    best_results = None

    for method in methods:
        for n_components in range(*n_components_range):
            print(f"Method: {method}, Components: {n_components}")
            try:
                X_train, X_test, X_TEST, MODEL_df, Y_train, Y_test = prepare_model_data(X, Y, method, n_components, outliers=True)

                model_results = {}
                models = {
                    "RandomForestRegressor": RandomForestRegressor(n_estimators=450, random_state=42, n_jobs=-1),
                    "DecisionTreeRegressor": DecisionTreeRegressor(random_state=42, max_depth=300),
                    # Add XGBoost if needed with proper GPU configuration
                    "XGBRegressor": XGBRegressor(objective='reg:squarederror', n_estimators=100, random_state=42, tree_method='hist', device='cuda')
                }

                for model_name, model in models.items():
                    print(f"Model: {model_name}")
                    result = train_and_evaluate_model(X_train, X_test, X_TEST, MODEL_df, Y_train, Y_test, model,
                                                      save_test_predictions=False, model_name=model_name)
                    model_results[model_name] = result

                    for target_idx, target_name in enumerate(Y.columns[:2]):
                        print(f"  Target: {target_name}")
                        print(f"    MSE: {result['MSE'][target_idx]:.4f}")
                        print(f"    R2: {result['R2'][target_idx]:.4f}")
                        print(f"    Spearman Correlation: {result['Spearman'][target_idx]:.4f}")

                    # Determine if this model is the best based on spearman correlation
                    # avg_mse = np.mean(result['MSE'])
                    avg_spearman = np.mean(result['Spearman'])
                    if avg_spearman > best_score:
                        best_score = avg_spearman
                        print(f"New best score: {best_score:.4f}")
                        best_model = result['Model']
                        best_results = result


                best_results[(method, n_components)] = model_results
            except Exception as e:
                print(f"Error with method {method} and components {n_components}: {e}")

    # Save the best model
    if best_model is not None:
        dump(best_model, 'best_model.joblib')
        print(f"Best model saved with average spearman correlation: {best_score:.4f}")

    return best_results

def train_model(X_path, Y_path):
    X = pd.read_csv(X_path)
    Y = pd.read_csv(Y_path)

    methods = ['UMAP']
    n_components_range = (37, 38)

    results = try_dimensionality_reduction_methods(X, Y, methods, n_components_range)
    return results

train_model("../data/Train_data_with_features.csv", "../data/Target_data.csv")
