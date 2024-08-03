import pandas as pd
import matplotlib.pyplot as plt


def calculate_differences(excel_file_path):
    # Load the data from the sheets
    without_dnt_df = pd.read_excel(excel_file_path, sheet_name='luminescence without DNT', engine='openpyxl')
    with_dnt_df = pd.read_excel(excel_file_path, sheet_name='luminescence with DNT', engine='openpyxl')

    # Check if both dataframes have the same shape
    if without_dnt_df.shape != with_dnt_df.shape:
        raise ValueError("The sheets do not have the same shape.")

    # Perform the subtraction
    difference_df = with_dnt_df.iloc[:, 1:] - without_dnt_df.iloc[:, 1:]

    # Calculate the mean and max difference for each variant
    mean_diff = difference_df.mean(axis=1)
    max_diff = difference_df.max(axis=1)

    y_data = pd.DataFrame({
        'name': without_dnt_df.iloc[:, 0],
        'dt_average': mean_diff,
        'DT_max': max_diff
    })

    return y_data


def plot_differences(ydata):
    # Plot scatter plot for mean difference values
    plt.figure(figsize=(12, 6))

    # Scatter plot for mean difference
    plt.subplot(1, 2, 1)
    plt.scatter(ydata['name'], ydata['dt_average'], color='blue')
    plt.title('Mean Difference per Variant')
    plt.xlabel('Variant number')
    plt.ylabel('Mean Difference')

    # Scatter plot for max difference
    plt.subplot(1, 2, 2)
    plt.scatter(ydata['name'], ydata['DT_max'], color='red')
    plt.title('Max Difference per Variant')
    plt.xlabel('Variant number')
    plt.ylabel('Max Difference')

    # Adjust layout to avoid overlap
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    scatter_plot = True
    save_target_data = False
    excel_file_path = "../data/Train_data.xlsx"
    y_data = calculate_differences(excel_file_path)

    # Print the resulting matrix of differences
    print(y_data.head(5))

    if scatter_plot:
        plot_differences(y_data)

    if save_target_data:
        # save target data (only in the first time)
        y_data.to_csv("../data/Target_data.csv", index=False)
        print("Saved data to ../data/Target_data.csv")