import pandas as pd


# function that reads from a text file dataframe of genomic data
def read_array(file_name):
    file = open(file_name, "r")  # first open the txt file
    file_content = file.readlines()

    # the first 2 lines are irrelevant for the dataframe
    header = file_content[2].split()  # the 3rd line is the columns names
    data = []
    for line in file_content[3:]:
        row = line.split()
        data.append(row)

    df = pd.DataFrame(data=data, columns=header)

    return df


# function to find the row that contains a specific gene by its name
def find_gene_name(df, gene_name, column):
    row_idx = df.index[df[column] == gene_name]  # returns the index of the row but not as int
    return row_idx.tolist()[0]  # convert to list and return the first (and only) index



# function that adds 2 columns containing the information of the indices of the gene location: "startCoordinate", "endCoordinate"
def create_coordinates_cols(df, column_name="Location"):
    df[["startCoordinate", "endCoordinate"]] = df[column_name].str.split('\.\.', expand=True)  # split the column by ".." between the 2 indices
    df[["startCoordinate", "endCoordinate"]].astype(int)
    return df


# function to delete a specific column
def drop_column(df, column):
    df = df.drop(columns=[column])
    return df



# check
# gene_df = read_array("array.txt")
# print(gene_df.shape)
# # print(df.columns)
# # print(df.head(5))
# # print(df[-1:])  # last line
#
# gene_row_idx = find_gene_name(gene_df, 'YAL067W-A', 'Synonym')
# print(gene_row_idx)
#
# gene_df_with_coordinates = create_coordinates_cols(gene_df, "Location")
# print(gene_df_with_coordinates.head(5))
#
# df_final = drop_column(gene_df_with_coordinates, "Location")
# print(df_final.head(5))
