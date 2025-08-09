import pandas as pd
from typing import List

def validate_summits_df(df: pd.DataFrame, required_columns: List[str] = None) -> pd.DataFrame:
    """Validate and standardize a summits DataFrame.

    This function checks that the input DataFrame contains the required columns
    (by default: 'chrom', 'start', 'end'), ensures that the data types of these 
    columns are appropriate, and verifies that the start positions are less than or 
    equal to the end positions.

    Args:
        df (pd.DataFrame): Input DataFrame containing summit information.
        required_columns (List[str], optional): List of required column names.
            Defaults to ['chrom', 'start', 'end'].

    Returns:
        pd.DataFrame: A copy of the input DataFrame with standardized column names.

    Raises:
        ValueError: If any required column is missing, if the data types are not as expected,
                    or if some rows have start greater than end.
    """
    if required_columns is None:
        required_columns = ['chrom', 'start', 'end']
    
    # check columns
    missing_cols = set(required_columns) - set(df.columns)
    if missing_cols:
        raise ValueError(f"Summits DataFrame is missing required columns: {missing_cols}")
    
    # copy dataframe, avoiding to modify original dataframe
    df_valid = df.copy()
    df_valid.columns = [col.lower() for col in df_valid.columns]

    # check if "chrom" is type of string type
    if not pd.api.types.is_string_dtype(df_valid['chrom']):
        raise ValueError("Column 'chrom' must be of string type")
    
    if not pd.api.types.is_numeric_dtype(df_valid['start']):
        raise ValueError("Column 'start' must be numeric.")
    
    if not pd.api.types.is_numeric_dtype(df_valid['end']):
        raise ValueError("Column 'end' must be numeric.")

    if (df_valid['start'] > df_valid['end']).any():
        raise ValueError("Some rows have start > end.")
    
    return df_valid