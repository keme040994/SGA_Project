# Created by: Dr. David John & Kenneth Meza.
# Created at: March, 2016.
# Updated at: March, 2016.


# ===========================
# AMALGAMATED MODEL FUNCTIONS
# ===========================

# FUNCTION: amalgamated_model_creator
def amalgamated_model_creator(composite_models_generated, cant_genes):
    """
    Creates the amalgamated model, taking a group of composite models and summing all the matrices
    values into one new matrix, and dividing every value in the new matrix by the amount of composite
    models summed.

    Args:
        composite_models_generated : LIST[Composite_Model, Composite_Model, ...]
            A list containing composite models, represented in matrix of floats
        cant_genes : INT
            The amount of genes for each chromosome

    Returns:
        MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            The amalgamated model
    """
    amalgamated_model = [[0 for i in range(cant_genes)] for i in range(cant_genes)]
    for composite_model in composite_models_generated:
        for i in range(0, cant_genes):
            for j in range(0, cant_genes):
                amalgamated_model[i][j] += composite_model[i][j]
    return element_divisor(amalgamated_model, len(composite_models_generated))


# FUNCTION: element_divisor
def element_divisor(matrix, value):
    """
    Divides every element on a matrix by a given value.

    Args:
        matrix : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            A given matrix
        value : INT
            A given number

    Returns:
        matrix : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            The matrix divided by the given value
    """
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            matrix[i][j] /= value
    return matrix
