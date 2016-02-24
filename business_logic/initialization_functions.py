# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: February, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from random import randint


# ========================
# INITIALIZATION FUNCTIONS
# ========================

# FUNCTION initial_population_creator
def initial_population_creator(pop_size, cant_genes, per_ones):
    """
    Returns a list of chromosomes, created by random methods.

    Args:
        pop_size : INT
            The desired population size
        cant_genes : INT
            The desired amount of genes for each chromosome
        per_ones : INT
            The percentage of ones based on matrix's size

    Returns:
        LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
    """
    initial_population = []
    for i in range(0, pop_size):
        initial_population.append(Chromosome(create_random_genes(cant_genes, per_ones)))
    return initial_population


# FUNCTION: create_random_genes
def create_random_genes(cant_genes, per_ones):
    """
    Creates a random matrix that will be used as genes for an initial population's chromosome.

    Args:
        cant_genes : INT
            The desired amount of genes for each chromosome
        per_ones : INT
            The percentage of ones based on matrix's size

    Returns:
        MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix represents the genes and has a random number of ones based on the 'per_ones'
    """
    matrix = [[0 for i in range(cant_genes)] for i in range(cant_genes)]
    matrix_size = len(matrix)*len(matrix)
    cant_ones = randint(0, (matrix_size*per_ones)//100)
    for i in range(0, cant_ones):
        x = randint(0, len(matrix)-1)
        y = randint(0, len(matrix)-1)
        matrix[x][y] = 1
    return matrix
