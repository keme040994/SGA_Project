# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: February, 2016.

# LIBRARIES
from business_logic.pnml_functions import *
from math import sqrt
import networkx as nx
import numpy as np


# =================
# GENERAL FUNCTIONS
# =================

# FUNCTION: view_population
def view_population(population):
    """
    Shows the genes of every chromosome in a given population.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
    """
    for i in range(0, len(population)):
        print(str(population[i]) + "\n")


# FUNCTION: likelihood_result_calculator
def likelihood_result_calculator(population, likelihood_function, rep):
    """
    Given a population and the replicates, it adds the likelihood result to every chromosome.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
        likelihood_function: INT
            Indicates the likelihood function that is going to be used
                1 = cotemporal, 2 = next_step_one, 3 = next_step_one_two
        rep : LIST[rep1, rep2, rep3, ...]
            A repN is a biological data used to calc the likelihood result
    """
    if likelihood_function == 1:
        for i in range(0, len(population)):
            di_graph = convert_matrix_to_digraph(population[i].get_genes())
            population[i].set_likelihood_result(pnml_cotemporal(di_graph, rep))
    elif likelihood_function == 2:
        for i in range(0, len(population)):
            di_graph = convert_matrix_to_digraph(population[i].get_genes())
            population[i].set_likelihood_result(pnml_next_step_one(di_graph, rep))
    elif likelihood_function == 3:
        for i in range(0, len(population)):
            di_graph = convert_matrix_to_digraph(population[i].get_genes())
            population[i].set_likelihood_result(pnml_next_step_one_two(di_graph, rep))


# FUNCTION: relative_likelihood_result_sorting
def relative_likelihood_result_sorting(population):
    """
    Sorts a population given, in descending order, based on the relative likelihood result of every chromosome.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
    """
    population.sort(key=lambda x: x.relative_likelihood_result, reverse=True)


# FUNCTION: relative_likelihood_result_calculator
def relative_likelihood_result_calculator(population):
    """
    Given a population, this function calculates the relative likelihood result (it is a way to make
    the likelihood result bigger) to every chromosome.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
    """
    total = sum_likelihood_result(population)
    for i in range(0, len(population)):
        population[i].set_relative_likelihood_result(population[i].get_likelihood_result()/total)


# FUNCTION: sum_likelihood_result
def sum_likelihood_result(population):
    """
    Sums all the likelihood results values on a population.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        INT
            The sum of all likelihood results on a population
    """
    total = 0
    for i in range(0, len(population)):
        total += population[i].get_likelihood_result()
    return total


# FUNCTION: convert_matrix_to_digraph
def convert_matrix_to_digraph(matrix):
    """
    Converts a python's matrix into a networkx's digraph.

    Args:
        matrix : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes

    Returns:
        nx.DiGraph()
            The converted matrix into a digraph
    """
    genes = np.asmatrix(matrix)
    di_graph = nx.to_networkx_graph(genes, create_using=nx.DiGraph())
    return di_graph


# FUNCTION: convert_digraph_to_matrix
def convert_digraph_to_matrix(digraph):
    """
    Converts a networkx's digraph into a python's matrix.

    Args:
        digraph : nx.DiGraph()
            It's a networkx's digraph representing genes

    Returns:
        MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            The converted diagraph into a mtrix
    """
    return [[int(i) for i in j] for j in nx.to_numpy_matrix(digraph).tolist()]


# FUNCTION: transpose_matrix
def transpose_matrix(matrix):
    """
    Given a matrix, this function returns its transpose.

    Args:
        matrix : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes

    Returns:
        MATRIX[[INT, INT, ...], [INT, INT, ...], ...]

    """
    t_matrix = []
    for i in range(0, len(matrix)):
        t_matrix.append(extract_column(matrix, i))
    return t_matrix


# FUNCTION: extract_column
def extract_column(matrix, column_number):  # column_number:
    """
    An auxiliary function that extracts a column given a matrix.

    Args:
        matrix : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes
        column_number : INT
            A number between 0 to N-1

    Returns:
        LIST[INT, INT, INT, ...]
            A list that contains the elements of the column
    """
    column = []
    for i in range(0, len(matrix)):
        column.append(matrix[i][column_number])
    return column


# FUNCTION: select_uniques_chromosomes
def select_uniques_chromosomes(population):
    """
    Returns a filtered population with only unique elements. The comparison criteria is the 'relative likelihood
    result'.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        LIST[Chromosome(), Chromosome(), ...]
            A new list filled with 'Chromosome' objects
    """
    unique_population = []
    for i in range(0, len(population)):
        if find_element_on_list(population[i].get_relative_likelihood_result(), unique_population):
            unique_population.append(population[i])
    return unique_population


# FUNCTION: find_element_on_list
def find_element_on_list(relative_likelihood_result, unique_population):
    """
    An auxiliary function that returns a boolean value if the element exists on the 'unique_population' parameter.

    Args:
        relative_likelihood_result : FLOAT
            A number representing the value of certain chromosome's relative likelihood result
        unique_population: LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects where each one is unique

    Returns:
        BOOLEAN
            True it the value is found or False if not
    """
    for i in range(0, len(unique_population)):
        # based on the 'relative likelihood  result'
        if relative_likelihood_result == unique_population[i].get_relative_likelihood_result():
            return False
    return True


# FUNCTION: fitness_average
def fitness_average(population):
    fitness_avg = 0
    for item in population:
        fitness_avg += item.get_fitness()
    return fitness_avg/len(population)


# FUNCTION: fitness_variance
def fitness_variance(population):
    fitness_avg = fitness_average(population)
    fitness_var = 0
    for i in range(0, len(population)):
        fitness_var += (population[i].get_fitness()-fitness_avg)**2
    return fitness_var/len(population)


# FUNCTION: rlr_average
def rlr_average(population):
    rlr_avg = 0
    for item in population:
        rlr_avg += item.get_relative_likelihood_result()
    return rlr_avg/len(population)


# FUNCTION: rlr_variance
def rlr_variance(population):
    rlr_avg = rlr_average(population)
    rlr_var = 0
    for i in range(0, len(population)):
        rlr_var += (population[i].get_relative_likelihood_result()-rlr_avg)**2
    return rlr_var/len(population)
