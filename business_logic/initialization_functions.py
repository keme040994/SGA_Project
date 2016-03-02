# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: February, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from business_logic.general_functions import *
from business_logic.mutation_functions import *
from business_logic.repair_functions import *
from business_logic.selection_functions import *
from random import randint


# ========================
# INITIALIZATION FUNCTIONS
# ========================

# FUNCTION: seed_population
def seed_population(pop_size, cant_genes, per_ones, likelihood_function, rep):
    part_a_percentage = [15, 30, 45, 60, 75, 90]
    part_b_percentage = [60, 70, 80]

    initial_population = []
    initial_population.extend(population_creator(pop_size, cant_genes, per_ones))
    initial_population = repair_population(initial_population)

    likelihood_result_calculator(initial_population, likelihood_function, rep)
    relative_likelihood_result_calculator(initial_population)
    relative_likelihood_result_sorting(initial_population)

    # Part A
    for percentage in part_a_percentage:
        unique_population = select_uniques_chromosomes(initial_population)

        keep_percentage = (len(initial_population)*percentage)//100
        complete_percentage = len(initial_population)-keep_percentage

        if len(unique_population) >= keep_percentage:
            new_population = unique_population[:keep_percentage]
        else:
            # there must be a fancier way to do this!
            new_population = unique_population

            index = 0
            while len(new_population) < keep_percentage:
                new_population.append(unique_population[index])
                if index >= len(unique_population):
                    index = 0
                else:
                    index += 1

        new_population.extend(population_creator(complete_percentage, cant_genes, per_ones))

        initial_population = new_population
        initial_population = repair_population(initial_population)

        likelihood_result_calculator(initial_population, likelihood_function, rep)
        relative_likelihood_result_calculator(initial_population)
        relative_likelihood_result_sorting(initial_population)

    # Part B
    for percentage in part_b_percentage:
        unique_population = select_uniques_chromosomes(initial_population)

        keep_percentage = (len(initial_population)*percentage)//100
        complete_percentage = len(initial_population)-keep_percentage

        if len(unique_population) >= keep_percentage:
            new_population = unique_population[:keep_percentage]
        else:
            # there must be a fancier way to do this!
            new_population = unique_population

            index = 0
            while len(new_population) < keep_percentage:
                new_population.append(unique_population[index])
                if index >= len(unique_population):
                    index = 0
                else:
                    index += 1

        aux_population = []
        for i in range(0, complete_percentage):
            random_index = randint(0, len(initial_population)-1)
            aux_population.append(initial_population[random_index])
        mutation_function(aux_population, 0.5, 3)  # this could be (more) aggressive!

        initial_population = new_population
        initial_population.extend(aux_population)
        initial_population = repair_population(initial_population)

        likelihood_result_calculator(initial_population, likelihood_function, rep)
        relative_likelihood_result_calculator(initial_population)
        relative_likelihood_result_sorting(initial_population)
    return initial_population


# FUNCTION population_creator
def population_creator(pop_size, cant_genes, per_ones):
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
    population = []
    for i in range(0, pop_size):
        population.append(Chromosome(create_random_genes(cant_genes, per_ones)))
    return population


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
