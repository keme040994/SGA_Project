# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: April, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from business_logic.general_functions import *
from copy import deepcopy
from random import randint
import random


# ===================
# SELECTION FUNCTIONS
# ===================

# FUNCTION: fitness_calculator
def fitness_calculator(ordered_population):
    """
    Calculates the 'fitness' that will be used for doing a rank based selection.

    Args:
        ordered_population : LIST[Chromosome]
            A list filled with 'Chromosome' objects sorted by the likelihood result
    """
    size = len(ordered_population)
    for i in range(0, size):
        ordered_population[i].fitness = (2*(size+1-(i+1)))/(size*(size+1))


# FUNCTION: selection_function
def selection_function(population, num_survivors, selection_prop):
    """
    Creates the new generation of chromosomes by doing the crossover on every two parents. The elitism happens here.

    Args:
        population : LIST[Chromosome]
            A list filled with 'Chromosome' objects
        num_survivors : INT
            The number of survivors, based on pre-calculated data using the percentage of elitism
        selection_prop : FLOAT
            The desired probability of selection
    Returns:
        LIST[Chromosome]
            A list filled with 'Chromosome' objects, containing the new population
    """
    match_list = match_list_creator(population)

    # Elitism
    unique_population = select_uniques_chromosomes(population)
    if len(unique_population) > num_survivors:
        new_population = unique_population[:num_survivors]
    else:
        new_population = unique_population
        num_survivors -= len(unique_population)
        if num_survivors % 2 != 0:
            num_survivors += 1

    # Selection
    for i in range(0, len(population)-num_survivors, 2):
        random1 = random.random()
        random2 = random.random()
        children = crossover_function(population[match_list_finder(match_list, random1)],
                                      population[match_list_finder(match_list, random2)],
                                      selection_prop)
        new_population.append(children[0])
        new_population.append(children[1])
    return new_population


# FUNCTION: calc_num_survivors
def calc_num_survivors(len_population, per_elitism):
    """
    Calculates the number of chromosomes who are going to survive.

    Args:
        len_population : INT
            A number that represents the size of the population
        per_elitism : INT
            The desired percentage of elitism

    Returns:
        INT
            A number indicating the number of survivors
    """
    num_survivors = (len_population*per_elitism)//100
    if num_survivors % 2 != 0:
        num_survivors += 1
    return num_survivors


# FUNCTION: crossover_function
def crossover_function(parent_a, parent_b, selection_prop):
    """
    Does the crossover given two parents and returns two offsprings. If the selection probability doesn't get matched,
    the two parents will be the new offsprings.

    Args:
        parent_a : Chromosome
            An object 'Chromosome' that will be used for mating
        parent_b : Chromosome
            An object 'Chromosome' that will be used for mating
        selection_prop : FLOAT
            The desired probability of selection

    Returns:
        TUPLE(Chromosome, Chromosome)
            A tuple formed by two offsprings, represented by a 'Chromosome' object
    """
    genes_a = deepcopy(parent_a.get_genes())
    genes_b = deepcopy(parent_b.get_genes())
    match_prop = random.random()
    if match_prop <= selection_prop:
        column_number = randint(0, len(genes_a)-1)
        for i in range(0, len(genes_a)):
            genes_a[i][column_number], genes_b[i][column_number] = genes_b[i][column_number], genes_a[i][column_number]
    return Chromosome(genes_a), Chromosome(genes_b)


# FUNCTION: match_list_creator
def match_list_creator(population):
    """
    Given a population this function returns a list of probabilities that will be used in the selection process.

    Args:
        population : LIST[Chromosome]
            A list filled with 'Chromosome' objects

    Returns:
        LIST[FLOAT]
            A list of probabilities in the range 0 to 1
    """
    match_list = []
    sum_ = 0
    for i in range(0, len(population)):
        sum_ += population[i].get_fitness()
        match_list.append(sum_)
    return match_list


# FUNCTION: match_list_finder
def match_list_finder(match_list, random_num):
    """
    Using the list created in the function 'match_list_creator', this function returns and index that represents the
    position of a parent in the ordered population that is being mating.

    Args:
        match_list : LIST[FLOAT]
            A list of probabilities in the range 0 to 1
        random_num: FLOAT
            A random generated number that is created in the selection function

    Returns:
        INT
            An index representing a parent position
    """
    for i in range(0, len(match_list)):
        if random_num <= match_list[i]:
            return i
    return len(match_list)
