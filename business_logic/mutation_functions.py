# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: February, 2016.

# LIBRARIES
from random import randint
import random


# ==================
# MUTATION FUNCTIONS
# ==================

# FUNCTION: mutation_function
def mutation_function(population, mutation_prob):
    """
    Given a probability of mutation, this function applies a mutation (change of a bit in the genes) to every
    chromosome in the population.

    Args:
        population : LIST[Chromosome]
            A list filled with 'Chromosome' objects
        mutation_prob : FLOAT
            The desired probability of mutation
    """
    for i in range(0, len(population)):
        if random.random() <= mutation_prob:
            random1 = randint(0, len(population[0].get_genes())-1)
            random2 = randint(0, len(population[0].get_genes())-1)
            while random1 == random2:
                random1 = randint(0, len(population[0].get_genes())-1)
                random2 = randint(0, len(population[0].get_genes())-1)
            population[i].bit_changer(random1, random2)
