# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: February, 2016.

# LIBRARIES
from business_logic.composite_model_functions import *
from business_logic.initialization_functions import *
from business_logic.mutation_functions import *
from business_logic.repair_functions import *
from business_logic.selection_functions import *
from data_logic.data_functions import *
from presentation_logic.GUI import *
import data_logic.data_input as data


# ========
#   MAIN
# ========

# FUNCTION: main
def main():
    """ The main function program, having the structure of the SGA. """

    # Initial Variables
    pop_size = 40  # Population size
    cant_genes = 12  # Amount of genes for each chromosome
    per_ones = 3  # Percentage of ones based on matrix's size
    per_elitism = 10  # Percentage of elitism
    selection_prop = 0.5  # Probability of selection
    mutation_prop = 0.3  # Probability of mutation
    cant_matings = 100  # Amount of matings (loops)

    rep = [data.a1, data.a2, data.a3]  # this data comes from 'data_logic.data_input.py'
    # rep = data_switcher(data.a_switch_log, data.a_switch_zscore, rep)  # log & zscore transforms

    # Creation of the initial population
    current_population = initial_population_creator(pop_size, cant_genes, per_ones)
    current_population = repair_population(current_population)

    for i in range(0, cant_matings):
        # 1) The applying of the likelihood on every chromosome to obtain the 'likelihood result'
        likelihood_result_calculator(current_population, rep)

        # 2) Finding the 'relative likelihood result' for every chromosome
        relative_likelihood_result_calculator(current_population)

        # 3) Sorting of the population based on the 'likelihood result' of every chromosome
        relative_likelihood_result_sorting(current_population)

        # 4) Given a sorted population, this function add ranks used on selection
        ranked_selection_calculator(current_population)

        print("FITNESS AVERAGE: " + str(fitness_average(current_population)))
        print("FITNESS VARIANCE: " + str(fitness_variance(current_population)))

        print("RLR AVERAGE: " + str(rlr_average(current_population)))
        print("RLR VARIANCE: " + str(rlr_variance(current_population)))

        # 5) Creation of the new population, by doing the selection process
        new_population = selection_function(current_population, per_elitism, selection_prop)

        # 6) Application of the mutation function to the population
        mutation_function(new_population, mutation_prop)

        # 7) Application of the repairing function to the population
        current_population = repair_population(new_population)

        cont = 0
        for j in range(0, len(current_population)):
            cont += current_population[j].counter_ones()
        print("Generation " + str(i+1) + " - TOTAL NUM. ONES: " + str(cont))
        print("Average: " + str(cont/pop_size))
        print()

    # Removing from the population the repeated chromosomes
    likelihood_result_calculator(current_population, rep)
    relative_likelihood_result_calculator(current_population)
    relative_likelihood_result_sorting(current_population)
    ranked_selection_calculator(current_population)

    current_population = select_uniques_chromosomes(current_population)

    print("NUM. OF UNIQUE ONES: " + str(len(current_population)))

    # Recalculation of 'likelihood' values and fitness on the current population
    likelihood_result_calculator(current_population, rep)
    relative_likelihood_result_calculator(current_population)
    relative_likelihood_result_sorting(current_population)
    ranked_selection_calculator(current_population)

    # Creating and displaying of the composite model
    composite_model = composite_model_calculator(current_population, cant_genes)
    view_composite_model(composite_model)
    show_matrix_gui(convert_composite_model_to_digraph(composite_model, data.a_protein_names), data.a_name)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

if __name__ == "__main__":
    main()
