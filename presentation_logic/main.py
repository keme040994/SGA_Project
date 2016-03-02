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
from presentation_logic.output_functions import *
import data_logic.data_input as data


# ========
#   MAIN
# ========

# FUNCTION: main
def main():
    """ The main function program, having the structure of the SGA. """

    # Initial Variables
    pop_size = 250  # Population size
    cant_genes = 12  # Amount of genes for each chromosome
    per_ones = 5  # Percentage of ones based on matrix's size
    likelihood_function = 1  # 1 = cotemporal, 2 = next_step_one, 3 = next_step_one_two
    per_elitism = 10  # Percentage of elitism
    selection_prop = 0.05  # Probability of selection
    mutation_prop = 0.3  # Probability of mutation
    cant_mutations = 1  # Amount of mutations per chromosome
    cant_matings = 25  # Amount of matings (loops)

    rep = [data.a1, data.a2, data.a3]  # this data comes from 'data_logic.data_input.py'
    # rep = data_switcher(data.a_switch_log, data.a_switch_zscore, rep)  # log & zscore transforms

    # Creation of the initial population using "seeding" method. See the function documentation
    print("* Creating the initial population...")
    current_population = seed_population(pop_size, cant_genes, per_ones, likelihood_function, rep)
    print("* Initial Population created, having " + str(len(select_uniques_chromosomes(current_population))) +
          " unique chromosomes.")

    for i in range(0, cant_matings):
        print("* Working on generation " + str(i+1) + "...")
        # 1) Given a sorted population, this function add ranks used on selection
        ranked_selection_calculator(current_population)

        # 2) Creation of the new population, by doing the selection process
        new_population = selection_function(current_population, per_elitism, selection_prop)

        # 3) Application of the mutation function to the population
        mutation_function(new_population, mutation_prop, cant_mutations)

        # 4) Application of the repairing function to the population
        current_population = repair_population(new_population)

        # 5) The applying of the likelihood on every chromosome to obtain the 'likelihood result'
        likelihood_result_calculator(current_population, likelihood_function, rep)

        # 6) Finding the 'relative likelihood result' for every chromosome
        relative_likelihood_result_calculator(current_population)

        # 7) Sorting of the population based on the 'likelihood result' of every chromosome
        relative_likelihood_result_sorting(current_population)

        print("\tcreated.")

    # Removing from the population the repeated chromosomes
    ranked_selection_calculator(current_population)
    current_population = select_uniques_chromosomes(current_population)
    print("* The last generation ends with " + str(len(current_population)) + " unique chromosomes.")

    # Recalculation of 'likelihood' values and fitness on the current population
    likelihood_result_calculator(current_population, likelihood_function, rep)
    relative_likelihood_result_calculator(current_population)
    relative_likelihood_result_sorting(current_population)
    ranked_selection_calculator(current_population)

    # Creating and displaying of the composite model
    print("* Creating the composite model...")
    composite_model = composite_model_creator(current_population, cant_genes)
    print("\tcreated.\n")
    view_composite_model(composite_model)
    create_composite_model_image(convert_composite_model_to_digraph(composite_model, data.a_protein_names), data.a_name)
    print("\n* Creating the image...\n\tcreated.")
    write_file(data.a_name, str(composite_model))
    print("* Creating the file...\n\tcreated.")
    print("* Done.")
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

if __name__ == "__main__":
    main()
