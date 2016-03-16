# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: March, 2016.

# LIBRARIES
from business_logic.amalgamated_model_functions import *
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
    per_filter_cm = 0.7  # For filtering values in the composite model
    per_filter_am = 0.7  # For filtering values in the amalgamated model
    cant_matings = 500  # Amount of matings (loops)
    cant_composite_model = 10  # Amount of composite models fro creation

    filter_likelihood_selection(likelihood_function)

    # The variable 'rep' can be found inside 'data_logic.data_input.py'. Refers to replications.
    # rep = data_switcher(data.a_switch_log, data.a_switch_zscore, rep)  # log & zscore transforms

    composite_models_generated = []
    for i in range(0, cant_composite_model):
        print("* GENERATING COMPOSITE MODEL " + str(i+1))

        # Creation of the initial population using "seeding" method. See the function documentation
        print("* Creating the initial population...")
        current_population = seed_population(pop_size, cant_genes, per_ones, likelihood_function, data.rep)
        print("* Initial Population created, having " + str(len(select_uniques_chromosomes(current_population))) +
              " unique chromosomes.")

        for j in range(0, cant_matings):
            print("* Working on generation " + str(j+1) + "...")
            # 1) Given a sorted population, this function add ranks used on selection
            ranked_selection_calculator(current_population)

            # 2) Creation of the new population, by doing the selection process
            new_population = selection_function(current_population, per_elitism, selection_prop)

            # 3) Application of the mutation function to the population
            mutation_function(new_population, mutation_prop, cant_mutations)

            # 4) Application of the repairing function to the population
            current_population = repair_population(new_population)

            # 5) The applying of the likelihood on every chromosome to obtain the 'likelihood result'
            likelihood_result_calculator(current_population, likelihood_function, data.rep)

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
        likelihood_result_calculator(current_population, likelihood_function, data.rep)
        relative_likelihood_result_calculator(current_population)
        relative_likelihood_result_sorting(current_population)
        ranked_selection_calculator(current_population)

        # Creating, storing and displaying of the composite model
        print("* Creating the composite model...")
        composite_model = composite_model_creator(current_population, cant_genes)
        print("\tcreated.\n")

        if likelihood_function == 1:
            composite_model = sum_matrix(composite_model, transpose_matrix(composite_model))

        composite_models_generated.append(composite_model)
        view_model(composite_model, "COMPOSITE MODEL (ROUNDED TO 3 DECIMALS)")

        if likelihood_function == 1:
            create_model_image_cotemporal(convert_model_to_digraph(composite_model, per_filter_cm, data.a_protein_names),
                                          data.a_name + " - " + str(i+1))
        else:
            create_model_image(convert_model_to_digraph(composite_model, per_filter_cm, data.a_protein_names),
                               data.a_name + " - " + str(i+1))

        print("\n* Creating the image...\n\tcreated.")
        write_file(data.a_name + " - " + str(i+1), str(composite_model))
        print("* Creating the file...\n\tcreated.")
        print("* Done.\n")

    print("...................................\n")

    print("* GENERATING AMALGAMATED MODEL")
    # Creating, storing and displaying of the amalgamated model
    print("* Creating the amalgamated model...")
    amalgamated_model = amalgamated_model_creator(composite_models_generated, cant_genes)
    print("\tcreated.\n")
    view_model(amalgamated_model, "AMALGAMATED MODEL (ROUNDED TO 3 DECIMALS)")

    if likelihood_function == 1:
        create_model_image_cotemporal(convert_model_to_digraph(amalgamated_model, per_filter_am, data.a_protein_names),
                                      data.a_name + " (AMALGAMATED MODEL)")
    else:
        create_model_image(convert_model_to_digraph(amalgamated_model, per_filter_am, data.a_protein_names),
                           data.a_name + " (AMALGAMATED MODEL)")

    print("\n* Creating the image...\n\tcreated.")
    write_file(data.a_name + " (AMALGAMATED MODEL)", str(amalgamated_model))
    print("* Creating the file...\n\tcreated.")
    print("* Done.\n")
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

if __name__ == "__main__":
    main()
