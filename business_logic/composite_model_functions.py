# Created by: Dr. David John & Kenneth Meza.
# Created at: February, 2016.
# Updated at: February, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from business_logic.general_functions import *
import networkx as nx


# =========================
# COMPOSITE MODEL FUNCTIONS
# =========================

# FUNCTION: composite_model_creator
def composite_model_creator(population, cant_genes):
    """
    Creates the composite model by multiplying every chromosome to its fitness and sum them all in a new matrix.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
        cant_genes : INT
            The desired amount of genes for each chromosome

    Returns:
        MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            The composite model
    """
    new_population = fitness_x_genes(population)
    composite_model = [[0 for i in range(cant_genes)] for i in range(cant_genes)]
    for i in range(0, len(new_population)):
        genes = population[i].get_genes()
        for j in range(0, len(genes)):
            for k in range(0, len(genes)):
                composite_model[j][k] += genes[j][k]
    return round_matrix(composite_model, 3)


# FUNCTION: round_matrix
def round_matrix(composite_model, cant_decimals):
    """
    Returns a composite model's matrix, where every element has being rounded by 3.

    Args:
        composite_model : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            A given composite model
        cant_decimals : INT
            The amount of decimals desired for the rounding

    Returns:
        A rounded composite model' matrix
    """
    for i in range(0, len(composite_model)):
        for j in range(0, len(composite_model)):
            composite_model[i][j] = round(composite_model[i][j], cant_decimals)
    return composite_model


# FUNCTION: fitness_x_genes
def fitness_x_genes(population):
    """
    This is an auxiliary function that multiply the fitness by the genes on every chromosome of a given population.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        LIST[Chromosome(), Chromosome(), ...]
            The new population with the genes multiplied by the fitness
    """
    new_population = []
    for i in range(0, len(population)):
        genes = population[i].get_genes()
        fitness = population[i].get_fitness()
        for j in range(0, len(genes)):
            for k in range(0, len(genes)):
                genes[j][k] *= fitness
        new_population.append(Chromosome(genes))
    return new_population


# FUNCTION: view_composite_model
def view_composite_model(composite_model):
    """
    Displays a composite model and all its content

    Args:
        composite_model : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            A given composite model
    """
    print("COMPOSITE MODEL (ROUNDED TO 3 DECIMALS)")
    for row in composite_model:
        for val in row:
            print("{:10.3f}".format(val), end="")
        print()


# FUNCTION: convert_composite_model_to_digraph
def convert_composite_model_to_digraph(composite_model, a_protein_names):
    """
    A special conversion function where given a composite model and the proteins names it combines all in a new digraph
    that will contains the information necessary for displaying the graph.

    Args:
        composite_model : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            A given composite model
        a_protein_names : LIST[STRING, STRING, STRING, ...]
            The name of all the proteins given in the variable 'a_protein_names' on 'data_logic.data_input.py'

    Returns:
        nx.DiGraph()
            The converted composite model into the true digraph format for displaying
    """
    cm_digraph = convert_matrix_to_digraph(composite_model)
    di_graph = nx.DiGraph()
    di_graph.add_nodes_from(a_protein_names)
    for i in range(0, len(cm_digraph.edges())):
        edge = cm_digraph.edges()[i]
        if cm_digraph[edge[0]][edge[1]]['weight'] >= 0.7:
            di_graph.add_edge(a_protein_names[edge[0]], a_protein_names[edge[1]],
                              weight=cm_digraph[edge[0]][edge[1]]['weight'])
    return di_graph
