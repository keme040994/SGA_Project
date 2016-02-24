# Created by: Dr. David John & Kenneth Meza.
# Created at: February, 2016.
# Updated at: February, 2016.

# LIBRARIES
import pydotplus as pydot


# =======
#   GUI
# =======

# FUNCTION: show_matrix_gui
def show_matrix_gui(di_graph, name):
    """
    Using a networkx's DAG (created based on the composite model), this function creates an image of the result of
    the Simple Genetic Algorithm using 'pydotplus' library.

    Args:
        di_graph : nx.DiGraph()
            It's a networkx's digraph representing the composite model
        name : STRING
            The name given in the variable 'a_name' on 'data_logic.data_input.py'

    Returns:
        PNG IMAGE
            An image with the result of the problem, saved on the root file
    """
    display_graph = pydot.Dot(graph_type='digraph')

    nodes = di_graph.nodes()
    for node in nodes:
        display_graph.add_node(pydot.Node(node))

    edges = di_graph.edges()
    for edge in edges:
        if di_graph[edge[0]][edge[1]]['weight'] == 1:
            display_graph.add_edge(pydot.Edge(edge[0], edge[1],
                                              label=str(di_graph[edge[0]][edge[1]]['weight']),
                                              color="green"))
        else:
            display_graph.add_edge(pydot.Edge(edge[0], edge[1],
                                              label=str(di_graph[edge[0]][edge[1]]['weight']),
                                              color="red"))
    display_graph.write_png("../" + name + ".png")
