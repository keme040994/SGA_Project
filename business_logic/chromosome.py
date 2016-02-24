# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: February, 2016.


# CLASS: Chromosome
class Chromosome:

    # CONSTRUCTOR
    def __init__(self, genes):
        self.genes = genes
        self.likelihood_result = None
        self.relative_likelihood_result = None
        self.fitness = None

    # ACCESSOR METHODS
    def get_genes(self):
        return self.genes

    def set_genes(self, genes):
        self.genes = genes

    def get_likelihood_result(self):
        return self.likelihood_result

    def set_likelihood_result(self, likelihood_result):
        self.likelihood_result = likelihood_result

    def get_relative_likelihood_result(self):
        return self.relative_likelihood_result

    def set_relative_likelihood_result(self, relative_likelihood_result):
        self.relative_likelihood_result = relative_likelihood_result

    def get_fitness(self):
        return self.fitness

    def set_fitness(self, fitness):
        self.fitness = fitness

    # METHODS
    def bit_changer(self, i, j):
        self.genes[i][j] = (1 - self.genes[i][j])

    def counter_ones(self):
        num_ones = 0
        for i in range(0, len(self.genes)):
            for j in range(0, len(self.genes)):
                num_ones += self.genes[i][j]
        return num_ones

    def __repr__(self):
        msg = ""
        for i in range(0, len(self.genes)):
            msg += str(self.genes[i]) + "\n"
        msg += "Likelihood Result: " + str(self.likelihood_result) + "\n"
        msg += "Relative Likelihood Result: " + str(self.relative_likelihood_result) + "\n"
        msg += "Fitness: " + str(self.fitness) + "\n"
        msg += "Num. Ones: " + str(self.counter_ones())
        return msg
