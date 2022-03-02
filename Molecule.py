#class for collecting and processing CMA benchmarking data, as well as generating LaTeX-formatted SI on a molecule by molecule basis.

class Molecule(object):

    def __init__(self,natty_obj,redundant_obj,ID):
        #nothing to init
        self.natty = natty_obj
        self.redundant = redundant_obj
        self.ID = ID
        #self.freq_data = freq_data
        #self.projection_mat = projection_mat
        #self.ted
        #self.LCIC 
    def run(self):
        print('Molecule Class!')
        print(self.ID)

         







