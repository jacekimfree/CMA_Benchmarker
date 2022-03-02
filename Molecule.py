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

         







    def freq_diff(F1, F2):
        return F1 - F2
    
    def averages(F_diff):
        return np.average(F_diff), np.average(abs(F_diff)) 
    def stdev(F_diff):
        return np.std(F_diff)
    
    def return_path_base(jobpath):
        name = os.path.basename(os.path.normpath(jobpath)) 
        return name


    def build_dataframe(d,basename):
        df = pd.DataFrame(data=d)
        df_i = df.to_csv('out.csv',index=False)
        # adds a space in between each molecule. 
        df.loc[df.shape[0]] = [None, None, None, None, None, None, None] 
        df.loc[df.shape[0]] = [None, None, None, None, None, None, None] 
        df.loc[df.shape[0]] = [None, None, None, None, None, None, None] 
    
        #df.loc[df.shape[0]] = [None, None, 'Signed Avg.', averages(F1diff)[0], 'Signed Avg.', averages(F1diff)[0], averages(F1F2diff)[0]]
        #df.loc[df.shape[0]] = [None, None, 'Average   .', averages(F1diff)[1], 'Average    ', averages(F1diff)[1], averages(F1F2diff)[1]]
        #df.loc[df.shape[0]] = [None, None, 'Std. Dev  .', stdev(F1diff),    'Std. Dev.  ', stdev(F1diff), stdev(F1F2diff)]
 
        df.loc[0,('Molecule')] = str(basename) 
        #df_i: dataframe, individual for each child directory
        frame.append(df)
        return frame, df_i 

#     def grab_geometries(self):
#         try:
#             with open("zmat", "r") as f:
#                 output = f.readlines()


