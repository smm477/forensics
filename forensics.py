import random
import subprocess
import os
import numpy as np
import pandas
import sys
import os

it = sys.argv[1]

'''Creates a random test set of specified size
Arguments:
    pop_file: a txt file with individuals from a specific population to sample from
    size: an integer specifying desired size of test set
        (must be <= number of individuals in given population)
Returns:
    indiv_list: a list of individual id's in test_set
'''
def create_test_set(pop_file, size, i):
    # create list of individuals from pop_file
    with open(pop_file) as f:
        indiv_list = f.readlines()

    # randomly chose individuals (without replacement) & create list of selected individuals
    test_set = []
    for ind in range(size):
        rand = random.randint(0, len(indiv_list)-1)
        test_set.append(indiv_list[rand])
        indiv_list.remove(indiv_list[rand])

    # turn test set list into txt file with an ID on each line
    filename = str(i) + 'test_ids.txt'
    with open(filename, 'w') as test_ids:
        test_ids.writelines(test_set)

    return test_set

'''Creates reference panel from given list/population of individuals of specified size
Arguments:
    pop_file: a txt file with individuals to sample from & test_set individauls
    size: an integer specifying desired size of reference panel
        (must be <= number of individuals in given file - number in test_set)
    test_set:a list of individual id's in the previously created test_set
'''
def create_ref_panel(pop_file, size, test_set, i):
    #creates list of individuals from pop_file
    with open(pop_file) as f:
        indiv_list = f.readlines()

    #removes individuals in the test_set from list of inidivuals to sample from
    for indiv in indiv_list:
        if indiv in test_set: indiv_list.remove(indiv)

    #randomly chose individuals (without replacement) & create list of selected individuals
    ref_panel = []
    for ind in range(size):
        rand = random.randint(0, len(indiv_list)-1)
        ref_panel.append(indiv_list[rand])
        indiv_list.remove(indiv_list[rand])

    #turn reference panel list into txt file with an ID on each line
    filename = str(i) + 'ref_panel_ids.txt'
    with open(filename, 'w') as ref_panel_ids:
        ref_panel_ids.writelines(ref_panel) 

'''Calculates r^2 value for imputation imputation_accuracy
Arguments:
    imp_file: file of imputed genotypes (.FORMAT)
    true_file: file of true genotypes (.FORMAT)
Returns:
    r_sq: the r^2 value (correlation coefficient) between the imputed and true
            genotypes for the individuals included in
'''
def imputation_accuracy(imp_file, true_file):
    input_imp = list(pandas.read_csv(imp_file, sep= '\t').loc[0].values.flatten())[2:]
    input_true = list(pandas.read_csv(true_file, sep= '\t').loc[0].values.flatten())[2:]

    list_imp, list_true = [], []
    for element in input_imp:
        num1 = int(element.partition('|')[0])
        num2 = int(element.partition('|')[2])
        list_imp.append([num1, num2])

    for element in input_true:
        num1 = int(element.partition('/')[0])
        num2 = int(element.partition('/')[2])
        list_true.append([num1, num2])

    imp_D, true_D = [], []
    for pair in list_imp:
        imp_D.append(pair[0] + pair[1])
    for pair in list_true:
        true_D.append(pair[0] + pair[1])

    cov = np.cov(imp_D, true_D)
    r_sq = cov[0][1]*cov[0][1] / (cov[0][0]*cov[1][1])
    return r_sq

'''Computes r^2 for each loci in test run
Arguments:
    filepath: the filepath for the directory of output files of 1 iteration
Returns:
    my_dict: a dictionary that maps imputation accuracy to CODIS loci for the given directory
'''
def files_to_imp_accuracy(filepath):
    codis = ['CSF1PO', 'D13S317', 'D18S51', 'D3S1358', 'D5S818', 'D7S820', 'D8S1179',
            'FGA', 'TH01', 'TPOX', 'vWA', 'D1S1656', 'D2S441', 'D2S1338', 'D10S1248', 'D12S391', 'D19S433', 'D22S1045']
    list_of_files = os.listdir(filepath)
    my_dict = {}
    for loci in codis:
        true = filepath+ '/' + str(it) +  '/' + loci+ 'codistest.GT.FORMAT'
        imp = filepath + '/' + str(it) + '/' + loci + '.GT.FORMAT'
        my_dict[loci] = imputation_accuracy(imp, true)
    return my_dict

# list of populations
# populations = ['ACB', 'CEU', 'FIN', 'IBS', 'LWK', 'PJL', 'YRI', 'ASW', 'CHB', 
#  'GBR', 'ITU', 'MSL', 'PUR', 'BEB', 'CLM', 'GIH', 'JPT', 'MXL', 'STU', 'CDX', 
#  'ESN', 'GWD', 'KHV', 'PEL', 'TSI']

# populations = ['GWD', 'PUR', 'CHS', 'IBS', 'GIH']

populations = ['GWD', 'PUR']

# ref panel sizes
# size = [10, 20, 30, 40, 50]
size = [10, 20]

# Computational Pipeline
def main():
  codis = ['CSF1PO', 'D13S317', 'D18S51', 'D3S1358', 'D5S818', 'D7S820', 'D8S1179',
            'FGA', 'TH01', 'TPOX', 'vWA', 'D1S1656', 'D2S441', 'D2S1338', 'D10S1248', 'D12S391', 'D19S433', 'D22S1045']
  columns = {'population': [], 'ref panel size':[] }
  for loci in codis:
    columns[loci] = []
  df = pandas.DataFrame(columns)
  counter = 0
  for s in size: # vary reference panel size
    for pop in populations: # Create a test set for each population
      counter+=1
      # create folder for files
      dir = '/workdir/forensics/' + pop + str(s)
      popfile = '/workdir/forensics/poplists/' + pop + '.txt'
      test_set = create_test_set(popfile, 5, it)
      create_ref_panel(popfile, s, test_set, it)
      # run bash script
      test = '/workdir/forensics/' + str(it) + 'test_ids.txt'
      ref = '/workdir/forensics/' + str(it) + 'ref_panel_ids.txt'
      subprocess.run(['/workdir/forensics/bash.sh', test, ref, str(it)])
      # rename tmp folder and move into population folder
      name = dir + '/' + str(it)
      os.rename('/workdir/forensics/tmp'+str(it), name)
      os.rename(test, name + '/test_ids.txt')
      os.rename(ref, name + '/ref_panel.txt')
      # calculate imputation accuracy
      accuracy_dict = files_to_imp_accuracy(dir)
      df.at[counter, 'population'] = pop
      df.at[counter,'ref panel size'] = s
      for element in codis:
        df.at[counter][element] = accuracy_dict[element]
 
  df.to_csv('/workdir/forensics/csv/'+str(it)+'.csv')

if __name__ == "__main__":
  main()
