# Adds and removes dummy atoms to fc.dat files in specified indicies
# python dummy.py filename [indicies] [indicies]
# use "-" and "+" to specify indicies and "," to separate them
# e.g. python dummy.py fc.dat -1,2 +10,11
# indicies can be in any order and there can be more than one specification

import sys 

import numpy as np
np.set_printoptions(precision=2, suppress=True, linewidth=300)

with open(sys.argv[1], 'r+') as f:
    lines = f.readlines()
    ndim = int(lines[0].split()[1])

    fcmat = []
    for line in lines[1:]:
        fcmat.append([float(i) for i in line.split()])
        
    fcnew = np.array(fcmat).reshape((ndim,ndim))  
    # print(fcnew)
    # print()
  
    # print(sys.argv[1:])

    sub, add = False, False
    sub_ind, add_ind = [],[]
    for arg in sys.argv[2:]:
        if "-" in arg:
            sub = True
            sub_ind += [int(i) for i in arg[1:].split(",")]
            # print(sub_ind)

        if "+" in arg:
            add = True
            add_ind += [int(i) for i in arg[1:].split(",")]
            # print(add_ind)

    sub_ind = sorted(sub_ind,reverse=True)
    add_ind = sorted(add_ind)

    if sub:
        for i in sub_ind:
            fcnew = np.delete(fcnew, slice((i-1)*3,i*3), axis=0)
            ndim -= 3
            fcnew = np.delete(fcnew, slice((i-1)*3,i*3), axis=1)
        # print(fcnew)
        # print()

    if add:
        for i in add_ind:
            fcnew = np.insert(fcnew, (i-1)*3, np.zeros((3,ndim)), axis=0)
            ndim += 3
            fcnew = np.insert(fcnew, [(i-1)*3], np.zeros((ndim,3)), axis=1)
        # print(fcnew)
        # print()

    fcnew = fcnew.reshape((-1, 3))

    new_txt = '    {0:<4}{1} \n'.format(int(ndim/3), ndim)
    for row in fcnew:
        new_txt += '       {0:>13.10f}       {1:>13.10f}       {2:>13.10f} \n'.format(row[0],row[1],row[2])

    f.seek(0)
    f.truncate()
    f.write(new_txt)
    
    f.close()
