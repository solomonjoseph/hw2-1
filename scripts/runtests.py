#%%
import os
import math
import numpy as np

TOTAL_CORES = 68
NUM_PARTICLES = 1.5e6
OUT_DIR = "../output"

#%%
def space(max_val, base=2, num_tests=-1):
    if num_tests == -1:
        num_tests = int(math.log(TOTAL_CORES, base)) + 1
    
    return [int(x) for x in np.logspace(0, math.log(max_val, base), num_tests, base=base)]

# print(space(TOTAL_CORES))
# print(space(NUM_PARTICLES))
sp = space(TOTAL_CORES)
print([int(NUM_PARTICLES * x / TOTAL_CORES) for x in sp])
print(sp)

#%%
os.chdir("build")
# os.system("ls")
# os.system("echo hi >> ../output/test.txt")

def strong_scaling(base=2, num_tests=-1):
    if num_tests == -1:
        num_tests = int(math.log(TOTAL_CORES, base)) + 1

    dest = OUT_DIR + "/strong_scaling.out"
    print(f"dest: {dest}")
    os.system("export OMP_PLACES=cores")
    os.system("export OMP_PROC_BIND=spread")
    os.system(f"echo NEW TEST >> {dest}")
    # os.system("echo NEW TEST ree")
    
    for x in np.logspace(0, math.log(TOTAL_CORES, base), num_tests, base=base):
        cores = int(x)
        print(f"running test on {cores} cores")
        os.system(f"echo cores: {cores} >> {dest}")
        os.system(f"export OMP_NUM_THREADS={cores}")
        os.system(f"./openmp -s 1 -n {NUM_PARTICLES} >> {dest}")
        print("done running test\n")


def weak_scaling(base=2, num_tests=-1):
    if num_tests == -1:
        num_tests = int(math.log(TOTAL_CORES, base)) + 1

    dest = OUT_DIR + "/weak_scaling.out"
    print(f"dest: {dest}")
    os.system("export OMP_PLACES=cores")
    os.system("export OMP_PROC_BIND=spread")
    os.system(f"echo NEW TEST >> {dest}")
    # os.system("echo NEW TEST ree")
    
    for x in np.logspace(0, math.log(TOTAL_CORES, base), num_tests, base=base):
        cores = int(x)
        print(f"running test on {cores} cores")
        os.system(f"echo cores: {cores} >> {dest}")
        os.system(f"export OMP_NUM_THREADS={cores}")
        os.system("echo $OMP_NUM_THREADS")
        os.system(f"./openmp -s 1 -n {NUM_PARTICLES} >> {dest}")
        print("done running test\n")



# strong_scaling()

# np.logspace(0, np.log(68), 7, base=np.e)
# export OMP_NUM_THREADS=68
# export OMP_PLACES=cores
# export OMP_PROC_BIND=spread
