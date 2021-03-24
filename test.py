from ArtinWedderburn import *

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def pass_fail(aw, defect_threshold=1.0e-10):
    if aw.total_defect < defect_threshold:
        print(bcolors.OKGREEN + "passed" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "failed" + bcolors.ENDC)

print("irreps for S2: ",end='')
aw2 = ArtinWedderburn(symmetric_group_algebra(2))
pass_fail(aw2)

print("irreps for S3: ",end='')
aw3 = ArtinWedderburn(symmetric_group_algebra(3))
pass_fail(aw3)

print("irreps for S4: ", end='')
aw4 = ArtinWedderburn(symmetric_group_algebra(4))
pass_fail(aw4)

print("irreps for sparse S2: ",end='')
aw2s = ArtinWedderburn(sparse_symmetric_group_algebra(2))
pass_fail(aw2)

print("irreps for sparse S3: ",end='')
aw3s = ArtinWedderburn(sparse_symmetric_group_algebra(3))
pass_fail(aw3)

print("irreps for sparse S4: ", end='')
aw4s = ArtinWedderburn(sparse_symmetric_group_algebra(4))
pass_fail(aw4)
