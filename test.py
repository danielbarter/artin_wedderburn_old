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
s2 = load_algebra_from_file('./symmetric_groups/2')
aw2 = ArtinWedderburn(s2)
pass_fail(aw2)

print("irreps for S3: ",end='')
s3 = load_algebra_from_file('./symmetric_groups/3')
aw3 = ArtinWedderburn(s3)
pass_fail(aw3)


print("irreps for S4: ", end='')
s4 = load_algebra_from_file('./symmetric_groups/4')
aw4 = ArtinWedderburn(s4)
pass_fail(aw4)
