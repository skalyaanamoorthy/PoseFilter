
#global res_min, res_max, olig_num, working_dir, toAlph, Testing
#Testing = 0
# Will be used as a switch- when it is on then we have lots of output to test, off we do not have any
global olig_num, res_min, res_max, toAlph
olig_num = 0
res_min = -1
res_max = -1
toAlph = []

def SetOligNum(num):
    global olig_num
    olig_num = num

def SetRes(rmin, rmax):
    global res_min, res_max
    res_min = rmin
    res_max = res_max