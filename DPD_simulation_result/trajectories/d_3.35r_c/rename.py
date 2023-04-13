import os
import glob


path = glob.glob("*\\*.mat")


for p in path:
    
    (a,b) = p.split("\\");
    
    os.rename(a,"L_"+a[2:])
    