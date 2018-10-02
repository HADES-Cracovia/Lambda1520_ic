
import subprocess
#print("run job for all channels")
import time
import os

#with open(fname) as f:
 #   content = f.readlines()
#content = [x.strip() for x in content]

#for k in content:   #take every name from vector content
#while():
#   time.sleep(30)
bashCommand = "squeue -u knowakow|wc -l"
output = subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)
print(bashCommand)
os.system("squeue -u knowakow")
a=os.system(bashCommand)
print(a)
print(output)
#print("coped following files:")
#for k in content:
 #   print(k)
