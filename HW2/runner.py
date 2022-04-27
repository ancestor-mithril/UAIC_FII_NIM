import subprocess
from multiprocessing import Pool
import multiprocessing
multiprocessing.cpu_count()

dimensions = 10


# swarms = 2
# swarms = 3

populationSizes = [50, 100, 300]
resetThresholds = [ 20, 100, 200, 1000]
inertia = [0.1, 0.3, 0.5]
# TODO: Add more
cognition = [0.5,  1.5, 2.0]
# TODO: Add more
social = [ 1.5,  2.0, 3.0]
# TODO: Add more
swarmAttraction = [ 0.1, 0.01, 0.001]
# TODO: Add more
chaosCoef = [ 0.0,  0.001, 0.01]
# TODO: Add more
topology = [ "Star", "Ring"]
selection = [ "true", "false"]
jitter = [ "true", "false"]



def create_experiments():
 file = open("experiment.txt", "w")
 for a in populationSizes:
  for b in resetThresholds:
   for c in inertia:
    for d in cognition:
     for e in social:
      for f in swarmAttraction:
       for g in chaosCoef:
        for h in topology:
        #  for i in selection:
        #   for j in jitter:
           file.write(f"{a} {b} {c} {d} {e} {f} {g} {h}\n")
  file.close()

lines = [x.rstrip('\n') for x in open("experiment.txt", "r").readlines()]

lines1 = lines[:len(lines)//2]
lines2 = lines[len(lines)//2:]




