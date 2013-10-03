#!/usr/bin/python
import sys
import subprocess
import random
import time

if len(sys.argv) != 7:
	print "Usage: " + sys.argv[0] + " <csogenerator> r m n p a"
	sys.exit(1)

csogenerator = sys.argv[1]
r = int(sys.argv[2])
list_m = map(int, sys.argv[3].split(","))
list_n = map(int, sys.argv[4].split(","))
list_p = map(float, sys.argv[5].split(","))
list_a = map(int, sys.argv[6].split(","))

old_seed = 0
for m in list_m:
	for n in list_n:
		for p in list_p:
			for a in list_a:
				if a <= m:
					for i in range (r):
					    #seed = int(time.time())
					    #r1 = random.Random(seed)
					    #actual_random = r1.random()
					    #actual_random =int(actual_random)
					    start = time.clock()
					    seed = random.randint(0,1000000)
					    if seed == old_seed:
						seed = random.randint(0,1000000)
						print "Found same seed. New seed is:", seed
					    print m, n, p, a, i, seed
					    command = csogenerator + " -s " + str(seed) + " -m "+ str(m) + " -n "+ str(n) +" -p "+ str(p) +" -a "+ str(a)
					    f = open("-s_" + str(seed) + "-m_" + str(m) + "-n_" + str(n) + "-p_" + str(p) + "-a_" + str(a)+ "-i_" + str(i) + ".xml", "w")
					    subprocess.call(command, shell=True, stdout=f)
					    old_seed = seed
					    f.close()