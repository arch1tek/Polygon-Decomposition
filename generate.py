# -----------------------------------------------------------------------------
#   Randomly generate points in a [0,1) x [0,1) square and write to generate_test1.txt
# -----------------------------------------------------------------------------

import os
import sys
import random

# prompt user to enter the number of vertices
n = int(sys.argv[1])

f = open("generate_test1.txt", "w+")
f.write("%i\n" % n)
for i in range(n):
    x = random.random()
    y = random.random()
    f.write("%f %f\n" % (x, y))
f.close()

os.system("python3 main_improved.py");
