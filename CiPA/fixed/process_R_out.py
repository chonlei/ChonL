#!/usr/bin/env python2

for i in xrange(10):
    f = open("control_%d/out.txt"%i, "r")
    g = open("control_%d/Chaste_out.txt"%i, "w")

    next(f)
    for line in f:
        if line.strip():
            g.write("\t".join(line.split()[2:-9]) + "\n")
        else:
            raise Exception("Not working...")

    f.close()
    g.close()
