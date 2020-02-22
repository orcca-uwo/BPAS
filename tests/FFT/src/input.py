import random
import sys
size=int(sys.argv[1])
file = open("input.h","w")
file.write("#define K %i\n"%size)
s = "#define tab {"
prime = 4179340454199820289
s+= "%i"%random.randint(0,prime-1)
rand = [random.randint(0,prime-1) for _ in range(0,size-1)]
file.write(s)
for k in range(0,size-1):
	file.write(",%i"%rand[k])
file.write("}\n")
file.close()
