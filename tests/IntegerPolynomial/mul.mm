
fd := fopen("degree.dat", READ):
d := readdata(fd, integer)[1]:
fclose(fd):

#print(d);

#d := 8192:

#print(d):

#while d < 4194305 do 

	a := readdata("a.dat", d):
	b := readdata("b.dat", d):
	Ax := add(op(1,a[i])*x^(i-1),i=1..d):
	Bx := add(op(1,b[i])*x^(i-1),i=1..d):

	t1 := time():
	Cx := expand(Ax*Bx):
	t2 := time():

	print(t2-t1):

	#d := d*2:

#end do:

#print(Ax):
#print(Bx):
#print(Cx):
