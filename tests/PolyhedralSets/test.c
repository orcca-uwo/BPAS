#include "PolyhedralSets/FME_Support_projection.h"
#include <stdlib.h>

int main(int argc, char ** argv)
{
	printf("\n\n*** Minimal Polyhedron Projection *** \n\n");


	
	project(argv[1], atoi(argv[2]), atoi(argv[3]));



	return 0;
}
