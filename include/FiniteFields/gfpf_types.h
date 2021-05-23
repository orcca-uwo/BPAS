#ifndef GFPF_TYPES_H_
#define GFPF_TYPES_H_

#include <stdio.h>
#ifdef SERIAL
#define cilk_spawn
#define cilk_sync
#define cilk_for for
#endif
/**************************************/
typedef unsigned short usfixn16;
typedef unsigned int usfixn32;
typedef unsigned long long int usfixn64;

#define BASE_U32 4294967296
#define U32_BASE 4294967296
/**************************************/

#define BASE_U64 18446744073709551616L
#define U64_BASE 18446744073709551616L

/**************************************/
//maximum value that can be set to a u64
#define U64_MAX 18446744073709551615ULL
#define U32_MAX 4294967295U

#define U64_MASK 18446744073709551615UL
#define U32_MASK 4294967295

/**************************************/
#define MAX_CMD_LEN 4096
#define MAX_PATH_LEN 4096

/**************************************/

typedef struct
{
  int k;
  usfixn64 radix;
} srgfn_prime;

/**************************************/
/**************************************/
#ifndef LINES
#define LINES
static char const *quadline = "------------------\n";
static char const *shortline = "----------------------------------\n";
static char const *longline =
    "---------------------------------------------------------------\n";
#endif

/**************************************/

static void
print_verification_msg (const char * msg, int status)
{
  if (msg == NULL)
    {
      printf ("ERROR: in print_verification_msg; msg is NULL!\n");
      exit (EXIT_FAILURE);
    }
  printf ("[verifying %s]...", msg);
  if (status == EXIT_FAILURE)
    {
      printf ("FAILED!\n");
      exit (EXIT_FAILURE);
    }
  else
    printf ("VERIFIED!\n");
  printf ("%s", quadline);
}

/**************************************/

static void
printVector (usfixn64 * vector, int vectorSize, int coefficientSize,
	     const char* vectorName)
{
  for (int i = 0; i < vectorSize; i++)
    {
      for (int j = 0; j < coefficientSize; j++)
//	cout
//	    << vectorName << "[" << i << "," << j << "]="
//	    << vector[i * coefficientSize + j] << endl;
	printf ("%s[%d,%d] = %llu \n", vectorName, i, j,
		vector[i * coefficientSize + j]);
      printf ("\n");
//      cout << endl;
    }
}

/**************************************/

static void
compare_vectors (usfixn64 * x, usfixn64 * y, int vectorSize,
		 int coefficientSize, const char* x_name, const char* y_name)
{
  int idx = 0;
  for (int i = 0; i < vectorSize; i++)
    {
      for (int j = 0; j < coefficientSize; j++)
	{
//	cout
//	    << vectorName << "[" << i << "," << j << "]="
//	    << vector[i * coefficientSize + j] << endl;
	  printf ("[%s,%s][%d,%d] = [%llu\t,\t%llu]", x_name, y_name, i, j,
		  x[idx], y[idx]);
	  if (x[idx] != y[idx])
	    printf ("... MISMATCH");
	  printf ("\n");
	  idx++;
//      cout << endl;
	}
    }
}

/**************************************/

static void
printVectorToFile (usfixn64 * vector, int vectorSize, int coefficientSize,
		   const char* fileName, int mode)
{
  FILE * writeFile;

  //write only
  if (mode == 0)
    writeFile = fopen (fileName, "w");
  //append to the end of file
  else if (mode == 1)
    {
      writeFile = fopen (fileName, "a");
    }

  usfixn32 idx = 0;
  for (int i = 0; i < vectorSize; i++)
    {
      for (int j = 0; j < coefficientSize; j++)
	{
	  fprintf (writeFile, "%llu\n", vector[idx++]);
	}
      fprintf (writeFile, "\n");
    }

  fclose (writeFile);
}

/**************************************/

static usfixn64*
readVectorFromFile (usfixn64 *vector, int vectorSize, int coefficientSize,
		    const char *fileName)
{
  FILE * readFile = fopen (fileName, "r");
  for (int i = 0; i < vectorSize; i++)
    for (int j = 0; j < coefficientSize; j++)
      {
	//ToDo: double check the following line.
	if (fscanf (readFile, "%llu", &vector[i * coefficientSize + j]) == EOF)
	  break;
      }

  fclose (readFile);
  return vector;
}

/**************************************/

#endif // end of GFPF_TYPES_H_
