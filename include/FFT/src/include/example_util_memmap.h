// -*- C++ -*-

#ifndef CILK_EXAMPLE_UTIL_MEMMAP_H_INCLUDED
#define CILK_EXAMPLE_UTIL_MEMMAP_H_INCLUDED

/* 
 * These simple utility functions hide Linux/Windows differences so that we can provide
 * examples with source code common to the different platforms.  Feel free
 * to use this unsupported file at your own risk.
 */

#ifdef _WIN32
#include <Windows.h>
#else
#include <stdlib.h>
#endif

#ifdef WIN32
extern "C++"
typedef struct { // Structure for mapped file
    void * pInFile;
    HANDLE hInMap;
    HANDLE hIn;
} MAPPED_FILE_HANDLE;

extern "C++"
inline
void * example_map_file(LPCSTR filename, unsigned int * pFsLow, int * error, MAPPED_FILE_HANDLE * pmFH)
{
    *error = 0;
    HANDLE hIn = CreateFile (filename, GENERIC_READ, 0, NULL,
                             OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hIn == INVALID_HANDLE_VALUE) {
        *error = 2;
        return NULL;
    }

    // Create a file mapping object on the input file. Use the file size. 
    HANDLE hInMap = CreateFileMapping (hIn, NULL, PAGE_READONLY, 0, 0, NULL);
    if (hInMap == INVALID_HANDLE_VALUE) {
        CloseHandle(hIn);
        *error = 3;
        return NULL;
    }

    // Map the input file 
    char * pInFile = (char *)MapViewOfFile (hInMap, FILE_MAP_READ, 0, 0, 0);
    if (pInFile == NULL) {
        CloseHandle (hInMap);
        CloseHandle (hIn);
        *error = 4;
        return NULL;
    }

    // Get the input file size. As the mapping succeeded, the file size is < 2 GB. 
    *pFsLow = GetFileSize (hIn, NULL);
    if (*pFsLow == 0xFFFFFFFF) {
        UnmapViewOfFile (pInFile);
        CloseHandle (hInMap);
        CloseHandle (hIn);
        *error = 5;
        return NULL;
    }

    pmFH->pInFile = pInFile;
    pmFH->hInMap  = hInMap;
    pmFH->hIn     = hIn;

    return pInFile;
}

extern "C++"
inline
void UnMapFile (MAPPED_FILE_HANDLE * pmFH)
{
     UnmapViewOfFile (pmFH->pInFile);
     CloseHandle (pmFH->hInMap);
     CloseHandle (pmFH->hIn);
}
#else  // Linux Memory mapping

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

extern "C++"
typedef struct { // Structure for mapped file
    void * pInFile;
    int hIn;        // file descriptor
    int fLen;
} MAPPED_FILE_HANDLE;

extern "C++"
inline
void * example_map_file(char* filename, unsigned int * pFsLow, int * error, MAPPED_FILE_HANDLE * pmFH)
{
    *error = 0;
    int hIn = open (filename, O_RDONLY);

    if (-1 == hIn) {
        *error = 2;
        return NULL;
    }

    struct stat fstat;
    int status = stat(filename, &fstat);
    if (-1 == status) {
        *error = 5;
        close (hIn);
        return NULL;
    }
    *pFsLow = (int)fstat.st_size;

    // Map the input file 
    char * pInFile = (char *)mmap(0, *pFsLow, PROT_READ, MAP_PRIVATE, hIn, 0);
    if (NULL == pInFile) {
        close (hIn);
        *error = 4;
        return NULL;
    }


    pmFH->pInFile = pInFile;
    pmFH->hIn     = hIn;
    pmFH->fLen    = *pFsLow;

    return pInFile;
}

extern "C++"
inline
void UnMapFile (MAPPED_FILE_HANDLE * pmFH)
{
     munmap (pmFH->pInFile, (off_t)pmFH->fLen);

     close (pmFH->hIn);
}

#endif //WIN32

#endif // CILK_EXAMPLE_UTIL_MEMMAP_H_INCLUDED
