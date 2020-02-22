
#ifndef _FAKE_MAPLE_EXTERNAL_C_H_
#define _FAKE_MAPLE_EXTERNAL_C_H_

#define M_INT long

typedef M_INT ***ALGEB;
typedef M_INT M_BOOL;
typedef unsigned M_INT M_UINT;
typedef int* MCallBackVector;

typedef struct MapleKernelVector {
  int i;
} *MKernelVector;


#define EXT_DECL
#define M_DECL
#define M_CDECL

#ifdef __cplusplus
extern "C" {
#endif

EXT_DECL ALGEB M_DECL EvalMapleStatement( MKernelVector kv, const char *statement );
EXT_DECL ALGEB M_DECL ToMapleName( MKernelVector kv, const char *n, M_BOOL is_global );
EXT_DECL ALGEB M_CDECL EvalMapleProc( MKernelVector kv, ALGEB fn, int nargs, ... );
EXT_DECL M_BOOL M_DECL MapleToM_BOOL( MKernelVector kv, ALGEB s );
EXT_DECL MKernelVector M_DECL StartMaple( int argc, char *argv[],
                MCallBackVector cb, void *user_data, void *info, char *errstr );
EXT_DECL void M_DECL StopMaple( MKernelVector kv );
EXT_DECL M_BOOL M_DECL RestartMaple( MKernelVector kv, char *errstr );
EXT_DECL M_BOOL M_DECL IsMapleUnnamedZero( MKernelVector kv, ALGEB s );
EXT_DECL char* M_DECL MapleToString( MKernelVector kv, ALGEB s );
EXT_DECL ALGEB M_DECL ToMapleInteger( MKernelVector kv, M_INT i );
EXT_DECL ALGEB M_CDECL ToMapleFunction( MKernelVector kv, ALGEB fn, int nargs, ... );

#ifdef __cplusplus
}
#endif

#endif /* _FAKE_MAPLE_EXTERNAL_C_H_ */
