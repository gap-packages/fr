/* C file produced by GAC */
#include "src/compiled.h"
#include <nag.h>
#include <nag_stdlib.h>
#include <nag_string.h>
#include <stdio.h>
#include <nagf12.h>
#include <nagf16.h>

/* global variables used in handlers */
static GVar G_FLOAT__STRING;
static Obj  GF_FLOAT__STRING;
static GVar G_F12FAC;
static GVar G_InfoFR;
static Obj  GC_InfoFR;

/* record names used in handlers */

/* information for the functions */
static Obj  NameFunc[3];
static Obj  NamsFunc[3];
static Int  NargFunc[3];
static Obj  DefaultName;

static void PrintInfo (char *s) {
  Obj t_1, t_2, t_3;

  t_1 = GC_InfoFR;
  CHECK_BOUND( t_1, "InfoFR" )
    t_3 = CALL_2ARGS( InfoDecision, t_1, INTOBJ_INT(2) );
  if ( t_3 == True ) {
    t_2 = NEW_PLIST( T_PLIST, 2 );
    SET_LEN_PLIST( t_2, 2 );
    C_NEW_STRING( t_3, strlen(s), s )
      SET_ELM_PLIST( t_2, 1, t_3 );
    CHANGED_BAG(t_2);
    SET_ELM_PLIST( t_2, 2, a_nev );
    CHANGED_BAG(t_2);
    CALL_1ARGS( InfoDoPrint, t_2 );
  }
}

/* handler for function 2 */
static Obj  HdlrFunc2 (
 Obj  self,
 Obj  a_pts,
 Obj  a_nev,
 Obj  a_ncv )
{
  Obj t_1 = 0;
  Obj t_2 = 0;
  Obj t_3 = 0;
  Obj t_4 = 0;
  Bag oldFrame;
  OLD_BRK_CURR_STAT;
    
    /* allocate new stack frame */
  SWITCH_TO_NEW_FRAME(self,0,0,oldFrame);
  REM_BRK_CURR_STAT();
  SET_BRK_CURR_STAT(0);

  /* pts[1][2] := 2; */
  C_ELM_LIST_FPL( t_1, a_pts, INTOBJ_INT(1) );
  C_ASS_LIST_FPL_INTOBJ( t_1, INTOBJ_INT(2), INTOBJ_INT(2) );
 
 
    /* return [ FLOAT_STRING( "3.14" ) ]; */
  t_1 = NEW_PLIST( T_PLIST, 1 );
  SET_LEN_PLIST( t_1, 1 );
  t_3 = NewBag(T_FLOAT,sizeof(double));
  *((double *) ADDR_OBJ(t_3)) = 3.14159265;
  SET_ELM_PLIST( t_1, 1, t_3 );
  CHANGED_BAG( t_1 );
  RES_BRK_CURR_STAT();
  SWITCH_TO_OLD_FRAME(oldFrame);
  return t_1;
}

/* handler for function 1 */
static Obj  HdlrFunc1 (
 Obj  self )
{
 Obj t_1 = 0;
 Obj t_2 = 0;
 Bag oldFrame;
 OLD_BRK_CURR_STAT
 
 /* allocate new stack frame */
 SWITCH_TO_NEW_FRAME(self,0,0,oldFrame);
 REM_BRK_CURR_STAT();
 SET_BRK_CURR_STAT(0);
 
 /* F12FAC := function ( pts, nev, ncv )
      pts[1][2] := 2;
      Info( InfoFR, 2, "Eigenvalues are", nev );
      return [ FLOAT_STRING( "3.14" ) ];
  end; */
 t_1 = NewFunction( NameFunc[2], NargFunc[2], NamsFunc[2], HdlrFunc2 );
 ENVI_FUNC( t_1 ) = CurrLVars;
 t_2 = NewBag( T_BODY, 0 );
 BODY_FUNC(t_1) = t_2;
 CHANGED_BAG( CurrLVars );
 AssGVar( G_F12FAC, t_1 );
 
 /* return; */
 RES_BRK_CURR_STAT();
 SWITCH_TO_OLD_FRAME(oldFrame);
 return 0;
 
 /* return; */
 RES_BRK_CURR_STAT();
 SWITCH_TO_OLD_FRAME(oldFrame);
 return 0;
}

/* 'InitKernel' sets up data structures, fopies, copies, handlers */
static Int InitKernel ( StructInitInfo * module )
{
 
 /* global variables used in handlers */
 InitFopyGVar( "FLOAT_STRING", &GF_FLOAT__STRING );
 InitCopyGVar( "InfoFR", &GC_InfoFR );
 
 /* information for the functions */
 InitGlobalBag( &DefaultName, "spectrum.g:DefaultName(117446121)" );
 InitHandlerFunc( HdlrFunc1, "spectrum.g:HdlrFunc1(117446121)" );
 InitGlobalBag( &(NameFunc[1]), "spectrum.g:NameFunc[1](117446121)" );
 InitHandlerFunc( HdlrFunc2, "spectrum.g:HdlrFunc2(117446121)" );
 InitGlobalBag( &(NameFunc[2]), "spectrum.g:NameFunc[2](117446121)" );
 
 /* return success */
 return 0;
 
}

/* 'InitLibrary' sets up gvars, rnams, functions */
static Int InitLibrary ( StructInitInfo * module )
{
 Obj func1;
 Obj body1;
 
 /* Complete Copy/Fopy registration */
 UpdateCopyFopyInfo();
 
 /* global variables used in handlers */
 G_FLOAT__STRING = GVarName( "FLOAT_STRING" );
 G_F12FAC = GVarName( "F12FAC" );
 G_InfoFR = GVarName( "InfoFR" );
 
 /* record names used in handlers */
 
 /* information for the functions */
 C_NEW_STRING( DefaultName, 14, "local function" )
 NameFunc[1] = DefaultName;
 NamsFunc[1] = 0;
 NargFunc[1] = 0;
 NameFunc[2] = DefaultName;
 NamsFunc[2] = 0;
 NargFunc[2] = 3;
 
 /* create all the functions defined in this module */
 func1 = NewFunction(NameFunc[1],NargFunc[1],NamsFunc[1],HdlrFunc1);
 ENVI_FUNC( func1 ) = CurrLVars;
 CHANGED_BAG( CurrLVars );
 body1 = NewBag( T_BODY, 0);
 BODY_FUNC( func1 ) = body1;
 CHANGED_BAG( func1 );
 CALL_0ARGS( func1 );
 
 /* return success */
 return 0;
 
}

/* 'PostRestore' restore gvars, rnams, functions */
static Int PostRestore ( StructInitInfo * module )
{
 
 /* global variables used in handlers */
 G_FLOAT__STRING = GVarName( "FLOAT_STRING" );
 G_F12FAC = GVarName( "F12FAC" );
 G_InfoFR = GVarName( "InfoFR" );
 
 /* record names used in handlers */
 
 /* information for the functions */
 NameFunc[1] = DefaultName;
 NamsFunc[1] = 0;
 NargFunc[1] = 0;
 NameFunc[2] = DefaultName;
 NamsFunc[2] = 0;
 NargFunc[2] = 3;
 
 /* return success */
 return 0;
 
}


/* <name> returns the description of this module */
static StructInitInfo module = {
 /* type        = */ 3,
 /* name        = */ "spectrum.g",
 /* revision_c  = */ 0,
 /* revision_h  = */ 0,
 /* version     = */ 0,
 /* crc         = */ 117446121,
 /* initKernel  = */ InitKernel,
 /* initLibrary = */ InitLibrary,
 /* checkInit   = */ 0,
 /* preSave     = */ 0,
 /* postSave    = */ 0,
 /* postRestore = */ PostRestore
};

StructInitInfo * Init__Dynamic ( void )
{
 return &module;
}

/* compiled code ends here */
