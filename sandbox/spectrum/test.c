/* C file produced by GAC */
#include "src/compiled.h"

#include "nag.h"

/* global variables used in handlers */
static GVar G_test;

/* record names used in handlers */

/* information for the functions */
static Obj  NameFunc[3];
static Obj  NamsFunc[3];
static Int  NargFunc[3];
static Obj  DefaultName;

/* handler for function 2 */
static Obj  HdlrFunc2 (
 Obj  self,
 Obj  a_x )
{
 Integer licomm=140, lcomm, irevcm, j, n, nconv, ncv, nev,
    nshift, nd, niter, *icomm = 0, *ix = 0, *iy = 0, exit_status;
  double *oldx = 0, *comm = 0, *eigv = 0, *eigest = 0, *mx = 0, *resid = 0,
    *v = 0, *x = 0, *y = 0, sigma;
  NagError fail;

  f12fac(n, nev, ncv, icomm, licomm, comm, lcomm, &fail);

 Obj t_1 = 0;
 Bag oldFrame;
 OLD_BRK_CURR_STAT
 
 /* allocate new stack frame */
 SWITCH_TO_NEW_FRAME(self,0,0,oldFrame);
 REM_BRK_CURR_STAT();
 SET_BRK_CURR_STAT(0);
 
 /* return x ^ 2; */
 t_1 = POW( a_x, INTOBJ_INT(2) );
 RES_BRK_CURR_STAT();
 SWITCH_TO_OLD_FRAME(oldFrame);
 return t_1;
 
 /* return; */
 RES_BRK_CURR_STAT();
 SWITCH_TO_OLD_FRAME(oldFrame);
 return 0;
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
 
 /* test := function ( x )
      return x ^ 2;
  end; */
 t_1 = NewFunction( NameFunc[2], NargFunc[2], NamsFunc[2], HdlrFunc2 );
 ENVI_FUNC( t_1 ) = CurrLVars;
 t_2 = NewBag( T_BODY, 0 );
 BODY_FUNC(t_1) = t_2;
 CHANGED_BAG( CurrLVars );
 AssGVar( G_test, t_1 );
 
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
 
 /* information for the functions */
 InitGlobalBag( &DefaultName, "test.g:DefaultName(-133654923)" );
 InitHandlerFunc( HdlrFunc1, "test.g:HdlrFunc1(-133654923)" );
 InitGlobalBag( &(NameFunc[1]), "test.g:NameFunc[1](-133654923)" );
 InitHandlerFunc( HdlrFunc2, "test.g:HdlrFunc2(-133654923)" );
 InitGlobalBag( &(NameFunc[2]), "test.g:NameFunc[2](-133654923)" );
 
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
 G_test = GVarName( "test" );
 
 /* record names used in handlers */
 
 /* information for the functions */
 C_NEW_STRING( DefaultName, 14, "local function" )
 NameFunc[1] = DefaultName;
 NamsFunc[1] = 0;
 NargFunc[1] = 0;
 NameFunc[2] = DefaultName;
 NamsFunc[2] = 0;
 NargFunc[2] = 1;
 
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
 G_test = GVarName( "test" );
 
 /* record names used in handlers */
 
 /* information for the functions */
 NameFunc[1] = DefaultName;
 NamsFunc[1] = 0;
 NargFunc[1] = 0;
 NameFunc[2] = DefaultName;
 NamsFunc[2] = 0;
 NargFunc[2] = 1;
 
 /* return success */
 return 0;
 
}


/* <name> returns the description of this module */
static StructInitInfo module = {
 /* type        = */ 3,
 /* name        = */ "test.g",
 /* revision_c  = */ 0,
 /* revision_h  = */ 0,
 /* version     = */ 0,
 /* crc         = */ -133654923,
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
