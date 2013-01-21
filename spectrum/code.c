/* C file produced by GAC */
#include "src/compiled.h"

/* global variables used in handlers */
static GVar G_Error;
static Obj  GF_Error;
static GVar G_f;

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
 Obj t_1 = 0;
 Obj t_2 = 0;
 Obj t_3 = 0;
 Bag oldFrame;
 OLD_BRK_CURR_STAT
 
 /* allocate new stack frame */
 SWITCH_TO_NEW_FRAME(self,0,0,oldFrame);
 REM_BRK_CURR_STAT();
 SET_BRK_CURR_STAT(0);
 
 /* if Length( x[1] ) = 2 then */
 C_ELM_LIST_FPL( t_3, a_x, INTOBJ_INT(1) )
 C_LEN_LIST_FPL( t_2, t_3 )
 t_1 = (Obj)(UInt)(((Int)t_2) == ((Int)INTOBJ_INT(2)));
 if ( t_1 ) {
  
  /* return [ x[10][1], x[10][2] ]; */
  t_1 = NEW_PLIST( T_PLIST, 2 );
  SET_LEN_PLIST( t_1, 2 );
  C_ELM_LIST_FPL( t_3, a_x, INTOBJ_INT(10) )
  C_ELM_LIST_FPL( t_2, t_3, INTOBJ_INT(1) )
  SET_ELM_PLIST( t_1, 1, t_2 );
  CHANGED_BAG( t_1 );
  C_ELM_LIST_FPL( t_3, a_x, INTOBJ_INT(10) )
  C_ELM_LIST_FPL( t_2, t_3, INTOBJ_INT(2) )
  SET_ELM_PLIST( t_1, 2, t_2 );
  CHANGED_BAG( t_1 );
  RES_BRK_CURR_STAT();
  SWITCH_TO_OLD_FRAME(oldFrame);
  return t_1;
  
 }
 /* fi */
 
 /* Error( "not ok" ); */
 t_1 = GF_Error;
 C_NEW_STRING( t_2, 6, "not ok" )
 CALL_1ARGS( t_1, t_2 );
 
 /* return; */
 RES_BRK_CURR_STAT();
 SWITCH_TO_OLD_FRAME(oldFrame);
 return 0;
 
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
 
 /* f := function ( x )
      if Length( x[1] ) = 2  then
          return [ x[10][1], x[10][2] ];
      fi;
      Error( "not ok" );
      return;
  end; */
 t_1 = NewFunction( NameFunc[2], NargFunc[2], NamsFunc[2], HdlrFunc2 );
 ENVI_FUNC( t_1 ) = CurrLVars;
 t_2 = NewBag( T_BODY, 0 );
 BODY_FUNC(t_1) = t_2;
 CHANGED_BAG( CurrLVars );
 AssGVar( G_f, t_1 );
 
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
 InitFopyGVar( "Error", &GF_Error );
 
 /* information for the functions */
 InitGlobalBag( &DefaultName, "code.g:DefaultName(-80113791)" );
 InitHandlerFunc( HdlrFunc1, "code.g:HdlrFunc1(-80113791)" );
 InitGlobalBag( &(NameFunc[1]), "code.g:NameFunc[1](-80113791)" );
 InitHandlerFunc( HdlrFunc2, "code.g:HdlrFunc2(-80113791)" );
 InitGlobalBag( &(NameFunc[2]), "code.g:NameFunc[2](-80113791)" );
 
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
 G_Error = GVarName( "Error" );
 G_f = GVarName( "f" );
 
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
 G_Error = GVarName( "Error" );
 G_f = GVarName( "f" );
 
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
 /* name        = */ "code.g",
 /* revision_c  = */ 0,
 /* revision_h  = */ 0,
 /* version     = */ 0,
 /* crc         = */ -80113791,
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
