! Common use statements and preprocessor macros
! should be included in the use statements

  USE cp_log_handling
  USE cp_error_handling

! The following macros are here to facilitate the use of error handling
! proposed in cp_error_handling.
! they assume at least
! 'use cp_error_handling, only: cp_assert, cp_a_l, cp_error_type'
! and 'use cp_log_handling, only: cp_to_string'
! They ere useful because they give a reference to the file and line
! number in the error message.


! this macro expands to a string that contains the filename.
! if the path is long the filename could make that some lines of code
! become too long and overlow (fortran compilers have a maximum line length)
! in this case substitute __FILE__ with "file" down here.
! obviously then the error messages will not give the filename.
! (otherwise make the file reference in the makefile relative vs. absolute)
#ifdef __SHORT_FILE__
#define CPSourceFileRef __SHORT_FILE__//' line '//TRIM(ADJUSTL(cp_to_string(__LINE__)))
#else
#define CPSourceFileRef __FILE__//' line '//TRIM(ADJUSTL(cp_to_string(__LINE__)))
#endif

! if the following macro is defined the longest form of macro
! expansions is used (and the error messages are more meaningful)

! inlines the test but does not write the file name, but that should
! be easily recovered from the routineP variable (that should contain
! the module name).
!
! We are trying to use a small amount of characters
! the test is not inlined (you have always a function call)

#define CPPrecondition(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error,failure)
#define CPPostcondition(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error,failure)
#define CPInvariant(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error,failure)
#define CPAssert(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error,failure)
#define CPErrorMessage(level,routineP,msg,error) \
CALL cp_error_message(level,routineP,msg,error)
#define CPPreconditionNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error)
#define CPPostconditionNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error)
#define CPInvariantNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error)
#define CPAssertNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error)
#define CPPreconditionNoErr(cond, level, routineN) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineN,__LINE__)
#define CPPostconditionNoErr(cond, level, routineN) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineN ,__LINE__)
#define CPInvariantNoErr(cond,level,routineP) \
IF(.NOT.(cond))CALL cp_a_l(0==1, level, routineP,__LINE__)
#define CPAssertNoErr(cond,level,routineP) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__)

#define CPAAssert(condition) \
CALL cp_simple_assert(condition, routineP, __LINE__)