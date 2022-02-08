MODULE DDE_SOLVER_M
  !
  ! Further change history is recorded in the git commit history.
  !____________________________________________________________________________
  !
  ! Last change:  ST 15 August 2011
  !               Changes to make g95 and Salford/Checkmate (and other
  !               compilers with array bounds checking) friendlier to
  !               the solver
  !____________________________________________________________________________
  !
  ! Note:
  ! If you supply a USER_TRIM_GET function to interpolate using discarded
  ! solution information, then depending on your compiler, you may need
  ! to change the file OPEN statements below.
  !____________________________________________________________________________
  !
  ! Note:
  ! Either of the Salford compile options Debug Win32, Release Win32, and
  ! Checkmate Win 32 compile may e used. In addition, either of the g95
  ! compile options work fine.
  !    g95 -ftrace=full dde_solver_m.f90 filename.f90 -o a.exe
  !    g95 -O3 -march=pentium4 dde_solver_m.f90 filename.f90 -o a.exe 
  ! The solver also has been used with the Lahey F90 and lf95 compilers
  ! and with the Compaq, Intel, Portland, Sun/Unix, and various other
  ! F90/F95 compilers. We will of course appreciate any feedback regarding
  ! compilation or run time issues you may encounter.
  !____________________________________________________________________________
  !
  ! This Fortran DDE solver was written by L.F. Shampine and S. Thompson.
  ! You are welcome to use it, but we assume no liability for its use.
  ! We have tested the program thoroughly, but bugs are always possible
  ! in a program of this size and complexity. If you believe that you
  ! have encountered a bug or have questions regarding usage of this
  ! solver, please provide one of us with details. You can get in touch
  ! with us at
  !
  !   L.F. Shampine                     S. Thompson
  !   Mathematics Department            Mathematics & Statistics Department
  !   Southern Methodist University     Radford University
  !   Dallas, TX 75275                  Radford, VA 24142
  !   shampine@mail.smu.edu             thompson@radford.edu
  !
  ! The latest release of the source code for the DDE_SOLVER module
  ! and a collection of example programs are available at
  !
  !   http://www.radford.edu/~thompson/ffddes/index.html
  !
  ! We recommend that you download and study the example programs since
  ! they illustrate usage for a variety of problems and may be used as
  ! templates for other problems. (See also the usage description below.)
 
  ! If your usage of this software requires a license, the following
  ! BSD license is applicable. If it is not appropriate for your
  ! usage, please contact one of the authors.
  !
  ! Copyright (c) 2009
  ! L.F. Shampine, Southern Methodist University, Dallas, TX
  ! S. Thompson, Radford University, Radford, VA
  ! All rights reserved.
  !
  ! Redistribution and use in source and binary forms, with or without
  ! modification, are permitted provided that the following conditions
  ! are met:
  !
  !    Redistributions of source code must retain the above copyright
  !    notice, this list of conditions, and the following disclaimer.
  !
  !    Redistributions in binary form must reproduce the above copyright
  !    notice, this list of conditions, and the following disclaimer
  !    in the documentation and/or other materials provided with the
  !    distribution.
  !
  !    Neither the name of the organizations nor the names of its
  !    contributors may be used to endorse or promote products derived
  !    from this software without specific prior written permission.
  !
  ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  ! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  ! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  ! COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  ! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  ! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  ! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  ! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  ! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  ! POSSIBILITY OF SUCH DAMAGE.

  !____________________________________________________________________________

  !                     Usage Instructions

  ! SECTION   I. General Information
  ! SECTION  II. Constant Delays and Constant History
  ! SECTION III. Non-constant Delays and/or History
  ! SECTION  IV. Event Functions
  ! SECTION   V. Options and Callable Functions
  !____________________________________________________________________________

  ! SECTION I. General Information

  ! DDE_SOLVER is a program for solving DDEs of the form
  !           y'(t) = f(t,y(t),y(beta(t,y(t))),y'(beta(t,y(t)))
  ! They are to be integrated from t_0 to t_f subject to
  !           y(t)  = history(t) for t<=t_0.
  ! The delay times must satisfy beta(t,y(t))<=t.
  ! DDE_SOLVER was designed to make solving easy problems as simple
  ! as possible. At the same time, it is capable of solving very
  ! complicated problems. The associated web site has a manuscript
  ! that explains the design and how it is realized. In simplest use
  ! just about all you have to do is define the problem. Even with this
  ! minimal call list of required arguments, there are simplifications.
  ! One is for the very common case of constant history function. For
  ! such a problem you can simply input the vector instead of writing a
  ! subroutine as for general history functions. Options are handled in
  ! two ways. Some are optional arguments to the solver that can be
  ! specified in any order in the call list after the required arguments
  ! using a keyword. This is the way that event functions are specified.
  ! One of the optional arguments is an options structure. An options
  ! structure is formed with the auxiliary function DDE_SET. All the
  ! arguments of this function are optional. They are set in any order
  ! using keywords. This is the way that error tolerances are set. If
  ! an option is not set in one of these two ways, it is given a
  ! default value.
  ! The web site has a collection of examples that illustrate the use
  ! of the solver. Each of these examples writes information to data
  ! files. A companion M-file reads this data and plots the solution
  ! in Matlab. The primary examples are taken from Chapter 4 of the
  ! book "Solving ODEs with Matlab", by L.F. Shampine, I. Gladwell, and
  ! S. Thompson published by Cambridge University Press.
  !   1. Example 4.4.1: Constant Delays and Constant History
  !   2. Example 4.4.2: Interpolation
  !   3. Example 4.4.3: Events and Change Routine
  !   4. Example 4.4.4: Events and use of the JUMPS Option
  !   5. Example 4.4.5: Events and Change Routine
  !   6. Neutral problem with state-dependent vanishing delay
  !____________________________________________________________________________

  ! SECTION II. Constant Delays and Constant History

  ! The most common type of problem is one for which the delays are
  ! constant so that the delay times are t-delay_1,...,t-delay_NEQN,
  ! and the history is a constant vector. DDE_SOLVER is designed to
  ! be particularly easy to use for such problems, at least if all
  ! the default values are acceptable. In this case, DDE_SOLVER can
  ! be invoked with
  ! SOL = DDE_SOLVER(NVAR,DDES,DELAYS,HISTORY,TSPAN)
  ! These arguments are required in every call to the solver. The only
  ! thing special here is the form of DELAYS and HISTORY. The arguments
  ! of this minimal list have the following meaning:
  ! NVAR is a vector of length 2 where
  !    NVAR(1) = NEQN is the number of DDEs
  !    NVAR(2) = NLAGS is the number of delays
  ! DDES is a subroutine of the form DDES(T,Y,Z,DY) that evaluates
  !    the differential equations for y'(t). The input arguments are
  !    T - the current time t
  !    Y - a vector of length NEQN containing an approximation
  !        to the current solution y(t)
  !    Z - an array with dimensions NEQN by NLAGS. Z(I,J) is an
  !        approximation to y_i(t-delay_j).
  !    The output argument is
  !    DY - a vector of length NEQN that contains y'(t).
  ! DELAYS is a vector of length NLAGS that contains the delays
  ! delay_1,...,delay_NLAGS.
  ! HISTORY is a vector of length NEQN containing the constant
  ! history HISTORY(i) = y_i(t) for t<=T0, where T0 is the
  ! initial time.
  ! TSPAN is a vector of length 2 where
  !    TSPAN(1) = T0, the initial time
  !    TSPAN(2) = TF, the final time
  ! DDE_SOLVER returns a solution structure SOL. In this usual mode of
  ! output, the solver selects the step size and returns the solution
  ! at all these points. Accordingly, there are fields
  !   SOL%NPTS         -- NPTS, number of output points
  !   SOL%T(NPTS)      -- values of independent variable, T
  !   SOL%Y(NPTS,NEQN) -- values of dependent variable, Y,
  !                       corresponding to values of SOL%T
  ! The example program ex4p4p1.f90 solves a problem with constant
  ! delays and constant history. It writes the solution to a data file.
  ! ex4p4p1.m reads the data file and plots the solution in Matlab.
  ! There are other possibilities for output. To specify in advance
  ! where you want the solution, just supply these times as the vector
  ! TSPAN. Just as with the usual mode of output, the first entry of
  ! TSPAN is T0 and the last is TF. All the entries of TSPAN must be in
  ! order. Options described below give you even more control. The
  ! INTERPOLATION option lets you use SOL to evaluate y(t) and y'(t)
  ! for any t in [t_0,t_f]. In effect, you have a function for y(t)
  ! and y'(t). To customize the output, you can write your
  ! own output function that the solver will call after each step.
  !____________________________________________________________________________

  ! SECTION III. Non-constant Delays and/or History

  ! If the delays and/or the history are not constant, you define
  ! them by means of subroutines.
  ! If the history is not constant, the HISTORY argument of the solver
  ! is a subroutine of the form HISTORY(T,Y). For given value of time
  ! T<=T0, the subroutine is to evaluate history(T) and return it as
  ! the vector Y of length NEQN.
  ! If the delays are not constant, the BETA argument of the solver
  ! is a subroutine of the form BETA(T,Y,BVAL). For input scalar T and
  ! vector Y of length NEQN, the subroutine is to evaluate the delays
  ! beta_i(T,Y) and return them as the vector BVAL of length NLAGS.
  !____________________________________________________________________________

  ! SECTION IV. Event Functions

  ! DDE_SOLVER provides for event functions. For i = 1,...,NEF, each
  ! event function has the form
  !             g_i(t,y(t),y(beta(t,y)),y'(t))
  ! An event is said to occur where one of these functions vanishes.
  ! The task is to locate events and report the solution there. This
  ! is done by means of the fields
  !   SOL%NE           -- NE, number of events
  !   SOL%TE(NE)       -- locations of events
  !   SOL%YE(NE,NEQN)  -- values of solution at events
  !   SOL%IE(NE)       -- identifies which event occurred
  ! You can define an event as terminal and have the integration stop
  ! then. Also, you can tell the solver that you are interested in
  ! roots of a particular event function only when it decreases through
  ! zero or increases through zero or both. These qualifiers are
  ! defined by means of an options structure and so are discussed in
  ! SECTION V.
  ! Sometimes you want to change the nature of the problem after an
  ! event occurs. This is done by means of a change function. The
  ! solver will call this function after each event. You could then,
  ! e.g., change problem parameters for the rest of the integration.
  ! To define the event functions, you need only inform DDE_SOLVER of
  ! how many event functions are present by changing NVARS to a vector
  ! of length 3 with NVARS(3) = NEF, the number of event functions, and
  ! supply a subroutine of the form EF(T,Y,DY,Z,G). The input arguments
  ! of EF are
  !    T  -  the current time t
  !    Y  -  a vector of length NEQN containing the solution y(t)
  !    DY -  a vector of length NEQN containing the derivative y'(t)
  !    Z -   an array with dimensions NEQN by NLAGS. The entry
  !          Z(I,J) contains y_i(beta_j(t,y)).
  ! EF is to evaluate the event functions and return these residuals
  ! in the array G of length NEF.
  ! If you wish to make changes at events, you need only supply a
  ! subroutine CHNG of the form
  ! CHNG(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT,DIRECTION,ISTERMINAL, &
  !      QUIT)
  ! The input arguments of CHNG are
  !    NEVENT-  the component of G that vanishes
  !    TEVENT-  the time t at which the event occurred
  !    YEVENT-  an array of length NEQN that contains
  !             the solution y(t)
  !    DYEVENT- an array of length NEQN that contains
  !             the derivative y'(t)
  !    QUIT-    a logical variable that is .FALSE. on input
  ! The output arguments are
  !    HINIT-   If you wish to supply a trial step size
  !             for the solver to attempt when the integration
  !             is re-started, you can define HINIT to be
  !             this step size.
  !    DIRECTION and ISTERMINAL - arrays of length NEF
  !             that you can supply to instruct DDE_SOLVER
  !             as to which events are reported and
  !             whether events are terminal. Each is
  !             described below in the options section.
  !    QUIT-    To terminate the integration, change the logical
  !             variable QUIT to .TRUE.
  ! The solver is informed about the presence of event functions by
  ! means of optional arguments. The name of the event function is
  ! supplied with the keyword EVENT_FCN. If there is a change function,
  ! its name is supplied with the keyword CHANGE_FCN. Optional
  ! arguments can appear anywhere in the call list after the required
  ! arguments. For example, if there are both an event function and
  ! a change function, the solver might be called with
  ! SOL = DDE_SOLVER(NVAR,DDES,BETA,HISTORY,TSPAN, &
  !                  EVENT_FCN=EF,CHANGE_FCN=CHNG)
  ! The program ex4p4p4.f90 is an example with an event function. The
  ! programs ex4p4p3.f90 and ex4p4p5.f90 are examples that have both
  ! an event function and a change function.
  !____________________________________________________________________________

  ! SECTION V. Options and Callable Functions

  ! DDE_SOLVER uses default values for most parameters. You can
  ! supply an options structure to override default values. It is an
  ! optional argument of DDE_SOLVER that is supplied with the keyword
  ! OPTIONS. An options structure is formed with DDE_SET. All of its
  ! arguments are optional and set with keywords. Following are the
  ! possibilities. The manuscript and example programs of the web site
  ! discuss them more fully.
  ! RE                    - A relative error tolerance applied to all
  !                         components. Default: RE = 1D-3.
  ! AE                    - An absolute error tolerance applied to all
  !                         components. Default: AE = 1D-6.
  ! AE_VECTOR             - Vector of absolute error tolerances applied
  !                         to corresponding components. If both AE and
  !                         AE_VECTOR are set accidentally, AE_VECTOR
  !                         is used. Default: All components = 1D-6.
  ! INTERPOLATION         - Set .TRUE. to save in SOL information needed
  !                         to evaluate y(t) and y'(t) with DDE_VAL.
  !                         Default: .FALSE.
  ! TERMINAL              - Logical vector of length NEF. Set an entry
  !                         .TRUE. if the integration is to terminate
  !                         when the corresponding event occurs.
  !                         Default: No event is terminal.
  ! DIRECTION             - Vector of length NEF. Set an entry to +1
  !                         if events are restricted to the function
  !                         increasing through 0; -1 if restricted to
  !                         decreasing through 0; and 0 to ignore the
  !                         direction. Default: Ignore direction
  ! JUMPS                 - Vector containing known discontinuities in
  !                         the history function and known discontinuities
  !                         in y'(t)at values of t>T0.
  !                         Default: JUMP in y'(t) at T0.
  ! TRACKING_LEVEL        - Integer indicating the maximum level for
  !                         tracking discontinuities. Default: infinity
  !                         for neutral problems, 7 otherwise.
  ! TRACK_DISCONTINUITIES - Logical value indicating whether to track
  !                         discontinuities. Default: .TRUE.
  ! HINIT                 - Trial initial step size. Default:
  !                         DDE_SOLVER determines the initial step size.
  ! HMAX                  - Maximum step size. Default: 0.1D0*(TF-T0).
  ! NEUTRAL               - Logical variable indicating whether the
  !                         problem is of neutral type, i.e., whether
  !                         the differential equations depend on
  !                         y'(beta(t,y(t))). Default: .FALSE.
  !                         The program neutvan.f90 illustrates solving
  !                         a neutral problem.
  ! THIT_EXACTLY          - Vector of times the solver must hit exactly.
  ! MAX_EVENTS            - Maximum allowable number of event occurrences.
  ! MAX_STEPS             - Maximum allowable number of integration steps.
  !
  ! The following two options are intended for large problems in which
  ! the size of the history queue may grow to be extremely large if the
  ! entire queue is saved. In the event they are used, it is further
  ! possible to provide an optional subroutine to perform interpolation
  ! (similar to the INTERPOLATION option). Refer to the accompanying
  ! programs by Steven White for examples of usage of these options.
  !
  ! MAX_DELAY             - Bound for maximum delay for use in trimming
  !                         the solution queue
  !                         Note:
  !                         Automatic discontinuity tracking is not
  !                         performed if MAX_DELAY>0 since necessary
  !                         solution history queue information will be
  !                         discarded.
  ! TRIM_FREQUENCY        - Frequency of trimming solution queue if
  !                         MAX_DELAY>0 is specified. Not used unless
  !                         MAX_DELAY>0.0D0.

  ! The following user callable functions are available.

  ! Set non-default integration options.
  !    OPTIONS = DDE_SET(...)
  ! Invoke the solver
  !    SOL = DDE_SOLVER(...)
  ! Print integration statistics:
  !    CALL PRINT_STATS(SOL)
  ! Approximate y(t) and y'(t) for an array of points T in [t_0,t_f].
  !   YINT = DDE_VAL(T,SOL)
  !   To evaluate only selected components, specify an integer vector
  ! COMPONENTS indicating which you want.
  !   YINT = DDE_VAL(T,SOL,COMPONENTS)
  !   If you do not want y'(t), call as
  !   YINT = DDE_VAL(T,SOL,DERIVATIVES=.FALSE.)
  !   To use DDE_VAL for the solution structure SOL, the solution must
  !   to be computed with INTERPOLATION set .TRUE. The program
  !   ex4p4p2.f90 illustrates the use of interpolation.
  ! To release all allocated arrays solution and option structures:
  !   CALL RELEASE_ARRAYS(SOL,OPTS)
  !   CALL RELEASE_INT(YINT)
  !____________________________________________________________________________

  ! Beginning of remaining DDE_SOLVER private variables and arrays.
  ! This private information is used throughout DDE_SOLVER.

  ! .. Generic Interface Blocks ..
  INTERFACE DDE_SOLVER
     !MODULE PROCEDURE DKL_1, DKL_2, DKL_3, DKL_4, ODEAVG
     MODULE PROCEDURE DKL_1, DKL_2, DKL_3, DKL_4
  END INTERFACE
  ! ..
  ! .. Derived Type Declarations ..
  TYPE, PUBLIC :: DDE_SOL
     INTEGER :: NPTS, FLAG, NE
     DOUBLE PRECISION, DIMENSION (:), POINTER :: T
     DOUBLE PRECISION, DIMENSION (:,:), POINTER :: Y
     DOUBLE PRECISION, DIMENSION (:), POINTER :: TE
     DOUBLE PRECISION, DIMENSION (:,:), POINTER :: YE
     DOUBLE PRECISION, DIMENSION (:,:), POINTER :: QUEUE
     DOUBLE PRECISION, DIMENSION (:), POINTER :: YOFT
     DOUBLE PRECISION, DIMENSION (:), POINTER :: TQUEUE
     INTEGER, DIMENSION (:), POINTER :: STATS
     INTEGER, DIMENSION (:), POINTER :: IE
     INTEGER, DIMENSION (:), POINTER :: IPOINT
     LOGICAL SHIFT
     DOUBLE PRECISION TSHIFT
  END TYPE DDE_SOL

  TYPE, PUBLIC :: DDE_OPTS
     LOGICAL, DIMENSION (:), POINTER :: ISTERMINAL
     INTEGER, DIMENSION (:), POINTER :: DIRECTION
     DOUBLE PRECISION :: HINIT, HMAX
     DOUBLE PRECISION, DIMENSION (:), POINTER :: ABSERR, RELERR, JUMPS
     LOGICAL :: NEUTRAL, TRACK_DISCONTINUITIES, INTERPOLATION
     INTEGER :: TRACKING_LEVEL, MAX_EVENTS, MAX_STEPS, MOVING_AVERAGE
     DOUBLE PRECISION, DIMENSION (:), POINTER :: THIT_EXACTLY
  END TYPE DDE_OPTS

  TYPE, PUBLIC :: DDE_INT
     DOUBLE PRECISION, DIMENSION (:), POINTER :: TVALS
     DOUBLE PRECISION, DIMENSION (:,:), POINTER :: YT, DT
     INTEGER, DIMENSION (:), POINTER :: COMPONENTS
  END TYPE DDE_INT

  TYPE, PUBLIC :: DDE_SOLVER_VERSION_TYPE
    INTEGER MAJOR, MINOR, POINT
    LOGICAL RELEASED
  END TYPE DDE_SOLVER_VERSION_TYPE

  ! ..
  ! .. Local Structures ..
  TYPE (DDE_SOL), POINTER, PRIVATE :: PASS_SOL
  ! ..
  ! .. Local Scalars ..
  DOUBLE PRECISION, PRIVATE :: EFAC, FACDN, FACUP, GRTOL, KCFAC, &
       MYTINIT, MY_MAX_DELAY, TQLAST
  DOUBLE PRECISION, PRIVATE :: REALMAX = HUGE(1D0)
  DOUBLE PRECISION, PRIVATE :: REMIN, RER, TOLFAC, U10, U13, U26, U65
  DOUBLE PRECISION, PRIVATE :: UROUND = EPSILON(1D0)
  DOUBLE PRECISION, PRIVATE, PARAMETER :: HUNDRED=1.0D2, PT05=0.05D0, &
       TEN=1.0D1
  INTEGER, PRIVATE :: ARRAY_STORAGE, IFIRST, IORDER, IPAIR, IRTMAX, &
       JN, JP, LENTREE, LIW, LQUEUE, LTQUEUE, MAXCOR, MAXTRY, MYMAX_EVENTS, &
       MYMAX_STEPS, MYMETHOD, MYN, MYNGUSER, MYNJUMPS, MYNTHIT, MYQCOLS, &
       NFAILC, NFAILS, NFEVAL, NGEMAX, NGEVAL, NLAGS, NSIG, NSPAN, NSTEPS, &
       ROOT_FUNCTIONS, NEQN_USER, NLAGS_USER, &
       MY_TRIM_FREQUENCY, MIN_DROP, UNUM655, UNUM656, ERROR_FLAG

  ! Rootfining arrays will be incremented by this amount.
  INTEGER, PRIVATE, PARAMETER :: ADD_TO_IW = 100
  ! Debugging prints will be written to this ouyput unit.
  INTEGER, PRIVATE, PARAMETER :: DBGUNIT = 10

  LOGICAL, PRIVATE :: CONSTANT_DELAYS, CONSTANT_HISTORY, HAVE_EVENT_FCN, &
       HAVE_OUT_FCN, MYNEUTRAL, DROPZ, SHIFTT, INFO_DROPPED, &
       SOLIPOINT_IS_ALLOCATED, HAVE_TRIM_GET, DEBUG, DONT_CALL_CHANGE
  ! ..
  ! .. Local Arrays ..
  DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: DARRAY(:,:), DTREE(:), GA(:), &
       GB(:), GTNEW(:), GTOLD(:), GTROOT(:), GTROOT2(:), KMAT(:,:), &
       MYABSER(:), MYBVAL(:), MYDYEXTP(:), MYDYNEW(:), MYDYOLD(:), &
       MYDYOUT(:), MYDYROOT(:), MYJUMPS(:), MYR(:), MYRELER(:), MYTHIT(:), &
       MYYEXTP(:), MYYNEW(:), MYYOLD(:), MYYOUT(:), MYYROOT(:), QUEUE(:,:), &
       TG(:), TQUEUE(:), TREE_TEMP(:), VTEMP3(:), ZANDD(:,:), ZARRAY(:,:), &
       ZANDD2(:,:), ZARRAY2(:,:), DARRAY2(:,:), Y0SAVE(:), F0SAVE(:)
  DOUBLE PRECISION, POINTER, PRIVATE :: PASS_DELAYS(:), PASS_HISTORY(:)
  INTEGER, ALLOCATABLE, PRIVATE :: INDEXG(:), LEVEL(:), MYDIRECTION(:), &
       MYJROOT(:), MYREPORT(:), STARTS(:)
  LOGICAL, ALLOCATABLE, PRIVATE :: MYISTERMINAL(:)

  ! Perstistent local variables in subroutine DDE_DRV2 .. 
  DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: PVARD(:)
  INTEGER, ALLOCATABLE, PRIVATE :: PVARI(:)
  LOGICAL, ALLOCATABLE, PRIVATE :: PVARL(:)

  INTEGER, PRIVATE :: MYIPOINT(27)
  INTEGER, PARAMETER :: LIPOINT=27

  ! .. Intrinsic Functions ..
  INTRINSIC EPSILON, HUGE
  ! ..
  ! .. Private Statements ..
  PRIVATE :: CD_BETA, CH_YINIT, DUMMY_CHANGE, DUMMY_GUSER, SOL_OUT, &
       TRIM_GET_DUMMY, DUMMY_BETA, DUMMY_YINIT
  ! ..
  ! .. DDE_SOLVER version
  INTEGER, PRIVATE, PARAMETER :: VERSION_MAJOR=2, VERSION_MINOR=0, VERSION_POINT=0
  ! .. Set VERSION_RELEASED to .TRUE. for a released version.
  ! .. Set it to .FALSE. during development.
  LOGICAL, PRIVATE, PARAMETER :: VERSION_RELEASED=.FALSE.

  ! ..
  ! End of DDE_SOLVER private private variables and arrays.
  !
  ! NB. IPOINT(8)-IPOINT(lipoint) is used to indicate if certain arrays
  ! have been allocated (used in RELEASE_ARRAYS cleanup subroutine):
  !   8 = SOL%T
  !   9 = SOL%Y
  !  10 = SOL%TE
  !  11 = SOL%YE
  !  12 = SOL%QUEUE
  !  13 = SOL%YOFT
  !  14 = SOL%TQUEUE
  !  15 = SOL%IE
  !  16 = SOL%IPOINT
  !  17 = SOL%STATS
  !  18 = YINT%TVALS
  !  19 = YINT%YT
  !  20 = YINT%DT
  !  21 = YINT%COMPONENTS
  !  22 = OPTS%ISTERMINAL
  !  23 = OPTS%DIRECTION
  !  24 = OPTS%ABSERR
  !  25 = OPTS%RELERR
  !  26 = OPTS%JUMPS
  !  27 = OPTS%THIS_EXACTLY
  !____________________________________________________________________________

CONTAINS

  !********Begin subroutines for generic DDE_SOLVER**************
  ! The appropriate DKL_* subroutine will be invoked when DDE_SOLVER
  ! is called.

  FUNCTION DKL_1(NVAR,DDES,BETA,HISTORY,TSPAN,OPTIONS,EVENT_FCN, &
       CHANGE_FCN,OUT_FCN,USER_TRIM_GET) RESULT (SOL)

    ! Both BETA and HISTORY are user supplied subroutines.

    ! .. Function Return Value ..
    TYPE (DDE_SOL), TARGET :: SOL
    ! ..
    ! .. Structure Arguments ..
    TYPE (DDE_OPTS), OPTIONAL :: OPTIONS
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION :: TSPAN(:)
    INTEGER :: NVAR(:)
    ! ..
    ! .. Subroutine Arguments ..
    OPTIONAL :: OUT_FCN, CHANGE_FCN, EVENT_FCN, USER_TRIM_GET
    ! ..
    INTERFACE
       SUBROUTINE DDES(T,Y,Z,DY)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         INTENT(IN):: T,Y,Z
         INTENT(OUT) :: DY
       END SUBROUTINE DDES
    END INTERFACE
   INTERFACE
       SUBROUTINE HISTORY(T,Y)
        DOUBLE PRECISION :: T
        DOUBLE PRECISION, DIMENSION(:) :: Y
        INTENT(IN):: T
        INTENT(OUT) :: Y
      END SUBROUTINE HISTORY
   END INTERFACE
   INTERFACE
      SUBROUTINE BETA(T,Y,BVAL)
        DOUBLE PRECISION :: T
        DOUBLE PRECISION, DIMENSION(:) :: Y
        DOUBLE PRECISION, DIMENSION(:) :: BVAL
        INTENT(IN):: T,Y
        INTENT(OUT) :: BVAL
      END SUBROUTINE BETA
   END INTERFACE
    INTERFACE
       SUBROUTINE CHANGE_FCN(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
            DIRECTION,ISTERMINAL,QUIT)
         INTEGER :: NEVENT
         INTEGER, DIMENSION(:) :: DIRECTION
         DOUBLE PRECISION :: TEVENT,HINIT
         DOUBLE PRECISION, DIMENSION(:) :: YEVENT,DYEVENT
         LOGICAL :: QUIT
         LOGICAL, DIMENSION(:) :: ISTERMINAL
         INTENT(IN) :: NEVENT,TEVENT
         INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,ISTERMINAL,QUIT
       END SUBROUTINE CHANGE_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE EVENT_FCN(T,Y,DYDT,Z,G)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         DOUBLE PRECISION, DIMENSION(:) :: G
         INTENT(IN):: T,Y,DYDT,Z
         INTENT(OUT) :: G
       END SUBROUTINE EVENT_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE OUT_FCN(T,Y,DY,N,NEVENT)
         INTEGER :: N,NEVENT
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
       END SUBROUTINE OUT_FCN
    END INTERFACE
   INTERFACE
      SUBROUTINE USER_TRIM_GET
      END SUBROUTINE USER_TRIM_GET
   END INTERFACE
    ! ..
    ! .. Local Structures ..
    TYPE (DDE_OPTS), TARGET :: OPTS
    ! ..
    ! .. Local Scalars ..
    INTEGER :: I, IQUEUE, IER, IFLAG, NEQN, NGUSER, NLAGS, NOUT, &
         IDROPM1, ISAVE
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC PRESENT, SIZE
    ! ..
    ! .. User subroutines ..
    ! EXTERNAL HISTORY,DDES,BETA,OUT_FCN,CHANGE_FCN,EVENT_FCN, &
    !          USER_TRIM_GET
    ! ..
    SOLIPOINT_IS_ALLOCATED = .FALSE.

    SHIFTT = .FALSE.
    SOL%SHIFT = SHIFTT
    SOL%TSHIFT = 0.0D0
    HAVE_TRIM_GET = PRESENT(USER_TRIM_GET)

    IF (PRESENT(OPTIONS)) THEN
       OPTS = OPTIONS
    ELSE
       OPTS = DDE_SET()
    END IF

    IF (SIZE(NVAR)<2) THEN
       PRINT *, ' You must input the number of equations, NVAR(1)'
       PRINT *, ' and the number of delays, NVAR(2).'
       STOP
    ELSE
       NEQN = NVAR(1)
       NLAGS = NVAR(2)
       NEQN_USER = NEQN
       NLAGS_USER = NLAGS
    END IF

    HAVE_EVENT_FCN = PRESENT(EVENT_FCN)
    IF (HAVE_EVENT_FCN) THEN
       IF (SIZE(NVAR)<3) THEN
          PRINT *, ' You must input the number of event functions, NVAR(3).'
          STOP
       ELSE
          NGUSER = NVAR(3)
       END IF

       ! Provide for 10 events.
       ALLOCATE (SOL%IE(10),SOL%TE(10),SOL%YE(10,NEQN),STAT=IER)
       CALL CHECK_STAT(IER,1)
       MYIPOINT(15) = 1
       MYIPOINT(10) = 1
       MYIPOINT(11) = 1
    ELSE
       NGUSER = 0
    END IF
    SOL%NE = 0

    ! IF (ALLOCATED(YOLD)) THEN
    !    DEALLOCATE(YOLD,STAT=IER)
    !    CALL CHECK_STAT(IER,1204)
    ! END IF
    ! IF (OPTS%NEUTRAL) THEN
    !    ALLOCATE(YOLD(2*NEQN),STAT=IER)
    !    CALL CHECK_STAT(IER,1205)
    ! ELSE
    !    ALLOCATE(YOLD(NEQN),STAT=IER)
    !    CALL CHECK_STAT(IER,1206)
    !    PRINT *, ' YOLD has been allocated.'
    ! END IF

    CALL EXPAND_OPTS(NEQN,NGUSER,OPTS)
    ! The following call was moved to DDE_DRV1.
    ! CALL HISTORY(TSPAN(1),YOLD)
    HAVE_OUT_FCN = PRESENT(OUT_FCN)
    NSPAN = SIZE(TSPAN)
    IF (HAVE_OUT_FCN .AND. NSPAN>2) THEN
       ! Change to output at every step.
       TSPAN = (/ TSPAN(1), TSPAN(NSPAN) /)
       NSPAN = 2
    END IF
    CALL CHECK_TSPAN(TSPAN,NOUT)
    ALLOCATE (SOL%T(NOUT),SOL%Y(NOUT,NEQN),STAT=IER)
    CALL CHECK_STAT(IER,2)
    MYIPOINT(8) = 1
    MYIPOINT(9) = 1
    SOL%NPTS = 0
    PASS_SOL => SOL

    NULLIFY (PASS_DELAYS,PASS_HISTORY)
    CONSTANT_DELAYS = .FALSE.
    CONSTANT_HISTORY = .FALSE.

    IF (PRESENT(CHANGE_FCN)) THEN
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  OUT_FCN,BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  OUT_FCN,BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  SOL_OUT,BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  SOL_OUT,BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    ELSE IF (HAVE_EVENT_FCN) THEN
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  OUT_FCN,BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  OUT_FCN,BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  SOL_OUT,BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  SOL_OUT,BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    ELSE
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  OUT_FCN,BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  OUT_FCN,BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  SOL_OUT,BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  SOL_OUT,BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    END IF

    IFLAG = ERROR_FLAG

    ! Generate an error message if the integration was not successful.
    IF (IFLAG/=0) THEN
       PRINT *, ' One or more errors occurred in DKL_1.'
    END IF

    ! Save the necessary interpolation structure if necessary.
    SOLIPOINT_IS_ALLOCATED = .FALSE.
    IF (IFLAG==0 .AND. OPTS%INTERPOLATION) THEN
       ! Add the interpolation arrays to SOL.
       ! Reset the dimensions of TQUEUE and QUEUE to return only
       ! the portions of TQUEUE and QUEUE which were actually used
       ! (okay to reset since the arrays are about to be deallocated).
       LTQUEUE = MYIPOINT(2)
       MYQCOLS = 10*LTQUEUE
       MYIPOINT(1) = MYQCOLS
       ALLOCATE (SOL%IPOINT(LIPOINT),SOL%YOFT(NEQN),SOL%QUEUE(NEQN,MYQCOLS), &
            SOL%TQUEUE(0:LTQUEUE),STAT=IER)
       CALL CHECK_STAT(IER,3)
       SOLIPOINT_IS_ALLOCATED = .TRUE.
       MYIPOINT(16) = 1
       MYIPOINT(13) = 1
       MYIPOINT(12) = 1
       MYIPOINT(14) = 1
       SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
       SOL%QUEUE(1:NEQN,1:MYQCOLS) = QUEUE(1:NEQN,1:MYQCOLS)
       SOL%TQUEUE(0:LTQUEUE) = TQUEUE(0:LTQUEUE)
    END IF

    IF (HAVE_TRIM_GET) THEN
       ! Write the remaining queue information.
       IQUEUE = MYIPOINT(2)
       DO I = 1, IQUEUE
          ISAVE = I
          IF (TQUEUE(I)>TQLAST) THEN
             IF (DEBUG) THEN
                WRITE(DBGUNIT,9999)
                WRITE(DBGUNIT,9998)
                9999 FORMAT(' About to call TRIM_SAVE following the')
                9998 FORMAT(' completion of the integration.')
             END IF
             CALL TRIM_SAVE(I,IQUEUE)
             CALL USER_TRIM_GET
             GOTO 10
          END IF
       END DO
10     CONTINUE

       IF (DEBUG) THEN
          IDROPM1 = IQUEUE - ISAVE + 1
          WRITE(DBGUNIT,9997) TQUEUE(IQUEUE)
          WRITE(DBGUNIT,9996) IDROPM1
          9997 FORMAT(' TNEW = ', D20.10)
          9996 FORMAT(' Dropped from end of queue = ', I10)
       END IF

       TQLAST = TQUEUE(IQUEUE)
    END IF

    ! Deallocate the remaining DDE_SOLVER arrays.
    DEALLOCATE (QUEUE,TQUEUE,STAT=IER)
    CALL CHECK_STAT(IER,4)

    ! IF (ALLOCATED(YOLD)) THEN
    !    DEALLOCATE(YOLD,STAT=IER)
    !    CALL CHECK_STAT(IER,1207)
    ! END IF

    ! Trim solution structure for output.
    CALL CONTRACT_SOL(PASS_SOL,NEQN)

    ALLOCATE (SOL%STATS(6),STAT=IER)
    CALL CHECK_STAT(IER,3)
    MYIPOINT(17) = 1
    IF (SOLIPOINT_IS_ALLOCATED) THEN
    ELSE
       ALLOCATE (SOL%IPOINT(LIPOINT),STAT=IER)
       CALL CHECK_STAT(IER,3)
       SOLIPOINT_IS_ALLOCATED = .TRUE.
       MYIPOINT(16) = 1
    END IF
    SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
    SOL%FLAG = IFLAG
    SOL%STATS(1) = NSTEPS
    SOL%STATS(2) = NFAILS
    SOL%STATS(3) = NFEVAL
    SOL%STATS(4) = NFAILC
    SOL%STATS(5) = ARRAY_STORAGE
    SOL%STATS(6) = ROOT_FUNCTIONS

    MYIPOINT(1:LIPOINT) = 0

    IF (ASSOCIATED(PASS_SOL)) NULLIFY(PASS_SOL)
    IF (ASSOCIATED(PASS_DELAYS)) NULLIFY(PASS_DELAYS)
    IF (ASSOCIATED(PASS_HISTORY)) NULLIFY(PASS_HISTORY)

    RETURN
  END FUNCTION DKL_1
  !____________________________________________________________________________

  FUNCTION DKL_2(NVAR,DDES,DELAYS,HISTORY,TSPAN,OPTIONS,EVENT_FCN, &
       CHANGE_FCN,OUT_FCN,USER_TRIM_GET) RESULT (SOL)

    ! HISTORY is a user supplied subroutine; and DELAYS (BETA) is a
    ! user supplied vector.

    ! .. Function Return Value ..
    TYPE (DDE_SOL), TARGET :: SOL
    ! ..
    ! .. Structure Arguments ..
    TYPE (DDE_OPTS), OPTIONAL :: OPTIONS
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, TARGET :: DELAYS(:)
    DOUBLE PRECISION :: TSPAN(:)
    INTEGER :: NVAR(:)
    ! ..
    ! .. Subroutine Arguments ..
    OPTIONAL :: OUT_FCN, CHANGE_FCN, EVENT_FCN, USER_TRIM_GET
    !..
    INTERFACE
       SUBROUTINE DDES(T,Y,Z,DY)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         INTENT(IN):: T,Y,Z
         INTENT(OUT) :: DY
       END SUBROUTINE DDES
    END INTERFACE
   INTERFACE
      SUBROUTINE HISTORY(T,Y)
        DOUBLE PRECISION :: T
        DOUBLE PRECISION, DIMENSION(:) :: Y
        INTENT(IN):: T
        INTENT(OUT) :: Y
      END SUBROUTINE HISTORY
   END INTERFACE
    INTERFACE
       SUBROUTINE CHANGE_FCN(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
            DIRECTION,ISTERMINAL,QUIT)
         INTEGER :: NEVENT
         INTEGER, DIMENSION(:) :: DIRECTION
         DOUBLE PRECISION :: TEVENT,HINIT
         DOUBLE PRECISION, DIMENSION(:) :: YEVENT,DYEVENT
         LOGICAL :: QUIT
         LOGICAL, DIMENSION(:) :: ISTERMINAL
         INTENT(IN) :: NEVENT,TEVENT
         INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,ISTERMINAL,QUIT
       END SUBROUTINE CHANGE_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE EVENT_FCN(T,Y,DYDT,Z,G)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         DOUBLE PRECISION, DIMENSION(:) :: G
         INTENT(IN):: T,Y,DYDT,Z
         INTENT(OUT) :: G
       END SUBROUTINE EVENT_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE OUT_FCN(T,Y,DY,N,NEVENT)
         INTEGER :: N,NEVENT
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
       END SUBROUTINE OUT_FCN
    END INTERFACE
   INTERFACE
      SUBROUTINE USER_TRIM_GET
      END SUBROUTINE USER_TRIM_GET
   END INTERFACE
    ! ..
    ! .. Local Structures ..
    TYPE (DDE_OPTS), TARGET :: OPTS
    ! .. Local Scalars ..
    INTEGER :: I, IQUEUE, IER, IFLAG, NEQN, NGUSER, NLAGS, NOUT, &
         IDROPM1, ISAVE
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC PRESENT, SIZE
    ! ..
    ! .. User subroutines ..
    ! EXTERNAL HISTORY,DDES,OUT_FCN,CHANGE_FCN,EVENT_FCN,USER_TRIM_GET
    ! ..
    SOLIPOINT_IS_ALLOCATED = .FALSE.

    SHIFTT = .FALSE.
    SOL%SHIFT = SHIFTT
    SOL%TSHIFT = 0.0D0

    HAVE_TRIM_GET = PRESENT(USER_TRIM_GET)

    IF (PRESENT(OPTIONS)) THEN
       OPTS = OPTIONS
    ELSE
       OPTS = DDE_SET()
    END IF

    IF (SIZE(NVAR)<2) THEN
       PRINT *, ' You must input the number of equations, NVAR(1)'
       PRINT *, ' and the number of delays, NVAR(2).'
       STOP
    ELSE
       NEQN = NVAR(1)
       NLAGS = NVAR(2)
       NEQN_USER = NEQN
       NLAGS_USER = NLAGS
    END IF

    HAVE_EVENT_FCN = PRESENT(EVENT_FCN)
    IF (HAVE_EVENT_FCN) THEN
       IF (SIZE(NVAR)<3) THEN
          PRINT *, ' You must input the number of event functions, NVAR(3).'
          STOP
       ELSE
          NGUSER = NVAR(3)
       END IF
       ! Provide for 10 events.
       ALLOCATE (SOL%IE(10),SOL%TE(10),SOL%YE(10,NEQN),STAT=IER)
       CALL CHECK_STAT(IER,5)
       MYIPOINT(15) = 1
       MYIPOINT(10) = 1
       MYIPOINT(11) = 1
    ELSE
       NGUSER = 0
    END IF
    SOL%NE = 0

    ! IF (ALLOCATED(YOLD)) THEN
    !    DEALLOCATE(YOLD,STAT=IER)
    !    CALL CHECK_STAT(IER,1204)
    ! END IF
    ! IF (OPTS%NEUTRAL) THEN
    !    ALLOCATE(YOLD(2*NEQN),STAT=IER)
    !    CALL CHECK_STAT(IER,1205)
    ! ELSE
    !    ALLOCATE(YOLD(NEQN),STAT=IER)
    !    CALL CHECK_STAT(IER,1206)
    !    PRINT *, ' YOLD has been allocated.'
    ! END IF

    CALL EXPAND_OPTS(NEQN,NGUSER,OPTS)
    ! The following call was moved to DDE_DRV1.
    ! CALL HISTORY(TSPAN(1),YOLD)
    HAVE_OUT_FCN = PRESENT(OUT_FCN)
    NSPAN = SIZE(TSPAN)
    IF (HAVE_OUT_FCN .AND. NSPAN>2) THEN
       ! Change to output at every step.
       TSPAN = (/ TSPAN(1), TSPAN(NSPAN) /)
       NSPAN = 2
    END IF
    CALL CHECK_TSPAN(TSPAN,NOUT)
    ALLOCATE (SOL%T(NOUT),SOL%Y(NOUT,NEQN),STAT=IER)
    CALL CHECK_STAT(IER,6)
    MYIPOINT(8) = 1
    MYIPOINT(9) = 1
    SOL%NPTS = 0
    PASS_SOL => SOL

    NULLIFY (PASS_DELAYS,PASS_HISTORY)
    PASS_DELAYS => DELAYS
    CONSTANT_DELAYS = .TRUE.
    CONSTANT_HISTORY = .FALSE.

    HAVE_OUT_FCN = PRESENT(OUT_FCN)
    IF (PRESENT(CHANGE_FCN)) THEN
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  OUT_FCN,DUMMY_BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
            DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  OUT_FCN,DUMMY_BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  SOL_OUT,DUMMY_BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  SOL_OUT,DUMMY_BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    ELSE IF (HAVE_EVENT_FCN) THEN
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    ELSE
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,HISTORY,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,HISTORY,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    END IF

    IFLAG = ERROR_FLAG

    ! Generate an error message if the integration was not successful.
    IF (IFLAG/=0) THEN
       PRINT *, ' One or more errors occurred in DKL_2.'
    END IF

    ! Save the necessary interpolation structure if necessary.
    SOLIPOINT_IS_ALLOCATED = .FALSE.
    IF (IFLAG==0 .AND. OPTS%INTERPOLATION) THEN
       ! Add the interpolation arrays to SOL.
       ! Reset the dimensions of TQUEUE and QUEUE to return only
       ! the portions of TQUEUE and QUEUE which were actually used
       ! (okay to reset since the arrays are about to be deallocated).
       LTQUEUE = MYIPOINT(2)
       MYQCOLS = 10*LTQUEUE
       MYIPOINT(1) = MYQCOLS
       ALLOCATE (SOL%IPOINT(LIPOINT),SOL%YOFT(NEQN),SOL%QUEUE(NEQN,MYQCOLS), &
            SOL%TQUEUE(0:LTQUEUE),STAT=IER)
       CALL CHECK_STAT(IER,7)
       SOLIPOINT_IS_ALLOCATED = .TRUE.
       MYIPOINT(16) = 1
       MYIPOINT(13) = 1
       MYIPOINT(12) = 1
       MYIPOINT(14) = 1
       SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
       SOL%QUEUE(1:NEQN,1:MYQCOLS) = QUEUE(1:NEQN,1:MYQCOLS)
       SOL%TQUEUE(0:LTQUEUE) = TQUEUE(0:LTQUEUE)
    END IF

    IF (HAVE_TRIM_GET) THEN
       ! Write the remaining queue information.
       IQUEUE = MYIPOINT(2)
       DO I = 1, IQUEUE
          ISAVE = I
          IF (TQUEUE(I)>TQLAST) THEN
                  
             IF (DEBUG) THEN
                WRITE(DBGUNIT,9999)
                WRITE(DBGUNIT,9998)
                9999 FORMAT(' About to call TRIM_SAVE following the')
                9998 FORMAT(' completion of the integration.')
             END IF
             
             CALL TRIM_SAVE(I,IQUEUE)
             CALL USER_TRIM_GET
             GOTO 10
          END IF
       END DO
10     CONTINUE

       IF (DEBUG) THEN
          IDROPM1 = IQUEUE - ISAVE + 1
          WRITE(DBGUNIT,9997) TQUEUE(IQUEUE)
          WRITE(DBGUNIT,9996) IDROPM1
          9997 FORMAT(' TNEW = ', D20.10)
          9996 FORMAT(' Dropped from end of queue = ', I10)
       END IF

       TQLAST = TQUEUE(IQUEUE)

    END IF

    ! Deallocate the remaining DDE_SOLVER arrays.
    DEALLOCATE (QUEUE,TQUEUE,STAT=IER)
    CALL CHECK_STAT(IER,8)

    ! IF (ALLOCATED(YOLD)) THEN
    !    DEALLOCATE(YOLD,STAT=IER)
    !    CALL CHECK_STAT(IER,1203)
    ! END IF

    ! Trim solution structure for output.
    CALL CONTRACT_SOL(PASS_SOL,NEQN)

    ALLOCATE (SOL%STATS(6),STAT=IER)
    CALL CHECK_STAT(IER,3)
    MYIPOINT(17) = 1
    IF (SOLIPOINT_IS_ALLOCATED) THEN
    ELSE
       ALLOCATE (SOL%IPOINT(LIPOINT),STAT=IER)
       CALL CHECK_STAT(IER,3)
       SOLIPOINT_IS_ALLOCATED = .TRUE.
       MYIPOINT(16) = 1
    END IF
    SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
    SOL%FLAG = IFLAG
    SOL%STATS(1) = NSTEPS
    SOL%STATS(2) = NFAILS
    SOL%STATS(3) = NFEVAL
    SOL%STATS(4) = NFAILC
    SOL%STATS(5) = ARRAY_STORAGE
    SOL%STATS(6) = ROOT_FUNCTIONS

    MYIPOINT(1:LIPOINT) = 0

    IF (ASSOCIATED(PASS_SOL)) NULLIFY(PASS_SOL)
    IF (ASSOCIATED(PASS_DELAYS)) NULLIFY(PASS_DELAYS)
    IF (ASSOCIATED(PASS_HISTORY)) NULLIFY(PASS_HISTORY)

    RETURN
  END FUNCTION DKL_2
  !____________________________________________________________________________
  FUNCTION DKL_3(NVAR,DDES,BETA,HISTORY,TSPAN,OPTIONS,EVENT_FCN, &
       CHANGE_FCN,OUT_FCN,USER_TRIM_GET) RESULT (SOL)

    ! BETA is a user supplied subroutine and HISTORY is a user
    ! supplied vector.

    ! .. Function Return Value ..
    TYPE (DDE_SOL), TARGET :: SOL
    ! ..
    ! .. Structure Arguments ..
    TYPE (DDE_OPTS), OPTIONAL :: OPTIONS
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, TARGET :: HISTORY(:)
    DOUBLE PRECISION :: TSPAN(:)
    INTEGER :: NVAR(:)
    ! ..
    ! .. Subroutine Arguments ..
    OPTIONAL :: OUT_FCN, CHANGE_FCN, EVENT_FCN, USER_TRIM_GET
    ! ..
   INTERFACE
      SUBROUTINE BETA(T,Y,BVAL)
        DOUBLE PRECISION :: T
        DOUBLE PRECISION, DIMENSION(:) :: Y
        DOUBLE PRECISION, DIMENSION(:) :: BVAL
        INTENT(IN):: T,Y
        INTENT(OUT) :: BVAL
      END SUBROUTINE BETA
   END INTERFACE
    INTERFACE
       SUBROUTINE DDES(T,Y,Z,DY)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         INTENT(IN):: T,Y,Z
         INTENT(OUT) :: DY
       END SUBROUTINE DDES
    END INTERFACE
    INTERFACE
       SUBROUTINE CHANGE_FCN(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
            DIRECTION,ISTERMINAL,QUIT)
         INTEGER :: NEVENT
         INTEGER, DIMENSION(:) :: DIRECTION
         DOUBLE PRECISION :: TEVENT,HINIT
         DOUBLE PRECISION, DIMENSION(:) :: YEVENT,DYEVENT
         LOGICAL :: QUIT
         LOGICAL, DIMENSION(:) :: ISTERMINAL
         INTENT(IN) :: NEVENT,TEVENT
         INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,ISTERMINAL,QUIT
       END SUBROUTINE CHANGE_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE EVENT_FCN(T,Y,DYDT,Z,G)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         DOUBLE PRECISION, DIMENSION(:) :: G
         INTENT(IN):: T,Y,DYDT,Z
         INTENT(OUT) :: G
       END SUBROUTINE EVENT_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE OUT_FCN(T,Y,DY,N,NEVENT)
         INTEGER :: N,NEVENT
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
       END SUBROUTINE OUT_FCN
    END INTERFACE
   INTERFACE
      SUBROUTINE USER_TRIM_GET
      END SUBROUTINE USER_TRIM_GET
   END INTERFACE
    ! ..
    ! .. Local Structures ..
    TYPE (DDE_OPTS), TARGET :: OPTS
    ! ..
    ! .. Local Scalars ..
    INTEGER :: I, IQUEUE, IER, IFLAG, NEQN, NGUSER, NLAGS, NOUT, &
         IDROPM1, ISAVE
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC PRESENT, SIZE
    ! ..
    ! .. User subroutines ..
    ! EXTERNAL DDES,BETA,OUT_FCN,CHANGE_FCN,EVENT_FCN,USER_TRIM_GET
    ! ..
    SOLIPOINT_IS_ALLOCATED = .FALSE.

    SHIFTT = .FALSE.
    SOL%SHIFT = SHIFTT
    SOL%TSHIFT = 0.0D0

    HAVE_TRIM_GET = PRESENT(USER_TRIM_GET)

    IF (PRESENT(OPTIONS)) THEN
       OPTS = OPTIONS
    ELSE
       OPTS = DDE_SET()
    END IF

    IF (SIZE(NVAR)<2) THEN
       PRINT *, ' You must input the number of equations, NVAR(1)'
       PRINT *, ' and the number of delays, NVAR(2).'
       STOP
    ELSE
       NEQN = NVAR(1)
       NLAGS = NVAR(2)
       NEQN_USER = NEQN
       NLAGS_USER = NLAGS
    END IF

    HAVE_EVENT_FCN = PRESENT(EVENT_FCN)
    IF (HAVE_EVENT_FCN) THEN
       IF (SIZE(NVAR)<3) THEN
          PRINT *, ' You must input the number of event functions, NVAR(3).'
          STOP
       ELSE
          NGUSER = NVAR(3)
       END IF
       ! Provide for 10 events.
       ALLOCATE (SOL%IE(10),SOL%TE(10),SOL%YE(10,NEQN),STAT=IER)
       CALL CHECK_STAT(IER,9)
       MYIPOINT(15) = 1
       MYIPOINT(10) = 1
       MYIPOINT(11) = 1
    ELSE
       NGUSER = 0
    END IF
    SOL%NE = 0

    ! IF (ALLOCATED(YOLD)) THEN
    !    DEALLOCATE(YOLD,STAT=IER)
    !    CALL CHECK_STAT(IER,1204)
    ! END IF
    ! IF (OPTS%NEUTRAL) THEN
    !    ALLOCATE(YOLD(2*NEQN),STAT=IER)
    !    CALL CHECK_STAT(IER,1205)
    ! ELSE
    !    ALLOCATE(YOLD(NEQN),STAT=IER)
    !    CALL CHECK_STAT(IER,1206)
    !    PRINT *, ' YOLD has been allocated.'
    ! END IF

    CALL EXPAND_OPTS(NEQN,NGUSER,OPTS)
    ! The following was moved to DDE_DRV1.
    ! YOLD(1:NVAR(1)) = HISTORY
    HAVE_OUT_FCN = PRESENT(OUT_FCN)
    NSPAN = SIZE(TSPAN)
    IF (HAVE_OUT_FCN .AND. NSPAN>2) THEN
       ! Change to output at every step.
       TSPAN = (/ TSPAN(1), TSPAN(NSPAN) /)
       NSPAN = 2
    END IF
    CALL CHECK_TSPAN(TSPAN,NOUT)
    ALLOCATE (SOL%T(NOUT),SOL%Y(NOUT,NEQN),STAT=IER)
    CALL CHECK_STAT(IER,10)
    MYIPOINT(8) = 1
    MYIPOINT(9) = 1
    SOL%NPTS = 0
    PASS_SOL => SOL

    NULLIFY (PASS_DELAYS,PASS_HISTORY)
    CONSTANT_DELAYS = .FALSE.
    PASS_HISTORY => HISTORY
    CONSTANT_HISTORY = .TRUE.

    HAVE_OUT_FCN = PRESENT(OUT_FCN)
    IF (PRESENT(CHANGE_FCN)) THEN
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  OUT_FCN,BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  OUT_FCN,BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  SOL_OUT,BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  SOL_OUT,BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    ELSE IF (HAVE_EVENT_FCN) THEN
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  OUT_FCN,BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  OUT_FCN,BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  SOL_OUT,BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  SOL_OUT,BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    ELSE
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  OUT_FCN,BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  OUT_FCN,BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  SOL_OUT,BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  SOL_OUT,BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    END IF

    IFLAG = ERROR_FLAG

    ! Generate an error message if the integration was not successful.
    IF (IFLAG/=0) THEN
       PRINT *, ' One or more errors occurred in DKL_3.'
    END IF

    ! Save the necessary interpolation structure if necessary.
    SOLIPOINT_IS_ALLOCATED = .FALSE.
    IF (IFLAG==0 .AND. OPTS%INTERPOLATION) THEN
       ! Add the interpolation arrays to SOL.
       ! Reset the dimensions of TQUEUE and QUEUE to return only
       ! the portions of TQUEUE and QUEUE which were actually used
       ! (okay to reset since the arrays are about to be deallocated).
       LTQUEUE = MYIPOINT(2)
       MYQCOLS = 10*LTQUEUE
       MYIPOINT(1) = MYQCOLS
       ALLOCATE (SOL%IPOINT(LIPOINT),SOL%YOFT(NEQN),SOL%QUEUE(NEQN,MYQCOLS), &
            SOL%TQUEUE(0:LTQUEUE),STAT=IER)
       CALL CHECK_STAT(IER,11)
       SOLIPOINT_IS_ALLOCATED = .TRUE.
       MYIPOINT(16) = 1
       MYIPOINT(15) = 1
       MYIPOINT(13) = 1
       MYIPOINT(12) = 1
       MYIPOINT(14) = 1
       SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
       SOL%QUEUE(1:NEQN,1:MYQCOLS) = QUEUE(1:NEQN,1:MYQCOLS)
       SOL%TQUEUE(0:LTQUEUE) = TQUEUE(0:LTQUEUE)
    END IF

    IF (HAVE_TRIM_GET) THEN
       ! Write the remaining queue information.
       IQUEUE = MYIPOINT(2)
       DO I = 1, IQUEUE
          ISAVE = I
          IF (TQUEUE(I)>TQLAST) THEN
             IF (DEBUG) THEN
                WRITE(DBGUNIT,9999)
                WRITE(DBGUNIT,9998)
                9999 FORMAT(' About to call TRIM_SAVE following the')
                9998 FORMAT(' completion of the integration.')
             END IF
             CALL TRIM_SAVE(I,IQUEUE)
             CALL USER_TRIM_GET
             GOTO 10
          END IF
       END DO
10     CONTINUE

       IF (DEBUG) THEN
          IDROPM1 = IQUEUE - ISAVE + 1
          WRITE(DBGUNIT,9997) TQUEUE(IQUEUE)
          WRITE(DBGUNIT,9996) IDROPM1
          9997 FORMAT(' TNEW = ', D20.10)
          9996 FORMAT(' Dropped from end of queue = ', I10)
       END IF

       TQLAST = TQUEUE(IQUEUE)

    END IF

    ! Deallocate the remaining DDE_SOLVER arrays.
    DEALLOCATE (QUEUE,TQUEUE,STAT=IER)
    CALL CHECK_STAT(IER,12)

    ! IF (ALLOCATED(YOLD)) THEN
    !    DEALLOCATE(YOLD,STAT=IER)
    !    CALL CHECK_STAT(IER,1211)
    ! END IF

    ! Trim solution structure for output.
    CALL CONTRACT_SOL(PASS_SOL,NEQN)

    ALLOCATE (SOL%STATS(6),STAT=IER)
    CALL CHECK_STAT(IER,3)
    MYIPOINT(17) = 1
    IF (SOLIPOINT_IS_ALLOCATED) THEN
    ELSE
       ALLOCATE (SOL%IPOINT(LIPOINT),STAT=IER)
       CALL CHECK_STAT(IER,3)
       SOLIPOINT_IS_ALLOCATED = .TRUE.
       MYIPOINT(16) = 1
    END IF
    SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
    SOL%FLAG = IFLAG
    SOL%STATS(1) = NSTEPS
    SOL%STATS(2) = NFAILS
    SOL%STATS(3) = NFEVAL
    SOL%STATS(4) = NFAILC
    SOL%STATS(5) = ARRAY_STORAGE
    SOL%STATS(6) = ROOT_FUNCTIONS

    MYIPOINT(1:LIPOINT) = 0

    IF (ASSOCIATED(PASS_SOL)) NULLIFY(PASS_SOL)
    IF (ASSOCIATED(PASS_DELAYS)) NULLIFY(PASS_DELAYS)
    IF (ASSOCIATED(PASS_HISTORY)) NULLIFY(PASS_HISTORY)

    RETURN
  END FUNCTION DKL_3
  !____________________________________________________________________________
  FUNCTION DKL_4(NVAR,DDES,DELAYS,HISTORY,TSPAN,OPTIONS,EVENT_FCN, &
       CHANGE_FCN,OUT_FCN,USER_TRIM_GET) RESULT (SOL)

    ! Both BETA and HISTORY are user supplied vectors.

    ! .. Function Return Value ..
    TYPE (DDE_SOL), TARGET :: SOL
    ! ..
    ! .. Structure Arguments ..
    TYPE (DDE_OPTS), OPTIONAL :: OPTIONS
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, TARGET :: DELAYS(:), HISTORY(:)
    DOUBLE PRECISION :: TSPAN(:)
    INTEGER :: NVAR(:)
    ! ..
    ! .. Subroutine Arguments ..
    OPTIONAL :: OUT_FCN, CHANGE_FCN, EVENT_FCN, USER_TRIM_GET
    ! ..
    INTERFACE
       SUBROUTINE DDES(T,Y,Z,DY)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         INTENT(IN):: T,Y,Z
         INTENT(OUT) :: DY
       END SUBROUTINE DDES
    END INTERFACE
    INTERFACE
       SUBROUTINE CHANGE_FCN(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
            DIRECTION,ISTERMINAL,QUIT)
         INTEGER :: NEVENT
         INTEGER, DIMENSION(:) :: DIRECTION
         DOUBLE PRECISION :: TEVENT,HINIT
         DOUBLE PRECISION, DIMENSION(:) :: YEVENT,DYEVENT
         LOGICAL :: QUIT
         LOGICAL, DIMENSION(:) :: ISTERMINAL
         INTENT(IN) :: NEVENT,TEVENT
         INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,ISTERMINAL,QUIT
       END SUBROUTINE CHANGE_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE EVENT_FCN(T,Y,DYDT,Z,G)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         DOUBLE PRECISION, DIMENSION(:) :: G
         INTENT(IN):: T,Y,DYDT,Z
         INTENT(OUT) :: G
       END SUBROUTINE EVENT_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE OUT_FCN(T,Y,DY,N,NEVENT)
         INTEGER :: N,NEVENT
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
       END SUBROUTINE OUT_FCN
    END INTERFACE
   INTERFACE
      SUBROUTINE USER_TRIM_GET
      END SUBROUTINE USER_TRIM_GET
   END INTERFACE
    ! ..
    ! .. Local Structures ..
    TYPE (DDE_OPTS), TARGET :: OPTS
    ! ..
    ! .. Local Scalars ..
    INTEGER :: I, IQUEUE, IER, IFLAG, NEQN, NGUSER, NLAGS, NOUT, &
         IDROPM1, ISAVE
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC PRESENT, SIZE
    ! ..
    ! .. User subroutines ..
    ! EXTERNAL DDES,OUT_FCN,CHANGE_FCN,EVENT_FCN,USER_TRIM_GET
    ! ..
    SOLIPOINT_IS_ALLOCATED = .FALSE.

    SHIFTT = .FALSE.
    SOL%SHIFT = SHIFTT
    SOL%TSHIFT = 0.0D0

    HAVE_TRIM_GET = PRESENT(USER_TRIM_GET)

    IF (PRESENT(OPTIONS)) THEN
       OPTS = OPTIONS
    ELSE
       OPTS = DDE_SET()
    END IF

    IF (SIZE(NVAR)<2) THEN
       PRINT *, ' You must input the number of equations, NVAR(1)'
       PRINT *, ' and the number of delays, NVAR(2).'
       STOP
    ELSE
       NEQN = NVAR(1)
       NLAGS = NVAR(2)
       NEQN_USER = NEQN
       NLAGS_USER = NLAGS
    END IF

    HAVE_EVENT_FCN = PRESENT(EVENT_FCN)
    IF (HAVE_EVENT_FCN) THEN
       IF (SIZE(NVAR)<3) THEN
          PRINT *, ' You must input the number of event functions, NVAR(3).'
          STOP
       ELSE
          NGUSER = NVAR(3)
       END IF
       ! Provide for 10 events.
       ALLOCATE (SOL%IE(10),SOL%TE(10),SOL%YE(10,NEQN),STAT=IER)
       CALL CHECK_STAT(IER,13)
       MYIPOINT(15) = 1
       MYIPOINT(10) = 1
       MYIPOINT(11) = 1
    ELSE
       NGUSER = 0
    END IF
    SOL%NE = 0
    CALL EXPAND_OPTS(NEQN,NGUSER,OPTS)

    ! IF (ALLOCATED(YOLD)) THEN
    !    DEALLOCATE(YOLD,STAT=IER)
    !    CALL CHECK_STAT(IER,1204)
    ! END IF
    ! IF (OPTS%NEUTRAL) THEN
    !    ALLOCATE(YOLD(2*NEQN),STAT=IER)
    !    CALL CHECK_STAT(IER,1205)
    ! ELSE
    !    ALLOCATE(YOLD(NEQN),STAT=IER)
    !    CALL CHECK_STAT(IER,1206)
    !    PRINT *, ' YOLD has been allocated.'
    ! END IF

    ! The following was moved to DDE_DRV1.
    ! YOLD(1:NVAR(1)) = HISTORY
    HAVE_OUT_FCN = PRESENT(OUT_FCN)
    NSPAN = SIZE(TSPAN)
    IF (HAVE_OUT_FCN .AND. NSPAN>2) THEN
       ! Change to output at every step.
       TSPAN = (/ TSPAN(1), TSPAN(NSPAN) /)
       NSPAN = 2
    END IF
    CALL CHECK_TSPAN(TSPAN,NOUT)
    ALLOCATE (SOL%T(NOUT),SOL%Y(NOUT,NEQN),STAT=IER)
    CALL CHECK_STAT(IER,14)
    MYIPOINT(8) = 1
    MYIPOINT(9) = 1
    SOL%NPTS = 0
    PASS_SOL => SOL

    NULLIFY (PASS_DELAYS,PASS_HISTORY)
    PASS_DELAYS => DELAYS
    CONSTANT_DELAYS = .TRUE.
    PASS_HISTORY => HISTORY
    CONSTANT_HISTORY = .TRUE.

    HAVE_OUT_FCN = PRESENT(OUT_FCN)
    IF (PRESENT(CHANGE_FCN)) THEN
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,CHANGE_FCN, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    ELSE IF (HAVE_EVENT_FCN) THEN
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
              DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,EVENT_FCN,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    ELSE
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,DDES,DUMMY_GUSER,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    END IF

    IFLAG = ERROR_FLAG

    ! Generate an error message if the integration was not successful.
    IF (IFLAG/=0) THEN
       PRINT *, ' One or more errors occurred in DKL_4.'
    END IF

    ! Save the necessary interpolation structure if necessary.
    SOLIPOINT_IS_ALLOCATED = .FALSE.
    IF (IFLAG==0 .AND. OPTS%INTERPOLATION) THEN
       ! Add the interpolation arrays to SOL.
       ! Reset the dimensions of TQUEUE and QUEUE to return only
       ! the portions of TQUEUE and QUEUE which were actually used
       ! (okay to reset since the arrays are about to be deallocated).
       LTQUEUE = MYIPOINT(2)
       MYQCOLS = 10*LTQUEUE
       MYIPOINT(1) = MYQCOLS
       ALLOCATE (SOL%IPOINT(LIPOINT),SOL%YOFT(NEQN),SOL%QUEUE(NEQN,MYQCOLS), &
            SOL%TQUEUE(0:LTQUEUE),STAT=IER)
       CALL CHECK_STAT(IER,15)
       SOLIPOINT_IS_ALLOCATED = .TRUE.
       MYIPOINT(16) = 1
       MYIPOINT(13) = 1
       MYIPOINT(12) = 1
       MYIPOINT(14) = 1
       SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
       SOL%QUEUE(1:NEQN,1:MYQCOLS) = QUEUE(1:NEQN,1:MYQCOLS)
       SOL%TQUEUE(0:LTQUEUE) = TQUEUE(0:LTQUEUE)
    END IF

    IF (HAVE_TRIM_GET) THEN
       ! Write the remaining queue information.
       IQUEUE = MYIPOINT(2)
       DO I = 1, IQUEUE
          ISAVE = I
          IF (TQUEUE(I)>TQLAST) THEN

             IF (DEBUG) THEN
                WRITE(DBGUNIT,9999)
                WRITE(DBGUNIT,9998)
                9999 FORMAT(' About to call TRIM_SAVE following the')
                9998 FORMAT(' completion of the integration.')
             END IF

             CALL TRIM_SAVE(I,IQUEUE)
             CALL USER_TRIM_GET
             GOTO 10
          END IF
       END DO
10     CONTINUE

       IF (DEBUG) THEN
          IDROPM1 = IQUEUE - ISAVE + 1
          WRITE(DBGUNIT,9997) TQUEUE(IQUEUE)
          WRITE(DBGUNIT,9996) IDROPM1
          9997 FORMAT(' TNEW = ', D20.10)
          9996 FORMAT(' Dropped from end of queue = ', I10)
       END IF
 
       TQLAST = TQUEUE(IQUEUE)

    END IF

    ! Deallocate the remaining DDE_SOLVER arrays.
    DEALLOCATE (QUEUE,TQUEUE,STAT=IER)
    CALL CHECK_STAT(IER,16)

    ! IF (ALLOCATED(YOLD)) THEN
    !    DEALLOCATE(YOLD,STAT=IER)
    !    CALL CHECK_STAT(IER,1216)
    ! END IF

    ! Trim solution structure for output.
    CALL CONTRACT_SOL(PASS_SOL,NEQN)
    ALLOCATE (SOL%STATS(6),STAT=IER)
    CALL CHECK_STAT(IER,3)
    MYIPOINT(17) = 1
    IF (SOLIPOINT_IS_ALLOCATED) THEN
    ELSE
       ALLOCATE (SOL%IPOINT(LIPOINT),STAT=IER)
       SOLIPOINT_IS_ALLOCATED = .TRUE.
       MYIPOINT(16) = 1
       CALL CHECK_STAT(IER,3)
    END IF
    SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
    SOL%FLAG = IFLAG
    SOL%STATS(1) = NSTEPS
    SOL%STATS(2) = NFAILS
    SOL%STATS(3) = NFEVAL
    SOL%STATS(4) = NFAILC
    SOL%STATS(5) = ARRAY_STORAGE
    SOL%STATS(6) = ROOT_FUNCTIONS

    MYIPOINT(1:LIPOINT) = 0

    IF (ASSOCIATED(PASS_SOL)) NULLIFY(PASS_SOL)
    IF (ASSOCIATED(PASS_DELAYS)) NULLIFY(PASS_DELAYS)
    IF (ASSOCIATED(PASS_HISTORY)) NULLIFY(PASS_HISTORY)

    RETURN
  END FUNCTION DKL_4
  !____________________________________________________________________________

  !**********END FUNCTIONS FOR GENERIC DDE_SOLVER**************

  !**********BEGIN PRIVATE AUXILIARY SUBROUTINES***************

  SUBROUTINE DUMMY_BETA(T,Y,BVAL)
     DOUBLE PRECISION, INTENT(IN) :: T
     DOUBLE PRECISION, INTENT(IN) :: Y(:)
     DOUBLE PRECISION, INTENT(OUT) :: BVAL(:)
     RETURN
  END SUBROUTINE DUMMY_BETA

  SUBROUTINE DUMMY_YINIT(T,Y)
     DOUBLE PRECISION, INTENT(IN) :: T
     DOUBLE PRECISION, INTENT(OUT) :: Y(:)
     RETURN
  END SUBROUTINE DUMMY_YINIT

  SUBROUTINE DUMMY_CHANGE(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT,DIRECTION, &
       ISTERMINAL,QUIT)

    ! .. Scalar Arguments ..
    DOUBLE PRECISION :: HINIT, TEVENT, A
    INTEGER :: NEVENT, I
    LOGICAL :: QUIT, GETRIDOF
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION :: DYEVENT(:), YEVENT(:)
    INTEGER :: DIRECTION(:)
    LOGICAL :: ISTERMINAL(:)
    INTENT(IN) :: NEVENT, TEVENT
    INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,ISTERMINAL,QUIT
    ! ..
    ! Get rid of some compiler warning messages about unused variables
    GETRIDOF = .FALSE.
    IF (GETRIDOF) THEN
       I = NEVENT
       QUIT = ISTERMINAL(1)
       A = TEVENT
       A = YEVENT(1)
       A = DYEVENT(1)
       A = HINIT
       A = DIRECTION(1)
       !NEVENT = I
       HINIT = A
    END IF

    RETURN
  END SUBROUTINE DUMMY_CHANGE
  !____________________________________________________________________________

  SUBROUTINE DUMMY_GUSER(T,Y,DYDT,Z,G)
    ! ..
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    DOUBLE PRECISION, DIMENSION(:) :: G
    LOGICAL DUMMY
    INTENT(IN):: T,Y,DYDT,Z
    INTENT(OUT) :: G

    ! Get rid of some needless compiler warning messages.
    DUMMY = .FALSE.
    IF (DUMMY) THEN
       G(1) = T
       G(1) = Y(1)
       G(1) = DYDT(1)
       G(1) = Z(1,1)
    END IF

    RETURN
  END SUBROUTINE DUMMY_GUSER
  !____________________________________________________________________________
  
  SUBROUTINE TRIM_GET_DUMMY

    ! Dummy routine called if the user does not specify a TRIM_GET subroutine.

    RETURN
  END SUBROUTINE TRIM_GET_DUMMY
  !____________________________________________________________________________

  SUBROUTINE CD_BETA(T,BVAL,NLAGS)

    ! .. Scalar Arguments ..
    DOUBLE PRECISION :: T
    INTEGER :: NLAGS
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION :: BVAL(NLAGS)
    ! ..
    BVAL(1:NLAGS) = T - PASS_DELAYS(1:NLAGS)

    RETURN
  END SUBROUTINE CD_BETA
  !____________________________________________________________________________
  SUBROUTINE CH_YINIT(Y,NEQN)

    ! .. Scalar Arguments ..
    INTEGER :: NEQN
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION :: Y(NEQN)
    ! ..
    Y = PASS_HISTORY

    RETURN
  END SUBROUTINE CH_YINIT
  !_______________________________________________________________________

  SUBROUTINE SOL_OUT(T,Y,DY,NEQN,NEVENT)

    ! .. Scalar Arguments ..
    DOUBLE PRECISION :: T
    INTEGER :: NEQN, NEVENT
    ! ..
    ! .. Array Arguments ..
    !DOUBLE PRECISION :: DY(NEQN), Y(NEQN)
    DOUBLE PRECISION :: DY(:), Y(:)
    ! ..
    ! .. Local Scalars ..
    INTEGER :: NE, NPTS
    LOGICAL GETRIDOF
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC SIZE
    ! ..
    IF (NEVENT>0) THEN
       NE = PASS_SOL%NE + 1
       IF (NE>SIZE(PASS_SOL%TE)) THEN
          CALL EXPAND_SOL_EVENTS(PASS_SOL,1000) !10)
       END IF
       PASS_SOL%NE = NE
       PASS_SOL%TE(NE) = T
       PASS_SOL%YE(NE,:) = Y
       PASS_SOL%IE(NE) = NEVENT
    ELSE
       NPTS = PASS_SOL%NPTS + 1
       IF (NPTS>SIZE(PASS_SOL%T)) THEN
          CALL EXPAND_SOL(PASS_SOL,1000) ! 100
       END IF
       PASS_SOL%NPTS = NPTS
       PASS_SOL%T(NPTS) = T
       PASS_SOL%Y(NPTS,:) = Y
    END IF
    ! Get rid of compiler warning message about DY not being used.
    GETRIDOF = .FALSE.
    IF (GETRIDOF) THEN
       DY(1) = T
    END IF

    RETURN
  END SUBROUTINE SOL_OUT
  !____________________________________________________________________________
  SUBROUTINE CHECK_STAT_WARN(IER,CALLED_FROM)

    ! .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: CALLED_FROM, IER
    ! ..
    IF (IER/=0) THEN
       PRINT *, ' A storage allocation error occurred at code location ', &
            CALLED_FROM
    END IF

    RETURN
  END SUBROUTINE CHECK_STAT_WARN
  !____________________________________________________________________________

  SUBROUTINE CHECK_STAT(IER,CALLED_FROM)

    ! .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: CALLED_FROM, IER
    ! ..
    IF (IER/=0) THEN
       PRINT *, ' A storage allocation error occurred at code location ', &
            CALLED_FROM
       STOP
    END IF

    RETURN
  END SUBROUTINE CHECK_STAT
  !____________________________________________________________________________

  SUBROUTINE CHECK_TSPAN(TSPAN,NOUT)

    ! .. Scalar Arguments ..
    INTEGER :: NOUT
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION :: TSPAN(:)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS
    ! ..
    IF (NSPAN<2) THEN
       PRINT *, ' TSPAN must have at least two entries.'
       STOP
    ELSE IF (ABS(TSPAN(NSPAN)-TSPAN(1))<=0.0D0) THEN
       PRINT *, ' The first and last entries of TSPAN must be different.'
       STOP
    END IF
    IF (NSPAN==2) THEN
       NOUT = 100
    ELSE
       NOUT = NSPAN
    END IF

    RETURN
  END SUBROUTINE CHECK_TSPAN
  !____________________________________________________________________________

  SUBROUTINE EXPAND_SOL(PASS_SOL,M)

    ! .. Structure Arguments ..
    TYPE (DDE_SOL), POINTER :: PASS_SOL
    ! ..
    ! .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: M
    ! ..
    ! .. Local Scalars ..
    INTEGER :: COLS, IER, ROWS
    ! ..
    ! .. Local Arrays ..
    DOUBLE PRECISION, ALLOCATABLE :: TEMP(:), YTEMP(:,:)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC SIZE
    ! ..
    ROWS = SIZE(PASS_SOL%Y,1)
    COLS = SIZE(PASS_SOL%Y,2)
    ALLOCATE (TEMP(ROWS+M),STAT=IER)
    CALL CHECK_STAT(IER,18)
    TEMP(1:ROWS) = PASS_SOL%T
    TEMP(ROWS+1:ROWS+M) = 0D0
    DEALLOCATE (PASS_SOL%T,STAT=IER)
    CALL CHECK_STAT(IER,19)
    ALLOCATE (PASS_SOL%T(ROWS+M),STAT=IER)
    MYIPOINT(8) = 1
    CALL CHECK_STAT(IER,20)
    PASS_SOL%T = TEMP
    DEALLOCATE (TEMP,STAT=IER)
    CALL CHECK_STAT(IER,21)

    ALLOCATE (YTEMP(ROWS+M,COLS),STAT=IER)
    CALL CHECK_STAT(IER,22)
    YTEMP(1:ROWS,1:COLS) = PASS_SOL%Y
    YTEMP(ROWS+1:ROWS+M,1:COLS) = 0D0
    DEALLOCATE (PASS_SOL%Y,STAT=IER)
    CALL CHECK_STAT(IER,23)
    ALLOCATE (PASS_SOL%Y(ROWS+M,COLS),STAT=IER)
    MYIPOINT(9) = 1
    CALL CHECK_STAT(IER,23)
    PASS_SOL%Y = YTEMP
    DEALLOCATE (YTEMP,STAT=IER)
    CALL CHECK_STAT(IER,24)

    RETURN
  END SUBROUTINE EXPAND_SOL
  !____________________________________________________________________________

  SUBROUTINE EXPAND_SOL_EVENTS(PASS_SOL,M)

    ! .. Structure Arguments ..
    TYPE (DDE_SOL), POINTER :: PASS_SOL
    ! ..
    ! .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: M
    ! ..
    ! .. Local Scalars ..
    INTEGER :: COLS, IER, ROWS
    ! ..
    ! .. Local Arrays ..
    DOUBLE PRECISION, ALLOCATABLE :: TEMP(:), YTEMP(:,:)
    INTEGER, ALLOCATABLE :: ITEMP(:)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC SIZE
    ! ..
    ROWS = SIZE(PASS_SOL%YE,1)
    COLS = SIZE(PASS_SOL%YE,2)

    ALLOCATE (ITEMP(ROWS+M),STAT=IER)
    CALL CHECK_STAT(IER,25)
    ITEMP(1:ROWS) = PASS_SOL%IE
    ITEMP(ROWS+1:ROWS+M) = 0
    DEALLOCATE (PASS_SOL%IE,STAT=IER)
    CALL CHECK_STAT(IER,26)
    ALLOCATE (PASS_SOL%IE(ROWS+M),STAT=IER)
    MYIPOINT(15) = 1
    CALL CHECK_STAT(IER,27)
    PASS_SOL%IE = ITEMP
    DEALLOCATE (ITEMP,STAT=IER)
    CALL CHECK_STAT(IER,28)

    ALLOCATE (TEMP(ROWS+M),STAT=IER)
    CALL CHECK_STAT(IER,29)
    TEMP(1:ROWS) = PASS_SOL%TE
    TEMP(ROWS+1:ROWS+M) = 0D0
    DEALLOCATE (PASS_SOL%TE,STAT=IER)
    CALL CHECK_STAT(IER,30)
    ALLOCATE (PASS_SOL%TE(ROWS+M),STAT=IER)
    MYIPOINT(10) = 1
    CALL CHECK_STAT(IER,31)
    PASS_SOL%TE = TEMP
    DEALLOCATE (TEMP,STAT=IER)
    CALL CHECK_STAT(IER,32)

    ALLOCATE (YTEMP(ROWS+M,COLS),STAT=IER)
    CALL CHECK_STAT(IER,33)
    YTEMP(1:ROWS,1:COLS) = PASS_SOL%YE
    YTEMP(ROWS+1:ROWS+M,1:COLS) = 0D0
    DEALLOCATE (PASS_SOL%YE,STAT=IER)
    CALL CHECK_STAT(IER,34)
    ALLOCATE (PASS_SOL%YE(ROWS+M,COLS),STAT=IER)
    CALL CHECK_STAT(IER,35)
    MYIPOINT(11) = 1
    PASS_SOL%YE = YTEMP
    DEALLOCATE (YTEMP,STAT=IER)
    CALL CHECK_STAT(IER,35)

    RETURN
  END SUBROUTINE EXPAND_SOL_EVENTS
  !____________________________________________________________________________

  SUBROUTINE CONTRACT_SOL(PASS_SOL,NEQN)

    ! .. Structure Arguments ..
    TYPE (DDE_SOL), POINTER :: PASS_SOL
    ! ..
    ! .. Local Scalars ..
    INTEGER :: COLS, IER, NE, NPTS, NEQN
    ! ..
    ! .. Local Arrays ..
    DOUBLE PRECISION, ALLOCATABLE :: TEMP(:), YTEMP(:,:)
    INTEGER, ALLOCATABLE :: ITEMP(:)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC SIZE
    ! ..
    IF (HAVE_OUT_FCN) THEN
       DEALLOCATE (PASS_SOL%T,PASS_SOL%Y,STAT=IER)
       CALL CHECK_STAT(IER,36)
       IF (HAVE_EVENT_FCN) THEN
          DEALLOCATE (PASS_SOL%TE,PASS_SOL%YE,PASS_SOL%IE,STAT=IER)
          CALL CHECK_STAT(IER,37)
       END IF
       RETURN
    END IF

    NPTS = PASS_SOL%NPTS
    NE = PASS_SOL%NE
    COLS = SIZE(PASS_SOL%Y,2)

    ALLOCATE (TEMP(NPTS),STAT=IER)
    CALL CHECK_STAT(IER,38)
    TEMP = PASS_SOL%T(1:NPTS)
    DEALLOCATE (PASS_SOL%T,STAT=IER)
    CALL CHECK_STAT(IER,39)
    ALLOCATE (PASS_SOL%T(NPTS),STAT=IER)
    CALL CHECK_STAT(IER,40)
    MYIPOINT(8) = 1
    PASS_SOL%T = TEMP
    DEALLOCATE (TEMP,STAT=IER)
    CALL CHECK_STAT(IER,41)

    ALLOCATE (YTEMP(NPTS,COLS),STAT=IER)
    CALL CHECK_STAT(IER,42)
    YTEMP = PASS_SOL%Y(1:NPTS,:)
    DEALLOCATE (PASS_SOL%Y,STAT=IER)
    CALL CHECK_STAT(IER,43)
    ALLOCATE (PASS_SOL%Y(NPTS,COLS),STAT=IER)
    CALL CHECK_STAT(IER,43)
    MYIPOINT(9) = 1
    PASS_SOL%Y = YTEMP
    CALL CHECK_STAT(IER,44)
    DEALLOCATE (YTEMP,STAT=IER)
    CALL CHECK_STAT(IER,44)

    IF (NE>0) THEN

       ALLOCATE (ITEMP(NE),STAT=IER)
       CALL CHECK_STAT(IER,45)
       ITEMP = PASS_SOL%IE(1:NE)
       DEALLOCATE (PASS_SOL%IE,STAT=IER)
       CALL CHECK_STAT(IER,46)
       ALLOCATE (PASS_SOL%IE(NE),STAT=IER)
       CALL CHECK_STAT(IER,47)
       MYIPOINT(15) = 1
       PASS_SOL%IE = ITEMP
       DEALLOCATE (ITEMP,STAT=IER)
       CALL CHECK_STAT(IER,48)

       ALLOCATE (TEMP(NE),STAT=IER)
       CALL CHECK_STAT(IER,49)
       TEMP = PASS_SOL%TE(1:NE)
       DEALLOCATE (PASS_SOL%TE,STAT=IER)
       CALL CHECK_STAT(IER,50)
       ALLOCATE (PASS_SOL%TE(NE),STAT=IER)
       CALL CHECK_STAT(IER,51)
       MYIPOINT(10) = 1
       PASS_SOL%TE = TEMP
       DEALLOCATE (TEMP,STAT=IER)
       CALL CHECK_STAT(IER,52)

       ALLOCATE (YTEMP(NE,COLS),STAT=IER)
       CALL CHECK_STAT(IER,53)
       YTEMP = PASS_SOL%YE(1:NE,:)
       DEALLOCATE (PASS_SOL%YE,STAT=IER)
       CALL CHECK_STAT(IER,54)
       ALLOCATE (PASS_SOL%YE(NE,COLS),STAT=IER)
       CALL CHECK_STAT(IER,55)
       MYIPOINT(11) = 1
       PASS_SOL%YE = YTEMP
       DEALLOCATE (YTEMP,STAT=IER)
       CALL CHECK_STAT(IER,55)

    END IF

    RETURN
  END SUBROUTINE CONTRACT_SOL
  !____________________________________________________________________________

  SUBROUTINE EXPAND_OPTS(NEQN,NGUSER,OPTS)

    ! .. Structure Arguments ..
    TYPE (DDE_OPTS) :: OPTS
    ! ..
    ! .. Scalar Arguments ..
    INTEGER :: NEQN, NGUSER
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: AE, RE
    INTEGER :: IER, NAE, NRE
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ASSOCIATED, SIZE
    ! ..
    NAE = SIZE(OPTS%ABSERR)
    IF (NAE==1) THEN
       AE = OPTS%ABSERR(1)
       DEALLOCATE (OPTS%ABSERR,STAT=IER)
       CALL CHECK_STAT(IER,56)
       ALLOCATE (OPTS%ABSERR(NAE),STAT=IER)
       MYIPOINT(24) = 1
       CALL CHECK_STAT(IER,57)
       OPTS%ABSERR(1) = AE
    ELSE IF (NAE/=NEQN) THEN
       PRINT *, ' AE_VECTOR must have NEQN components.'
       STOP
    END IF
    NRE = SIZE(OPTS%RELERR)
    IF (NRE==1) THEN
       RE = OPTS%RELERR(1)
       DEALLOCATE (OPTS%RELERR,STAT=IER)
       CALL CHECK_STAT(IER,56)
       ALLOCATE (OPTS%RELERR(NRE),STAT=IER)
       MYIPOINT(25) = 1
       CALL CHECK_STAT(IER,57)
       OPTS%RELERR(1) = RE
    ELSE IF (NRE/=NEQN) THEN
       PRINT *, ' RE_VECTOR must have NEQN components.'
       STOP
    END IF

    IF (.NOT. ASSOCIATED(OPTS%ISTERMINAL)) THEN
       ALLOCATE (OPTS%ISTERMINAL(NGUSER),STAT=IER)
       CALL CHECK_STAT(IER,58)
       IF (NGUSER /= 0) MYIPOINT(22) = 1
       OPTS%ISTERMINAL(1:NGUSER) = .FALSE.
    END IF
    IF (.NOT. ASSOCIATED(OPTS%DIRECTION)) THEN
       ALLOCATE (OPTS%DIRECTION(NGUSER),STAT=IER)
       CALL CHECK_STAT(IER,59)
       IF (NGUSER /= 0) MYIPOINT(23) = 1
       OPTS%DIRECTION(1:NGUSER) = 0
    END IF

    RETURN
  END SUBROUTINE EXPAND_OPTS
  !____________________________________________________________________________

  SUBROUTINE TRIM_QUEUE(TRIM_GET)

    ! Trim unneeded solution history when MAX_DELAY option is used.

    ! .. Subroutine Arguments ..
    !EXTERNAL TRIM_GET
    INTERFACE
       SUBROUTINE USER_TRIM_GET
       END SUBROUTINE USER_TRIM_GET
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: TTEMP
    INTEGER :: I, IDROP, IDROPM1, IQUEUE
    ! ..
    IF (TQUEUE(1)<TQUEUE(0)) THEN
       PRINT *, ' Integration right to left is not permitted.'
       STOP
    END IF
    IQUEUE = MYIPOINT(2)
    TTEMP = TQUEUE(IQUEUE) - MY_MAX_DELAY
    IF (TTEMP>=TQUEUE(1)) THEN
       ! Bracket TTEMP:
       DO I = 1, IQUEUE - 1
          IDROP = I
          IF (TQUEUE(I)<=TTEMP .AND. TTEMP<TQUEUE(I+1)) THEN
             GO TO 10
          END IF
       END DO
       PRINT *, ' MAX_DELAY is not valid in subroutine TRIM_QUEUE.'
       STOP
10     CONTINUE
       ! Do not drop fewer than MIN_DROP points:
       IF (IDROP>=MIN_DROP) THEN
          ! Drop this many points:
          IDROPM1 = IDROP - 1

          IF (DEBUG) THEN
             WRITE(DBGUNIT,9999)
             9999 FORMAT(' About to call TRIM_SAVE')
          END IF

          ! Save the information to be dropped.
          IF (HAVE_TRIM_GET) CALL TRIM_SAVE(1,IDROPM1)
          TQLAST = TQUEUE(IDROPM1)
          IF (HAVE_TRIM_GET) CALL TRIM_GET
          ! Move kept tqueue points left:
          TQUEUE(1:IQUEUE-IDROPM1) = TQUEUE(IDROP:IQUEUE)
          ! Move kept queue columns left:
          QUEUE(1:MYN,1:10*(IQUEUE-IDROPM1)) = QUEUE(1:MYN, &
               10*IDROPM1+1:10*IQUEUE)
          ! Indicate that queue information has been dropped:
          INFO_DROPPED = .TRUE.
          ! Reset the bracketing search interval flags:
          IF (NLAGS>0) THEN
             DO I = 1, NLAGS
                STARTS(I) = MAX(STARTS(I)-IDROPM1,1)
             END DO
          END IF
          ! Reset the number of points in the queue:
          MYIPOINT(2) = IQUEUE - IDROPM1
          IQUEUE = IQUEUE - IDROPM1

          !IF (DEBUG) THEN
          !   PRINT *, ' TNEW = ', TQUEUE(IQUEUE), ' Dropped: ', IDROPM1, &
          ! ' Kept: ', MYIPOINT(2)
          ! END IF

          IF (DEBUG) THEN
             WRITE(DBGUNIT,9998) TQUEUE(IQUEUE), IDROPM1, MYIPOINT(2)
             9998 FORMAT(' TNEW = ', D20.10, ' Dropped = ', I10, &
                         ' KEPT = ', I10) 
          END IF

       END IF
    END IF

    RETURN
  END SUBROUTINE TRIM_QUEUE
  !____________________________________________________________________________

  SUBROUTINE TRIM_SAVE(JSTART,JSTOP)

    ! TRIM_SAVE saves the discarded solution when the history queue
    ! is trimmed. UNUM655=655 and UNUM656=656 specify the units to
    ! which te data are written. The information is written to files
    ! with names 'dropped_tqueue.dat' and 'dropped_queue.dat'
    ! The information can be retrieved from these files with TRIM_GET.

    ! .. Scalar Arguments ..
    INTEGER :: JSTART, JSTOP
    ! INTEGER I,J,K
    ! ..
    ! .. Local Scalars ..
    INTEGER :: IDROPM1

    OPEN (UNIT=UNUM655,FILE='dropped_tqueue.dat',FORM='UNFORMATTED', &
         STATUS='OLD',POSITION='REWIND')
    OPEN (UNIT=UNUM656,FILE='dropped_queue.dat',FORM='UNFORMATTED', &
         STATUS='OLD',POSITION='REWIND')

    WRITE (UNUM655) TQLAST

    ! Dropping IDROPM1 points:
    IDROPM1 = JSTOP - JSTART + 1
    WRITE (UNUM655) IDROPM1

    ! Dropped queue information to be dropped:
    WRITE (UNUM655) TQUEUE(JSTART:JSTOP)
    WRITE (UNUM656) QUEUE(1:MYN,10*(JSTART-1)+1:10*JSTOP)
    ! DO I = JSTART, JSTOP
    !    WRITE(UNUM655) TQUEUE(I)
    ! END DO
    ! DO I = JSTART, JSTOP
    !    ISTART = 10*(I-1) + 1
    !    DO J = 1, MYN
    !       DO K = ISTART, ISTART+9
    !          WRITE(UNUM656) QUEUE(J,K)
    !       END DO
    !    END DO
    ! END DO

    CLOSE (UNUM655)
    CLOSE (UNUM656)

    RETURN
  END SUBROUTINE TRIM_SAVE
  !____________________________________________________________________________

  SUBROUTINE PREP_OUT(NGUSER,I,JROOT,INDEXG,NROOT,ISTERMINAL,TEVENT)

    ! .. Scalar Arguments ..
    INTEGER :: I, NGUSER, NROOT
    LOGICAL, OPTIONAL :: TEVENT
    ! ..
    ! .. Array Arguments ..
    INTEGER :: INDEXG(NGUSER), JROOT(NGUSER)
    LOGICAL, OPTIONAL :: ISTERMINAL(NGUSER)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC PRESENT
    ! ..
    ! Note:
    ! MYREPORT is set in DDE_GRT5 and DDE_GRT2. If it is, the root
    ! is not reported according to how direction is defined.
    NROOT = 0
    IF ((JROOT(I)==1) .AND. (INDEXG(I)==0) .AND. (MYREPORT(I)==1)) THEN
       NROOT = I
       IF (PRESENT(TEVENT)) TEVENT = ISTERMINAL(I)
       RETURN
    END IF

    RETURN
  END SUBROUTINE PREP_OUT
  !____________________________________________________________________________

  !************END PRIVATE AUXILIARY SUBROUTINES*****************

  !************BEGIN PUBLIC AUXILIARY SUBROUTINES****************

  FUNCTION DDE_SET(RE,RE_VECTOR,AE,AE_VECTOR,HINIT,HMAX,ISTERMINAL, &
       DIRECTION,JUMPS,NEUTRAL,TRACK_DISCONTINUITIES,TRACKING_LEVEL, &
       INTERPOLATION,THIT_EXACTLY,MAX_EVENTS,MAX_STEPS, &
       MOVING_AVERAGE,MAX_DELAY,TRIM_FREQUENCY) RESULT (OPTS)

    ! .. Function Return Value ..
    TYPE (DDE_OPTS) :: OPTS
    ! ..
    ! .. Scalar Arguments ..
    DOUBLE PRECISION, OPTIONAL :: AE, HINIT, HMAX, RE, MAX_DELAY
    INTEGER, OPTIONAL :: MAX_EVENTS, MAX_STEPS, TRACKING_LEVEL, &
         MOVING_AVERAGE, TRIM_FREQUENCY
    LOGICAL, OPTIONAL :: INTERPOLATION, NEUTRAL, TRACK_DISCONTINUITIES
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, OPTIONAL :: AE_VECTOR(:), RE_VECTOR(:), JUMPS(:), &
         THIT_EXACTLY(:)
    INTEGER, OPTIONAL :: DIRECTION(:)
    LOGICAL, OPTIONAL :: ISTERMINAL(:)
    ! ..
    ! .. Local Scalars ..
    INTEGER :: IER, NAE, NRE
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC MINVAL, PRESENT, SIZE
    ! ..
    MYIPOINT(1:LIPOINT) = 0

    ! Temporarily block the use of the MOVING_AVERAGE option.
    IF (PRESENT(MOVING_AVERAGE)) THEN
       PRINT *, ' This option is not yet activated.'
       STOP
    END IF
    ! ..
    ! Determine whether averaging will be used and if so which type.
    DROPZ = .FALSE.

    IF (PRESENT(NEUTRAL)) THEN
       OPTS%NEUTRAL = NEUTRAL
    ELSE
       OPTS%NEUTRAL = .FALSE.
    END IF

    IF (PRESENT(TRACK_DISCONTINUITIES)) THEN
       OPTS%TRACK_DISCONTINUITIES = TRACK_DISCONTINUITIES
    ELSE
       OPTS%TRACK_DISCONTINUITIES = .TRUE.
    END IF

    IF (PRESENT(TRACKING_LEVEL)) THEN
       OPTS%TRACKING_LEVEL = TRACKING_LEVEL
    ELSE
       OPTS%TRACKING_LEVEL = 7
       IF (PRESENT(NEUTRAL)) OPTS%TRACKING_LEVEL = 0
    END IF

    IF (PRESENT(MAX_EVENTS)) THEN
       OPTS%MAX_EVENTS = MAX_EVENTS
    ELSE
       OPTS%MAX_EVENTS = 200
    END IF

    IF (PRESENT(MAX_STEPS)) THEN
       OPTS%MAX_STEPS = MAX_STEPS
    ELSE
       OPTS%MAX_STEPS = 10000
    END IF

    IF (PRESENT(INTERPOLATION)) THEN
       OPTS%INTERPOLATION = INTERPOLATION
    ELSE
       OPTS%INTERPOLATION = .FALSE.
    END IF

    IF (PRESENT(RE_VECTOR)) THEN
       IF (MINVAL(RE_VECTOR)<0D0) THEN
          PRINT *, ' All components of RE_VECTOR must be non-negative.'
          STOP
       END IF
       NRE = SIZE(RE_VECTOR)
    ELSE
       NRE = 1
    END IF
    ALLOCATE (OPTS%RELERR(NRE),STAT=IER)
    CALL CHECK_STAT(IER,60)
    MYIPOINT(25) = 1
    IF (PRESENT(RE_VECTOR)) THEN
       OPTS%RELERR = RE_VECTOR
    ELSE IF (PRESENT(RE)) THEN
       IF (RE<=0D0) THEN
          PRINT *, ' RE must be non-negative.'
          STOP
       END IF
       OPTS%RELERR(1) = RE
    ELSE
       OPTS%RELERR(1) = 1D-3
    END IF

    IF (PRESENT(AE_VECTOR)) THEN
       IF (MINVAL(AE_VECTOR)<0D0) THEN
          PRINT *, ' All components of AE_VECTOR must be non-negative.'
          STOP
       END IF
       NAE = SIZE(AE_VECTOR)
    ELSE
       NAE = 1
    END IF
    ALLOCATE (OPTS%ABSERR(NAE),STAT=IER)
    CALL CHECK_STAT(IER,60)
    MYIPOINT(24) = 1
    IF (PRESENT(AE_VECTOR)) THEN
       OPTS%ABSERR = AE_VECTOR
    ELSE IF (PRESENT(AE)) THEN
       IF (AE<=0D0) THEN
          PRINT *, ' AE must be non-negative.'
          STOP
       END IF
       OPTS%ABSERR(1) = AE
    ELSE
       OPTS%ABSERR(1) = 1D-6
    END IF

    IF (PRESENT(JUMPS)) THEN
       ALLOCATE (OPTS%JUMPS(SIZE(JUMPS)),STAT=IER)
       CALL CHECK_STAT(IER,61)
       MYIPOINT(26) = 1
       OPTS%JUMPS = JUMPS
    ELSE
       NULLIFY (OPTS%JUMPS)
       MYIPOINT(26) = 0
    END IF

    IF (PRESENT(HINIT)) THEN
       OPTS%HINIT = HINIT
    ELSE
       OPTS%HINIT = 0D0
    END IF

    IF (PRESENT(HMAX)) THEN
       OPTS%HMAX = HMAX
    ELSE
       OPTS%HMAX = -1D0
    END IF

    IF (PRESENT(ISTERMINAL)) THEN
       ALLOCATE (OPTS%ISTERMINAL(SIZE(ISTERMINAL)),STAT=IER)
       CALL CHECK_STAT(IER,62)
       MYIPOINT(22) = 1
       OPTS%ISTERMINAL = ISTERMINAL
    ELSE
       NULLIFY (OPTS%ISTERMINAL)
       MYIPOINT(22) = 0
    END IF

    IF (PRESENT(DIRECTION)) THEN
       ALLOCATE (OPTS%DIRECTION(SIZE(DIRECTION)),STAT=IER)
       CALL CHECK_STAT(IER,63)
       MYIPOINT(23) = 1
       OPTS%DIRECTION = DIRECTION
    ELSE
       NULLIFY (OPTS%DIRECTION)
       MYIPOINT(23) = 0
    END IF

    IF (PRESENT(THIT_EXACTLY)) THEN
       ALLOCATE (OPTS%THIT_EXACTLY(SIZE(THIT_EXACTLY)),STAT=IER)
       CALL CHECK_STAT(IER,64)
       MYIPOINT(27) = 1
       OPTS%THIT_EXACTLY = THIT_EXACTLY
    ELSE
       NULLIFY (OPTS%THIT_EXACTLY)
       MYIPOINT(27) = 0
    END IF

    ! Trim queue parameters.
    MY_MAX_DELAY = 0.0D0
    MY_TRIM_FREQUENCY = 0
    MIN_DROP = 0
    IF (PRESENT(MAX_DELAY)) THEN
       IF (MAX_DELAY>0.0D0) THEN
          IF (OPTS%INTERPOLATION) THEN
             PRINT *, ' The INTERPOLATION option may not be used if'
             PRINT *, ' the solution queue is to be trimmed during'
             PRINT *, ' the integration. Stopping.'
             STOP
          END IF
          MY_MAX_DELAY = MAX(MY_MAX_DELAY,MAX_DELAY)
          MY_TRIM_FREQUENCY = 500
          MIN_DROP = 100
          IF (PRESENT(TRIM_FREQUENCY)) THEN
             IF (TRIM_FREQUENCY>0) THEN
                !MY_TRIM_FREQUENCY=MAX(MY_TRIM_FREQUENCY,TRIM_FREQUENCY)
                MY_TRIM_FREQUENCY = TRIM_FREQUENCY
                MIN_DROP = MIN(MIN_DROP,MY_TRIM_FREQUENCY/2)
             END IF
          END IF
       END IF
    END IF

    RETURN
  END FUNCTION DDE_SET
  !____________________________________________________________________________

  SUBROUTINE PRINT_STATS(SOL)

    ! .. Structure Arguments ..
    TYPE (DDE_SOL) :: SOL
    ! ..
    PRINT *, ' '
    PRINT *, ' Integration statistics: '
    PRINT *, ' '
    PRINT *, ' Number of successful steps                  = ', &
         SOL%STATS(1)
    PRINT *, ' Number of failed steps                      = ', &
         SOL%STATS(2)
    PRINT *, ' Number of derivative evaluations            = ', &
         SOL%STATS(3)
    PRINT *, ' Number of failed corrector iterations       = ', &
         SOL%STATS(4)
    PRINT *, ' Total array storage used                    = ', &
         SOL%STATS(5)
    PRINT *, ' Total number of root finding functions used = ', &
         SOL%STATS(6)
    PRINT *, ' '

    RETURN
  END SUBROUTINE PRINT_STATS
  !____________________________________________________________________________

  FUNCTION DDE_VAL(T,SOL,COMPONENTS,DERIVATIVES) RESULT (YINT)

    ! Interpolate the solution and derivative using the solution
    ! structure SOL produced by a call to DDE_SOLVER with
    ! INTERPOLATION = .TRUE.

    !     Usage:
    !     YINT = DDE_VAL(T,SOL,[COMPONENTS],[DERIVATIVES])
    !            [ ] indicates the argument is optional.
    !     On Entry:
    !     T              -  Array of length NT>0 at which interpolation
    !                       of the solution is desired.
    !     SOL            -  Solution structure produced by the a previous
    !                       call to DDE_SOLVER.
    !     COMPONENTS     -  Optional vector specifying which components of
    !                       the solution are to be interpolated. If not
    !                       present, all components will be interpolated.
    !     DERIVATIVES    -  Optional logical flag which indicates whether
    !                       derivatives are to be calculated.  Derivatives
    !                       will be calculated if DERIVATIVES is present
    !                       and .TRUE.
    !     ON RETURN:
    !     YINT%COMPONENTS - Vector to indicate which components were
    !                       interpolated. If COMPONENTS is not present
    !                       in the argument list for DDE_VAL, all
    !                       components will be interpolated and
    !                       YINT%COMPONENTS will be the vector with
    !                       components 1,...,NEQN. If COMPONENTS is
    !                       present in the argument list,
    !                       YINT%COMPONENTS = COMPONENTS.
    !     YINT%TVALS      - Copy of the T array.
    !     YINT%YT         - Interpolated solution. Array with dimensions
    !                       SIZE(T) by SIZE(COMPONENTS).
    !                       YINT%YT(I,J) is an approximation to solution
    !                       component J at time T(I).
    !     YINT%DT         - If DERIVATIVES = .TRUE., the interpolated
    !                       derivatives array with dimensions
    !                       SIZE(T) by SIZE(COMPONENTS).
    !                       YINT%YT(I,J) is an approximation to derivative
    !                       component J at time T(I).

    ! .. Function Return Value ..
    TYPE (DDE_INT), TARGET :: YINT
    ! ..
    ! .. Structure Arguments ..
    !TYPE (DDE_SOL), TARGET :: SOL
    TYPE (DDE_SOL) :: SOL
    ! ..
    ! .. Scalar Arguments ..
    LOGICAL, OPTIONAL :: DERIVATIVES
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION :: T(:)
    INTEGER, OPTIONAL :: COMPONENTS(:)
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: C, DELINT, TN, TO, TVAL
    INTEGER :: I, IBEGIN, IDIR, IER, INDEXO, INEW, IOLD, IQUEUE, ISTART, &
         J, JQUEUE, N, NC, NT
    LOGICAL :: POLYD, POLYS
    ! ..
    ! .. Local Arrays ..
    DOUBLE PRECISION :: W(0:9), WD(0:9)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC PRESENT, SIZE
    ! ..
    ! Check dimensions.

    N = SIZE(SOL%QUEUE,1)

    NT = SIZE(T)
    IF (NT<1) THEN
       PRINT *, ' The length of T must be positive in DDE_VAL.'
       STOP
    END IF
    ! Process the components array, YINT%COMPONENTS.
    IF (PRESENT(COMPONENTS)) THEN
       NC = SIZE(COMPONENTS)
       IF (NC<1 .OR. NC>N) THEN
          PRINT *, ' The size of COMPONENTS must be between 1 and NEQN.'
          GOTO 30
       ELSE
          DO I = 1, NC
             IF (COMPONENTS(I)<1 .OR. COMPONENTS(I)>N) THEN
                PRINT *, ' COMPONENTS(I) must be between 1 and NEQN'
                GOTO 30
             END IF
          END DO
       END IF
       ALLOCATE (YINT%COMPONENTS(NC),STAT=IER)
       CALL CHECK_STAT(IER,65)
       !SOL%IPOINT(21) = 1
       YINT%COMPONENTS = COMPONENTS
    ELSE
       NC = N
       ALLOCATE (YINT%COMPONENTS(NC),STAT=IER)
       CALL CHECK_STAT(IER,66)
       !SOL%IPOINT(21) = 1
       DO I = 1, NC
          YINT%COMPONENTS(I) = I
       END DO
    END IF
    ! Set the local calculate derivatives flag.
    IF (PRESENT(DERIVATIVES)) THEN
       POLYD = DERIVATIVES
    ELSE
       POLYD = .FALSE.
    END IF
    ! Process T array, YINT%TVALS.
    ALLOCATE (YINT%TVALS(NT),STAT=IER)
    CALL CHECK_STAT(IER,67)
    !SOL%IPOINT(18) = 1
    YINT%TVALS(1:NT) = T(1:NT)
    ! Process solution array, YINT%YT.
    ALLOCATE (YINT%YT(NT,NC),STAT=IER)
    CALL CHECK_STAT(IER,68)
    !SOL%IPOINT(19) = 1
    POLYS = .TRUE.
    ! Process derivative array, YINT%DT.
    IF (POLYD) THEN
       ALLOCATE (YINT%DT(NT,NC),STAT=IER)
       CALL CHECK_STAT(IER,69)
       !SOL%IPOINT(20) = 1
    ELSE
       ALLOCATE (YINT%DT(0,0),STAT=IER)
       CALL CHECK_STAT(IER,69)
       !SOL%IPOINT(20) = 1
    END IF

    ! IQUEUE points to the most recent addition to the queue.
    ! JQUEUE points to the oldest addition to the queue.
    IQUEUE = SOL%IPOINT(2)
    JQUEUE = SOL%IPOINT(3)

    ! Check that an integration step has been completed.
    IF (IQUEUE==0) THEN
       PRINT *, ' The first integration step has not been completed.'
       PRINT *, ' There is insufficient information in the solution'
       PRINT *, ' queue to perform the interpolation requested in'
       PRINT *, ' DDE_VAL.'
       GOTO 30
    END IF
    ! We land here if T is spanned by the solution queue. Interpolate
    ! using the solution history queue.

    ! Set some method flags.
    IDIR = SOL%IPOINT(5)
    IPAIR = 1
    IORDER = 6
    MYMETHOD = 4
    ! The interpolation loop begins here.
    ISTART = 0
    DO 20 I = 1, NT

       TVAL = T(I)

       ! Check if TVAL falls before the initial time.
       IF ((IDIR==1 .AND. TVAL<SOL%TQUEUE(0)) .OR. &
          (IDIR==-1 .AND. TVAL>SOL%TQUEUE(0))) THEN
          PRINT *, ' T falls before the initial time.'
          PRINT *, ' There is insufficient information in the solution'
          PRINT *, ' queue to perform the interpolation requested in'
          PRINT *, ' DDE_VAL. Use your history function to obtain the'
          PRINT *, ' interpolated solution requested in DDE_VAL.'
          GOTO 30
       END IF

       ! Check if T falls beyond the last point in the queue.
       IF ((IDIR==1 .AND. TVAL>SOL%TQUEUE(IQUEUE)) .OR. &
          (IDIR==-1 .AND. TVAL<SOL%TQUEUE(IQUEUE))) THEN
          PRINT *, ' The requested value of T is not spanned by the'
          PRINT *, ' solution queue. It is not possible to perform'
          PRINT *, ' interpolation requested in DDE_VAL.'
          GOTO 30
       END IF

       ! Bracket TVAL in the queue.

       IOLD = JQUEUE
       INEW = IQUEUE
       IF (IOLD/=1) THEN
          PRINT *, ' IOLD is not 1 in DDE_VAL.'
       END IF

       ! The bracketing search begins here.
       ! Linear search; remember last bracketing interval.

       IF (IDIR==1) THEN
          ! SOL%TQUEUE is increasing.
          IF (TVAL<=SOL%TQUEUE(ISTART)) THEN
             DO J = ISTART, 0, -1
                IF (TVAL>=SOL%TQUEUE(J)) THEN
                   ! SOL%TQUEUE(J)<=TVAL<=SOL%TQUEUE(J+1)
                   ISTART = J
                   INDEXO = J + 1
                   TO = SOL%TQUEUE(J)
                   TN = SOL%TQUEUE(INDEXO)
                   GOTO 10
                END IF
             END DO
          ELSE
             ! TVAL>SOL%TQUEUE(ISTART)
             DO J = ISTART + 1, INEW
                IF (TVAL<=SOL%TQUEUE(J)) THEN
                   ! SOL%TQUEUE(J-1)<=TVAL<=SOL%TQUEUE(J)
                   ISTART = J - 1
                   INDEXO = J
                   TO = SOL%TQUEUE(ISTART)
                   TN = SOL%TQUEUE(INDEXO)
                   GOTO 10
                END IF
             END DO
          END IF
       ELSE
          ! SOL%TQUEUE is decreasing.
          IF (TVAL>=SOL%TQUEUE(ISTART)) THEN
             DO J = ISTART, 0, -1
                IF (TVAL<=SOL%TQUEUE(J)) THEN
                   ! SOL%TQUEUE(J)>=TVAL>=SOL%TQUEUE(J+1)
                   ISTART = J
                   INDEXO = J + 1
                   TO = SOL%TQUEUE(J)
                   TN = SOL%TQUEUE(INDEXO)
                   GOTO 10
                END IF
             END DO
          ELSE
             ! TVAL<SOL%TQUEUE(ISTART)
             DO J = ISTART + 1, INEW
                IF (TVAL>=SOL%TQUEUE(J)) THEN
                   ! SOL%TQUEUE(J-1)>=TVAL>=SOL%TQUEUE(J)
                   ISTART = J - 1
                   INDEXO = J
                   TO = SOL%TQUEUE(ISTART)
                   TN = SOL%TQUEUE(INDEXO)
                   GOTO 10
                END IF
             END DO
          END IF
       END IF

       ! The bracketing search ends here.
10     CONTINUE
       IF (((IDIR==1) .AND. (TVAL<TO .OR. TVAL>TN)) .OR. &
          ((IDIR==-1) .AND. (TVAL>TO .OR. TVAL<TN))) THEN
          PRINT *, ' A search error occurred in DDE_VAL.'
          STOP
       END IF

       ! The data to be interpolated corresponds to TO = TOLD, TN = TNEW.
       ! INDEXO is the pointer for the corresponding slug of data. The
       ! data for TOLD to be interpolated begins in queue location
       ! (1,IBEGIN+1).
       IBEGIN = (INDEXO-1)*10
       ! Define the step size corresponding to the interpolation data.
       DELINT = TN - TO

       ! Determine the value of C at which the solution or the
       ! derivative is to be interpolated.
       IF (ABS(DELINT)<=0.0D0) THEN
          C = 0.0D0
       ELSE
          C = (TVAL-TO)/DELINT
       END IF
       ! Calculate the interpolation coefficients.
       CALL DDE_WSET(C,W,WD,POLYS,POLYD,MYMETHOD)

       ! Interpolate the solution.
       SOL%YOFT(1:N) = SOL%QUEUE(1:N,IBEGIN+1) + DELINT * &
           (W(0)*SOL%QUEUE(1:N,IBEGIN+2)+W(2)*SOL%QUEUE(1:N,IBEGIN+3)+ &
            W(3)*SOL%QUEUE(1:N,IBEGIN+4)+W(4)*SOL%QUEUE(1:N,IBEGIN+5)+ &
            W(5)*SOL%QUEUE(1:N,IBEGIN+6)+W(6)*SOL%QUEUE(1:N,IBEGIN+7)+ &
            W(7)*SOL%QUEUE(1:N,IBEGIN+8)+W(8)*SOL%QUEUE(1:N,IBEGIN+9)+ &
            W(9)*SOL%QUEUE(1:N,IBEGIN+10))
       ! Load the solution.
       YINT%YT(I,:) = SOL%YOFT(YINT%COMPONENTS)
       IF (POLYD) THEN
          ! Interpolate the derivative.
          SOL%YOFT(1:N) = WD(0)*SOL%QUEUE(1:N,IBEGIN+2) + &
               WD(2)*SOL%QUEUE(1:N,IBEGIN+3) + WD(3)*SOL%QUEUE(1:N,IBEGIN+4) + &
               WD(4)*SOL%QUEUE(1:N,IBEGIN+5) + WD(5)*SOL%QUEUE(1:N,IBEGIN+6) + &
               WD(6)*SOL%QUEUE(1:N,IBEGIN+7) + WD(7)*SOL%QUEUE(1:N,IBEGIN+8) + &
               WD(8)*SOL%QUEUE(1:N,IBEGIN+9) + WD(9)*SOL%QUEUE(1:N,IBEGIN+10)
          ! Load the derivative.
          YINT%DT(I,:) = SOL%YOFT(YINT%COMPONENTS)
       END IF

       ! The interpolation loop ends here.
20  END DO
30  CONTINUE

    RETURN
  END FUNCTION DDE_VAL
  !____________________________________________________________________________
  SUBROUTINE DDE_DRV1(N,NLAG,OPTIONS,DERIVS,GUSER,CHANGE,OUT_FCN,BETA, &
       YINIT,TSPAN,NGUSER,TRIM_GET)

    ! This first level driver must be called by DKL_1, DKL_2, DKL_3, or DKL_4.
    ! PRIVATE global arrays are used throughout DDE_SOLVER.

    !                        Major Subroutines

    !  The major tasks performed by DDE_DRV1 include the following.
    !     1. First level driver called by the user.
    !     2. Allocate the initial storage and partition the work arrays.
    !     3. Call DDE_DRV2 to perform the integration.
    !     4. Allocate additional storage when requested by DDE_DRV2.
    !     5. Deallocate the work arrays upon completion of the integration.
    !     6. If the INTERPOLATION option iflag is .TRUE., function
    !        DDE_VAL may be used following the call to DDE_SOLVER to
    !        interpolate the solution and derivative.

    !  The major tasks performed by DDE_DRV2 include the following.
    !     1. Second level driver called by DDE_DRV1.
    !     2. Set the integration parameters and flags.
    !     3. Calculate the initial step size.
    !     4. Set up and manage the root finding related arrays.
    !     5. Request additional storage from DDE_DRV1 as necessary.
    !     6. Call DDE_DRV3 to perform the step integrations.
    !     7. Interpolate the solution and derivative at output points.
    !     8. Call CHANGE to allow problem changes when events are located.

    !  The major tasks performed by DDE_DRV3 include the following.
    !     1. Third level driver called by DDE_DRV2.
    !     2. Manage the extrapolatory root finding before each step.
    !     3. Call DDE_STEP to take the step.
    !     4. Call DDE_ERR1 to estimate the local error.
    !     5. Call DDE_HST3 to calculate step sizes.
    !     6. Manage the interpolatory root finding following each step.
    !     7. Force event times to be included as integration mesh points.
    !     8. Call DDE_QUE1 to update the solution history queue.

    !  The major tasks performed by DDE_STEP include the following.
    !     1. Core step integrator called by DDE_DRV3.
    !     2. Perform the corrector iterations and apply the Runge-Kutta
    !        formulas.
    !     3. Call DDE_ZSET to set the Z array (and DZ array for neutral
    !        problems) for use in the user derivative subroutine when
    !        called by DDE_STEP.

    !  The major tasks performed by DDE_ZSET include the following.
    !     1. Core interpolation subroutine called by DDE_STEP.
    !     2. Calculate the Z array of delayed solution values.
    !     3. Calculate the DZ array of delayed derivative values for
    !        neutral problems.

    !                       Auxiliary Subroutines

    !     DDE_VAL:   Interpolate the solution (user callable).
    !     DDE_HST1:  Calculate the initial step size.
    !     DDE_HST2:  Interpolation subroutine called by DDE_HST1.
    !     DDE_HST3:  Calculate step sizes for DDE_DRV3.
    !     DDE_HST4:  Limit the step size to be no larger than HMAX.
    !     DDE_QUE1:  Update the solution history queue.
    !     DDE_SRC2:  Perform searches of the solution history queue.
    !     DDE_SRC1:  Utility subroutine called by DDE_SRC2.
    !     DDE_GRT2:  Rootfinder called by DDE_GRT5.
    !     DDE_GRT3:  Evaluate the event residuals for DDE_GRT5.
    !     DDE_GRT4:  Manage the root finding array storage.
    !     DDE_GRT5:  Direct the interpolatory root finding for DDE_DRV3.
    !     DDE_GRT6:  Set up the initial event rootfinding functions.
    !     DDE_ERR1:  Calculate the local error estimate.
    !     DDE_POLY:  Core subroutine to evaluate the solution and
    !                derivative polynomials.
    !     DDE_WSET:  Define the polynomial coefficients.
    !     DDE_DERV:  Wrapper subroutine to evaluate the user defined
    !                derivatives (DERIVS).
    !     DDE_ZSET2: Utility subroutine for DDE_ZSET.
    !     DDE_ZSET3: Utility subroutine for DDE_ZSET (and DDE_HST2).
    !     DDE_INTP:  Interpolate the solution and derivative for DDE_DRV2
    !                at output points for use in OUT_FCN.
    !     DDE_BETA:  Wrapper subroutine to evaluate the user defined
    !                delay times (BETA).
    !     DDE_YINIT: Wrapper subroutine to evaluate the user defined
    !                history function (YINIT).
    !     DDE_TREE:  Construct the discontinuity tree for constant
    !                delay problems.
    !     DDE_SORT:  Sort arrays for DDE_TREE.

    !     SIZES OF Q, W, and INTEGER WORK ARRAYS

    INTEGER STEPS
    ! ..
    ! .. Structure Arguments ..
    TYPE (DDE_OPTS), TARGET :: OPTIONS
    ! ..
    ! .. Scalar Arguments ..
    INTEGER :: IFLAG, N, NGUSER, NLAG
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION :: TSPAN(:)
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL BETA, CHANGE, DERIVS, GUSER, OUT_FCN, YINIT, TRIM_GET
    INTERFACE
       SUBROUTINE BETA(T,Y,BVAL)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y
          DOUBLE PRECISION, DIMENSION(:) :: BVAL
          INTENT(IN):: T,Y
          INTENT(OUT) :: BVAL
       END SUBROUTINE BETA
    END INTERFACE
    INTERFACE
       SUBROUTINE CHANGE(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
            DIRECTION,ISTERMINAL,QUIT)
         INTEGER :: NEVENT
         INTEGER, DIMENSION(:) :: DIRECTION
         DOUBLE PRECISION :: TEVENT,HINIT
         DOUBLE PRECISION, DIMENSION(:) :: YEVENT,DYEVENT
         LOGICAL :: QUIT
         LOGICAL, DIMENSION(:) :: ISTERMINAL
         INTENT(IN) :: NEVENT,TEVENT
         INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,ISTERMINAL,QUIT
       END SUBROUTINE CHANGE
    END INTERFACE
    INTERFACE
       SUBROUTINE DERIVS(T,Y,Z,DY)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y,DY
          DOUBLE PRECISION, DIMENSION(:,:) :: Z
          INTENT(IN):: T,Y,Z
          INTENT(OUT) :: DY
       END SUBROUTINE DERIVS
    END INTERFACE
    INTERFACE
       SUBROUTINE GUSER(T,Y,DYDT,Z,G)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         DOUBLE PRECISION, DIMENSION(:) :: G
         INTENT(IN):: T,Y,DYDT,Z
         INTENT(OUT) :: G
       END SUBROUTINE GUSER
    END INTERFACE
    INTERFACE
       SUBROUTINE OUT_FCN(T,Y,DY,N,NEVENT)
         INTEGER :: N,NEVENT
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
       END SUBROUTINE OUT_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE YINIT(T,Y)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y
         INTENT(IN):: T
         INTENT(OUT) :: Y
      END SUBROUTINE YINIT
    END INTERFACE
    INTERFACE
       SUBROUTINE TRIM_GET
       END SUBROUTINE TRIM_GET
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    INTEGER :: I, IER, LG, LG_NEW, LIW_NEW, LQUEUE_NEW, LTQUEUE_NEW, QCOLS, &
         QCOLS_NEW
    ! ..
    ! .. Local Arrays ..
    DOUBLE PRECISION, ALLOCATABLE :: G_TEMP(:), Q_TEMP(:,:), TQUEUE_TEMP(:)
    INTEGER, ALLOCATABLE :: ITEMP(:), JTEMP(:), LTEMP(:)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ASSOCIATED, SIZE
    ! ..

    ! Define the integration status flag.
    IFLAG = 0

    DEBUG = .FALSE.
    ! DEBUG = .TRUE. ! Produces tons of output that is usually of no interest

    MYNEUTRAL = OPTIONS%NEUTRAL

    MYNGUSER = NGUSER
    IF (MYNGUSER>0) THEN
       ALLOCATE (MYISTERMINAL(MYNGUSER),MYDIRECTION(MYNGUSER), &
            MYREPORT(MYNGUSER),STAT=IER)
       CALL CHECK_STAT(IER,70)
       IF (ASSOCIATED(OPTIONS%ISTERMINAL)) THEN
          MYISTERMINAL(1:MYNGUSER) = OPTIONS%ISTERMINAL(1:MYNGUSER)
       ELSE
          MYISTERMINAL(1:MYNGUSER) = .FALSE.
       END IF
       IF (ASSOCIATED(OPTIONS%DIRECTION)) THEN
          MYDIRECTION(1:MYNGUSER) = OPTIONS%DIRECTION(1:MYNGUSER)
       ELSE
          MYDIRECTION(1:MYNGUSER) = 0           
       END IF
    END IF

    IF (ASSOCIATED(OPTIONS%JUMPS)) THEN
       MYNJUMPS = SIZE(OPTIONS%JUMPS)
       ALLOCATE (MYJUMPS(MYNJUMPS),STAT=IER)
       CALL CHECK_STAT(IER,71)
       MYJUMPS(1:MYNJUMPS) = OPTIONS%JUMPS(1:MYNJUMPS)
    ELSE
       MYNJUMPS = 0
    END IF
    MYMAX_EVENTS = OPTIONS%MAX_EVENTS
    MYMAX_STEPS = OPTIONS%MAX_STEPS
    IF (ASSOCIATED(OPTIONS%THIT_EXACTLY)) THEN
       MYNTHIT = SIZE(OPTIONS%THIT_EXACTLY)
       ALLOCATE (MYTHIT(MYNTHIT),STAT=IER)
       CALL CHECK_STAT(IER,72)
       MYTHIT(1:MYNTHIT) = OPTIONS%THIT_EXACTLY(1:MYNTHIT)
    ELSE
       MYNTHIT = 0
    END IF
    NLAGS = NLAG
    MYN = N
    UNUM655 = 655
    UNUM656 = 656

    IF (HAVE_TRIM_GET) THEN
       OPEN(UNIT=UNUM655,FILE='dropped_tqueue.dat',FORM="UNFORMATTED",&
            STATUS="REPLACE",POSITION="REWIND")
       OPEN(UNIT=UNUM656,FILE='dropped_queue.dat',FORM="UNFORMATTED",&
            STATUS="REPLACE",POSITION="REWIND")
    END IF

    ! Allocate the work arrays.

    ! History arrays, ZARRAY and DARRAY.
    IF (NLAGS>0) THEN
       ALLOCATE (ZARRAY(N,NLAGS),DARRAY(N,NLAGS),ZANDD(N,2*NLAGS), &
            STARTS(NLAGS),MYBVAL(NLAGS),DARRAY2(NEQN_USER,NLAGS_USER), &
            ZANDD2(NEQN_USER,2*NLAGS_USER),ZARRAY2(NEQN_USER,NLAGS_USER), &
            STAT=IER)
       CALL CHECK_STAT(IER,74)
    END IF

    ! Solution queue arrays Q(N,QCOLS), TQUEUE(0:LTQUEUE).
    ! The queue starts with enough space to hold this many steps and
    ! increments by this amount as necessary.
    STEPS = MY_TRIM_FREQUENCY + 100
    LTQUEUE = STEPS
    QCOLS = 10*STEPS
    MYQCOLS = QCOLS
    LQUEUE = N*QCOLS

    ! Root finding arrays, MYJROOT(LIW), INDEXG(LIW), LEVEL(LIW),
    ! GTOLD, GTNEW, GTROOT, GA, GB, GTROOT2, TG.
    LIW = 200 + NLAGS
    LG = LIW
    ALLOCATE (QUEUE(N,QCOLS),TQUEUE(0:LTQUEUE),MYJROOT(LIW),INDEXG(LIW), &
         LEVEL(LIW),VTEMP3(2*N),KMAT(N,20),GTOLD(LG),GTNEW(LG),GTROOT(LG), &
         GA(LG),GB(LG),GTROOT2(LG),TG(LG),STAT=IER)
    CALL CHECK_STAT(IER,75)
    ! For array bounds checking when reallocate and copy to new.
    QUEUE(1:N,1:QCOLS) = 0D0
    TQUEUE(0:LTQUEUE) = 0D0
    MYJROOT(1:LIW) = 0
    INDEXG(1:LIW) = 0
    LEVEL(1:LIW) = 0
    VTEMP3(1:2*N) = 0D0
    KMAT(1:N,1:20)= 0D0
    GTOLD(1:LG) = 0D0
    GTNEW(1:LG) = 0D0
    GTROOT(1:LG) = 0D0
    GA(1:LG) = 0D0
    GB(1:LG) = 0D0
    GTROOT2(1:LG) = 0D0
    TG(1:LG) = 0D0

    ! Work arrays.
    ALLOCATE (MYYOLD(N),MYDYOLD(N),MYYNEW(N),MYDYNEW(N),MYYROOT(N), &
         MYDYROOT(N),MYYOUT(N),MYDYOUT(N),MYABSER(N),MYRELER(N),MYR(N), &
         MYYEXTP(N),MYDYEXTP(N),Y0SAVE(NEQN_USER),F0SAVE(NEQN_USER), &
         PVARD(16), PVARI(27), PVARL(5),STAT=IER)
    CALL CHECK_STAT(IER,76)
    PVARD(1:16) = 0D0
    PVARI(1:27) = 0
    PVARL(1:5) = .FALSE.

    ! Define the integration status flag.
    ! IFLAG = 0

    ! Load the initial solution.
    ! CALL DDE_YINIT(TSPAN(1),YINIT)
    ! YOLD(1:NEQN_USER) = VTEMP3(1:NEQN_USER)
    ! MYYOLD(1:NEQN_USER) = YOLD(1:NEQN_USER)
    ! Y0SAVE(1:NEQN_USER) = YOLD(1:NEQN_USER)
    CALL DDE_YINIT(TSPAN(1),YINIT)
    MYYOLD(1:NEQN_USER) = VTEMP3(1:NEQN_USER)
    Y0SAVE(1:NEQN_USER) = MYYOLD(1:NEQN_USER)

    ! Load YOLD into the first slot of the queue for use until
    ! the first integration step has been completed.
    QUEUE(1:N,1) = MYYOLD(1:N)

10  CONTINUE

    IF (DEBUG) THEN
       WRITE(DBGUNIT, 9999)
       9999 FORMAT(' Call level 2 driver DDE_DRV2 from DDE_DRV1.')
    END IF

    CALL DDE_DRV2(DERIVS,GUSER,CHANGE,OUT_FCN,BETA,YINIT,TSPAN,NSPAN, &
         IFLAG,TRIM_GET,OPTIONS)

    ! Check if DDE_DRV2 is asking for more storage for the work
    ! arrays and return to the calling program if it is not.

    IF (IFLAG<10) THEN

    IF (DEBUG) THEN
       WRITE(DBGUNIT, 9996)
       WRITE(DBGUNIT, 9997)
       WRITE(DBGUNIT, 9998)
       9996 FORMAT(' ')
       9997 FORMAT(' In DDE_DRV1 deallocate work arrays and return ')
       9998 FORMAT(' to calling programc.')
    END IF

       ! Deallocate the work arrays and return to the calling program.

       ! Approximate total length of the arrays used.
       ARRAY_STORAGE = 10*MYIPOINT(6) + LQUEUE + LTQUEUE + 35*N + &
            (7*N+3)*NLAGS + LENTREE + 48
       ! Total number of root functions used.
       ROOT_FUNCTIONS = MYIPOINT(6)
       ! Note: Some arrays are deallocated in DKL_*.
       DEALLOCATE (MYYOLD,MYDYOLD,MYYNEW,MYDYNEW,MYYROOT,MYDYROOT,MYYOUT, &
            MYDYOUT,MYABSER,MYRELER,MYR,MYYEXTP,MYDYEXTP,Y0SAVE,F0SAVE, &
            PVARD,PVARI,PVARL,STAT=IER)
       CALL CHECK_STAT(IER,77)
       DEALLOCATE (MYJROOT,INDEXG,LEVEL,VTEMP3,KMAT,GTOLD,GTNEW,GTROOT,GA, &
            GB,GTROOT2,TG,STAT=IER)
       CALL CHECK_STAT(IER,78)
       IF (NLAGS>0) THEN
          DEALLOCATE (ZARRAY,DARRAY,ZANDD,STARTS,MYBVAL,DARRAY2,ZARRAY2, &
               ZANDD2,STAT=IER)
          CALL CHECK_STAT(IER,79)
       END IF
       IF (MYNGUSER>0) THEN
          DEALLOCATE (MYISTERMINAL,MYDIRECTION,MYREPORT,STAT=IER)
          CALL CHECK_STAT(IER,80)
       END IF
       IF (ASSOCIATED(OPTIONS%JUMPS)) THEN
          DEALLOCATE (MYJUMPS,STAT=IER)
          CALL CHECK_STAT(IER,81)
       END IF
       IF (ASSOCIATED(OPTIONS%THIT_EXACTLY)) THEN
          DEALLOCATE (MYTHIT,STAT=IER)
          CALL CHECK_STAT(IER,82)
       END IF
       IF (LENTREE>0) THEN
          DEALLOCATE (DTREE,STAT=IER)
          CALL CHECK_STAT(IER,83)
       END IF
       ! Note: We need to add the interpolation arrays to SOL
       ! in DKL_* before deallocating them.
       ! DEALLOCATE(QUEUE,TQUEUE,STAT=IER)
       ! CALL CHECK_STAT(IER,84)
    ELSE
       ! Or, allocate the requested storage if it is.
       ! IFLAG = 10 means more storage for solution history queue.
       ! IFLAG = 11 means more storage for root finding.
       IF (IFLAG==10) THEN
          QCOLS_NEW = QCOLS + 10*STEPS
          LTQUEUE_NEW = LTQUEUE + STEPS
          ! LQUEUE_NEW = LQUEUE + N*QCOLS_NEW
          LQUEUE_NEW = N*QCOLS_NEW

          IF (DEBUG) THEN
             WRITE(DBGUNIT,9993)
             WRITE(DBGUNIT,9994)
             WRITE(DBGUNIT,9995) LQUEUE_NEW
             9993 FORMAT(' ')
             9994 FORMAT(' In DDE_DRV1. Allocate more storage for the')
             9995 FORMAT(' solution queue: LTQUEUE_NEW = ', I0)
          END IF

          ALLOCATE (Q_TEMP(N,QCOLS_NEW),TQUEUE_TEMP(0:LTQUEUE_NEW),STAT=IER)
          CALL CHECK_STAT(IER,85)
          ! For array bounds checking when reallocate and copy to new.
          Q_TEMP(1:N,1:QCOLS_NEW) = 0D0
          TQUEUE_TEMP(0:LTQUEUE_NEW) = 0D0
          Q_TEMP(1:N,1:QCOLS) = QUEUE(1:N,1:QCOLS)
          TQUEUE_TEMP(0:LTQUEUE) = TQUEUE(0:LTQUEUE)
          DEALLOCATE (QUEUE,TQUEUE,STAT=IER)
          CALL CHECK_STAT(IER,86)
          ALLOCATE (QUEUE(N,QCOLS_NEW),TQUEUE(0:LTQUEUE_NEW),STAT=IER)
          CALL CHECK_STAT(IER,87)
          ! For array bounds checking when reallocate and copy to new.
          QUEUE(1:N,1:QCOLS_NEW) = 0D0
          TQUEUE(0:LTQUEUE_NEW) = 0D0 
          QUEUE(1:N,1:QCOLS) = Q_TEMP(1:N,1:QCOLS)
          TQUEUE(0:LTQUEUE) = TQUEUE_TEMP(0:LTQUEUE)
          DEALLOCATE (Q_TEMP,TQUEUE_TEMP,STAT=IER)
          CALL CHECK_STAT(IER,88)
          ! Update the dimensions.
          LQUEUE = LQUEUE_NEW
          QCOLS = QCOLS_NEW
          MYQCOLS = QCOLS
          LTQUEUE = LTQUEUE_NEW

          ! End of IF IFLAG = 10.
       END IF

       IF (IFLAG==11) THEN               

          ! The size of the MYJROOT, INDEXG, and LEVEL ARRAYS will
          ! be increased by ADD_TO_IW WORDS. Potentially, every one
          ! of the LIW root functions could have a root at the same
          ! point.

          LIW_NEW = LIW + ADD_TO_IW

          IF (DEBUG) THEN
             WRITE(DBGUNIT, 9993)
             WRITE(DBGUNIT, 9991)
             WRITE(DBGUNIT, 9992) LIW_NEW
             9991 FORMAT(' IFLAG = 11. Allocate root finding arrays')
             9992 FORMAT(' in DDE_DRV1: LIW_NEW = ', I0)
          END IF

          ALLOCATE (JTEMP(LIW_NEW),ITEMP(LIW_NEW),LTEMP(LIW_NEW),STAT=IER)
          CALL CHECK_STAT(IER,89)
          ! For array bounds checking when reallocate and copy to new.
          JTEMP(1:LIW_NEW) = 0
          ITEMP(1:LIW_NEW) = 0
          LTEMP(1:LIW_NEW) = 0
          JTEMP(1:LIW) = MYJROOT(1:LIW)
          ITEMP(1:LIW) = INDEXG(1:LIW)
          LTEMP(1:LIW) = LEVEL(1:LIW)
          DEALLOCATE (MYJROOT,INDEXG,LEVEL,STAT=IER)
          CALL CHECK_STAT(IER,90)
          ALLOCATE (MYJROOT(LIW_NEW),INDEXG(LIW_NEW),LEVEL(LIW_NEW), &
               STAT=IER)
          CALL CHECK_STAT(IER,91)
          ! For array bounds checking when reallocate and copy to new.
          MYJROOT(1:LIW_NEW) = 0
          INDEXG(1:LIW_NEW) = 0
          LEVEL(1:LIW_NEW) = 0
          MYJROOT(1:LIW) = JTEMP(1:LIW)
          INDEXG(1:LIW) = ITEMP(1:LIW)
          LEVEL(1:LIW) = LTEMP(1:LIW)
          DEALLOCATE (JTEMP,ITEMP,LTEMP,STAT=IER)
          CALL CHECK_STAT(IER,92)
          ! Update the dimensions.
          LIW = LIW_NEW
          ! The sizes of the residual arrays will be increased
          ! by ADD_TO_IW WORDS.
          LG_NEW = LG + ADD_TO_IW         

          ALLOCATE (G_TEMP(LG_NEW),STAT=IER)
          CALL CHECK_STAT(IER,93)
          ! For array bounds checking when reallocate and copy to new.
          G_TEMP(1:LG_NEW) = 0D0
          G_TEMP(1:LG) = GTOLD(1:LG)
          DEALLOCATE (GTOLD,STAT=IER)
          CALL CHECK_STAT(IER,94)
          ALLOCATE (GTOLD(LG_NEW),STAT=IER)
          CALL CHECK_STAT(IER,95)
          ! For array bounds checking when reallocate and copy to new.
          GTOLD(1:LG_NEW) = 0D0
          GTOLD(1:LG) = G_TEMP(1:LG)
          G_TEMP(1:LG) = GTNEW(1:LG)
          DEALLOCATE (GTNEW,STAT=IER)
          CALL CHECK_STAT(IER,96)
          ALLOCATE (GTNEW(LG_NEW),STAT=IER)
          CALL CHECK_STAT(IER,97)
          ! For array bounds checking when reallocate and copy to new.
          GTNEW(1:LG_NEW) = 0D0
          GTNEW(1:LG) = G_TEMP(1:LG)
          G_TEMP(1:LG) = GTROOT(1:LG)
          DEALLOCATE (GTROOT,STAT=IER)
          CALL CHECK_STAT(IER,98)
          ALLOCATE (GTROOT(LG_NEW),STAT=IER)
          CALL CHECK_STAT(IER,99)
          ! For array bounds checking when reallocate and copy to new.
          GTROOT(1:LG_NEW) = 0D0
          GTROOT(1:LG) = G_TEMP(1:LG)
          G_TEMP(1:LG) = GA(1:LG)
          DEALLOCATE (GA,STAT=IER)
          CALL CHECK_STAT(IER,100)
          ALLOCATE (GA(LG_NEW),STAT=IER)
          CALL CHECK_STAT(IER,101)
          ! For array bounds checking when reallocate and copy to new.
          GA(1:LG_NEW) = 0D0
          GA(1:LG) = G_TEMP(1:LG)
          G_TEMP(1:LG) = GB(1:LG)
          DEALLOCATE (GB,STAT=IER)
          CALL CHECK_STAT(IER,102)
          ALLOCATE (GB(LG_NEW),STAT=IER)
          CALL CHECK_STAT(IER,103)
          ! For array bounds checking when reallocate and copy to new.
          GB(1:LG_NEW) = 0
          GB(1:LG) = G_TEMP(1:LG)
          G_TEMP(1:LG) = GTROOT2(1:LG)
          DEALLOCATE (GTROOT2,STAT=IER)
          CALL CHECK_STAT(IER,104)
          ALLOCATE (GTROOT2(LG_NEW),STAT=IER)
          CALL CHECK_STAT(IER,105)
          ! For array bounds checking when reallocate and copy to new.
          GTROOT2(1:LG_NEW) = 0D0
          GTROOT2(1:LG) = G_TEMP(1:LG)
          G_TEMP(1:LG) = TG(1:LG)
          DEALLOCATE (TG,STAT=IER)
          CALL CHECK_STAT(IER,106)
          ALLOCATE (TG(LG_NEW),STAT=IER)
          CALL CHECK_STAT(IER,107)
          ! For array bounds checking when reallocate and copy to new.
          TG(1:LG_NEW) = 0D0
          TG(1:LG) = G_TEMP(1:LG)
          DEALLOCATE (G_TEMP,STAT=IER)
          CALL CHECK_STAT(IER,108)

          ! Update the dimensions.
          LG = LG_NEW

          IF (DEBUG) THEN
             WRITE(DBGUNIT, 9990) LG_NEW
          END IF
          9990 FORMAT(' IFLAG = 11. LG_NEW = ', I10)

          ! End of IF IFLAG = 11.
       END IF

       GOTO 10

    END IF

    ERROR_FLAG = IFLAG

    RETURN
  END SUBROUTINE DDE_DRV1
  !____________________________________________________________________________
  SUBROUTINE DDE_DRV2(DERIVS,GUSER,CHANGE,OUT_FCN,BETA,YINIT,TSPAN,NSPAN, &
       IFLAG,TRIM_GET,OPTIONS)

    ! This second level driver must be called by DDE_DRV1.

    ! .. Structure Arguments ..
    TYPE (DDE_OPTS), TARGET :: OPTIONS
    ! ..
    ! .. Scalar Arguments ..
    INTEGER, INTENT (INOUT) :: IFLAG
    INTEGER :: NSPAN
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION :: TSPAN(:)
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL BETA, CHANGE, DERIVS, GUSER, OUT_FCN, YINIT, TRIM_GET
    INTERFACE
       SUBROUTINE BETA(T,Y,BVAL)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y
          DOUBLE PRECISION, DIMENSION(:) :: BVAL
          INTENT(IN):: T,Y
          INTENT(OUT) :: BVAL
       END SUBROUTINE BETA
    END INTERFACE
    INTERFACE
       SUBROUTINE CHANGE(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
            DIRECTION,ISTERMINAL,QUIT)
         INTEGER :: NEVENT
         INTEGER, DIMENSION(:) :: DIRECTION
         DOUBLE PRECISION :: TEVENT,HINIT
         DOUBLE PRECISION, DIMENSION(:) :: YEVENT,DYEVENT
         LOGICAL :: QUIT
         LOGICAL, DIMENSION(:) :: ISTERMINAL
         INTENT(IN) :: NEVENT,TEVENT
         INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,ISTERMINAL,QUIT
       END SUBROUTINE CHANGE
    END INTERFACE
    INTERFACE
       SUBROUTINE DERIVS(T,Y,Z,DY)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y,DY
          DOUBLE PRECISION, DIMENSION(:,:) :: Z
          INTENT(IN):: T,Y,Z
          INTENT(OUT) :: DY
       END SUBROUTINE DERIVS
    END INTERFACE
    INTERFACE
       SUBROUTINE GUSER(T,Y,DYDT,Z,G)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         DOUBLE PRECISION, DIMENSION(:) :: G
         INTENT(IN):: T,Y,DYDT,Z
         INTENT(OUT) :: G
       END SUBROUTINE GUSER
    END INTERFACE
    INTERFACE
       SUBROUTINE OUT_FCN(T,Y,DY,N,NEVENT)
         INTEGER :: N,NEVENT
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
       END SUBROUTINE OUT_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE YINIT(T,Y)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y
         INTENT(IN):: T
         INTENT(OUT) :: Y
      END SUBROUTINE YINIT
    END INTERFACE
    INTERFACE
       SUBROUTINE TRIM_GET
       END SUBROUTINE TRIM_GET
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    !DOUBLE PRECISION, SAVE :: BIG, HINIT, HINITM, HLEFT, HMAX, HNEXT, &
    !     LAGMIN, TEMP, TFINAL, TINC, TNEW, TOLD, TOUT, TROOT, TROOT2, TVAL
    !INTEGER, SAVE :: ANCESL, DFLAG, IA, I, IDIR, IER, IHFLAG, INDEX, INEED, &
    !     IQUEUE, ISPAN, ITREE, J, JNDEX, JQUEUE, KROOT, LEVMAX, LEXTRA, &
    !     MORDER, N, NFOUND, NG, NGUSER, NGUSP1, NROOT, NUMPT, OFFRT
    !LOGICAL, SAVE :: AUTORF, ITERATE, QUIT, TEVENT, USETREE
    DOUBLE PRECISION :: BIG, HINIT, HINITM, HLEFT, HMAX, HNEXT, LAGMIN, &
         TEMP, TFINAL, TINC, TNEW, TOLD, TOUT, TROOT, TROOT2, TVAL, &
         HBEFORE
    INTEGER :: ANCESL, DFLAG, IA, I, IDIR, IER, IHFLAG, INDEX, INEED, &
         IQUEUE, ISPAN, ITREE, J, JNDEX, JQUEUE, KROOT, LEVMAX, LEXTRA, &
         MORDER, N, NFOUND, NG, NGUSER, NGUSP1, NROOT, NUMPT, OFFRT, ITASK
    LOGICAL :: AUTORF, ITERATE, QUIT, TEVENT, USETREE
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN, SIGN, SQRT

    ! ..
    ! If this is a call following an IFLAG = 10 or 11 return, restore
    ! the values of the local variables.

    IF (IFLAG==10 .OR. IFLAG==11) THEN
       ITASK = 1
       CALL DDE_PVAR(BIG, HINIT, HINITM, HLEFT, HMAX, HNEXT, LAGMIN, &
            TEMP, TFINAL, TINC, TNEW, TOLD, TOUT, TROOT, TROOT2, TVAL, &
            ANCESL, DFLAG, IA, I, IDIR, IER, IHFLAG, INDEX, INEED, &
            IQUEUE, ISPAN, ITREE, J, JNDEX, JQUEUE, KROOT, LEVMAX, LEXTRA, &
            MORDER, N, NFOUND, NG, NGUSER, NGUSP1, NROOT, NUMPT, OFFRT, &
            AUTORF, ITERATE, QUIT, TEVENT, USETREE, ITASK)
    END IF
    ! ..
    ! If this is the first initialize the persistant local variables
    ! to avoid problems assigning a variable to an undefined variable
    ! when DDE_PAR is called (the values used are set in DDE_DRV1 and
    ! are irrelevant).
    IF (IFLAG==0) THEN
       ITASK = 1            
       CALL DDE_PVAR(BIG, HINIT, HINITM, HLEFT, HMAX, HNEXT, LAGMIN, &
            TEMP, TFINAL, TINC, TNEW, TOLD, TOUT, TROOT, TROOT2, TVAL, &
            ANCESL, DFLAG, IA, I, IDIR, IER, IHFLAG, INDEX, INEED, &
            IQUEUE, ISPAN, ITREE, J, JNDEX, JQUEUE, KROOT, LEVMAX, LEXTRA, &
            MORDER, N, NFOUND, NG, NGUSER, NGUSP1, NROOT, NUMPT, OFFRT, &
            AUTORF, ITERATE, QUIT, TEVENT, USETREE, ITASK)

       N = MYN
       NGUSER = MYNGUSER

       ! Define parameters and flags used throughout the code.

       ! Initialize the derivative counter.
       NFEVAL = 0
       ! Initialize the residual evaluation counter.
       NGEVAL = 0
       ! Initialize the safety counter for NEGVAL.
       NGEMAX = 500
       ! Initialize the successful steps counter.
       NSTEPS = 0
       ! Initialize the failed steps counter.
       NFAILS = 0
       ! Initialize the failed corrector iteration step counter.
       NFAILC = 0
       ! Maximum number of corrections.
       MAXCOR = 5
       ! Convergence factor for the predictor corrector iteration.
       KCFAC = 1.0D0
       ! Maximum number of step halvings at a discontinuity.
       ! (safety flag).
       MAXTRY = 100
       ! Maximum factor by which the step size will be decreased.
       ! FACDN = 10.0D0
       FACDN = 2.0D0
       ! Maximum factor by which the step size will be increased.
       ! FACUP = 5.0D0
       FACUP = 2.0D0
       ! Attempt to obtain a local error of efac times the error
       ! tolerance vector.
       EFAC = 0.5D0
       ! Decrease the next estimated step size by a factor of
       ! TOLFAC**(1/6).
       TOLFAC = 15.0D0
       ! Roundoff factors for termination criteria.
       U10 = 10.0D0*UROUND
       U13 = 13.0D0*UROUND
       U26 = 26.0D0*UROUND
       U65 = 65.0D0*UROUND
       ! Floor value RER for the relative error tolerances.
       REMIN = 1.0D-12
       RER = 2.0D0*UROUND + REMIN
       ! Convergence parameters for the root finder.
       NSIG = 14
       GRTOL = 0.0D0
       IF (NSIG>0) GRTOL = 10.0D0**(-NSIG)
       GRTOL = MAX(GRTOL,U10)
       ! Array storage.
       ARRAY_STORAGE = 0
       ! Root finding functions.
       ROOT_FUNCTIONS = 0
       ! Use the sixth order method for the basic integration.
       IORDER = 6
       ! Use the (5,6) embedded pair.
       IPAIR = 1
       MYMETHOD = 4
       ! If IRTMAX = 2 and a root TROOT is loacted in (TOLD,TROOT),
       ! second pass root finding will be done for the interval
       ! (TROOT,TNEW).
       IRTMAX = 2
       ! Trim solution queue information.
       INFO_DROPPED = .FALSE.

       ! Initialize the problem.
       TOLD = TSPAN(1)
       TNEW = TOLD
       TFINAL = TSPAN(NSPAN)
       HINIT = OPTIONS%HINIT
       HMAX = OPTIONS%HMAX

       ! Define the absolute error tolerances. The tolerances are
       ! relaxed by a factor of 100 for all but the most averaged
       ! components.
       IA = SIZE(OPTIONS%ABSERR)
       IF (IA==1) THEN
          MYABSER(1:N) = OPTIONS%ABSERR(1)
       ELSE
          MYABSER(1:N) = OPTIONS%ABSERR(1:N)
       END IF

       ! Define the relative error tolerances. The tolerances
       ! are relaxed to the smaller of 10*RELERR and 0.05
       ! for all but the most averaged components.
       IA = SIZE(OPTIONS%RELERR)
       IF (IA==1) THEN
          MYRELER(1:N) = OPTIONS%RELERR(1)
       ELSE
          MYRELER(1:N) = OPTIONS%RELERR(1:N)
       END IF

       AUTORF = OPTIONS%TRACK_DISCONTINUITIES
       IF (AUTORF .AND. MY_MAX_DELAY>0.0D0) THEN
          ! Do not perform discontinuity tracking if the
          ! trim queue  option is used.
          AUTORF = .FALSE.
       END IF
       IF (NLAGS<1) NLAGS = 0
       IF (NLAGS==0) AUTORF = .FALSE.
       IF (HMAX<0D0) THEN
          HMAX = 0.1D0*(TFINAL-TOLD)
       END IF
       MYNEUTRAL = OPTIONS%NEUTRAL
       ISPAN = 1
       TINC = TSPAN(2) - TSPAN(1)
       ! TINC will be reset following the call to D6HST1.
       IF (OPTIONS%INTERPOLATION) TINC = TSPAN(NSPAN) - TSPAN(ISPAN)
       MYTINIT = TOLD
       TQUEUE(0) = MYTINIT
       TQLAST = MYTINIT
       IQUEUE = 0
       MYIPOINT(2) = IQUEUE
       JQUEUE = 1
       MYIPOINT(3) = JQUEUE
       IFIRST = 0
       IF (NLAGS>0) THEN
          STARTS(1:NLAGS) = 0
       END IF
       LEVMAX = OPTIONS%TRACKING_LEVEL
       LEVMAX = MAX(LEVMAX,0)
       ! Do an extra level of tracking if TINIT is one of the
       ! jump points since this indicates the solution itself
       ! may have a jump at TINIT.
       LEXTRA = 0
       IF (LEVMAX==7) THEN
          IF (MYNJUMPS>0) THEN
             DO I = 1, MYNJUMPS
                ! IF (MYJUMPS(I)==MYTINIT) LEXTRA = 1
                IF (ABS(MYJUMPS(I)-MYTINIT)<=0.0D0) LEXTRA = 1
             END DO
          END IF
          LEVMAX = LEVMAX + LEXTRA
       END IF

       ! Check for some common errors. Also reset the relative
       ! error tolerances to the floor value if necessary.
       IF (N<1) THEN
          PRINT *, ' N is less than one.'
          GOTO 70
       END IF
       IF (ABS(TFINAL-TOLD)<=0.0D0) THEN
          PRINT *, ' TFINAL and TOLD must be different.'
          GOTO 70
       END IF
       DO I = 1, N
          IF (ABS(MYRELER(I))>0.0D0 .AND. MYRELER(I)<RER) THEN
             PRINT *, ' You have specified a nonzero relative error'
             PRINT *, ' which is smaller than the floor value. The'
             PRINT *, ' value has been reset and execution will'
             PRINT *, ' continue.'
             MYRELER(I) = RER
          END IF
          IF (MYABSER(I)<0.0D0) THEN
             PRINT *, ' ABSER is negative.'
             GOTO 70
          END IF
          IF (MYRELER(I)<0.0D0) THEN
             PRINT *, ' RELER is negative.'
             GOTO 70
          END IF
          IF (ABS(MYABSER(I))<=0.0D0 .AND. ABS(MYRELER(I))<=0.0D0) THEN
             PRINT *, ' Both RELER and ABSER are zero.'
             GOTO 70
          END IF
       END DO
       IF (ABS(TINC)<=0.0D0) THEN
          PRINT *, ' The output print increment TINC must be nonzero.'
          GOTO 70
       END IF

       ! Calculate the initial derivatives.
       IF (NLAGS>0) THEN
          CALL DDE_BETA(TOLD,MYYOLD,BETA)
          DO I = 1, NLAGS
             TVAL = MYBVAL(I)
             IF (ABS(TVAL-MYTINIT)<=0.0D0) THEN
                ZARRAY(1:N,I) = MYYOLD(1:N)
             ELSE
                CALL DDE_YINIT(TVAL,YINIT)
                ZARRAY(1:N,I) = VTEMP3(1:N)
             END IF
             IF (MYNEUTRAL) THEN
                DARRAY(1:N,I) = VTEMP3(N+1:2*N)
             END IF
          END DO
       END IF
       CALL DDE_DERV(TOLD,MYYOLD,MYDYOLD,DERIVS)

       ! Determine the direction of the integration.
       ! IDIR =  1 means integration to the right.
       ! IDIR = -1 means integration to the left.
       IDIR = 1
       IF (TFINAL<TOLD) IDIR = -1
       MYIPOINT(5) = IDIR

       ! Check that the exact hit times are in (TINIT,TFINAL).
       IF (MYNTHIT>0) THEN
          IDIR = 1
          IF (TFINAL<TOLD) IDIR = -1
          MYIPOINT(5) = IDIR
          IF (IDIR==1) THEN
             DO I = 1, MYNTHIT
                IF (MYTHIT(I)<=TOLD .OR. MYTHIT(I)>=TFINAL) THEN
                   PRINT *, ' Specified exact hit times must be'
                   PRINT *, ' strictly between the initial and'
                   PRINT *, ' final integration times.'
                   GOTO 70
                END IF
             END DO
          ELSE
             DO I = 1, MYNTHIT
                IF (MYTHIT(I)>=TOLD .OR. MYTHIT(I)<=TFINAL) THEN
                   PRINT *, ' Specified exact hit times must be'
                   PRINT *, ' strictly between the initial and'
                   PRINT *, ' final integration times.'
                   GOTO 70
                END IF
             END DO
          END IF
       END IF

       ! If the delays are constant, build the discontinuity
       ! tree directly rather than use root finding.
       LENTREE = 0
       ITREE = 0
       USETREE = .FALSE.
       IF (CONSTANT_DELAYS) THEN
          IF (NLAGS>0) THEN
             CALL DDE_TREE(MYTINIT,TFINAL,LEVMAX,BETA,MYYOLD,DFLAG)
             IF (DFLAG/=0) THEN
                PRINT *, ' An error occurred in DDE_TREE.'
                IFLAG = DFLAG
                RETURN
             END IF
             AUTORF = .FALSE.
             USETREE = .TRUE.
             ITREE = 2
             ! Recalculate the initial delays.
             IF (NLAGS>0) CALL DDE_BETA(TOLD,MYYOLD,BETA)
          END IF
       ELSE
          ! If the delays are not constant and exact hit points
          ! have been specified, the tree will consist of the
          ! exact hit times.
          IF (MYNTHIT>0) THEN
             USETREE = .TRUE.
             LENTREE = MYNTHIT + 2
             ITREE = LENTREE
             ! Allocation of DTREE.
             ALLOCATE (DTREE(LENTREE),STAT=IER)
             CALL CHECK_STAT(IER,109)
             ! Sort the exact hit times.
             CALL DDE_SORT(MYTHIT,MYTHIT,MYNTHIT,IDIR,IFLAG)
             IF (IFLAG/=0) RETURN
             DTREE(1) = TOLD
             DTREE(2:ITREE-1) = MYTHIT(1:MYNTHIT)
             DTREE(ITREE) = TFINAL
          END IF
       END IF

       ! Set up the root finding functions.
       JN = 0
       JP = 0
       CALL DDE_GRT6(AUTORF,NG)
       IF (NG>0 .AND. LIW<NG) THEN
          PRINT *, ' The length of the integer work array is too small.'
          GOTO 70
       END IF

       ! Calculate the initial delays.
       IF (NLAGS>0) CALL DDE_BETA(TOLD,MYYOLD,BETA)

       ! Give TINC the correct sign.
       TINC = SIGN(ABS(TINC),TFINAL-TOLD)

       ! Finish the initialization.
       TOUT = TOLD + TINC
       INDEX = 1

       ! IDIR =  1 means integration to the right.
       ! IDIR = -1 means integration to the left.
       IDIR = 1
       IF (TFINAL<TOLD) IDIR = -1
       MYIPOINT(5) = IDIR

       ! Output the initial solution.
       NROOT = 0
       CALL OUT_FCN(TOLD,MYYOLD,MYDYOLD,N,NROOT)

       ! Turn off the root finding until the second step is completed.
       IFIRST = 0

       ! Calculate the initial step size if the user did not supply it.
       IF (ABS(HINIT)<=0.0D0) THEN
          DO I = 1, N
             KMAT(I,1) = 0.5D0*(MYABSER(I)+MYRELER(I)*ABS(MYYOLD(I)))
             IF (ABS(KMAT(I,1))<=0.0D0) THEN
                PRINT *, ' A component of the solution vanished making it'
                PRINT *, ' impossible to continue with a pure relative'
                PRINT *, ' error test.'
                GOTO 70
             END IF
          END DO
          MORDER = 6
          BIG = SQRT(REALMAX)
          CALL DDE_HST1(DERIVS,TOLD,TOUT,MYYOLD,MYDYOLD,KMAT(:,1),MORDER, &
               BIG,KMAT(:,2),KMAT(:,3),KMAT(:,4),KMAT(:,5),BETA,YINIT,IHFLAG, &
               HINIT)
          IF (IHFLAG/=0) GOTO 70
          ! Recalculate the initial delays.
          IF (NLAGS>0) CALL DDE_BETA(TOLD,MYYOLD,BETA)
       END IF

       ! Limit the initial step size to be no larger than HMAX, no
       ! larger than the distance to the first output point, no
       ! no larger than 0.1D0 time the smallest nonzero delay, and
       ! no larger than the distance to the first point in the
       ! tree. But use only the first and last entries of TSPAN if
       ! the interpolation option is set.
       IF (OPTIONS%INTERPOLATION) THEN
          LAGMIN = ABS(TSPAN(NSPAN)-TSPAN(1))
       ELSE
          LAGMIN = MIN(ABS(HMAX),ABS(TSPAN(2)-TSPAN(1)))
       END IF
       IF (NLAGS>0) THEN
          DO I = 1, NLAGS
             TEMP = ABS(MYTINIT-MYBVAL(I))
             IF (TEMP>U65) LAGMIN = MIN(LAGMIN,TEMP)
          END DO
          LAGMIN = 0.1D0*LAGMIN
          HINITM = MIN(ABS(HINIT),ABS(HMAX),LAGMIN)
          ! If the discontinuity tree is being used, limit the step
          ! size to not step past the first discontinuity time.
          IF (USETREE) HINITM = MIN(HINITM,ABS(DTREE(2)-MYTINIT))
          HINIT = SIGN(HINITM,TFINAL-TOLD)
          HMAX = SIGN(ABS(HMAX),TFINAL-TOLD)
          HNEXT = HINIT
       END IF

       ! If this is an ode, limit the step size in a manner
       ! consistent with the above.
       IF (NLAGS==0) THEN
          LAGMIN = 0.1D0*LAGMIN
          HINITM = MIN(ABS(HINIT),ABS(HMAX),LAGMIN)
          HINIT = SIGN(HINITM,TFINAL-TOLD)
          HMAX = SIGN(ABS(HMAX),TFINAL-TOLD)
          HNEXT = HINIT
       END IF

       ! Check if the initial step size is too small.
       TEMP = U13*ABS(TOLD)
       IF (ABS(HINIT)<=TEMP) THEN
          PRINT *, ' The required step size is too small.'
          GOTO 70
       END IF

       ! Limit the initial step size to be no larger than
       ! ABS(TFINAL-TOLD) and no larger than HMAX.
       HINITM = MIN(ABS(HINIT),ABS(TFINAL-TOLD))
       HINIT = SIGN(HINITM,TFINAL-TOLD)
       HMAX = SIGN(ABS(HMAX),TFINAL-TOLD)
       HNEXT = HINIT
       IF (DEBUG) THEN
          WRITE(DBGUNIT,9999) HNEXT, HINIT
          9999 FORMAT(' Spot 1 (START). HNEXT, HINIT = ', 2E20.10)
       END IF              

       ! ISPAN = 1 and TINC = TSPAN(2)-TSPAN(1) at this point.
       ! If the interpolation option is set, update ISPAN and
       ! TINC if necessary.
       IF (OPTIONS%INTERPOLATION) THEN
          DO I = 1, NSPAN - 1
             IF (IDIR==1) THEN
                IF (TOLD+HINIT>TSPAN(ISPAN+1)) THEN
                   ISPAN = ISPAN + 1
                ELSE
                   GOTO 10
                END IF
             ELSE
                IF (TOLD+HINIT<TSPAN(ISPAN+1)) THEN
                   ISPAN = ISPAN + 1
                ELSE
                   GOTO 10
                END IF
             END IF
          END DO
10        CONTINUE
          IF (ISPAN<NSPAN) TINC = TSPAN(ISPAN+1) - TSPAN(ISPAN)
          TOUT = TOLD + TINC
       END IF

       ! End of IFLAG = 0 initialization.
    END IF

    ! Call the integrator to advance the integration one step.
20  CONTINUE

    ! NOTE:
    ! A call will land here if an IFLAG>9 return was made to
    ! DDE_DRV1 to request more storage for the Q, W, or IW arrays.

    ! Check if the size of the queue needs to be increased.
    NUMPT = LQUEUE/(10*N)
    IQUEUE = MYIPOINT(2)
    IF (IQUEUE>=NUMPT-1) THEN
       IFLAG = 10
       ! First save the local variables.
       ITASK = 2
       CALL DDE_PVAR(BIG, HINIT, HINITM, HLEFT, HMAX, HNEXT, LAGMIN, &
            TEMP, TFINAL, TINC, TNEW, TOLD, TOUT, TROOT, TROOT2, TVAL, &
            ANCESL, DFLAG, IA, I, IDIR, IER, IHFLAG, INDEX, INEED, &
            IQUEUE, ISPAN, ITREE, J, JNDEX, JQUEUE, KROOT, LEVMAX, LEXTRA, &
            MORDER, N, NFOUND, NG, NGUSER, NGUSP1, NROOT, NUMPT, OFFRT, &
            AUTORF, ITERATE, QUIT, TEVENT, USETREE, ITASK)
       ! Get more storage.
       RETURN
    END IF

    ! Check if the sizes of the root finding arrays need to be
    ! increased.
    IF (LIW<NG+ADD_TO_IW/5) THEN

       IF (DEBUG) THEN
          WRITE(DBGUNIT,9998)
          WRITE(DBGUNIT,9997)
          9998 FORMAT(' Return to DDE_RV2 to DDE_DRV2 to allocate more')
          9997 FORMAT(' storage for the root finding arrays.')
       END IF

       IFLAG = 11
       ! First save the local variables.
       ITASK = 2
       CALL DDE_PVAR(BIG, HINIT, HINITM, HLEFT, HMAX, HNEXT, LAGMIN, &
            TEMP, TFINAL, TINC, TNEW, TOLD, TOUT, TROOT, TROOT2, TVAL, &
            ANCESL, DFLAG, IA, I, IDIR, IER, IHFLAG, INDEX, INEED, &
            IQUEUE, ISPAN, ITREE, J, JNDEX, JQUEUE, KROOT, LEVMAX, LEXTRA, &
            MORDER, N, NFOUND, NG, NGUSER, NGUSP1, NROOT, NUMPT, OFFRT, &
            AUTORF, ITERATE, QUIT, TEVENT, USETREE, ITASK)
       ! Get more storage.
       RETURN
    END IF

    IFLAG = 0

    ! Limit the step size not to step past the next discontinuity
    ! in the tree.
    IF (USETREE) THEN
       HLEFT = DTREE(ITREE) - TOLD
       IF (HLEFT>U26*MAX(ABS(TNEW),ABS(TFINAL))) THEN
          HBEFORE = HNEXT
          HNEXT = SIGN(MIN(ABS(HNEXT),ABS(HLEFT)),HLEFT)
          IF (DEBUG) THEN
             WRITE(DBGUNIT,9996) HNEXT, HBEFORE
             9996 FORMAT(' Spot 2 (TREE). HNEXT, HBEFORE = ', 2E20.10)
          END IF
       END IF
    END IF

    ! Continue the integration.

    IF (DEBUG) THEN
       WRITE(DBGUNIT,9994)
       WRITE(DBGUNIT,9995) TOLD, HNEXT
       9995 FORMAT(' Perform new step with DDE_DRV3. T, H =', 2D15.5)
       9994 FORMAT(' ')
    END IF

    NGEVAL = 0
    CALL DDE_DRV3(TOLD,TNEW,TROOT,TFINAL,DERIVS,NG,GUSER,HNEXT,HMAX,BETA, &
         YINIT,INDEX,TROOT2,TRIM_GET)

    IF (INDEX==-14) GOTO 50
    IF (INDEX==2) THEN
       TROOT = TNEW
       TROOT2 = TNEW
       MYYROOT(1:N) = MYYNEW(1:N)
       MYDYROOT(1:N) = MYDYNEW(1:N)
    END IF

    ! Terminate the integration if an abnormal return was made from
    ! the integrator.
    IF (INDEX<2) GOTO 70

    ! Check if this return was due to a root at the initial point
    ! and return if it is.
    IF (INDEX==5) THEN
       GOTO 20
    END IF

    ! Turn on the root finding.
    IFIRST = 1

    ! Have we reached a point in the discontinuity tree?
    IF (USETREE .AND. ITREE<LENTREE) THEN

       IF (ABS(TNEW-DTREE(ITREE))<U26*MAX(ABS(TNEW),ABS(TFINAL))) THEN
          ITREE = ITREE + 1

       IF (DEBUG) THEN
          WRITE(DBGUNIT,9993) TNEW
          9993 FORMAT(' We have reached a point in the discontinuity', &
               ' tree. TNEW = ', D20.10)
       END IF

       END IF
    END IF

    ! Is this a root return?

    IF (INDEX==3 .OR. INDEX==4) THEN

       ! Yes, this is a root return.
30     CONTINUE

       IF (DEBUG) THEN
          WRITE(DBGUNIT,9992) TROOT
          9992 FORMAT(' Root return from DDE_DRV3. TROOT = ', D20.10)
       END IF

       ! Are we in one step mode?
       IF (NSPAN==2) THEN
          ! Call PREP_OUT and OUT_FCN in a loop to account
          ! for the possibility that more than one of the
          ! event functions may have a zero at TROOT.
          DO I = 1, NGUSER
             CALL PREP_OUT(NGUSER,I,MYJROOT,INDEXG,NROOT)
             IF (NROOT/=0) THEN
                CALL OUT_FCN(TNEW,MYYNEW,MYDYNEW,N,NROOT)
             END IF
          END DO
       END IF

       ! Have we passed an output point?
       IDIR = MYIPOINT(5)
       IF ((((TROOT>TOUT) .AND. (IDIR==1)) .OR. &
            ((TROOT<TOUT) .AND. (IDIR==-1))) .AND. (NSPAN/=2)) THEN
          ! We have passed the next output point.
          ! TROOT is beyond tout.
          ! Interpolate at TOUT.

            IF (DEBUG) THEN
               WRITE(DBGUNIT,9989)
               WRITE(DBGUNIT,9990) TROOT
               WRITE(DBGUNIT,9991) TOUT
               9989 FORMAT(' We have passed the next output point with')
               9990 FORMAT(' TROOT = ', D20.10, ' beyond TOUT so')
               9991 FORMAT(' interpolate at TOUT = ', D20.10)
            END IF

          CALL DDE_INTP(TOLD,TNEW,TOUT,MYYOUT,MYDYOUT)
          NROOT = 0
          ! Make a print return at TOUT.

          IF (DEBUG) THEN
             WRITE(DBGUNIT,9988) TOUT
             9988 FORMAT(' We have reached an output point TOUT = ', D20.10)
          END IF

          CALL OUT_FCN(TOUT,MYYOUT,MYDYOUT,N,NROOT)
          ! Update the print return output point.
          ISPAN = ISPAN + 1
          IF (ISPAN<NSPAN) TINC = TSPAN(ISPAN+1) - TSPAN(ISPAN)
          TOUT = TOUT + TINC
          ! Loop until the output point TOUT lies beyond the
          ! root, TROOT.
          GOTO 30
       END IF

       ! The root (TROOT) does not lie beyond the next output point.
       ! Make a print return at the root unless we are in one step
       ! mode and the return was made above.

         IF (DEBUG) THEN
            WRITE(DBGUNIT,9986)
            WRITE(DBGUNIT,9987) TROOT, TOUT
            9986 FORMAT(' TROOT does not lie beyond the next output')
            9987 FORMAT(' point. TROOT, TOUT = ', 2D20.10)
         END IF

       IF (NSPAN/=2) THEN
          ! Call PREP_OUT and OUT_FCN in a loop to account
          ! for the possibility that more than one of the
          ! event functions may have a zero at TROOT.
          DO I = 1, NGUSER
             CALL PREP_OUT(NGUSER,I,MYJROOT,INDEXG,NROOT)
             IF (NROOT/=0) THEN
                CALL OUT_FCN(TROOT,MYYROOT,MYDYROOT,N,NROOT)
             END IF
          END DO
       END IF

       ! Restart the integration at the root or continue the
       ! integration from TNEW.

       ! First check if we are close enough to the final
       ! integration time to terminate the integration.
       HLEFT = TFINAL - TROOT
       IF (ABS(HLEFT)<=U26*MAX(ABS(TROOT),ABS(TFINAL))) THEN
          MYYOUT(1:N) = MYYROOT(1:N)
          MYDYOUT(1:N) = MYDYROOT(1:N)
          ! CALL DDE_DERV(TFINAL,MYYOUT,MYDYOUT,DERIVS)
          ! Make a print return for the final integration time.
          IF (NSPAN/=2) THEN
             NROOT = 0
             CALL OUT_FCN(TFINAL,MYYOUT,MYDYOUT,N,NROOT)
          END IF
          ! Adjust the final mesh time in the history queue.
          IQUEUE = MYIPOINT(2)
          TQUEUE(IQUEUE) = TFINAL
          ! Terminate the integration.
          GOTO 60
       END IF

       ! It is not yet time to terminate the integration so prepare
       ! for the restart at the root.

       ! Reset the integration flag.
       INDEX = 1

       ! Turn off the root finding until the second step is completed.
       IFIRST = 0
       OFFRT = 1

       IF (NGUSER>0) THEN
          ! Call PREP_OUT and OUT_FCN in a loop to account for the
          ! possibility that more than one of the event functions
          ! may have a zero at TROOT.
          DO I = 1, NGUSER
             CALL PREP_OUT(NGUSER,I,MYJROOT,INDEXG,NROOT,MYISTERMINAL,TEVENT)
             IF (NROOT>0) THEN
                HINIT = HNEXT
                HINIT = SIGN(MAX(ABS(HINIT),U65),HINIT)
                QUIT = .FALSE.
                IF (DONT_CALL_CHANGE) THEN
                ELSE
                   CALL CHANGE(NROOT,TROOT,MYYROOT,MYDYROOT,HINIT, &
                               MYDIRECTION, MYISTERMINAL,QUIT)
                   IF (QUIT) GOTO 50
                END IF
             END IF
          END DO
       END IF

       ! Check if a terminal event occurred. If so make a print return
       ! and terminate the integration.
       IF (NGUSER>0) THEN
          IF (TEVENT) THEN
             ! Call to OUT_FCN already made above.
             ! CALL OUT_FCN(TROOT,MYYROOT,MYDYROOT,N,NROOT)
             ! Adjust the final time reached in the history queue.
             IQUEUE = MYIPOINT(2)
             TQUEUE(IQUEUE) = TROOT
             GOTO 60
          END IF
       END IF
       HINIT = SIGN(MAX(ABS(HINIT),U65),HINIT)
       HINITM = MIN(ABS(HNEXT),ABS(HINIT),LAGMIN)
       HINIT = SIGN(MAX(HINITM,U65),HINIT)

       ! Accommodate the additional root functions. Determine
       ! how many of the root functions have a zero at TROOT.
       IF (OFFRT/=0) THEN
          KROOT = 0
          ANCESL = -1
          IF (NG>NGUSER) THEN
             ANCESL = 10000
             NGUSP1 = NGUSER + 1
             DO I = NGUSP1, NG
                IF (MYJROOT(I)==1) THEN
                   KROOT = KROOT + 1
                   ! Minimum level from which this root came.
                   IF (LEVMAX>0) ANCESL = MIN(ANCESL,LEVEL(I))
                END IF
             END DO
             ! Reset KROOT:
             IF (KROOT/=0) KROOT = NLAGS*(1+JP+JN)
             IF (LEVMAX>0) THEN
                ! Do not add this root if have reached the maximum
                ! tracking level.
                IF (ANCESL>=LEVMAX) THEN
                   KROOT = 0
                END IF
             END IF
             NFOUND = NFOUND + 1
          END IF

          !Check if there is enough storage to continue.
          !INEED = NG + KROOT
          !INEED = 2*NG + ADD_TO_IW
          !IF (INEED>LIW) THEN
          !   PRINT *, ' The length of the integer work array is too small.'
          !   GOTO 70
          !END IF

          ! Check if the maximum allowable number of events have occurred.
          !IF (MYMAX_EVENTS>0 .AND. NFOUND>=MYMAX_EVENTS) THEN
          IF (MYMAX_EVENTS>0 .AND. NG>=MYMAX_EVENTS) THEN
             !PRINT *, ' The integration will be terminated since the'
             !PRINT *, ' maximum number of events, MAX_EVENTS, has'
             !PRINT *, ' been found. You can increase the default value'
             !PRINT *, ' of 200 by adding MAX_EVENTS=desired number'
             !PRINT *, ' in your options call to DDE_SET. If this error'
             !PRINT *, ' message occurs your problem may hace a'
             !PRINT *, ' vanishing delay or the discontinuity'
             !PRINT *, ' rootfinding may have generated too many'
             !PRINT *, ' too many root functions. You may need to disable'
             !PRINT *, ' discontinuity tracking for this problem by'
             !PRINT *, ' adding TRACK_DISCONTINUITUES=.FALSE. in your'
             !PRINT *, ' options call to DDE_SET.'
             !GOTO 70
             PRINT *, ' Discontinuity root finding will be turned off'
             PRINT *, ' because the maximum allowable number (200) of'
             PRINT *, ' root functions has been reachded.'
             OFFRT = 0
          END IF

          ! Set up the added root finding functions.
          IF (KROOT>0 .AND. OFFRT==1) THEN
             CALL DDE_GRT4(KROOT,TROOT,ANCESL)
             IF (DEBUG) THEN
                WRITE(DBGUNIT,9985) NG, NG+KROOT
                9985 FORMAT('Increase NG from ', I5, ' to ', I5)
             END IF
          END IF

          ! Define the new number of root functions.

          NG = NG + KROOT
          IF (OFFRT==0) NG = NGUSER
          MYIPOINT(6) = NG
       END IF

       ! Re-set the initial time.
       ! TOLD = TROOT

       ! Define the initial step size for the restart.
       HNEXT = HINIT ! Comment out if brave
       IF (DEBUG) THEN
          WRITE(DBGUNIT,9984) HNEXT, HINIT
          9984 FORMAT(' Spot 3 (RESTART). HNEXT, HINIT = ', 2E20.10)
       END IF

       ! If a second root was located in (TROOT,TNEW) limit the
       ! step size so as not to miss it on restarting.
       ! IF (TROOT2 /= TROOT) THEN
       IF (ABS(TROOT2-TROOT)>0.0D0) THEN
          HLEFT = 0.25D0*(TROOT2-TROOT)
          HNEXT = SIGN(MIN(ABS(HNEXT),ABS(HLEFT)),TFINAL-TNEW)
          IF (DEBUG) THEN
             WRITE(DBGUNIT,9983) HNEXT, HINIT
             9983 FORMAT(' Spot 4 (TROOT2). HNEXT, HINIT = ', 2E20.10)
          END IF
       END IF

       ! Use the interpolated solution at the root to update the
       ! left solution.
       ! MYYOLD(1:N) = MYYROOT(1:N)

       ! Update the left derivatives with the system derivatives for
       ! the new left solution (not the interpolated derivatives in
       ! MYDYROOT).
       ! Note: The error vector R will be used as scratch storage.
       IF (NLAGS>0) THEN
          CALL DDE_BETA(TROOT,MYYROOT,BETA)
          CALL DDE_ZSET(TROOT,MYYROOT,HNEXT,YINIT,JNDEX,1,KMAT(:,1), &
               KMAT(:,3),KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8), &
               KMAT(:,9),KMAT(:,10),TOLD,MYYOLD,ITERATE)
          IF (JNDEX/=2) THEN
             INDEX = JNDEX
             GOTO 70
          END IF
       END IF
       CALL DDE_DERV(TROOT,MYYROOT,MYDYROOT,DERIVS)

       ! Re-set the initial time.
       TOLD = TROOT
       TNEW = TROOT

       ! Use the interpolated solution at the root to update the
       ! left solution.
       MYYOLD(1:N) = MYYROOT(1:N)
       MYDYOLD(1:N) = MYDYROOT(1:N)
       MYYNEW(1:N) = MYYROOT(1:N)
       MYDYNEW(1:N) = MYDYROOT(1:N)

       ! Reset the step size if it will throw us beyond the final
       ! integration time on the next step.
       HLEFT = TFINAL - TROOT
       HNEXT = SIGN(MIN(ABS(HNEXT),ABS(HLEFT)),TFINAL-TROOT)

       ! Call the integrator again.
       GOTO 20

    ELSE

       ! No, this is not a root return.
40     CONTINUE

    ! Are we in one step mode?
       IF (NSPAN==2) THEN
          NROOT = 0
          CALL OUT_FCN(TNEW,MYYNEW,MYDYNEW,N,NROOT)
       END IF

       ! Have we reached an output point?
       IDIR = MYIPOINT(5)
       IF ((((TROOT>TOUT) .AND. (IDIR==1)) .OR. &
            ((TROOT<TOUT) .AND. (IDIR==-1))) .AND. (NSPAN/=2)) THEN
          ! We have passed the next output point.
          ! (TNEW is beyond TOUT.)
          ! Interpolate at TOUT.
          CALL DDE_INTP(TOLD,TNEW,TOUT,MYYOUT,MYDYOUT)
          ! Make a print return at TOUT.
          NROOT = 0
          CALL OUT_FCN(TOUT,MYYOUT,MYDYOUT,N,NROOT)
          ! Update the print return output point.
          ISPAN = ISPAN + 1
          IF (ISPAN<NSPAN) TINC = TSPAN(ISPAN+1) - TSPAN(ISPAN)
          TOUT = TOUT + TINC
          ! Loop until the output point TOUT lies beyond the
          ! root, TROOT.
          GOTO 40
       END IF

       ! TNEW does not lie beyond the next output point. Prepare for
       ! another call to the integrator or terminate the integration.

       ! Extrapolate the solution to TFINAL and terminate the
       ! integration if close enough to TFINAL.
       HLEFT = TFINAL - TNEW
       IF (ABS(HLEFT)<=U26*MAX(ABS(TNEW),ABS(TFINAL))) THEN
          MYYOUT(1:N) = MYYNEW(1:N)
          MYDYOUT(1:N) = MYDYNEW(1:N)
          ! CALL DDE_DERV(TFINAL,MYYOUT,MYDYOUT,DERIVS)
          ! Make a print return for the final integration time.
          IF (NSPAN/=2) THEN
             NROOT = 0
             CALL OUT_FCN(TFINAL,MYYOUT,MYDYOUT,N,NROOT)
          END IF
          ! Adjust the final mesh time in the history queue.
          IQUEUE = MYIPOINT(2)
          TQUEUE(IQUEUE) = TFINAL
          ! TERMINATE THE INTEGRATION.
          GOTO 60
       END IF

       ! It is not yet time to terminate the integration. Prepare to
       ! call the integrator again to continue the integration from
       ! TNEW.

       ! Update the independent variable.
       TOLD = TNEW

       ! Update the left solution.
       MYYOLD(1:N) = MYYNEW(1:N)

       ! Update the left derivatives with the system derivatives
       ! for the new left solution.
       MYDYOLD(1:N) = MYDYNEW(1:N)

       ! Reset the step size if it will throw us beyond the final
       ! integration time on the next step.
       HLEFT = TFINAL - TNEW
       HNEXT = SIGN(MIN(ABS(HNEXT),ABS(HLEFT)),TFINAL-TNEW)

       ! Call the integrator again.
       GOTO 20

       ! End of is this a root return block.
    END IF

50  CONTINUE
    ! User termination.
    PRINT *, ' The integration was terminated by the user.'
    IFLAG = 0
    GOTO 80

60  CONTINUE
    ! Normal termination.
    IFLAG = 0
    GOTO 80

70  CONTINUE
    ! Abnormal termination.
    PRINT *, ' An abnormal error was encountered in DDE_DRV2.'
    IFLAG = 1

80  CONTINUE

    RETURN
  END SUBROUTINE DDE_DRV2
  !____________________________________________________________________________

  SUBROUTINE DDE_PVAR(BIG, HINIT, HINITM, HLEFT, HMAX, HNEXT, LAGMIN, &
     TEMP, TFINAL, TINC, TNEW, TOLD, TOUT, TROOT, TROOT2, TVAL, &
     ANCESL, DFLAG, IA, I, IDIR, IER, IHFLAG, INDEX, INEED, &
     IQUEUE, ISPAN, ITREE, J, JNDEX, JQUEUE, KROOT, LEVMAX, LEXTRA, &
     MORDER, N, NFOUND, NG, NGUSER, NGUSP1, NROOT, NUMPT, OFFRT, &
     AUTORF, ITERATE, QUIT, TEVENT, USETREE, ITASK)

  ! The purpose of DDE_PVAR is to save and restore the persistent local
  ! used in subroutine DDE_DRV2. This routine is called only by DDE_DRV2.

  ! Caution: ITASK must have a value of 1 or 2.
      INTEGER ITASK
      DOUBLE PRECISION :: BIG, HINIT, HINITM, HLEFT, HMAX, HNEXT, LAGMIN, &
         TEMP, TFINAL, TINC, TNEW, TOLD, TOUT, TROOT, TROOT2, TVAL
      INTEGER :: ANCESL, DFLAG, IA, I, IDIR, IER, IHFLAG, INDEX, INEED, &
         IQUEUE, ISPAN, ITREE, J, JNDEX, JQUEUE, KROOT, LEVMAX, LEXTRA, &
         MORDER, N, NFOUND, NG, NGUSER, NGUSP1, NROOT, NUMPT, OFFRT
      LOGICAL :: AUTORF, ITERATE, QUIT, TEVENT, USETREE

      IF (ITASK==1) THEN

         ! RESTORE ...

         ! DOUBLE PRECISION 
         BIG     = PVARD(1)
         HINIT   = PVARD(2)
         HINITM  = PVARD(3)
         HLEFT   = PVARD(4)
         HMAX    = PVARD(5)
         HNEXT   = PVARD(6)
         LAGMIN  = PVARD(7)
         TEMP    = PVARD(8)
         TFINAL  = PVARD(9)
         TINC    = PVARD(10)
         TNEW    = PVARD(11)
         TOLD    = PVARD(12)
         TOUT    = PVARD(13)
         TROOT   = PVARD(14)
         TROOT2  = PVARD(15)
         TVAL    = PVARD(16)

         ! INTEGER
         ANCESL  = PVARI(1)
         DFLAG   = PVARI(2)
         IA      = PVARI(3)
         I       = PVARI(4)
         IDIR    = PVARI(5)
         IER     = PVARI(6)
         IHFLAG  = PVARI(7)
         INDEX   = PVARI(8)
         INEED   = PVARI(9)
         IQUEUE  = PVARI(10)
         ISPAN   = PVARI(11)
         ITREE   = PVARI(12)
         J       = PVARI(13)
         JNDEX   = PVARI(14)
         JQUEUE  = PVARI(15)
         KROOT   = PVARI(16)
         LEVMAX  = PVARI(17)
         LEXTRA  = PVARI(18)
         MORDER  = PVARI(19)
         N       = PVARI(20)
         NFOUND  = PVARI(21)
         NG      = PVARI(22)
         NGUSER  = PVARI(23)
         NGUSP1  = PVARI(24)
         NROOT   = PVARI(25)
         NUMPT   = PVARI(26)
         OFFRT   = PVARI(27)

         ! LOGICAL
         AUTORF  = PVARL(1)
         ITERATE = PVARL(2)
         QUIT    = PVARL(3)
         TEVENT  = PVARL(4)
         USETREE = PVARL(5)

      ELSE

         ! SAVE ...

         ! DOUBLE PRECISION 
         PVARD(1)  = BIG
         PVARD(2)  = HINIT
         PVARD(3)  = HINITM
         PVARD(4)  = HLEFT
         PVARD(5)  = HMAX
         PVARD(6)  = HNEXT
         PVARD(7)  = LAGMIN
         PVARD(8)  = TEMP
         PVARD(9)  = TFINAL
         PVARD(10) = TINC
         PVARD(11) = TNEW
         PVARD(12) = TOLD
         PVARD(13) = TOUT
         PVARD(14) = TROOT
         PVARD(15) = TROOT2
         PVARD(16) = TVAL

         ! INTEGER
         PVARI(1)  = ANCESL
         PVARI(2)  = DFLAG
         PVARI(3)  = IA
         PVARI(4)  = I
         PVARI(5)  = IDIR
         PVARI(6)  = IER
         PVARI(7)  = IHFLAG
         PVARI(8)  = INDEX
         PVARI(9)  = INEED
         PVARI(10) = IQUEUE
         PVARI(11) = ISPAN
         PVARI(12) = ITREE
         PVARI(13) = J
         PVARI(14) = JNDEX
         PVARI(15) = JQUEUE
         PVARI(16) = KROOT
         PVARI(17) = LEVMAX
         PVARI(18) = LEXTRA
         PVARI(19) = MORDER
         PVARI(20) = N
         PVARI(21) = NFOUND
         PVARI(22) = NG
         PVARI(23) = NGUSER
         PVARI(24) = NGUSP1
         PVARI(25) = NROOT
         PVARI(26) = NUMPT
         PVARI(27) = OFFRT

         ! LOGICAL
         PVARL(1)  = AUTORF
         PVARL(2)  = ITERATE
         PVARL(3)  = QUIT
         PVARL(4)  = TEVENT
         PVARL(5)  = USETREE

      END IF

    RETURN
  END SUBROUTINE DDE_PVAR
  !____________________________________________________________________________

  SUBROUTINE DDE_DRV3(TOLD,TNEW,TROOT,TFINAL,DERIVS,NG,GUSER,HNEXT,HMAX, &
       BETA,YINIT,INDEX,TROOT2,TRIM_GET)

    ! This third level driver must be called by DDE_DRV2.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: HMAX, TFINAL, TOLD
    DOUBLE PRECISION, INTENT (INOUT) :: HNEXT, TNEW, TROOT
    DOUBLE PRECISION, INTENT (OUT) :: TROOT2
    INTEGER, INTENT (INOUT) :: INDEX
    INTEGER, INTENT (IN) :: NG
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL BETA, DERIVS, GUSER, YINIT, TRIM_GET
    INTERFACE
       SUBROUTINE BETA(T,Y,BVAL)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y
          DOUBLE PRECISION, DIMENSION(:) :: BVAL
          INTENT(IN):: T,Y
          INTENT(OUT) :: BVAL
       END SUBROUTINE BETA
    END INTERFACE
    INTERFACE
       SUBROUTINE DERIVS(T,Y,Z,DY)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y,DY
          DOUBLE PRECISION, DIMENSION(:,:) :: Z
          INTENT(IN):: T,Y,Z
          INTENT(OUT) :: DY
       END SUBROUTINE DERIVS
    END INTERFACE
    INTERFACE
       SUBROUTINE GUSER(T,Y,DYDT,Z,G)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         DOUBLE PRECISION, DIMENSION(:) :: G
         INTENT(IN):: T,Y,DYDT,Z
         INTENT(OUT) :: G
       END SUBROUTINE GUSER
    END INTERFACE
    INTERFACE
       SUBROUTINE YINIT(T,Y)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y
         INTENT(IN):: T
         INTENT(OUT) :: Y
      END SUBROUTINE YINIT
    END INTERFACE
   INTERFACE
      SUBROUTINE TRIM_GET
      END SUBROUTINE TRIM_GET
   END INTERFACE
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: H, HZSET, RATIO, TEMP
    INTEGER :: I, IERTST, IMETH, ITRY, JNDEX, N, NGUSER
    LOGICAL :: CONVERGED, ITERATE
    ! ..
    ! .. Local Arrays ..
    DOUBLE PRECISION :: TVALS(4)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, MIN, SIGN
    ! ..
    N = MYN
    NGUSER = MYNGUSER

    ! Check the termination flag.
    IF (INDEX==0) INDEX = 2
    IF (INDEX<1) INDEX = -7
    IF (INDEX>5) INDEX = -7
    IF (INDEX==-7) THEN
       PRINT *, ' INDEX is less than 1 or greater than 5.'
       GOTO 70
    END IF

    ! Initialization.
    TNEW = TOLD
    TROOT = TOLD
    TROOT2 = TOLD
    RATIO = 0.0D0
    HNEXT = SIGN(MIN(ABS(HNEXT),ABS(HMAX)),TFINAL-TOLD)
    H = HNEXT
    IF (INDEX<0) GOTO 70
    ! If this is the first call or an integration restart.
    IF (INDEX==1) GOTO 10
    GOTO 30

10  CONTINUE

    ! Calculate the initial residuals if necessary.
    IF (NG==0) GOTO 30
    IF (NGUSER>0) THEN
       IF (NLAGS>0) CALL DDE_BETA(TOLD,MYYOLD,BETA)
       IMETH = 1
       ITERATE = .FALSE.
       CALL DDE_ZSET(TOLD,MYYOLD,H,YINIT,INDEX,IMETH,KMAT(:,1),KMAT(:,3), &
            KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9), &
            KMAT(:,10),TOLD,MYYOLD,ITERATE)
       ITERATE = .FALSE.
       IF (INDEX==-12 .OR. INDEX==-4) THEN
          PRINT *, ' An error occurred in the first call to'
          PRINT *, ' DDE_ZSET from DDE_DRV3.'
          GOTO 70
       END IF
    END IF
    CALL DDE_GRT3(TOLD,MYYOLD,MYDYOLD,GUSER,GTOLD,BETA)
    NGEVAL = NGEVAL + 1
    GTROOT(1:NG) = GTOLD(1:NG)

    ! Check for a zero at the left endpoint.
    DO I = 1, NG
       MYJROOT(I) = 0
       IF (ABS(GTOLD(I))<=0.0D0) MYJROOT(I) = 1
    END DO
    DO I = 1, NG
       IF (MYJROOT(I)==1) GOTO 20
    END DO
    GOTO 30
20  CONTINUE
    INDEX = 5
    MYYROOT(1:N) = MYYOLD(1:N)
    GOTO 80
30  CONTINUE

    ! The step integration section begins here.

    ! Integration step from TOLD to TNEW = TOLD + H.

    ! Call DDE_STEP To calculate the new approximation MYYNEW.
    ! Note: The error vector R will be used as scratch storage.

    ITRY = 0
40  CONTINUE
    ITRY = ITRY + 1

    CALL DDE_STEP(TOLD,H,DERIVS,TNEW,INDEX,BETA,YINIT,CONVERGED)

    IF (.NOT. (CONVERGED)) THEN
       IF (ITRY>MAXTRY) THEN
          PRINT *, ' The step size was halved MAXTRY times in DDE_STEP'
          PRINT *, ' due to convergence failures. The integration'
          PRINT *, ' will be terminated.'
          INDEX = -15
          GOTO 70
       END IF
       H = 0.5D0*H
       GOTO 40
    END IF

    ! Check for error return.
    IF (INDEX<0) GOTO 70

    ! Error testing and step size selection section.

    ! Check if it is advisable to continue the integration
    ! with a pure relative error test.
    DO I = 1, N
       IF (MYYNEW(I)*MYYOLD(I)<=0.0D0 .AND. MYRELER(I)<=0.0D0) THEN
          PRINT *, ' A component of the solution vanished making it'
          PRINT *, ' impossible to continue with a pure relative'
          PRINT *, ' error test.'
          INDEX = -8
          GOTO 70
       END IF
    END DO

    ! Calculate the maximum local error estimate RATIO.
    CALL DDE_ERR1(INDEX,RATIO)

    ! Check if the estimated error is acceptable and calculate an
    ! estimated step size with which to redo the step or with which
    ! to take the next step.
    CALL DDE_HST3(IERTST,RATIO,H,HNEXT)
    IF (NLAGS>0) CALL DDE_BETA(TNEW,MYYNEW,BETA)

    ! Limit the step size to be no larger than the maximum step size.
    CALL DDE_HST4(HNEXT,HMAX)

    ! If the estimated error is small enough, continue.
    IF (IERTST==0) GOTO 60

    ! If the step needs to be redone, check if the estimated
    ! step size is too small.

    ! Update the failed step counter.
    NFAILS = NFAILS + 1

    TEMP = MAX(ABS(TOLD),ABS(TNEW))
    TEMP = U13*TEMP
    IF (ABS(HNEXT)<=TEMP) THEN
       PRINT *, ' The required step size is too small.'
       INDEX = -4
       GOTO 70
    END IF

    ! If the step size is not too small, redo the step with the
    ! estimated step size as the new step size.
    H = HNEXT
    GOTO 30

60  CONTINUE

    ! K8 contains the derivative at T+H so there is no need to
    ! calculate the derivatives for the approximate solution
    ! at this point.
    MYDYNEW(1:N) = KMAT(1:N,9)
    ! The step integration section ends here.

    ! The interpolatory root finding section begins here.

    IF (DEBUG) THEN
       WRITE(DBGUNIT,9999)
       9999 FORMAT(' Begin interpolatory root finding section.')
    END IF

    IF (NG>0) THEN
       ! Calculate the residuals at TNEW, using exact derivatives.
       IF (NGUSER>0) THEN
          IF (NLAGS>0) CALL DDE_BETA(TNEW,MYYNEW,BETA)
          IMETH = 1
          HZSET = TNEW - TOLD
          ITERATE = .FALSE.
          CALL DDE_ZSET(TNEW,MYYNEW,HZSET,YINIT,INDEX,IMETH,KMAT(:,1), &
               KMAT(:,3),KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8), &
               KMAT(:,9),KMAT(:,10),TOLD,MYYOLD,ITERATE)
          ITERATE = .FALSE.
          IF (INDEX==-12 .OR. INDEX==-4) THEN
             PRINT *, ' An error occurred in the third call to'
             PRINT *, ' DDE_ZSET from DDE_DRV3.'
             GOTO 70
          END IF
       END IF
       CALL DDE_GRT3(TNEW,MYYNEW,MYDYNEW,GUSER,GTNEW,BETA)
       NGEVAL = NGEVAL + 1

       ! Do not make another return for a previously reported root.
       ! Do the interpolatory root finding.
       INDEX = 2
       TVALS(1) = TOLD
       TVALS(2) = TNEW
       TVALS(3) = TOLD
       TVALS(4) = TNEW
       IF (IFIRST==1) CALL DDE_GRT5(TVALS,MYYOLD,GUSER,H,INDEX,KMAT(:,1), &
            KMAT(:,3),KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8), &
            KMAT(:,9),KMAT(:,10),BETA,TOLD,YINIT)
       TROOT = TVALS(3)
       TROOT2 = TVALS(4)

       ! Check if the root finder was successful.
       IF (INDEX<0) THEN
          PRINT *, ' The root finder was not successful.'
          GOTO 70
       END IF

       ! If the root occurred at TNEW then do not increase
       ! the step size for the next step.
       IF (INDEX==4) THEN
          HNEXT = SIGN(MIN(ABS(HNEXT),ABS(H)),H)
          IF (DEBUG) THEN
             WRITE(DBGUNIT,9998) HNEXT, H
             9998 FORMAT(' Spot 5 (TNEW ROOT). HNEXT, H = ', 2E20.10)
          END IF
       END IF
       ! If the root did not occur at TNEW then ...
       ! The interpolatory root finding section ends here.

       ! The integration step to the root begins here.

       IF (DEBUG) THEN
          WRITE(DBGUNIT,9997)
          9997 FORMAT(' Perform integration step to the root.')
       END IF

       IF (INDEX==3) THEN

          ! DDE_DRV2 will restart the integration at TROOT.
          H = TROOT - TOLD
          ! HNEXT = H
          ! For the next step choose the larger of the distances from
          ! the root to the endpoints.
          HNEXT = SIGN(MAX(ABS(TROOT-TOLD),ABS(TNEW-TROOT)),H)
          IF (DEBUG) THEN
             WRITE(DBGUNIT,9996) HNEXT, H
             9996 FORMAT(' Spot 6 (TNEW). HNEXT, HLEFT = ', 2E20.10)
          END IF
          ! Call DDE_STEP to calculate the new approximation MYYNEW
          ! at TROOT. Use an heuristic and assume that the error test
          ! will pass since it passed for the original TNEW. The
          ! purpose of this call is to get new Runge-Kutta weights
          ! that correspond to the interval (TOLD,TROOT) rather than
          ! to the original interval (TOLD,TNEW).
          ! NOTE:
          ! The error vector R will be used as scratch storage.
          ITRY = 1

          CALL DDE_STEP(TOLD,H,DERIVS,TNEW,JNDEX,BETA,YINIT,CONVERGED)

          IF (.NOT. (CONVERGED)) THEN
             PRINT *, ' After a root of an event function was located'
             PRINT *, ' and the integration step was repeated up to'
             PRINT *, ' the root, the corrector iteration in DDE_STEP'
             PRINT *, ' did not converge.'
             INDEX = -15
             GOTO 70
          END IF

          ! Check for error return.
          IF (JNDEX<0) THEN
             INDEX = JNDEX
             GOTO 70
          END IF
          TNEW = TROOT

          ! K8 contains the derivative at T+H so there is no need to
          ! recalculate the derivatives for the approximate solution
          ! at this point.
          MYDYNEW(1:N) = KMAT(1:N,9)

          ! Load MYYNEW and MYDYNEW into MYYROOT and MYDYROOT.
          ! (MYYNEW and MYDYNEW have changed since the root
          ! finding was done.)
          MYYROOT(1:N) = MYYNEW(1:N)
          MYDYROOT(1:N) = MYDYNEW(1:N)

          ! NOTE:
          ! Do not recalculate the residuals since that would overwrite
          ! the exact zeroes loaded into G(TROOT) by DDE_GRT1.

       END IF

       ! Update the residuals.
       GTOLD(1:NG) = GTNEW(1:NG)

    END IF
    ! The integration step to the root ends here.

    IF (DEBUG) THEN
       WRITE(DBGUNIT,9995)
       9995 FORMAT(' End integration step to the root.')
    END IF

    ! Update the history queue before returning.
    ! Note:
    ! TNEW is now the first of either the root or the last successful
    ! integration point.
    CALL DDE_QUE1(TNEW,TRIM_GET)

    ! Call DAFTER to allow any end of successful step bookkeeping.
    ! CALL DAFTER(N,MYYNEW,MYDYNEW,TNEW,H)

    ! Update the successful step counter and check if we have
    ! reached the maximum allowable number of steps.
    NSTEPS = NSTEPS + 1
    IF (NSTEPS<0) THEN
       PRINT *, ' NSTEPS wrapped around in subroutine DDE_DRV3.'
       STOP
    END IF
    IF (NSTEPS>=MYMAX_STEPS) THEN
       PRINT *, ' The integration will be terminated since the'
       PRINT *, ' maximum allowable number of steps, MAX_STEPS,'
       PRINT *, ' has been taken without completing the integration.'
       PRINT *, ' The default value of 10,000 steps may be changed'
       PRINT *, ' by adding MAX_STEPS=desired number in your options'
       PRINT *, ' call to DDE_SOLVER.'
       INDEX = -15
       GOTO 70
    END IF

    GOTO 80

70  CONTINUE
    PRINT *, ' A terminal error occurred in DDE_DRV3.'

80  CONTINUE

    RETURN
  END SUBROUTINE DDE_DRV3
  !____________________________________________________________________________

  SUBROUTINE DDE_STEP(TOLD,H,DERIVS,TNEW,INDEX,BETA,YINIT,CONVERGED)

    ! This step integrator must be called by DDE_DRV3.

    ! The purpose of DDE_STEP is to perform a single step of the
    ! integration for DDE_SOLVER. DDE_STEP is called by the interval
    ! oriented driver DDE_DRV3. DDE_STEP extrapolates the most recent
    ! solution polynomial to obtain an initial solution if a delay
    ! time falls beyond the solution history queue. It then performs
    ! a predictor-corrector like iteration of the Runge-Kutta-Sarafyan
    ! solution polynomials.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (INOUT) :: H
    DOUBLE PRECISION, INTENT (OUT) :: TNEW
    DOUBLE PRECISION, INTENT (IN) :: TOLD
    INTEGER, INTENT (OUT) :: INDEX
    LOGICAL, INTENT (OUT) :: CONVERGED
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL BETA, DERIVS, YINIT
    INTERFACE
       SUBROUTINE BETA(T,Y,BVAL)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y
         DOUBLE PRECISION, DIMENSION(:) :: BVAL
         INTENT(IN):: T,Y
         INTENT(OUT) :: BVAL
       END SUBROUTINE BETA
    END INTERFACE
    INTERFACE
       SUBROUTINE DERIVS(T,Y,Z,DY)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y,DY
          DOUBLE PRECISION, DIMENSION(:,:) :: Z
          INTENT(IN):: T,Y,Z
          INTENT(OUT) :: DY
       END SUBROUTINE DERIVS
    END INTERFACE
    INTERFACE
       SUBROUTINE YINIT(T,Y)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y
          INTENT(IN):: T
          INTENT(OUT) :: Y
       END SUBROUTINE YINIT
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: HNM1, KDIFF, KFAC, KRATIO, TNM1, TR, TVAL, XDENOM, XNUMER
    INTEGER :: I, IBEGIN, ICOR, IFAM, IMETH, INDEXO, IQUEUE, N
    LOGICAL :: ITERATE
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    ! ..
    ! Columns 1-10 of KMAT are K0,...,K9.
    ! Columns 11-20 of KMAT are KP0,...KP9.
    INDEX = 2
    N = MYN
    TNEW = TOLD
    ITERATE = .FALSE.
    ICOR = 0
    CONVERGED = .TRUE.

    ! If IFAM = 0, Forward Euler will be used for extrapolations
    ! during the prediction (ICOR = 0), if the delay time falls beyond
    ! the solution history queue. Otherwise, the most recent successful
    ! Y((6,10)) polynomial will be used for the extrapolations during
    ! the prediction, if the delay falls beyond the solution history
    ! queue. In both cases, the most recent Y((6,10)) polynomial iterate
    ! will be used during the corrector iteration (ICOR>0).
    IFAM = 1

    ! If the first step has not been completed, use Euler extrapolations
    ! even if IFAM = 1 above.
    IQUEUE = MYIPOINT(2)
    IF (IQUEUE==0) IFAM = 0

    ! The information for the most recent successful Y((6,10))
    ! is needed if IFAM = 0.
    IF (IQUEUE/=0) THEN
       INDEXO = IQUEUE
       IBEGIN = (INDEXO-1)*10
       HNM1 = TQUEUE(INDEXO) - TQUEUE(INDEXO-1)
       TNM1 = TQUEUE(INDEXO-1)
    END IF

    ! Calculate K0.
    ! Do not re-compute the derivative at TVAL = TOLD since it is
    ! available in MYDYOLD.
    MYYNEW(1:N) = MYYOLD(1:N)
    KMAT(1:N,1) = MYDYOLD(1:N)

10  CONTINUE

    ! Calculate K1.
    TVAL = TOLD + H/9.0D0
    MYYNEW(1:N) = MYYOLD(1:N) + H*(KMAT(1:N,1)/9.0D0)
    IF (NLAGS>0) THEN
       CALL DDE_BETA(TVAL,MYYNEW,BETA)
       IF (ICOR==0 .AND. IFAM/=0) THEN
          ! Most recent successful Y((6,10)) polynomial to extrapolate.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,HNM1,YINIT,INDEX,IMETH, &
               QUEUE(:,IBEGIN+2),QUEUE(:,IBEGIN+3),QUEUE(:,IBEGIN+4), &
               QUEUE(:,IBEGIN+5),QUEUE(:,IBEGIN+6),QUEUE(:,IBEGIN+7), &
               QUEUE(:,IBEGIN+8),QUEUE(:,IBEGIN+9),QUEUE(:,IBEGIN+10),TNM1, &
               QUEUE(:,IBEGIN+1),ITERATE)
       END IF
       IF (ICOR==0 .AND. IFAM==0) THEN
          ! Use Euler extrapolations.
          IMETH = 9
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1),KMAT(:,3), &
               KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9), &
               KMAT(:,10),TOLD,MYYOLD,ITERATE)
       END IF
       IF (ICOR>0) THEN
          ! Most recent Y((6,10)) iterate to correct.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1), &
               KMAT(:,13),KMAT(:,14),KMAT(:,15),KMAT(:,16),KMAT(:,17), &
               KMAT(:,18),KMAT(:,19),KMAT(:,20),TOLD,MYYOLD,ITERATE)
          ITERATE = .FALSE.
       END IF
       IF (INDEX==-12 .OR. INDEX==-4) GOTO 40
    END IF
    CALL DDE_DERV(TVAL,MYYNEW,KMAT(:,2),DERIVS)
    NFEVAL = NFEVAL + 1
    IF (ICOR>0) THEN
       ! Force convergence of H*K1.
       KRATIO = 0.0D0
       DO I = 1, N
          KDIFF = KMAT(I,2) - KMAT(I,12)
          XNUMER = ABS(H)*ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)*ABS(H)*MAX(ABS(KMAT(I, &
               2)),ABS(KMAT(I,12)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
    END IF

    ! Calculate K2.
    TVAL = TOLD + H/6.0D0
    MYYNEW(1:N) = MYYOLD(1:N) + H*((1.0D0/24.0D0)*KMAT(1:N,1)+(3.0D0/ &
         24.0D0)*KMAT(1:N,2))

    IF (NLAGS>0) THEN
       CALL DDE_BETA(TVAL,MYYNEW,BETA)
       IF (ICOR==0 .AND. IFAM/=0) THEN
          ! Most recent successful Y((6,10)) polynomial to extrapolate.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,HNM1,YINIT,INDEX,IMETH, &
               QUEUE(:,IBEGIN+2),QUEUE(:,IBEGIN+3),QUEUE(:,IBEGIN+4), &
               QUEUE(:,IBEGIN+5),QUEUE(:,IBEGIN+6),QUEUE(:,IBEGIN+7), &
               QUEUE(:,IBEGIN+8),QUEUE(:,IBEGIN+9),QUEUE(:,IBEGIN+10),TNM1, &
               QUEUE(:,IBEGIN+1),ITERATE)
       END IF
       IF (ICOR==0 .AND. IFAM==0) THEN
          ! Use Euler extrapolations.
          IMETH = 9
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1),KMAT(:,3), &
               KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9), &
               KMAT(:,10),TOLD,MYYOLD,ITERATE)
       END IF
       IF (ICOR>0) THEN
          ! Most recent Y((6,10)) iterate to correct.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1), &
               KMAT(:,13),KMAT(:,14),KMAT(:,15),KMAT(:,16),KMAT(:,17), &
               KMAT(:,18),KMAT(:,19),KMAT(:,20),TOLD,MYYOLD,ITERATE)
          ITERATE = .FALSE.
       END IF
       IF (INDEX==-12 .OR. INDEX==-4) GOTO 40
    END IF
    CALL DDE_DERV(TVAL,MYYNEW,KMAT(:,3),DERIVS)
    NFEVAL = NFEVAL + 1
    IF (ICOR>0) THEN
       ! Force convergence of H*K2.
       DO I = 1, N
          KDIFF = KMAT(I,3) - KMAT(I,13)
          XNUMER = ABS(H)*ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)*ABS(H)*MAX(ABS(KMAT(I, &
               3)),ABS(KMAT(I,13)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
    END IF

    ! Calculate K3.
    TVAL = TOLD + H/3.0D0
    MYYNEW(1:N) = MYYOLD(1:N) + H*((1.0D0/6.0D0)*KMAT(1:N,1)-(3.0D0/6.0D0) &
         *KMAT(1:N,2)+(4.0D0/6.0D0)*KMAT(1:N,3))
    IF (NLAGS>0) THEN
       CALL DDE_BETA(TVAL,MYYNEW,BETA)
       IF (ICOR==0 .AND. IFAM/=0) THEN
          ! Most recent successful Y((6,10)) polynomial to extrapolate.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,HNM1,YINIT,INDEX,IMETH, &
               QUEUE(:,IBEGIN+2),QUEUE(:,IBEGIN+3),QUEUE(:,IBEGIN+4), &
               QUEUE(:,IBEGIN+5),QUEUE(:,IBEGIN+6),QUEUE(:,IBEGIN+7), &
               QUEUE(:,IBEGIN+8),QUEUE(:,IBEGIN+9),QUEUE(:,IBEGIN+10),TNM1, &
               QUEUE(:,IBEGIN+1),ITERATE)
       END IF
       IF (ICOR==0 .AND. IFAM==0) THEN
          ! Use Euler extrapolations.
          IMETH = 9
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1),KMAT(:,3), &
               KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9), &
               KMAT(:,10),TOLD,MYYOLD,ITERATE)
       END IF
       IF (ICOR>0) THEN
          ! Most recent Y((6,10)) iterate to correct.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1), &
               KMAT(:,13),KMAT(:,14),KMAT(:,15),KMAT(:,16),KMAT(:,17), &
               KMAT(:,18),KMAT(:,19),KMAT(:,20),TOLD,MYYOLD,ITERATE)
          ITERATE = .FALSE.
       END IF
       IF (INDEX==-12 .OR. INDEX==-4) GOTO 40
    END IF
    CALL DDE_DERV(TVAL,MYYNEW,KMAT(:,4),DERIVS)
    NFEVAL = NFEVAL + 1
    IF (ICOR>0) THEN
       ! Force convergence of H*K3.
       DO I = 1, N
          KDIFF = KMAT(I,4) - KMAT(I,14)
          XNUMER = ABS(H)*ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)*ABS(H)*MAX(ABS(KMAT(I, &
               4)),ABS(KMAT(I,14)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
    END IF

    ! Calculate K4.
    TVAL = TOLD + H/2.0D0
    MYYNEW(1:N) = MYYOLD(1:N) + H*((1.0D0/8.0D0)*KMAT(1:N,1)+(3.0D0/8.0D0) &
         *KMAT(1:N,4))
    IF (NLAGS>0) THEN
       CALL DDE_BETA(TVAL,MYYNEW,BETA)
       IF (ICOR==0 .AND. IFAM/=0) THEN
          ! Most recent successful Y((6,10)) polynomial to extrapolate.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,HNM1,YINIT,INDEX,IMETH, &
               QUEUE(:,IBEGIN+2),QUEUE(:,IBEGIN+3),QUEUE(:,IBEGIN+4), &
               QUEUE(:,IBEGIN+5),QUEUE(:,IBEGIN+6),QUEUE(:,IBEGIN+7), &
               QUEUE(:,IBEGIN+8),QUEUE(:,IBEGIN+9),QUEUE(:,IBEGIN+10),TNM1, &
               QUEUE(:,IBEGIN+1),ITERATE)
       END IF
       IF (ICOR==0 .AND. IFAM==0) THEN
          ! Use Euler extrapolations.
          IMETH = 9
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1),KMAT(:,3), &
               KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9), &
               KMAT(:,10),TOLD,MYYOLD,ITERATE)
       END IF
       IF (ICOR>0) THEN
          ! Most recent Y((6,10)) iterate to correct.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1), &
               KMAT(:,13),KMAT(:,14),KMAT(:,15),KMAT(:,16),KMAT(:,17), &
               KMAT(:,18),KMAT(:,19),KMAT(:,20),TOLD,MYYOLD,ITERATE)
          ITERATE = .FALSE.
       END IF
       IF (INDEX==-12 .OR. INDEX==-4) GOTO 40
    END IF
    CALL DDE_DERV(TVAL,MYYNEW,KMAT(:,5),DERIVS)
    NFEVAL = NFEVAL + 1
    IF (ICOR>0) THEN
       ! Force convergence of H*K4.
       DO I = 1, N
          KDIFF = KMAT(I,5) - KMAT(I,15)
          XNUMER = ABS(H)*ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)*ABS(H)*MAX(ABS(KMAT(I, &
               5)),ABS(KMAT(I,15)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
    END IF

    ! Calculate K5.
    TVAL = TOLD + H*(2.0D0/3.0D0)
    KFAC = 2.0D0/21.0D0
    MYYNEW(1:N) = MYYOLD(1:N) + H*(-4.0D0*KFAC*KMAT(1:N,1)+3.0D0*KFAC*KMAT &
         (1:N,2)+12.0D0*KFAC*KMAT(1:N,3)-12.0D0*KFAC*KMAT(1:N,4)+ &
         8.0D0*KFAC*KMAT(1:N,5))
    IF (NLAGS>0) THEN
       CALL DDE_BETA(TVAL,MYYNEW,BETA)
       IF (ICOR==0 .AND. IFAM/=0) THEN
          ! Most recent successful Y((6,10)) polynomial to extrapolate.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,HNM1,YINIT,INDEX,IMETH, &
               QUEUE(:,IBEGIN+2),QUEUE(:,IBEGIN+3),QUEUE(:,IBEGIN+4), &
               QUEUE(:,IBEGIN+5),QUEUE(:,IBEGIN+6),QUEUE(:,IBEGIN+7), &
               QUEUE(:,IBEGIN+8),QUEUE(:,IBEGIN+9),QUEUE(:,IBEGIN+10),TNM1, &
               QUEUE(:,IBEGIN+1),ITERATE)
       END IF
       IF (ICOR==0 .AND. IFAM==0) THEN
          ! Use Euler extrapolations.
          IMETH = 9
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1),KMAT(:,3), &
               KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9), &
               KMAT(:,10),TOLD,MYYOLD,ITERATE)
       END IF
       IF (ICOR>0) THEN
          ! Most recent Y((6,10)) iterate to correct.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1), &
               KMAT(:,13),KMAT(:,14),KMAT(:,15),KMAT(:,16),KMAT(:,17), &
               KMAT(:,18),KMAT(:,19),KMAT(:,20),TOLD,MYYOLD,ITERATE)
          ITERATE = .FALSE.
       END IF
       IF (INDEX==-12 .OR. INDEX==-4) GOTO 40
    END IF
    CALL DDE_DERV(TVAL,MYYNEW,KMAT(:,6),DERIVS)
    NFEVAL = NFEVAL + 1
    IF (ICOR>0) THEN
       ! Force convergence of H*K5.
       DO I = 1, N
          KDIFF = KMAT(I,6) - KMAT(I,16)
          XNUMER = ABS(H)*ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)*ABS(H)*MAX(ABS(KMAT(I, &
               6)),ABS(KMAT(I,16)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
    END IF

    ! Calculate K6.
    TVAL = TOLD + H*(5.0D0/6.0D0)
    KFAC = (5.0D0/456.0D0)
    MYYNEW(1:N) = MYYOLD(1:N) + H*(13.0D0*KFAC*KMAT(1:N,1)+3.0D0*KFAC*KMAT &
         (1:N,2)-16.0D0*KFAC*KMAT(1:N,3)+66.0D0*KFAC*KMAT(1:N,4)- &
         32.0D0*KFAC*KMAT(1:N,5)+42.0D0*KFAC*KMAT(1:N,6))
    IF (NLAGS>0) THEN
       CALL DDE_BETA(TVAL,MYYNEW,BETA)
       IF (ICOR==0 .AND. IFAM/=0) THEN
          ! Most recent successful Y((6,10)) polynomial to extrapolate.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,HNM1,YINIT,INDEX,IMETH, &
               QUEUE(:,IBEGIN+2),QUEUE(:,IBEGIN+3),QUEUE(:,IBEGIN+4), &
               QUEUE(:,IBEGIN+5),QUEUE(:,IBEGIN+6),QUEUE(:,IBEGIN+7), &
               QUEUE(:,IBEGIN+8),QUEUE(:,IBEGIN+9),QUEUE(:,IBEGIN+10),TNM1, &
               QUEUE(:,IBEGIN+1),ITERATE)
       END IF
       IF (ICOR==0 .AND. IFAM==0) THEN
          ! Use Euler extrapolations.
          IMETH = 9
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1),KMAT(:,3), &
               KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9), &
               KMAT(:,10),TOLD,MYYOLD,ITERATE)
       END IF
       IF (ICOR>0) THEN
          ! Most recent Y((6,10)) iterate to correct.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1), &
               KMAT(:,13),KMAT(:,14),KMAT(:,15),KMAT(:,16),KMAT(:,17), &
               KMAT(:,18),KMAT(:,19),KMAT(:,20),TOLD,MYYOLD,ITERATE)
          ITERATE = .FALSE.
       END IF
       IF (INDEX==-12 .OR. INDEX==-4) GOTO 40
    END IF
    CALL DDE_DERV(TVAL,MYYNEW,KMAT(:,7),DERIVS)
    NFEVAL = NFEVAL + 1
    IF (ICOR>0) THEN
       ! Force convergence of H*K6.
       DO I = 1, N
          KDIFF = KMAT(I,7) - KMAT(I,17)
          XNUMER = ABS(H)*ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)*ABS(H)*MAX(ABS(KMAT(I, &
               7)),ABS(KMAT(I,17)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
    END IF

    ! Calculate K7.
    TVAL = TOLD + H
    MYYNEW(1:N) = MYYOLD(1:N) + H*((35.0D0/322.0D0)*KMAT(1:N,1)-(81.0D0/ &
         322.0D0)*KMAT(1:N,2)+(240.0D0/322.0D0)*KMAT(1:N,3)- &
         (360.0D0/322.0D0)*KMAT(1:N,4)+(680.0D0/322.0D0)*KMAT(1:N,5)- &
         (420.0D0/322.0D0)*KMAT(1:N,6)+(228.0D0/322.0D0)*KMAT(1:N,7))
    IF (NLAGS>0) THEN
       CALL DDE_BETA(TVAL,MYYNEW,BETA)
       IF (ICOR==0 .AND. IFAM/=0) THEN
          ! Most recent successful Y((6,10)) polynomial to extrapolate.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,HNM1,YINIT,INDEX,IMETH, &
               QUEUE(:,IBEGIN+2),QUEUE(:,IBEGIN+3),QUEUE(:,IBEGIN+4), &
               QUEUE(:,IBEGIN+5),QUEUE(:,IBEGIN+6),QUEUE(:,IBEGIN+7), &
               QUEUE(:,IBEGIN+8),QUEUE(:,IBEGIN+9),QUEUE(:,IBEGIN+10),TNM1, &
               QUEUE(:,IBEGIN+1),ITERATE)
       END IF
       IF (ICOR==0 .AND. IFAM==0) THEN
          ! Use Euler extrapolations.
          IMETH = 9
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1),KMAT(:,3), &
               KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9), &
               KMAT(:,10),TOLD,MYYOLD,ITERATE)
       END IF
       IF (ICOR>0) THEN
          ! Most recent Y((6,10)) iterate to correct.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1), &
               KMAT(:,13),KMAT(:,14),KMAT(:,15),KMAT(:,16),KMAT(:,17), &
               KMAT(:,18),KMAT(:,19),KMAT(:,20),TOLD,MYYOLD,ITERATE)
          ITERATE = .FALSE.
       END IF
       IF (INDEX==-12 .OR. INDEX==-4) GOTO 40
    END IF
    CALL DDE_DERV(TVAL,MYYNEW,KMAT(:,8),DERIVS)
    NFEVAL = NFEVAL + 1
    IF (ICOR>0) THEN
       ! Force convergence of H*K7.
       DO I = 1, N
          KDIFF = KMAT(I,8) - KMAT(I,18)
          XNUMER = ABS(H)*ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)*ABS(H)*MAX(ABS(KMAT(I, &
               8)),ABS(KMAT(I,18)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
    END IF

    ! Calculate K8.
    TVAL = TOLD + H
    MYYNEW(1:N) = MYYOLD(1:N) + H*((161.0D0/3000.0D0)*KMAT(1:N,1)+(57.0D0/ &
         250.0D0)*KMAT(1:N,3)+(21.0D0/200.0D0)*KMAT(1:N,4)+ &
         (17.0D0/75.0D0)*KMAT(1:N,5)+(21.0D0/200.0D0)*KMAT(1:N,6)+ &
         (57.0D0/250.0D0)*KMAT(1:N,7)+(161.0D0/3000.0D0)*KMAT(1:N,8))
    IF (NLAGS>0) THEN
       CALL DDE_BETA(TVAL,MYYNEW,BETA)
       IF (ICOR==0 .AND. IFAM/=0) THEN
          ! Most recent successful Y((6,10)) polynomial to extrapolate.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,HNM1,YINIT,INDEX,IMETH, &
               QUEUE(:,IBEGIN+2),QUEUE(:,IBEGIN+3),QUEUE(:,IBEGIN+4), &
               QUEUE(:,IBEGIN+5),QUEUE(:,IBEGIN+6),QUEUE(:,IBEGIN+7), &
               QUEUE(:,IBEGIN+8),QUEUE(:,IBEGIN+9),QUEUE(:,IBEGIN+10),TNM1, &
               QUEUE(:,IBEGIN+1),ITERATE)
       END IF
       IF (ICOR==0 .AND. IFAM==0) THEN
          ! Use EULER extrapolations.
          IMETH = 9
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1),KMAT(:,3), &
               KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9), &
               KMAT(:,10),TOLD,MYYOLD,ITERATE)
       END IF
       IF (ICOR>0) THEN
          ! Most recent Y((6,10)) iterate to correct.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1), &
               KMAT(:,13),KMAT(:,14),KMAT(:,15),KMAT(:,16),KMAT(:,17), &
               KMAT(:,18),KMAT(:,19),KMAT(:,20),TOLD,MYYOLD,ITERATE)
          ITERATE = .FALSE.
       END IF
       IF (INDEX==-12 .OR. INDEX==-4) GOTO 40
    END IF
    CALL DDE_DERV(TVAL,MYYNEW,KMAT(:,9),DERIVS)
    NFEVAL = NFEVAL + 1
    IF (ICOR>0) THEN
       ! Force convergence of H*K8.
       DO I = 1, N
          KDIFF = KMAT(I,9) - KMAT(I,19)
          XNUMER = ABS(H)*ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)*ABS(H)*MAX(ABS(KMAT(I, &
               9)),ABS(KMAT(I,19)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
    END IF

    ! Calculate K9.
    TVAL = TOLD + H*(2.0D0/3.0D0)
    MYYNEW(1:N) = MYYOLD(1:N) + H*((7.0D0/135.0D0)*KMAT(1:N,1)+(32.0D0/ &
         135.0D0)*KMAT(1:N,3)+(12.0D0/135.0D0)*KMAT(1:N,4)+ &
         (32.0D0/135.0D0)*KMAT(1:N,5)+(7.0D0/135.0D0)*KMAT(1:N,6))
    IF (NLAGS>0) THEN
       CALL DDE_BETA(TVAL,MYYNEW,BETA)
       IF (ICOR==0 .AND. IFAM/=0) THEN
          ! Most recent successful Y((6,10)) polynomial to extrapolate.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,HNM1,YINIT,INDEX,IMETH, &
               QUEUE(:,IBEGIN+2),QUEUE(:,IBEGIN+3),QUEUE(:,IBEGIN+4), &
               QUEUE(:,IBEGIN+5),QUEUE(:,IBEGIN+6),QUEUE(:,IBEGIN+7), &
               QUEUE(:,IBEGIN+8),QUEUE(:,IBEGIN+9),QUEUE(:,IBEGIN+10),TNM1, &
               QUEUE(:,IBEGIN+1),ITERATE)
       END IF
       IF (ICOR==0 .AND. IFAM==0) THEN
          ! Use Euler extrapolations.
          IMETH = 9
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1),KMAT(:,3), &
               KMAT(:,4),KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9), &
               KMAT(:,10),TOLD,MYYOLD,ITERATE)
       END IF
       IF (ICOR>0) THEN
          ! Most recent Y((6,10)) iterate to correct.
          IMETH = 1
          CALL DDE_ZSET(TVAL,MYYNEW,H,YINIT,INDEX,IMETH,KMAT(:,1), &
               KMAT(:,13),KMAT(:,14),KMAT(:,15),KMAT(:,16),KMAT(:,17), &
               KMAT(:,18),KMAT(:,19),KMAT(:,20),TOLD,MYYOLD,ITERATE)
          ITERATE = .FALSE.
       END IF
       IF (INDEX==-12 .OR. INDEX==-4) GOTO 40
    END IF
    CALL DDE_DERV(TVAL,MYYNEW,KMAT(:,10),DERIVS)
    NFEVAL = NFEVAL + 1
    IF (ICOR>0) THEN
       ! Force convergence of H*K9.
       DO I = 1, N
          KDIFF = KMAT(I,10) - KMAT(I,20)
          XNUMER = ABS(H)*ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)*ABS(H)*MAX(ABS(KMAT(I, &
               10)),ABS(KMAT(I,20)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
    END IF

    ! Calculate MYYNEW.
    TVAL = TOLD + H
    CALL DDE_POLY(TOLD,MYYOLD,H,KMAT(:,1),KMAT(:,3),KMAT(:,4),KMAT(:,5), &
         KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9),KMAT(:,10),MYYNEW, &
         MYYNEW,TVAL,1)

    IF (ICOR>0) THEN
       ! Force convergence of MYYNEW.
       DO I = 1, N
          KDIFF = MYYNEW(I) - KMAT(I,11)
          XNUMER = ABS(H)*ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)* &
                   MAX(ABS(MYYNEW(I)),ABS(KMAT(I,11)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
       ! Force convergence of MYDYNEW if this is a neutral problem.
       IF (.NOT. MYNEUTRAL) GOTO 20
       DO I = 1, N
          KDIFF = KMAT(I,9) - KMAT(I,19)
          XNUMER = ABS(KDIFF)
          XDENOM = MYABSER(I) + MYRELER(I)*MAX(ABS(KMAT(I,9)), &
                                               ABS(KMAT(I,19)))
          IF (XDENOM<=0.0D0) XDENOM = 1.0D0
          TR = XNUMER/XDENOM
          KRATIO = MAX(KRATIO,TR)
       END DO
20     CONTINUE
    END IF

    IF (ICOR==0 .AND. (.NOT. ITERATE)) THEN
       ! A correction will not be performed since no extrapolations
       ! were performed.
       GOTO 30
    END IF

    ! Test for convergence of the derivative approximations.
    IF (ICOR>0) THEN
       IF (KRATIO<=KCFAC) THEN
          ! Do not perform any further corrections.
          GOTO 30
       ELSE
          IF (ICOR==MAXCOR) THEN
             ! Update the failed iteration counter.
             NFAILC = NFAILC + 1
             ! Reduce the step size and try again.
             CONVERGED = .FALSE.
             GOTO 40
          END IF
       END IF
    END IF

    ! Prepare to do a correction.
    ICOR = ICOR + 1
    IF (ICOR<=MAXCOR) THEN
       KMAT(1:N,11) = MYYNEW(1:N)
       KMAT(1:N,12:20) = KMAT(1:N,2:10)
       ! And do it
       GOTO 10
    END IF

30  CONTINUE

    ! Update the independent variable.
    TNEW = TOLD + H

    ! The step was successful.

40  CONTINUE

    RETURN
  END SUBROUTINE DDE_STEP
  !____________________________________________________________________________
  SUBROUTINE DDE_HST1(DERIVS,A,B,Y,YPRIME,ETOL,MORDER,BIG,SPY,PV,YP,SF, &
       BETA,YINIT,INDEX,H)

    ! The function of DDE_HST1 is to compute a starting step size to
    ! be used in solving ode initial value problems. It is based on a
    ! modification of Buddy Watt's HSTART from the SLATEC library. It
    ! is called by DDE_DRV2 if HINIT = 0.0D0.

    ! DDE_HST1 computes a starting step size to be used by DDE_SOLVER.
    ! DDE_HST1 is based on an estimate of the local Lipschitz constant
    ! for the differential equation (lower bound on a norm of the
    ! Jacobian), a bound on the differential equation first derivative),
    ! and a bound on the partial derivative of the equation with respect
    ! to the independent variable (all approximated near the initial
    ! point A).

    ! On input you must provide the following:
    !      DERIVS -- Derivative subroutine.

    !      A -- This is the initial point of integration.

    !      B -- This is a value of the independent variable used to define
    !             the direction of integration. A reasonable choice is to
    !             set B to the first point at which a solution is desired.
    !             You can also use B, if necessary, to restrict the length
    !             of the first integration step because the algorithm will
    !             not compute a starting step length which is bigger than
    !             ABS(B-A), unless B has been chosen too close to A.
    !             (It is presumed that DDE_HST1 has been called with B
    !             different from A on the machine being used. Also see the
    !             discussion about the parameter SMALL.)

    !      Y(*) -- This is the vector of initial values of the solution.

    !      YPRIME(*) -- This is the vector of derivatives.
    !      ETOL -- This is the vector of error tolerances corresponding to
    !             the solution components. It is assumed that all
    !             elements are positive. Following the first integration
    !             step, the tolerances are expected to be used by the
    !             integrator in an error test which roughly requires that
    !                       ABS(LOCAL ERROR)<=ETOL
    !             for each vector component. For consistency with the error
    !             control in DDE_DRV1, It is recommended that ETOL be
    !                ETOL(I) = 0.5D0*(MYABSER + MYRELER*ABS(Y(I)))
    !             where MYABSER and MYRELER are the absolute and relative
    !             error tolerances for DDE_DRV1. This version uses the
    !             first column of the KMAT matrix for ETOL.
    !      MORDER -- This is the order of the formula which will be used by
    !             the initial value method for taking the first integration
    !             step. MORDER should be 6 for DDE_DRV1.
    !      SMALL HAS BEEN SET TO UROUND.
    !      SMALL -- This is a small positive machine dependent constant
    !             which is used for protecting against computations with
    !             numbers which are too small relative to the precision
    !             of floating point arithmetic. Small should be set to
    !             (approximately) the smallest positive double precision
    !             number such that (1.+SMALL)>1. on the machine being
    !             used. The quantity SMALL**(3/8) is used in computing
    !             increments of variables for approximating derivatives
    !             by differences. Also the algorithm will not compute
    !             a starting step length which is smaller than
    !             100*SMALL*ABS(A).
    !      BIG -- This is a large positive machine dependent constant which
    !             is used for preventing machine overflows. A reasonable
    !             choice is to set BIG to (approximately) the square root
    !             of the largest double precision number which can be
    !             held in the machine.
    !      SPY(*),PV(*),YP(*),SF(*)-- These are double precision work
    !             arrays of length MYN which provide the routine with
    !             needed storage space. This version uses columns 2-5
    !             of the KMAT array for these arrays.
    !      BETA, YINIT -- These are the corresponding three external
    !             subroutines for DDE_DRV1. Refer to the documentation
    !             for DDE_DRV1 for a description of these subroutines.
    ! On output  (after the return from DDE_HST1),
    !      H -- is an appropriate starting step size to be attempted
    !             by the differential equation method if INDEX = 0.
    !      INDEX -- Error flag. If INDEX = 0, DDE_HST1 was successful.
    !           Otherwise an error was encountered. In the latter
    !           case, H does not contains an appropriate starting
    !           step size.
    !           Possible error returns:
    !           INDEX = -1 means a delay time was encountered which
    !                      was beyond the corresponding time
    !           INDEX = -2 means B - A is equal to 0.0D0. (This will
    !                      not happen if DDE_HST1 is called from
    !                      DDE_DRV2.)
    !           INDEX =-14 means the integration was terminated
    !                      in one of the user subroutines called
    !                      by DDE_HST1.
    !           All input parameters in the call list remain unchanged
    !           except for the working arrays SPY(*),PV(*),YP(*),SF(*).
    !     Note:
    !     If a delay time falls beyond the initial point, DDE_HST1 will
    !     use forward Euler extrapolations methods to approximate
    !     the delayed solution and delayed derivative. Otherwise it
    !     will use subroutine YINIT.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: A, B, BIG
    DOUBLE PRECISION, INTENT (OUT) :: H
    INTEGER, INTENT (OUT) :: INDEX
    INTEGER, INTENT (IN) :: MORDER
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: ETOL(:)
    DOUBLE PRECISION, INTENT (INOUT) :: PV(:), SF(:), SPY(:), Y(:), YP(:), &
         YPRIME(:)
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL BETA, DERIVS, YINIT
    INTERFACE
       SUBROUTINE BETA(T,Y,BVAL)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y
          DOUBLE PRECISION, DIMENSION(:) :: BVAL
          INTENT(IN):: T,Y
          INTENT(OUT) :: BVAL
       END SUBROUTINE BETA
    END INTERFACE
    INTERFACE
       SUBROUTINE DERIVS(T,Y,Z,DY)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y,DY
          DOUBLE PRECISION, DIMENSION(:,:) :: Z
          INTENT(IN):: T,Y,Z
          INTENT(OUT) :: DY
       END SUBROUTINE DERIVS
    END INTERFACE
    INTERFACE
       SUBROUTINE YINIT(T,Y)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y
         INTENT(IN):: T
         INTENT(OUT) :: Y
      END SUBROUTINE YINIT
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: ABSDX, DA, DELF, DELY, DFDUB, DFDXB, DX, DY, FBND, &
         RELPER, SRYDPB, TOLEXP, TOLMIN, TOLP, TOLSUM, TVAL, YDPB
    INTEGER :: J, K, LK
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, DBLE, FLOAT, LOG10, MAX, MAXVAL, MIN, MIN0, SIGN, SQRT
    ! ..
    ! Begin block permitting; exits to 160.
    INDEX = 0
    DX = B - A
    IF (ABS(DX)<=0.0D0) THEN
       INDEX = -2
       GOTO 120
    END IF
    ABSDX = ABS(DX)
    RELPER = UROUND**0.375D0

    ! Compute an approximate bound(DFDXB) on the partial derivative of
    ! the equation with respect to the independent variable. Protect
    ! against an overflow. Also compute a bound (FBND) on the first
    ! derivative locally.

    DA = SIGN(MAX(MIN(RELPER*ABS(A),ABSDX),100.0D0*UROUND*ABS(A)),DX)
    IF (ABS(DA)<=0.0D0) DA = RELPER*DX
    IF (NLAGS>0) CALL DDE_BETA(A+DA,Y,BETA)
    TVAL = A + DA
    ! IMETH = 9
    ! CALL DDE_HST2(Y,DA,INDEX,YINIT,A,YPRIME,YPRIME,YPRIME,YPRIME, &
    ! YPRIME,YPRIME,YPRIME,YPRIME,YPRIME,IMETH)
    CALL DDE_HST2(Y,INDEX,YINIT,A,YPRIME)
    INDEX = 0
    CALL DDE_DERV(TVAL,Y,SF,DERIVS)
    YP(1:MYN) = SF(1:MYN) - YPRIME(1:MYN)
    DELF = MAXVAL(ABS(YP(1:MYN)))
    DFDXB = BIG
    IF (DELF<BIG*ABS(DA)) DFDXB = DELF/ABS(DA)
    FBND = MAXVAL(ABS(SF(1:MYN)))

    ! Compute an estimate (DFDUB) of the local Lipschitz constant for
    ! the system of differential equations. This also represents an
    ! estimate of the norm of the Jacobian locally. Three iterations
    ! (two when MYN=1) are used to estimate the Lipschitz constant by
    ! numerical differences. The first perturbation vector is based on
    ! the initial derivatives and direction of integration. The second
    ! perturbation vector is formed using another evaluation of the
    ! differential equation. The third perturbation vector is formed
    ! using perturbations based only on the initial values. Components
    ! that are zero are always changed to non-zero values (except on
    ! the first iteration). When information is available, care is
    ! taken to ensure that components of the perturbation vector have
    ! signs which are consistent with the slopes of local solution
    ! curves. Also choose the largest bound (FBND) for the first
    ! derivative. The perturbation vector size is held constant for
    ! all iterations.

    ! Compute the change from the size of the vector of initial
    ! values.
    DELY = RELPER*MAXVAL(ABS(Y(1:MYN)))
    IF (ABS(DELY)<=0.0D0) DELY = RELPER
    DELY = SIGN(DELY,DX)
    DELF = MAXVAL(ABS(YPRIME(1:MYN)))
    FBND = MAX(FBND,DELF)
    IF (ABS(DELF)<=0.0D0) GOTO 10
    ! Use initial derivatives for first perturbation.
    SPY(1:MYN) = YPRIME(1:MYN)
    YP(1:MYN) = YPRIME(1:MYN)
    GOTO 20
10  CONTINUE
    ! Cannot have a null perturbation vector.
    SPY(1:MYN) = 0.0D0
    YP(1:MYN) = 1.0D0
    DELF = MAXVAL(ABS(YP(1:MYN)))
20  CONTINUE
    DFDUB = 0.0D0
    LK = MIN0(MYN+1,3)
    DO K = 1, LK

       ! Define the perturbed vector of initial values.
       PV(1:MYN) = Y(1:MYN) + DELY*(YP(1:MYN)/DELF)
       IF (K==2) GOTO 30

       ! Evaluate derivatives associated with the perturbed vector
       ! and compute the corresponding differences.
       IF (NLAGS>0) THEN
          CALL DDE_BETA(A,PV,BETA)
          DO J = 1, NLAGS
             TVAL = MYBVAL(J)
             IF ((TVAL<=A .AND. B>A) .OR. (TVAL>=A .AND. B<A)) THEN
                CALL DDE_YINIT(TVAL,YINIT)
                ZARRAY(1:MYN,J) = VTEMP3(1:MYN)
                IF (MYNEUTRAL) DARRAY(1:MYN,J) = VTEMP3(MYN+1:2*MYN)
             ELSE
                ZARRAY(1:MYN,J) = Y(1:MYN)
                IF (MYNEUTRAL) DARRAY(1:MYN,J) = YPRIME(1:MYN)
             END IF
          END DO
       END IF
       INDEX = 0
       CALL DDE_DERV(A,PV,YP,DERIVS)
       PV(1:MYN) = YP(1:MYN) - YPRIME(1:MYN)
       GOTO 40
30     CONTINUE

       ! Use a shifted value of the independent variable in computing
       ! one estimate.
       IF (NLAGS>0) CALL DDE_BETA(A+DA,PV,BETA)
       TVAL = A + DA
       ! IMETH = 9
       ! CALL DDE_HST2(Y,DA,INDEX,YINIT,A,YPRIME,YPRIME,YPRIME,YPRIME, &
       ! YPRIME,YPRIME,YPRIME,YPRIME,YPRIME,IMETH)
       CALL DDE_HST2(Y,INDEX,YINIT,A,YPRIME)
       INDEX = 0
       CALL DDE_DERV(TVAL,PV,YP,DERIVS)
       PV(1:MYN) = YP(1:MYN) - SF(1:MYN)
40     CONTINUE

       ! Choose largest bounds on the first derivative and a local
       ! Lipschitz constant.
       FBND = MAX(FBND,MAXVAL(ABS(YP(1:MYN))))
       DELF = MAXVAL(ABS(PV(1:MYN)))
       ! Exit
       IF (DELF>=BIG*ABS(DELY)) GOTO 70
       DFDUB = MAX(DFDUB,DELF/ABS(DELY))
       ! Exit
       IF (K==LK) GOTO 80

       ! Choose the next perturbation vector.
       IF (ABS(DELF)<=0.0D0) DELF = 1.0D0
       DO J = 1, MYN
          IF (K==2) GOTO 50
          DY = ABS(PV(J))
          IF (ABS(DY)<=0.0D0) DY = DELF
          GOTO 60
50        CONTINUE
          DY = Y(J)
          IF (ABS(DY)<=0.0D0) DY = DELY/RELPER
60        CONTINUE
          IF (ABS(SPY(J))<=0.0D0) SPY(J) = YP(J)
          IF (ABS(SPY(J))>0.0D0) DY = SIGN(DY,SPY(J))
          YP(J) = DY
       END DO
       DELF = MAXVAL(ABS(YP(1:MYN)))
    END DO
70  CONTINUE

    ! Protect against an overflow.
    DFDUB = BIG
80  CONTINUE

    ! Compute a bound (YDPB) on the norm of the second derivative.
    YDPB = DFDXB + DFDUB*FBND

    ! Define the tolerance parameter upon which the starting step
    ! size is to be based. A value in the middle of the error
    ! tolerance range is selected.
    TOLMIN = BIG
    TOLSUM = 0.0D0
    DO K = 1, MYN
       TOLEXP = LOG10(ETOL(K))
       TOLMIN = MIN(TOLMIN,TOLEXP)
       TOLSUM = TOLSUM + TOLEXP
    END DO
    TOLP = 10.0D0**(0.5D0*(TOLSUM/DBLE(FLOAT(MYN))+TOLMIN)/ &
         DBLE(FLOAT(MORDER+1)))

    ! Compute a starting step size based on the above first and
    ! second derivative information.
    ! Restrict the step length to be not bigger than ABS(B-A)
    ! (unless B is too close to A).
    H = ABSDX
    IF (ABS(YDPB)>0.0D0 .OR. ABS(FBND)>0.0D0) GOTO 90

    ! Both first derivative term (FBND) and second derivative term
    ! (YDPB) are zero.
    IF (TOLP<1.0D0) H = ABSDX*TOLP
    GOTO 110
90  CONTINUE

    IF (ABS(YDPB)>0.0D0) GOTO 100

    ! Only the second derivative term (YDPB) is zero.
    IF (TOLP<FBND*ABSDX) H = TOLP/FBND
    GOTO 110
100 CONTINUE

    ! The second derivative term (YDPB) is non-zero.
    SRYDPB = SQRT(0.5D0*YDPB)
    IF (TOLP<SRYDPB*ABSDX) H = TOLP/SRYDPB
110 CONTINUE

    ! Further restrict the step length to be not bigger
    ! than 1/DFDUB.
    IF (H*DFDUB>1.0D0) H = 1.0D0/DFDUB

    ! Finally, restrict the step length to be not smaller than
    ! 100*UROUND*ABS(A). However, if A=0. and the computed H
    ! underflowed to zero, the algorithm returns UROUND*ABS(B)
    ! for the step length.
    H = MAX(H,100.0D0*UROUND*ABS(A))
    IF (ABS(H)<=0.0D0) H = UROUND*ABS(B)

    ! Now set the direction of integration.
    H = SIGN(H,DX)

    ! Normal returns.
    RETURN

    ! Error returns.
120 CONTINUE

    PRINT *, ' A - B is equal to zero in DDE_HST1.'

    RETURN
  END SUBROUTINE DDE_HST1
  !____________________________________________________________________________
  SUBROUTINE DDE_HST2(YOLD,INDEX,YINIT,TOLD,K0)

    ! The purpose of DDE_HST2 is to interpolate the delayed solution
    ! and derivative during the first step of the integration for the
    ! DDE_HST1 initial step size selection subroutine for DDE_SOLVER.
    ! It is called only by DDE_HST1 if a delay time falls beyond the
    ! initial time on the first integration step. It calculates
    ! ZARRAY and DARRAY if MYNEUTRAL is .TRUE. ZARRAY(J,I) is the
    ! approximate solution for component J evaluated at the Ith delay
    ! time, MYBVAL(I). DARRAY(J,I) is the approximate derivative for
    ! component J evaluated at the Ith delay time, MYBVAL(I). Euler
    ! extrapolations are used if necessary.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: TOLD
    INTEGER, INTENT (OUT) :: INDEX
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: K0(:), YOLD(:)
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL YINIT
    INTERFACE
       SUBROUTINE YINIT(T,Y)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y
         INTENT(IN):: T
         INTENT(OUT) :: Y
      END SUBROUTINE YINIT
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: TVAL
    INTEGER :: I
    ! ..
    INDEX = 2

    ! Do nothing if NLAGS = 0.
    IF (NLAGS<=0) RETURN

    DO I = 1, NLAGS
       ! Define the Ith delay time.
       TVAL = MYBVAL(I)
       IF (TVAL<=MYTINIT) THEN
          ! Use the history function.
          CALL DDE_YINIT(TVAL,YINIT)
          ZARRAY(1:MYN,I) = VTEMP3(1:MYN)
          IF (MYNEUTRAL) DARRAY(1:MYN,I) = VTEMP3(MYN+1:2*MYN)
       ELSE
          ! Calculate Euler extrapolations.
          ZARRAY(1:MYN,I) = YOLD(1:MYN) + (TVAL-TOLD)*K0(1:MYN)
          DARRAY(1:MYN,I) = K0(1:MYN)
       END IF
    END DO

    RETURN
  END SUBROUTINE DDE_HST2
  !____________________________________________________________________________

  SUBROUTINE DDE_HST3(IERTST,RATIO,H,HNEXT)

    ! The purpose of DDE_HST3 is to calculate the step size for
    ! DDE_SOLVER.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: H
    DOUBLE PRECISION, INTENT (OUT) :: HNEXT
    DOUBLE PRECISION, INTENT (INOUT) :: RATIO
    INTEGER, INTENT (OUT) :: IERTST
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: DORDER, HEXP, RATMAX, RATMIN
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, DBLE, FLOAT, MAX, MIN, SIGN
    ! ..
    HNEXT = H
    IERTST = 0

    RATIO = RATIO*ABS(H)
    ! If RATIO is greater than 1.0D0, the error is too large.
    IF (RATIO>1.0D0) GOTO 10

    ! The estimated error for this step is small enough.
    ! Estimate the step size for the next step. Limit the
    ! increase to a factor of at most FACUP and a decrease
    ! to a factor of at most FACDN.

    IERTST = 0
    DORDER = DBLE(FLOAT(IORDER))
    RATMAX = EFAC*(FACDN**DORDER)
    RATMIN = EFAC/(FACUP**DORDER)
    RATIO = MAX(RATIO,RATMIN)
    RATIO = MIN(RATIO,RATMAX)
    HEXP = 1.0D0/DORDER
    HNEXT = ((EFAC/RATIO)**HEXP)*ABS(H)
    HNEXT = SIGN(HNEXT,H)
    IF (DEBUG) THEN
       WRITE(DBGUNIT,9999) HNEXT, H
       9999 FORMAT(' Spot 7 (HST3). HNEXT, H = ', 2E20.10)
       END IF
    GOTO 20

10  CONTINUE

    ! The estimated error for this step is too large.
    ! Estimate a smaller step size with which to redo
    ! the step. Limit the decrease to a factor of at
    ! most FACDN.
    IERTST = 1
    DORDER = DBLE(FLOAT(IORDER))
    RATMAX = EFAC*(FACDN**DORDER)
    RATIO = MIN(RATIO,RATMAX)
    HEXP = 1.0D0/DORDER
    HNEXT = ((EFAC/RATIO)**HEXP)*ABS(H)
    HNEXT = SIGN(HNEXT,H)

20  CONTINUE

    RETURN
  END SUBROUTINE DDE_HST3
  !____________________________________________________________________________

  SUBROUTINE DDE_HST4(HNEXT,HMAX)

    ! The purpose of DDE_HST4 is to limit the step size to be no
    ! larger than the maximum step size for DDE_SOLVER.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: HMAX
    DOUBLE PRECISION, INTENT (INOUT) :: HNEXT
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: HMAXN
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, MIN, SIGN
    ! ..
    HMAXN = ABS(HMAX)
    HNEXT = SIGN(MIN(ABS(HNEXT),HMAXN),HNEXT)

    RETURN
  END SUBROUTINE DDE_HST4
  !____________________________________________________________________________
  
  SUBROUTINE DDE_QUE1(TNEW,TRIM_GET)

    ! The purpose of DDE_QUE1 is to update the solution history queue
    ! for DDE_SOLVER. It is called by DDE_DRV3 following each
    ! successful integration step.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: TNEW
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: HTEMP
    INTEGER :: IBEGIN, IQUEUE, JQUEUE, NUMPT

    !EXTERNAL TRIM_GET
    INTERFACE
       SUBROUTINE TRIM_GET
       END SUBROUTINE TRIM_GET
    END INTERFACE
    ! ..
    ! IQUEUE points to the most recent addition to the queue.
    ! JQUEUE points to the oldest.
    ! MYN is the number of ddes in the system.
    ! LQUEUE = length of queue.
    ! NUMPT = LQUEUE/(10*MYN) is the maximum number of points for
    ! which the solution can be saved.
    IQUEUE = MYIPOINT(2)
    JQUEUE = MYIPOINT(3)
    NUMPT = LQUEUE/(10*MYN)

    ! Do not duplicate a time in the queue.
    HTEMP = TNEW - TQUEUE(IQUEUE)
    IF (ABS(HTEMP)<=0D0) THEN
       RETURN
    END IF

    ! Update the pointer to the most recent information in the
    ! circular queue.
    IQUEUE = IQUEUE + 1

    IF (IQUEUE>NUMPT) THEN
       ! Cannot get to the next statement if everything is
       ! working correctly:
       PRINT *, ' The QUEUE wrapped around. This is not permitted'
       PRINT *, ' in this version of DDE_SOLVER.'
       STOP
    END IF

    ! Update the pointer information with the new queue pointer
    ! information.
    MYIPOINT(2) = IQUEUE
    MYIPOINT(3) = JQUEUE

    ! Update the independent variable information.
    TQUEUE(IQUEUE) = TNEW

    ! Store the integration information for (TOLD,TNEW) in the
    ! queue. There is no need to save DYOLD since K0 = DYOLD.
    ! There is no need to store K1 since it is not used in the
    ! interpolation subroutines.

    ! Update the solution queue.
    IBEGIN = (IQUEUE-1)*10
    QUEUE(1:MYN,IBEGIN+1) = MYYOLD(1:MYN)
    QUEUE(1:MYN,IBEGIN+2) = KMAT(1:MYN,1)
    QUEUE(1:MYN,IBEGIN+3) = KMAT(1:MYN,3)
    QUEUE(1:MYN,IBEGIN+4) = KMAT(1:MYN,4)
    QUEUE(1:MYN,IBEGIN+5) = KMAT(1:MYN,5)
    QUEUE(1:MYN,IBEGIN+6) = KMAT(1:MYN,6)
    QUEUE(1:MYN,IBEGIN+7) = KMAT(1:MYN,7)
    QUEUE(1:MYN,IBEGIN+8) = KMAT(1:MYN,8)
    QUEUE(1:MYN,IBEGIN+9) = KMAT(1:MYN,9)
    QUEUE(1:MYN,IBEGIN+10) = KMAT(1:MYN,10)

    ! Trim the solution queue.
    IF (NSTEPS > 0 .AND. MY_TRIM_FREQUENCY > 0) THEN
       IF (MY_MAX_DELAY > 0.0D0 .AND. &
             (NSTEPS/MY_TRIM_FREQUENCY)*MY_TRIM_FREQUENCY == NSTEPS) THEN
          CALL TRIM_QUEUE(TRIM_GET)
          IQUEUE = MYIPOINT(2)
       END IF
    END IF

    RETURN
  END SUBROUTINE DDE_QUE1
  !____________________________________________________________________________

  SUBROUTINE DDE_SRC1(HNEXT,I,INEW,INDEXO,TVAL)

    ! The purpose of DDE_SRC1 is to perform linear searches
    ! of the solution history queue for DDE_SOLVER.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: HNEXT, TVAL
    INTEGER, INTENT (IN) :: I, INEW
    INTEGER, INTENT (OUT) :: INDEXO
    ! ..
    ! .. Local Scalars ..
    INTEGER :: ISTART, J
    ! ..
    ! Note: The queue is ordered linearly in this version.
    ! The bracketing search begins here.
    ! Linear search; remember last bracketing interval.
    ! Start the search at the interval last located.
    ISTART = STARTS(I)

    IF (HNEXT>0.0D0) THEN
       ! TQUEUE is increasing.
       IF (TVAL<=TQUEUE(ISTART)) THEN
          DO J = ISTART, 0, -1
             IF (TVAL>=TQUEUE(J)) THEN
                ! TQUEUE(J)<=TVAL<=TQUEUE(J+1)
                ISTART = J
                INDEXO = J + 1
                GOTO 10
             END IF
          END DO
       ELSE
          ! TVAL>TQUEUE(ISTART)
          DO J = ISTART + 1, INEW
             IF (TVAL<=TQUEUE(J)) THEN
                ! TQUEUE(J-1)<=TVAL<=TQUEUE(J)
                ISTART = J - 1
                INDEXO = J
                GOTO 10
             END IF
          END DO
       END IF
    ELSE
       ! TQUEUE is decreasing.
       IF (TVAL>=TQUEUE(ISTART)) THEN
          DO J = ISTART, 0, -1
             IF (TVAL<=TQUEUE(J)) THEN
                ! TQUEUE(J)>=TVAL>=TQUEUE(J+1)
                ISTART = J
                INDEXO = J + 1
                GOTO 10
             END IF
          END DO
       ELSE
          ! TVAL<TQUEUE(ISTART)
          DO J = ISTART + 1, INEW
             IF (TVAL>=TQUEUE(J)) THEN
                ! TQUEUE(J-1)>=TVAL>=TQUEUE(J)
                ISTART = J - 1
                INDEXO = J
                GOTO 10
             END IF
          END DO
       END IF
    END IF

    ! The bracketing search ends here.
10  CONTINUE

    ! Save the search interval index for the next call to DDE_SRC1.
    STARTS(I) = ISTART
    IF (INFO_DROPPED) THEN
       IF (ISTART<1) THEN
          PRINT *, ' One of the delays is largerr than MAX_DELAY.'
          STOP
       END IF
    END IF

    RETURN
  END SUBROUTINE DDE_SRC1
  !____________________________________________________________________________
  SUBROUTINE DDE_SRC2(TVAL,NUMPT,IOLD,INEW,HNEXT,INDEXO,TOLD,TNEW,INDEX,I)

    ! The purpose of DDE_SRC2 is to assist in the queue search and
    ! interpolation for DDE_SOLVER.

    ! Note:
    ! This routine is not intended to be called with HNEXT=0.0D0
    ! or INEW=0. (The latter case would correspond to the first
    ! step of the integration and is handled separately in
    ! DDE_ZSET. The former should never occur and should be
    ! treated as an error.)

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: HNEXT, TVAL
    DOUBLE PRECISION, INTENT (OUT) :: TNEW, TOLD
    INTEGER, INTENT (IN) :: I, INEW, IOLD, NUMPT
    INTEGER, INTENT (OUT) :: INDEX, INDEXO
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    ! ..
    INDEX = 2

    ! INEW points to the most recent addition to the queue.
    ! IOLD points to the oldest.

    ! Check for extrapolation outside the queue. Also handle
    ! the case in which TVAL = TNEW.
    IF (HNEXT>0.0D0) THEN
       IF ((TVAL-TQUEUE(IOLD-1))<(-1.0D3*UROUND*MAX(ABS(TVAL), &
            ABS(TQUEUE(IOLD-1))))) THEN
          PRINT *, ' The integration was terminated to avoid extrapolating'
          PRINT *, ' past the solution history queue.'
          INDEX = -12
          GOTO 10
       END IF
       IF (TVAL>=TQUEUE(INEW)) THEN
          IF ((TVAL-TQUEUE(INEW))>(1.0D3*UROUND*MAX(ABS(TVAL), &
               ABS(TQUEUE(INEW))))) THEN
             INDEXO = INEW
             TOLD = TQUEUE(INEW-1)
             TNEW = TQUEUE(INEW)
             INDEX = -11
             GOTO 10
          END IF
       END IF
    END IF

    IF (HNEXT<0.0D0) THEN
       IF ((TVAL-TQUEUE(IOLD-1))>(1.0D3*UROUND*MAX(ABS(TVAL), &
          ABS(TQUEUE(IOLD-1))))) THEN
          PRINT *, ' The integration was terminated to avoid extrapolating'
          PRINT *, ' past the solution history queue.'
          INDEX = -12
          GOTO 10
       END IF
       IF (TVAL<=TQUEUE(INEW)) THEN
          IF ((TVAL-TQUEUE(INEW))<(-1.0D3*UROUND*MAX(ABS(TVAL), &
             ABS(TQUEUE(INEW))))) THEN
             INDEXO = INEW
             TOLD = TQUEUE(INEW-1)
             TNEW = TQUEUE(INEW)
             INDEX = -11
             GOTO 10
          END IF
       END IF
    END IF

    ! Handle the case in which there is only one point in the queue.
    IF ((NUMPT==1) .OR. ((IOLD==1) .AND. (INEW==1))) THEN
       INDEXO = 1
       TOLD = TQUEUE(0)
       TNEW = TQUEUE(1)
       GOTO 10
    END IF

    ! We get here only if there is more than one point in the
    ! queue and TQUEUE(IOLD-1)<=TVAL<TQUEUE(1,INEW-1).

    ! Handle the case in which 1<=IOLD<INEW (so that the
    ! data is already ordered). Note that the queue may not
    ! yet be full.
    CALL DDE_SRC1(HNEXT,I,INEW,INDEXO,TVAL)
    TOLD = TQUEUE(INDEXO-1)
    TNEW = TQUEUE(INDEXO)

10  CONTINUE

    RETURN
  END SUBROUTINE DDE_SRC2
  !____________________________________________________________________________

  SUBROUTINE DDE_ERR1(INDEX,RATIO)

    ! The purpose of DDE_ERR1 is to calculate error estimates
    ! for DDE_SOLVER. DDE_ERR1 is called by DDE_DRV3.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (OUT) :: RATIO
    INTEGER, INTENT (OUT) :: INDEX
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: ALPHA, TR, W9, XDENOM, XNUMER
    INTEGER :: I, METH5
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, MAX
    ! ..
    INDEX = 2

    ! Use the difference between the fifth and sixth order
    ! methods to estimate the error. The estimate is
    ! (Y(6,8) - Y(5,10))(C) with C = 1.
    ALPHA = -184.0D0/3375.0D0
    MYR(1:MYN) = ALPHA*KMAT(1:MYN,1)
    ALPHA = 564.0D0/3375.0D0
    MYR(1:MYN) = MYR(1:MYN) + ALPHA*KMAT(1:MYN,3)
    ALPHA = -60.0D0/3375.0D0
    MYR(1:MYN) = MYR(1:MYN) + ALPHA*KMAT(1:MYN,4)
    ALPHA = -1720.0D0/3375.0D0
    MYR(1:MYN) = MYR(1:MYN) + ALPHA*KMAT(1:MYN,5)
    ALPHA = -735.0D0/3375.0D0
    MYR(1:MYN) = MYR(1:MYN) + ALPHA*KMAT(1:MYN,6)
    ALPHA = -1596.0D0/3375.0D0
    MYR(1:MYN) = MYR(1:MYN) + ALPHA*KMAT(1:MYN,7)
    ALPHA = -644.0D0/3375.0D0
    MYR(1:MYN) = MYR(1:MYN) + ALPHA*KMAT(1:MYN,8)
    ALPHA = 1000.0D0/3375.0D0
    MYR(1:MYN) = MYR(1:MYN) + ALPHA*KMAT(1:MYN,9)
    ALPHA = 3375.0D0/3375.0D0
    MYR(1:MYN) = MYR(1:MYN) + ALPHA*KMAT(1:MYN,10)
    ! Multiply the contents of R by
    ! (Y(6,8) - Y(5,10))(1) = -Y(5,10)(1).
    IF (IPAIR==1) METH5 = 17
    IF (IPAIR==2) METH5 = 25
    IF (IPAIR==3) METH5 = 27
    IF (METH5==17) THEN
       W9 = 0.25D0
    ELSE
       IF (METH5==18) THEN
          W9 = -1.0D0
       ELSE
          IF (METH5==25) THEN
             W9 = 2.0D0
          ELSE
             IF (METH5==26) THEN
                W9 = -2.0D0
             ELSE
                ! Can make another substitution at this point.
                W9 = 1.0D0
             END IF
          END IF
       END IF
    END IF
    MYR(1:MYN) = W9*MYR(1:MYN)

    ! Calculate the error ratio (maximum estimated error per
    ! step divided by the weight).
    RATIO = 0.0D0
    DO I = 1, MYN
       XNUMER = ABS(MYR(I))
       XDENOM = MYRELER(I)*ABS(MYYNEW(I)) + MYABSER(I)
       TR = XNUMER/XDENOM
       RATIO = MAX(RATIO,TR)
    END DO
    !   DO I = 1, MYN
    !     XNUMER = ABS(MYR(I))
    !     XDENOM = MYRELER(I)*ABS(MYYNEW(I)) + MYABSER(I)
    !     TR = XNUMER/XDENOM
    !     RATIO = MAX(RATIO,TR)
    !   END DO

    ! The following value of TOLFAC will decrease the next
    ! estimated step size by a factor of(1/15)**(1/6) and
    ! therefore may increase the number of steps by about 60%.
    RATIO = TOLFAC*RATIO

    RETURN
  END SUBROUTINE DDE_ERR1
  !____________________________________________________________________________
  SUBROUTINE DDE_POLY(TOLD,POLD,H,K0,K2,K3,K4,K5,K6,K7,K8,K9,PVAL,DPVAL, &
       TVAL,ITASK)

    ! The purpose of DDE_POLY is to evaluate the solution polynomial
    ! and/or the derivative polynomial for DDE_SOLVER. DDE_POLY is
    ! the core interpolation subroutine.

    !     Parameters:
    !     TOLD     = Old integration time.
    !     POLD(*)  = Y(T0).
    !     H        = Step size.
    !     K0, K2,  = RKS derivative stage approximations.
    !     K3, K4,
    !     K5, K6,
    !     K7, K8,
    !     K9
    !     PVAL(*)  = Array of length N. If requested, the
    !                interpolated solution will be returned
    !                in this array. Not used if ITASK = 2,5.
    !     DPVAL(*) = Array of length N. If requested, the
    !                interpolated derivative will be returned
    !                in this array. Not used if ITASK = 1,4.
    !     TVAL     = Time at which the solution and/or
    !                derivative is to be approximated.
    !     ITASK    = Task flag, with the following meanings:
    !              =  1 Return the interpolated solution in
    !                   PVAL.
    !              =  2 return the interpolated derivative
    !                   IN DPVAL.
    !              =  3 Return the interpolated solution in
    !                   PVAL and return the interpolated
    !                   derivative in DPVAL.
    !              =  4 Add the input values of PVAL to the
    !                   interpolated solution and return
    !                   the results in PVAL.
    !              =  5 Add the input values of DPVAL
    !                   to the interpolated derivative
    !                   and return the results in DPVAL.
    !              =  6 Add the input values of PVAL to the
    !                   interpolated solution and add the
    !                   input values of DPVAL to the
    !                   interpolated derivative; and return
    !                   the results in PVAL and DPVAL,
    !                   respectively.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: H, TOLD, TVAL
    INTEGER, INTENT (IN) :: ITASK
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (OUT) :: DPVAL(:), PVAL(:)
    DOUBLE PRECISION, INTENT (IN) :: K0(:), K2(:), K3(:), K4(:), K5(:), &
         K6(:), K7(:), K8(:), K9(:), POLD(:)
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: C
    LOGICAL :: POLYD, POLYS
    ! ..
    ! .. Local Arrays ..
    DOUBLE PRECISION :: W(0:9), WD(0:9)
    ! ..
    ! Check for illegal method flags and zero step size.
    IF (ITASK<1 .OR. ITASK>3) THEN
       PRINT *, ' An illegal value of ITASK was encountered in DDE_POLY.'
       STOP
    END IF
    IF (IPAIR<1 .OR. IPAIR>3) THEN
       PRINT *, ' An illegal value of IPAIR was encountered in DDE_POLY.'
       STOP
    END IF
    ! IF (H==0.0D0) THEN
    !     PRINT *, ' The step size cannot be zero in DDE_POLY.'
    !     STOP
    ! END IF

    ! Determine what is to be done.
    IF (ITASK==1) THEN
       POLYS = .TRUE.
       POLYD = .FALSE.
    END IF
    IF (ITASK==2) THEN
       POLYS = .FALSE.
       POLYD = .TRUE.
    END IF
    IF (ITASK==3) THEN
       POLYS = .TRUE.
       POLYD = .TRUE.
    END IF

    ! Determine the value of C at which the solution or derivative
    ! is to be interpolated.
    IF (ABS(H)<=0.0D0) THEN
       C = 0.0D0
    ELSE
       C = (TVAL-TOLD)/H
    END IF

    ! Calculate the interpolation coefficients and/or their
    ! derivatives.
    CALL DDE_WSET(C,W,WD,POLYS,POLYD,MYMETHOD)

    IF (POLYS) THEN
       ! Evaluate the interpolation polynomial.
       PVAL(1:MYN) = POLD(1:MYN) + H*(W(0)*K0(1:MYN)+W(2)*K2(1:MYN)+W(3)*K3 &
            (1:MYN)+W(4)*K4(1:MYN)+W(5)*K5(1:MYN)+W(6)*K6(1:MYN)+ &
            W(7)*K7(1:MYN)+W(8)*K8(1:MYN)+W(9)*K9(1:MYN))
    END IF

    IF (POLYD) THEN
       ! Evaluate the derivative of the interpolation polynomial derivative.
       DPVAL(1:MYN) = WD(0)*K0(1:MYN) + WD(2)*K2(1:MYN) + WD(3)*K3(1:MYN) + &
            WD(4)*K4(1:MYN) + WD(5)*K5(1:MYN) + WD(6)*K6(1:MYN) + &
            WD(7)*K7(1:MYN) + WD(8)*K8(1:MYN) + WD(9)*K9(1:MYN)
    END IF

    RETURN
  END SUBROUTINE DDE_POLY
  !____________________________________________________________________________

  SUBROUTINE DDE_WSET(C,W,WD,POLYS,POLYD,METHOD)

    ! The purpose of DDE_WSET is to calculate the interpolation
    ! coefficients and/or their derivatives for DDE_SOLVER. It
    ! is called by the core interpolation subroutine DDE_POLY.

    !     PARAMETERS:
    !     C      = Value for which the coefficients will be formed.
    !     W(*)   = Interpolation coefficients. Dimension W(0:9).
    !     WD(*)  = Derivative coefficients. Dimension WD(0:9).
    !     POLYS  = Logical variable. If .TRUE., the interpolation
    !              coefficients will be formed.
    !     POLYD  = Logical variable. If .TRUE., the derivative
    !              coefficients will be formed.
    !     METHOD = Which method is to be used.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: C
    INTEGER, INTENT (IN) :: METHOD
    LOGICAL, INTENT (IN) :: POLYD, POLYS
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (INOUT) :: W(0:9), WD(0:9)
    ! ..
    ! SOLUTION.
    IF (POLYS) THEN

       ! Calculate W_9.

       ! Primary integration methods.

       IF (IORDER==6) THEN
          IF (METHOD==4) THEN
             W(9) = 0.0D0
          ELSE
             IF (METHOD==5) THEN
                W(9) = 0.75D0*(C*C*(22.0D0+C*(-137.0D0+ &
                     C*(298.0D0+C*(-273.0D0+C*(90.0D0))))))
             ELSE
                IF (METHOD==6) THEN
                   W(9) = (27.0D0/4.0D0)*(C*C*(-2.0D0+C*(7.0D0+ &
                        C*(-8.0D0+C*(3.0D0)))))
                ELSE
                   IF (METHOD==7) THEN
                      W(9) = (9.0D0/4.0D0)*(C*C*(-4.0D0+C*(11.0D0+ &
                           C*(-10.0D0+C*(3.0D0)))))
                   ELSE
                      ! Can make substitution for another
                      ! method at this point.
                      W(9) = 0.0D0
                   END IF
                END IF
             END IF
          END IF
       END IF

       ! Subsidiary integration methods.

       IF (IORDER==5) THEN
          IF (METHOD==17) THEN
             W(9) = 0.25D0*C**5
          ELSE
             IF (METHOD==18) THEN
                W(9) = -(C**5)
             ELSE
                IF (METHOD==25) THEN
                   W(9) = 2.0D0*(C**5)
                ELSE
                   IF (METHOD==26) THEN
                      W(9) = -2.0D0*(C**5)
                   ELSE
                      ! Can make another substitution at this point.
                      W(9) = C**5
                   END IF
                END IF
             END IF
          END IF
       END IF

       ! Define the remaining interpolation coefficients.

       ! Calculate W_0.
       W(0) = (C*(1.0D0+C*(-2427.0D0/500.0D0+C*(3586.0D0/375.0D0+C*(- &
            1659.0D0/200.0D0+C*(66.0D0/25.0D0)))))) - (184.0D0/3375.0D0)*W(9)

       ! Calculate W_1.
       W(1) = 0.0D0

       ! Calculate W_2.
       W(2) = (C*C*(708.0D0/125.0D0+C*(-1977.0D0/125.0D0+ &
            C*(789.0D0/50.0D0+C*(-27.0D0/5.0D0))))) + (188.0D0/1125.0D0)*W(9)

       ! Calculate W_3.
       W(3) = (C*C*(87.0D0/50.0D0+C*(-243.0D0/50.0D0+ &
            C*(201.0D0/40.0D0+C*(-9.0D0/5.0D0))))) - (4.0D0/225.0D0)*W(9)

       ! Calculate W_4.
       W(4) = (C*C*(-88.0D0/25.0D0+C*(1226.0D0/75.0D0+ &
            C*(-21.0D0+C*(42.0D0/5.0D0))))) - (344.0D0/675.0D0)*W(9)

       ! Calculate W_5.
       W(5) = (C*C*(-21.0D0/100.0D0+C*(21.0D0/25.0D0+ &
            C*(-21.0D0/40.0D0)))) - (49.0D0/225.0D0)*W(9)

       ! Calculate W_6.
       W(6) = (C*C*(228.0D0/125.0D0+C*(-1197.0D0/125.0D0+ &
            C*(741.0D0/50.0D0+C*(-171.0D0/25.0D0))))) - &
            (532.0D0/1125.0D0)*W(9)

       ! Calculate W_7.
       W(7) = (C*C*(-161.0D0/250.0D0+C*(1127.0D0/750.0D0+ &
            C*(-161.0D0/200.0D0)))) - (644.0D0/3375.0D0)*W(9)

       ! Calculate W_8.
       W(8) = (C*C*C*(2.0D0+C*(-5.0D0+C*(3.0D0)))) + (8.0D0/27.0D0)*W(9)

    END IF

    ! Derivative.

    IF (POLYD) THEN

       ! Calculate WD_9.

       ! Primary integration methods.

       IF (IORDER==6) THEN
          IF (METHOD==4) THEN
             WD(9) = 0.0D0
          ELSE
             IF (METHOD==5) THEN
                WD(9) = 0.75D0*(C*(44.0D0+C*(-411.0D0+ &
                     C*(1192.0D0+C*(-1365.0D0+C*(540.0D0))))))
             ELSE
                IF (METHOD==6) THEN
                   WD(9) = (27.0D0/4.0D0)*(C*(-4.0D0+C*(21.0D0+ &
                        C*(-32.0D0+C*(15.0D0)))))
                ELSE
                   IF (METHOD==7) THEN
                      WD(9) = (9.0D0/4.0D0)*(C*(-8.0D0+C*(33.0D0+ &
                           C*(-40.0D0+C*(15.0D0)))))
                   ELSE
                      ! Can make substitution for another
                      ! method at this point.
                      WD(9) = 0.0D0
                   END IF
                END IF
             END IF
          END IF
       END IF

       ! Subsidiary integration methods.

       IF (IORDER==5) THEN
          IF (METHOD==17) THEN
             WD(9) = 1.25D0*(C**4)
          ELSE
             IF (METHOD==18) THEN
                WD(9) = -5.0D0*(C**4)
             ELSE
                IF (METHOD==25) THEN
                   WD(9) = 10.0D0*(C**4)
                ELSE
                   IF (METHOD==26) THEN
                      WD(9) = -10.0D0*(C**4)
                   ELSE
                      ! Can make another substitution at this point.
                      WD(9) = 5.0D0*(C**4)
                   END IF
                END IF
             END IF
          END IF
       END IF

       ! Define the remaining derivative coefficients.

       ! Calculate WD_0.
       WD(0) = (1.0D0+C*(-2427.0D0/250.0D0+C*(3586.0D0/125.0D0+ &
            C*(-1659.0D0/50.0D0+C*(66.0D0/5.0D0))))) - &
            (184.0D0/3375.0D0)*WD(9)

       ! Calculate WD_1.
       WD(1) = 0.0D0

       ! Calculate WD_2.
       WD(2) = (C*(1416.0D0/125.0D0+C*(-5931.0D0/125.0D0+ &
            C*(1578.0D0/25.0D0+C*(-27.0D0))))) + (188.0D0/1125.0D0)*WD(9)

       ! Calculate WD_3.
       WD(3) = (C*(87.0D0/25.0D0+C*(-729.0D0/50.0D0+ &
            C*(201.0D0/10.0D0+C*(-9.0D0))))) - (4.0D0/225.0D0)*WD(9)

       ! Calculate WD_4.
       WD(4) = (C*(-176.0D0/25.0D0+C*(1226.0D0/25.0D0+ &
            C*(-84.0D0+C*(42.0D0))))) - (344.0D0/675.0D0)*WD(9)

       ! Calculate WD_5.
       WD(5) = (C*(-21.0D0/50.0D0+C*(63.0D0/25.0D0+C*(-21.0D0/10.0D0)))) - &
            (49.0D0/225.0D0)*WD(9)

       ! Calculate WD_6.
       WD(6) = (C*(456.0D0/125.0D0+C*(-3591.0D0/125.0D0+ &
            C*(1482.0D0/25.0D0+C*(-171.0D0/5.0D0))))) - &
            (532.0D0/1125.0D0)*WD(9)

       ! Calculate WD_7.
       WD(7) = (C*(-161.0D0/125.0D0+C*(1127.0D0/250.0D0+ &
            C*(-161.0D0/50.0D0)))) - (644.0D0/3375.0D0)*WD(9)

       ! Calculate WD_8.
       WD(8) = (C*C*(6.0D0+C*(-20.0D0+C*(15.0D0)))) + (8.0D0/27.0D0)*WD(9)

    END IF

    RETURN
  END SUBROUTINE DDE_WSET
  !____________________________________________________________________________

  SUBROUTINE DDE_DERV(T,YSOL,DY,DERIVS)

    ! DDE_DERV is an internal subroutine to calculate derivatives.

    ! .. Scalar Arguments ..
    REAL (KIND(1D0)), INTENT (IN) :: T
    ! ..
    ! .. Array Arguments ..
    REAL (KIND(1D0)), INTENT (INOUT) :: DY(:)
    REAL (KIND(1D0)), INTENT (IN) :: YSOL(:)
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL DERIVS
    INTERFACE
       SUBROUTINE DERIVS(T,Y,Z,DY)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         INTENT(IN):: T,Y,Z
         INTENT(OUT) :: DY
       END SUBROUTINE DERIVS
    END INTERFACE
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC KIND
    ! ..
    IF (MYNEUTRAL) THEN
        IF (NLAGS>0) THEN
            ZANDD(1:MYN,1:NLAGS) = ZARRAY(1:MYN,1:NLAGS)
            ZANDD(1:MYN,NLAGS+1:2*NLAGS) = DARRAY(1:MYN,1:NLAGS)
        END IF
        CALL DERIVS(T,YSOL,ZANDD,DY)
    ELSE
        CALL DERIVS(T,YSOL,ZARRAY,DY)
    END IF

    RETURN
  END SUBROUTINE DDE_DERV
  !____________________________________________________________________________

  SUBROUTINE DDE_BETA(T,Y,BETA)

    ! DDE_BETA is an internal subroutine to calculate delays.

    ! .. Scalar Arguments ..
    INTEGER I
    DOUBLE PRECISION, INTENT (IN) :: T
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: Y(:)
    ! ..
    ! .. Subroutine Arguments ..
    INTERFACE
       SUBROUTINE BETA(T, Y, BVAL)
          DOUBLE PRECISION, INTENT(IN) :: T
          DOUBLE PRECISION, INTENT(IN) :: Y(:)
          DOUBLE PRECISION, INTENT(OUT) :: BVAL(:)
       END SUBROUTINE BETA
    END INTERFACE
    ! ..
    ! Do nothing if NLAGS = 0.
    IF (NLAGS>0) THEN
        IF (CONSTANT_DELAYS) THEN
            CALL CD_BETA(T,MYBVAL,NLAGS)
        ELSE
            CALL BETA(T,Y,MYBVAL)
        END IF
    END IF

    IF (INFO_DROPPED) THEN
       DO I = 1, NLAGS
          IF (MYBVAL(I)>TQUEUE(0) .AND. MYBVAL(I)<TQUEUE(1)) THEN
             PRINT *, ' A delay is larger than the value of MAX_DELAY'
             PRINT *, ' you specified. Stopping.'
             STOP
          END IF
       END DO
    END IF

    RETURN
  END SUBROUTINE DDE_BETA
  !____________________________________________________________________________

  SUBROUTINE DDE_YINIT(T,YINIT)

    ! DDE_YINIT is an internal subroutine to calculate the initial history.

    ! VTEMP3(1:MYN) is the initial solution. If this is a neutral
    ! problem, VTEMP3(MYN+1:2*N) is the initial derivative.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: T
    ! ..
    ! .. Subroutine Arguments ..
    INTERFACE
       SUBROUTINE YINIT(T,Y)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y
          INTENT(IN):: T
          INTENT(OUT) :: Y
       END SUBROUTINE YINIT
    END INTERFACE
    ! ..
    IF (CONSTANT_HISTORY) THEN
        CALL CH_YINIT(VTEMP3(1:MYN),MYN)
        IF (MYNEUTRAL) VTEMP3(MYN+1:2*MYN) = 0D0            
    ELSE
        CALL YINIT(T,VTEMP3)
    END IF

    RETURN
  END SUBROUTINE DDE_YINIT
  !____________________________________________________________________________

  SUBROUTINE DDE_ZSET(T,Y,HNEXT,YINIT,INDEX,IMETH,K0,K2,K3,K4,K5,K6,K7,K8, &
       K9,TOLD,YOLD,ITERATE)

    ! The purpose of DDE_ZSET is to interpolate the past solution
    ! history queue for DDE_SOLVER. DDE_ZSET calculates
    ! ZARRAY and DARRAY if MYNEUTRAL is .TRUE. ZARRAY(J,I) is the
    ! approximate solution for component J evaluated at the Ith
    ! delay time, MYBVAL(I). DARRAY(J,I) is the approximate derivative
    ! for component J evaluated at the Ith delay time, MYBVAL(I).

    !     The following parameters must be input to DDE_ZSET.
    !     T       = Current value of the independent variable.
    !     Y       = Solution at time T.
    !     HNEXT   =  1.0D0 if the integration goes left to right.
    !             = -1.0D0 if the integration goes right to left.
    !     YINIT   = External function to calculate the initial
    !               delayed solution.
    !     IMETH   = Flag to indicate whether to extrapolate the old
    !               solution polynomials or perform euler extrapolations
    !               (near the initial time only).
    !     K0,K2,K3= RKS derivative approximations for the integration
    !     K4,K5,K6, interval beginning at TOLD.
    !     K7,K8,K9
    !     TOLD    = Last successful integration time.
    !     YOLD    = Solution at TOLD.
    !     The following will be output by DDE_ZSET.
    !     ITERATE     = .TRUE. if subroutine DDE_ZSET3 is called
    !     INDEX    - =   2 if DDE_ZSET was successful.
    !                   -4 if HNEXT = 0.
    !                  -11 if an error occurred in DDE_SRC2.
    !                  -12 if an error occurred in DDE_SRC2.
    !     The private array, ZARRAY(N,NLAGS), will be calculated, as
    !     DDE_ZSET assumes BETA was called before it is called.

    !     Interpolation subroutine call information.

    !               DDE_DRV2 calls
    !              /      \
    !          DDE_INTP    DDE_ZSET
    !                      /   \
    !               DDE_ZSET2   DDE_ZSET3
    !                  |
    !               DDE_ZSET3

    !               DDE_STEP calls
    !                   |
    !                DDE_ZSET
    !               /        \
    !        DDE_ZSET2        DDE_ZSET3
    !           |
    !        DDE_ZSET3
    !               DDE_HST1 calls
    !                  |
    !               DDE_HST2
    !                  |
    !               DDE_ZSET3

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: HNEXT, T, TOLD
    INTEGER, INTENT (IN) :: IMETH
    INTEGER, INTENT (OUT) :: INDEX
    LOGICAL, INTENT (INOUT) :: ITERATE
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: K0(:), K2(:), K3(:), K4(:), K5(:), &
         K6(:), K7(:), K8(:), K9(:), Y(:), YOLD(:)
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL YINIT
    INTERFACE
       SUBROUTINE YINIT(T,Y)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y
         INTENT(IN):: T
         INTENT(OUT) :: Y
      END SUBROUTINE YINIT
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: TVAL
    INTEGER :: I, IQUEUE
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS
    ! ..
    INDEX = 2

    ! Do nothing if NLAGS = 0.
    IF (NLAGS<=0) RETURN

    IF (ABS(HNEXT)<=0.0D0) THEN
       PRINT *, ' The required step size is too small.'
       INDEX = -4
       GOTO 20
    END IF

    ! IQUEUE Points to the most recent addition to the queue.
    IQUEUE = MYIPOINT(2)

    ! If the first integration step has not been completed, use the
    ! initial delay interval information using euler extrapolations
    ! if necessary.
    IF (IQUEUE==0) THEN

       DO 10 I = 1, NLAGS

          ! Define the Ith delay time.
          TVAL = MYBVAL(I)

          ! IF (TVAL==T) THEN
          IF (ABS(TVAL-T)<=0.0D0) THEN
             ! Extrapolate for MYDYEXTP; do not interpolate for MYYEXTP.
             IF (MYNEUTRAL) THEN
                CALL DDE_ZSET3(TOLD,YOLD,HNEXT,K0,K2,K3,K4,K5,K6,K7,K8,K9,I, &
                     IMETH)
                DARRAY(1:MYN,I) = MYDYEXTP(1:MYN)
                ITERATE = .TRUE.
                ! Do not interpolate for MYYEXTP.
                ZARRAY(1:MYN,I) = Y(1:MYN)
             ELSE
                ! Do not interpolate for MYYEXTP.
                ZARRAY(1:MYN,I) = Y(1:MYN)
             END IF
             GOTO 10
          END IF

          ! If the delay does not fall in the initial interval and is
          ! not equal to current time T, extrapolate.
          IF ((TVAL>MYTINIT) .AND. (HNEXT>0.0D0) .OR. &
              (TVAL<MYTINIT) .AND. (HNEXT<0.0D0)) THEN
             CALL DDE_ZSET3(TOLD,YOLD,HNEXT,K0,K2,K3,K4,K5,K6,K7,K8,K9,I, &
                  IMETH)
             ITERATE = .TRUE.
             ZARRAY(1:MYN,I) = MYYEXTP(1:MYN)
             IF (MYNEUTRAL) DARRAY(1:MYN,I) = MYDYEXTP(1:MYN)
          ELSE
             ! Otherwise use the history function.
             CALL DDE_YINIT(TVAL,YINIT)
             ZARRAY(1:MYN,I) = VTEMP3(1:MYN)
             IF (MYNEUTRAL) THEN
                DARRAY(1:MYN,I) = VTEMP3(MYN+1:2*MYN)
             END IF
          END IF
10     END DO

       GOTO 20

    END IF

    ! Else, at least one integration step has been completed
    ! (IQUEUE>0) so interpolate using the history queue.
    CALL DDE_ZSET2(T,Y,YINIT,HNEXT,INDEX,TOLD,YOLD,K0,K2,K3,K4,K5,K6,K7, &
         K8,K9,IMETH,ITERATE)

20  CONTINUE

    RETURN
  END SUBROUTINE DDE_ZSET
  !____________________________________________________________________________
  SUBROUTINE DDE_ZSET2(T,Y,YINIT,HNEXT,INDEX,TOLD,YOLD,K0,K2,K3,K4,K5,K6, &
       K7,K8,K9,IMETH,ITERATE)

    ! The purpose of DDE_ZSET2 is to search and interpolate the past
    ! solution history queue for DDE_SOLVER. DDE_ZSET2 is called only
    ! by DDE_ZSET. If necessary, it calls DDE_ZSET3 to calculate Euler
    ! extrapolations. DDE_ZSET2 calculates ZARRAY and DARRAY if
    ! MYNEUTRAL is true. ZARRAY(J,I) is the approximate solution for
    ! component J evaluated at the Ith delay time, MYBVAL(I).
    ! DARRAY(J,I) is the approximate derivative for component J
    ! evaluated at the Ith delay time, MYBVAL(I).

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: HNEXT, T, TOLD
    INTEGER, INTENT (IN) :: IMETH
    INTEGER, INTENT (OUT) :: INDEX
    LOGICAL, INTENT (INOUT) :: ITERATE
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: K0(:), K2(:), K3(:), K4(:), K5(:), &
         K6(:), K7(:), K8(:), K9(:), Y(:), YOLD(:)
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL YINIT
    INTERFACE
       SUBROUTINE YINIT(T,Y)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y
         INTENT(IN):: T
         INTENT(OUT) :: Y
      END SUBROUTINE YINIT
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: DELINT, TN, TO, TVAL
    INTEGER :: I, IBEGIN, INDEXO, INEW, IOLD, IQUEUE, ITASK, JQUEUE, NUMPT
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS
    ! ..
    INDEX = 2

    ! Do nothing if NLAGS = 0.
    IF (NLAGS<=0) RETURN

    ! IQUEUE points to the most recent addition to the queue.
    ! JQUEUE points to the oldest.
    IQUEUE = MYIPOINT(2)
    JQUEUE = MYIPOINT(3)
    ! NUMPT is the maximum number of points for which history data
    ! can be stored.
    NUMPT = LQUEUE/(10*MYN)

    DO 20 I = 1, NLAGS

       ! Define the Ith delay time.
       TVAL = MYBVAL(I)
       ! IF (TVAL==T) THEN
       IF (ABS(TVAL-T)<=0.0D0) THEN
          IF (MYNEUTRAL) THEN
             ! Extrapolate for MYDYEXTP.
             CALL DDE_ZSET3(TOLD,YOLD,HNEXT,K0,K2,K3,K4,K5,K6,K7,K8,K9,I, &
                  IMETH)
             DARRAY(1:MYN,I) = MYDYEXTP(1:MYN)
             ITERATE = .TRUE.
             ! Do not interpolate for MYYEXTP.
             ZARRAY(1:MYN,I) = Y(1:MYN)
          ELSE
             ZARRAY(1:MYN,I) = Y(1:MYN)
          END IF
          GOTO 20
       END IF

       ! Extrapolate if the step size exceeds the delay.
       IF ((TVAL>TOLD) .AND. (HNEXT>0.0D0) .OR. &
           (TVAL<TOLD) .AND. (HNEXT<0.0D0)) THEN
          CALL DDE_ZSET3(TOLD,YOLD,HNEXT,K0,K2,K3,K4,K5,K6,K7,K8,K9,I,IMETH)
          ITERATE = .TRUE.
          ZARRAY(1:MYN,I) = MYYEXTP(1:MYN)
          IF (MYNEUTRAL) DARRAY(1:MYN,I) = MYDYEXTP(1:MYN)
          GOTO 20
       END IF

       ! Interpolate using the solution history queue.
       IOLD = JQUEUE
       INEW = IQUEUE

       ! Use YINIT or call the search routine to locate the old slug
       ! of data to be interpolated for each value of MYBVAL(I).
       ! Use YINIT if TVAL falls in the initial delay interval.
       IF ((HNEXT>0.0D0) .AND. (TVAL>MYTINIT)) GOTO 10
       IF ((HNEXT<0.0D0) .AND. (TVAL<MYTINIT)) GOTO 10
       CALL DDE_YINIT(TVAL,YINIT)
       ZARRAY(1:MYN,I) = VTEMP3(1:MYN)
       IF (MYNEUTRAL) THEN
          DARRAY(1:MYN,I) = VTEMP3(MYN+1:2*MYN)
       END IF
       GOTO 20

       ! Otherwise use the history queue.
10     CONTINUE
       CALL DDE_SRC2(TVAL,NUMPT,IOLD,INEW,HNEXT,INDEXO,TO,TN,INDEX,I)

       ! Check for error return.
       IF (INDEX==-11) GOTO 30
       IF (INDEX==-12) GOTO 30

       ! The data to be interpolated corresponds to TO = TOLD,
       ! TN = TNEW. INDEXO is the pointer for the corresponding
       ! slug of data. The data for told to be interpolated
       ! begins in queue location (1,IBEGIN+1).
       IBEGIN = (INDEXO-1)*10

       ! Define the step size corresponding to the interpolation data.
       DELINT = TN - TO

       ! Interpolate the data.
       ITASK = 1
       IF (MYNEUTRAL) ITASK = 3
       CALL DDE_POLY(TO,QUEUE(:,IBEGIN+1),DELINT,QUEUE(:,IBEGIN+2), &
            QUEUE(:,IBEGIN+3),QUEUE(:,IBEGIN+4),QUEUE(:,IBEGIN+5), &
            QUEUE(:,IBEGIN+6),QUEUE(:,IBEGIN+7),QUEUE(:,IBEGIN+8), &
            QUEUE(:,IBEGIN+9),QUEUE(:,IBEGIN+10),MYR,VTEMP3(MYN+1:2*MYN), &
            TVAL,ITASK)
       ZARRAY(1:MYN,I) = MYR(1:MYN)
       IF (MYNEUTRAL) DARRAY(1:MYN,I) = VTEMP3(MYN+1:2*MYN)

20  END DO

30  CONTINUE

    RETURN
  END SUBROUTINE DDE_ZSET2
  !____________________________________________________________________________

  SUBROUTINE DDE_INTP(TOLD,TNEW,TOUT,YOUT,DYOUT)

    ! The purpose of DDE_INTP is to evaluate the solution polynomial
    ! and the derivative polynomial for DDE_SOLVER. It does not
    ! calculate ZARRAY or DARRAY.

    ! The input parameters TOLD, DYOLD, and TNEW are output from
    ! the last call to DDE_DRV3. TOUT is the value of the
    ! independent variable at which the solution and derivatives are
    ! to be interpolated. The interpolated solution is returned in
    ! YOUT. The interpolated derivatives are returned in DYOUT.
    ! NOTE:
    ! DDE_INTP is called only by DDE_DRV2 to obtain the interpolated
    ! solution and derivative to pass to OUT_FCN.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: TNEW, TOLD, TOUT
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (OUT) :: DYOUT(:), YOUT(:)
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: DELINT
    ! ..
    ! Define the step size that was used.
    DELINT = TNEW - TOLD

    ! Do the interpolation.
    CALL DDE_POLY(TOLD,MYYOLD,DELINT,KMAT(:,1),KMAT(:,3),KMAT(:,4), &
         KMAT(:,5),KMAT(:,6),KMAT(:,7),KMAT(:,8),KMAT(:,9),KMAT(:,10),YOUT, &
         DYOUT,TOUT,3)

    RETURN
  END SUBROUTINE DDE_INTP
  !____________________________________________________________________________

  SUBROUTINE DDE_ZSET3(TOLD,YOLD,H,K0,K2,K3,K4,K5,K6,K7,K8,K9,I,IMETH)

    ! The purpose of DDE_ZSET3 is to evaluate the solution polynomial
    ! and the derivative polynomial for dde_solver using Euler
    ! extrapolations near the initial time if necessary. It calculates
    ! ZARRAY and DARRAY if MYNEUTRAL is .TRUE. ZARRAY(J,I) is the
    ! approximate solution for component J evaluated at the Ith delay
    ! time, MYBVAL(I). DARRAY(J,I) is the approximate derivative for
    ! component J evaluated at the Ith delay time, MYBVAL(I).

    ! Note:
    ! DDE_ZSET3 is called by DDE_ZSET and by DDE_ZSET2.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: H, TOLD
    INTEGER, INTENT (IN) :: I, IMETH
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: K0(:), K2(:), K3(:), K4(:), K5(:), &
         K6(:), K7(:), K8(:), K9(:), YOLD(:)
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: TVAL, Z
    INTEGER :: ITASK
    ! ..
    ! Do nothing if NLAGS = 0.
    IF (NLAGS<=0) RETURN

    ! Define the Ith delay time.
    TVAL = MYBVAL(I)

    IF (IMETH==9) GOTO 10

    ! Extrapolate the solution polynomial.
    ITASK = 1
    IF (MYNEUTRAL) ITASK = 3
    CALL DDE_POLY(TOLD,YOLD,H,K0,K2,K3,K4,K5,K6,K7,K8,K9,MYYEXTP,MYDYEXTP, &
         TVAL,ITASK)
    ZARRAY(1:MYN,I) = MYYEXTP(1:MYN)
    IF (MYNEUTRAL) DARRAY(1:MYN,I) = MYDYEXTP(1:MYN)
    GOTO 20

10  CONTINUE

    ! Perform Euler extrapolations.
    Z = TVAL - TOLD
    ! MYYEXTP(1:MYN) = YOLD(1:MYN)
    MYYEXTP(1:MYN) = YOLD(1:MYN) + Z*K0(1:MYN)
    ZARRAY(1:MYN,I) = YOLD(1:MYN) + Z*K0(1:MYN)
    IF (MYNEUTRAL) THEN
       MYDYEXTP(1:MYN) = K0(1:MYN)
       DARRAY(1:MYN,I) = K0(1:MYN)
    END IF

20  CONTINUE

    RETURN
  END SUBROUTINE DDE_ZSET3
  !____________________________________________________________________________
  SUBROUTINE DDE_GRT2(HMIN,JFLAG,XA,XB,GX,X,ALPHA,X2,IMAX,LAST,IRCALL)

    ! The purpose of DDE_GRT2 is to perform root finding for DDE_SOLVER.
    ! DDE_GRT2 is a modification of the root finder used in the LSODAR
    ! ode solver. It is called by DDE_GRT5.

    ! This routine finds the leftmost root of a set of arbitrary
    ! event functions.

    ! Here the sign of XB - XA is arbitrary, but is constant for a given
    ! problem, and leftmost means nearest to XA. The values of the vector
    ! valued function G(X) = (G(I), I=1,...,NG) are communicated through
    ! the call sequence of DDE_GRT2. The method used is the Illinois
    ! algorithm.
    ! Note:
    ! This routine is a modification of the Petzold/Shampine root finder
    ! used in the LSODAR ode solver. One modification is the use of the
    ! added parameter IRCALL. DDE_GRT2 is called with the values 0, 1, 2
    ! corresponding respectively to extrapolatory root finding,
    ! interpolatory root finding on an interval (TOLD,TNEW), and
    ! interpolatory root finding on an interval (TROOT,TNEW) where TROOT
    ! is the root located when the IRCALL = 1 call is made. When called
    ! with IRCALL = 2, DDE_GRT2 does not alter the MYJROOT(*) array.

    ! Description of parameters.
    ! HMIN   = Resolution parameter in X. Input only. When a root is
    !          found, it is located only to within an error of HMIN in X.
    !          Typically, HMIN should be set to something on the order of
    !                10.0D0*UROUND*MAX(ABS(XA),ABS(XB)),
    !          where UROUND is the unit roundoff of the machine.
    ! JFLAG  = Integer flag for input and output communication.
    !          On input, set JFLAG = 0 on the first call for the problem,
    !          and leave it unchanged until the problem is completed.
    !          (The problem is completed when JFLAG>=2 on return.)
    !          On output, JFLAG has the following values and meanings:
    !          JFLAG = 1 means DDE_GRT2 needs a value of G(X). Set
    !                    GX = G(X) and call again.
    !          JFLAG = 2 means a root has been found. The root is at,
    !                    X and GX contains G(X). (Actually, X is the
    !                    rightmost approximation to the root on an
    !                    interval (XA,XB) of size HMIN or less.)
    !          JFLAG = 3 means X = XB is a root, with one or more
    !                    of the GI being zero at XB and no sign
    !                    changes in (XA,XB). GX contains G(X) on output.
    !          JFLAG = 4 means no roots (of odd multiplicity) were
    !                    found in (XA,XB) (no sign changes).
    ! XA,XB  = Endpoints of the interval where roots are sought.
    !          XB and XA are input when JFLAG = 0 (first call), and
    !          must be left unchanged between calls until the problem is
    !          completed. XA and XB must be distinct, but XB - XA may be
    !          of either sign. However, the notion of left and right
    !          will be used to mean nearer to XA or XB, respectively.
    !          When JFLAG>=2 on return, XA and XB are output, and
    !          are the endpoints of the relevant interval.
    ! GX     = Array of length NG containing G(X). GX is input
    !          when JFLAG = 1, and output when JFLAG>=2.
    ! X      = Independent variable value. Output only.
    !          When JFLAG = 1 on output, X is the point at which G(X)
    !          is to be evaluated and loaded into GX.
    !          When JFLAG = 2 or 3, X is the root.
    !          When JFLAG = 4, X is the right endpoint, XB, of the interval.
    ! MYJROOT= Integer array of length NG. Output only.
    !          When JFLAG = 2 or 3, MYJROOT indicates which components
    !          of G(X) have a root at X. MYJROOT(I) is 1 if the Ith
    !          component has a root, and MYJROOT(I) = 0 otherwise.

    !    Note: The GA and GB arrays are private (not arguments).

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (INOUT) :: ALPHA, X, X2, XA, XB
    DOUBLE PRECISION, INTENT (IN) :: HMIN
    INTEGER, INTENT (INOUT) :: IMAX, JFLAG, LAST
    INTEGER, INTENT (IN) :: IRCALL
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (INOUT) :: GX(:)
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: T2, TMAX
    INTEGER :: I, IMXOLD, NG, NXLAST
    LOGICAL :: BXROOT, BZROOT, SGNCHG
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, SIGN
    ! ..
    NG = MYIPOINT(6)

    IF (JFLAG==1) GOTO 90

    ! JFLAG /= 1. Check for change in sign of G or zero at XB.
    IMAX = 0
    TMAX = 0.0D0
    BZROOT = .FALSE.
    DO 20 I = 1, NG
       IF (ABS(GB(I))>0.0D0) GOTO 10
       BZROOT = .TRUE.
       GOTO 20

       ! At this point, GA(I) has been checked and cannot be zero.
10     IF (ABS(SIGN(1.0D0,GA(I))-SIGN(1.0D0,GB(I)))<=0.0D0) GOTO 20

       ! Do not process the root if MYDIRECTION says not to.
       IF (MYNGUSER>0) THEN
          IF (I<=MYNGUSER .AND. SIGN(1.0D0,GA(I))>0 .AND. MYDIRECTION(I)==1) &
               THEN
             ! MYREPORT(I) = 0
             GOTO 20
          END IF
          IF (I<=MYNGUSER .AND. SIGN(1.0D0,GA(I))<0 .AND. MYDIRECTION(I)==-1 &
               ) THEN
             ! MYREPORT(I) = 0
             GOTO 20
          END IF
       END IF

       T2 = ABS(GB(I)/(GB(I)-GA(I)))
       IF (T2<=TMAX) GOTO 20
       TMAX = T2
       IMAX = I
20  END DO
    IF (IMAX>0) GOTO 30
    SGNCHG = .FALSE.
    GOTO 40
30  SGNCHG = .TRUE.
40  IF (.NOT. SGNCHG) GOTO 200

    ! There is a sign change. Find the first root in the interval.
    BXROOT = .FALSE.
    LAST = 1
    NXLAST = 0

    ! Repeat until the first root in the interval is found. Loop point.
50  IF (BXROOT) GOTO 170
    IF (NXLAST==LAST) GOTO 60
    ALPHA = 1.0D0
    GOTO 80
60  IF (LAST==0) GOTO 70
    ALPHA = 0.5D0*ALPHA
    GOTO 80
70  ALPHA = 2.0D0*ALPHA
80  X2 = XB - (XB-XA)*GB(IMAX)/(GB(IMAX)-ALPHA*GA(IMAX))
    IF ((ABS(X2-XA)<HMIN) .AND. (ABS(XB-XA)>10.0D0*HMIN)) &
       X2 = XA + 0.1D0*(XB-XA)
    JFLAG = 1
    X = X2

    ! Return to the calling routine to get a value of GX = G(X).
    RETURN
90  IMXOLD = IMAX

    ! Check to see in which interval G changes sign.
    IMAX = 0
    TMAX = 0.0D0
    BZROOT = .FALSE.
    DO 110 I = 1, NG
       IF (ABS(GX(I))>0.0D0) GOTO 100
       BZROOT = .TRUE.
       GOTO 110

       ! Neither GA(I) nor GX(I) can be zero at this point.
100    IF (ABS(SIGN(1.0D0,GA(I))-SIGN(1.0D0,GX(I)))<=0.0D0) GOTO 110

       ! Do not process the root if MYDIRECTION says not to.
       IF (MYNGUSER>0) THEN
          IF (I<=MYNGUSER .AND. SIGN(1.0D0,GA(I))>0 .AND. MYDIRECTION(I)==1) &
               THEN
             ! MYREPORT(I) = 0
             GOTO 110
          END IF
          IF (I<=MYNGUSER .AND. SIGN(1.0D0,GA(I))<0 .AND. MYDIRECTION(I)==-1 &
               ) THEN
             ! MYREPORT(I) = 0
             GOTO 110
          END IF
       END IF

       T2 = ABS(GX(I)/(GX(I)-GA(I)))
       IF (T2<=TMAX) GOTO 110
       TMAX = T2
       IMAX = I
110 END DO
    IF (IMAX>0) GOTO 120
    SGNCHG = .FALSE.
    IMAX = IMXOLD
    GOTO 130
120 SGNCHG = .TRUE.
130 NXLAST = LAST
    IF (.NOT. SGNCHG) GOTO 140

    ! Sign change between XA and X2, so replace XB with X2.
    XB = X2
    GB(1:NG) = GX(1:NG)
    LAST = 1
    BXROOT = .FALSE.
    GOTO 160
140 IF (.NOT. BZROOT) GOTO 150

    ! Zero value at X2 and no sign change in (XA,X2), so X2 is a root.
    XB = X2
    GB(1:NG) = GX(1:NG)
    BXROOT = .TRUE.
    GOTO 160

    ! No sign change between XA and X2. Replace XA with X2.
150 CONTINUE
    GA(1:NG) = GX(1:NG)
    XA = X2
    LAST = 0
    BXROOT = .FALSE.
160 IF (ABS(XB-XA)<=HMIN) BXROOT = .TRUE.
    GOTO 50

    ! Return with XB as the root. Set MYJROOT. Set X = XB and GX = GB.
170 JFLAG = 2
    X = XB
    GX(1:NG) = GB(1:NG)
    DO 190 I = 1, NG
       IF (IRCALL/=2) MYJROOT(I) = 0
       IF (ABS(GB(I))>0.0D0) GOTO 180
       IF (IRCALL/=2) MYJROOT(I) = 1
       GOTO 190
180    IF ((ABS(SIGN(1.0D0,GA(I))-SIGN(1.0D0,GB(I)))>0.0D0) .AND. &
          (IRCALL/=2)) MYJROOT(I) = 1
190 END DO
    RETURN

    ! No sign change in the interval. Check for zero at right endpoint.
200 IF (.NOT. BZROOT) GOTO 210

    ! Zero value at XB and no sign change in (XA,XB). Return JFLAG = 3.
    X = XB
    GX(1:NG) = GB(1:NG)
    DO I = 1, NG
       IF (IRCALL/=2) MYJROOT(I) = 0
       IF (ABS(GB(I))<=0.0D0 .AND. IRCALL/=2) MYJROOT(I) = 1
    END DO
    JFLAG = 3
    RETURN

    ! No sign changes in this interval. Set X = XB. Return JFLAG = 4.
210 GX(1:NG) = GB(1:NG)
    X = XB
    JFLAG = 4

    RETURN
  END SUBROUTINE DDE_GRT2
  !____________________________________________________________________________

  SUBROUTINE DDE_GRT3(T,Y,DY,GUSER,G,BETA)

    ! The purpose of DDE_GRT3 is to evaluate the residuals for the
    ! root finder for DDE_SOLVER. It is called by subroutine
    ! DDE_GRT5 in conjunction with the DDE_GRT2 root finder.

    !     Input -
    !     T         = Independent variable (time).
    !     Y(*)      = Solution at time T.
    !     DY(*)     = Approximation to derivative at time T
    !                 (not exact derivative from DERIVS).
    !     GUSER     = User provided external subroutine to
    !                 evaluate the residuals for the user root
    !                 functions.
    !     BETA      = Name of subroutine in which to evaluate delays.
    !                 This subroutine must define the event function
    !                 residuals, G(I), I=1,...,NG.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: T
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: DY(:), Y(:)
    DOUBLE PRECISION, INTENT (OUT) :: G(:)
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL BETA, GUSER
    INTERFACE
       SUBROUTINE BETA(T,Y,BVAL)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y
          DOUBLE PRECISION, DIMENSION(:) :: BVAL
          INTENT(IN):: T,Y
          INTENT(OUT) :: BVAL
       END SUBROUTINE BETA
    END INTERFACE
    INTERFACE
       SUBROUTINE GUSER(T,Y,DYDT,Z,G)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
          DOUBLE PRECISION, DIMENSION(:,:) :: Z
          DOUBLE PRECISION, DIMENSION(:) :: G
          INTENT(IN):: T,Y,DYDT,Z
          INTENT(OUT) :: G
       END SUBROUTINE GUSER
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    INTEGER :: I, NG, NGUSP1
    ! ..
    NG = MYIPOINT(6)

    ! Evaluate the residuals for the user root functions.
    IF (MYNGUSER>0) THEN
       IF (MYNEUTRAL) THEN
          IF (NLAGS>0) THEN
             ZANDD(1:MYN,1:NLAGS) = ZARRAY(1:MYN,1:NLAGS)
             ZANDD(1:MYN,NLAGS+1:2*NLAGS) = DARRAY(1:MYN,1:NLAGS)
          END IF
          CALL GUSER(T,Y,DY,ZANDD,G)
       ELSE
          CALL GUSER(T,Y,DY,ZARRAY,G)
       END IF
    END IF

    ! Evaluate the residuals for the discontinuity root functions.
    IF (NG>MYNGUSER) THEN
       IF (JP>0) THEN
          G(MYNGUSER+1:MYNGUSER+JP) = T - TG(MYNGUSER+1:MYNGUSER+JP)
       END IF
       CALL DDE_BETA(T,Y,BETA)
       NGUSP1 = MYNGUSER + JP + 1
       DO I = NGUSP1, NG
          G(I) = MYBVAL(INDEXG(I)) - TG(I)
       END DO
    END IF

    RETURN
  END SUBROUTINE DDE_GRT3
  !____________________________________________________________________________

  SUBROUTINE DDE_GRT4(K,TROOT,ANCESL)

    ! The purpose of DDE_GRT4 is to do bookkeeping for DDE_SOLVER when
    ! new residual functions are added.

    ! Note:
    ! DDE_GRT4 assumes that before it is called (like in DDE_DRV2), K is
    ! set to a positive number and that storage was checked to make sure
    ! there is enough space to make the necessary changes in the
    ! MYJROOT and INDEXG vectors.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: TROOT
    INTEGER, INTENT (IN) :: ANCESL, K
    ! ..
    ! .. Local Scalars ..
    INTEGER :: I, IG, J, NG, NGNEW
    ! ..
    NG = MYIPOINT(6)

    ! The new number of root functions will be NGNEW.
    NGNEW = NG + K

    ! Load TROOT into the new TG slots.
    TG(NG+1:NGNEW) = TROOT

    ! Load exact zeroes into the new residual slots.
    GTOLD(NG+1:NGNEW) = 0.0D0
    GTNEW(NG+1:NGNEW) = 0.0D0
    GTROOT(NG+1:NGNEW) = 0.0D0
    GA(NG+1:NGNEW) = 0.0D0
    GB(NG+1:NGNEW) = 0.0D0
    GTROOT2(NG+1:NGNEW) = 0.0D0

    ! Flag the ode component numbers corresponding to the new root
    ! functions, i.e., those with a root at TROOT.
    ! Note: K = NLAGS*(1+JP+JN).
    DO I = 1, K
       MYJROOT(NG+I) = 1
       LEVEL(NG+I) = ANCESL + 1
    END DO

    ! Set INDEXG(I) for event functions G = MYBVAL - MYJUMPS(*).
    IG = NG
    DO J = 1, 1 + JP + JN
       DO I = 1, NLAGS
          IG = IG + 1
          INDEXG(IG) = I
       END DO
    END DO

    RETURN
  END SUBROUTINE DDE_GRT4
  !____________________________________________________________________________

  SUBROUTINE DDE_GRT5(TVALS,YOLD,GUSER,HINTP,INDEX1,K0,K2,K3,K4,K5,K6,K7, &
       K8,K9,BETA,TINTP,YINIT)

    ! The purpose of DDE_GRT5 is to direct the root finding for
    ! DDE_SOLVER. It is called by subroutine DDE_DRV3 to perform
    ! interpolatory root finding.

    ! In this version two root finding passes are made. The first
    ! searches for roots in (TOLD,TNEW). If a root TROOT1 is
    ! found, a second pass searches for roots in (TROOT1,TNEW).
    ! If a second root TROOT2 is found, the step size is limited
    ! by TROOT2 - TROOT1 when the integration is restarted to
    ! avoid skipping TROOT2 when the integration is continued.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: HINTP, TINTP
    INTEGER, INTENT (INOUT) :: INDEX1
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: K0(:), K2(:), K3(:), K4(:), K5(:), &
         K6(:), K7(:), K8(:), K9(:), YOLD(:)
    DOUBLE PRECISION, INTENT (INOUT) :: TVALS(4)
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL BETA, GUSER, YINIT
    INTERFACE
       SUBROUTINE BETA(T,Y,BVAL)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y
          DOUBLE PRECISION, DIMENSION(:) :: BVAL
          INTENT(IN):: T,Y
          INTENT(OUT) :: BVAL
       END SUBROUTINE BETA
    END INTERFACE
    INTERFACE
       SUBROUTINE GUSER(T,Y,DYDT,Z,G)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
          DOUBLE PRECISION, DIMENSION(:,:) :: Z
          DOUBLE PRECISION, DIMENSION(:) :: G
          INTENT(IN):: T,Y,DYDT,Z
          INTENT(OUT) :: G
       END SUBROUTINE GUSER
    END INTERFACE
    INTERFACE
       SUBROUTINE YINIT(T,Y)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y
         INTENT(IN):: T
         INTENT(OUT) :: Y
      END SUBROUTINE YINIT
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: BNEW, BOLD, GTOL, HDUM1, HDUM2, TNEW, TOLD, TROOT1, &
         TROOT2
    INTEGER :: I, IDUM1, IDUM2, IMETH, IRCALL, JFLAG, NG
    LOGICAL :: ITERATE
    ! ..
    ! .. Local Arrays ..
    DOUBLE PRECISION :: TROOT(2)
    INTEGER :: INDEX(2)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, SIGN
    ! ..
    NG = MYIPOINT(6)

    TOLD = TVALS(1)
    TNEW = TVALS(2)

    INDEX(1) = 2
    INDEX(2) = 2
    TROOT(1) = TOLD
    TROOT(2) = TOLD

    ! IRTMAX = 1 means root finding on (TOLD,TNEW). If a root is located
    ! and IRTMAX = 2, root finding will then be done on the interval
    ! (TROOT,TNEW). If a second root is found, it will be used to limit
    ! the step size when the integration is restarted at TROOT.
    IF (IRTMAX/=2) IRTMAX = 1

    DO 30 IRCALL = 1, IRTMAX

       IF (IRCALL==2) THEN
          ! The second pass is made only if the first pass found
          ! a root with INDEX(1) = 3.
          IF (INDEX(1)/=3) GOTO 40
       END IF

       ! Initialize the root finder.
       IF (IRCALL==1) THEN
          BOLD = TOLD
          BNEW = TNEW
          MYJROOT(1:NG) = 0
          IF (MYNGUSER>0) THEN
             MYREPORT(1:MYNGUSER) = 1
          END IF
       ELSE
          BOLD = TROOT(1)
          BNEW = TNEW
          ! Note: Do not reset MYJROOT!
       END IF
       JFLAG = 0
       GTOL = GRTOL*MAX(ABS(TOLD),ABS(TNEW))
       GTOL = MAX(GTOL,U10)

       ! Load the residuals at the endpoints.
       IF (IRCALL==1) THEN
          GA(1:NG) = GTOLD(1:NG)
          GB(1:NG) = GTNEW(1:NG)
       ELSE
          GA(1:NG) = GTROOT2(1:NG)
          GB(1:NG) = GTNEW(1:NG)
       END IF

       ! Ignore roots not consistent with the direction vector.
       IF (MYNGUSER>0 .AND. IRCALL==1) THEN
          DO I = 1, MYNGUSER
             IF (SIGN(1.0D0,GA(I))>0 .AND. MYDIRECTION(I)==1) THEN
                MYREPORT(I) = 0
             END IF
             IF (SIGN(1.0D0,GA(I))<0 .AND. MYDIRECTION(I)==-1) THEN
                MYREPORT(I) = 0
             END IF
          END DO
       END IF

       ! Call the root finder.
10     CONTINUE
       CALL DDE_GRT2(GTOL,JFLAG,BOLD,BNEW,GTROOT2,TROOT(IRCALL),HDUM1, &
            HDUM2,IDUM1,IDUM2,IRCALL)

       ! Follow the instructions from the root finder.
       IF (JFLAG==1) THEN

          ! Evaluate the residuals for the root finder.

          ! First approximate the solution and the derivative.
          CALL DDE_POLY(TINTP,YOLD,HINTP,K0,K2,K3,K4,K5,K6,K7,K8,K9,MYYROOT, &
               MYDYROOT,TROOT(IRCALL),3)

          ! If user root finding is being done, set the delayed values.
          IF (MYNGUSER>0) THEN
             CALL DDE_BETA(TROOT(IRCALL),MYYROOT,BETA)
             ITERATE = .FALSE.
             IMETH = 1
             CALL DDE_ZSET(TROOT(IRCALL),MYYROOT,HINTP,YINIT,INDEX(IRCALL), &
                  IMETH,K0,K2,K3,K4,K5,K6,K7,K8,K9,TINTP,YOLD,ITERATE)
             ITERATE = .FALSE.
             IF (INDEX(IRCALL)==-12 .OR. INDEX(IRCALL)==-4) THEN
                INDEX(1) = -12
                INDEX(2) = -12
                PRINT *, ' An error occurred when DDE_GRT5 called DDE_ZSET.'
                GOTO 20
             END IF
          END IF

          ! Now calculate the residuals.
          CALL DDE_GRT3(TROOT(IRCALL),MYYROOT,MYDYROOT,GUSER,GTROOT2,BETA)
          NGEVAL = NGEVAL + 1

          ! Check for too many residual evaluations for this call
          ! to DDE_GRT5.
          IF (NGEVAL>NGEMAX) THEN
             INDEX(1) = -10
             INDEX(2) = -10
             GOTO 20
          END IF

          ! Call the root finder again.
          GOTO 10

       END IF

       IF (JFLAG==2) THEN

          ! A root was found.
          ! Approximate the solution and the derivative at the root.
          IF (IRCALL==1) THEN
             CALL DDE_POLY(TINTP,YOLD,HINTP,K0,K2,K3,K4,K5,K6,K7,K8,K9, &
                  MYYROOT,MYDYROOT,TROOT(1),3)
             ! Flag the components which have a root at TROOT(1).
             IF (IRTMAX==2) MYJROOT(1:NG) = -MYJROOT(1:NG)
          ELSE
             ! Restore MYJROOT to its state after the IRCALL = 1 call.
             DO I = 1, NG
                IF (MYJROOT(I)==-1) THEN
                   MYJROOT(I) = 1
                ELSE
                   MYJROOT(I) = 0
                END IF
             END DO
          END IF

          ! Set the flag to indicate a root was found.
          INDEX(IRCALL) = 3

          ! Load an exact zero for the residual.
          IF (IRCALL==1) THEN
             DO I = 1, NG
                IF (MYJROOT(I)==1) GTROOT2(I) = 0.0D0
             END DO
          END IF

          ! Restore the solution for the first root.
          IF (IRCALL==2) THEN
             CALL DDE_POLY(TINTP,YOLD,HINTP,K0,K2,K3,K4,K5,K6,K7,K8,K9, &
                  MYYROOT,MYDYROOT,TROOT(1),3)
          END IF

          ! Return to DDE_DRV3.
          GOTO 20

       END IF

       IF (JFLAG==3) THEN

          ! The right endpoint is a root.
          ! Load the solution and derivative into YROOT and DYROOT.
          MYYROOT(1:MYN) = MYYNEW(1:MYN)
          MYDYROOT(1:MYN) = MYDYNEW(1:MYN)

          ! Set the flag to indicate a root was found.
          INDEX(IRCALL) = 4

          ! Load an exact zero for the residual.
          DO I = 1, NG
             IF (MYJROOT(I)==1) GTROOT2(I) = 0.0D0
             IF (MYJROOT(I)==1) THEN
                MYJROOT(I) = 1
             ELSE
                MYJROOT(I) = 0
             END IF
          END DO

          ! Return to DDE_DRV3.
          GOTO 20

       END IF

       IF (JFLAG==4) THEN

          ! No sign changes were detected.

          ! Set the flag to indicate no sign changes were detected.
          INDEX(IRCALL) = 2

          ! Restore the solution for the first root.
          IF (IRCALL==2) THEN
             CALL DDE_POLY(TINTP,YOLD,HINTP,K0,K2,K3,K4,K5,K6,K7,K8,K9, &
                  MYYROOT,MYDYROOT,TROOT(1),3)
             ! Restore the solution for the first root.
             DO I = 1, NG
                IF (MYJROOT(I)==-1) THEN
                   MYJROOT(I) = 1
                ELSE
                   MYJROOT(I) = 0
                END IF
             END DO
          END IF

          ! Return to DDE_DRV3.

       END IF

20     CONTINUE

30  END DO

40  CONTINUE

    TROOT1 = TROOT(1)
    TROOT2 = TROOT(2)
    IF (INDEX(2)/=3) TROOT2 = TROOT1
    INDEX1 = INDEX(1)
    TVALS(3) = TROOT1
    TVALS(4) = TROOT2

    RETURN
  END SUBROUTINE DDE_GRT5
  !____________________________________________________________________________

  SUBROUTINE DDE_GRT6(AUTORF,NG)

    ! DDE_GRT6 sets up the initial event root finding system for DDE_DRV2.

    ! .. Scalar Arguments ..
    INTEGER, INTENT (OUT) :: NG
    LOGICAL, INTENT (IN) :: AUTORF
    ! ..
    ! .. Local Scalars ..
    INTEGER :: I, IDIR, IG, J
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS
    ! ..
    ! User event functions GUSER(1),...,GUSER(MYNGUSER):
    NG = 0
    MYIPOINT(6) = NG

    IF (MYNGUSER>0) THEN
       NG = MYNGUSER
       MYIPOINT(6) = NG
       TG(1:MYNGUSER) = MYTINIT
       INDEXG(1:MYNGUSER) = 0
       MYJROOT(1:MYNGUSER) = 0
       LEVEL(1:MYNGUSER) = -1
    END IF

    IF (.NOT. (AUTORF)) RETURN

    IG = MYNGUSER
    IDIR = MYIPOINT(5)

    ! JN is the number of jumps before MYTINIT.
    ! JP is the number of jumps after MYTINIT.
    JN = 0
    JP = 0
    IF (MYNJUMPS>0) THEN
       DO I = 1, MYNJUMPS
          IF (IDIR==1) THEN
             IF (MYJUMPS(I)<MYTINIT) JN = JN + 1
             IF (MYJUMPS(I)>MYTINIT) JP = JP + 1
          ELSE
             IF (MYJUMPS(I)>MYTINIT) JN = JN + 1
             IF (MYJUMPS(I)<MYTINIT) JP = JP + 1
          END IF
       END DO
    END IF

    ! Event function G = T - MYJUMPS(*) for MJJUMPS(*)
    ! after MYTINIT.
    IF (JP>0) THEN
       DO I = 1, MYNJUMPS
          IF (IDIR==1 .AND. MYJUMPS(I)>MYTINIT .OR. &
               IDIR==-1 .AND. MYJUMPS(I)<MYTINIT) THEN
             IG = IG + 1
             TG(IG) = MYJUMPS(I)
             INDEXG(IG) = 0
             MYJROOT(IG) = 0
             LEVEL(IG) = -1
          END IF
       END DO
    END IF

    IF (IG/=MYNGUSER+JP) THEN
       PRINT *, ' An error was encountered in DDE_GRT6.'
       STOP
    END IF

    ! MYTINIT is always treated as a jump point.

    ! Event functions G = MYBVAL - MYTINIT.
    IF (NLAGS>0) THEN
       DO I = 1, NLAGS
          IG = IG + 1
          TG(IG) = MYTINIT
          INDEXG(IG) = I
          MYJROOT(IG) = 0
          LEVEL(IG) = 0
       END DO
    END IF

    ! Event functions G = MYBVAL - MYJUMPS(*).
    IF (JP+JN>0) THEN
       DO J = 1, MYNJUMPS
          ! IF (MYJUMPS(J) /= MYTINIT) THEN
          IF (ABS(MYJUMPS(J)-MYTINIT)>0.0D0) THEN
             IF (NLAGS>0) THEN
                DO I = 1, NLAGS
                   IG = IG + 1
                   TG(IG) = MYJUMPS(J)
                   INDEXG(IG) = I
                   MYJROOT(IG) = 0
                   LEVEL(IG) = 0
                END DO
             END IF
          END IF
       END DO
    END IF

    IF (IG/=MYNGUSER+JP+NLAGS+(JP+JN)*NLAGS) THEN
       PRINT *, ' An error was encountered in DDE_GRT6.'
       STOP
    END IF

    NG = MYNGUSER + JP + (JP+JN+1)*NLAGS
    MYIPOINT(6) = NG

    RETURN
  END SUBROUTINE DDE_GRT6
  !____________________________________________________________________________

  SUBROUTINE DDE_SORT(DX,DY,N,KFLAG,IFLAG)

    ! DDE_SORT Sorts array DX and optionally makes the same interchanges
    ! in array DY. The array DX may be sorted in increasing order or
    ! decreasing order. A slightly modified quicksort algorithm is used.

    !     Description of parameters:
    !     DX - array of values to be sorted.
    !     DY - Array to be (optionally) carried along.
    !     N  - Number of values in array DX to be sorted.
    !     KFLAG - Control parameter.
    !           =  2  means sort DX in increasing order and carry DY along.
    !           =  1  means sort DX in increasing order (ignoring DY).
    !           = -1  means sort DX in decreasing order (ignoring DY).
    !           = -2  means sort DX in decreasing order and carry DY along.
    !     Note: This routine is a modified SLATEC routine.

    ! .. Scalar Arguments ..
    INTEGER, INTENT (OUT) :: IFLAG
    INTEGER, INTENT (IN) :: KFLAG, N
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (INOUT) :: DX(:), DY(:)
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: R, T, TT, TTY, TY
    INTEGER :: I, IJ, J, K, KK, L, M, NN
    ! ..
    ! .. Local Arrays ..
    INTEGER :: IL(21), IU(21)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, INT
    ! ..
    IFLAG = 0
    NN = N
    IF (NN<1) THEN
       IFLAG = 6
       PRINT *, ' An array to be sorted in DDE_SORT contains'
       PRINT *, ' only one element.'
       RETURN
    END IF

    KK = ABS(KFLAG)
    IF (KK/=1 .AND. KK/=2) THEN
       IFLAG = 7
       PRINT *, ' An illegal sorting option was encountered in DDE_SORT.'
       RETURN
    END IF

    ! Alter array DX to get decreasing order if needed.
    IF (KFLAG<=-1) THEN
       DX(1:NN) = -DX(1:NN)
    END IF

    IF (KK==2) GOTO 90

    ! Sort DX only.
    M = 1
    I = 1
    J = NN
    R = 0.375D0

10  IF (I==J) GOTO 50
    IF (R<=0.5898437D0) THEN
       R = R + 3.90625D-2
    ELSE
       R = R - 0.21875D0
    END IF

20  K = I
    ! Select a central element of the array and save it in
    ! location T.
    IJ = I + INT((J-I)*R)
    T = DX(IJ)

    ! If first element of array is greater than T, interchange
    ! with T.
    IF (DX(I)>T) THEN
       DX(IJ) = DX(I)
       DX(I) = T
       T = DX(IJ)
    END IF
    L = J

    ! If last element of array is less than than T, interchange
    ! with T.
    IF (DX(J)<T) THEN
       DX(IJ) = DX(J)
       DX(J) = T
       T = DX(IJ)
       ! If first element of array is greater than T, interchange
       ! with T.
       IF (DX(I)>T) THEN
          DX(IJ) = DX(I)
          DX(I) = T
          T = DX(IJ)
       END IF
    END IF

    ! Find an element in the second half of the array which is
    ! smaller than T.
30  L = L - 1
    IF (DX(L)>T) GOTO 30

    ! Find an element in the first half of the array which is
    ! greater than T.
40  K = K + 1
    IF (DX(K)<T) GOTO 40

    ! Interchange these elements.
    IF (K<=L) THEN
       TT = DX(L)
       DX(L) = DX(K)
       DX(K) = TT
       GOTO 30
    END IF

    ! Save upper and lower subscripts of the array yet to
    ! be sorted.
    IF (L-I>J-K) THEN
       IL(M) = I
       IU(M) = L
       I = K
       M = M + 1
    ELSE
       IL(M) = K
       IU(M) = J
       J = L
       M = M + 1
    END IF
    GOTO 60

    ! Begin again on another portion of the unsorted array.
50  M = M - 1
    IF (M==0) GOTO 180
    I = IL(M)
    J = IU(M)

60  IF (J-I>=1) GOTO 20
    IF (I==1) GOTO 10
    I = I - 1

70  I = I + 1
    IF (I==J) GOTO 50
    T = DX(I+1)
    IF (DX(I)<=T) GOTO 70
    K = I

80  DX(K+1) = DX(K)
    K = K - 1
    IF (T<DX(K)) GOTO 80
    DX(K+1) = T
    GOTO 70

    ! Sort DX and carry DY along.
90  M = 1
    I = 1
    J = NN
    R = 0.375D0

100 IF (I==J) GOTO 140
    IF (R<=0.5898437D0) THEN
       R = R + 3.90625D-2
    ELSE
       R = R - 0.21875D0
    END IF

110 K = I
    ! Select a central element of the array and save it in
    ! location T.
    IJ = I + INT((J-I)*R)
    T = DX(IJ)
    TY = DY(IJ)

    ! If first element of array is greater than T, interchange
    ! with T.
    IF (DX(I)>T) THEN
       DX(IJ) = DX(I)
       DX(I) = T
       T = DX(IJ)
       DY(IJ) = DY(I)
       DY(I) = TY
       TY = DY(IJ)
    END IF
    L = J

    ! If last element of array is less than T, interchange
    ! with T.
    IF (DX(J)<T) THEN
       DX(IJ) = DX(J)
       DX(J) = T
       T = DX(IJ)
       DY(IJ) = DY(J)
       DY(J) = TY
       TY = DY(IJ)

       ! If first element of array is greater than T, interchange
       ! with T.
       IF (DX(I)>T) THEN
          DX(IJ) = DX(I)
          DX(I) = T
          T = DX(IJ)
          DY(IJ) = DY(I)
          DY(I) = TY
          TY = DY(IJ)
       END IF
    END IF

    ! Find an element in the second half of the array which is
    ! smaller than T.
120 L = L - 1
    IF (DX(L)>T) GOTO 120

    ! Find an element in the first half of the array which is
    ! greater than T.
130 K = K + 1
    IF (DX(K)<T) GOTO 130
    ! Interchange these elements.
    IF (K<=L) THEN
       TT = DX(L)
       DX(L) = DX(K)
       DX(K) = TT
       TTY = DY(L)
       DY(L) = DY(K)
       DY(K) = TTY
       GOTO 120
    END IF

    ! Save upper and lower subscripts of the array yet to be sorted.
    IF (L-I>J-K) THEN
       IL(M) = I
       IU(M) = L
       I = K
       M = M + 1
    ELSE
       IL(M) = K
       IU(M) = J
       J = L
       M = M + 1
    END IF
    GOTO 150

    ! Begin again on another portion of the unsorted array.
140 M = M - 1
    IF (M==0) GOTO 180
    I = IL(M)
    J = IU(M)

150 IF (J-I>=1) GOTO 110
    IF (I==1) GOTO 100
    I = I - 1

160 I = I + 1
    IF (I==J) GOTO 140
    T = DX(I+1)
    TY = DY(I+1)
    IF (DX(I)<=T) GOTO 160
    K = I

170 DX(K+1) = DX(K)
    DY(K+1) = DY(K)
    K = K - 1
    IF (T<DX(K)) GOTO 170
    DX(K+1) = T
    DY(K+1) = TY
    GOTO 160

    ! Clean up.
180 IF (KFLAG<=-1) THEN
       DX(1:NN) = -DX(1:NN)
    END IF

    RETURN
  END SUBROUTINE DDE_SORT
  !____________________________________________________________________________

  SUBROUTINE DDE_TREE(TINIT,TFINAL,LEVMAX,BETA,Y,IFLAG)

    ! If all delays are constant, DDE_TREE is used to compute the
    ! array of times at which discontinuities may occur. The
    ! automatic root finding normally used in DDE_SOLVER is then
    ! bypassed by having DDE_SOLVER step exactly to each of the
    ! potential discontinuity times TDISC(2),..., TDISC(LENTREE).

    !     Input:
    !     TFINAL  - final integration time (> TINIT).
    !     LEVMAX  - Maximum derivative level to flag for possible
    !               derivative discontinuities
    !               (0<=LEVMAX<=order of integration method).
    !     BETA    - External user subroutine to evaluate delay times.
    !     Y       - Work array of length NEQN (passed as Y to DDE_BETA).

    !     Output:
    !     Note: The contents of MYBVAL(*) will be altered.
    !     Note: DTREE and LENTREE are private and declared elsewhere.
    !     DTREE   - The discontinuity tree.
    !     LENTREE - Length of the discontinuity tree.
    !               If TINIT<TFINAL, the discontinuity times are
    !               DTREE(1) = TINIT<DTREE(2)<...<DTREE(LENTREE) = TFINAL.
    !               If TINIT>TFINAL, the discontinuity times are
    !               DTREE(1) = TINIT>DTREE(2)>...>DTREE(LENTREE) = TFINAL.
    !               Duplicates in the tree are discarded. Values are
    !               considered duplicates if they differ by less than
    !                     TOL = U13*MAX(ABS(TINIT),ABS(TFINAL))
    !               where U13 IS 13 units of roundoff.
    !     IFLAG   - Error flag.
    !               An error occurred if IFLAG is not 0; refer to the
    !               printed error message. The contents of LENTREE and
    !               TREE(*) are valid only if IFLAG = 0.

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: TFINAL, TINIT
    INTEGER, INTENT (OUT) :: IFLAG
    INTEGER, INTENT (IN) :: LEVMAX
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: Y(:)
    ! ..
    ! .. Subroutine Arguments ..
    !EXTERNAL BETA
    INTERFACE
       SUBROUTINE BETA(T,Y,BVAL)
          DOUBLE PRECISION :: T
          DOUBLE PRECISION, DIMENSION(:) :: Y
          DOUBLE PRECISION, DIMENSION(:) :: BVAL
          INTENT(IN):: T,Y
          INTENT(OUT) :: BVAL
       END SUBROUTINE BETA
    END INTERFACE
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: TC, TOL, SPACINGM, SPACINGI
    INTEGER :: ADDED, ADD_TO_TREE, C, I, IDIR, IER, ITEND, J, JDIR, K, &
         LAST, MAXLEN, NUML, NUMLAGS, TAKEOUT
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ABS, DBLE, MAX
    ! ..
    ! Note: This subroutine may change the value of NUMLAGS.
    NUMLAGS = NLAGS

    LENTREE = 0

    ! Check for illegal parameters.
    IF (NUMLAGS<1) THEN
       IFLAG = 2
       PRINT *, ' NUMLAGS must be constant to use DDE_TREE.'
       RETURN
    END IF
    IF (.NOT. (CONSTANT_DELAYS)) THEN
       IFLAG = 2
       PRINT *, ' The delays must be constant in DDE_TREE.'
       RETURN
    END IF
    ! ST 06-02-04:
    ! IF (LEVMAX<0 .OR. LEVMAX>8) THEN
    IF (LEVMAX<0) THEN
       IFLAG = 2
       PRINT *, ' The maximum tracking level must be nonnegative'
       PRINT *, ' to use DDE_TREE.'
       RETURN
    END IF
    IF (ABS(TINIT-TFINAL)<=0.0D0) THEN
       IFLAG = 2
       PRINT *, ' TINIT and TFINAL must be different in DDE_TREE.'
    END IF

    ! Determine the delays and check them.
    ! Calculate the delay times.
    CALL DDE_BETA(TINIT,Y,BETA)
    ! Convert the delay times to the actual delays.
    MYBVAL(1:NUMLAGS) = TINIT - MYBVAL(1:NUMLAGS)

    ! Note: JDIR is not the same as IDIR in DDE_DRV2.
    ! Integration left to right.
    JDIR = -1

    ! Integration right to left.
    IF (TFINAL<TINIT) JDIR = 1
    DO I = 1, NUMLAGS
       IF (DBLE(JDIR)*MYBVAL(I)>0.0D0) THEN
          IFLAG = 2
          PRINT *, ' Illegal delays were encountered in DDE_TREE.'
          RETURN
       END IF
    END DO

    ! Initial allocation of DTREE.
    MAXLEN = 200 !
    ADD_TO_TREE = 100
    ALLOCATE (DTREE(MAXLEN),STAT=IER)
    CALL CHECK_STAT(IER,110)
    GOTO 20

    ! Increase the length of DTREE and start over, if necessary.
10  CONTINUE
    MAXLEN = MAXLEN + ADD_TO_TREE
    DEALLOCATE (DTREE,STAT=IER)
    ALLOCATE (DTREE(MAXLEN),STAT=IER)
20  CONTINUE
    IFLAG = 0
    TOL = U13*MAX(ABS(TINIT),ABS(TFINAL))
    IF (MYNTHIT>0) THEN
       ! Sort the exact hit times.
       CALL DDE_SORT(MYTHIT,MYTHIT,MYNTHIT,JDIR,IFLAG)
       IF (IFLAG/=0) RETURN
       LENTREE = MYNTHIT + 2
       DTREE(1) = TINIT
       DTREE(2:LENTREE-1) = MYTHIT(1:MYNTHIT)
       DTREE(LENTREE) = TFINAL
    ELSE
       LENTREE = 2
       DTREE(1) = TINIT
       DTREE(2) = TFINAL
    END IF

    ! Sort the delays (decreasing if JDIR = -1 and increasing if
    ! JDIR = 1) and discard unnecessary values and those too
    ! close together.
    NUML = NUMLAGS
    IF (NUMLAGS>1) THEN
       CALL DDE_SORT(MYBVAL,MYBVAL,NUMLAGS,JDIR,IFLAG)
       IF (IFLAG/=0) RETURN
       DO I = NUMLAGS, 2, -1
          IF ((ABS(MYBVAL(I-1)-MYBVAL(I))<=TOL) .OR. &
             (ABS(MYBVAL(I)-TINIT)<=0.0D0)) THEN
             TC = MYBVAL(I-1) - MYBVAL(I)
             NUML = NUML - 1
          END IF
       END DO
    END IF
    IF (NUML<1) THEN
       IFLAG = 2
       PRINT *, ' An illegal number of delays was encountered in DDE_TREE.'
       RETURN
    END IF
    NUMLAGS = NUML

    ! Generate the discontinuity tree.
    C = 1
    DTREE(C) = TINIT
    LAST = 1
    ! If it is present, add MYJUMPS to the tree.
    IF (MYNJUMPS>0) THEN
       DO I = 1, MYNJUMPS
          IF (ABS(MYJUMPS(I)-TINIT)>0.0D0) THEN
             C = C + 1
             DTREE(C) = MYJUMPS(I)
             LAST = LAST + 1
          END IF
       END DO
    END IF

    DO K = 1, NUMLAGS
       DO I = 1, LAST
          ADDED = 0
          J = 1
30        CONTINUE
          C = C + 1
          ! Is DTREE big enough to continue?
          IF (MAXLEN<C) GOTO 10
          ! Define the next candidate for DTREE.
          TC = DTREE(I) + DBLE(J)*MYBVAL(K)
          ! Add it if we have not reached TFINAL.
          IF (((TC<=TFINAL) .AND. (JDIR==-1)) .OR. &
              ((TC>=TFINAL) .AND. (JDIR==1))) THEN
             DTREE(C) = TC
             ADDED = ADDED + 1
          ELSE
             ! We reached TFINAL.
             C = C - 1
             GOTO 40
          END IF
          J = J + 1
          ! Have we reached the maximum tracking level?
          IF (LEVMAX>0 .AND. J>LEVMAX) GOTO 40
          GOTO 30
40        CONTINUE
       END DO
       LAST = LAST + ADDED
    END DO
    LENTREE = C

    ! Add TFINAL to the tree, if necessary.
    IF (ABS(DTREE(LENTREE)-TFINAL)>0.0D0) THEN
       C = C + 1
       IF (MAXLEN<C) GOTO 10
       LENTREE = C
       DTREE(LENTREE) = TFINAL
    END IF

    ! Sort the tree (increasing if JDIR = -1 and decreasing if JDIR = 1)
    ! and drop values beyond TFINAL or too close together. Delete entries
    ! that are too close together.
    IF (LENTREE>1) THEN
       CALL DDE_SORT(DTREE,DTREE,LENTREE,-JDIR,IFLAG)
       IF (IFLAG/=0) RETURN
       ITEND = C
       LENTREE = C
       DO I = LENTREE, 1, -1
          IF (((DTREE(I)>TFINAL) .AND. (JDIR==-1)) .OR. &
              ((DTREE(I)<TFINAL) .AND. (JDIR==1))) THEN
             ITEND = ITEND - 1
          END IF
       END DO
       LENTREE = ITEND
       DO I = LENTREE, 2, -1
          IF (ABS(DTREE(I)-DTREE(I-1))<=TOL) THEN
             ! Drop DTREE(I) and move everything below it up.
             IF (I/=LENTREE) THEN
                DO J = I, ITEND-1 ! ST 08/02/2011
                   DTREE(J) = DTREE(J+1)
                END DO
             END IF
             ITEND = ITEND - 1
          END IF
       END DO
       LENTREE = ITEND
    END IF

    ! Ensure that TINIT was not deleted from the tree.
    DO I = 1, LENTREE
       IF (ABS(TINIT-DTREE(I))<=TOL) THEN
          DTREE(I) = TINIT
          GOTO 50
       END IF
    END DO
50  CONTINUE

    ! Delete the jump times and multiple of delays that occur
    ! before TINIT from the tree.
    ! JN is the number of jumps before MYTINIT.
    JN = 0
    IDIR = MYIPOINT(5)
    IF (MYNJUMPS>0) THEN
       DO I = 1, MYNJUMPS
          IF (IDIR==1) THEN
             IF (MYJUMPS(I)<TINIT) JN = JN + 1
          ELSE
             IF (MYJUMPS(I)>TINIT) JN = JN + 1
          END IF
       END DO
    END IF
    TAKEOUT = 0
    IF (JN>0) THEN
       DO I = 1, LENTREE
          IF (ABS(DTREE(I)-TINIT)<=0.0D0) THEN
             TAKEOUT = I - 1
             GOTO 60
          END IF
       END DO
60     CONTINUE
       IF (TAKEOUT>0) THEN
          DTREE(1:LENTREE-TAKEOUT) = DTREE(TAKEOUT+1:LENTREE)
          LENTREE = LENTREE - TAKEOUT
       END IF
    END IF

    ! Trim the discontinuity tree.
    IF (IFLAG==0) THEN
       ALLOCATE (TREE_TEMP(LENTREE),STAT=IER)
       CALL CHECK_STAT(IER,111)
       TREE_TEMP(1:LENTREE) = DTREE(1:LENTREE)
       DEALLOCATE (DTREE,STAT=IER)
       CALL CHECK_STAT(IER,112)
       ALLOCATE (DTREE(LENTREE),STAT=IER)
       CALL CHECK_STAT(IER,113)
       DTREE(1:LENTREE) = TREE_TEMP(1:LENTREE)
       DEALLOCATE (TREE_TEMP,STAT=IER)
       CALL CHECK_STAT(IER,114)

       IF (DEBUG) THEN
          WRITE(DBGUNIT,9997)
          WRITE(DBGUNIT,9999) LENTREE
          SPACINGM = 1D20
          SPACINGI = 0D0
          WRITE(DBGUNIT,9996)
          WRITE(DBGUNIT,9998) I, DTREE(1), SPACINGI
          DO I=2,LENTREE
             SPACINGI = DTREE(I) - DTREE(I-1)
             SPACINGM = MIN(SPACINGM, SPACINGI)
             WRITE(DBGUNIT,9998) I, DTREE(I), SPACINGI 
          END DO
          WRITE(DBGUNIT,9995) SPACINGM
          9999 FORMAT(' The discontinuity tree follows. Length = ', I10)
          9998 FORMAT(I10, 2D20.10)
          9997 FORMAT(' ')
          9996 FORMAT(5X, ' Component', 6X, 'Value', 13X, 'Spacing')
          9995 FORMAT(2X, 'Minimum spacing = ', D20.10)
       END IF

    END IF

    RETURN
  END SUBROUTINE DDE_TREE
  !____________________________________________________________________________

  SUBROUTINE DDE_USER(TVAL,YOFT,DYOFT)

    ! User callable subroutine to interpolate the solution and the
    ! derivative. This subroutine is not called by DDE_SOLVER.

    !     Usage:
    !     CALL DDE_USER(TVAL,YOFT,DYOFT)
    !     Input:
    !     TVAL           -  Value of T at which interpolation of the
    !                       solution and the derivative is desired
    !                       is desired. TVAL must be greater then
    !                       the initial time.
    !     YOFT           -  Array of length NEQN (= the number of ddes).
    !     DYOFT          -  Array of length NEQN.
    !     Output:
    !     YOFT           -  Interpolated solution Y(TVAL).
    !     DYOFT          -  Interpolated derivative Y(TVAL).

    ! .. Scalar Arguments ..
    DOUBLE PRECISION, INTENT (IN) :: TVAL
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, INTENT (OUT) :: DYOFT(:), YOFT(:)
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: C, DELINT, TN, TO
    INTEGER :: IBEGIN, IDIR, INDEXO, INEW, IOLD, IQUEUE, ISTART, J, &
         JQUEUE, N
    LOGICAL :: POLYD, POLYS
    ! ..
    ! .. Local Arrays ..
    DOUBLE PRECISION :: W(0:9), WD(0:9)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC ALLOCATED
    ! ..
    ! Check if this is a legal call.
    IF (.NOT. (ALLOCATED(QUEUE))) THEN
       PRINT *, ' DDE_USER may only be called by the user routines'
       PRINT *, ' DDES, BETA, HISTORY, EF, and CHNG. Use DDE_VAL'
       PRINT *, ' to perform the requested interpolation.'
       RETURN
    END IF

    N = MYN

    ! IQUEUE points to the most recent addition to the queue.
    ! JQUEUE points to the oldest addition to the queue.
    IQUEUE = MYIPOINT(2)
    JQUEUE = MYIPOINT(3)

    ! Set some method flags.
    IDIR = MYIPOINT(5)
    IPAIR = 1
    IORDER = 6
    MYMETHOD = 4
    POLYS = .TRUE.
    POLYD = .TRUE.

    ! Check if TVAL falls before the initial time.
    IF ((IDIR==1 .AND. TVAL<TQUEUE(0)) .OR. &
       (IDIR==-1 .AND. TVAL>TQUEUE(0))) THEN
       PRINT *, ' T falls before the initial time.'
       PRINT *, ' There is insufficient information in the solution'
       PRINT *, ' queue to perform the interpolation requested in'
       PRINT *, ' DDE_USER. Use your history function to obtain the'
       PRINT *, ' interpolated solution requested in DDE_USER.'
       GOTO 20
    END IF

    ! Check whether an integration step has been completed. Extrapolate
    ! the initial solution if it hasn't. (The initial solution is
    ! loaded in the queue in DDE_DRV1.)
    IF (IQUEUE==0) THEN
       YOFT(1:N) = QUEUE(1:N,1)
       DYOFT(1:N) = 0
       GOTO 20
    END IF

    ! Check if T falls beyond the last point in the queue. Extrapolate
    ! the most recent solution polynomial if it does.
    IF ((IDIR==1 .AND. TVAL>TQUEUE(IQUEUE)) .OR. &
       (IDIR==-1 .AND. TVAL<TQUEUE(IQUEUE))) THEN
       TO = TQUEUE(IQUEUE-1)
       TN = TQUEUE(IQUEUE)
       INDEXO = IQUEUE
    ELSE
       ! Otherwise, TVAL IS SPANNED by the solution queue;
       ! so interpolate using the solution history queue.
       ! First bracket TVAL in the queue.
       ISTART = 0
       IOLD = JQUEUE
       INEW = IQUEUE
       IF (IOLD/=1) THEN
          PRINT *, ' IOLD is not 1 in DDE_USER.'
          GOTO 20
       END IF
       ! The bracketing search begins here.
       ! Linear search; remember last bracketing interval.
       IF (IDIR==1) THEN
          ! TQUEUE is increasing.
          IF (TVAL<=TQUEUE(ISTART)) THEN
             DO J = ISTART, 0, -1
                IF (TVAL>=TQUEUE(J)) THEN
                   ! TQUEUE(J)<=TVAL<=TQUEUE(J+1)
                   ISTART = J
                   INDEXO = J + 1
                   TO = TQUEUE(J)
                   TN = TQUEUE(INDEXO)
                   GOTO 10
                END IF
             END DO
          ELSE
             ! TVAL>TQUEUE(ISTART)
             DO J = ISTART + 1, INEW
                IF (TVAL<=TQUEUE(J)) THEN
                   ! TQUEUE(J-1)<=TVAL<=TQUEUE(J)
                   ISTART = J - 1
                   INDEXO = J
                   TO = TQUEUE(ISTART)
                   TN = TQUEUE(INDEXO)
                   GOTO 10
                END IF
             END DO
          END IF
       ELSE
          ! TQUEUE is decreasing.
          IF (TVAL>=TQUEUE(ISTART)) THEN
             DO J = ISTART, 0, -1
                IF (TVAL<=TQUEUE(J)) THEN
                   ! TQUEUE(J)>=TVAL>=TQUEUE(J+1)
                   ISTART = J
                   INDEXO = J + 1
                   TO = TQUEUE(J)
                   TN = TQUEUE(INDEXO)
                   GOTO 10
                END IF
             END DO
          ELSE
             ! TVAL<TQUEUE(ISTART)
             DO J = ISTART + 1, INEW
                IF (TVAL>=TQUEUE(J)) THEN
                   ! TQUEUE(J-1)>=TVAL>=TQUEUE(J)
                   ISTART = J - 1
                   INDEXO = J
                   TO = TQUEUE(ISTART)
                   TN = TQUEUE(INDEXO)
                   GOTO 10
                END IF
             END DO
          END IF
       END IF
    END IF
    ! The bracketing search ends here.

10  CONTINUE
    ! The data to be interpolated corresponds to TO = TOLD, TN = TNEW.
    ! INDEXO is the pointer for the corresponding slug of data. The
    ! data for TOLD to be interpolated begins in queue location
    ! (1,IBEGIN+1).
    IBEGIN = (INDEXO-1)*10

    ! Define the step size corresponding to the interpolation data.
    DELINT = TN - TO

    ! Determine the value of C at which the solution or the
    ! derivative is to be interpolated.
    IF (ABS(DELINT)<=0.0D0) THEN
       C = 0.0D0
    ELSE
       C = (TVAL-TO)/DELINT
    END IF

    ! Calculate the interpolation coefficients.
    CALL DDE_WSET(C,W,WD,POLYS,POLYD,MYMETHOD)

    ! Interpolate the solution.
    YOFT(1:N) = QUEUE(1:N,IBEGIN+1) + DELINT * &
        (W(0)*QUEUE(1:N,IBEGIN+2)+W(2)*QUEUE(1:N,IBEGIN+3)+ &
         W(3)*QUEUE(1:N,IBEGIN+4)+W(4)*QUEUE(1:N,IBEGIN+5)+ &
         W(5)*QUEUE(1:N,IBEGIN+6)+W(6)*QUEUE(1:N,IBEGIN+7)+ &
         W(7)*QUEUE(1:N,IBEGIN+8)+ W(8)*QUEUE(1:N,IBEGIN+9)+ &
         W(9)*QUEUE(1:N,IBEGIN+10))

    ! Interpolate the derivative.
    DYOFT(1:N) = &
         WD(0)*QUEUE(1:N,IBEGIN+2) + WD(2)*QUEUE(1:N,IBEGIN+3) + &
         WD(3)*QUEUE(1:N,IBEGIN+4) + WD(4)*QUEUE(1:N,IBEGIN+5) + &
         WD(5)*QUEUE(1:N,IBEGIN+6) + WD(6)*QUEUE(1:N,IBEGIN+7) + &
         WD(7)*QUEUE(1:N,IBEGIN+8) + WD(8)*QUEUE(1:N,IBEGIN+9) + &
         WD(9)*QUEUE(1:N,IBEGIN+10)

20  CONTINUE

    RETURN
  END SUBROUTINE DDE_USER
  !____________________________________________________________________________

  SUBROUTINE DDE_IQUEUE(IQ)

    ! User callable companion subroutine for DDE_USER to retrieve
    ! the most recent history queue index.

    ! Usage:
    ! CALL DDE_IQUEUE(IQ)
    ! Output:
    ! IQ - Index for the most recent queue entry.

    ! .. Scalar Arguments ..
    INTEGER, INTENT (OUT) :: IQ
    ! ..
    ! IQ points to the most recent addition to the queue.
    IQ = MYIPOINT(2)

    RETURN
  END SUBROUTINE DDE_IQUEUE
  !____________________________________________________________________________

  SUBROUTINE DDE_VAL_USER(TVALS,NT,SOL,DERIVATIVES,YT,DT)

    ! Interpolate the solution and derivative using the solution
    ! structure sol produced by a call to DDE_SOLVER with
    ! INTERPOLATION = .TRUE. This subroutine is not called by
    ! DDE_SOLVER nor is it intended for general use. Its intended
    ! use is for testing when the SOL solution structure from a
    ! previous integration must be used a large number of times
    ! in which event the linear bracket searches and the
    ! allocation and deallocation of arrays in DDE_VAL is too
    ! time consuming.

    !     Usage:
    !     CALL DDE_VAL_USER(TVALS,NT,SOL,DERIVATIVES,YT,DT)
    !     On Entry:
    !     TVALS       -  Array of length NT>0 at which interpolation
    !                    of the solution is desired.
    !     NT          -  Length of TVALS array.
    !     SOL         -  Solution structure produced by the a previous
    !                    call to DDE_SOLVER.
    !     DERIVATIVES -  Logical flag which indicates whether
    !                    derivatives are to be calculated.
    !     YT,DT       -  Arrays of size (NT,N) where N=SIZE(SOL%QUEUE,1)
    !                    is the number of ddes in the call to DDE_SOLVER
    !                    that produced SOL.
    !     On Return:
    !     YT - Interpolated solution. Array with dimensions NT by N.
    !                       YT(I,J) is an approximation to solution
    !                       component J at time TVALS(I).
    !     DT - Interpolated derivatives. Array with dimensions
    !                       NT by N. DT(I,J) is an approximation to
    !                       derivative component J at time TVALS(I).

    ! .. Structure Arguments ..
    TYPE (DDE_SOL), TARGET :: SOL
    ! ..
    ! .. Scalar Arguments ..
    INTEGER :: NT
    LOGICAL :: DERIVATIVES
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION :: DT(NT,*), TVALS(NT), YT(NT,*)
    ! ..
    ! .. Local Scalars ..
    DOUBLE PRECISION :: C, DELINT, TN, TO, TVAL
    INTEGER :: I, IBEGIN, IDIR, INDEXO, INEW, IOLD, IQ, IQUEUE, JPQ, JQ, &
             JQUEUE, N
    LOGICAL :: POLYD, POLYS
    ! ..
    ! .. Local Arrays ..
    DOUBLE PRECISION :: W(0:9), WD(0:9)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC SIZE
    ! ..
    ! Check dimensions.
    N = SIZE(SOL%QUEUE,1)
    IF (N<0) THEN
       PRINT *, ' The number of DDEs must be positive in DDE_VAL_USER.'
    END IF

    IF (NT<1) THEN
       PRINT *, ' The length of TVALS must be positive in DDE_VAL_USER.'
       STOP
    END IF

    POLYS = .TRUE.
    ! Set the local calculate derivatives flag.
    IF (DERIVATIVES) THEN
       POLYD = .TRUE.
    ELSE
       POLYD = .FALSE.
    END IF

    ! IQUEUE points to the most recent addition to the queue.
    ! JQUEUE points to the oldest addition to the queue.
    IQUEUE = SOL%IPOINT(2)
    JQUEUE = SOL%IPOINT(3)

    ! Check that an integration step has been completed.
    IF (IQUEUE==0) THEN
       PRINT *, ' The first integration step has not been completed.'
       PRINT *, ' There is insufficient information in the solution'
       PRINT *, ' history queue to perform the interpolation'
       PRINT *, ' requested in DDE_VAL_USER.'
       GOTO 90
    END IF

    ! We land here if T is spanned by the solution queue.
    ! Interpolate  using the solution history queue.

    ! Set some method flags.
    IDIR = SOL%IPOINT(5)
    ! These are the same values used in DDE_SOLVER, so no harm.
    IPAIR = 1
    IORDER = 6
    MYMETHOD = 4

    ! The interpolation loop begins here.
    DO 80 I = 1, NT

       TVAL = TVALS(I)

       ! Check if TVAL falls before the initial time.
       IF ((IDIR==1 .AND. TVAL<SOL%TQUEUE(0)) .OR. &
          (IDIR==-1 .AND. TVAL>SOL%TQUEUE(0))) THEN
          PRINT *, ' T falls before the initial time.'
          PRINT *, ' There is insufficient information in the solution'
          PRINT *, ' queue to perform the interpolation requested in'
          PRINT *, ' DDE_VAL. Use your history function to obtain the'
          PRINT *, ' interpolated solution requested in DDE_VAL_USER.'
          GOTO 90
       END IF

       ! Check if T falls beyond the last point in the queue.
       IF ((IDIR==1 .AND. TVAL>SOL%TQUEUE(IQUEUE)) .OR. &
          (IDIR==-1 .AND. TVAL<SOL%TQUEUE(IQUEUE))) THEN
          PRINT *, ' The requested value of T is not spanned by the'
          PRINT *, ' solution queue. It is not possible to perform'
          PRINT *, ' interpolation requested in DDE_VAL_USER.'
          GOTO 90
       END IF

       ! Bracket TVAL in the queue.
       IOLD = JQUEUE
       INEW = IQUEUE
       IF (IOLD/=1) THEN
          PRINT *, ' IOLD is not 1 in DDE_VAL_USER.'
       END IF

       ! The binary bracketing search begins here.
       IF (IDIR==1) THEN
          ! SOL%TQUEUE is increasing.
          IF (TVAL>=SOL%TQUEUE(INEW-1)) THEN
             INDEXO = INEW
             GOTO 70
          END IF
          IQ = IOLD
          JQ = INEW
10        CONTINUE
          IF (JQ-IQ==1) GOTO 30
          JPQ = (JQ+IQ)/2
          IF (TVAL<SOL%TQUEUE(JPQ-1)) GOTO 20
          IQ = JPQ
          GOTO 10
20        CONTINUE
          JQ = JPQ
          GOTO 10
30        CONTINUE
          INDEXO = IQ
          GOTO 70
       ELSE
          ! SOL%TQUEUE is decreasing.
          IF (TVAL<=SOL%TQUEUE(INEW-1)) THEN
             INDEXO = INEW
             GOTO 70
          END IF
          IQ = IOLD
          JQ = INEW
40        CONTINUE
          IF (JQ-IQ==1) GOTO 60
          JPQ = (JQ+IQ)/2
          IF (TVAL>SOL%TQUEUE(JPQ-1)) GOTO 50
          IQ = JPQ
          GOTO 40
50        CONTINUE
          JQ = JPQ
          GOTO 40
60        CONTINUE
          INDEXO = IQ
          GOTO 70
       END IF

       ! The bracketing search ends here.
70     CONTINUE
       TO = SOL%TQUEUE(INDEXO-1)
       TN = SOL%TQUEUE(INDEXO)
       IF (((IDIR==1) .AND. (TVAL<TO .OR. TVAL>TN)) .OR. &
          ((IDIR==-1) .AND. (TVAL>TO .OR. TVAL<TN))) THEN
          PRINT *, ' A search error occurred in DDE_VAL_USER.'
          STOP
       END IF

       ! THE data to be interpolated corresponds to TO = TOLD, TN = TNEW.
       ! INDEXO is the pointer for the corresponding slug of data. The
       ! data for TOLD to be interpolated begins in queue location
       ! (1,IBEGIN+1).
       IBEGIN = (INDEXO-1)*10

       ! Define the step size corresponding to the interpolation data.
       DELINT = TN - TO

       ! Determine the value of C at which the solution or the
       ! derivative is to be interpolated.
       IF (ABS(DELINT)<=0.0D0) THEN
          C = 0.0D0
       ELSE
          C = (TVAL-TO)/DELINT
       END IF

       ! Calculate the interpolation coefficients.
       CALL DDE_WSET(C,W,WD,POLYS,POLYD,MYMETHOD)

       ! Interpolate the solution.
       YT(I,1:N) = SOL%QUEUE(1:N,IBEGIN+1) + DELINT*(W(0)*SOL%QUEUE(1:N, &
            IBEGIN+2)+W(2)*SOL%QUEUE(1:N,IBEGIN+3)+W(3)*SOL%QUEUE(1:N,IBEGIN+4 &
            )+W(4)*SOL%QUEUE(1:N,IBEGIN+5)+W(5)*SOL%QUEUE(1:N,IBEGIN+6)+ &
            W(6)*SOL%QUEUE(1:N,IBEGIN+7)+W(7)*SOL%QUEUE(1:N,IBEGIN+8)+ &
            W(8)*SOL%QUEUE(1:N,IBEGIN+9)+W(9)*SOL%QUEUE(1:N,IBEGIN+10))

       IF (POLYD) THEN
          ! Interpolate the derivative.
          DT(I,1:N) = WD(0)*SOL%QUEUE(1:N,IBEGIN+2) + &
               WD(2)*SOL%QUEUE(1:N,IBEGIN+3) + WD(3)*SOL%QUEUE(1:N,IBEGIN+4) + &
               WD(4)*SOL%QUEUE(1:N,IBEGIN+5) + WD(5)*SOL%QUEUE(1:N,IBEGIN+6) + &
               WD(6)*SOL%QUEUE(1:N,IBEGIN+7) + WD(7)*SOL%QUEUE(1:N,IBEGIN+8) + &
               WD(8)*SOL%QUEUE(1:N,IBEGIN+9) + WD(9)*SOL%QUEUE(1:N,IBEGIN+10)
       END IF

       ! The interpolation loop ends here.
80  END DO

90  CONTINUE

    RETURN
  END SUBROUTINE DDE_VAL_USER
  !____________________________________________________________________________
        SUBROUTINE DAFTER(N,Y,DY,T,H)

  ! Perform end of successful step bookkeeping (not called, as stands).

  ! .. Scalar Arguments ..
         INTEGER :: N
  ! ..
  ! .. Array Arguments ..
          DOUBLE PRECISION :: Y(N), DY(N)
  ! ..
  ! .. Local Scalars ..
          DOUBLE PRECISION :: T, H, TEMPDUMMY
          INTEGER IDUMMY
  ! ..
        IDUMMY = 0
        IF (IDUMMY /= 0) THEN
           TEMPDUMMY = Y(1)
           TEMPDUMMY = DY(1)
           TEMPDUMMY = T
           TEMPDUMMY = H
           T = TEMPDUMMY
        END IF
        IF (DEBUG) THEN
           WRITE(DBGUNIT,9999) T, H
           9999 FORMAT(2X, ' Successful step taken. T and H = ', 2D20.10)
        END IF
  ! ..
        RETURN
      END SUBROUTINE DAFTER
  !____________________________________________________________________________
  SUBROUTINE RELEASE_ARRAYS(SOL,OPTS)

    ! Deallocate work arrays and solution structure fields.
    ! NB: Do not call this subroutine before a return has been made
    ! by DDE_SOLVER to the calling program.

    ! This subroutine is not called by DDE_SOLVER. To prevent a memory
    ! leak from occurring, it should be called after the solution has
    ! been processed by the calling program in order the deallocate the
    ! fields in the solution structure SOL, and in the options structure
    ! OPTS.

    ! PLEASE NOTE:
    ! If an error occurs when this subroutine is called, it was not
    ! possible to DEALLOCATE an array. Please contact Skip Thompson
    ! in this case and provide the exact error message. To continue
    ! using the solver until the problem is corrected, comment out
    ! the call to CHECK_STAT in the line following the attempted
    ! deallocation that caused the error message.

    ! .. Structure Arguments ..
    TYPE (DDE_OPTS) :: OPTS
    TYPE (DDE_SOL) :: SOL
    ! ..
    ! .. Local Scalars ..
    INTEGER :: IER, L1, L2, TOTAL
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC SIZE
    ! ..

    TOTAL = 0

    IF (SOL%IPOINT(8)==1) THEN
       IER = 0
       L1 = SIZE(SOL%T)
       TOTAL = TOTAL + L1
       IF (L1>0) DEALLOCATE (SOL%T,STAT=IER)
       CALL CHECK_STAT(IER,1150)
    END IF

    IF (SOL%IPOINT(9)==1) THEN
       IER = 0
       L1 = SIZE(SOL%Y,1)
       L2 = SIZE(SOL%Y,2)
       TOTAL = TOTAL + L1*L2
       IF (L1*L2>0) DEALLOCATE (SOL%Y,STAT=IER)
       CALL CHECK_STAT(IER,1151)
    END IF

    IF (SOL%IPOINT(10)==1) THEN
       IER = 0
       L1 = SIZE(SOL%TE)
       TOTAL = TOTAL + L1
       IF (L1>0) DEALLOCATE (SOL%TE,STAT=IER)
       CALL CHECK_STAT(IER,1152)
    END IF

    IF (SOL%IPOINT(11)==1) THEN
       IER = 0
       L1 = SIZE(SOL%YE,1)
       L2 = SIZE(SOL%YE,2)
       TOTAL = TOTAL + L1*L2
       IF (L1*L2>0) DEALLOCATE (SOL%YE,STAT=IER)
       CALL CHECK_STAT(IER,1153)
    END IF

    IF (SOL%IPOINT(12)==1) THEN
       IER = 0
       L1 = SIZE(SOL%QUEUE,1)
       L2 = SIZE(SOL%QUEUE,2)
       TOTAL = TOTAL + L1*L2
       IF (L1*L2>0) DEALLOCATE (SOL%QUEUE,STAT=IER)
       CALL CHECK_STAT(IER,1154)
    END IF

    IF (SOL%IPOINT(13)==1) THEN
       IER = 0
       L1 = SIZE(SOL%YOFT)
       TOTAL = TOTAL + L1
       IF (L1>0) DEALLOCATE (SOL%YOFT,STAT=IER)
       CALL CHECK_STAT(IER,1155)
    END IF

    IF (SOL%IPOINT(14)==1) THEN
       IER = 0
       L1 = SIZE(SOL%TQUEUE)
       TOTAL = TOTAL + L1
       IF (L1>0) DEALLOCATE (SOL%TQUEUE,STAT=IER)
       CALL CHECK_STAT(IER,1156)
    END IF

    IF (SOL%IPOINT(15)==1) THEN
       IER = 0
       L1 = SIZE(SOL%IE)
       TOTAL = TOTAL + L1
       IF (L1>0) DEALLOCATE (SOL%IE,STAT=IER)
       CALL CHECK_STAT(IER,1157)
    END IF

    IF (SOL%IPOINT(17)==1) THEN
       IER = 0
       L1 = SIZE(SOL%STATS)
       TOTAL = TOTAL + L1
       IF (L1>0) DEALLOCATE (SOL%STATS,STAT=IER)
       CALL CHECK_STAT(IER,1158)
    END IF

    IF (SOL%IPOINT(22)==1) THEN
       IF (ASSOCIATED(OPTS%ISTERMINAL)) THEN
          IER = 0
          L1 = SIZE(OPTS%ISTERMINAL)
          TOTAL = TOTAL + L1
          IF (L1>0) DEALLOCATE (OPTS%ISTERMINAL,STAT=IER)
          CALL CHECK_STAT(IER,1159)
       END IF
    END IF

    IF (SOL%IPOINT(23)==1) THEN
       IF (ASSOCIATED(OPTS%DIRECTION)) THEN
          IER = 0
          L1 = SIZE(OPTS%DIRECTION)
          TOTAL = TOTAL + L1
          IF (L1>0) DEALLOCATE (OPTS%DIRECTION,STAT=IER)
          CALL CHECK_STAT(IER,1160)
       END IF
    END IF

    !IF (SOL%IPOINT(24)==1) THEN
    !   IER = 0
    !   L1 = SIZE(OPTS%ABSERR)
    !   TOTAL = TOTAL + L1
    !   DEALLOCATE (OPTS%ABSERR,STAT=IER)
    !   CALL CHECK_STAT(IER,1161)
    !END IF

    !IF (SOL%IPOINT(25)==1) THEN
    !   IF (ASSOCIATED(OPTS%RELERR)) THEN
    !      IER = 0
    !      L1 = SIZE(OPTS%RELERR)
    !      TOTAL = TOTAL + L1
    !      DEALLOCATE (OPTS%RELERR,STAT=IER)
    !      CALL CHECK_STAT(IER,1162)
    !   END IF
    !END IF

    IF (SOL%IPOINT(26)==1) THEN
       IF (ASSOCIATED(OPTS%JUMPS)) THEN
          IER = 0
          L1 = SIZE(OPTS%JUMPS)
          TOTAL = TOTAL + L1
          DEALLOCATE (OPTS%JUMPS,STAT=IER)
          CALL CHECK_STAT(IER,1163)
       END IF
    END IF

    IF (SOL%IPOINT(27)==1) THEN
       IF (ASSOCIATED(OPTS%THIT_EXACTLY)) THEN
          IER = 0
          L1 = SIZE(OPTS%THIT_EXACTLY)
          TOTAL = TOTAL + L1
          DEALLOCATE (OPTS%THIT_EXACTLY,STAT=IER)
          CALL CHECK_STAT(IER,1164)
       END IF
    END IF

    IER = 0
    L1 = SIZE(SOL%IPOINT)
    TOTAL = TOTAL + L1
    DEALLOCATE (SOL%IPOINT,STAT=IER)
    CALL CHECK_STAT(IER,1165)

    PRINT *, ' In subroutine RELEASE_ARRAYS, total number of words'
    PRINT *, ' released = ', TOTAL

    RETURN
  END SUBROUTINE RELEASE_ARRAYS
  !____________________________________________________________________________

  SUBROUTINE RELEASE_INT(YINT)

    ! Deallocate interpolation structure fields.

    ! This subroutine is not called by DDE_SOLVER. To prevent a memory
    ! leak from occurring, it should be called after the solution has
    ! been processed by the calling program in order the deallocate
    ! the fields in any interpolation structure YINT produced by the
    ! interpolation module DDE_VAL. A deallocation error will occur
    ! if YINT is not a valid solution structure containing fields
    ! YINT%YINT, YINT%YT, YINT%DT and YINT%COMPONENTS.

    ! .. Structure Arguments ..
    TYPE (DDE_INT) :: YINT
    ! ..
    ! .. Local Scalars ..
    INTEGER IER, L1, L2, TOTAL
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC SIZE
    ! ..

    TOTAL = 0

    L1 = SIZE(YINT%TVALS)
    TOTAL = TOTAL + L1
    DEALLOCATE (YINT%TVALS,STAT=IER)
    CALL CHECK_STAT(IER,1200)

    L1 = SIZE(YINT%YT,1)
    L2 = SIZE(YINT%YT,2)
    TOTAL = TOTAL + L1*L2
    DEALLOCATE (YINT%YT,STAT=IER)
    CALL CHECK_STAT(IER,1201)

    L1 = SIZE(YINT%DT,1)
    L2 = SIZE(YINT%DT,2)
    TOTAL = TOTAL + L1*L2
    DEALLOCATE (YINT%DT,STAT=IER)
    CALL CHECK_STAT(IER,1202)

    L1 = SIZE(YINT%COMPONENTS)
    TOTAL = TOTAL + L1
    DEALLOCATE (YINT%COMPONENTS,STAT=IER)
    CALL CHECK_STAT(IER,1203)
    PRINT *, ' In subroutine RELEASE_INT, the total number'
    PRINT *, ' of words released was ', TOTAL, '.'
    RETURN
  END SUBROUTINE RELEASE_INT

  FUNCTION DDE_SOLVER_VERSION() RESULT(VERSION)
    TYPE(DDE_SOLVER_VERSION_TYPE) :: VERSION
    VERSION%MAJOR = VERSION_MAJOR
    VERSION%MINOR = VERSION_MINOR
    VERSION%POINT = VERSION_POINT
    VERSION%RELEASED = VERSION_RELEASED
    RETURN
  END FUNCTION DDE_SOLVER_VERSION

!************END PUBLIC AUXILIARY SUBROUTINES*********************

END MODULE DDE_SOLVER_M
!____________________________________________________________________________
