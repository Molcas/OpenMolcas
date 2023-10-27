!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine CHO_SCF_RDINP(DFonly,LuSpool)
!***********************************************************************
!
!  Purpose:   If DFonly, use defaults only.
!             Else, read and process input for Cholesky section in SCF
!
!***********************************************************************
      Use Fock_util_global, only: Deco, DensityCheck, Estimate, Update
      Use Cholesky, only: ChFracMem, timings
      Use ChoSCF, only: ALGO, dmpk, nScreen, ReOrd
      use Constants, only: Zero, Half
      Implicit None
      Logical  DFonly
      Integer LuSpool

#include "print.fh"
      Character(LEN=180) KWord, Key
      Character(LEN=180), External :: Get_Ln
      character(LEN=13), parameter:: SECNAM = 'CHO_SCF_RDINP'
      Real*8 dmpk_dfl
      Integer i, iChrct, iPrint, iRout, jRout, Last, n
      Integer, External::iCLast
!
      iRout=1
      iPrint=nPrint(iRout)
!                                                                      *
!**** Algorithms for using Cholesky vectors in SCF *********************
!
!   ALGO:
!          0  --->  Integrals are regenerated on the fly
!                   from a set of Cholesky vectors resorted on disk
!                   (Used only for debugging! Amazingly slow!)
!
!          1  --->  The resorted Cholesky vectors are used directly
!                   by the Fock matrix builder routines and contracted
!                   with the proper density matrices. Uses
!                   vectors resorted either on disk or on the fly
!                   (Used only for debugging! Very slow!)
!
!          2  --->  As in option 1 but using the MO-basis transformed
!                   vectors for computing the exchange term
!
!          3  --->  As in option 2 but using the MO-basis vectors
!                   transformed directly in reduced sets
!
!          4  --->  Local-exchange (LK) algorithm for the exchange term
!
!***********************************************************************
!                                                                      *
!     Default  parameters

#ifdef _MOLCAS_MPP_
      ChFracMem=0.3d0
#else
      ChFracMem=half
#endif

      IF (DFonly) THEN
         ALGO  = 4
         REORD =.false.
         DECO  =.true.
         DensityCheck=.false.
         timings=.false.
         NSCREEN = 10    ! default screening interval (# of red sets)
         dmpk = 0.045d0   ! default damping of the screening threshold
         Estimate = .false.
         Update = .true.
         goto 999  !return flag
      ENDIF
!
!    set some parameters if not specified in ChoInput section
         ALGO  = 4
         REORD =.false.
         DECO  =.true.
         DensityCheck=.false.
         timings=.false.
         NSCREEN = 10
         dmpk = 0.045d0
         Estimate = .false.
         Update = .true.

         dmpk_dfl = dmpk
!***********************************************************************
!                                                                      *
      iPrint=5
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Process the input
!
!-------------------------------------------------------------------*
! The big turning point.                                            *
!-------------------------------------------------------------------*
1000  Continue
!-------------------------------------------------------------------*
! Use Get_Ln to read the lines.                                     *
!-------------------------------------------------------------------*
      Key=Get_Ln(LuSpool)
      Kword=Key
      Call UpCase(Kword)
!-------------------------------------------------------------------*
! The keywords and their labels.                                    *
!-------------------------------------------------------------------*

      If (KWord(1:1).eq.'*')    Go To 1000
      If (KWord.eq.'')       Go To 1000
      If (KWord(1:4).eq.'ALGO') Go To 900
      If (KWord(1:4).eq.'REOR') Go To 800
      If (KWord(1:4).eq.'NODE') Go To 810
      If (KWord(1:4).eq.'DCHK') Go To 820
      If (KWord(1:4).eq.'TIME') Go To 830
      If (KWord(1:4).eq.'SCRN') Go To 840
      If (KWord(1:4).eq.'DMPK') Go To 850
      If (KWord(1:4).eq.'UPDA') Go To 860
      If (KWord(1:4).eq.'ESTI') Go To 870
      If (KWord(1:4).eq.'LOCK') Go To 880
      If (KWord(1:4).eq.'LK  ') Go To 880
      If (KWord(1:4).eq.'NOLK') Go To 881
      If (KWord(1:4).eq.'MEMF') Go To 890
      If (KWord(1:4).eq.'PRIN') Go To 700
      If (KWord(1:4).eq.'ENDC') Go To 998
      If (KWord(1:4).eq.'END ') Go To 998
      If (KWord(1:4).eq.'ENDO') Go To 998

!-------------------------------------------------------------------*
! Control section
!-------------------------------------------------------------------*
      iChrct=Len(KWord)
      Last=iCLast(KWord,iChrct)
      Write (6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      Write (6,*) SECNAM, ' Error in keyword.'
      Call Quit_OnUserError()
!                                                                      *
!***** ALGO ************************************************************
!                                                                      *
!-----Read Cholesky algorithm parameters
!
 900  Continue
!
       READ(LuSpool,*) ALGO
!
#ifdef _DEBUGPRINT_
      if(ALGO.eq.0)then
      Write(6,*) 'Integral regeneration from Cholesky vectors reordered on disk'
      Write(6,*)
      elseif(ALGO.eq.1)then
      Write(6,*) 'Density-based Cholesky. Default reorder: on the fly'
      Write(6,*)
      elseif(ALGO.eq.2)then
      Write(6,*) 'MO-based-Exchange Cholesky. Default reorder: on the fly'
      Write(6,*)
      elseif(ALGO.eq.3)then
      Write(6,*) 'MO-based-Exchange Cholesky. MO-transformation in reduced sets'
      Write(6,*)
      elseif(ALGO.eq.4)then
      Write(6,*) 'Local-Exchange (LK) algorithm.'
      Write(6,*)
      endif
#endif

      Go To 1000
!                                                                      *
!***** REOR ************************************************************
!                                                                      *
 800   Continue
       REORD=.true.
#ifdef _DEBUGPRINT_
      Write(6,*) 'Vectors reordered on DISK'
      Write(6,*)
#endif
!
      Go To 1000
!                                                                      *
!***** NODE ************************************************************
!                                                                      *
 810   Continue
       DECO=.false.
#ifdef _DEBUGPRINT_
      Write(6,*) 'Not-Using Decomposed density matrix'
      Write(6,*)
#endif
!
      Go To 1000
!                                                                      *
!***** DCHK ************************************************************
!                                                                      *
 820   Continue
       DensityCheck=.true.
!
      Go To 1000
!                                                                      *
!***** TIME ************************************************************
!                                                                      *
 830   Continue
       timings=.true.
!
      Go To 1000
!                                                                      *
!***** SCRN ************************************************************
!                                                                      *
 840   Continue
       READ(LuSpool,*) NSCREEN
!
      Go To 1000
!***** DMPK ************************************************************
!                                                                      *
 850   Continue
       READ(LuSpool,*) dmpk
       If (dmpk .lt. Zero) Then
        write(6,*)'OBS! Specified Negative DMPK value. Restore Defaults'
        dmpk=dmpk_dfl
       EndIf
!
      Go To 1000
!                                                                      *
!***** UPDA ************************************************************
!                                                                      *
 860   Continue
       Update=.true.
#ifdef _DEBUGPRINT_
      Write(6,*) 'Local-K with updating of the true diagonals'
      Write(6,*)
#endif
!
      Go To 1000
!                                                                      *
!***** ESTI ************************************************************
!                                                                      *
 870   Continue
       Estimate=.true.
#ifdef _DEBUGPRINT_
      Write(6,*) 'Local-K with evaluation of the diagonals from the current vec '
      Write(6,*)
#endif
!
      Go To 1000
!                                                                      *
!***** LOCK or LK ******************************************************
!                                                                      *
 880   Continue
       algo=4
#ifdef _DEBUGPRINT_
       Write(6,*) 'Local-Exchange (LK) algorithm.'
       Write(6,*)
#endif
!
      Go To 1000
!                                                                      *
!***** NoLK ************************************************************
!                                                                      *
 881   Continue
       algo=3
#ifdef _DEBUGPRINT_
       Write(6,*) 'Local-Exchange (LK) screening turned off! '
       Write(6,*)
#endif
!
      Go To 1000
!                                                                      *
!***** MemF ************************************************************
!                                                                      *
 890   Continue
       READ(LuSpool,*) ChFracMem
!
      Go To 1000
!                                                                      *
!***** PRIN ************************************************************
!                                                                      *
!-----Print level

 700  Key=Get_Ln(LuSpool)
      KWord=Key
      Call Get_I1(1,n)
      Do i = 1, n
         KWord=Get_Ln(LuSpool)
         Call Get_I1(1,jRout)
         Call Get_I1(2,iPrint)
         nPrint(jRout)=iPrint
      End Do
      Go To 1000
!                                                                      *
!***** ENDOFchoinput  **************************************************
!                                                                      *
!-----EndofChoinput
!
 998  Continue
!                                                                      *
!***********************************************************************
!                                                                      *
999   Continue
      Return
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Error handling
!
      Write (6,*) SECNAM, ' Premature end of input file.'
      Call Quit_OnUserError()
      Write (6,*) SECNAM, ' Error while reading input file.'
      Call Quit_OnUserError()
!
      End SubRoutine CHO_SCF_RDINP
