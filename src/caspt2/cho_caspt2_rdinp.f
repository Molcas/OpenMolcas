************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine CHO_CASPT2_RDINP(DFonly,LuSpool)
************************************************************************
*
*  Purpose:   If DFonly, use defaults only.
*             Else, read and process input for Cholesky section
*             in CASPt2
*
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "output.fh"
#include "WrkSpc.fh"
      Character(180) KWord, Key, Get_Ln
      External Get_Ln
      Logical  DFonly,REORD,DECO,timings,DensityCheck
      character(16) SECNAM
      parameter (SECNAM = 'CHO_CASPT2_RDINP')
      Integer  ALGO
*
      Common /CHORAS / REORD,DECO,ALGO
      Common /CHOTIME/ timings
      COMMON /CHODENSITY/ DensityCheck

#include "chocaspt2.fh"

*
      iRout=1
      iPrint=nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
*     Algorithms for generating MO integrals in CASPT2
*
*        iALGO :
*               0  --> MOLINT file is generated from the
*                      transformed Cholesky vectors. Both "Coulomb" and
*                      "Exchange(1,2)" integrals are computed and stored
*                      on disk
*
*               1  --> Only the "Exchange" integrals are computed and
*                      combined directly in order to compute the RHS
*                      of the caspt2 equations. The latter is then
*                      stored on disk. The AO Fock matrix is
*                      computed during the MO transformation of the
*                      vectors and it is stored on disk
*
************************************************************************
*
***** Algorithms for using Cholesky vectors in Fock matrix generation **
*
*   ALGO:
*          0  --->  Integrals are regenerated on the fly
*                   from a set of Cholesky vectors resorted on disk
*
*          1  --->  The resorted Cholesky vectors are used directly
*                   by the Fock matrix builder routines and contracted
*                   with the proper density matrices. Uses
*                   vectors resorted either on disk or on the fly
*
*          2  --->  As in option 1 but using the MO-basis transformed
*                   vectors for computing the exchange term
*
*                                                                      *
************************************************************************
*                                                                      *
*     Default  parameters

      IF (DFonly) THEN
         iAlGO = 1
         ALGO  = 2
         REORD =.false.
         DECO  =.true.
         DensityCheck=.false.
         timings=.false.
         goto 999  !return flag
      ENDIF
*
*    set some parameters if not specified in ChoInput section
         iAlGO = 1
         ALGO  = 2
         REORD =.false.
         DECO  =.true.
         DensityCheck=.false.
         timings=.false.

*                                                                      *
************************************************************************
*     Define Blank lines                                               *
*
      Do i = 1, 80
         BLine(i:i) = ' '
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*-----Process the input
*
*-------------------------------------------------------------------*
* The big turning point.                                            *
*-------------------------------------------------------------------*
1000  Continue
*-------------------------------------------------------------------*
* Use Get_Ln to read the lines.                                     *
*-------------------------------------------------------------------*
      Key=Get_Ln(LuSpool)
      Kword=Key
      Call UpCase(Kword)
*-------------------------------------------------------------------*
* The keywords and their labels.                                    *
*-------------------------------------------------------------------*

      If (KWord(1:1).eq.'*')    Go To 1000
      If (KWord.eq.BLine)       Go To 1000
      If (KWord(1:4).eq.'ALGO') Go To 900
      If (KWord(1:4).eq.'IALG') Go To 950
      If (KWord(1:4).eq.'REOR') Go To 800
      If (KWord(1:4).eq.'DECO') Go To 810
      If (KWord(1:4).eq.'TIME') Go To 820
      If (KWord(1:4).eq.'DCHK') Go To 830
      If (KWord(1:4).eq.'PRIN') Go To 700
      If (KWord(1:4).eq.'ENDC') Go To 998
      If (KWord(1:4).eq.'END ') Go To 998
      If (KWord(1:4).eq.'ENDO') Go To 998

*-------------------------------------------------------------------*
* Control section                                                   *
*-------------------------------------------------------------------*
      iChrct=Len(KWord)
      Last=iCLast(KWord,iChrct)
      WRITE(6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      Call ErrTra
      WRITE(6,*) SECNAM, ' Error in keyword.'
      CALL ABEND()
*                                                                      *
****** ALGO ************************************************************
*                                                                      *
*-----Read Cholesky algorithm parameters
*
 900  Continue
*
       READ(LuSpool,*) ALGO
*
*
      Go To 1000
*                                                                      *
****** IALG ************************************************************
*                                                                      *
 950   Continue
       READ(LuSpool,*) iALGO
*
      Go To 1000
*                                                                      *
****** REOR ************************************************************
*                                                                      *
 800   Continue
       REORD=.true.
      WRITE(6,*)
     &'Vectors reordered on FILE'
      WRITE(6,*)
*
      Go To 1000
*                                                                      *
****** DECO ************************************************************
*                                                                      *
 810   Continue
       DECO=.true.
      WRITE(6,*)
     &'Decomposed densty matrix'
      WRITE(6,*)
*
      Go To 1000
*                                                                      *
****** TIME ************************************************************
*                                                                      *
 820   Continue
       timings=.true.
*
      Go To 1000
*                                                                      *
****** DCHK ************************************************************
*                                                                      *
 830   Continue
       DensityCheck=.true.
*
      Go To 1000

*                                                                      *
****** PRIN ************************************************************
*                                                                      *
*-----Print level
*
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
*                                                                      *
****** END  ************************************************************
*                                                                      *
*-----End of input
*
 998  Continue
*                                                                      *
************************************************************************
*                                                                      *
999   Continue
* SB: this printout is misleading if one does not use the Cholesky
* approximation and superfluous otherwise, since the algorithms
* are not documented. The user ignores it anyway.
      ! If (IPRGLB.ge.TERSE) Then
      !    WRITE(6,'(1X,A,I4)') 'Cholesky algorithm in CASPT2 = ',iALGO
      !    WRITE(6,*)
      ! End If

      Return
*                                                                      *
************************************************************************
*                                                                      *
*-----Error handling
*
*
      End
