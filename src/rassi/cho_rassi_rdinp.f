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
      SubRoutine CHO_RASSI_RDINP(DFonly,LuSpool)
************************************************************************
*
*  Purpose:   If DFonly, use defaults only.
*             Else, read and process input for Cholesky section
*             in RASSI
*
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "rassi.fh"
#include "print.fh"
      Character*180 KWord, Key, Get_Ln
      External Get_Ln
      Logical  DFonly,timings
      Logical  Update,Estimate,Deco,PseudoChoMOs
      character*15 SECNAM
      parameter (SECNAM = 'CHO_RASSI_RDINP')
      Integer  ALGO,Nscreen
      Real*8   dmpk
*
      Common /CHORASSI / ALGO,Nscreen,dmpk
      COMMON /CHOTIME / timings
      COMMON /LKSCREEN / Estimate, Update, Deco, PseudoChoMOs
*
      iRout=1
      iPrint=nPrint(iRout)
*
***** Algorithms for using Cholesky vectors in RASSI ******************
*
*   ALGO:
*
*          1  --->  compute fock matrices in AO-basis from vectors
*                   read in reduced sets and transformed on the fly.
*                   (TU|VX) integrals are also returned
*
*          2  --->  compute fock matrices in AO-basis from vectors
*                   read in reduced sets and transformed on the fly
*                   (TU|VX) integrals are also returned.
*                   Exchange contributions are computed using the
*                   "Local K" scheme
*                                                                      *
************************************************************************
*                                                                      *
*     Default  parameters
#if defined (_MOLCAS_MPP_)
      ChFracMem=0.3d0
#else
      ChFracMem=0.0d0
#endif
      IF (DFonly) THEN
         ALGO  = 2
         timings=.false.
         dmpk=1.0d-1
         Nscreen=10
         Deco=.true.
         Update=.true.
         Estimate=.false.
         goto 999  !return flag
      ENDIF
*
*    set some parameters if not specified in ChoInput section
         ALGO  = 2
         timings=.false.
         dmpk=1.0d-1
         Nscreen=10  ! to be reset to 0 if doing RI calculation
         Deco=.true.
         Update=.true.
         Estimate=.false.
         PseudoChoMOs=.false.

         dmpk_dfl=1.0d-1
************************************************************************
*                                                                      *
      iPrint=5
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
      If (KWord.eq.'')       Go To 1000
      If (KWord(1:4).eq.'ALGO') Go To 900
      If (KWord(1:4).eq.'LOCK') Go To 810
      If (KWord(1:4).eq.'NOLK') Go To 811
      If (KWord(1:4).eq.'SCRN') Go To 820
      If (KWord(1:4).eq.'DMPK') Go To 825
      If (KWord(1:4).eq.'TIME') Go To 830
      If (KWord(1:4).eq.'UPDA') Go To 840
      If (KWord(1:4).eq.'ESTI') Go To 850
      If (KWord(1:4).eq.'MEMF') Go To 860
      If (KWord(1:4).eq.'NODE') Go To 870
      If (KWord(1:4).eq.'PSEU') Go To 880
      If (KWord(1:4).eq.'PRIN') Go To 700
      If (KWord(1:4).eq.'ENDC') Go To 998
      If (KWord(1:4).eq.'END ') Go To 998
      If (KWord(1:4).eq.'ENDO') Go To 998

*-------------------------------------------------------------------*
* Control section
*-------------------------------------------------------------------*
      iChrct=Len(KWord)
      Last=iCLast(KWord,iChrct)
      Write (6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      Call ErrTra
      Write (6,*) SECNAM, ' Error in keyword.'
      Call AbEnd()
*                                                                      *
****** ALGO ************************************************************
*                                                                      *
*-----Read Cholesky algorithm parameters
*
 900  Continue
c      Call Get_F1(1,Eps)
c      Call Get_F1(2,rds)
c      Call Get_I1(1,ALGO)
C      If (nToken(KWord).gt.1) goto 988
*
       READ(LuSpool,*) ALGO
*
      if(ALGO.eq.1)then
      Write(6,*)
     &'Default CD-RASSI algorithm reset to  ',ALGO
      Write(6,*)
      elseif(ALGO.eq.2)then
      Write(6,*)
     &'Default CD-RASSI algorithm reset to  ',ALGO
      Write(6,*)
      else
      Write(6,*)
     &'The specified algorithm is not implemented. Option Ignored '
      Write(6,*)
      endif
*
      Go To 1000
*                                                                      *
****** LOCK ************************************************************
*                                                                      *
 810   Continue
       algo=2
      Write(6,*)
     &'Default CD-RASSI algorithm reset to  ',ALGO
      Write(6,*)
     &'Using Local K scheme for Exchange matrices '

*
      Go To 1000
*                                                                      *
*                                                                      *
****** NOLK ************************************************************
*                                                                      *
 811   Continue
       algo=1
      Write(6,*)
     &'Default CD-RASSI algorithm reset to  ',ALGO
      Write(6,*)
     &'Local K scheme for Exchange matrices turned off ! '

*
      Go To 1000
*                                                                      *
*                                                                      *
****** SCRN ************************************************************
*                                                                      *
 820   Continue
       READ(LuSpool,*) Nscreen
*
      Go To 1000
*                                                                      *
****** DMPK ************************************************************
*                                                                      *
 825   Continue
       READ(LuSpool,*) DMPK
       If (dmpk .lt. 0.0d0) Then
        write(6,*)'OBS! Specified Negative DMPK value. Restore Defaults'
        dmpk=dmpk_dfl
      EndIf
*
      Go To 1000
*                                                                      *
****** TIME ************************************************************
*                                                                      *
 830   Continue
       timings=.true.
*
      Go To 1000
*                                                                      *
****** UPDA ************************************************************
*                                                                      *
 840   Continue
       Update=.true.
      Write(6,*)
     &'Local-K with updating of the true diagonals'
      Write(6,*)
*
      Go To 1000
*                                                                      *
****** ESTI ************************************************************
*                                                                      *
 850   Continue
       Estimate=.true.
      Write(6,*)
     &'Local-K with evaluation of the diagonals from the current vec '
      Write(6,*)
*
      Go To 1000
*                                                                      *
****** MEMF ************************************************************
*                                                                      *
 860   Continue
       READ(LuSpool,*) ChFracMem
*
      Go To 1000
*                                                                      *
****** NODE ************************************************************
*                                                                      *
 870   Continue
       If (PseudoChoMOs) Then
          Write(6,*) ' The keyword NODEcompose is incompatible with'//
     &               ' the previously specified keyword PSEUdo.'//
     &               ' NODEcompose will be ignored'
       Else
          Deco=.false.
          Write(6,*)'Canonical inactive orbitals used in LK CD-RASSI.'
       EndIf
       Write(6,*)
*
      Go To 1000
*                                                                      *
****** PSEU ************************************************************
*                                                                      *
 880   Continue
       If (.not.Deco) Then
          Write(6,*) ' The keyword PSEUdo is incompatible with'//
     &               ' the previously specified keyword NODEco'//
     &               'mpose. PSEUdo will be ignored'
       Else
          PseudoChoMOs=.true.
          Write(6,*)'Pseudo Cholesky orbitals used in LK CD-RASSI.'
       EndIf
       Write(6,*)
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
999   Return
*                                                                      *
************************************************************************
*                                                                      *
*-----Error handling
*
*977  Call ErrTra
*     Write (6,*) SECNAM, ' Premature end of input file.'
*     Call AbEnd()
*988  Call ErrTra
*     Write (6,*) SECNAM, ' Error while reading input file.'
*     Call AbEnd()
*
      End
