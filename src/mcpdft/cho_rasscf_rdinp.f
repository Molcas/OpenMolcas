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
      SubRoutine CHO_RASSCF_RDINP_m(DFonly,LuInput)
************************************************************************
*                                                                      *
*  Purpose:   If DFonly, use defaults only.                            *
*             Else, read and process input for Cholesky section        *
*             in RASSCF                                                *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='CHO_RASS')
#include "WrkSpc.fh"
      Character*180 KWord, Key, Get_Ln
      External Get_Ln
      Logical  DFonly,DensityCheck,timings,DoLock,Deco
      Logical  DoCholesky,Estimate,Update
      character*16 SECNAM
      parameter (SECNAM = 'CHO_RASSCF_RDINP')
      Integer  ALGO,Nscreen
      Real*8   dmpk
*
      Common /CHLCAS / DoCholesky,ALGO
      COMMON /CHODENSITY/ DensityCheck
      COMMON /CHOTIME / timings
      Common /CHOLK / DoLocK,Deco,dmpk,Nscreen
      COMMON /CHOSCREEN/ Estimate,Update
      COMMON /CHOPAR/ ChFracMem

*
*
***** Algorithms for using Cholesky vectors in RASSCF ******************
*
*   ALGO:
*
*          1  --->  compute fock matrices in AO-basis from vectors
*                   read in reduced sets and transformed on the fly
*                   (PU|VX) integrals are also returned
*
*                   If DoLock=.true. then the Exchange terms are
*                   computed by emplpying the "Local K" scheme
*
*
*          2  --->  compute fock matrices in AO-basis from vectors
*                   read in reduced sets and transformed on the fly
*                   Only the (TU|VX) integrals are returned
*                   while the vectors are directly contracted
*                   with the 2-body density in order to construct
*                   the Q-matrix
*                                                                      *
************************************************************************
*                                                                      *
*     Default  parameters
#if defined (_MOLCAS_MPP_)
      ChFracMem=0.3D0
#else
      ChFracMem=0.0D0
#endif
*
*     set some parameters if not specified in ChoInput section
      ALGO  = 1
      DensityCheck=.false.
      Deco=.true.
      timings=.false.
      DoLock=.true.
      Nscreen=10
      dmpk=1.0d-1
      Update=.true.
      Estimate=.false.
      IF (DFonly) goto 999  !return flag

      dmpk_dfl=1.0d-1
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
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
      Key=Get_Ln(LuInput)
      Kword=Key
      Call UpCase(Kword)
*-------------------------------------------------------------------*
* The keywords and their labels.                                    *
*-------------------------------------------------------------------*

      If (KWord(1:1).eq.'*')    Go To 1000
      If (KWord.eq.'')       Go To 1000
      If (KWord(1:4).eq.'ALGO') Go To 900
      If (KWord(1:4).eq.'LOCK') Go To 910
      If (KWord(1:4).eq.'LK  ') Go To 910
      If (KWord(1:4).eq.'NOLK') Go To 915
      If (KWord(1:4).eq.'DMPK') Go To 920
      If (KWord(1:4).eq.'NODE') Go To 930
      If (KWord(1:4).eq.'SCRN') Go To 940
      If (KWord(1:4).eq.'MEMF') Go To 950
      If (KWord(1:4).eq.'DCHK') Go To 820
      If (KWord(1:4).eq.'TIME') Go To 830
      If (KWord(1:4).eq.'ESTI') Go To 840
      If (KWord(1:4).eq.'UPDA') Go To 850
      If (KWord(1:4).eq.'PRIN') Go To 700
      If (KWord(1:4).eq.'ENDC') Go To 998
      If (KWord(1:4).eq.'END ') Go To 998
      If (KWord(1:4).eq.'ENDO') Go To 998

*-------------------------------------------------------------------*
* Control section
*-------------------------------------------------------------------*
      iChrct=Len(KWord)
      Last=iCLast(KWord,iChrct)
      Write(LF,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      Call ErrTra
      Write(LF,*) SECNAM, ' Error in keyword.'
      Call Quit_OnUserError()
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
       READ(LuInput,*) ALGO
*
      if(ALGO.eq.1)then
      Write(LF,*)
     &'Default RASSCF algorithm reset to  ',ALGO
      Write(LF,*)
      elseif(ALGO.eq.2)then
      Write(LF,*)
     &'Default RASSCF algorithm reset to  ',ALGO
      Write(LF,*)
      Write(LF,*)' !!! STILL UNDER DEBUGGING !!! '
      else
      Write(LF,*)
     &'The specified algorithm is not implemented. Option Ignored '
      Write(LF,*)
      endif
*
      Go To 1000
*                                                                      *
****** LOCK or LK ******************************************************
*                                                                      *
 910   Continue
      DoLocK=.true.
c      Write(LF,*)
c     &'Using Local K scheme for Exchange matrices '
*
      Go To 1000
*                                                                      *
****** NOLK ************************************************************
*                                                                      *
 915   Continue
      DoLocK=.false.
c      Write(LF,*)
c     &'LK screening for Exchange matrices turned off !'
*
      Go To 1000
*                                                                      *
****** DMPK ************************************************************
*                                                                      *
 920   Continue
       READ(LuInput,*) dmpk
       If (dmpk .lt. 0.0D0) Then
        write(6,*)'OBS! Specified Negative DMPK value. Restore Defaults'
        dmpk=dmpk_dfl
       EndIf
*
      Go To 1000
*                                                                      *
****** NODE ************************************************************
*                                                                      *
 930   Continue
      Deco=.false.
      Write(LF,*)
     &'Not-Using Cholesky decomposed Inactive density '
*
      Go To 1000
*                                                                      *
****** SCRN ************************************************************
*                                                                      *
 940   Continue
       READ(LuInput,*) Nscreen
*
      Go To 1000
*                                                                      *
****** MEMF ************************************************************
*                                                                      *
 950   Continue
       READ(LuInput,*) ChFracMem
*
      Go To 1000
*                                                                      *
****** DCHK ************************************************************
*                                                                      *
 820   Continue
       DensityCheck=.true.
      Write(LF,*)
     &'Non-valid option. IGNORED !! '
*
      Go To 1000
*                                                                      *
****** ESTI ************************************************************
*                                                                      *
 840   Continue
       Estimate=.true.
      Write(LF,*)
     &'Diagonal integrals estimated from the current Cholesky vectors'
*
      Go To 1000
*                                                                      *
****** UPDA ************************************************************
*                                                                      *
 850   Continue
       Update=.true.
      Write(LF,*)
     &'Updating of the true diagonal integrals'
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
****** PRIN ************************************************************
*                                                                      *
*-----Print level
*
 700  Key=Get_Ln(LuInput)
      KWord=Key
      Call Get_I1(1,n)
      Do i = 1, n
         KWord=Get_Ln(LuInput)
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
999      Write(LF,'(1X,A,I4)')
c     &'Default Cholesky algorithm in RASSCF = ',ALGO
c       Write(LF,*)
       If(ALGO.eq.2)Then
         Write(LF,*)'Local K scheme not implemented for the chosen algo'
     &'rithm. LocK keyword ignored !'
         DoLocK=.false.
       EndIf

      Return
*
      End
