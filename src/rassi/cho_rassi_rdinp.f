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
      Use Fock_util_global, only: Deco, Estimate, PseudoChoMOs, Update
      Use Cholesky, only: timings
      Implicit Real*8 (A-H,O-Z)
#include "rassi.fh"
#include "print.fh"
      Character(len=180) KWord, Key, Get_Ln
      External Get_Ln
      Logical  DFonly
*
#include "chorassi.fh"
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
#ifdef _MOLCAS_MPP_
      ChFracMem=0.3d0
#else
      ChFracMem=0.0d0
#endif
      ALGO  = 2
      timings=.false.
      dmpk=1.0d-1
      Nscreen=10  ! to be reset to 0 if doing RI calculation
      Deco=.true.
      Update=.true.
      Estimate=.false.

      IF (.not. DFonly) THEN
*
*    set some parameters if not specified in ChoInput section
        PseudoChoMOs=.false.
        dmpk_dfl=1.0d-1
        iPrint=5
*                                                                      *
************************************************************************
*                                                                      *
*-----Process the input
*
*-------------------------------------------------------------------*
      do
*-------------------------------------------------------------------*
* Use Get_Ln to read the lines.                                     *
*-------------------------------------------------------------------*
      Key=Get_Ln(LuSpool)
      Kword=Key
      Call UpCase(Kword)
*-------------------------------------------------------------------*
* The keywords and their labels.                                    *
*-------------------------------------------------------------------*
      If (KWord(1:1).eq.'*') cycle
      If (KWord.eq.'') cycle
      select case (Kword(1:4))
      case ('ALGO')
****** ALGO ************************************************************
*                                                                      *
*-----Read Cholesky algorithm parameters
*
c      Call Get_F1(1,Eps)
c      Call Get_F1(2,rds)
c      Call Get_I1(1,ALGO)
C      If (nToken(KWord).gt.1) call abend()
*
       READ(LuSpool,*) ALGO
      if(ALGO.eq.1)then
        Write(6,*) 'Default CD-RASSI algorithm reset to  ',ALGO
        Write(6,*)
      elseif (ALGO.eq.2) then
        Write(6,*) 'Default CD-RASSI algorithm reset to  ',ALGO
        Write(6,*)
      else
        Write(6,*) 'The specified algorithm is not implemented.'
        Write(6,*) 'Option Ignored '
        Write(6,*)
      endif
*
      case ('LOCK')
        algo=2
        Write(6,*) 'Default CD-RASSI algorithm reset to  ',ALGO
        Write(6,*) 'Using Local K scheme for Exchange matrices '

      case ('NOLK')
        algo=1
        Write(6,*) 'Default CD-RASSI algorithm reset to  ',ALGO
        Write(6,*) 'Local K scheme for Exchange matrices turned off ! '

      case ('SCRN')
       READ(LuSpool,*) Nscreen

      case ('DMPK')
      READ(LuSpool,*) DMPK
      If (dmpk .lt. 0.0d0) Then
      write(6,*)'OBS! Specified Negative DMPK value. Restore Defaults'
      dmpk=dmpk_dfl
      End If

      case ('TIME')
            timings=.true.

      case ('UPDA')
        Update=.true.
        Write(6,*) 'Local-K with updating of the true diagonals'
        Write(6,*)

      case ('ESTI')
        Estimate=.true.
        Write(6,*)
     &'Local-K with evaluation of the diagonals from the current vec '
        Write(6,*)
      case ('MEMF')
        read(LuSpool,*) ChFracMem

      case ('NODE')
        If (PseudoChoMOs) Then
          Write(6,*) ' The keyword NODEcompose is incompatible with'//
     &               ' the previously specified keyword PSEUdo.'//
     &               ' NODEcompose will be ignored'
        Else
          Deco=.false.
          Write(6,*)'Canonical inactive orbitals used in LK CD-RASSI.'
        End If
        Write(6,*)

      case ('PSEU')
        If (.not.Deco) Then
          Write(6,*) ' The keyword PSEUdo is incompatible with'//
     &               ' the previously specified keyword NODEco'//
     &               'mpose. PSEUdo will be ignored'
        Else
          PseudoChoMOs=.true.
          Write(6,*)'Pseudo Cholesky orbitals used in LK CD-RASSI.'
        EndIf
      Write(6,*)

      case ('PRIN')
        Key=Get_Ln(LuSpool)
        KWord=Key
        Call Get_I1(1,n)
        Do i = 1, n
          KWord=Get_Ln(LuSpool)
          Call Get_I1(1,jRout)
          Call Get_I1(2,iPrint)
          nPrint(jRout)=iPrint
        End Do

      case ('ENDC')
        exit
      case ('END ')
        exit
      case ('ENDO')
        exit

      case default
        iChrct=Len(KWord)
        Last=iCLast(KWord,iChrct)
        Write (6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
        Write (6,*) 'CHO_RASSI_RDINP: Error in keyword.'
        Call AbEnd()

      end select

      end do

      end if

      end subroutine CHO_RASSI_RDINP
