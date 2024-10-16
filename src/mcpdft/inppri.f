************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine InpPri_m()
************************************************************************
*                                                                      *
*     Echo input                                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use OneDat, only: sNoOri
      Use Functionals, only: Init_Funcs, Print_Info
      Use KSDFT_Info, only: CoefR, CoefX
      use printlevel, only: silent, usual
      use mcpdft_output, only: lf, iPrLoc
      use Fock_util_global, only: docholesky
      use rctfld_module, only: lRF
      use mcpdft_input, only: mcpdft_options
      use stdalloc, only: mma_allocate, mma_deallocate

      Implicit Real*8 (A-H,O-Z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
      Character(LEN=8)   Fmt1,Fmt2, Label
      Character(LEN=120)  Line,BlLine,StLine
      Character(LEN=3) lIrrep(8)

      Real*8, Allocatable:: Tmp0(:)

* Print level:
      IPRLEV=IPRLOC(1)
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
      lPaper=132
      Zero = 0.0D0
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      lLine=Len(Line)
      Do i=1,lLine
        BlLine(i:i)=' '
        StLine(i:i)='*'
      End Do
      lPaper=132
      left=(lPaper-lLine)/2
      Write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
      IF (IPRLEV.EQ.SILENT) GOTO 900
*----------------------------------------------------------------------*
*     Print the project title                                          *
*----------------------------------------------------------------------*
      IF(IPRLEV >= USUAL) THEN
       If ( nTit > 0 ) then
         Write(LF,*)
         nLine=nTit+5
         Do i=1,nLine
           Line=BlLine
           If ( i.eq.1 .or. i.eq.nLine )
     &     Line=StLine
           If ( i.eq.3 )
     &     Line='Project:'
           If ( i.ge.4 .and. i.le.nLine-2 )
     &     Write(Line,'(A72)')Title(i-3)
           Call Center_Text(Line)
           Write(LF,Fmt1) '*'//Line//'*'
         End Do
       Write(LF,*)
       End If
      END IF
*----------------------------------------------------------------------*
*     Print the ONEINT file identifier                                 *
*----------------------------------------------------------------------*
      IF(IPRLEV.GE.USUAL) THEN
       Write(LF,*)
       Write(LF,Fmt1) 'Header of the ONEINT file:'
       Write(LF,Fmt1) '--------------------------'
       Write(Line,'(36A2)') (Header(i),i=1,36)
       Write(LF,Fmt1)  trim(adjustl(Line))
       Write(Line,'(36A2)') (Header(i),i=37,72)
       Write(LF,Fmt1)  trim(adjustl(Line))
       Write(LF,*)
*----------------------------------------------------------------------*
*     Print the status of ORDINT                                       *
*----------------------------------------------------------------------*
       Write(LF,*)
       If (lSquare) Then
         Write(LF,Fmt1) 'OrdInt status: squared'
       Else
         Write(LF,Fmt1) 'OrdInt status: non-squared'
       End If
       Write(LF,*)
!----------------------------------------------------------------------*
!     Print cartesian coordinates of the system                        *
!----------------------------------------------------------------------*
       Call PrCoor
      END IF
*----------------------------------------------------------------------*
*     Print orbital and wavefunction specifications                    *
*----------------------------------------------------------------------*
      IF(IPRLEV.GE.USUAL) THEN
      Write(LF,*)
      Line=' '
      Write(Line(left-2:),'(A)') 'Wave function specifications:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt1)'-----------------------------'
      Write(LF,*)
      If (NFR.gt.0)
     &Write(LF,Fmt2//'A,T45,I6)')'Number of frozen shell electrons',
     &                           2*NFR
      Write(LF,Fmt2//'A,T45,I6)')'Number of closed shell electrons',
     &                           2*NIN
      Write(LF,Fmt2//'A,T45,I6)')'Number of electrons in active shells',
     &                           NACTEL
C.. for RAS
      Write(LF,Fmt2//'A,T45,I6)')'Max number of holes in RAS1 space',
     &                           NHOLE1
      Write(LF,Fmt2//'A,T45,I6)')'Max nr of electrons in RAS3 space',
     &                           NELEC3


      If (NFR.gt.0)
     &Write(LF,Fmt2//'A,T45,I6)')'Number of frozen orbitals',
     &                           NFR
      Write(LF,Fmt2//'A,T45,I6)')'Number of inactive orbitals',
     &                           NIN
      Write(LF,Fmt2//'A,T45,I6)')'Number of active orbitals',
     &                           NAC
      Write(LF,Fmt2//'A,T45,I6)')'Number of secondary orbitals',
     &                           NSEC
      Write(LF,Fmt2//'A,T45,F6.1)')'Spin quantum number',
     &                           (DBLE(ISPIN-1))/2.0d0
      Write(LF,Fmt2//'A,T45,I6)')'State symmetry',
     &                           STSYM
      Call CollapseOutput(0,'Wave function specifications:')
*
      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym = 1, nSym
         lIrrep(iSym) = adjustr(lIrrep(iSym))
      End Do
*
      Write(LF,*)
      Line=' '
      Write(Line(left-2:),'(A)') 'Orbital specifications:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt1)'-----------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A,T47,8I4)') 'Symmetry species',
     &                            (iSym,iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8(1X,A))') '                ',
     &                            (lIrrep(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Frozen orbitals',
     &                            (nFro(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Inactive orbitals',
     &                            (nIsh(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Active orbitals',
     &                            (nAsh(iSym),iSym=1,nSym)
        Write(LF,Fmt2//'A,T47,8I4)') 'RAS1 orbitals',
     &                            (nRs1(iSym),iSym=1,nSym)
        Write(LF,Fmt2//'A,T47,8I4)') 'RAS2 orbitals',
     &                            (nRs2(iSym),iSym=1,nSym)
        Write(LF,Fmt2//'A,T47,8I4)') 'RAS3 orbitals',
     &                            (nRs3(iSym),iSym=1,nSym)

      Write(LF,Fmt2//'A,T47,8I4)') 'Secondary orbitals',
     &                            (nSsh(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Deleted orbitals',
     &                            (nDel(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Number of basis functions',
     &                            (nBas(iSym),iSym=1,nSym)
      Call CollapseOutput(0,'Orbital specifications:')
      Write(LF,*)

      Write(LF,Fmt2//'A,T45,I6)')'Number of root(s) required',
     &                             NROOTS
*TRS
      Call Get_iScalar('Relax CASSCF root',iRlxRoot)
*TRS
      If (irlxroot.ne.0)
     &       Write(LF,Fmt2//'A,T45,I6)')'Root chosen for geometry opt.',
     &                             IRLXROOT

      Write(LF,Fmt2//'A,T45,I6)')'highest root included in the CI',
     &                           LROOTS

      Call CollapseOutput(0,'CI expansion specifications:')

      END IF

      IF(IPRLEV.GE.USUAL) THEN
       Write(LF,*)
       Line=' '
       Write(Line(left-2:),'(A)') 'Optimization specifications:'
       Call CollapseOutput(1,Line)
       Write(LF,Fmt1)'----------------------------'
       Write(LF,*)
       If (DoCholesky) Then
        Call Get_iScalar('System BitSwitch',iDoRI)
        if (Iand(iDoRI,1024).Eq.1024) then
             Write(LF,Fmt2//'A,T45,I6)')'RASSCF algorithm: LK RI/DF'

        else
             Write(LF,Fmt2//'A,T45,I6)')'RASSCF algorithm: LK Cholesky'
        endif
       Else
        Write(LF,Fmt2//'A,T45,I6)')'RASSCF algorithm: Conventional'
       EndIf
        Write(LF,Fmt2//'A)') 'This is a MC-PDFT calculation '//
     &   'with functional: '//mcpdft_options%otfnal%otxc
        Write(LF,Fmt2//'A,T45,ES10.3)')'Exchange scaling factor',CoefX
        Write(LF,Fmt2//'A,T45,ES10.3)')'Correlation scaling factor',
     &                                 CoefR
       If (mcpdft_options%grad) then
        Write(LF,Fmt1) 'Potentials are computed for gradients'
       end if
       If ( lRF ) then
         Call mma_allocate(Tmp0,nTot1+4,Label='Tmp0')
         iRc=-1
         iOpt=ibset(0,sNoOri)
         iComp=1
         iSyLbl=1
         Label='Mltpl  0'
         Call RdOne(iRc,iOpt,Label,iComp,Tmp0,iSyLbl)
         Tot_Nuc_Charge=Tmp0(nTot1+4)
         If ( iRc.ne.0 ) then
            Write(LF,*) 'InpPri: iRc from Call RdOne not 0'
            Write(LF,*) 'Label = ',Label
            Write(LF,*) 'iRc = ',iRc
            Call Abend
         Endif
         Call mma_deallocate(Tmp0)
         Tot_El_Charge=Zero
         Do iSym=1,nSym
            Tot_El_Charge=Tot_El_Charge
     &                   -2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
         End Do
         Tot_El_Charge=Tot_El_Charge-DBLE(nActEl)
         Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
         iCharge=Int(Tot_Charge)
         Call PrRF(.False.,NonEq,iCharge,2)
       End If
       Call CollapseOutput(0,'Optimization specifications:')
       Write(LF,Fmt1)
     &  'Starting CI array(s) will be read from file'
      END IF
      Write(LF,*)

!---- Print out grid information
       Call Put_dScalar('DFT exch coeff',CoefX)
       Call Put_dScalar('DFT corr coeff',CoefR)
       Call Funi_Print()
       IF(IPRLEV.GE.USUAL) THEN
          Write(6,*)
          Write(6,'(6X,A)') 'DFT functional specifications'
          Write(6,'(6X,A)') '-----------------------------'
          Call libxc_version()
          Call Init_Funcs(mcpdft_options%otfnal%xc)
          Call Print_Info()
          Write(6,*)
       END IF
  900 CONTINUE
      Call XFlush(LF)

      Return
      End
