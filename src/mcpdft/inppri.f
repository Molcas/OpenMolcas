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
      use mspdft, only: dogradmspd
      use mcpdft_output, only: silent, usual, lf, iPrLoc
      use Fock_util_global, only: docholesky

      Implicit Real*8 (A-H,O-Z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "ciinfo.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
      Character*8   Fmt1,Fmt2, Label
      Character*120  Line,BlLine,StLine
      Character*3 lIrrep(8)
      Character*80 KSDFT2

* Print level:
      IPRLEV=IPRLOC(1)
      DoDMRG=.False.
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
*----------------------------------------------------------------------*
*     Print cartesian coordinates of the system                        *
*----------------------------------------------------------------------*
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

      If(.Not.DoDMRG) GoTo 113

      Line=' '
      Write(Line(left-2:),'(A)') 'DMRG sweep specifications:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt1)'--------------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A,T45,I6)')'Number of renormalized basis',
     &                           MxDMRG
      Write(LF,Fmt2//'A,T45,I6)')'Number of root(s) required',
     &                           NROOTS
* NN.14 FIXME: haven't yet checked whether geometry opt. works correctly with DMRG
      Write(LF,Fmt2//'A,T45,I6)')'Root chosen for geometry opt.',
     &                           IRLXROOT
      Call CollapseOutput(0,'DMRG sweep specifications:')

*     Skip printing CI specifications in DMRG-CASSCF
      GoTo 114

 113  Continue

      Line=' '
      Write(Line(left-2:),'(A)') 'CI expansion specifications:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt1)'----------------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A,T40,I11)')'Number of CSFs',
     &                           NCSASM(STSYM)
      Write(LF,Fmt2//'A,T40,I11)')'Number of determinants',
     &                           NDTASM(STSYM)
        n_Det=2
        n_unpaired_elec=(iSpin-1)
        n_paired_elec=nActEl-n_unpaired_elec
        If(n_unpaired_elec+n_paired_elec/2.eq.nac.or.
     &     NDTASM(STSYM).eq.1) n_Det = 1
        If(KSDFT.eq.'DIFF')   n_Det = 1
        If(KSDFT.eq.'ROKS')   n_Det = 1

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

 114  Continue

      END IF

* Check that the user doesn't try to calculate more roots than it's possible
* NN.14 FIXME: in DMRG-CASSCF, skip this check for the time
*              since Block DMRG code will check this internally
*     If (NROOTS .GT. NCSASM(STSYM)) Then
      If (.false.) Then
!      If (.NOT.DoDMRG .AND. NROOTS .GT. NCSASM(STSYM)) Then
         Write(LF,*) '************ ERROR ***********'
         Write(LF,*) ' You can''t ask for more roots'
         Write(LF,*) ' than there are configurations '
         Write(LF,*) '******************************'
         Write(LF,*)
!         Call Quit_OnUserError()
      End If
* If the calculation will be too big:
      call GetMem('ChkMx','Max','Real',iDum,MaxRem)
      WillNeedMB=(8.0D0*1.50D0*6.0D0*NDTASM(STSYM)/1.048D6)
      AvailMB=(8.0D0*MaxRem/1.048D6)
      if (WillNeedMB .gt. AvailMB) then
        write(6,*)
        write(6,*)' *************************************************'
        write(6,*)' Sorry, but your calculation will probably be too'
        write(6,*)' large for the available memory.'
        write(6,*)' The number of determinants is ',NDTASM(STSYM)
        write(6,*)' During CI equation solution, there will be'
        write(6,*)' up to six vectors of this size in memory.'
        write(6,*)' We estimate an additional 50% for other stuff.'
        write(6,*)
        write(6,'(A,F9.1,A)')' This alone will need at least ',
     &                                            WillNeedMB,' MB.'
        write(6,'(A,F9.1,A)')' Available at this point is ',
     &                                               AvailMB,' MB.'
        write(6,*)' Please increase MOLCAS_MEM, and try again.'
        write(6,*)' *************************************************'
        write(6,*)
        Call Quit_OnUserError()
      end if

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
       KSDFT2 = KSDFT
       IF(KSDFT(1:2).eq.'T:'.or.KSDFT(1:3).eq.'FT:') Then
        KSDFT2 = KSDFT(index(KSDFT,'T:')+2:)
        Write(LF,Fmt2//'A)') 'This is a MC-PDFT calculation '//
     &   'with functional: '//KSDFT
        Write(LF,Fmt2//'A,T45,E10.3)')'Exchange scaling factor',CoefX
        Write(LF,Fmt2//'A,T45,E10.3)')'Correlation scaling factor',
     &                                 CoefR
       end if
       If (dogradPDFT.or.dogradMSPD) then
        Write(LF,Fmt1) 'Potentials are computed for gradients'
       end if
       If ( lRF ) then
         Call GetMem('Ovrlp','Allo','Real',iTmp0,nTot1+4)
         iRc=-1
         iOpt=ibset(0,sNoOri)
         iComp=1
         iSyLbl=1
         Label='Mltpl  0'
         Call RdOne(iRc,iOpt,Label,iComp,Work(iTmp0),iSyLbl)
         Tot_Nuc_Charge=Work(iTmp0+nTot1+3)
         If ( iRc.ne.0 ) then
            Write(LF,*) 'InpPri: iRc from Call RdOne not 0'
            Write(LF,*) 'Label = ',Label
            Write(LF,*) 'iRc = ',iRc
            Call Abend
         Endif
         Call GetMem('Ovrlp','Free','Real',iTmp0,nTot1+4)
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

       Call Put_dScalar('EThr',ThrE)
*
*---- Print out grid information in case of DFT
*
       If (KSDFT.ne.'SCF') Then
         Call Put_dScalar('DFT exch coeff',CoefX)
         Call Put_dScalar('DFT corr coeff',CoefR)
         Call Funi_Print()
         IF(IPRLEV.GE.USUAL) THEN
            Write(6,*)
            Write(6,'(6X,A)') 'DFT functional specifications'
            Write(6,'(6X,A)') '-----------------------------'
            Call libxc_version()
            Call Init_Funcs(KSDFT2)
            Call Print_Info()
            Write(6,*)
         END IF
       End If
  900 CONTINUE
      Call XFlush(LF)
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
