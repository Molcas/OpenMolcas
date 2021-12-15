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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2003, Valera Veryazov                                  *
*               2017, Roland Lindh                                     *
************************************************************************
      SubRoutine PrFin(OneHam,Ovlp,Dens,TwoHam,nDT,
     &                 EOrb,OccNo,nEO,CMO,nCMO,note,iCase,MssVlc,Darwin)
************************************************************************
*                                                                      *
*     purpose: Final printout                                          *
*                                                                      *
*     input:                                                           *
*       Dens    : the last total density                               *
*       TwoHam  : the last two-electron ham.                           *
*       OneHam  : one-electron hamiltonian of length nOOK              *
*       Ovlp    : overlap in AO basis of length nOOK                   *
*       EOrb    : orbital energies of length nEO                       *
*       OccNo   : occupation numbers of length nEO                     *
*       CMO     : molecular orbital coefficients of length nCMO        *
*                                                                      *
*     called from: Final                                               *
*                                                                      *
*     calls to: R1IntB, RelEny, SetUp, Ortho                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: UHF - V.Veryazov, 2003                                  *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      External EFP_ON
*
      Real*8 Dens(nDT),TwoHam(nDT),OneHam(nDT),Ovlp(nDT),
     &       EOrb(nEO),OccNo(nEO),CMO(nCMO),MssVlc(nDT),Darwin(nDT)
      Character*80 Note
*
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "file.fh"
#include "stdalloc.fh"
#include "rctfld.fh"
#include "oneswi.fh"
      Logical Do_SpinAV
      COMMON  / SPAVE_L  / Do_SpinAV
*
*---- Define local variables
      Character Fmt*60
      Logical PrEne,PrOcc, get_BasisType
      Logical DeBug, First, NonEq, Dff, Do_DFT, FullMlk, Reduce_Prt
      External Reduce_Prt
      Real*8, Dimension(:), Allocatable:: RFfld, Scr2, Scr3
cnf
      Logical Do_ESPF, EFP_On
cnf
*     Save ipScr1
      Character AlphaLabel*30
      Dimension Dumm0(1),Dumm1(1)
*
*----------------------------------------------------------------------*
*     Warning!
*  this routine uses a dirty trick in UHF case
*  ipScr1 is used to keep a temporary array in memory
*  which is released during a second call.
*----------------------------------------------------------------------*
*


#include "SysDef.fh"

*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      jPrint=iPrint
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0
      If (iPL.le.1) jPrint=1
*
*----------------------------------------------------------------------*
*
#ifdef _DEBUGPRINT_
      Debug=.true.
#else
      Debug=.false.
      If (jPrint.ge.4) Debug=.True.
#endif
*
      Fmt = '(6X,A,T50,F17.10)'
      AlphaLabel=' '
      if(iUHF.eq.1.and.iCase.eq.0) AlphaLabel=' (alpha) '
      if(iUHF.eq.1.and.iCase.eq.1) AlphaLabel=' (beta)  '

      if (Do_SpinAV) AlphaLabel=Alphalabel(1:9)//'and (spin-averaged)'

*---- Compute relativistic corrections
      If (lRel) Then
         Call RelEny(ERelMV,ERelDC,Dens,
     &               MssVlc,Darwin,nBT)
         If (jPrint.ge.2) Then
            Write(6,*)
            Write(6,'(6X,A)') '1st order relativistic corrections'
            Write(6,Fmt)'Total energy',EneV+ERelMV+ERelDC
            Write(6,Fmt)'Mass-velocity correction',ERelMV
            Write(6,Fmt)'1-el Darwin correction',ERelDC
            Write(6,Fmt)'Sum of relatvity corrections',ERelDC+ERelMV
            Write(6,*)
         End If
      End If
*
*---- Print numerical quadrature information
      iSpin=1
      if(iUHF.eq.1) iSpin=2
      If (KSDFT.ne.'SCF'.and.iCase.eq.0) Call Print_NQ_Info(iSpin)
*
*---- Write out last density matrix to output
      If (DeBug) Then
         Write(6,'(6x,A)')'Last density matrix (interpolated) '
     &                  //'in AO basis'
         ij = 1
         Do iSym = 1, nSym
            Write(6,*)' symmetry',iSym
            Call TriPrt(' ',' ',Dens(ij),nBas(iSym))
            ij=ij+nBas(iSym)*(nBas(iSym)+1)/2
         End Do
         Write(6,*)
         Write(6,'(6x,A)')'Last 2-el. Hamiltonian (interpolated) '
     &                  //'in AO basis'
         ij = 1
         Do iSym = 1, nSym
            Write(6,*)' symmetry',iSym
            Call TriPrt(' ',' ',TwoHam(ij),nBas(iSym))
            ij=ij+nBas(iSym)*(nBas(iSym)+1)/2
         End Do
         Write(6,*)
         Write(6,'(6x,A)')'Last 1-el. Hamiltonian (interpolated) '
     &                  //'in AO basis'
         ij = 1
         Do iSym = 1, nSym
            Write(6,*)' symmetry',iSym
            Call TriPrt(' ',' ',OneHam(ij),nBas(iSym))
            ij=ij+nBas(iSym)*(nBas(iSym)+1)/2
         End Do
         Write(6,*)
      End If
*
*...  Process the external potential with the total electronic density
*
cnf
      Call DecideOnESPF(Do_ESPF)
      If ( (Do_ESPF .or. lRF .or. KSDFT.ne.'SCF' .or. EFP_On()) .and.
     &      .Not.NDDO .and. iCase.eq.0) Then
cnf      If ( (lRF .or. KSDFT.ne.'SCF') .and.  .Not.NDDO .and.
         iCharge=Int(Tot_Charge)
         NonEq=.False.
*        Call Get_PotNuc(PotNuc)
*        Call Get_dScalar('PotNuc',PotNuc)
         Call Peek_dScalar('PotNuc',PotNuc)
         Call mma_allocate(RFfld,nBT,Label='RFfld')
         call dcopy_(nBT,[Zero],0,RFfld,1)
         First=.True.
         Dff = .False.
         Do_DFT=.False. ! We do not need to redo the DFT!
         iDumm=1
         Call DrvXV(RFfld,RFfld,Dens,
     &             PotNuc,nBT,First,Dff,NonEq,lRF,KSDFT,ExFac,
     &             iCharge,iSpin,
     &             Dumm0,Dumm1,iDumm,'SCF ',Do_DFT)
         Call mma_deallocate(RFfld)
*
*------- Print multipole analysis of the reaction field contributions
*
         Call RFmltp()
*
      End If
*
*
*---- Print orbitals (the case InVec=3 and nIter=0 is set up in RdInp)
      Fmt = '(6X,A)'
      If (iPrOrb.ge.1) Then
         FullMlk=.True.
         PrEne=.True.
         PrOcc=.True.
         If (iPrOrb.eq.1) Then
            EHomo = -99999.0d+00
            ELumo =  99999.0d+00
            ij = 1
            Do iSym = 1, nSym
              Do iv=1,nOrb(iSym)
              if(OccNo(ij+iv-1).gt.0.001) then
               EHomo = Max(EHomo,EOrb(ij + iv - 1))
              else
               ELumo = Min(ELumo,EOrb(ij + iv -1 ))
              endif
              End Do
               ij    = ij + nOrb(iSym)
            End Do
cvv            ThrEne = Abs(Two*EHomo - ELumo)
            ThrEne = ELumo+0.5
            If (jPrint.ge.2) Then
               Write(6,*)
               Write(6,Fmt)'All orbitals with orbital energies smaller '
     &                //'than  E(LUMO)+0.5 are printed'
            End If
c         print *,'Elumo',Elumo
         Else
            If (jPrint.ge.2) Then
               Write(6,*)
               Write(6,'(6X,2A,E11.4,A)')'All orbitals with orbital',
     &              ' energies smaller than',ThrEne,' are printed'
            End If
         End If
         ThrOcc = -99999.0d+00
         If (KSDFT.eq.'SCF') Then
            If(iUHF.eq.0) Then
               Note='SCF orbitals'//AlphaLabel
               If (kIvo.ne.0)  Note='SCF orbitals + IVO'
               If (iCoCo.ne.0)
     &            Note='SCF orbitals + arbitrary occupations'
            Else
               Note='UHF orbitals'//AlphaLabel
               If (kIvo.ne.0)  Note='UHF orbitals + IVO'
               If (iCoCo.ne.0)
     &            Note='UHF orbitals + arbitrary occupations'
            End If
         Else
            If(iUHF.eq.0) Then
               Note='RKS-DFT orbitals'//AlphaLabel
               If (kIvo.ne.0)  Note='RKS-DFT orbitals + IVO'
               If (iCoCo.ne.0)
     &            Note='RKS-DFT orbitals + arbitrary occupations'
            Else
               Note='UKS-DFT orbitals'//AlphaLabel
               If (kIvo.ne.0)  Note='UKS-DFT orbitals + IVO'
               If (iCoCo.ne.0)
     &            Note='UKS-DFT orbitals + arbitrary occupations'
            End If
         End If
         If (jPrint.ge.2)
     &      Call PriMO(Note,PrOcc,PrEne,ThrOcc,ThrEne,
     &              nSym,nBas,nOrb,Name,EOrb,OccNo,CMO,iPrForm)
      Else
         If (jPrint.ge.2)
     &      Write(6,Fmt)'No orbitals printed'
         FullMlk=.False.
      End If
*
      If (InVec.ne.3 . or . nIter(nIterP).gt.0) Then
         Call mma_allocate(Scr2,nBB,Label='Scr2')
         Call mma_allocate(Scr3,nnB,Label='Scr3')
*
*------- Prepare CMO in symmetry blocks nBas x nBas
         iVec  = 0
         iCMO  = 0
         Do iSym = 1, nSym
            Do i = 1, nBas(iSym)*nOrb(iSym)
               Scr2(iVec + i) = CMO(iCMO + i)
            End Do
            iVec = iVec + nBas(iSym)*nOrb(iSym)
            Do i = 1, nBas(iSym)*(nBas(iSym) - nOrb(iSym))
               Scr2(iVec + i) = Zero
            End Do
            iVec = iVec + nBas(iSym)*(nBas(iSym) - nOrb(iSym))
            iCMO = iCMO + nOrb(iSym)*nBas(iSym)
         End Do
*
*------- Prepare occupation numbers
*
         iOr = 0
         iBs = 0
         Do iSym = 1, nSym
            Do j = 1, nOrb(iSym)
               Scr3(iBs + j) = OccNo(iOr + j)
            End Do
            Do j = nOrb(iSym) + 1, nBas(iSym)
               Scr3(iBs + j) = Zero
            End Do
            iOr = iOr + nOrb(iSym)
            iBs = iBs + nBas(iSym)
         End Do

         If(.not.NoProp) Then
*
*------- Population analysis
*
         jCase=iCase
         if(iUHF.eq.0) jCase=2
         Call Charge(nSym,nBas,Name,Scr2,Scr3,
     &               Ovlp,jCase,FullMlk,.True.)

         If (get_BasisType('ANO')) Then
            iRc = 0
            Call LoProp(iRc)
*
*------- NBO analysis
*
            Call Nat_Bond_Order(nSym,nBas,Name,jCase)
         End If
         End If
*
*------- ESPF analysis
*
         If (Do_ESPF) Call espf_analysis(.True.)
*
         Call mma_deallocate(Scr3)
         Call mma_deallocate(Scr2)
*
      End If

      End subroutine PrFin
