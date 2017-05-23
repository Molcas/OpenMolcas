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
*               1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine WrInp_SCF(SIntTh)
************************************************************************
*                                                                      *
*     purpose: Write input                                             *
*                                                                      *
*     called from: Scf                                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*     modified by M. Schuetz @teokem.lu.se, 1995                       *
*                                                                      *
*     input: SIntTh: Threshold for Integral prescreening               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
      Real*8 SIntTh
*

#include "mxdm.fh"
#include "infscf.fh"
#include "infso.fh"
#include "rctfld.fh"
#include "ldfscf.fh"
*
*---- Define local variables
      Character*60 Fmt, FmtR, FmtI
      Character*72 Line
      Character*3 lIrrep(8)
      Logical NonEq
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _DEBUG_
      Call QEnter('WrInp')
#endif
      If (jPrint.ge.2) Then
         Call CollapseOutput(1,'   Input section:')
         Write(6,'(3X,A)')     '   --------------'
         Write(6,*)
      End If
*
*---- Print out header of the integral file
      If (jPrint.ge.2) Then
         Write(6,'(6X,A)')  'Header of the integral files:'
         Write(Line,'(A72)') Header(1)
         Write(6,'(6X,A)') Line(:mylen(Line))
         Write(Line,'(A72)') Header(2)
         Write(6,'(6X,A)') Line(:mylen(Line))
         Write(6,*)
      End If
*
*---- Print out title (if any)
      If (nTit.gt.0) Then
         If (jPrint.ge.3) Then
             Call Banner(Title,nTit,lPaper-7)
             Write(6,*)
         Else If (jPrint.ge.2) Then
             Write (6,'(6X,A)') ' Title:'
             Do iTit = 1, nTit
                Write (6,'(8X,A)') Title(iTit)(:mylen(Title(iTit)))
             End Do
         End If
      End If
*
*---- Print the coordinates of the system
      If (jPrint.ge.2) Call PrCoor
*
*---- Print out Print level
*     Write(*,'(1X,A,I3)') 'Print level=',iPrint
*     Write(*,*)
*
*
      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym = 1, nSym
         Call RightAd(lIrrep(iSym))
      End Do
*
      If (jPrint.ge.2) Then
        Call CollapseOutput(0,'   Input section:')
        Write(6,*)
      End If
*---- Print out orbital informations
      Fmt='(6X,A,T35,8I4)'
      If (jPrint.ge.2) Then
         Call CollapseOutput(1,'   Orbital specifications:')
         Write(6,'(3X,A)')     '   -----------------------'
         Write(6,*)
         Write(6,Fmt)'Symmetry species',         (i,i=1,nSym)
         Write(6,'(6X,A,T35,8(1X,A))')
     &         '                ',         (lIrrep(i),i=1,nSym)
         Write(6,Fmt)'Frozen orbitals',          (nFro(i),i=1,nSym)
      End If
*
      If (Aufb) Then
       if(iUHF.eq.0) Then
        If (nAufb(1).eq.-1) Then
           Tot_El_Charge=Tot_Charge-Tot_Nuc_Charge
*--------- Check that Tot_El_Charge is a multiple of two!
           mtmp=Int(-Tot_El_Charge+0.1D0)
           If (Mod(mtmp,2).ne.0) Then
              Write (6,*)
     &          'WrInp: Error in definition of molecular charge!'
              Write (6,*)
     &          'Current implementation only allows double occupations.'
              Write (6,*) 'Tot_Charge    :',Tot_Charge
              Write (6,*) 'Tot_El_Charge :',Tot_El_Charge
              Write (6,*) 'Tot_Nuc_Charge:',Tot_Nuc_Charge
              Call Abend()
           End If
           nAufb(1)=mtmp/2
        endif
       else
        If (nAufb(1).eq.-1) Then
           Tot_El_Charge=Tot_Charge-Tot_Nuc_Charge
*--------- Check that Tot_El_Charge is a multiple of two!
           mtmp=Int(-Tot_El_Charge+0.1D0)
c if ZSPIN is not set - make difference alpha-beta = 0 or 1
            nAufb(2)=(mtmp-iAu_ab)/2
            nAufb(1)=Int(-Tot_El_Charge-nAufb(2))
        endif
c          Write(6,*) ' CHARGE + UHF is un'
c           Call Abend()
        End If
        if(iUHF.eq.0.and.jPrint.ge.2) then
        Write(6,Fmt)'Aufbau',                 nAufb(1)
        else if (jPrint.ge.3) Then
        Write(6,Fmt)'Aufbau alpha',                 nAufb(1)
        Write(6,Fmt)'Aufbau beta ',                 nAufb(2)
        endif
        If (Teee.and.jPrint.ge.2) Then
           Write (6,'(a,f6.3)') '      Start temperature =',RTemp
           Write (6,'(a,f6.3)') '      End temperature   =',TStop
           Write (6,'(a,f6.3)') '      Temperature Factor=',TemFac
        End If
      Else
        if(iUHF.eq.0.and.jPrint.ge.2) then
        Write(6,Fmt)'Occupied orbitals',    (nOcc(i,1),i=1,nSym)
        Write(6,Fmt)'Secondary orbitals',   (nOrb(i)-nOcc(i,1),i=1,nSym)
        else if (jPrint.ge.2) Then
        Write(6,Fmt)'Occupied orbitals alpha',      (nOcc(i,1),i=1,nSym)
        Write(6,Fmt)'Occupied orbitals beta ',      (nOcc(i,2),i=1,nSym)
        Write(6,Fmt)'Secondary orbitals alpha',
     &              (nOrb(i)-nOcc(i,1),i=1,nSym)
        Write(6,Fmt)'Secondary orbitals beta',
     &              (nOrb(i)-nOcc(i,2),i=1,nSym)

        endif
      End If
*
      Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
      Call put_dscalar('Total Charge    ',Tot_Charge)
      If (jPrint.ge.2) Then
         Write(6,Fmt)'Deleted orbitals         ',(nDel(i),i=1,nSym)
         Write(6,Fmt)'Total number of orbitals ',(nOrb(i),i=1,nSym)
         Write(6,Fmt)'Number of basis functions',(nBas(i),i=1,nSym)
         Call CollapseOutput(0,'   Orbital specifications:')
         Write(6,*)
         Write(6,'(6X,A,T45,F10.3)') 'Molecular charge',Tot_Charge
         Write(6,*)
      End If
*
*---- Print out reaction field specifications
*
      iCharge=Int(Tot_Charge)
      NonEq=.False.
      Call PrRF(.False.,NonEq,iCharge,jPrint)
      If ( RFpert.and.jPrint.ge.2 ) then
         Write(6,'(6X,A)')'Reaction field specifications:'
         Write(6,'(6X,A)')'------------------------------'
         Write(6,*)
         Write(6,'(6X,A)')'The Reaction field is added as a '//
     &                    'perturbation and has been determined '//
     &                    'in a previos calculation'
         Write(6,*)
      End If
*
*---- Print out grid information in case of DFT
*
      If (KSDFT.ne.'SCF') Then
         Call Put_dScalar('EThr',EThr)
         Call Funi_Print
         If (jPrint.ge.2) Then
            Write(6,*)
            If (One_Grid) Then
               Write (6,'(6X,A)') 'The same grid will be used for all'
     &                          //' iterations.'
            Else
               Write (6,'(6X,A)') 'A smaller intermediate grid will b'
     &                          //'e used the first few iterations.'
            End If
            Write(6,*)
         End If
      End If
*
*---- Print out informations concerning Direct/Conventional scheeme
      FmtI= '(6X,A,T50,I8)'
      FmtR= '(6X,A,T50,E8.2)'
      If (jPrint.ge.2) Then
         Call CollapseOutput(1,'   Optimization specifications:')
         Write(6,'(3X,A)')     '   ----------------------------'
         Write(6,*)
      End If
      If (DSCF.and.jPrint.ge.2) Then
         If (nDisc.eq.0.and.nCore.eq.0) Then
            Write(6,'(6X,A)')'SCF Algorithm: Direct'
         Else If (nDisc*1024.ge.nCore) Then
            Write(6,'(6X,A)')'SCF Algorithm: Semi-direct'
*
*---------- The threshold to be used in the direct SCF procedure is
*           defined by the energy threshold.
            Write (6,FmtI)
     &        'Max MByte of integrals on disk/process:',nDisc
            Write(6,FmtR) 'Threshold for saving integrals on disc',
     &                    Thize
         Else
            Write(6,'(6X,A)')'SCF Algorithm: Semi-direct in-core'
*
*---------- The threshold to be used in the direct SCF procedure is
*           defined by the energy threshold.
            Write (6,FmtI)
     &        'Max kByte of integrals in memory/process:',nCore
            Write(6,FmtR) 'Threshold for saving integrals in memory',
     &                    Thize
         End If
         If (PreSch) Then
            Write(6,'(6X,A)')'Prescreening Scheme: '//
     &                       'Only Integral value'
         Else
            Write(6,'(6X,A)')'Prescreening Scheme: '//
     &                       'Integral*Density value'
         End If
      Else If (jPrint.ge.2) Then
         if(iUHF.eq.0) then
            if(.not.DoCholesky)then
              Write(6,'(6X,A)')'SCF Algorithm: Conventional'
            else
              If (DoLDF) Then
                 If (LDF_IntegralMode.eq.0) Then
                    Write(6,'(6X,A)')
     &          'SCF Algorithm: Direct LDF using conventional integrals'
                 Else If (LDF_IntegralMode.eq.1) Then
                    Write(6,'(6X,A)')
     &              'SCF Algorithm: Direct robust LDF'
                 Else If (LDF_IntegralMode.eq.2) Then
                    Write(6,'(6X,A)')
     &              'SCF Algorithm: Direct nonrobust LDF'
                 Else If (LDF_IntegralMode.eq.3) Then
                    Write(6,'(6X,A)')
     &              'SCF Algorithm: Direct half-and-half LDF'
                 Else
                    Write(6,'(6X,A,I6)') 'Unknown LDF integral mode:',
     &                                   LDF_IntegralMode
                    Call LDF_NotImplemented()
                 End If
                 If (LDF_IntegralPrescreening.lt.0.0d0) Then
                    Write(6,'(6X,A,A)')
     &              'Integral prescreening threshold determined from',
     &              ' target accuracy'
                 Else
                    Write(6,FmtR)
     &              'Integral prescreening threshold',
     &              LDF_IntegralPrescreening
                 End If
                 If (LDF_ContributionPrescreening.lt.0.0d0) Then
                    Write(6,'(6X,A,A)')
     &              'Contribution prescreening threshold determined',
     &              ' from target accuracy'
                 Else
                    Write(6,FmtR)
     &              'Contribution prescreening threshold',
     &              LDF_ContributionPrescreening
                 End If
              Else
                 Call Get_iScalar('System BitSwitch',iDoRI)
                 if (Iand(iDoRI,1024).Eq.1024) then
                    if (LKon) then
                       Write(6,'(6X,A)')'SCF Algorithm: LK-RI/DF'
                    else
                       Write(6,'(6X,A)')'SCF Algorithm: RI/DF'
                    endif
                 else
                    if (LKon) then
                       Write(6,'(6X,A)')'SCF Algorithm: LK-Cholesky'
                    else
                       Write(6,'(6X,A)')'SCF Algorithm: Cholesky'
                    endif
                 endif
              End If
            endif
         else
            if(.not.DoCholesky)then
              Write(6,'(6X,A)')'SCF Algorithm: Conventional USCF'
            else
              If (DoLDF) Then
                 If (LDF_IntegralMode.eq.0) Then
                    Write(6,'(6X,A)')
     &          'SCF Algorithm: Direct LDF using conventional integrals'
                 Else If (LDF_IntegralMode.eq.1) Then
                    Write(6,'(6X,A)')
     &              'SCF Algorithm: Direct robust LDF USCF'
                 Else If (LDF_IntegralMode.eq.2) Then
                    Write(6,'(6X,A)')
     &              'SCF Algorithm: Direct nonrobust LDF USCF'
                 Else If (LDF_IntegralMode.eq.3) Then
                    Write(6,'(6X,A)')
     &              'SCF Algorithm: Direct half-and-half LDF USCF'
                 Else
                    Write(6,'(6X,A,I6)') 'Unknown LDF integral mode:',
     &                                   LDF_IntegralMode
                    Call LDF_NotImplemented()
                 End If
                 If (LDF_IntegralPrescreening.lt.0.0d0) Then
                    Write(6,'(6X,A,A)')
     &              'Integral prescreening threshold determined from',
     &              ' target accuracy'
                 Else
                    Write(6,FmtR)
     &              'Integral prescreening threshold',
     &              LDF_IntegralPrescreening
                 End If
                 If (LDF_ContributionPrescreening.lt.0.0d0) Then
                    Write(6,'(6X,A,A)')
     &              'Contribution prescreening threshold determined',
     &              ' from target accuracy'
                 Else
                    Write(6,FmtR)
     &              'Contribution prescreening threshold',
     &              LDF_ContributionPrescreening
                 End If
              Else
                 Call Get_iScalar('System BitSwitch',iDoRI)
                 if (Iand(iDoRI,1024).Eq.1024) then
                    if (LKon) then
                       Write(6,'(6X,A)')'SCF Algorithm: LK-RI/DF USCF'
                    else
                       Write(6,'(6X,A)')'SCF Algorithm: RI/DF USCF'
                    endif
                 else
                    if (LKon) then
                       Write(6,'(6X,A)')
     &                               'SCF Algorithm: LK-Cholesky USCF'
                    else
                       Write(6,'(6X,A)')'SCF Algorithm: Cholesky USCF'
                    endif
                 endif
              End If
            endif
         endif
      End If
      If (NoExchange) Then
         Write(6,'(6X,A)')
     &   'NOTE: exchange contributions will not be computed'
      End If
*
      If (jPrint.ge.2) Then
*
*---- Print out informations concerning difference scheme used
      If (MiniDn) Then
         Write(6,'(6X,A)')'Minimized density differences are used'
      Else
        if(.not.DDnOFF)then
          Write(6,'(6X,A)')'D(i)-D(i-1) density differences are used'
        else
          Write(6,'(6X,A)')'The actual AO density is used'
        endif
      End If
      Write(6,FmtI)'Number of density matrices in core',nMem
*
*---- Print out number of iterations
      Write(6,FmtI) 'Maximum number of NDDO SCF iterations',nIter(0)
      Write(6,FmtI) 'Maximum number of HF  SCF iterations',nIter(1)
*
*---- Print out thresholds for SCF
      Write(6,FmtR) 'Threshold for SCF energy change',EThr
      Write(6,FmtR) 'Threshold for density matrix',DThr
      Write(6,FmtR) 'Threshold for Fock matrix',FThr
      Write(6,FmtR) 'Threshold for linear dependence',DelThr
      If (Diis) Then
         Write(6,FmtR) 'Threshold at which DIIS is turned on',
     &                 DiisTh
         Write(6,FmtR) 'Threshold at which QNR/C2DIIS is turned on',
     &                 QNRTh
         Write(6,FmtR) 'Threshold for Norm(delta) (QNR/C2DIIS)',
     &                 DltNTh
      End If
      If (DSCF) Then
         Write(6,FmtR) 'Threshold for contribution to Fock matrix',
     &                 SIntTh
      End If
*
*---- Print out acceleration scheeme used in SCF
      Fmt = '(6x,A,A)'
      If (Diis) Then
         Write(6,Fmt) 'DIIS extrapolation of the SCF procedure'
      Else If (iDKeep.eq.0) Then
         Write(6,Fmt) 'No damping of the SCF procedure'
      Else If (iDKeep.ge.3) Then
         Write(6,Fmt) 'Extended dynamic damping of the SCF procedure'
      Else
         Write(6,Fmt) 'Dynamic damping of the SCF procedure'
      End If
*
*---- Print out IVO information (if any)
      If (kIvo.ne.0) Write(6,Fmt) 'Improved virtual orbitals.'
*
*---- Print out information about orbitals punched on the file
      If (jVOut.le.0) Then
         Write(6,Fmt) 'No vectors punched'
      Else If (jVOut.eq.1) Then
         If(iUHF.eq.0) Then
            Write(6,Fmt) 'All non deleted orbitals punched on: SCFORB'
         Else
            Write(6,Fmt) 'All non deleted orbitals punched on: UHFORB'
         End If
      Else
         If(iUHF.eq.0) Then
            Write(6,Fmt) 'All orbitals punched on: SCFORB'
         Else
            Write(6,Fmt) 'All orbitals punched on: UHFORB'
         End If
      End If
      Call CollapseOutput(0,'   Optimization specifications:')
      Write(6,*)
*
*---- Print out
      If (InVec.eq.0) Then
         Write(6,Fmt) 'Starting vectors from core diagonalization'
      Else If (InVec.eq.1) Then
         Write(6,Fmt)
     &         'NDDO MOs are generated before actual HF SCF computation'
      Else If (InVec.eq.2) Then
         Write(6,Fmt) 'Input vectors read from INPORB'
         Write(6,Fmt) 'Orbital file label: ',VTitle(:mylen(VTitle))
      Else If (InVec.eq.3) Then
         Write(6,Fmt) 'Input density matrix read from RUNFILE'
      Else If (InVec.eq.4) Then
         Write(6,Fmt) 'Restart...'
      Else If (InVec.eq.5) Then
         Write(6,Fmt) 'Input vectors from NDDO calculation'
      Else
         Write(6,Fmt) StVec
      End If
      If (Scrmbl) Write(6,Fmt) 'Start orbitals are scrambled in order '
     &      //'to introduce symmetry breaking'
      Write(6,*)
      Write(6,*)
*
      End If
*
*     LDF: print a warning message if nonrobust or half-n-half is used
*
      If (DoLDF) Then
         If (LDF_IntegralMode.eq.2 .or. LDF_IntegralMode.eq.3) Then
            Call WarningMessage(0,
     &      'WARNING: nonrobust and half-n-half LDF integral '//
     &      'representations may lead to large errors and '//
     &      'provide inferior numerical stability compared to '//
     &      'the robust representation')
         End If
      End If
*
#ifdef _DEBUG_
      Call QExit('WrInp')
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Call XFlush(6)
      Return
      End
