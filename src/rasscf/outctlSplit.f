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
      Subroutine OutCtlSplit(CMO,OCCN,SMAT,lOPTO)
************************************************************************
*                                                                      *
*     Control section for the final output                             *
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
************* GLMJ ************
      use OneDat, only: sNoOri, sOpSiz
      use rctfld_module, only: lRF
      use general_data, only: CleanMask
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_global, only: CBLBM, CMax, DE, ECAS, ESX, FDIAG,
     &                         HALFQ, IBLBM, iPCMRoot, iSPDen, iSupSM,
     &                         iSymBB, ITER, jBLBM, KSDFT, NAC, NACPAR,
     &                         NACPR2, BName, NIN, NONEQ, NSEC, OutFmt1,
     &                         RFPert, RlxGrd, RotMax, Tot_Charge,
     &                         Tot_El_Charge, Tot_Nuc_Charge, Via_DFT,
     &                         ixSym, iADR15, IPT2, iRLXRoot, Ener
      use SplitCas_Data, only: MxIterSplit,ThrSplit,
     &                         lRootSplit,EnerSplit,GapSpli,PerSplit,
     &                         PerCSpli,iDimBlockA
      use printlevel, only: DEBUG,USUAL,TERSE,VERBOSE
      use output_ras, only: LF,IPRLOC
      use general_data, only: NACTEL,NHOLE1,NELEC3,ISPIN,STSYM,NSYM,
     &                        NTOT1,NCONF,JOBIPH,NASH,NBAS,NDEL,NFRO,
     &                        NISH,NRS1,NRS2,NRS3,NSSH,NTOT,NTOT2
      use spinfo, only: NCSASM,NDTASM

      Implicit None

#include "rasdim.fh"
      Character(LEN=16), Parameter:: ROUTINE='OUTCTL  '
      Real*8, Allocatable:: DSave(:)
      Character(LEN=8)  Fmt2, Label
      Character(LEN=3) lIrrep(8)
      Character(LEN=80) Note
      Character(LEN=120) Line
      Logical FullMlk, get_BasisType
cnf
      Logical Do_ESPF,lSave, lOPTO
cnf
      Real*8 CMO(*),OCCN(*),SMAT(*)
      Real*8 Temp(2,mxRoot)
      Real*8 Dum(1)
      Integer iDum(56)

** (SVC) added for new supsym vector input
*     Integer NXSYM(mxOrb),nUND(mxOrb)

      Integer, External::  Cho_X_GetTol
      Real*8, Allocatable:: Tmp0(:), X1(:), X2(:), X3(:), X4(:),
     &                      CMON(:), D(:), X6(:), CMOSO(:)

      Real*8 CASDFT_Funct, EAV, Edc, Emv, EneTmp, Erel, percent
      Integer i, IAD03, IAD12, IAD14, IAD15, iCharge, iComp, iDimN,
     &        iDimO, iDImV, iEnd, Ind, iOpt, iPrLev, iRC, iRC1, iRC2,
     &        iStart, iSyLbl, iSym, iTemp, iTol, kRoot, left, LuTmp,
     &        NAO, nDCInt, nMVInt, NO
      Integer, External:: IsFreeUnit

*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
C Local print level (if any)
      IPRLEV=IPRLOC(6)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

* Additional DFT correlation energy, if any
      CASDFT_Funct=0.0d0
      If (KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM') Then
          Call Get_dScalar('CASDFT energy',CASDFT_Funct)
      End If
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      left=6
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'

      IF (IPRLEV.GE.USUAL .AND..NOT.lOPTO) THEN
* Start of long if-block A over IPRLEV
*----------------------------------------------------------------------*
*     Print orbital and wavefunction specifications                    *
*----------------------------------------------------------------------*
        Write(LF,*)
        Line=''
        Write(Line(left-2:),'(A)') 'Wave function specifications:'
        Call CollapseOutput(1,Line)
        Write(LF,Fmt2//'A)')'-----------------------------'
        Write(LF,*)
        Write(LF,Fmt2//'A,T45,I6)')'Number of closed shell electrons',
     &                           2*NIN
      Write(LF,Fmt2//'A,T45,I6)')'Number of electrons in active shells',
     &                           NACTEL
        Write(LF,Fmt2//'A,T45,I6)')'Max number of holes in RAS1 space',
     &                           NHOLE1
        Write(LF,Fmt2//'A,T45,I6)')'Max nr of electrons in RAS3 space',
     &                           NELEC3
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
        Line=''
        Write(Line(left-2:),'(A)') 'Orbital specifications:'
        Call CollapseOutput(1,Line)
        Write(LF,Fmt2//'A)')'-----------------------'
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
        Line=''
        Write(Line(left-2:),'(A)') 'CI expansion specifications:'
        Call CollapseOutput(1,Line)
        Write(LF,Fmt2//'A)')'----------------------------'
        Write(LF,*)
        Write(LF,Fmt2//'A,T40,I11)')'Number of CSFs',
     &                           NCSASM(STSYM)
        Write(LF,Fmt2//'A,T40,I11)')'Number of determinants',
     &                           NDTASM(STSYM)
        write(LF,Fmt2//'A,T45,I6)')  'Root required ', lrootSplit
        if (EnerSplit)
     &    write(LF,Fmt2//'A,T44,F7.2)')'Energy Gap (eV) in SplitCAS',
     &                            GapSpli
        if (PerSplit)
     &    write(LF,Fmt2//'A,T44,F7.1)')'Percentage sought in SplitCAS',
     &                            PercSpli
        percent = Real(iDimBlockA)/Real(NCSASM(STSYM))*100.0d0
        write(LF,Fmt2//'A,T42,I9,A,F5.1,A)')'A-Block Size in '//
     &                          'SplitCAS (CSFs)',iDimBlockA,' (',
     &                          percent,' %)'
        Write(LF,Fmt2//'A,T45,ES10.3)')'Threshold for SplitCAS',
     &                            ThrSplit
        write(LF,Fmt2//'A,T47,I4)') 'Maximum number of SplitCAS '//
     &                     'iterations', MxIterSplit
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
              Write(LF,*) 'OutCtl: iRc from Call RdOne not 0'
              Write(LF,*) 'Label = ',Label
              Write(LF,*) 'iRc = ',iRc
              Call Abend
           Endif
           Call mma_deallocate(Tmp0)
           Tot_El_Charge=0.0D0
           Do iSym=1,nSym
              Tot_El_Charge=Tot_El_Charge
     &                     -2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
           End Do
           Tot_El_Charge=Tot_El_Charge-DBLE(nActEl)
           Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
           iCharge=Int(Tot_Charge)
           Call PrRF(.False.,NonEq,iCharge,2)
           Write(LF,Fmt2//'A,T45,I2)')' Reaction field from state:',
     &                                IPCMROOT
        End If
        If ( RFpert ) then
           Write(LF,*)
           Write(LF,*)
           Write(LF,Fmt2//'A)')'Reaction field specifications:'
           Write(LF,Fmt2//'A)')'------------------------------'
           Write(LF,*)
           Write(LF,'(6X,A)')'The Reaction field has been added as a '//
     &                      'perturbation and has been determined '//
     &                      'in a previous calculation'
           Write(LF,*)
        End If
        If (KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM') Call Print_NQ_Info()
        Call CollapseOutput(0,'CI expansion specifications:')

* End of long if-block A over IPRLEV
      END IF

      EAV = ENER(lRootSplit,ITER)

      CASDFT_Funct=0.0d0
      If (KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM') Then
         If(nConf.eq.1.and.nActEl.eq.nac) Then
           Call Get_dScalar('CASDFT energy',CASDFT_Funct)
           EAV = EAV + CASDFT_Funct
         Else
           Call Get_dScalar('CASDFT energy',CASDFT_Funct)
           EAV=EAV+CASDFT_Funct - VIA_DFT - HALFQ
         End If
      End If

      IF(IPRLEV.GE.USUAL .AND..NOT.lOPTO) THEN
* Start of long if-block B over IPRLEV
        Write(LF,*)
        Line=''
        Write(Line(left-2:),'(A)') 'Final optimization specifications:'
        Call CollapseOutput(1,Line)
        Write(LF,Fmt2//'A)')'------------------------------'
        Write(LF,*)

        Write(LF,Fmt2//'A,T45,F20.8)') 'SplitCAS CI-energy',EAV
        Write(LF,Fmt2//'A,T45,F20.8)') 'SplitCAS RAS-energy',ECAS
        Write(LF,Fmt2//'A,T45,F20.8)') 'Super-CI energy',ESX
        Write(LF,Fmt2//'A,T45,F20.8)') 'RASSCF energy change',DE
        Write(LF,Fmt2//'A,T50,ES10.3)')'Max change in MO coeff.',CMAX
        Write(LF,Fmt2//'A,T50,ES10.3)')
     &          'Max non-diagonal density matrix element',ROTMAX
        Write(LF,Fmt2//'A,T50,ES10.3)') 'Max. BLB matrix element',CBLBM
        Write(LF,Fmt2//'A,I4,A,I4,A,I4,A)')
     &     '(orbital pair',IBLBM,',',JBLBM,' in symmetry',ISYMBB,')'
        If (irlxRoot.ne.0)
     &    Write(LF,Fmt2//'A,T45,ES10.3)')
     &         'Norm of electronic gradient',RLXGRD

        If ( ISUPSM.ne.0 ) then
          Write(LF,Fmt2//'A)')
     &    'Supersymmetry has been used to disable selected orbital'//
     &    ' rotations'
          iEnd=0
          Do iSym=1,nSym
            iStart=iEnd+1
            iEnd=iEnd+nBas(iSym)
            iTemp=0
            Do i=iStart,iEnd
              iTemp=iTemp+IXSYM(i)
            End Do
            If ( iTemp.gt.0 ) then
               Write(LF,Fmt2//'A,I3)')
     &         'Supersymmetry vector for symmetry species',iSym
               Write(LF,Fmt2//'30I3)') (IXSYM(i),i=iStart,iEnd)
            End If
          End Do
        End If

        If (Allocated(CleanMask)) then
          Write(LF,Fmt2//'A)')
     &    'The cleanup option has been used to set MO-coefficients'//
     &    ' explicitly to zero'
        End If
        Call CollapseOutput(0,'Final optimization conditions:')

* End of long if-block B over IPRLEV
      END IF

      call dcopy_(2*mxRoot,[0.0d0],0,Temp,1)
      iRc1=0
      iRc2=0
      iOpt=ibset(0,sOpSiz)
      iComp=1
      iSyLbl=1
      nMVInt=0
      nDCInt=0
      Label='MassVel'
      Call iRdOne(iRc1,iOpt,Label,iComp,iDum,iSyLbl)
      If (iRc1.eq.0) nMVInt=iDum(1)
      Label='Darwin'
      Call iRdOne(iRc2,iOpt,Label,iComp,iDum,iSyLbl)
      If (iRc2.eq.0) nDCInt=iDum(1)
      If ( (nMVInt+nDCInt).ne.0 ) Then
        IAD12=IADR15(12)
        CALL mma_allocate(X1,NTOT1,Label='X1')
        CALL mma_allocate(X2,NTOT1,Label='X2')
        CALL mma_allocate(X3,NTOT ,Label='X3')
        CALL mma_allocate(X4,NTOT2,Label='X4')

        kRoot = lRootSplit
        CALL DDAFILE(JOBIPH,2,X4,NTOT2,IAD12)
        CALL DDAFILE(JOBIPH,2,X3,NTOT,IAD12)
        Call RelEne(Temp(1,kRoot),Temp(2,kRoot),nSym,nBas,
     &              X4,X3,X2,X1)

        CALL mma_deallocate(X4)
        CALL mma_deallocate(X3)
        CALL mma_deallocate(X2)
        CALL mma_deallocate(X1)
      End If

      If ( (nMVInt+nDCInt).ne.0) then
        Emv=Temp(1,lRootSplit)
        Edc=Temp(2,lRootSplit)
        Erel=ENER(lRootSplit,ITER)+Emv+Edc
        EneTmp = Erel
      Else
        EneTmp = ENER(lRootSplit,ITER)+CASDFT_Funct-VIA_DFT-HALFQ
      End If

      IF(IPRLEV.GE.USUAL .AND..NOT.lOPTO) THEN
* Start of long if-block C over IPRLEV
*----------------------------------------------------------------------*
*     Print final state energies (including relativistic correction)   *
*----------------------------------------------------------------------*
        Write(LF,*)
        Write(LF,*)
        Write(LF,Fmt2//'A)')'Final state energy(ies):'
        Write(LF,Fmt2//'A)')'------------------------'
        Write(LF,*)
        If ( (nMVInt+nDCInt).ne.0) then
         Write(LF,Fmt2//'A)')
     &        'root     nonrelativistic        mass-velocity       '//
     &        'Darwin-contact         relativistic'
         Write(LF,Fmt2//'A)')
     &        '             energy                term             '//
     &        '     term                 energy'
          i=lRootSplit
          Emv=Temp(1,i)
          Edc=Temp(2,i)
          Erel=ENER(I,ITER)+Emv+Edc
          Write(LF,Fmt2//'I3,4(1X,F20.8))')
     &          I,ENER(I,ITER),Emv,Edc,Erel
        Else
          Write(LF,Fmt2//'A,I3,A,F20.8,A)')
     &          'root number',lrootSplit,' E =',
     &            EneTmp,' a.u.'
        End If

      ELSE IF(IPRLEV.GE.TERSE .AND..NOT.lOPTO) THEN
        Write(LF,Fmt2//'A,I3,A,F20.8,A)')
     &    'root number',lRootSplit,' E =',
     &    EneTmp,' a.u.'

* End of long if-block C over IPRLEV
      END IF

      call Store_Energies(1,[EneTmp],1)

      Call Put_iScalar('NumGradRoot',irlxroot)
      Call Put_dScalar('Last energy',ECAS)
      Call Put_dScalar('Average energy',EAV)

      iTol = Cho_X_GetTol(8)
      Call Add_Info('E_RASSCF',ENER(1,ITER),1,iTol)

*---------------------------------------------------------------
* New JOBIPH layout: Also write hamiltonian matrix at IADR15(17):
*---------------------------------------------------------------
      IAD15=IADR15(17)
      CALL DDAFILE(JOBIPH,1,ENER(1,ITER),1,IAD15)

*----------------------------------------------------------------------*
*     Print the multipole analysis of the solvation energy             *
*----------------------------------------------------------------------*
      Call RFMltp()

*----------------------------------------------------------------------*
*     Print the final orbitals                                         *
*----------------------------------------------------------------------*
      FullMlk=(OutFmt1.NE.'NOTHING ')
      IF (IPRLEV.GE.USUAL .AND..NOT.lOPTO) THEN
* Start of if-block D over IPRLEV
        If (OutFmt1.NE.'NOTHING ') then
          If(IPT2.EQ.0) THEN
            CALL PRIMO_RASSCF('Pseudonatural active'//
     &                  ' orbitals and approximate occupation numbers',
     &                  FDIAG,OCCN,CMO)
          Else
            CALL PRIMO_RASSCF('All orbitals are'//
     &                  ' eigenfunctions of the PT2 Fock matrix',
     &                  FDIAG,OCCN,CMO)
          EndIf
        EndIf
* End of if-block D over IPRLEV
      END IF

*----------------------------------------------------------------------*
*     Also put on RUNFILE (in the future...):
*----------------------------------------------------------------------*
      Call Put_dArray('RASSCF orbitals',CMO,NTOT2)

*----------------------------------------------------------------------*
*     compute properties and Mulliken's orbital populations            *
*----------------------------------------------------------------------*
      IAD12=IADR15(12)
      IAD03=IADR15(3)
      IAD14=IADR15(14)
*BOR0511
*     Save original orbitals for the spin density matrices
      Call mma_allocate(cmon,nTot2,Label='CMON')
      call dcopy_(ntot2,cmo,1,cmon,1)
*BOR0511
      FullMlk=(OutFmt1.NE.'NOTHING ')

************************************************************************
*     save the 1st order density for gradients                         *
************************************************************************
      Call mma_allocate(DSave,nTot1,Label='DSave')
      Call Get_dArray_chk('D1ao',DSave,NTOT1)

* Read natural orbitals
      If ( NAC.GT.0 ) then
        CALL DDAFILE(JOBIPH,2,CMO,NTOT2,IAD12)
        CALL DDAFILE(JOBIPH,2,OCCN,NTOT,IAD12)
      End If
*
* Put the density matrix of this state on the runfile for
*  LoProp utility
      Call mma_allocate(D,nTot1,Label='D')
      D(:)=0.0D0
      Call DONE_RASSCF(CMO,OCCN,D)
      Call Put_dArray('D1ao',D,NTOT1)
      Call mma_deallocate(D)

      IF (IPRLEV.GE.USUAL) THEN
*       Start of if-block E over IPRLEV
*
*       Compute Mulliken's population analysis
        Write(LF,'(/6X,A,I3)')
     *  'Mulliken population Analysis for root number:',lRootSplit
        Write(LF,'(6X,A)')
     *  '-----------------------------------------------'
        Write(LF,*)
        lSave = lRootSplit.eq.iRlxRoot
        CALL CHARGE(nsym,nbas,BName,CMO,OCCN,SMAT,2,FullMlk,lSave)
        Write(LF,*)
*
*       Compute properties
        Write(LF,'(/6X,2A,I3)')
     *  'Expectation values of various properties ',
     *  'for root number:',lRootSplit
        Write(LF,'(6X,2A)')
     *  '-----------------------------------------',
     *  '------------------'
        Write(LF,*)
      END IF
* End of if-block E over IPRLEV

*     Write out info on a temporary vector file. Prpt
*     will need to read this.
      Note='Temporary orbital file used by prpt3.'
      LuTmp=50
      LuTmp=IsFreeUnit(LuTmp)
      Call WrVec('TMPORB',LuTmp,'CO',nSym,nBas,nBas,
     &          CMO,OCCN,Dum,iDum,Note)
      CALL PRPT()
*
*     Compute spin orbitals and spin population
*     (Note: this section overwrites the pseudo natural orbitals
*            with the spin orbitals).
*
      CALL mma_allocate(X6,NACPAR,Label='X6')
      CALL DDAFILE(JOBIPH,0,X6,NACPAR,IAD03)
      CALL DDAFILE(JOBIPH,2,X6,NACPAR,IAD03)
      CALL DDAFILE(JOBIPH,0,X6,NACPR2,IAD03)
      CALL DDAFILE(JOBIPH,0,X6,NACPR2,IAD03)
      CALL DBLOCK(X6)

      IF (IPRLEV.GE.VERBOSE) THEN
* Start of long if-block F over IPRLEV
        IF(ISPDEN .EQ. 1) THEN
*         Print spin density matrix
          Write(LF,'(/6X,A,I3)')
     &    'Spin density matrix for root number:',lRootSplit
          Write(LF,'(6X,A)')
     &    '--------------------------------------'
          Write(LF,*)
          IND=1
          IDIMV=0
          IDIMO=0
          IDIMN=0
          DO ISYM=1,NSYM
            NAO=NASH(ISYM)
            IF(NAO .EQ. 0) GO TO 50
            NO=NBAS(ISYM)
            IDIMV=IDIMV+NAO*NAO
            IDIMO=IDIMO+NAO
            IDIMN=IDIMN+NO*NAO
            If(iprlev.ge.VERBOSE) then
              Write(LF,'(/6X,A,I2)')
     &        'symmetry species',ISYM
              Write(LF,*)
              CALL TRIPRT(' ',' ',X6(IND),NASH(ISYM))
            Endif
            IND=IND+NASH(ISYM)*(NASH(ISYM)+1)/2
  50        CONTINUE
          END DO
        END IF

* End of long if-block F over IPRLEV
      END IF

*     Compute spin orbitals and spin population
      CALL DCOPY_(NTOT,[0.0D0],0,OCCN,1)
*SVC-11-01-2007 store original cmon in cmoso, which gets changed
      CALL mma_allocate(CMOSO,NTOT2,Label='CMOSO')
      CALL DCOPY_(NTOT2,CMON,1,CMOSO,1)

      CALL SPINORB(X6,CMOSO,OCCN,lRootSplit)
      CALL mma_deallocate(X6)
      CALL DDAFILE(JOBIPH,1,CMOSO,NTOT2,IAD14)
      CALL DDAFILE(JOBIPH,1,OCCN,NTOT,IAD14)

      IF (IPRLEV.GE.USUAL) THEN
* Start of long if-block G over IPRLEV
        IF(ISPDEN .EQ. 1) THEN
          Write(LF,'(/6X,A,I3)')
     &   'Mulliken spin population Analysis for root number:',lRootSplit
          Write(LF,'(6X,A)')
     &    '---------------------------------------------------'
          Write(LF,*)
          CALL CHARGE(nsym,nbas,BName,cmoso,OCCN,SMAT,3,FullMlk,.False.)
          Write(LF,*)
        ENDIF

*BOR0511
*       Do LoProp Charge analysis
        If (get_BasisType('ANO')) Then
          Write(LF,'(/6X,A,I3)')
     *    'LoProp population Analysis for root number:',lRootSplit
          Write(LF,'(6X,A)')
     *    '-----------------------------------------------'
          Write(LF,*)
          Write(LF,*)
          Line=''
          Write(Line(left-2:),'(A)') 'LoProp Analysis:'
          Call CollapseOutput(1,Line)
          Write(LF,Fmt2//'A)') '----------------'
          Write(LF,*)
          Write(LF,*)
* Ugly trick: use IRC to pass some info to LoProp ...
          iRC = Min((Abs(lRootSplit-iRlxRoot)),1)
          Call LoProp(iRC)
          Write(LF,*)
          Write(LF,'(6X,A,I2)')
     &    'Natural Bond Order Analysis for root number:',lRootSplit
          Call Nat_Bond_Order(nSym,nBas,BName,2)
          Call CollapseOutput(0,'LoProp Analysis:')
          Write(6,*)
        End If

* End of long if-block G over IPRLEV
      END IF
      CALL mma_deallocate(CMOSO)

*----- ESPF analysis
      Call DecideOnESPF(Do_ESPF)
      If (Do_ESPF) Call espf_analysis(.False.)


*     Restore the correct 1st order density for gradient calculations.
*
      Call Put_dArray('D1ao',DSave,NTOT1)
      Call mma_deallocate(DSave)
*                                                                      *
************************************************************************
*                                                                      *
      CALL mma_deallocate(CMON)
*----------------------------------------------------------------------*
      End Subroutine OutCtlSplit
