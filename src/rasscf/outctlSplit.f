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
      Implicit Real*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='OUTCTL  ')
#include "ciinfo.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
#include "splitcas.fh"
#include "SysDef.fh"

      Character*8  Fmt2, Label
      Character*3 lIrrep(8)
      Character*80 Note
      Character*120 Line
      Logical FullMlk, get_BasisType
cnf
      Logical Do_ESPF,lSave, lOPTO
cnf
      DIMENSION CMO(*),OCCN(*),SMAT(*)
      Dimension Temp(2,mxRoot)
      Dimension Dum(1),iDum(56)

** (SVC) added for new supsym vector input
*      DIMENSION NXSYM(mxOrb),nUND(mxOrb)

      Integer  Cho_X_GetTol
      External Cho_X_GetTol

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
     &                           LSYM
        Call CollapseOutput(0,'Wave function specifications:')
*
        Call Get_cArray('Irreps',lIrrep,24)
        Do iSym = 1, nSym
           Call RightAd(lIrrep(iSym))
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
     &                           NCSASM(LSYM)
        Write(LF,Fmt2//'A,T40,I11)')'Number of determinants',
     &                           NDTASM(LSYM)
        write(LF,Fmt2//'A,T45,I6)')  'Root required ', lrootSplit
        if (EnerSplit)
     &    write(LF,Fmt2//'A,T44,F7.2)')'Energy Gap (eV) in SplitCAS',
     &                            GapSpli
        if (PerSplit)
     &    write(LF,Fmt2//'A,T44,F7.1)')'Percentage sought in SplitCAS',
     &                            PercSpli
        percent = Real(iDimBlockA)/Real(NCSASM(LSYM))*100.0d0
        write(LF,Fmt2//'A,T42,I9,A,F5.1,A)')'A-Block Size in '//
     &                          'SplitCAS (CSFs)',iDimBlockA,' (',
     &                          percent,' %)'
        Write(LF,Fmt2//'A,T45,E10.3)')'Threshold for SplitCAS',
     &                            ThrSplit
        write(LF,Fmt2//'A,T47,I4)') 'Maximum number of SplitCAS '//
     &                     'iterations', MxIterSplit
        If ( lRF ) then
           Call GetMem('Ovrlp','Allo','Real',iTmp0,nTot1+4)
           iRc=-1
           iOpt=2
           iComp=1
           iSyLbl=1
           Label='Mltpl  0'
           Call RdOne(iRc,iOpt,Label,iComp,Work(iTmp0),iSyLbl)
           Tot_Nuc_Charge=Work(iTmp0+nTot1+3)
           If ( iRc.ne.0 ) then
              Write(LF,*) 'OutCtl: iRc from Call RdOne not 0'
              Write(LF,*) 'Label = ',Label
              Write(LF,*) 'iRc = ',iRc
              Call QTrace
              Call Abend
           Endif
           Call GetMem('Ovrlp','Free','Real',iTmp0,nTot1+4)
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
        If (KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM') Call Print_NQ_Info(iSpin)
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
        Write(LF,Fmt2//'A,T50,E10.3)')'Max change in MO coeff.',CMAX
        Write(LF,Fmt2//'A,T50,E10.3)')
     &          'Max non-diagonal density matrix element',ROTMAX
        Write(LF,Fmt2//'A,T50,E10.3)') 'Max. BLB matrix element',CBLBM
        Write(LF,Fmt2//'A,I4,A,I4,A,I4,A)')
     &     '(orbital pair',IBLBM,',',JBLBM,' in symmetry',ISYMBB,')'
        If (irlxRoot.ne.0)
     &    Write(LF,Fmt2//'A,T45,E10.3)')
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

        If ( ICLEAN.ne.0 ) then
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
      iOpt=1
      iComp=1
      iSyLbl=1
      nMVInt=0
      nDCInt=0
      Call iRdOne(iRc1,iOpt,'MassVel ',iComp,iDum,iSyLbl)
      If (iRc1.eq.0) nMVInt=iDum(1)
      Call iRdOne(iRc2,iOpt,'Darwin  ',iComp,iDum,iSyLbl)
      If (iRc2.eq.0) nDCInt=iDum(1)
      If ( (nMVInt+nDCInt).ne.0 ) Then
        IAD12=IADR15(12)
        CALL GETMEM('OPER','ALLO','REAL',LX1,NTOT1)
        CALL GETMEM('DEN1','ALLO','REAL',LX2,NTOT1)
        CALL GETMEM('OCCN','ALLO','REAL',LX3,NTOT )
        CALL GETMEM('CMON','ALLO','REAL',LX4,NTOT2)

        kRoot = lRootSplit
        CALL DDAFILE(JOBIPH,2,work(LX4),NTOT2,IAD12)
        CALL DDAFILE(JOBIPH,2,work(LX3),NTOT,IAD12)
        Call RelEne(Temp(1,kRoot),Temp(2,kRoot),nSym,nBas,
     &              work(LX4),work(LX3),work(LX2),work(LX1))

        CALL GETMEM('CMON','FREE','REAL',LX4,NTOT2)
        CALL GETMEM('OCCN','FREE','REAL',LX3,NTOT )
        CALL GETMEM('DEN1','FREE','REAL',LX2,NTOT1)
        CALL GETMEM('OPER','FREE','REAL',LX1,NTOT1)
      End If

      If ( (nMVInt+nDCInt).ne.0) then
        Emv=Temp(1,lRootSplit)
        Edc=Temp(2,lRootSplit)
        Erel=ENER(lRootSplit,ITER)+Emv+Edc
        EneTmp = Erel
*          Work(ipEneTmp+i-1) = eRel
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
      Call GetMem('cmon','ALLO','REAL',icmon,nTot2)
      call dcopy_(ntot2,cmo,1,work(icmon),1)
*BOR0511
      FullMlk=(OutFmt1.NE.'NOTHING ')

************************************************************************
*     save the 1st order density for gradients                         *
************************************************************************
      Call Get_D1AO(ipDSave,NTOT1)

* Read natural orbitals
      If ( NAC.GT.0 ) then
        CALL DDAFILE(JOBIPH,2,CMO,NTOT2,IAD12)
        CALL DDAFILE(JOBIPH,2,OCCN,NTOT,IAD12)
      End If
*
* Put the density matrix of this state on the runfile for
*  LoProp utility
      Call GetMem('DState','ALLO','REAL',ipD,nTot1)
      call dcopy_(nTot1,[0.0D0],0,Work(ipD),1)
      Call DONE_RASSCF(CMO,OCCN,Work(ipD))
      Call Put_D1AO(Work(ipD),NTOT1)
      Call Free_Work(ipD)

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
        CALL CHARGE(nsym,nbas,name,CMO,OCCN,SMAT,2,FullMlk,lSave)
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
      CALL GETMEM('RHO1S','ALLO','REAL',LX6,NACPAR)
      CALL DDAFILE(JOBIPH,0,WORK(LX6),NACPAR,IAD03)
      CALL DDAFILE(JOBIPH,2,WORK(LX6),NACPAR,IAD03)
      CALL DDAFILE(JOBIPH,0,WORK(LX6),NACPR2,IAD03)
      CALL DDAFILE(JOBIPH,0,WORK(LX6),NACPR2,IAD03)
      CALL DBLOCK(WORK(LX6))

      IF (IPRLEV.GE.VERBOSE) THEN
* Start of long if-block F over IPRLEV
        IF(ISPDEN .EQ. 1) THEN
*         Print spin density matrix
          Write(LF,'(/6X,A,I3)')
     &    'Spin density matrix for root number:',lRootSplit
          Write(LF,'(6X,A)')
     &    '--------------------------------------'
          Write(LF,*)
          IND=0
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
              CALL TRIPRT(' ',' ',WORK(LX6+IND),NASH(ISYM))
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
      CALL GETMEM('CMOSO','ALLO','REAL',ICMOSO,NTOT2)
      CALL DCOPY_(NTOT2,WORK(ICMON),1,WORK(ICMOSO),1)

      CALL SPINORB(WORK(LX6),WORK(ICMOSO),OCCN,lRootSplit)
      CALL GETMEM('RHO1S','FREE','REAL',LX6,NTOTSP)
      CALL DDAFILE(JOBIPH,1,WORK(ICMOSO),NTOT2,IAD14)
      CALL DDAFILE(JOBIPH,1,OCCN,NTOT,IAD14)

      IF (IPRLEV.GE.USUAL) THEN
* Start of long if-block G over IPRLEV
        IF(ISPDEN .EQ. 1) THEN
          Write(LF,'(/6X,A,I3)')
     &   'Mulliken spin population Analysis for root number:',lRootSplit
          Write(LF,'(6X,A)')
     &    '---------------------------------------------------'
          Write(LF,*)
          CALL CHARGE(nsym,nbas,name,Work(icmoso),OCCN,SMAT,3,FullMlk,
     &               .False.)
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
          Call Nat_Bond_Order(nSym,nBas,Name,2)
          Call CollapseOutput(0,'LoProp Analysis:')
          Write(6,*)
        End If

* End of long if-block G over IPRLEV
      END IF
      CALL GETMEM('CMOSO','FREE','REAL',ICMOSO,NTOT2)

*----- ESPF analysis
      Call DecideOnESPF(Do_ESPF)
      If (Do_ESPF) Call espf_analysis(.False.)


*     Restore the correct 1st order density for gradient calculations.
*
      Call Put_D1AO(Work(ipDSave),NTOT1)
      Call GetMem('DSave','Free','REAL',ipDSave,nTot1)
*                                                                      *
************************************************************************
*                                                                      *
      CALL GETMEM('CMON','FREE','REAL',icmon,NTOT2)
*----------------------------------------------------------------------*
      Return
      End
