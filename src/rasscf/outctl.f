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
      Subroutine OutCtl(CMO,OCCN,SMAT,lOPTO)
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
#ifdef _DMRG_
!     module dependencies
      use qcmaquis_interface_cfg
      use qcmaquis_interface_utility_routines, only: print_dmrg_info
#endif
      use OneDat, only: sNoOri, sOpSiz
      use rctfld_module, only: lRF
      use general_data, only: CleanMask
      use stdalloc, only: mma_allocate, mma_deallocate
      use gas_data, only: iDoGAS, NGAS, NGSSH
      use input_ras, only: KeyCION
      use rasscf_global, only: CBLBM, CMAX, DE, DoDMRG, ECAS, ESX,
     &                         FDIAG, HalfQ, IBLBM, ICICH, iPCMRoot,
     &                         iPT2, iRLXRoot, iSPDen, iSupSM, iSymBB,
     &                         ITER, JBLBM, kIVO, KSDFT, lRoots,
     &                         MaxOrbOut, NAC, NACPAR, NACPR2, BName,
     &                         NIN, NONEQ, nRoots, NSEC, OutFmt1,
     &                         RFPert, RLXGrd, RotMax,       Tot_Charge,
     &                         Tot_El_Charge, Tot_Nuc_Charge, via_DFT,
     &                         iRoot, Weight, iCI, cCI, ixSym, iADR15,
     &                         Ener
#if defined (_ENABLE_CHEMPS2_DMRG_)
      use rasscf_global, only: ThrE
#endif

#ifdef _ENABLE_DICE_SHCI_
      use rasscf_global, only: dice_eps1, dice_eps2, dice_iter,
     &                         dice_restart, dice_SampleN, dice_stoc,
     &                         nRef_dice, diceOcc
#endif
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
      use rasscf_global, only: DoBlockDMRG, MxDMRG
#endif
#ifdef _ENABLE_CHEMPS2_DMRG_
      use rasscf_global, only: ChemPS2_blb, ChemPS2_lrestart,
     &                         ChemPS2_Noise, ChemPS2_restart,
     &                         Davidson_Tol, Do3RDM, HFOcc,
     &                         Max_canonical, Max_Sweep
#endif
      use PrintLevel, only: DEBUG,USUAL,TERSE,VERBOSE
      use output_ras, only: LF,IPRLOC
      use general_data, only: NACTEL,NHOLE1,NELEC3,ISPIN,STSYM,NSYM,
     &                        NTOT1,NCONF,NTOT,JOBIPH,NASH,NBAS,NDEL,
     &                        NFRO,NISH,NRS1,NRS2,NRS3,NSSH,NSSH,NTOT2
      use spinfo, only: NCSASM,NDTASM
      use DWSol, only: DWSolv, DWSol_fixed, W_SOLV
      use Molcas, only: MxRoot
      use RASDim, only: MxRef

      Implicit None

      Real*8 CMO(*),OCCN(*),SMAT(*)
      Logical lOPTO
      Character(LEN=16), Parameter :: ROUTINE='OUTCTL  '

      Character(LEN=8)  Fmt2, Label
      Character(LEN=3) lIrrep(8)
      Character(LEN=80) Note
      Character(LEN=120) Line
#ifdef _ENABLE_CHEMPS2_DMRG_
      Character(LEN=3) SNAC
      Integer iHFOcc
#endif
      Logical FullMlk, get_BasisType
      Logical Do_ESPF,lSave, Do_DM
      Real*8 Temp(2,mxRoot)
      Real*8, Allocatable:: DSave(:), Tmp0(:), X1(:), X2(:), X3(:),
     &                      X4(:), HEFF(:,:), CMON(:), DM(:), DMs(:,:),
     &                      DState(:), X6(:), CMOSO(:), EneTmp(:)

** (SVC) added for new supsym vector input
*      DIMENSION NXSYM(mxOrb),nUND(mxOrb)

      Integer, External :: Cho_X_GetTol
      Real*8 Dum(1)
      Integer iDum(56)

      REAL*8 CASDFT_Funct, EAV, EDC, Emv, Erel, vNentropy, xnu
      Integer i, iAd03, iAd12, iAd14, iAd15, iCharge, iComp, iDimN,
     &        iDimO, iDimV, iEnd, iGAS, Ind, iOpt, iPrLev, iRC, iRC1,
     &        iRC2, iRef, iStart, iSyLbl, iSym, iTemp, iTol, j, kRoot,
     &        left, luTmp, NAO, nDCInt, nMVInt, NO
#ifdef _ENABLE_DICE_SHCI_
      Integer iref_dice
#endif
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

      IF (IPRLEV.GE.DEBUG .AND..NOT.lOPTO) THEN
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
      If(.not.iDoGas) then
        Write(LF,Fmt2//'A,T45,I6)')'Max number of holes in RAS1 space',
     &                           NHOLE1
        Write(LF,Fmt2//'A,T45,I6)')'Max nr of electrons in RAS3 space',
     &                           NELEC3
      End If
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
      IF(.not.iDoGas)then
        Write(LF,Fmt2//'A,T47,8I4)') 'RAS1 orbitals',
     &                            (nRs1(iSym),iSym=1,nSym)
        Write(LF,Fmt2//'A,T47,8I4)') 'RAS2 orbitals',
     &                            (nRs2(iSym),iSym=1,nSym)
        Write(LF,Fmt2//'A,T47,8I4)') 'RAS3 orbitals',
     &                            (nRs3(iSym),iSym=1,nSym)
      Else
        DO IGAS=1,NGAS
          Write(LF,Fmt2//'A,I1,A,T47,8I4)') 'GAS',IGAS,' orbitals',
     &                            (ngssh(igas,iSym),iSym=1,nSym)
        END DO
      End If

      Write(LF,Fmt2//'A,T47,8I4)') 'Secondary orbitals',
     &                            (nSsh(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Deleted orbitals',
     &                            (nDel(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Number of basis functions',
     &                            (nBas(iSym),iSym=1,nSym)
      If (kIVO) Then
        Write(LF,Fmt2//'A,T47)') 'Improved Virtual Orbitals '//
     &                           'option is used'
        Write(LF,Fmt2//'A,T47)') 'Molecular Orbitals are NOT '//
     &                           'suitable for CASPT2 & MRCI!'
      End If
      Call CollapseOutput(0,'Orbital specifications:')
      Write(LF,*)

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
      If(.Not.DoBlockDMRG) GoTo 113

#ifdef _ENABLE_DICE_SHCI_
      Line=' '
      Write(Line(left-2:),'(A)') 'DICE specifications:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt2//'A)')'--------------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A,T70,L6)')'Heat-bath configuration interaction
     &(JCTC, 2017, 13, 1595)', DoBlockDMRG
      Write(LF,Fmt2//'A,T45,L6)')'Semistochastic algorithm',Dice_stoc
      Write(LF,Fmt2//'A,T45,L6)')'Full restart',dice_restart
      Write(LF,Fmt2//'A,T45,I6)')'Max iterations',dice_iter
      Write(LF,Fmt2//'A,T45,ES10.3)')'Epsilon1',
     &                           dice_eps1
      Write(LF,Fmt2//'A,T45,ES10.3)')'Epsilon2',
     &                           dice_eps2
      Write(LF,Fmt2//'A,T45,I6)')'SampleN',
     &                           dice_sampleN
      Write(LF,Fmt2//'A,T45)')'Occupation guess'
      do iref_dice=1,nref_dice
         write(LF,Fmt2//'A)') trim(diceocc(iref_dice))
      enddo
      Call CollapseOutput(0,'DICE specifications:')

*     Skip printing CI specifications in DICE
      GoTo 114
#endif

      Line=''
      Write(Line(left-2:),'(A)') 'DMRG sweep specifications:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt2//'A)')'--------------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A,T45,I6)')'Number of renormalized basis',
     &                           MxDMRG
      Write(LF,Fmt2//'A,T45,I6)')'Number of root(s) required',
     &                           NROOTS
#ifdef _ENABLE_CHEMPS2_DMRG_
      Write(LF,Fmt2//'A,T45,I6)')'Maximum number of sweeps',
     &                           max_sweep
      Write(LF,Fmt2//'A,T45,I6)')'Maximum number of sweeps in RDM',
     &                           max_canonical
      Write(LF,Fmt2//'A,T45,ES10.3)')'Threshold for restarting',
     &                           chemps2_blb
      Write(LF,Fmt2//'A,T45,ES10.3)')'Minimum Davidson tolerance',
     &                           davidson_tol
      Write(LF,Fmt2//'A,T45,ES10.3)')'DMRG convergence threshold',
     &                           THRE/2.0
      Write(LF,Fmt2//'A,T45,ES10.3)')'Noise prefactor',
     &                           chemps2_noise
      Write(LF,Fmt2//'A,T45,L6)')'Restart from previous calculation',
     &                           chemps2_restart
      Write(LF,Fmt2//'A,T45,L6)')'Calculate 3-RDM and F.4-RDM',
     &                           Do3RDM
      Write(LF,Fmt2//'A,T45,I6)')'Restart scheme in 3-RDM and F.4-RDM',
     &                           chemps2_lrestart
      write(SNAC, '(I3)') NAC
      Write(LF,Fmt2//'A,T45,'//trim(adjustl(SNAC))//'I2)')
     &                           'Occupation guess',
     &                           (HFOCC(ihfocc), ihfocc=1,NAC)
#endif

* NN.14 FIXME: haven't yet checked whether geometry opt. works correctly with DMRG
      Write(LF,Fmt2//'A,T45,I6)')'Root chosen for geometry opt.',
     &                           IRLXROOT
      Call CollapseOutput(0,'DMRG sweep specifications:')

*     Skip printing CI specifications in DMRG-CASSCF
      GoTo 114

 113  Continue
#endif

      if (.not. doDMRG) then
        Line=''
        Write(Line(left-2:),'(A)') 'CI expansion specifications:'
        Call CollapseOutput(1,Line)
        Write(LF,Fmt2//'A)')'----------------------------'
        Write(LF,*)
        Write(LF,Fmt2//'A,T40,I11)')'Number of CSFs',
     &                           NCSASM(STSYM)
        Write(LF,Fmt2//'A,T40,I11)')'Number of determinants',
     &                           NDTASM(STSYM)
      end if
      Write(LF,Fmt2//'A,T45,I6)')'Number of root(s) required',
     &                           NROOTS
      If ( ICICH.eq.0 ) then
         If ( nRoots.eq.1 ) then
           if(doDMRG)then
            Write(LF,Fmt2//'A,(T45,10I6))')'DMRG root used',
     &                                  IROOT(1)
           else
            Write(LF,Fmt2//'A,(T45,10I6))')'CI root used',
     &                                  IROOT(1)
           end if
         Else
           if(doDMRG)then
            Write(LF,Fmt2//'A,(T45,10I6))')'DMRG roots used',
     &                                    (IROOT(i),i=1,nRoots)
           else
            Write(LF,Fmt2//'A,(T45,10I6))')'CI roots used',
     &                                    (IROOT(i),i=1,nRoots)
           end if
            Write(LF,Fmt2//'A,(T45,10F6.3))')'weights',
     &                                      (Weight(i),i=1,nRoots)
         End If
      Else
         Do i=1,nRoots
            Write(LF,Fmt2//'A,T45,I6)')'root',i
            Write(LF,Fmt2//'A,T45,10I6)')'Reference configuartions',
     &                                  (iCI(iRef,i),iRef=1,mxRef)
            Write(LF,Fmt2//'A,T45,10F6.3)')'CI-coeff',
     &                                    (cCI(iRef,i),iRef=1,mxRef)
          End Do
      End If
      if(doDMRG)then
        Write(LF,Fmt2//'A,T45,I6)')'highest root included in the DMRG',
     &                           LROOTS
      else
        Write(LF,Fmt2//'A,T45,I6)')'highest root included in the CI',
     &                           LROOTS
      end if
      If (irlxroot.ne.0)
     &Write(LF,Fmt2//'A,T45,I6)')'Root passed to geometry opt.',
     &                           iRlxRoot
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
     &                   -2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
         End Do
         Tot_El_Charge=Tot_El_Charge-DBLE(nActEl)
         Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
         iCharge=Int(Tot_Charge)
         Call PrRF(.False.,NonEq,iCharge,2)
         if (DWSolv%DWZeta == -12345d+00) then
             Write(LF,Fmt2//'A)')
     &         'Weights of the reaction field are specified by RFROOT'
             Write(LF,Fmt2//'(T45,10F6.3))') (W_SOLV(i),i=1,nRoots)
         else if (DWSolv%DWZeta < 0.0d+00) then
           Call DWSol_fixed(i,j)
           if (i==0 .and. j==0) then
             Write(LF,Fmt2//'A)') 'Unrecognized negative DWZeta (DWSOl)'
             Write(LF,Fmt2//'A,T51,A)')
     &         'Dynamically weighted solvation is ',
     &         'automatically turned off!'
           else
             Write(LF,Fmt2//'A,T45,I2,X,I2)')
     &         'Reaction field from states:', i, j
             if (max(i,j) > nRoots) then
               Write(LF,Fmt2//'A)')
     &           'The specified state is too high! Cannot proceed...'
               Call Quit_OnUserError()
             end if
           end if
         else if (IPCMROOT <= 0) then
           Write(LF,Fmt2//'A,T45,T15)')' Reaction field from state:',
     &                                 ' State-Averaged'
           if (DWSolv%DWZeta /= 0.0d+00) then
             Write(LF,Fmt2//'A,T51,A)')
     &         'Dynamically weighted solvation is ',
     &         'automatically turned off!'
             DWSolv%DWZeta = 0.0d+00
           end if
         else
           Write(LF,Fmt2//'A,T45,I2)')' Reaction field from state:',
     &                                IPCMROOT
           if (DWSolv%DWZeta > 0.0d+00) then
             Write(LF,Fmt2//'A,ES10.3,A,I1,A)')
     &         'Dynamically weighted solvation is used with DWSOlv = ',
     &         DWSolv%DWZeta," (DWTYpe = ",DWSolv%DWType,")"
             Write(LF,Fmt2//'A,(T45,10F6.3))')
     &         'Final weights for the reaction field',
     &                                   (W_SOLV(i),i=1,nRoots)
           end if
         end if
      End If
      If ( RFpert ) then
         Write(LF,*)
         Write(LF,*)
         Write(LF,Fmt2//'A)')'Reaction field specifications:'
         Write(LF,Fmt2//'A)')'------------------------------'
         Write(LF,*)
         Write(LF,'(6X,A)')'The Reaction field has been added as a '//
     &                    'perturbation and has been determined '//
     &                    'in a previous calculation'
         Write(LF,*)
      End If
      If (KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM') Call Print_NQ_Info()
      Call CollapseOutput(0,'CI expansion specifications:')

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
 114  Continue
#endif

* End of long if-block A over IPRLEV
      END IF

      EAV = 0.0d0

      Do kRoot = 1,nRoots
        EAV=EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
      End Do
      CASDFT_Funct=0.0d0
      If (KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM') Then
         If(nConf.eq.1.and.nActEl.eq.nac) Then
           Call Get_dScalar('CASDFT energy',CASDFT_Funct)
           ECAS = ECAS + CASDFT_Funct
           EAV=ECAS
         Else
           Call Get_dScalar('CASDFT energy',CASDFT_Funct)
           EAV=EAV+CASDFT_Funct - VIA_DFT - HALFQ
           ECAS = ECAS + CASDFT_Funct
         End If
      End If

      IF(IPRLEV.GE.USUAL .AND..NOT.lOPTO) THEN
* Start of long if-block B over IPRLEV

*----------------------------------------------------------------------*
*     Print final optimization conditions                              *
*----------------------------------------------------------------------*

      Write(LF,*)
      Line=''
      Write(Line(left-2:),'(A)') 'Final optimization conditions:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt2//'A)')'------------------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A,T45,F20.8)')
     &     'Average CI energy',EAV


      If (KeyCION.and.doDMRG)then
*        write(*,*)"KeyCION.and.doDMRG",KeyCION,doDMRG -- yma
*        If DMRG, no 2'-DMs thus omit ECAS and other properties
      else
        If (irlxroot.eq.0) Then
          Write(LF,Fmt2//'A,T45,F20.8)')
     &         'RASSCF energy',ECAS
        Else
          Write(LF,Fmt2//'A,I2,T45,F20.8)')
     &       'RASSCF energy for state ',irlxroot,ECAS
        End If

        Write(LF,Fmt2//'A,T45,F20.8)')
     &     'Super-CI energy',ESX
        Write(LF,Fmt2//'A,T45,F20.8)')
     &     'RASSCF energy change',DE
        Write(LF,Fmt2//'A,T50,ES10.3)')
     &     'Max change in MO coefficients',CMAX
        Write(LF,Fmt2//'A,T50,ES10.3)')
     &     'Max non-diagonal density matrix element',ROTMAX
        Write(LF,Fmt2//'A,T50,ES10.3)')
     &     'Maximum BLB matrix element',CBLBM
        Write(LF,Fmt2//'A,I4,A,I4,A,I4,A)')
     &     '(orbital pair',IBLBM,',',JBLBM,' in symmetry',ISYMBB,')'
        If (irlxRoot.ne.0)
     &  Write(LF,Fmt2//'A,T45,ES10.3)')
     &     'Norm of electronic gradient',RLXGRD

      End if

      If ( ISUPSM.ne.0 ) then
         Write(LF,Fmt2//'A)')
     &   'Supersymmetry has been used to disable selected orbital'//
     &   ' rotations'
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
     &   'The cleanup option has been used to set MO-coefficients'//
     &   ' explicitly to zero'
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
        Do kRoot = 1,lRoots
          CALL DDAFILE(JOBIPH,2,X4,NTOT2,IAD12)
          CALL DDAFILE(JOBIPH,2,X3,NTOT,IAD12)
          Call RelEne(Temp(1,kRoot),Temp(2,kRoot),nSym,nBas,
     &                X4,X3,X2,X1)
        End Do
        CALL mma_deallocate(X4)
        CALL mma_deallocate(X3)
        CALL mma_deallocate(X2)
        CALL mma_deallocate(X1)
      End If

      Call mma_allocate(EneTmp,lRoots,Label='EneTmp')
      If ( (nMVInt+nDCInt).ne.0) then
         Do i=1,lRoots
            Emv=Temp(1,i)
            Edc=Temp(2,i)
            Erel=ENER(I,ITER)+Emv+Edc
            EneTmp(i) = eRel
         End Do
      Else
         Do i=1,lRoots
            EneTmp(i) = ENER(I,ITER)+CASDFT_Funct-VIA_DFT-HALFQ
         End Do
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
         Do i=1,lRoots
            Emv=Temp(1,i)
            Edc=Temp(2,i)
            Erel=ENER(I,ITER)+Emv+Edc
            Write(LF,Fmt2//'I3,4(1X,F20.8))')
     &            I,ENER(I,ITER),Emv,Edc,Erel
         End Do
       Else
         Do i=1,lRoots
           Call PrintResult(LF,Fmt2//'A,I3,A,F16.8)',
     & 'RASSCF root number',I,' Total energy:',EneTmp(i),1)
         End Do
       End If
      ELSE IF(IPRLEV.GE.TERSE .AND..NOT.lOPTO) THEN
        Do i=1,lRoots
           Call PrintResult(LF,Fmt2//'A,I3,A,F16.8)',
     & 'RASSCF root number',I,' Total energy:',EneTmp(i),1)
        End Do
* End of long if-block C over IPRLEV
      END IF

      Call Store_Energies(lRoots,EneTmp,irlxroot)
      Call Put_iScalar('NumGradRoot',irlxroot)
      Call Put_dScalar('Average energy',EAV)
      Call mma_deallocate(EneTmp)

      iTol = Cho_X_GetTol(8)
      if(doDMRG) iTol = 6
*
      Line(1:8)='E_RASSCF'
      j=8
      Do i=1, nRoots
         If (nRoots.gt.1) Then
            If (i-1.lt.10) Then
               j=11
               Write (Line(9:j),'(A,I1,A)') '[',i-1,']'
            Else If (i-1.lt.100) Then
               j=12
               Write (Line(9:j),'(A,I2,A)') '[',i-1,']'
            Else If (i-1.lt.1000) Then
               j=13
               Write (Line(9:j),'(A,I3,A)') '[',i-1,']'
            Else If (i-1.lt.10000) Then
               j=14
               Write (Line(9:j),'(A,I4,A)') '[',i-1,']'
            End If
         Else
            j=8
         End If
         Call Add_Info(Line(1:j),ENER(iRoot(i),ITER),1,iTol)
      End Do
*---------------------------------------------------------------
* New JOBIPH layout: Also write hamiltonian matrix at IADR15(17):
*---------------------------------------------------------------
      CALL mma_allocate(HEFF,LROOTS,LROOTS,Label='HEFF')
      HEFF(:,:)=0.0D0
      DO J=1,LROOTS
       HEFF(J,J)=ENER(J,ITER)
      END DO
      IAD15=IADR15(17)
      CALL DDAFILE(JOBIPH,1,HEFF,LROOTS**2,IAD15)
      Call mma_deallocate(HEFF)
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
          If(kIvo) Then
            Call ivogen_rasscf(nSym,nBas,nFro,nIsh,nAsh,nTot2,nTot,
     &                         CMO,FDIAG)
            CALL PRIMO_RASSCF('Pseudonatural active'//
     &       ' orbitals and approximate occupation numbers + IVO,'//
     &       ' not suitable for CASPT2',
     &                         FDIAG,OCCN,CMO)
          Else
            CALL PRIMO_RASSCF('Pseudonatural active'//
     &                  ' orbitals and approximate occupation numbers',
     &                         FDIAG,OCCN,CMO)
          End If
        Else
           CALL PRIMO_RASSCF('All orbitals are'//
     &                  ' eigenfunctions of the PT2 Fock matrix',
     &                         FDIAG,OCCN,CMO)
        EndIf
      EndIf
* End of if-block D over IPRLEV
      END IF
*----------------------------------------------------------------------*
*     Also put on RUNFILE (in the future...):
*----------------------------------------------------------------------*
      Call Put_dArray('RASSCF orbitals',CMO,NTOT2)
      !! Fix https://molcasforum.univie.ac.at/viewtopic.php?id=1009
      IF (IPRLEV.LT.USUAL) Call Put_dArray('RASSCF OrbE',FDIAG,NTOT)
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

*                                                                      *
************************************************************************
*                                                                      *
*     Here follows a very long loop over KROOT:
*
*     But first save the 1st order density for gradients
*
      Call mma_allocate(DSave,nTot1,Label='DSave')
      Call Get_dArray_chk('D1AO',DSave,NTOT1)
*
*     The dipole moments will also be stored over all kroot states.
*
      Call mma_allocate(DM,3,Label='DM')
      Call mma_allocate(DMs,3,LROOTS,Label='DMs')
      DMs(:,:)=0.0D0
      Do_DM=.False.
*
      DO KROOT=1,LROOTS
*
* Read natural orbitals
        If ( NAC.GT.0 ) then
          CALL DDAFILE(JOBIPH,2,CMO,NTOT2,IAD12)
          CALL DDAFILE(JOBIPH,2,OCCN,NTOT,IAD12)
        End If
*
* Put the density matrix of this state on the runfile for
*  LoProp utility
        Call mma_allocate(DState,nTot1,Label='DState')
        DState(:)=0.0D0
        Call DONE_RASSCF(CMO,OCCN,DState)
        Call Put_dArray('D1ao',DState,NTOT1)
        Call mma_deallocate(DState)

        IF (IPRLEV.GE.USUAL) THEN
*
*       Start of if-block E over IPRLEV
*
        vNentropy=0.0d0
        Do i=1,NTOT ! very "inefficient" but the simplest
           xnu=OCCN(i)/2.0d0
           If (xnu.gt.0.0d0) vNentropy = vNentropy + xnu*log(xnu)
        End Do
        vNentropy = -vNentropy/log(2.0d0)
*
        Write(LF,'(6X,A,I3,A,F8.5)')
     *  'Von Neumann Entropy (Root ',KROOT,') = ',vNentropy
        Write(LF,*)
*
*       Compute Mulliken's population analysis
*
        Write(LF,'(/6X,A,I3)')
     *  'Mulliken population analysis for root number:',KROOT
        Write(LF,'(6X,A)')
     *  '-----------------------------------------------'
        Write(LF,*)
        lSave = KRoot.eq.iRlxRoot
        CALL CHARGE(nsym,nbas,BName,CMO,OCCN,SMAT,2,FullMlk,lSave)
        Write(LF,*)
*
*       Compute properties
*
        Write(LF,'(/6X,2A,I3)')
     *  'Expectation values of various properties ',
     *  'for root number:',KROOT
        Write(LF,'(6X,2A)')
     *  '-----------------------------------------',
     *  '------------------'
        Write(LF,*)
        END IF
* End of if-block E over IPRLEV

*
*     Write out info on a temporary vector file. Prpt
*     will need to read this.
        Note='Temporary orbital file used by prpt3.'
        LuTmp=50
        LuTmp=IsFreeUnit(LuTmp)
        Call WrVec('TMPORB',LuTmp,'CO',nSym,nBas,nBas,
     &            CMO,OCCN,Dum,iDum,Note)
        CALL PRPT()
*                                                                      *
************************************************************************
*       Store away the dipole moment of this state                     *
*
        Call Qpg_dArray('Dipole Moment',Do_DM,iDum(1))
        If (Do_DM) Then
*          Write (6,*) 'iRoot=',kRoot
           Call Get_dArray('Dipole Moment',DM,3)
*          Call RecPrt('Dipole Moment',' ',DM,1,3)
           Call DCopy_(3,DM,1,DMs(:,KROOT),1)
        End If
*                                                                      *
************************************************************************
*                                                                      *
*
*       Compute spin orbitals and spin population
*       (Note: this section overwrites the pseudo natural orbitals
*              with the spin orbitals).
*
* PAM2008: Only for at most MAXORBOUT orbitals. Default 10, reset by
*       keyword.
        IF(KROOT.LE.MAXORBOUT) THEN

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
     &    'Spin density matrix for root number:',KROOT
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
            Write(LF,'(/6X,A,I2)') 'symmetry species',ISYM
            Write(LF,*)
            CALL TRIPRT(' ',' ',X6(IND),NASH(ISYM))
            IND=IND+NASH(ISYM)*(NASH(ISYM)+1)/2
  50        CONTINUE
          END DO
         END IF
* End of long if-block F over IPRLEV
        END IF

*       Compute spin orbitals and spin population
        CALL DCOPY_(NTOT,[0.0D0],0,OCCN,1)
*SVC-11-01-2007 store original cmon in cmoso, which gets changed
        CALL mma_allocate(CMOSO,NTOT2,Label='CMOSO')
        CALL DCOPY_(NTOT2,CMON,1,CMOSO,1)

        if(.not.doDMRG)then
          CALL SPINORB(X6,CMOSO,OCCN,kroot)
        end if
        CALL mma_deallocate(X6)
        CALL DDAFILE(JOBIPH,1,CMOSO,NTOT2,IAD14)
        CALL DDAFILE(JOBIPH,1,OCCN,NTOT,IAD14)

        IF (IPRLEV.GE.USUAL) THEN
* Start of long if-block G over IPRLEV
         IF(ISPDEN .EQ. 1) THEN
          Write(LF,'(/6X,A,I3)')
     &    'Mulliken spin population analysis for root number:',KROOT
          Write(LF,'(6X,A)')
     &    '---------------------------------------------------'
          Write(LF,*)
          CALL CHARGE(nsym,nbas,BName,cmoso,OCCN,SMAT,3,FullMlk,
     &                .False.)
          Write(LF,*)
         ENDIF

*BOR0511
*
*       Do LoProp Charge analysis
*
         If (get_BasisType('ANO')) Then
           Write(LF,'(/6X,A,I3)')
     *     'LoProp population analysis for root number:',KROOT
           Write(LF,'(6X,A)')
     *     '-----------------------------------------------'
           Write(LF,*)
           Write(LF,*)
           Line=''
           Write(Line(left-2:),'(A)') 'LoProp analysis:'
           Call CollapseOutput(1,Line)
           Write(LF,Fmt2//'A)') '----------------'
           Write(LF,*)
           Write(LF,*)
*
* Ugly trick: use IRC to pass some info to LoProp ...
*
           iRC = Min((Abs(KRoot-iRlxRoot)),1)
           Call LoProp(iRC)
           Write(LF,*)
           Write(LF,'(6X,A,I3)')
     &     'Natural Bond Order analysis for root number:',KROOT
           Call Nat_Bond_Order(nSym,nBas,BName,2)
           Call CollapseOutput(0,'LoProp analysis:')
           Write(6,*)
         End If
* End of long if-block G over IPRLEV
         END IF
         CALL mma_deallocate(CMOSO)
cnf
*
*------- ESPF analysis
         Call DecideOnESPF(Do_ESPF)
         If (Do_ESPF) Call espf_analysis(lSave)
cnf

*PAM2008: End of restriction KROOT.LE.MAXORBOUT, see above.
        END IF
*
      END DO
*
*     Restore the correct 1st order density for gradient calculations.
*
      Call Put_dArray('D1ao',DSave,NTOT1)
      Call mma_deallocate(DSave)
*
*     Save the list of dipole moments on the run file.
*
      If (Do_DM)
     &   Call Put_dArray('Last Dipole Moments',DMs,3*LROOTS)
*     Call RecPrt('Last Dipole Moments',' ',DM),3,LROOTS)
      Call mma_deallocate(DM)
      Call mma_deallocate(DMs)
*                                                                      *
************************************************************************
*                                                                      *
      CALL mma_deallocate(CMON)
*----------------------------------------------------------------------*

      End Subroutine OutCtl
