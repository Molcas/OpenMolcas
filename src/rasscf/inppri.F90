!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine InpPri(lOPTO)
!***********************************************************************
!                                                                      *
!     Echo input                                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use OneDat, only: sNoOri
use fcidump, only: DumpOnly
use fciqmc, only: DoNECI
use CC_CI_mod, only: Do_CC_CI
use Fock_util_global, only: DoLocK
use Functionals, only: Init_Funcs, Print_Info
use KSDFT_Info, only: CoefR, CoefX
use rctfld_module, only: lRF
use gas_data, only: iDoGAS, NGAS, NGSSH, IGSOCCX
use rasscf_global, only: KSDFT, DoBlockDMRG, DoDMRG, ICICH, ICIRST, iPCMRoot, iRLXRoot, iSupSM, ITMAX, l_casdft, lRoots, lSquare, &
                         LvShft, MAXIT, n_Det, NAC, NFR, NIN, NONEQ, NROOTS, NSEC, nTit, RFPert, ThrE, ThrSX, ThrTE, Tot_Charge, &
                         Tot_El_Charge, Tot_Nuc_Charge, Title, Header, iRoot, Weight, iCI, cCI, ixSym, NQUNE
#ifdef _DMRG_
use qcmaquis_interface_cfg
use qcmaquis_interface_utility_routines, only: print_dmrg_info
#endif
#ifdef _ENABLE_DICE_SHCI_
use rasscf_global, only: dice_eps1, dice_eps2, dice_iter, dice_Restart, dice_SampleN, Dice_Stoc, nRef_Dice, diceocc
#endif
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
use rasscf_global, only: MXDMRG, ChemPS2_blb, ChemPS2_lreStart, ChemPS2_Noise, ChemPS2_Restart, Davidson_tol, Do3RDM, HFOcc, &
                         Max_canonical, Max_Sweep
#endif
use SplitCas_Data, only: DoSPlitCas, MxIterSplit, ThrSplit, lRootSplit, NumSplit, iDimBlockA, EnerSplit, GapSpli, PerSplit, &
                         PerCSpli, fOrdSplit
use PrintLevel, only: USUAL, SILENT
use output_ras, only: LF, IPRLOC
use general_data, only: NACTEL, NHOLE1, NELEC3, ISPIN, STSYM, NSYM, NSEL, NTOT1, NASH, NBAS, NDEL, NFRO, NISH, NRS1, NRS2, NRS3, &
                        NSSH
use spinfo, only: DoComb, NCNFTP, NCSASM, NDTASM, NDTFTP, I_ELIMINATE_GAS_MOLCAS, NCSF_HEXS
use DWSol, only: DWSolv, DWSol_fixed, W_SOLV
use RASDim, only: MxRef
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero

implicit none
logical lOPTO
character(len=8) Fmt1, Fmt2, Label
character(len=120) Line, BlLine, StLine
character(len=3) lIrrep(8)
character(len=2) GASidx
character(len=80) KSDFT2
logical DoCholesky
real*8, allocatable :: Tmp0(:)
real*8 AvailMB, WillNeedMB
integer i, iCharge, iComp, iDoRI, iEnd, iGAS, iOpt, iPrLev, iRC, iRef, iStart, iSyLbl, iSym, iTemp, j, left, lLine, lPaper, &
        MaxRem, n_paired_elec, n_unpaired_elec, nLine
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
character(len=3) SNAC
integer iHFOcc
#endif
#ifdef _DMRG_
character(len=100) :: dmrg_start_guess
#endif
#ifdef _ENABLE_DICE_SHCI_
integer iRef_Dice
#endif

! Print level:
IPRLEV = IPRLOC(1)
!----------------------------------------------------------------------*
!     Start and define the paper width                                 *
!----------------------------------------------------------------------*
lPaper = 132
!----------------------------------------------------------------------*
!     Initialize blank and header lines                                *
!----------------------------------------------------------------------*
lLine = len(Line)
do i=1,lLine
  BlLine(i:i) = ' '
  StLine(i:i) = '*'
end do
lPaper = 132
left = (lPaper-lLine)/2
write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
if (IPRLEV == SILENT) goto 900
!----------------------------------------------------------------------*
!     Print the project title                                          *
!----------------------------------------------------------------------*
if (IPRLEV >= USUAL) then
  if (nTit > 0) then
    write(LF,*)
    nLine = nTit+5
    do i=1,nLine
      Line = BlLine
      if ((i == 1) .or. (i == nLine)) Line = StLine
      if (i == 3) Line = 'Project:'
      if ((i >= 4) .and. (i <= nLine-2)) write(Line,'(A72)') Title(i-3)
      call Center_Text(Line)
      write(LF,Fmt1) '*'//Line//'*'
    end do
    write(LF,*)
  end if
end if
!----------------------------------------------------------------------*
!     Print the ONEINT file identifier                                 *
!----------------------------------------------------------------------*
if ((IPRLEV >= USUAL) .and. (.not. lOPTO)) then
  write(LF,*)
  write(LF,Fmt1) 'Header of the ONEINT file:'
  write(LF,Fmt1) '--------------------------'
  write(Line,'(36A2)') (Header(i),i=1,36)
  write(LF,Fmt1) trim(adjustl(Line))
  write(Line,'(36A2)') (Header(i),i=37,72)
  write(LF,Fmt1) trim(adjustl(Line))
  write(LF,*)
  !--------------------------------------------------------------------*
  !     Print the status of ORDINT                                     *
  !--------------------------------------------------------------------*
  write(LF,*)
  if (lSquare) then
    write(LF,Fmt1) 'OrdInt status: squared'
  else
    write(LF,Fmt1) 'OrdInt status: non-squared'
  end if
  write(LF,*)
  !--------------------------------------------------------------------*
  !     Print cartesian coordinates of the system                      *
  !--------------------------------------------------------------------*
  call PrCoor()
end if
!----------------------------------------------------------------------*
!     Print orbital and wavefunction specifications                    *
!----------------------------------------------------------------------*
if ((IPRLEV >= USUAL) .and. (.not. lOPTO)) then
  write(LF,*)
  Line = ' '
  write(Line(left-2:),'(A)') 'Wave function specifications:'
  call CollapseOutput(1,Line)
  write(LF,Fmt1) '-----------------------------'
  write(LF,*)
  if (NFR > 0) write(LF,Fmt2//'A,T45,I6)') 'Number of frozen shell electrons',2*NFR
  write(LF,Fmt2//'A,T45,I6)') 'Number of closed shell electrons',2*NIN
  write(LF,Fmt2//'A,T45,I6)') 'Number of electrons in active shells',NACTEL
  if (.not. idogas) then
    !.. for RAS
    write(LF,Fmt2//'A,T45,I6)') 'Max number of holes in RAS1 space',NHOLE1
    write(LF,Fmt2//'A,T45,I6)') 'Max nr of electrons in RAS3 space',NELEC3
  else
    !.. for GAS
    do IGAS=1,NGAS
      write(GASidx,'(i0)') IGAS
      write(LF,Fmt2//'A,A,A,T45,2I6)') 'Min/Max nr of electrons up to GAS',GASidx,' sp.',igsoccx(igas,1),igsoccx(igas,2)
    end do
  end if

  if (NFR > 0) write(LF,Fmt2//'A,T45,I6)') 'Number of frozen orbitals',NFR
  write(LF,Fmt2//'A,T45,I6)') 'Number of inactive orbitals',NIN
  write(LF,Fmt2//'A,T45,I6)') 'Number of active orbitals',NAC
  write(LF,Fmt2//'A,T45,I6)') 'Number of secondary orbitals',NSEC
  write(LF,Fmt2//'A,T45,F6.1)') 'Spin quantum number',(dble(ISPIN-1))/2.0d0
  write(LF,Fmt2//'A,T45,I6)') 'State symmetry',STSYM
  call CollapseOutput(0,'Wave function specifications:')

  call Get_cArray('Irreps',lIrrep,24)
  do iSym=1,nSym
    lIrrep(iSym) = adjustr(lIrrep(iSym))
  end do

  write(LF,*)
  Line = ' '
  write(Line(left-2:),'(A)') 'Orbital specifications:'
  call CollapseOutput(1,Line)
  write(LF,Fmt1) '-----------------------'
  write(LF,*)
  write(LF,Fmt2//'A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
  write(LF,Fmt2//'A,T47,8(1X,A))') '                ',(lIrrep(iSym),iSym=1,nSym)
  write(LF,Fmt2//'A,T47,8I4)') 'Frozen orbitals',(nFro(iSym),iSym=1,nSym)
  write(LF,Fmt2//'A,T47,8I4)') 'Inactive orbitals',(nIsh(iSym),iSym=1,nSym)
  write(LF,Fmt2//'A,T47,8I4)') 'Active orbitals',(nAsh(iSym),iSym=1,nSym)
  if (.not. iDoGas) then
    write(LF,Fmt2//'A,T47,8I4)') 'RAS1 orbitals',(nRs1(iSym),iSym=1,nSym)
    write(LF,Fmt2//'A,T47,8I4)') 'RAS2 orbitals',(nRs2(iSym),iSym=1,nSym)
    write(LF,Fmt2//'A,T47,8I4)') 'RAS3 orbitals',(nRs3(iSym),iSym=1,nSym)
  else
    do IGAS=1,NGAS
      write(LF,Fmt2//'A,I1,A,T47,8I4)') 'GAS',IGAS,' orbitals',(ngssh(igas,iSym),iSym=1,nSym)
    end do
  end if

  write(LF,Fmt2//'A,T47,8I4)') 'Secondary orbitals',(nSsh(iSym),iSym=1,nSym)
  write(LF,Fmt2//'A,T47,8I4)') 'Deleted orbitals',(nDel(iSym),iSym=1,nSym)
  write(LF,Fmt2//'A,T47,8I4)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
  call CollapseOutput(0,'Orbital specifications:')
  write(LF,*)

# if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
  if (.not. DoBlockDMRG) goto 113

# ifdef _ENABLE_DICE_SHCI_
  Line = ' '
  write(Line(left-2:),'(A)') 'DICE specifications:'
  call CollapseOutput(1,Line)
  write(LF,Fmt1) '--------------------------'
  write(LF,*)
  write(LF,Fmt2//'A,T70,L6)') 'Heat-bath configuration interaction (JCTC, 2017, 13, 1595)',DoBlockDMRG
  write(LF,Fmt2//'A,T45,L6)') 'Semistochastic algorithm',Dice_stoc
  write(LF,Fmt2//'A,T45,L6)') 'Full restart',dice_restart
  write(LF,Fmt2//'A,T45,I6)') 'Max iterations',dice_iter
  write(LF,Fmt2//'A,T45,ES10.3)') 'Epsilon1',dice_eps1
  write(LF,Fmt2//'A,T45,ES10.3)') 'Epsilon2',dice_eps2
  write(LF,Fmt2//'A,T45,I6)') 'SampleN',dice_sampleN
  write(LF,Fmt2//'A,T45)') 'Occupation guess'
  do iref_dice=1,nref_dice
    write(LF,Fmt2//'A)') trim(diceocc(iref_dice))
  end do
  call CollapseOutput(0,'DICE specifications:')

  ! Skip printing CI specifications in DICE
  goto 114
# endif

# if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
  Line = ' '
  write(Line(left-2:),'(A)') 'DMRG sweep specifications:'
  call CollapseOutput(1,Line)
  write(LF,Fmt1) '--------------------------'
  write(LF,*)
  write(LF,Fmt2//'A,T45,I6)') 'Number of renormalized basis',MxDMRG
  write(LF,Fmt2//'A,T45,I6)') 'Number of root(s) required',NROOTS

  write(LF,Fmt2//'A,T45,I6)') 'Maximum number of sweeps',max_sweep
  write(LF,Fmt2//'A,T45,I6)') 'Maximum number of sweeps in RDM',max_canonical
  write(LF,Fmt2//'A,T45,ES10.3)') 'Threshold for restarting',chemps2_blb
  write(LF,Fmt2//'A,T45,ES10.3)') 'Minimum Davidson tolerance',davidson_tol
  write(LF,Fmt2//'A,T45,ES10.3)') 'DMRG convergence threshold',THRE/2.0
  write(LF,Fmt2//'A,T45,ES10.3)') 'Noise prefactor',chemps2_noise
  write(LF,Fmt2//'A,T45,L6)') 'Restart from previous calculation',chemps2_restart
  write(LF,Fmt2//'A,T45,L6)') 'Calculate 3-RDM and F.4-RDM',Do3RDM
  write(LF,Fmt2//'A,T45,I6)') 'Restart scheme in 3-RDM and F.4-RDM',chemps2_lrestart
  write(SNAC,'(I3)') NAC
  write(LF,Fmt2//'A,T45,'//trim(adjustl(SNAC))//'I2)') 'Occupation guess',(HFOCC(ihfocc),ihfocc=1,NAC)

  ! NN.14 FIXME: haven't yet checked whether geometry opt. works correctly with DMRG
  write(LF,Fmt2//'A,T45,I6)') 'Root chosen for geometry opt.',IRLXROOT
  call CollapseOutput(0,'DMRG sweep specifications:')

  ! Skip printing CI specifications in DMRG-CASSCF
  goto 114

113 continue
# endif
# endif

  Line = ' '
  if (doDMRG) then
    write(Line(left-2:),'(A)') 'DMRG specifications:'
  else
    write(Line(left-2:),'(A)') 'CI expansion specifications:'
  end if
  call CollapseOutput(1,Line)
  write(LF,Fmt1) '----------------------------'
  write(LF,*)

  if (doDMRG) then  !> Information for QCMaquis-DMRG
#   ifdef _DMRG_
    if (dmrg_orbital_space%initial_occ(1,1) > 0) then
      dmrg_start_guess = 'Single determinant'
    else if (dmrg_warmup%doCIDEAS) then
      dmrg_start_guess = 'CI-DEAS'
    else
      dmrg_start_guess = 'Random numbers (default)'
    end if
    call print_dmrg_info(lf,fmt2,1,dmrg_start_guess,nroots,thre)
#   endif
  else
    write(LF,Fmt2//'A,T40,I11)') 'Number of CSFs',NCSASM(STSYM)
    if (I_ELIMINATE_GAS_MOLCAS > 0) write(LF,Fmt2//'A,T40,I11)') 'Number of highly excited CSFs',nCSF_HEXS
    if (DoComb) then
      write(LF,Fmt2//'A,T40,I11)') 'Number of spin combinations',NDTASM(STSYM)
      write(LF,Fmt2//'A,T40,I11)') 'Number of determinants',2*NDTASM(STSYM)-NDTFTP(1)*NCNFTP(1,STSYM)
    else
      write(LF,Fmt2//'A,T40,I11)') 'Number of determinants',NDTASM(STSYM)
    end if
  end if
  n_Det = 2
  n_unpaired_elec = (iSpin-1)
  n_paired_elec = nActEl-n_unpaired_elec
  if ((n_unpaired_elec+n_paired_elec/2 == nac) .or. (NDTASM(STSYM) == 1)) n_Det = 1
  if (KSDFT == 'DIFF') n_Det = 1
  if (KSDFT == 'ROKS') n_Det = 1

  if (.not. DoSplitCAS) then  ! GLMJ
    write(LF,Fmt2//'A,T45,I6)') 'Number of root(s) required',NROOTS
    if (irlxroot /= 0) write(LF,Fmt2//'A,T45,I6)') 'Root chosen for geometry opt.',IRLXROOT
    if (ICICH == 0) then

      if (nRoots == 1) then

        if (doDMRG) then
          write(LF,Fmt2//'A,(T45,10I6))') 'DMRG root used',IROOT(1)
        else
          write(LF,Fmt2//'A,(T45,10I6))') 'CI root used',IROOT(1)
        end if
        write(LF,Fmt1) '   '

      else

        if (doDMRG) then
          write(LF,Fmt2//'A,(T45,10I6))') 'DMRG roots used',(IROOT(i),i=1,nRoots)
        else
          write(LF,Fmt2//'A,(T45,10I6))') 'CI roots used',(IROOT(i),i=1,nRoots)
        end if
        write(LF,Fmt2//'A,(T45,10F6.3))') 'weights',(Weight(i),i=1,nRoots)
      end if

    else
      do i=1,nRoots
        write(LF,Fmt2//'A,T45,I6)') 'selected root',iRoot(i)
        write(LF,Fmt2//'A,T45,10I6)') 'Reference configurations',(iCI(i,iRef),iRef=1,mxRef)
        write(LF,Fmt2//'A,T45,10F6.3)') 'CI-coeff',(cCI(i,iRef),iRef=1,mxRef)
      end do
    end if
    if (.not. doDMRG) then
      write(LF,Fmt2//'A,T45,I6)') 'highest root included in the CI',LROOTS
      write(LF,Fmt2//'A,T45,I6)') 'max. size of the explicit Hamiltonian',NSEL
    end if
  else
    write(LF,Fmt2//'A,T45,I6)') 'Root required ',lrootSplit
    if (NumSplit) write(LF,Fmt2//'A,T45,I6)') 'A-Block Size in SplitCAS',iDimBlockA
    if (EnerSplit) write(LF,Fmt2//'A,T44,F7.2)') 'Energy Gap (eV) in SplitCAS',GapSpli
    if (PerSplit) write(LF,Fmt2//'A,T44,F7.1)') 'Percentage sought in SplitCAS',PercSpli
    write(LF,Fmt2//'A,T45,ES10.3)') 'Threshold for SplitCAS',ThrSplit
    !write(LF,Fmt2//'A,T49,G10.3)')'Thrs over the root to be opt in SplitCAS', ThrSplit
    write(LF,Fmt2//'A,T47,I4)') 'Maximum number of SplitCAS iterations',MxIterSplit
    if (FordSplit) then
      write(LF,Fmt2//'A,T47)') 'CI coeff: 1st-order approximation'
    else
      write(LF,Fmt2//'A,T47)') 'CI coeff: 0th-order approximation'
    end if
  end if
  call CollapseOutput(0,'CI expansion specifications:')

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
114 continue
#endif

end if
if (lOPTO) then
  write(LF,*)
  Line = ' '
  write(Line(left-2:),'(A)') 'RASSCF input specifications:'
  call CollapseOutput(1,Line)
  write(LF,Fmt1) '----------------------------'
  if (.not. DoSplitCAS) then
    write(LF,Fmt2//'A,T45,I6)') 'Number of root(s) required',NROOTS
    if (irlxroot /= 0) write(LF,Fmt2//'A,T45,I6)') 'Root chosen for geometry opt.',IRLXROOT
    if (ICICH == 0) then
      if (nRoots == 1) then
        write(LF,Fmt2//'A,(T45,10I6))') 'CI root used',IROOT(1)
      else
        write(LF,Fmt2//'A,(T45,10I6))') 'CI roots used',(IROOT(i),i=1,nRoots)
        write(LF,Fmt2//'A,(T45,10F6.3))') 'weights',(Weight(i),i=1,nRoots)
      end if
    else
      do i=1,nRoots
        write(LF,Fmt2//'A,T45,I6)') 'selected root',iRoot(i)
        write(LF,Fmt2//'A,T45,10I6)') 'Reference configurations',(iCI(i,iRef),iRef=1,mxRef)
        write(LF,Fmt2//'A,T45,10F6.3)') 'CI-coeff',(cCI(i,iRef),iRef=1,mxRef)
      end do
    end if
  end if
  call CollapseOutput(0,'RASSCF input specifications:')
end if
! Check that the user doesn't try to calculate more roots than it's possible
! NN.14 FIXME: in DMRG-CASSCF, skip this check for the time
!              since Block DMRG code will check this internally
!if (NROOTS > NCSASM(STSYM)) then
if ((.not. any([DoNECI,Do_CC_CI,DumpOnly,doDMRG,doBlockDMRG])) .and. (NROOTS > NCSASM(STSYM))) then
  write(LF,*) '************ ERROR ***********'
  write(LF,*) ' You can''t ask for more roots'
  write(LF,*) ' than there are configurations '
  write(LF,*) '******************************'
  write(LF,*)
  call Quit_OnUserError()
end if
! If the calculation will be too big:
call mma_MaxDBLE(MaxRem)
WillNeedMB = (8.0d0*1.50d0*6.0d0*NDTASM(STSYM)/1.048d6)
AvailMB = (8.0d0*MaxRem/1.048d6)
if (WillNeedMB > AvailMB) then
  write(6,*)
  write(6,*) ' *************************************************'
  write(6,*) ' Sorry, but your calculation will probably be too'
  write(6,*) ' large for the available memory.'
  write(6,*) ' The number of determinants is ',NDTASM(STSYM)
  write(6,*) ' During CI equation solution, there will be'
  write(6,*) ' up to six vectors of this size in memory.'
  write(6,*) ' We estimate an additional 50% for other stuff.'
  write(6,*)
  write(6,'(A,F9.1,A)') ' This alone will need at least ',WillNeedMB,' MB.'
  write(6,'(A,F9.1,A)') ' Available at this point is ',AvailMB,' MB.'
  write(6,*) ' Please increase MOLCAS_MEM, and try again.'
  write(6,*) ' *************************************************'
  write(6,*)
  call Quit_OnUserError()
end if

if ((IPRLEV >= USUAL) .and. (.not. lOPTO)) then
  write(LF,*)
  Line = ' '
  write(Line(left-2:),'(A)') 'Optimization specifications:'
  call CollapseOutput(1,Line)
  write(LF,Fmt1) '----------------------------'
  write(LF,*)
  call DecideOnCholesky(DoCholesky)
  if (DoCholesky) then
    call Get_iScalar('System BitSwitch',iDoRI)
    if (iand(iDoRI,1024) == 1024) then
      if (DoLocK) then
        write(LF,Fmt2//'A,T45,I6)') 'RASSCF algorithm: LK RI/DF'
      else
        write(LF,Fmt2//'A,T45,I6)') 'RASSCF algorithm: RI/DF'
      end if
    else
      if (DoLocK) then
        write(LF,Fmt2//'A,T45,I6)') 'RASSCF algorithm: LK Cholesky'
      else
        write(LF,Fmt2//'A,T45,I6)') 'RASSCF algorithm: Cholesky'
      end if
    end if
  else
    write(LF,Fmt2//'A,T45,I6)') 'RASSCF algorithm: Conventional'
  end if
  !*********************************************************************
  ! Some printout for mcpdft method
  !*********************************************************************
  KSDFT2 = KSDFT
  if (l_casdft) then
    KSDFT2 = KSDFT(index(KSDFT,'T:')+2:)
    write(LF,Fmt2//'A)') 'This is a MC-PDFT calculation with functional: '//KSDFT
    write(LF,Fmt2//'A,T45,ES10.3)') 'Exchange scaling factor',CoefX
    write(LF,Fmt2//'A,T45,ES10.3)') 'Correlation scaling factor',CoefR
  end if
  !*********************************************************************

  write(LF,Fmt2//'A,T45,I6)') 'Maximum number of macro iterations',MAXIT
  write(LF,Fmt2//'A,T45,I6)') 'Maximum number of SX iterations',ITMAX
  write(LF,Fmt2//'A,T45,ES10.3)') 'Threshold for RASSCF energy',THRE
  call Put_dScalar('EThr',ThrE)
  write(LF,Fmt2//'A,T45,ES10.3)') 'Threshold for max MO rotation',THRTE
  write(LF,Fmt2//'A,T45,ES10.3)') 'Threshold for max BLB element',THRSX
  write(LF,Fmt2//'A,T45,ES10.3)') 'Level shift parameter',LVSHFT
  if (NQUNE /= 0) write(LF,Fmt1) 'Make Quasi-Newton update'
  if (ISUPSM /= 0) then
    write(LF,Fmt1) 'Supersymmetry is used to disable selected orbital rotations'
    iEnd = 0
    do iSym=1,nSym
      iStart = iEnd+1
      iEnd = iEnd+nBas(iSym)
      iTemp = 0
      do i=iStart,iEnd
        iTemp = iTemp+IXSYM(i)
      end do
      if (iTemp > 0) then
        write(LF,Fmt2//'A,I3)') 'Supersymmetry vector for symmetry species',iSym
        write(LF,Fmt2//'30I3)') (IXSYM(i),i=iStart,iEnd)
      end if
    end do
  end if
  if (lRF) then
    call mma_allocate(Tmp0,nTot1+4,Label='Tmp0')
    iRc = -1
    iOpt = ibset(0,sNoOri)
    iComp = 1
    iSyLbl = 1
    Label = 'Mltpl  0'
    call RdOne(iRc,iOpt,Label,iComp,Tmp0,iSyLbl)
    Tot_Nuc_Charge = Tmp0(nTot1+4)
    if (iRc /= 0) then
      write(LF,*) 'InpPri: iRc from Call RdOne not 0'
      write(LF,*) 'Label = ',Label
      write(LF,*) 'iRc = ',iRc
      call Abend()
    end if
    call mma_deallocate(Tmp0)
    Tot_El_Charge = Zero
    do iSym=1,nSym
      Tot_El_Charge = Tot_El_Charge-2.0d0*dble(nFro(iSym)+nIsh(iSym))
    end do
    Tot_El_Charge = Tot_El_Charge-dble(nActEl)
    Tot_Charge = Tot_Nuc_Charge+Tot_El_Charge
    iCharge = int(Tot_Charge)
    call PrRF(.false.,NonEq,iCharge,2)

    if (DWSolv%DWZeta == -12345d+00) then
      write(LF,Fmt2//'A)') 'Weights of the reaction field are specified by RFROOT'
      write(LF,Fmt2//'(T45,10F6.3))') (W_SOLV(i),i=1,nRoots)
    else if (DWSolv%DWZeta < 0.0d+00) then
      call DWSol_fixed(i,j)
      if ((i == 0) .and. (j == 0)) then
        write(LF,Fmt2//'A)') 'Unrecognized negative DWZeta (DWSOl)'
        write(LF,Fmt2//'A,T51,A)') 'Dynamically weighted solvation is ','automatically turned off!'
      else
        write(LF,Fmt2//'A,T45,I2,X,I2)') 'Reaction field from states:',i,j
        if (max(i,j) > nRoots) then
          write(LF,Fmt2//'A)') 'The specified state is too high! Cannot proceed...'
          call Quit_OnUserError()
        end if
      end if
    else if (IPCMROOT <= 0) then
      write(LF,Fmt2//'A,T44,A)') 'Reaction field from state:',' State-Averaged'
      if (DWSolv%DWZeta /= 0.0d+00) then
        write(LF,Fmt2//'A,T51,A)') 'Dynamically weighted solvation is ','automatically turned off!'
        DWSolv%DWZeta = 0.0d+00
      end if
    else
      write(LF,Fmt2//'A,T44,I2)') 'Reaction field from state:',IPCMROOT
      if (DWSolv%DWZeta > 0.0d+00) write(LF,Fmt2//'A,ES10.3,A,I1,A)') 'Dynamically weighted solvation is used with DWSOlv = ', &
                                                                      DWSolv%DWZeta,' (DWTYpe = ',DWSolv%DWType,')'
    end if
  end if
  !if (DWSCF%do_DW) then
  !  write(LF,Fmt2//'A)') 'Dynamically weighted MCSCF is enabled'
  !  write(LF,Fmt2//'A,T44,I2)') 'Target state:', DWSCF%DWRoot
  !  write(LF,Fmt2//'A,ES10.3,A,I1,A)') 'Dynamically weighted MCSCF is used with DWZEta = ',DWSCF%DWZeta,' (DWTYpe = ', &
  !                                     DWSCF%DWType,')'
  !end if
  call CollapseOutput(0,'Optimization specifications:')
  if (RFpert) then
    write(LF,*)
    write(LF,Fmt1) 'Reaction field specifications:'
    write(LF,Fmt1) '------------------------------'
    write(LF,*)
    write(LF,'(6X,A)') 'The Reaction field is added as a perturbation and has been determined in a previous calculation'
    write(LF,*)
  end if
  if (ICIRST == 1) write(LF,Fmt1) 'Starting CI array(s) will be read from file'
  call Put_dScalar('EThr',ThrE)

  ! Print out grid information in case of DFT

  if (KSDFT /= 'SCF') then
    call Put_dScalar('DFT exch coeff',CoefX)
    call Put_dScalar('DFT corr coeff',CoefR)
    call Funi_Print()
    if (IPRLEV >= USUAL) then
      write(6,*)
      write(6,'(6X,A)') 'DFT functional specifications'
      write(6,'(6X,A)') '-----------------------------'
      call libxc_version()
      call Init_Funcs(KSDFT2)
      call Print_Info()
      write(6,*)
    end if
  end if
end if
write(LF,*)

900 continue
call XFlush(LF)
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine InpPri
