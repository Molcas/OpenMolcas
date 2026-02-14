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

subroutine OutCtlSplit(CMO,OCCN,SMAT,lOPTO)
!***********************************************************************
!                                                                      *
!     Control section for the final output                             *
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

use OneDat, only: sNoOri, sOpSiz
use rctfld_module, only: lRF
use rasscf_global, only: BName, CBLBM, CMax, DE, ECAS, Ener, ESX, FDIAG, HALFQ, iADR15, IBLBM, iPCMRoot, IPT2, iRLXRoot, iSPDen, &
                         iSupSM, iSymBB, ITER, ixSym, jBLBM, KSDFT, NAC, NACPAR, NACPR2, NIN, NONEQ, nRoots, NSEC, OutFmt1, &
                         RFPert, RlxGrd, RotMax, Tot_Charge, Tot_El_Charge, Tot_Nuc_Charge, Via_DFT
use SplitCas_Data, only: EnerSplit, GapSpli, iDimBlockA, lRootSplit, MxIterSplit, PerCSpli, PerSplit, ThrSplit
use PrintLevel, only: DEBUG, TERSE, USUAL, VERBOSE
use output_ras, only: IPRLOC
use general_data, only: CleanMask, ISPIN, JOBIPH, NACTEL, NASH, NBAS, NCONF, NDEL, NELEC3, NFRO, NHOLE1, NISH, NRS1, NRS2, NRS3, &
                        NSSH, NSYM, NTOT, NTOT1, NTOT2, STSYM
use spinfo, only: NCSASM, NDTASM
use DWSol, only: DWSol_fixed, DWSolv, W_SOLV
use Molcas, only: MxRoot
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CMO(*), OCCN(*), SMAT(*)
logical(kind=iwp) :: lOPTO
integer(kind=iwp) :: i, IAD03, IAD12, IAD14, IAD15, iCharge, iComp, iDimN, iDimO, iDImV, iDum(56), iEnd, Ind, iOpt, iPrLev, iRC, &
                     iRC1, iRC2, iStart, iSyLbl, iSym, iTemp, iTol, j, kRoot, left, LuTmp, NAO, nDCInt, nMVInt, NO
real(kind=wp) :: CASDFT_Funct, Dum(1), EAV, Edc, Emv, EneTmp, Erel, percent, Temp(2,mxRoot)
logical(kind=iwp) :: Do_ESPF, FullMlk, get_BasisType, lSave
character(len=120) Line
character(len=80) Note
character(len=8) Fmt2, Label
character(len=3) lIrrep(8)
real(kind=wp), allocatable :: CMON(:), CMOSO(:), D(:), DSave(:), Tmp0(:), X1(:), X2(:), X3(:), X4(:), X6(:)
integer, external :: Cho_X_GetTol, IsFreeUnit

!----------------------------------------------------------------------*
!     Start and define the paper width                                 *
!----------------------------------------------------------------------*
! Local print level (if any)
IPRLEV = IPRLOC(6)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering OUTCTLSPLIT'

! Additional DFT correlation energy, if any
CASDFT_Funct = Zero
if ((KSDFT /= 'SCF') .and. (KSDFT /= 'PAM')) call Get_dScalar('CASDFT energy',CASDFT_Funct)
!----------------------------------------------------------------------*
!     Initialize blank and header lines                                *
!----------------------------------------------------------------------*
left = 6
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'

if ((IPRLEV >= USUAL) .and. (.not. lOPTO)) then
  ! Start of long if-block A over IPRLEV
  !--------------------------------------------------------------------*
  !     Print orbital and wavefunction specifications                  *
  !--------------------------------------------------------------------*
  write(u6,*)
  Line = ''
  write(Line(left-2:),'(A)') 'Wave function specifications:'
  call CollapseOutput(1,Line)
  write(u6,Fmt2//'A)') '-----------------------------'
  write(u6,*)
  write(u6,Fmt2//'A,T45,I6)') 'Number of closed shell electrons',2*NIN
  write(u6,Fmt2//'A,T45,I6)') 'Number of electrons in active shells',NACTEL
  write(u6,Fmt2//'A,T45,I6)') 'Max number of holes in RAS1 space',NHOLE1
  write(u6,Fmt2//'A,T45,I6)') 'Max nr of electrons in RAS3 space',NELEC3
  write(u6,Fmt2//'A,T45,I6)') 'Number of inactive orbitals',NIN
  write(u6,Fmt2//'A,T45,I6)') 'Number of active orbitals',NAC
  write(u6,Fmt2//'A,T45,I6)') 'Number of secondary orbitals',NSEC
  write(u6,Fmt2//'A,T45,F6.1)') 'Spin quantum number',Half*real(ISPIN-1,kind=wp)
  write(u6,Fmt2//'A,T45,I6)') 'State symmetry',STSYM
  call CollapseOutput(0,'Wave function specifications:')

  call Get_cArray('Irreps',lIrrep,24)
  do iSym=1,nSym
    lIrrep(iSym) = adjustr(lIrrep(iSym))
  end do

  write(u6,*)
  Line = ''
  write(Line(left-2:),'(A)') 'Orbital specifications:'
  call CollapseOutput(1,Line)
  write(u6,Fmt2//'A)') '-----------------------'
  write(u6,*)
  write(u6,Fmt2//'A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8(1X,A))') '                ',(lIrrep(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Frozen orbitals',(nFro(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Inactive orbitals',(nIsh(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Active orbitals',(nAsh(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'RAS1 orbitals',(nRs1(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'RAS2 orbitals',(nRs2(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'RAS3 orbitals',(nRs3(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Secondary orbitals',(nSsh(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Deleted orbitals',(nDel(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
  call CollapseOutput(0,'Orbital specifications:')
  write(u6,*)
  Line = ''
  write(Line(left-2:),'(A)') 'CI expansion specifications:'
  call CollapseOutput(1,Line)
  write(u6,Fmt2//'A)') '----------------------------'
  write(u6,*)
  write(u6,Fmt2//'A,T40,I11)') 'Number of CSFs',NCSASM(STSYM)
  write(u6,Fmt2//'A,T40,I11)') 'Number of determinants',NDTASM(STSYM)
  write(u6,Fmt2//'A,T45,I6)') 'Root required ',lrootSplit
  if (EnerSplit) write(u6,Fmt2//'A,T44,F7.2)') 'Energy Gap (eV) in SplitCAS',GapSpli
  if (PerSplit) write(u6,Fmt2//'A,T44,F7.1)') 'Percentage sought in SplitCAS',PercSpli
  percent = real(iDimBlockA,kind=wp)/real(NCSASM(STSYM),kind=wp)*100.0_wp
  write(u6,Fmt2//'A,T42,I9,A,F5.1,A)') 'A-Block Size in SplitCAS (CSFs)',iDimBlockA,' (',percent,' %)'
  write(u6,Fmt2//'A,T45,ES10.3)') 'Threshold for SplitCAS',ThrSplit
  write(u6,Fmt2//'A,T47,I4)') 'Maximum number of SplitCAS iterations',MxIterSplit
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
      write(u6,*) 'OutCtl: iRc from Call RdOne not 0'
      write(u6,*) 'Label = ',Label
      write(u6,*) 'iRc = ',iRc
      call Abend()
    end if
    call mma_deallocate(Tmp0)
    Tot_El_Charge = Zero
    do iSym=1,nSym
      Tot_El_Charge = Tot_El_Charge-Two*real(nFro(iSym)+nIsh(iSym),kind=wp)
    end do
    Tot_El_Charge = Tot_El_Charge-real(nActEl,kind=wp)
    Tot_Charge = Tot_Nuc_Charge+Tot_El_Charge
    iCharge = int(Tot_Charge)
    call PrRF(.false.,NonEq,iCharge,2)
    if (DWSolv%DWZeta == -12345.0_wp) then
      write(u6,Fmt2//'A)') 'Weights of the reaction field are specified by RFROOT'
      write(u6,Fmt2//'(T45,10F6.3))') (W_SOLV(i),i=1,nRoots)
    else if (DWSolv%DWZeta < Zero) then
      call DWSol_fixed(i,j)
      if ((i == 0) .and. (j == 0)) then
        write(u6,Fmt2//'A)') 'Unrecognized negative DWZeta (DWSOl)'
        write(u6,Fmt2//'A,T51,A)') 'Dynamically weighted solvation is ','automatically turned off!'
      else
        write(u6,Fmt2//'A,T45,I2,X,I2)') 'Reaction field from states:',i,j
        if (max(i,j) > nRoots) then
          write(u6,Fmt2//'A)') 'The specified state is too high! Cannot proceed...'
          call Quit_OnUserError()
        end if
      end if
    else if (IPCMROOT <= 0) then
      write(u6,Fmt2//'A,T45,T15)') ' Reaction field from state:',' State-Averaged'
      if (DWSolv%DWZeta /= Zero) then
        write(u6,Fmt2//'A,T51,A)') 'Dynamically weighted solvation is ','automatically turned off!'
        DWSolv%DWZeta = Zero
      end if
    else
      write(u6,Fmt2//'A,T45,I2)') ' Reaction field from state:',IPCMROOT
      if (DWSolv%DWZeta > Zero) then
        write(u6,Fmt2//'A,ES10.3,A,I1,A)') 'Dynamically weighted solvation is used with DWSOlv = ',DWSolv%DWZeta,' (DWTYpe = ', &
                                           DWSolv%DWType,')'
        write(u6,Fmt2//'A,(T45,10F6.3))') 'Final weights for the reaction field:',(W_SOLV(i),i=1,nRoots)
      end if
    end if
  end if
  if (RFpert) then
    write(u6,*)
    write(u6,*)
    write(u6,Fmt2//'A)') 'Reaction field specifications:'
    write(u6,Fmt2//'A)') '------------------------------'
    write(u6,*)
    write(u6,'(6X,A)') 'The Reaction field has been added as a perturbation and has been determined in a previous calculation'
    write(u6,*)
  end if
  if ((KSDFT /= 'SCF') .and. (KSDFT /= 'PAM')) call Print_NQ_Info()
  call CollapseOutput(0,'CI expansion specifications:')

! End of long if-block A over IPRLEV
end if

EAV = ENER(lRootSplit,ITER)

CASDFT_Funct = Zero
if ((KSDFT /= 'SCF') .and. (KSDFT /= 'PAM')) then
  if ((nConf == 1) .and. (nActEl == nac)) then
    call Get_dScalar('CASDFT energy',CASDFT_Funct)
    EAV = EAV+CASDFT_Funct
  else
    call Get_dScalar('CASDFT energy',CASDFT_Funct)
    EAV = EAV+CASDFT_Funct-VIA_DFT-HALFQ
  end if
end if

if ((IPRLEV >= USUAL) .and. (.not. lOPTO)) then
  ! Start of long if-block B over IPRLEV
  write(u6,*)
  Line = ''
  write(Line(left-2:),'(A)') 'Final optimization specifications:'
  call CollapseOutput(1,Line)
  write(u6,Fmt2//'A)') '------------------------------'
  write(u6,*)

  write(u6,Fmt2//'A,T45,F20.8)') 'SplitCAS CI-energy',EAV
  write(u6,Fmt2//'A,T45,F20.8)') 'SplitCAS RAS-energy',ECAS
  write(u6,Fmt2//'A,T45,F20.8)') 'Super-CI energy',ESX
  write(u6,Fmt2//'A,T45,F20.8)') 'RASSCF energy change',DE
  write(u6,Fmt2//'A,T50,ES10.3)') 'Max change in MO coeff.',CMAX
  write(u6,Fmt2//'A,T50,ES10.3)') 'Max non-diagonal density matrix element',ROTMAX
  write(u6,Fmt2//'A,T50,ES10.3)') 'Max. BLB matrix element',CBLBM
  write(u6,Fmt2//'A,I4,A,I4,A,I4,A)') '(orbital pair',IBLBM,',',JBLBM,' in symmetry',ISYMBB,')'
  if (irlxRoot /= 0) write(u6,Fmt2//'A,T45,ES10.3)') 'Norm of electronic gradient',RLXGRD

  if (ISUPSM /= 0) then
    write(u6,Fmt2//'A)') 'Supersymmetry has been used to disable selected orbital rotations'
    iEnd = 0
    do iSym=1,nSym
      iStart = iEnd+1
      iEnd = iEnd+nBas(iSym)
      iTemp = 0
      do i=iStart,iEnd
        iTemp = iTemp+IXSYM(i)
      end do
      if (iTemp > 0) then
        write(u6,Fmt2//'A,I3)') 'Supersymmetry vector for symmetry species',iSym
        write(u6,Fmt2//'30I3)') (IXSYM(i),i=iStart,iEnd)
      end if
    end do
  end if

  if (allocated(CleanMask)) write(u6,Fmt2//'A)') 'The cleanup option has been used to set MO-coefficients explicitly to zero'
  call CollapseOutput(0,'Final optimization conditions:')

  ! End of long if-block B over IPRLEV
end if

call dcopy_(2*mxRoot,[Zero],0,Temp,1)
iRc1 = 0
iRc2 = 0
iOpt = ibset(0,sOpSiz)
iComp = 1
iSyLbl = 1
nMVInt = 0
nDCInt = 0
Label = 'MassVel'
call iRdOne(iRc1,iOpt,Label,iComp,iDum,iSyLbl)
if (iRc1 == 0) nMVInt = iDum(1)
Label = 'Darwin'
call iRdOne(iRc2,iOpt,Label,iComp,iDum,iSyLbl)
if (iRc2 == 0) nDCInt = iDum(1)
if ((nMVInt+nDCInt) /= 0) then
  IAD12 = IADR15(12)
  call mma_allocate(X1,NTOT1,Label='X1')
  call mma_allocate(X2,NTOT1,Label='X2')
  call mma_allocate(X3,NTOT,Label='X3')
  call mma_allocate(X4,NTOT2,Label='X4')

  kRoot = lRootSplit
  call DDAFILE(JOBIPH,2,X4,NTOT2,IAD12)
  call DDAFILE(JOBIPH,2,X3,NTOT,IAD12)
  call RelEne(Temp(1,kRoot),Temp(2,kRoot),nSym,nBas,X4,X3,X2,X1)

  call mma_deallocate(X4)
  call mma_deallocate(X3)
  call mma_deallocate(X2)
  call mma_deallocate(X1)
end if

if ((nMVInt+nDCInt) /= 0) then
  Emv = Temp(1,lRootSplit)
  Edc = Temp(2,lRootSplit)
  Erel = ENER(lRootSplit,ITER)+Emv+Edc
  EneTmp = Erel
else
  EneTmp = ENER(lRootSplit,ITER)+CASDFT_Funct-VIA_DFT-HALFQ
end if

if ((IPRLEV >= USUAL) .and. (.not. lOPTO)) then
  ! Start of long if-block C over IPRLEV
  !--------------------------------------------------------------------*
  !     Print final state energies (including relativistic correction) *
  !--------------------------------------------------------------------*
  write(u6,*)
  write(u6,*)
  write(u6,Fmt2//'A)') 'Final state energy(ies):'
  write(u6,Fmt2//'A)') '------------------------'
  write(u6,*)
  if ((nMVInt+nDCInt) /= 0) then
    write(u6,Fmt2//'A)') 'root     nonrelativistic        mass-velocity       Darwin-contact         relativistic'
    write(u6,Fmt2//'A)') '             energy                term                  term                 energy'
    i = lRootSplit
    Emv = Temp(1,i)
    Edc = Temp(2,i)
    Erel = ENER(I,ITER)+Emv+Edc
    write(u6,Fmt2//'I3,4(1X,F20.8))') I,ENER(I,ITER),Emv,Edc,Erel
  else
    write(u6,Fmt2//'A,I3,A,F20.8,A)') 'root number',lrootSplit,' E =',EneTmp,' a.u.'
  end if

else if ((IPRLEV >= TERSE) .and. (.not. lOPTO)) then
  write(u6,Fmt2//'A,I3,A,F20.8,A)') 'root number',lRootSplit,' E =',EneTmp,' a.u.'

  ! End of long if-block C over IPRLEV
end if

call Store_Energies(1,[EneTmp],1)

call Put_iScalar('NumGradRoot',irlxroot)
call Put_dScalar('Last energy',ECAS)
call Put_dScalar('Average energy',EAV)

iTol = Cho_X_GetTol(8)
call Add_Info('E_RASSCF',ENER(1,ITER),1,iTol)

!----------------------------------------------------------------
! New JOBIPH layout: Also write hamiltonian matrix at IADR15(17):
!----------------------------------------------------------------
IAD15 = IADR15(17)
call DDAFILE(JOBIPH,1,ENER(1,ITER),1,IAD15)

!----------------------------------------------------------------------*
!     Print the multipole analysis of the solvation energy             *
!----------------------------------------------------------------------*
call RFMltp()

!----------------------------------------------------------------------*
!     Print the final orbitals                                         *
!----------------------------------------------------------------------*
FullMlk = (OutFmt1 /= 'NOTHING ')
if ((IPRLEV >= USUAL) .and. (.not. lOPTO)) then
  ! Start of if-block D over IPRLEV
  if (OutFmt1 /= 'NOTHING ') then
    if (IPT2 == 0) then
      call PRIMO_RASSCF('Pseudonatural active orbitals and approximate occupation numbers',FDIAG,OCCN,CMO)
    else
      call PRIMO_RASSCF('All orbitals are eigenfunctions of the PT2 Fock matrix',FDIAG,OCCN,CMO)
    end if
  end if
  ! End of if-block D over IPRLEV
end if

!----------------------------------------------------------------------*
!     Also put on RUNFILE (in the future...):                          *
!----------------------------------------------------------------------*
call Put_dArray('RASSCF orbitals',CMO,NTOT2)

!----------------------------------------------------------------------*
!     compute properties and Mulliken's orbital populations            *
!----------------------------------------------------------------------*
IAD12 = IADR15(12)
IAD03 = IADR15(3)
IAD14 = IADR15(14)
!BOR0511
! Save original orbitals for the spin density matrices
call mma_allocate(cmon,nTot2,Label='CMON')
call dcopy_(ntot2,cmo,1,cmon,1)
!BOR0511
FullMlk = (OutFmt1 /= 'NOTHING ')

!***********************************************************************
!     save the 1st order density for gradients                         *
!***********************************************************************
call mma_allocate(DSave,nTot1,Label='DSave')
call Get_dArray_chk('D1ao',DSave,NTOT1)

! Read natural orbitals
if (NAC > 0) then
  call DDAFILE(JOBIPH,2,CMO,NTOT2,IAD12)
  call DDAFILE(JOBIPH,2,OCCN,NTOT,IAD12)
end if

! Put the density matrix of this state on the runfile for
!  LoProp utility
call mma_allocate(D,nTot1,Label='D')
D(:) = Zero
call DONE_RASSCF(CMO,OCCN,D)
call Put_dArray('D1ao',D,NTOT1)
call mma_deallocate(D)

if (IPRLEV >= USUAL) then
  ! Start of if-block E over IPRLEV

  ! Compute Mulliken's population analysis
  write(u6,'(/6X,A,I3)') 'Mulliken population Analysis for root number:',lRootSplit
  write(u6,'(6X,A)') '-----------------------------------------------'
  write(u6,*)
  lSave = lRootSplit == iRlxRoot
  call CHARGE(nsym,nbas,BName,CMO,OCCN,SMAT,2,FullMlk,lSave)
  write(u6,*)

  ! Compute properties
  write(u6,'(/6X,A,I3)') 'Expectation values of various properties for root number:',lRootSplit
  write(u6,'(6X,A)') '-----------------------------------------------------------'
  write(u6,*)
end if
! End of if-block E over IPRLEV

! Write out info on a temporary vector file. Prpt
! will need to read this.
Note = 'Temporary orbital file used by prpt3.'
LuTmp = 50
LuTmp = IsFreeUnit(LuTmp)
call WrVec('TMPORB',LuTmp,'CO',nSym,nBas,nBas,CMO,OCCN,Dum,iDum,Note)
call PRPT()

! Compute spin orbitals and spin population
! (Note: this section overwrites the pseudo natural orbitals
!        with the spin orbitals).

call mma_allocate(X6,NACPAR,Label='X6')
call DDAFILE(JOBIPH,0,X6,NACPAR,IAD03)
call DDAFILE(JOBIPH,2,X6,NACPAR,IAD03)
call DDAFILE(JOBIPH,0,X6,NACPR2,IAD03)
call DDAFILE(JOBIPH,0,X6,NACPR2,IAD03)
call DBLOCK(X6)

if (IPRLEV >= VERBOSE) then
  ! Start of long if-block F over IPRLEV
  if (ISPDEN == 1) then
    ! Print spin density matrix
    write(u6,'(/6X,A,I3)') 'Spin density matrix for root number:',lRootSplit
    write(u6,'(6X,A)') '--------------------------------------'
    write(u6,*)
    IND = 1
    IDIMV = 0
    IDIMO = 0
    IDIMN = 0
    do ISYM=1,NSYM
      NAO = NASH(ISYM)
      if (NAO == 0) GO TO 50
      NO = NBAS(ISYM)
      IDIMV = IDIMV+NAO*NAO
      IDIMO = IDIMO+NAO
      IDIMN = IDIMN+NO*NAO
      if (iprlev >= VERBOSE) then
        write(u6,'(/6X,A,I2)') 'symmetry species',ISYM
        write(u6,*)
        call TRIPRT(' ',' ',X6(IND),NASH(ISYM))
      end if
      IND = IND+NASH(ISYM)*(NASH(ISYM)+1)/2
50    continue
    end do
  end if

! End of long if-block F over IPRLEV
end if

! Compute spin orbitals and spin population
call DCOPY_(NTOT,[Zero],0,OCCN,1)
!SVC-11-01-2007 store original cmon in cmoso, which gets changed
call mma_allocate(CMOSO,NTOT2,Label='CMOSO')
call DCOPY_(NTOT2,CMON,1,CMOSO,1)

call SPINORB(X6,CMOSO,OCCN)
call mma_deallocate(X6)
call DDAFILE(JOBIPH,1,CMOSO,NTOT2,IAD14)
call DDAFILE(JOBIPH,1,OCCN,NTOT,IAD14)

if (IPRLEV >= USUAL) then
  ! Start of long if-block G over IPRLEV
  if (ISPDEN == 1) then
    write(u6,'(/6X,A,I3)') 'Mulliken spin population Analysis for root number:',lRootSplit
    write(u6,'(6X,A)') '---------------------------------------------------'
    write(u6,*)
    call CHARGE(nsym,nbas,BName,cmoso,OCCN,SMAT,3,FullMlk,.false.)
    write(u6,*)
  end if

  !BOR0511
  ! Do LoProp Charge analysis
  if (get_BasisType('ANO')) then
    write(u6,'(/6X,A,I3)') 'LoProp population Analysis for root number:',lRootSplit
    write(u6,'(6X,A)') '-----------------------------------------------'
    write(u6,*)
    write(u6,*)
    Line = ''
    write(Line(left-2:),'(A)') 'LoProp Analysis:'
    call CollapseOutput(1,Line)
    write(u6,Fmt2//'A)') '----------------'
    write(u6,*)
    write(u6,*)
    ! Ugly trick: use IRC to pass some info to LoProp ...
    iRC = min((abs(lRootSplit-iRlxRoot)),1)
    call LoProp(iRC)
    write(u6,*)
    write(u6,'(6X,A,I2)') 'Natural Bond Order Analysis for root number:',lRootSplit
    call Nat_Bond_Order(nSym,nBas,BName,2)
    call CollapseOutput(0,'LoProp Analysis:')
    write(u6,*)
  end if

  ! End of long if-block G over IPRLEV
end if
call mma_deallocate(CMOSO)

! ESPF analysis
call DecideOnESPF(Do_ESPF)
if (Do_ESPF) call espf_analysis(.false.)

! Restore the correct 1st order density for gradient calculations.

call Put_dArray('D1ao',DSave,NTOT1)
call mma_deallocate(DSave)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(CMON)
!----------------------------------------------------------------------*

end subroutine OutCtlSplit
