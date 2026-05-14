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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!               2017, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PrFin(OneHam,Ovlp,Dens,TwoHam,nDT,EOrb,OccNo,nEO,CMO,nCMO,note,iCase,MssVlc,Darwin)
!***********************************************************************
!                                                                      *
!     purpose: Final printout                                          *
!                                                                      *
!     input:                                                           *
!       Dens    : the last total density                               *
!       TwoHam  : the last two-electron ham.                           *
!       OneHam  : one-electron hamiltonian of length nOOK              *
!       Ovlp    : overlap in AO basis of length nOOK                   *
!       EOrb    : orbital energies of length nEO                       *
!       OccNo   : occupation numbers of length nEO                     *
!       CMO     : molecular orbital coefficients of length nCMO        *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use SpinAV, only: Do_SpinAV
use InfSCF, only: BName, EneV, ExFac, iCoCo, InVec, iPrForm, iPrint, iPrOrb, jPrint, kIvo, KSDFT, lRel, nBas, nBB, nBT, nD, nIter, &
                  nIterP, nnB, NoProp, nOrb, nSYm, PotNuc, ThrEne, Tot_Charge
use rctfld_module, only: lRF
use NDDO, only: oneel_NDDO
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDT, nEO, nCMO, iCase
real(kind=wp), intent(in) :: OneHam(nDT), Ovlp(nDT), Dens(nDT), TwoHam(nDT), EOrb(nEO), OccNo(nEO), CMO(nCMO), MssVlc(nDT), &
                             Darwin(nDT)
character(len=80), intent(out) :: Note
integer(kind=iwp) :: iBs, iCharge, iCMO, ij, i_Or, iPL, iRC, iSpin, iSym, iv, iVec, jCase
real(kind=wp) :: EHomo, ELumo, ERelDC, ERelMV
logical(kind=iwp) :: DeBug, Dff, Do_DFT, Do_ESPF, First, FullMlk, get_BasisType, NonEq, PrEne, PrOcc
character(len=60) :: Frmt
character(len=30) :: AlphaLabel
real(kind=wp), allocatable :: RFfld(:), Scr2(:), Scr3(:)
real(kind=wp), parameter :: ThrOcc = -99999.0_wp
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: EFP_On, Reduce_Prt

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

jPrint = iPrint
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
if (iPL <= 1) jPrint = 1

!----------------------------------------------------------------------*

#ifdef _DEBUGPRINT_
Debug = .true.
#else
Debug = .false.
if (jPrint >= 4) Debug = .true.
#endif

Frmt = '(6X,A,T50,F17.10)'
AlphaLabel = ' '
if (nD == 2) then
  if (iCase == 0) AlphaLabel = ' (alpha) '
  if (iCase == 1) AlphaLabel = ' (beta)  '
end if

if (Do_SpinAV) AlphaLabel = Alphalabel(1:9)//'and (spin-averaged)'

! Compute relativistic corrections
if (lRel) then
  call RelEny(ERelMV,ERelDC,Dens,MssVlc,Darwin,nBT)
  if (jPrint >= 2) then
    write(u6,*)
    write(u6,'(6X,A)') '1st order relativistic corrections'
    write(u6,Frmt) 'Total energy',EneV+ERelMV+ERelDC
    write(u6,Frmt) 'Mass-velocity correction',ERelMV
    write(u6,Frmt) '1-el Darwin correction',ERelDC
    write(u6,Frmt) 'Sum of relatvity corrections',ERelDC+ERelMV
    write(u6,*)
  end if
end if

! Print numerical quadrature information
iSpin = 1
if (nD == 2) iSpin = 2
if ((KSDFT /= 'SCF') .and. (iCase == 0)) call Print_NQ_Info()

! Write out last density matrix to output
if (DeBug) then
  write(u6,'(6x,A)') 'Last density matrix (interpolated) in AO basis'
  ij = 1
  do iSym=1,nSym
    write(u6,*) ' symmetry',iSym
    call TriPrt(' ',' ',Dens(ij),nBas(iSym))
    ij = ij+nTri_Elem(nBas(iSym))
  end do
  write(u6,*)
  write(u6,'(6x,A)') 'Last 2-el. Hamiltonian (interpolated) in AO basis'
  ij = 1
  do iSym=1,nSym
    write(u6,*) ' symmetry',iSym
    call TriPrt(' ',' ',TwoHam(ij),nBas(iSym))
    ij = ij+nTri_Elem(nBas(iSym))
  end do
  write(u6,*)
  write(u6,'(6x,A)') 'Last 1-el. Hamiltonian (interpolated) in AO basis'
  ij = 1
  do iSym=1,nSym
    write(u6,*) ' symmetry',iSym
    call TriPrt(' ',' ',OneHam(ij),nBas(iSym))
    ij = ij+nTri_Elem(nBas(iSym))
  end do
  write(u6,*)
end if

! Process the external potential with the total electronic density

!nf
call DecideOnESPF(Do_ESPF)
!nf if ((lRF .or. (KSDFT /= 'SCF')) .and. (.not. oneel_NDDO) .and.
if ((Do_ESPF .or. lRF .or. (KSDFT /= 'SCF') .or. EFP_On()) .and. (.not. oneel_NDDO) .and. (iCase == 0)) then
  iCharge = int(Tot_Charge)
  NonEq = .false.
  !call Get_PotNuc(PotNuc)
  !call Get_dScalar('PotNuc',PotNuc)
  call Peek_dScalar('PotNuc',PotNuc)
  call mma_allocate(RFfld,nBT,Label='RFfld')
  RFfld(:) = Zero
  First = .true.
  Dff = .false.
  Do_DFT = .false. ! We do not need to redo the DFT!
  call DrvXV(RFfld,RFfld,Dens,PotNuc,nBT,First,Dff,NonEq,lRF,KSDFT,ExFac,iCharge,iSpin,'SCF ',Do_DFT)
  call mma_deallocate(RFfld)

  ! Print multipole analysis of the reaction field contributions

  call RFmltp()

end if

! Print orbitals (the case InVec=3 and nIter=0 is set up in RdInp)
Frmt = '(6X,A)'
if (iPrOrb >= 1) then
  FullMlk = .true.
  PrEne = .true.
  PrOcc = .true.
  if (iPrOrb == 1) then
    EHomo = -99999.0_wp
    ELumo = 99999.0_wp
    ij = 1
    do iSym=1,nSym
      do iv=1,nOrb(iSym)
        if (OccNo(ij+iv-1) > 0.001_wp) then
          EHomo = max(EHomo,EOrb(ij+iv-1))
        else
          ELumo = min(ELumo,EOrb(ij+iv-1))
        end if
      end do
      ij = ij+nOrb(iSym)
    end do
    ThrEne = ELumo+Half
    if (jPrint >= 2) then
      write(u6,*)
      write(u6,Frmt) 'All orbitals with orbital energies smaller than  E(LUMO)+0.5 are printed'
    end if
  else
    if (jPrint >= 2) then
      write(u6,*)
      write(u6,'(6X,A,ES11.4,A)') 'All orbitals with orbital energies smaller than',ThrEne,' are printed'
    end if
  end if
  if (KSDFT == 'SCF') then
    if (nD == 1) then
      Note = 'SCF orbitals'//AlphaLabel
      if (kIvo /= 0) Note = 'SCF orbitals + IVO'
      if (iCoCo /= 0) Note = 'SCF orbitals + arbitrary occupations'
    else
      Note = 'UHF orbitals'//AlphaLabel
      if (kIvo /= 0) Note = 'UHF orbitals + IVO'
      if (iCoCo /= 0) Note = 'UHF orbitals + arbitrary occupations'
    end if
  else
    if (nD == 1) then
      Note = 'RKS-DFT orbitals'//AlphaLabel
      if (kIvo /= 0) Note = 'RKS-DFT orbitals + IVO'
      if (iCoCo /= 0) Note = 'RKS-DFT orbitals + arbitrary occupations'
    else
      Note = 'UKS-DFT orbitals'//AlphaLabel
      if (kIvo /= 0) Note = 'UKS-DFT orbitals + IVO'
      if (iCoCo /= 0) Note = 'UKS-DFT orbitals + arbitrary occupations'
    end if
  end if
  if (jPrint >= 2) call PriMO(Note,PrOcc,PrEne,ThrOcc,ThrEne,nSym,nBas,nOrb,BName,EOrb,OccNo,CMO,iPrForm)
else
  if (jPrint >= 2) write(u6,Frmt) 'No orbitals printed'
  FullMlk = .false.
end if

if ((InVec /= 3) .or. (nIter(nIterP) > 0)) then
  call mma_allocate(Scr2,nBB,Label='Scr2')
  call mma_allocate(Scr3,nnB,Label='Scr3')

  ! Prepare CMO in symmetry blocks nBas x nBas
  iVec = 0
  iCMO = 0
  do iSym=1,nSym
    Scr2(iVec+1:iVec+nBas(iSym)*nOrb(iSym)) = CMO(iCMO+1:iCMO+nBas(iSym)*nOrb(iSym))
    iVec = iVec+nBas(iSym)*nOrb(iSym)
    Scr2(iVec+1:iVec+nBas(iSym)*(nBas(iSym)-nOrb(iSym))) = Zero
    iVec = iVec+nBas(iSym)*(nBas(iSym)-nOrb(iSym))
    iCMO = iCMO+nOrb(iSym)*nBas(iSym)
  end do

  ! Prepare occupation numbers

  i_Or = 0
  iBs = 0
  do iSym=1,nSym
    Scr3(iBs+1:iBs+nOrb(iSym)) = OccNo(i_Or+1:i_Or+nOrb(iSym))
    Scr3(iBs+nOrb(iSym)+1:iBs+nBas(iSym)) = Zero
    i_Or = i_Or+nOrb(iSym)
    iBs = iBs+nBas(iSym)
  end do

  if (.not. NoProp) then

    ! Population analysis

    jCase = iCase
    if (nD == 1) jCase = 2
    call Charge(nSym,nBas,BName,Scr2,Scr3,Ovlp,jCase,FullMlk,.true.)

    if (get_BasisType('ANO')) then
      iRc = 0
      call LoProp(iRc)

      ! NBO analysis

      call Nat_Bond_Order(nSym,nBas,BName,jCase)
    end if
  end if

  ! ESPF analysis

  if (Do_ESPF) call espf_analysis(.true.)

  call mma_deallocate(Scr3)
  call mma_deallocate(Scr2)

end if

end subroutine PrFin
