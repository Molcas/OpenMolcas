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

use SpinAV, only: Do_SpinAV
use InfSCF, only: EneV, ExFac, iCoCo, InVec, iPrForm, iPrint, iPrOrb, nD, jPrint, kIvo, KSDFT, lRel, Name, nBB, nBT, nIterP, nnB, &
                  NoProp, nSYm, PotNuc, ThrEne, ThrOcc, Tot_Charge, nBas, nOrb, nIter
use Constants, only: Zero
use stdalloc, only: mma_allocate, mma_deallocate
use rctfld_module, only: lRF
use NDDO, only: oneel_NDDO

implicit none
integer nDT, nEO, nCMO
real*8 Dens(nDT), TwoHam(nDT), OneHam(nDT), Ovlp(nDT), EOrb(nEO), OccNo(nEO), CMO(nCMO), MssVlc(nDT), Darwin(nDT)
character(len=80) Note
external EFP_ON
! Define local variables
integer iCase, i, iBs, iCharge, iCMO, ij, iOr, iPL, iRC, iSpin, iSym, iv, iVec, j, jCase
integer, external :: iPrintLevel
real*8 EHomo, ELumo, ERelMV, ERelDC
character(len=60) Fmt
logical PrEne, PrOcc, get_BasisType
logical DeBug, First, NonEq, Dff, Do_DFT, FullMlk, Reduce_Prt
external Reduce_Prt
real*8, dimension(:), allocatable :: RFfld, Scr2, Scr3
!nf
logical Do_ESPF, EFP_On
!nf
character AlphaLabel*30
#include "SysDef.fh"

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

Fmt = '(6X,A,T50,F17.10)'
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
    write(6,*)
    write(6,'(6X,A)') '1st order relativistic corrections'
    write(6,Fmt) 'Total energy',EneV+ERelMV+ERelDC
    write(6,Fmt) 'Mass-velocity correction',ERelMV
    write(6,Fmt) '1-el Darwin correction',ERelDC
    write(6,Fmt) 'Sum of relatvity corrections',ERelDC+ERelMV
    write(6,*)
  end if
end if

! Print numerical quadrature information
iSpin = 1
if (nD == 2) iSpin = 2
if ((KSDFT /= 'SCF') .and. (iCase == 0)) call Print_NQ_Info()

! Write out last density matrix to output
if (DeBug) then
  write(6,'(6x,A)') 'Last density matrix (interpolated) in AO basis'
  ij = 1
  do iSym=1,nSym
    write(6,*) ' symmetry',iSym
    call TriPrt(' ',' ',Dens(ij),nBas(iSym))
    ij = ij+nBas(iSym)*(nBas(iSym)+1)/2
  end do
  write(6,*)
  write(6,'(6x,A)') 'Last 2-el. Hamiltonian (interpolated) in AO basis'
  ij = 1
  do iSym=1,nSym
    write(6,*) ' symmetry',iSym
    call TriPrt(' ',' ',TwoHam(ij),nBas(iSym))
    ij = ij+nBas(iSym)*(nBas(iSym)+1)/2
  end do
  write(6,*)
  write(6,'(6x,A)') 'Last 1-el. Hamiltonian (interpolated) in AO basis'
  ij = 1
  do iSym=1,nSym
    write(6,*) ' symmetry',iSym
    call TriPrt(' ',' ',OneHam(ij),nBas(iSym))
    ij = ij+nBas(iSym)*(nBas(iSym)+1)/2
  end do
  write(6,*)
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
  call dcopy_(nBT,[Zero],0,RFfld,1)
  First = .true.
  Dff = .false.
  Do_DFT = .false. ! We do not need to redo the DFT!
  call DrvXV(RFfld,RFfld,Dens,PotNuc,nBT,First,Dff,NonEq,lRF,KSDFT,ExFac,iCharge,iSpin,'SCF ',Do_DFT)
  call mma_deallocate(RFfld)

  ! Print multipole analysis of the reaction field contributions

  call RFmltp()

end if

! Print orbitals (the case InVec=3 and nIter=0 is set up in RdInp)
Fmt = '(6X,A)'
if (iPrOrb >= 1) then
  FullMlk = .true.
  PrEne = .true.
  PrOcc = .true.
  if (iPrOrb == 1) then
    EHomo = -99999.0d+00
    ELumo = 99999.0d+00
    ij = 1
    do iSym=1,nSym
      do iv=1,nOrb(iSym)
        if (OccNo(ij+iv-1) > 0.001) then
          EHomo = max(EHomo,EOrb(ij+iv-1))
        else
          ELumo = min(ELumo,EOrb(ij+iv-1))
        end if
      end do
      ij = ij+nOrb(iSym)
    end do
    ThrEne = ELumo+0.5
    if (jPrint >= 2) then
      write(6,*)
      write(6,Fmt) 'All orbitals with orbital energies smaller than  E(LUMO)+0.5 are printed'
    end if
  else
    if (jPrint >= 2) then
      write(6,*)
      write(6,'(6X,A,ES11.4,A)') 'All orbitals with orbital energies smaller than',ThrEne,' are printed'
    end if
  end if
  ThrOcc = -99999.0d+00
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
  if (jPrint >= 2) call PriMO(Note,PrOcc,PrEne,ThrOcc,ThrEne,nSym,nBas,nOrb,Name,EOrb,OccNo,CMO,iPrForm)
else
  if (jPrint >= 2) write(6,Fmt) 'No orbitals printed'
  FullMlk = .false.
end if

if ((InVec /= 3) .or. (nIter(nIterP) > 0)) then
  call mma_allocate(Scr2,nBB,Label='Scr2')
  call mma_allocate(Scr3,nnB,Label='Scr3')

  ! Prepare CMO in symmetry blocks nBas x nBas
  iVec = 0
  iCMO = 0
  do iSym=1,nSym
    do i=1,nBas(iSym)*nOrb(iSym)
      Scr2(iVec+i) = CMO(iCMO+i)
    end do
    iVec = iVec+nBas(iSym)*nOrb(iSym)
    do i=1,nBas(iSym)*(nBas(iSym)-nOrb(iSym))
      Scr2(iVec+i) = Zero
    end do
    iVec = iVec+nBas(iSym)*(nBas(iSym)-nOrb(iSym))
    iCMO = iCMO+nOrb(iSym)*nBas(iSym)
  end do

  ! Prepare occupation numbers

  iOr = 0
  iBs = 0
  do iSym=1,nSym
    do j=1,nOrb(iSym)
      Scr3(iBs+j) = OccNo(iOr+j)
    end do
    do j=nOrb(iSym)+1,nBas(iSym)
      Scr3(iBs+j) = Zero
    end do
    iOr = iOr+nOrb(iSym)
    iBs = iBs+nBas(iSym)
  end do

  if (.not. NoProp) then

    ! Population analysis

    jCase = iCase
    if (nD == 1) jCase = 2
    call Charge(nSym,nBas,Name,Scr2,Scr3,Ovlp,jCase,FullMlk,.true.)

    if (get_BasisType('ANO')) then
      iRc = 0
      call LoProp(iRc)

      ! NBO analysis

      call Nat_Bond_Order(nSym,nBas,Name,jCase)
    end if
  end if

  ! ESPF analysis

  if (Do_ESPF) call espf_analysis(.true.)

  call mma_deallocate(Scr3)
  call mma_deallocate(Scr2)

end if

end subroutine PrFin
