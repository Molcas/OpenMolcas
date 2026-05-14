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
! Copyright (C) 1990, Markus P. Fuelscher                              *
!               2013, Giovanni Li Manni                                *
!               2016, Andrew M. Sand                                   *
!               2024, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine compute_mcpdft_energy(cmo,e_mcscf,e_states)

use mcpdft_input, only: mcpdft_options
use KSDFT_Info, only: do_pdftpot
use mspdftgrad, only: D1aoMS, D1SaoMS, DIDA, P2MOT
use libxc_parameters, only: FuncExtParams
use rasscf_global, only: IADR15, lRoots, NAC, NACPAR, NACPR2, PotNuc
use general_data, only: ispin, jobiph, jobold, nash, nbas, nish, norb, nsym, ntot1, ntot2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: iwp, wp

implicit none
real(kind=wp), intent(in) :: cmo(ntot2), e_mcscf(lroots)
real(kind=wp), intent(out) :: e_states(lroots)
integer(kind=iwp) :: charge, dmDisk, IAD19, IADR19(30), iJOB, iSA, isym, jroot, niaia, NQ
real(kind=wp) :: e_ot, e_state, e_wfn
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: casdm1(:), casdm1s(:), Coul(:), dm1_cas(:), dm1_core(:), dm1s(:), folded_dm1(:), folded_dm1_cas(:), &
                              folded_dm1s(:), hcore(:), P2D(:), P2t(:)
integer(kind=iwp), external :: get_charge
real(kind=wp), external :: energy_mcwfn

! Molecular charge
charge = get_charge()

! Load h_mu,nu (1e in AO basis)
call mma_allocate(hcore,nTot1,label='hcore')
call get_hcore(hcore)

! Here we calculate the D1 Inactive matrix (in AO basis).
call mma_allocate(dm1_core,ntot2,label='dm1_core')
call Get_D1I_RASSCF(CMO,dm1_core)

if (mcpdft_options%grad .and. mcpdft_options%mspdft) call Fold2(nSym,nBas,dm1_core,DIDA(:,lroots+1))

iJOB = 0
IAD19 = 0
call f_Inquire('JOBOLD',Found)
if (.not. found) then
  call f_Inquire('JOBIPH',Found)
  if (Found) JOBOLD = JOBIPH
end if
if (Found) iJOB = 1
if (iJOB == 1) then
  if (JOBOLD <= 0) then
    JOBOLD = 20
    call DaName(JOBOLD,'JOBOLD')
  end if
end if
IADR19(:) = 0
call IDaFile(JOBOLD,2,IADR19,15,IAD19)
IADR15 = IADR19
dmDisk = IADR19(3)

call mma_allocate(casdm1,NACPAR,Label='casdm1')
call mma_allocate(dm1_cas,NTOT2,Label='dm1_cas')
call mma_allocate(casdm1s,NACPAR,Label='casdm1s')
call mma_allocate(dm1s,NTOT2,Label='dm1s')

call mma_allocate(folded_dm1,nTot1,Label='folded_dm1')
call mma_allocate(folded_dm1_cas,nTot1,Label='folded_dm1_cas')
call mma_allocate(folded_dm1s,nTot1,Label='folded_dm1s')
call mma_allocate(P2d,NACPR2,Label='P2D')

call mma_allocate(coul,ntot1,Label='coul')

! This iSA is used to control gradient calculations.  Analytic gradients
! (in ALASKA) will only run if iSA=1, and iSA will only be set to one if
! the on-top potentials are computed as part of this calculation.
iSA = 99
call Put_iScalar('SA ready',iSA)

! We loop over the number of states.  For each state, we read the density
! matrices from the JOBIPH file.
do jroot=1,lroots
  ! Read in the density matrices for <jroot>.
  casdm1(:) = Zero
  dm1_cas(:) = Zero
  casdm1s(:) = Zero
  dm1s(:) = Zero
  folded_dm1(:) = Zero
  folded_dm1s(:) = Zero
  folded_dm1_cas(:) = Zero
  P2D(:) = Zero

  ! Get the D1 Active matrix for this state.  These should probably be
  ! most easily read from the previous JOBIPH file.  Then, convert D1A from
  ! the MO to the AO basis.

  call DDaFile(JOBOLD,2,casdm1,NACPAR,dmDisk)
  call DDaFile(JOBOLD,2,casdm1s,NACPAR,dmDisk)
  call DDaFile(JOBOLD,2,P2d,NACPR2,dmDisk)
  ! This dummy read is to cycle the dmDisk so next
  ! time we encounter this block, we load the next
  ! states densities
  call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)

  ! This should be in energy_ot() but it gets
  ! modified in dblock. But also, why do we need
  ! this in the DFT section???
  call Put_dArray('D1mo',casdm1,NACPAR)

  ! Generate total density
  if (NASH(1) /= NAC) then
    call dblock(casdm1)
    call dblock(casdm1s)
  end if
  call Get_D1A_RASSCF(CMO,casdm1,dm1_cas)
  call Get_D1A_RASSCF(CMO,casdm1s,dm1s)

  call Fold(nSym,nBas,dm1_core,folded_dm1)
  call Fold(nSym,nBas,dm1_cas,folded_dm1_cas)
  folded_dm1(:) = folded_dm1(:)+folded_dm1_cas

  call Fold(nSym,nBas,dm1s,folded_dm1s)

  ! Maybe I can write all of these matrices to file, then modify stuff in
  ! the nq code to read in the needed density.  In other words, I need to
  ! replace the next call with something that supports multistates.

  ! save things for MSPDFT gradients
  if (mcpdft_options%grad .and. mcpdft_options%mspdft) then
    call fold2(nsym,nbas,dm1_cas,DIDA(:,jroot))
    d1aoms(:,jroot) = folded_dm1(:)
    call P2_contraction(casdm1,P2MOt(:,jroot))
    if (ispin /= 1) d1saoms(:,jroot) = folded_dm1s(:)
  end if

  do_pdftpot = (mcpdft_options%grad .and. (mcpdft_options%mspdft .or. (jroot == mcpdft_options%rlxroot)))
  e_ot = mcpdft_options%otfnal%energy_ot(folded_dm1,folded_dm1s,P2d,charge)

  call get_coulomb(cmo,dm1_core,dm1_cas,coul)

  e_wfn = energy_mcwfn(folded_dm1,hcore,coul,PotNuc,ntot1)

  e_state = e_wfn+e_ot

  if (mcpdft_options%otfnal%is_hybrid()) &
    e_state = mcpdft_options%otfnal%lambda*e_mcscf(jRoot)+(One-mcpdft_options%otfnal%lambda)*e_state

  call Print_MCPDFT_2(PotNuc,e_wfn,e_ot,jroot,e_mcscf(jroot))

  ! JB replacing e_mcscf with MC-PDFT energy for MS-PDFT use
  e_states(jroot) = e_state

  ! At this point, the energy calculation is done.  Now I need to build the
  ! fock matrix if this root corresponds to the relaxation root.
  if (mcpdft_options%grad) then
    ! Determine size of Q matrix
    NQ = 0
    NIAIA = 0
    do ISYM=1,NSYM
      NQ = max(NQ,NASH(ISYM)*NORB(ISYM))
      NIAIA = NIAIA+(NASH(ISYM)+NISH(ISYM))**2
    end do
    if (NQ < NIAIA) NQ = NIAIA

    if ((.not. mcpdft_options%mspdft) .and. (jroot == mcpdft_options%rlxroot)) &
      call savefock_pdft(cmo,hcore(:)+coul(:),casdm1,nq,p2d)
    if (mcpdft_options%mspdft) call savefock_mspdft(CMO,hcore(:)+coul(:),casdm1,NQ,p2d,jroot)
  end if
end do !loop over roots

if (mcpdft_options%grad) then
  dmDisk = IADR19(3)
  do jroot=1,mcpdft_options%rlxroot-1
    call DDaFile(JOBOLD,0,casdm1,NACPAR,dmDisk)
    call DDaFile(JOBOLD,0,casdm1s,NACPAR,dmDisk)
    call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
    call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
  end do
  call DDaFile(JOBOLD,2,casdm1,NACPAR,dmDisk)
  ! Andrew added this line to fix heh2plus
  call DDaFile(JOBOLD,2,casdm1s,NACPAR,dmDisk)
  call Put_dArray('D1mo',casdm1,NACPAR)
  call DDaFile(JOBOLD,2,P2d,NACPR2,dmDisk)
  call Put_dArray('P2mo',P2d,NACPR2)
  ! from before???
  call mma_allocate(P2t,NACPR2,Label='P2t')
  P2t(:) = Zero
  call P2_contraction(casdm1,P2t)
  call Put_dArray('P2MOt',P2t,NACPR2)
  call mma_deallocate(P2t)

  if (NASH(1) /= NAC) call DBLOCK(casdm1)
  call Get_D1A_RASSCF(CMO,casdm1,dm1_cas)

  call Fold(nSym,nBas,dm1_core,folded_dm1)
  call Fold(nSym,nBas,dm1_cas,folded_dm1_cas)
  folded_dm1(:) = folded_dm1(:)+folded_dm1_cas(:)
  call Put_dArray('D1ao',folded_dm1,nTot1)

  !  Generate spin-density
  dm1s(:) = Zero
  if (NASH(1) /= NAC) call DBLOCK(casdm1s)
  call Get_D1A_RASSCF(CMO,casdm1s,dm1s)
  call Fold(nSym,nBas,dm1s,folded_dm1s)
  call Put_dArray('D1Sao',folded_dm1s,nTot1)

  call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
end if

call mma_deallocate(folded_dm1_cas)
call mma_deallocate(folded_dm1)
call mma_deallocate(folded_dm1s)
call mma_deallocate(casdm1)
call mma_deallocate(dm1_cas)
call mma_deallocate(casdm1s)
call mma_deallocate(dm1s)
call mma_deallocate(hcore)
call mma_deallocate(coul)
call mma_deallocate(P2D)
call mma_deallocate(dm1_core)

! This deallocation SHOULD NOT BE HERE
! but it has to do with how the FuncExtParams is loaded...
if (allocated(FuncExtParams)) call mma_deallocate(FuncExtParams)

end subroutine compute_mcpdft_energy
