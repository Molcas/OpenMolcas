!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
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
  use definitions,only:iwp,wp
  use constants,only:zero,one
  use mcpdft_input,only:mcpdft_options
  Use KSDFT_Info,only:do_pdftpot
  use mspdftgrad,only:P2MOT,D1aoMS,DIDA,D1SaoMS
  use stdalloc,only:mma_allocate,mma_deallocate
  use libxc_parameters,only:FuncExtParams
  use rasscf_global,only:IADR15,lRoots,NAC,NACPAR,NACPR2,PotNuc
  use general_data,only:nash,norb,nsym,ntot2,ntot1,jobiph,ispin,jobold,nbas,nish
  implicit none

  real(kind=wp),intent(in) :: e_mcscf(*), cmo(*)
  real(kind=wp),intent(out) :: e_states(*)

  real(kind=wp),allocatable :: folded_dm1(:),folded_dm1_cas(:),folded_dm1s(:), &
                               dm1_core(:),casdm1(:),dm1_cas(:),dm1s(:),casdm1s(:),P2D(:),P2t(:),Coul(:),hcore(:)
  logical(kind=iwp) :: Found
  integer(kind=iwp),external :: get_charge
  integer(kind=iwp) :: IAD19,iJOB,dmDisk,IADR19(1:30),jroot,NQ,isym,charge,iSA,niaia
  real(kind=wp),external :: energy_mcwfn
  real(kind=wp) :: e_state,e_ot,e_wfn


  ! Molecular charge
  charge = get_charge()

  ! Load h_mu,nu (1e in AO basis)
  call mma_allocate(hcore,nTot1,label='hcore')
  call get_hcore(hcore)

  ! Here we calculate the D1 Inactive matrix (in AO basis).
  call mma_allocate(dm1_core,ntot2,label="dm1_core")
  Call Get_D1I_RASSCF(CMO,dm1_core)

  if(mcpdft_options%grad .and. mcpdft_options%mspdft) then
    Call Fold2(nSym,nBas,dm1_core,DIDA(:,lroots+1))
  endif

  iJOB = 0
  IAD19 = 0
  Call f_Inquire('JOBOLD',Found)
  If(.not. found) then
    Call f_Inquire('JOBIPH',Found)
    if(Found) JOBOLD = JOBIPH
  endif
  If(Found) iJOB = 1
  If(iJOB == 1) Then
    if(JOBOLD <= 0) Then
      JOBOLD = 20
      Call DaName(JOBOLD,'JOBOLD')
    endif
  endif
  IADR19(:) = 0
  Call IDaFile(JOBOLD,2,IADR19,15,IAD19)
  IADR15 = IADR19
  dmDisk = IADR19(3)

  Call mma_allocate(casdm1,NACPAR,Label='casdm1')
  Call mma_allocate(dm1_cas,NTOT2,Label='dm1_cas')
  Call mma_allocate(casdm1s,NACPAR,Label='casdm1s')
  Call mma_allocate(dm1s,NTOT2,Label='dm1s')

  Call mma_allocate(folded_dm1,nTot1,Label='folded_dm1')
  Call mma_allocate(folded_dm1_cas,nTot1,Label='folded_dm1_cas')
  Call mma_allocate(folded_dm1s,nTot1,Label='folded_dm1s')
  Call mma_allocate(P2d,NACPR2,Label='P2D')

  Call mma_allocate(coul,ntot1,Label='coul')

  ! This iSA is used to control gradient calculations.  Analytic gradients
  ! (in ALASKA) will only run if iSA=1, and iSA will only be set to one if
  ! the on-top potentials are computed as part of this calculation.
  iSA = 99
  Call Put_iScalar('SA ready',iSA)

  ! We loop over the number of states.  For each state, we read the density
  ! matrices from the JOBIPH file.
  do jroot = 1,lroots
    ! Read in the density matrices for <jroot>.
    casdm1(:) = zero
    dm1_cas(:) = zero
    casdm1s(:) = zero
    dm1s(:) = zero
    folded_dm1(:) = zero
    folded_dm1s(:) = zero
    folded_dm1_cas(:) = zero
    P2D(:) = zero

    ! Get the D1 Active matrix for this state.  These should probably be
    ! most easily read from the previous JOBIPH file.  Then, convert D1A from
    ! the MO to the AO basis.

    Call DDaFile(JOBOLD,2,casdm1,NACPAR,dmDisk)
    Call DDaFile(JOBOLD,2,casdm1s,NACPAR,dmDisk)
    Call DDaFile(JOBOLD,2,P2d,NACPR2,dmDisk)
    ! This dummy read is to cycle the dmDisk so next
    ! time we encounter this block, we load the next
    ! states densities
    Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)

    ! This should be in energy_ot() but it gets
    ! modified in dblock. But also, why do we need
    ! this in the DFT section???
    Call Put_dArray('D1mo',casdm1,NACPAR)

    ! Generate total density
    If(NASH(1) /= NAC) then
      Call dblock(casdm1)
      call dblock(casdm1s)
    endif
    Call Get_D1A_RASSCF(CMO,casdm1,dm1_cas)
    Call Get_D1A_RASSCF(CMO,casdm1s,dm1s)

    Call Fold(nSym,nBas,dm1_core,folded_dm1)
    Call Fold(nSym,nBas,dm1_cas,folded_dm1_cas)
    folded_dm1(:) = folded_dm1(:)+folded_dm1_cas

    Call Fold(nSym,nBas,dm1s,folded_dm1s)

!Maybe I can write all of these matrices to file, then modify stuff in
!the nq code to read in the needed density.  In other words, I need to
!replace the next call with something that supports multistates.

    ! save things for MSPDFT gradients
    IF(mcpdft_options%grad .and. mcpdft_options%mspdft) THEN
      call fold2(nsym,nbas,dm1_cas,DIDA(:,jroot))
      d1aoms(:,jroot) = folded_dm1(:)
      Call P2_contraction(casdm1,P2MOt(:,jroot))
      if(ispin /= 1) then
        d1saoms(:,jroot) = folded_dm1s(:)
      endif
    ENDIF

    do_pdftpot = (mcpdft_options%grad .and. (mcpdft_options%mspdft .or. (jroot == mcpdft_options%rlxroot)))
    e_ot = mcpdft_options%otfnal%energy_ot(folded_dm1,folded_dm1s,casdm1,P2d,charge)

    call get_coulomb(cmo,dm1_core,dm1_cas,coul)

    e_wfn = energy_mcwfn(folded_dm1,hcore,coul,PotNuc,ntot1)

    e_state = e_wfn+e_ot

    IF(mcpdft_options%otfnal%is_hybrid()) THEN
      e_state = mcpdft_options%otfnal%lambda*e_mcscf(jRoot)+(one-mcpdft_options%otfnal%lambda)*e_state
    ENDIF

    Call Print_MCPDFT_2(PotNuc,e_wfn,e_ot,jroot,e_mcscf(jroot))

    ! JB replacing e_mcscf with MC-PDFT energy for MS-PDFT use
    e_states(jroot) = e_state

    ! At this point, the energy calculation is done.  Now I need to build the
    ! fock matrix if this root corresponds to the relaxation root.
    if(mcpdft_options%grad) then
      ! Determine size of Q matrix
      NQ = 0
      NIAIA = 0
      do ISYM = 1,NSYM
        NQ = MAX(NQ,NASH(ISYM)*NORB(ISYM))
        NIAIA = NIAIA+(NASH(ISYM)+NISH(ISYM))**2
      enddo
      if(NQ < NIAIA) NQ = NIAIA

      if((.not. mcpdft_options%mspdft) .and. jroot == mcpdft_options%rlxroot) then
        call savefock_pdft(cmo,hcore(:)+coul(:),casdm1,nq,p2d)
      endif
      if(mcpdft_options%mspdft) then
        call savefock_mspdft(CMO,hcore(:)+coul(:),casdm1,NQ,p2d,jroot)
      endif
    endif
  enddo !loop over roots

  if(mcpdft_options%grad) then
    dmDisk = IADR19(3)
    do jroot = 1,mcpdft_options%rlxroot-1
      Call DDaFile(JOBOLD,0,casdm1,NACPAR,dmDisk)
      Call DDaFile(JOBOLD,0,casdm1s,NACPAR,dmDisk)
      Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
      Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
    enddo
    Call DDaFile(JOBOLD,2,casdm1,NACPAR,dmDisk)
!        Andrew added this line to fix heh2plus
    Call DDaFile(JOBOLD,2,casdm1s,NACPAR,dmDisk)
    Call Put_dArray('D1mo',casdm1,NACPAR)
    Call DDaFile(JOBOLD,2,P2d,NACPR2,dmDisk)
    Call Put_dArray('P2mo',P2d,NACPR2)
    ! from before???
    Call mma_allocate(P2t,NACPR2,Label='P2t')
    P2t(:) = zero
    Call P2_contraction(casdm1,P2t)
    Call Put_dArray('P2MOt',P2t,NACPR2)
    Call mma_deallocate(P2t)

    If(NASH(1) /= NAC) Call DBLOCK(casdm1)
    Call Get_D1A_RASSCF(CMO,casdm1,dm1_cas)

    Call Fold(nSym,nBas,dm1_core,folded_dm1)
    Call Fold(nSym,nBas,dm1_cas,folded_dm1_cas)
    folded_dm1(:) = folded_dm1(:)+folded_dm1_cas(:)
    Call Put_dArray('D1ao',folded_dm1,nTot1)

    !  Generate spin-density
    dm1s(:) = zero
    IF(NASH(1) /= NAC) CALL DBLOCK(casdm1s)
    Call Get_D1A_RASSCF(CMO,casdm1s,dm1s)
    Call Fold(nSym,nBas,dm1s,folded_dm1s)
    Call Put_dArray('D1Sao',folded_dm1s,nTot1)

    Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
  endif

  Call mma_deallocate(folded_dm1_cas)
  Call mma_deallocate(folded_dm1)
  Call mma_deallocate(folded_dm1s)
  Call mma_deallocate(casdm1)
  Call mma_deallocate(dm1_cas)
  Call mma_deallocate(casdm1s)
  Call mma_deallocate(dm1s)
  call mma_deallocate(hcore)
  Call mma_deallocate(coul)
  Call mma_deallocate(P2D)
  Call mma_deallocate(dm1_core)

  ! This deallocation SHOULD NOT BE HERE
  ! but it has to do with how the FuncExtParams is loaded...
  If(Allocated(FuncExtParams)) Call mma_deallocate(FuncExtParams)

ENDSubroutine compute_mcpdft_energy
