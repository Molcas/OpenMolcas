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
! Copyright (C) 2018, Andrew M. Sand                                   *
!               2019, Thais R. Scott                                   *
!               2021, Jie J. Bao                                       *
!               2024, Matthew R. Hennefarth                            *
!***********************************************************************

Subroutine SaveFock_PDFT(cmo,hcore,coul,D1Act,NQ,p2d,state)
  use definitions,only:iwp,wp,u6
  use constants,only:zero
  use printlevel,only:debug
  use mcpdft_output,only:iPrLoc
  use stdalloc,only:mma_allocate,mma_deallocate
  use wadr,only:fockocc
  use rasscf_global,only:nFint,ISTORP,ntot4
  use mspdftgrad,only:F1MS,F2MS,FocMS,FxyMS
  use general_data,only:ntot1,nbas,nfro,norb,nsym
  implicit none

  real(kind=wp),intent(in) :: cmo(*),D1Act(*),hcore(*),coul(*),P2D(*)
  integer(kind=iwp),intent(in) :: NQ,state

  integer(kind=iwp) :: iSA,iprlev
  real(kind=wp),allocatable :: ONTOPT(:),ONTOPO(:)
  real(kind=wp),allocatable :: FA_V(:),FI_V(:)
  real(kind=wp),allocatable :: Q(:),dm2(:),fock(:),h1e(:)

  call mma_allocate(fock,ntot4,label='fock')
  call mma_allocate(h1e,ntot1,label='h1e')

  fock(:) = zero

  write(u6,'(2X,A)') 'Calculating potentials for analytic gradients for MS-PDFT'

  IPRLEV = IPRLOC(3)

  h1e(:) = hcore(:ntot1)+coul(:ntot1)
  call ao2mo_1e(cmo,h1e,h1e,nsym,nbas,norb,nfro)

  ! loading one-electron potential and two-electron potential
  ! Used as F1 and F2 in equations 58 and 59 in Ref1.
  Call mma_allocate(ONTOPT,nfint,Label='OnTopT')
  Call mma_allocate(ONTOPO,ntot1,Label='OnTopO')
  Call Get_dArray('ONTOPT',OnTopT,NFINT)
  Call Get_dArray('ONTOPO',OnTopO,NTOT1)

  ! Store for later...
  Call Get_TUVX(OnTopT,f2ms(:,state))

  CALL mma_allocate(FI_V,Ntot1,Label='FI_V')
  CALL mma_allocate(FA_V,Ntot1,Label='FA_V')
  ! Note that these are stored in MO basis
  Call Get_dArray('FI_V',FI_V,NTOT1)
  Call Get_dArray('FA_V',FA_V,NTOT1)

  fi_v(:) = h1e(:)+OnTopO(:)+FI_V(:)
  F1MS(:,state) = fi_v(:)

  ! ____________________________________________________________
  ! This next part is to generate the MC-PDFT generalized fock operator.
  ! The corrections (from the potentials) to FI and FA are built in the NQ
  ! part of the code, for efficiency's sake.  It still needs to be
  ! debugged.

  !Reordering of the two-body density matrix.
  IF(ISTORP(NSYM+1) > 0) THEN
    call mma_allocate(dm2,ISTORP(NSYM+1),label="dm2")
    CALL PMAT_RASSCF(P2d,dm2)
  else
    call mma_allocate(dm2,1,label="dm2")
    dm2(:) = zero
  ENDIF

  !Must add to existing fock operator (occ/act).
  CALL mma_allocate(Q,NQ,Label='Q') ! q-matrix(1symmblock)
  CALL fock_update(fock,fi_v,fa_v,D1Act,dm2,Q,OnTopT,CMO)
  Call mma_deallocate(Q)
  call mma_deallocate(dm2)
  Call mma_deallocate(OnTopO)
  Call mma_deallocate(OnTopT)
  CALL mma_deallocate(FI_V)
  CALL mma_deallocate(FA_V)

  focms(:,state) = fockocc(:)
  FxyMS(:,state) = fock(:)
  IF(IPRLEV >= DEBUG) THEN
    write(u6,*) 'FOCC_OCC'
    call wrtmat(fockocc,1,ntot1,1,ntot1)
  ENDIF

  iSA = 1
  Call Put_iScalar('SA ready',iSA)

  call mma_deallocate(fock)
  call mma_deallocate(h1e)

EndSubroutine SaveFock_PDFT
