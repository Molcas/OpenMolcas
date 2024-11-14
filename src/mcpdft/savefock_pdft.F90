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
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine savefock_pdft(cmo,h1e,d1act,nq,p2d)
  use definitions,only:iwp,wp,u6
  use constants,only:zero
  use stdalloc,only:mma_allocate,mma_deallocate
  use wadr,only:fockocc
  use rasscf_global,only:istorp,ntot4,nfint,nacpr2
  use general_data,only:ntot1,nbas,nfro,norb,nsym

  implicit none

  real(kind=wp),intent(in) :: cmo(*),h1e(*),d1act(*),p2d(*)
  integer(kind=iwp),intent(in) :: nq

  real(kind=wp),allocatable :: ontopo(:),ontopt(:),h1e_mo(:),tuvx_tmp(:),fi_v(:),fa_v(:),dm2(:),fock(:),q(:)
  call mma_allocate(h1e_mo,ntot1,label='h1e_mo')
  call mma_allocate(fock,ntot4,label='fock')

  fock(:) = zero

  write(u6,'(2X,A)') 'Calculating potentials for analytical gradients for MC-PDFT'

  call ao2mo_1e(cmo,h1e,h1e_mo,nsym,nbas,norb,nfro)

  ! Loading 1e and 2e potentials
  call mma_allocate(ontopt,nfint,label='ontopt')
  call mma_allocate(ontopo,ntot1,label='ontopo')
  call get_darray('ONTOPT',ontopt,nfint)
  call get_darray('ONTOPO',ontopo,ntot1)

  ! store for latter..
  call mma_allocate(tuvx_tmp,nacpr2,label='tuvx_tmp')
  call get_tuvx(ontopt,tuvx_tmp)
  call put_darray('F2_PDFT         ',tuvx_tmp(:),nacpr2)
  call mma_deallocate(tuvx_tmp)

  call mma_allocate(fi_v,ntot1,label='fi_v')
  call mma_allocate(fa_v,ntot1,label='fa_v')
  ! Note that these are stored in MO basis
  call get_darray('FI_V',fi_v,size(fi_v))
  call get_darray('FA_V',fa_v,size(fa_v))

  fi_v(:) = fi_v(:)+ontopo(:)+h1e_mo(:)
  call put_darray('F1_PDFT         ',fi_v(:),ntot1)

  ! Now we generate generalized fock operator and fockocc

  ! Reordering of the 2body density matrix
  if(istorp(nsym+1) > 0) then
    call mma_allocate(dm2,istorp(nsym+1),label='dm2')
    call pmat_rasscf(p2d,dm2)
  else
    call mma_allocate(dm2,1,label='dm2')
    dm2(:) = zero
  endif

  call mma_allocate(q,nq,label='q')
  call fock_update(fock,fi_v,fa_v,d1act,dm2,q,ontopt,cmo)

  call put_darray('FockOcc',fockocc,ntot1)
  call put_darray('Fock_PDFT',fock,ntot4)

  call mma_deallocate(q)
  call mma_deallocate(fock)
  call mma_deallocate(h1e_mo)
  call mma_deallocate(fi_v)
  call mma_deallocate(fa_v)
  call mma_deallocate(ontopo)
  call mma_deallocate(ontopt)
  call mma_deallocate(dm2)

  call put_iscalar('SA ready',1)

endsubroutine savefock_pdft
