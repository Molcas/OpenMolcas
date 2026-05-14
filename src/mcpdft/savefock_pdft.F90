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

use wadr, only: fockocc
use rasscf_global, only: istorp, nacpar, nacpr2, nfint, ntot4
use general_data, only: nbas, nfro, norb, nsym, ntot1, ntot2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: cmo(ntot2), h1e(ntot1), d1act(nacpar), p2d(nacpr2)
integer(kind=iwp), intent(in) :: nq
real(kind=wp), allocatable :: dm2(:), fa_v(:), fi_v(:), fock(:), h1e_mo(:), ontopo(:), ontopt(:), q(:), tuvx_tmp(:)

call mma_allocate(h1e_mo,ntot1,label='h1e_mo')
call mma_allocate(fock,ntot4,label='fock')

fock(:) = Zero

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
call put_darray('F2_PDFT',tuvx_tmp(:),nacpr2)
call mma_deallocate(tuvx_tmp)

call mma_allocate(fi_v,ntot1,label='fi_v')
call mma_allocate(fa_v,ntot1,label='fa_v')
! Note that these are stored in MO basis
call get_darray('FI_V',fi_v,ntot1)
call get_darray('FA_V',fa_v,ntot1)

fi_v(:) = fi_v(:)+ontopo(:)+h1e_mo(:)
call put_darray('F1_PDFT',fi_v(:),ntot1)

! Now we generate generalized fock operator and fockocc

! Reordering of the 2body density matrix
if (istorp(nsym+1) > 0) then
  call mma_allocate(dm2,istorp(nsym+1),label='dm2')
  call pmat_rasscf(p2d,dm2)
else
  call mma_allocate(dm2,1,label='dm2')
  dm2(:) = Zero
end if

call mma_allocate(q,nq,label='q')
call fock_update(fock,fi_v,fa_v,d1act,dm2,q,nq,ontopt,cmo)

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

end subroutine savefock_pdft
