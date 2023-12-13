!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine gugaci(ireturn)

use gugaci_global, only: iw_downwei, iw_sta, lenvec, logic_calpro, logic_grad, LuCiDia, max_node, max_orb, max_vector, mroot, &
                         nci_dim, nci_h0, norb_all, vcm, vector1, vector2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: maxplcon, mxvec, nc, nc0, nc1, npl
real(kind=wp) :: sc0, sc1, sc2
real(kind=wp), external :: c_time

iw_downwei(:,:) = 0
iw_sta(:,:) = 0

ireturn = 100
call version_info()
! open scratch files and read from strach file
call gugaciinit()

! user input file
logic_grad = .false.
call mole_inf()

call paras_calculate()
call arrange_orbital()

call allocate_casrst()
call dbl_upwalk()         ! add by wyb 01.9.5
call ext_downwalk()       ! add by wyb 01.9.5
call active_drt()         ! add by wyb 01.9.5
call value_of_pl_in_dbl()

nc0 = norb_all*(norb_all+1)/2
nc1 = nc0*(nc0+1)/2
if (nc1 > 2*max_vector) then
  write(u6,*) 'Not enough space to store MO integrals! number of orbitals should be less than ',max_orb
# ifndef MOLPRO
  call abend()
# endif
end if

call mem_intinnindex_alloc()
lenvec = nc1
call mma_allocate(vector1,nc1,label='vector1')
vector1(1:nc1) = Zero
call int_sort()
call mma_deallocate(vector1)

mxvec = 60000000*10
if (mroot*nci_dim <= mxvec) then
  nc1 = mroot*nci_dim
else
  nc1 = nci_dim
end if
call mma_allocate(vector1,nc1,label='vector1')
!call mma_allocate(vector2,nc1,label='vector2')
lenvec = nc1
vector1(1:nc1) = Zero

!write(u6,*) '==================================================='

sc1 = c_time()

call allocate_subdrt(1,1)
call allocate_subdrtl(1,1)

call memcidiag_alloc()
call diagonal_loop_wyb()
call memcidiag_dealloc()
sc2 = c_time()
!write(u6,*) vector1(585)
!call abend()

write(u6,*)
write(u6,*) '==================================================='
write(u6,'(a30,i10,f14.2,a1)') '   end of h_diagonal, nci_dim=',nci_dim,sc2-sc1,'s'
write(u6,*) '==================================================='
write(u6,*)

call write_ml(lucidia,vector1,nci_dim,1)
!do i=1,74902
!  write(21,'(i8,f18.8)') i,vector1(i)
!end do
!call abend()

call allocate_vplp_memory()
call allocate_int_memory()

nc = nci_h0 !iw_sta(2,1)
nc1 = nc*(nc+1)/2

call mma_allocate(vcm,mroot*nci_h0,label='vcm')
if (nc1 <= lenvec) then
  call mma_allocate(vector2,lenvec,label='vector2')
else
  call mma_deallocate(vector1)
  call mma_allocate(vector1,nc*(nc+1)/2,label='vector1')
  call mma_allocate(vector2,nc*(nc+1)/2,label='vector2')
  vector1(:) = Zero
  call read_ml(Lucidia,vector1,nci_dim,1)
end if

vector2(1:nc1) = Zero

call geth0()
call xflush(u6)

if (nc1 > lenvec) then
  call mma_deallocate(vector1)
  call mma_deallocate(vector2)
  call mma_allocate(vector1,lenvec,label='vector1')
  call mma_allocate(vector2,lenvec,label='vector2')
end if

sc0 = c_time()
call guga_ploop(npl,maxplcon)

call deallocate_subdrt()
call deallocate_subdrtl()

sc1 = c_time()
call xflush(u6)

write(u6,'(a25,i10,f14.2,a1)') '  end of pl_serach, n_pl=',npl,sc1-sc0,'s'
write(u6,*) '=============================================='

if (maxplcon < max_node) maxplcon = max_node
call allocate_subdrt(2,maxplcon)
call allocate_subdrtl(2,maxplcon)

call cidiagonalize(mxvec)

sc2 = c_time()

call xflush(u6)

write(u6,910) sc2-sc0
write(u6,*)
call deallocate_int_memory()
if (allocated(vcm)) call mma_deallocate(vcm)
call mma_deallocate(vector1)
call mma_deallocate(vector2)

if (logic_calpro) then
  logic_grad = .true.
  call memdengrad_alloc()

  nc0 = norb_all*(norb_all+1)/2
  nc1 = nc0*(nc0+1)/2
  call mma_allocate(vector1,nci_dim,label='vector1')
  call mma_allocate(vector2,nc1,label='vector2')
  vector1(:) = Zero
  vector2(:) = Zero

  call cidenmat()

  ! calculate properties
# ifndef MOLPRO
  call cipro()
# endif
  call mma_deallocate(vector1)
  call mma_deallocate(vector2)
  call memdengrad_free()
end if

call deallocate_vplp_memory()
call deallocate_subdrt()
call deallocate_subdrtl()
call deallocate_casrst()
call mem_intinnindex_dealloc()
call gugafinalize()

ireturn = 0

return

910 format(/,1x,'end of ci energy calculation, takes ',f10.2,1x,'seconds'/)

end subroutine gugaci
