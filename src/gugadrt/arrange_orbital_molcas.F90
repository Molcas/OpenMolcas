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

subroutine arrange_orbital_molcas()
!****************************************************
! arrange orbital for ci calculation, in meld and
! molcas program, the orbital are arranged as the symmetry
! block. we transfer them to ci order
! -----
! map_order_orbital    ab ---> ci

use gugadrt_global, only: lsm_inn, max_orb, ng_sm, nlsm_all, norb_all, norb_dz, norb_inn
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i, iccount, im, j, la, lr, lr_scf, lr_scf0, lra, lrd, lsmid, lsmorbcount(ng_sm), lsmr, ms, nim
logical(kind=iwp) :: logi_norb_inn(norb_all), logic_assign_actorb
integer(kind=iwp), allocatable :: map_orb_order(:), map_tmp(:)

call mma_allocate(map_orb_order,max_orb,label='map_orb_order')

logi_norb_inn(1:norb_all) = .false.

nim = 0
lsmorbcount(1) = nim
do im=2,ng_sm
  nim = nim+nlsm_all(im-1)
  lsmorbcount(im) = nim
end do

! FIXME: logic_assign_actorb was undefined
logic_assign_actorb = .false.
if (logic_assign_actorb) then
  do lr=1,norb_inn
    lr_scf = map_orb_order(lr)
    logi_norb_inn(lr_scf) = .true.
  end do
else
  do lr=1,norb_inn
    lsmr = lsm_inn(lr)
    lsmorbcount(lsmr) = lsmorbcount(lsmr)+1
    lr_scf = lsmorbcount(lsmr)
    map_orb_order(lr) = lr_scf
    logi_norb_inn(lr_scf) = .true.
  end do
end if
lr_scf0 = norb_all
la = norb_inn+1
do ms=ng_sm,1,-1
  lr_scf0 = lr_scf0-nlsm_all(ms)
  lr_scf = lr_scf0
  do lra=1,nlsm_all(ms)
    lr_scf = lr_scf+1
    if (logi_norb_inn(lr_scf)) cycle
    map_orb_order(la) = lr_scf
    la = la+1
  end do
end do

iccount = 1
do lrd=1,norb_inn
  !ipwt(lrd) = iccount
  iccount = iccount+2
end do
lsmorbcount = 0
do lrd=norb_dz,1,-1
  lsmid = lsm_inn(lrd)
  lsmorbcount(lsmid) = lsmorbcount(lsmid)+1
  !ipws(lrd)=(lsmorbcount(lsmid)-1)*3+1
end do

call mma_allocate(map_tmp,max_orb,label='map_tmp')

map_tmp(1:norb_all) = map_orb_order(1:norb_all)
do i=1,norb_all
  do j=1,norb_all
    if (map_tmp(j) == i) then
      map_orb_order(i) = j
      exit
    end if
  end do
end do

call mma_deallocate(map_tmp)

call mma_deallocate(map_orb_order)

!write(u6,*) 'map_order_orbit'
!write(u6,1001) map_orb_order(1:norb_all)
return

!1001 format(20(1x,i3))

end subroutine arrange_orbital_molcas
