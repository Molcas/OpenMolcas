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

subroutine arrange_orbital()
!****************************************************
!  arrange orbital for ci calculation, in meld and
!  molcas program, the orbitals are arranged as the symmetry
!  block. we transfer them to ci order
! -----
! map_order_orbital    ab ---> ci

use gugaci_global, only: jp2, jp3, logic_assign_actorb, lsm_inn, map_orb_order, max_orb, ng_sm, nlsm_all, norb_all, norb_inn, &
                         norb_number !, norb_dz
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i, im, iorb, isum2, isum3, j, la, lr, lr_scf, lr_scf0, lra, lsmorbcount(ng_sm), lsmr, map_tmp(max_orb), ms, nim
logical(kind=iwp) :: logi_norb_inn(norb_all)

logi_norb_inn(1:norb_all) = .false.
iorb = norb_all
do la=1,norb_all
  norb_number(la) = iorb
  iorb = iorb-1
end do

nim = 0
lsmorbcount(1) = nim
do im=2,ng_sm
  nim = nim+nlsm_all(im-1)
  lsmorbcount(im) = nim
end do

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

isum2 = 0
isum3 = 0
do i=1,ng_sm
  jp2(i) = isum2
  isum2 = isum2+i
  jp3(i) = isum3
  isum3 = isum3+isum2
end do

!iccount = 1
!do lrd=1,norb_inn
!  ipwt(lrd) = iccount
!  iccount = iccount+2
!end do
!lsmorbcount = 0
!do lrd=norb_dz,1,-1
!  lsmid = lsm_inn(lrd)
!  lsmorbcount(lsmid) = lsmorbcount(lsmid)+1
!  ipws(lrd) = (lsmorbcount(lsmid)-1)*3+1
!end do

map_tmp(1:norb_all) = map_orb_order(1:norb_all)
do i=1,norb_all
  do j=1,norb_all
    if (map_tmp(j) == i) then
      map_orb_order(i) = j
      exit
    end if
  end do
end do

return

end subroutine arrange_orbital
