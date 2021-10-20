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

subroutine seg_drt()

use gugaci_global, only: ihy, ipae, iy, jj, jj_sub, jpad, jpae, jphy, max_wei, ndim, no, nohy, norb_act, norb_dz, norb_inn
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i, ihypos, ij1, in_, iw, j(4), j1, j2, j3, j4, jde, jp, jp0, jpe, jpend, jpmax, jps, jpsta, lr
integer(kind=iwp), allocatable :: jpihy(:)

jpmax = no(norb_inn+1)
iy(:,:) = 0
ihy(:) = 0
jj_sub(:,:) = 0
jphy(:) = 0

iy(1,jpad) = 1
ndim = 0
if (norb_act == 0) then
  if (jpad /= ipae) return
  ndim = 1
  jj_sub(1,jpad) = jpae
  return
end if
do i=1,4
  j(i) = jj(i,jpad)
  if (j(i) /= 0) then
    iy(1,j(i)) = 1
    jj_sub(i,jpad) = j(i)
  end if
end do

jps = 1
jpe = no(norb_inn+1)
do jp=jps,jpe
  jj_sub(1:4,jp) = 0
  if (iy(1,jp) == 0) cycle
  do i=1,4
    j(i) = jj(i,jp)
    if (j(i) == 0) cycle
    iy(1,j(i)) = 1
    jj_sub(i,jp) = j(i)
  end do
end do

if (iy(1,jpae) == 0) then
  ndim = 0
  nohy = 0
  return
end if

iy(1,1:jpmax) = 0
iy(1,jpae) = 1

do jp0=no(norb_inn-1)+1,no(norb_inn)
  if (jj_sub(1,jp0) /= jpae) jj_sub(1,jp0) = 0
  if (jj_sub(2,jp0) /= jpae) jj_sub(2,jp0) = 0
  if (jj_sub(3,jp0) /= jpae) jj_sub(3,jp0) = 0
  if (jj_sub(4,jp0) /= jpae) jj_sub(4,jp0) = 0
end do
do lr=norb_inn-1,norb_dz+1,-1
  jps = no(lr)+1
  jpe = no(lr+1)
  do jde=jps,jpe
    j1 = jj_sub(1,jde)
    j2 = jj_sub(2,jde)
    j3 = jj_sub(3,jde)
    j4 = jj_sub(4,jde)
    iy(1,jde) = iy(1,j1)+iy(1,j2)+iy(1,j3)+iy(1,j4)
    if (iy(1,jde) /= 0) then
      if ((j2 /= 0) .and. (iy(1,j2) /= 0)) iy(2,jde) = iy(1,j1)
      if ((j3 /= 0) .and. (iy(1,j3) /= 0)) iy(3,jde) = iy(1,j1)+iy(1,j2)
      if ((j4 /= 0) .and. (iy(1,j4) /= 0)) iy(4,jde) = iy(1,jde)-iy(1,j4)
      !if (jde == 1) continue
    else
      do jp0=no(lr-1)+1,no(lr)
        if (jj_sub(1,jp0) == jde) jj_sub(1,jp0) = 0
        if (jj_sub(2,jp0) == jde) jj_sub(2,jp0) = 0
        if (jj_sub(3,jp0) == jde) jj_sub(3,jp0) = 0
        if (jj_sub(4,jp0) == jde) jj_sub(4,jp0) = 0
      end do
    end if
  end do
end do
lr = norb_dz+1
jde = jpad
j1 = jj_sub(1,jde)
j2 = jj_sub(2,jde)
j3 = jj_sub(3,jde)
j4 = jj_sub(4,jde)
iy(1,jde) = iy(1,j1)+iy(1,j2)+iy(1,j3)+iy(1,j4)
if (iy(1,jde) /= 0) then
  if ((j2 /= 0) .and. (iy(1,j2) /= 0)) iy(2,jde) = iy(1,j1)
  if ((j3 /= 0) .and. (iy(1,j3) /= 0)) iy(3,jde) = iy(1,j1)+iy(1,j2)
  if ((j4 /= 0) .and. (iy(1,j4) /= 0)) iy(4,jde) = iy(1,jde)-iy(1,j4)
  !if (jde == 1) continue
else
  do jp0=no(lr)+1,no(lr+1)
    if (jj_sub(1,jp0) == jde) jj_sub(1,jp0) = 0
    if (jj_sub(2,jp0) == jde) jj_sub(2,jp0) = 0
    if (jj_sub(3,jp0) == jde) jj_sub(3,jp0) = 0
    if (jj_sub(4,jp0) == jde) jj_sub(4,jp0) = 0
  end do
end if

ndim = iy(1,jpad)
!=======================================================================
jphy(jpad) = 1
ihy(1) = 1
ihy(2) = 0

ihypos = 3
do i=1,4
  j(i) = jj_sub(i,jpad)
  if (j(i) == 0) cycle
  iw = 0
  if (i /= 1) iw = iy(i,jpad)
  jphy(j(i)) = 1
  ihy(ihypos) = 1
  ihy(ihypos+1) = iw
  ihypos = ihypos+2
end do

call mma_allocate(jpihy,max_wei,label='jpihy')
lr = norb_dz+1
jpsta = no(lr)+1
jpend = jpae
do jp=jpsta,jpend
  if (iy(1,jp) == 0) cycle
  call ajphy(jp,in_,jpihy)
  jphy(jp) = ihypos
  ihy(ihypos) = in_
  do ij1=1,in_
    ihypos = ihypos+1
    ihy(ihypos) = jpihy(ij1)
  end do
  ihypos = ihypos+1
end do
nohy = ihypos-1
!write(u6,'(3x,4i10)') jpad,jpae,ndim,nohy
call mma_deallocate(jpihy)

return

end subroutine seg_drt

function iwalk_ad(jdbl,jext,iwa,iwd)

use gugaci_global, only: iseg_downwei, iw_sta, jpad_upwei, log_prod
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iwalk_ad
integer(kind=iwp), intent(in) :: jdbl, jext, iwa, iwd
integer(kind=iwp) :: ioff, isegdown, iwup

if (log_prod == 3) then
  ! mrpt2
  ioff = 0
  if (jext == 1) ioff = iw_sta(jdbl,jext)
  iwup = jpad_upwei(jdbl)
  isegdown = iseg_downwei(jext)
  iwalk_ad = (iwa*iwup+iwd)*isegdown+ioff
else
  ! mrci
  iwup = jpad_upwei(jdbl)
  isegdown = iseg_downwei(jext)
  iwalk_ad = (iwa*iwup+iwd)*isegdown+iw_sta(jdbl,jext)
end if

return

end function iwalk_ad

!function iwalk_ad_mrso(jdbl,jext,iwa,iwd)
!
!use gugaci_global, only: iseg_downwei, jpad_upwei
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp) :: iwalk_ad_mrso
!integer(kind=iwp), intent(in) :: jdbl, jext, iwa, iwd
!integer(kind=iwp) :: isegdown, iwup
!
!iwup = jpad_upwei(jdbl)
!isegdown = iseg_downwei(jext)
!iwalk_ad_mrso = (iwa*iwup+iwd)*isegdown+isub_sta(jdbl,jext)
!
!return
!
!end function iwalk_ad_mrso

subroutine ajphy(jp,in_,jpihy)

use gugaci_global, only: iy, jj, jj_sub, jpad, kk, max_node, max_wei, no, norb_dz
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: jp
integer(kind=iwp), intent(inout) :: in_, jpihy(max_wei)
integer(kind=iwp) :: i, idr, j, jpe, jpn, jy, l, lr, nn
integer(kind=iwp), allocatable :: iin(:)

if (jp == jpad) then
  in_ = 1
  jpihy(1) = 0
  return
end if
call mma_allocate(iin,[0,max_node],label='iin')
iin(0) = 0
lr = kk(jp)
!write(u6,*) '  ajphy,jp,start,end',jp,no(nst-lr)+1,no(nst-lr+1)
jpe = no(lr+1)
iin(jpad:jpe) = 0
iin(jp) = 1
do jpn=no(lr-1),jpad,-1
  iin(jpn) = iin(jj_sub(1,jpn))+iin(jj_sub(2,jpn))+iin(jj_sub(3,jpn))+iin(jj_sub(4,jpn))
end do
in_ = iin(jpad)
idr = 0
do l=1,in_
  jpihy(l) = 0
  jy = l
  nn = jpad
  do i=norb_dz+1,lr-1
    do j=1,4
      if (jj(j,nn) == 0) cycle
      if (jy <= iin(jj_sub(j,nn))) then
        idr = j
        exit
      end if
      jy = jy-iin(jj_sub(j,nn))
    end do
    if (idr /= 1) jpihy(l) = jpihy(l)+iy(idr,nn)
    nn = jj_sub(idr,nn)
  end do
end do
call mma_deallocate(iin)

return

end subroutine ajphy

function K_COE(JBL,JBR,IDDL,IDDR)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: k_coe
integer(kind=iwp), intent(in) :: jbl, jbr, iddl, iddr
integer(kind=iwp) :: idl, idr

K_COE = 0
IDL = IDDL-1
IDR = IDDR-1
if ((IDL == 3) .and. (IDR == 3)) K_COE = 200
if ((IDL == 1) .and. (IDR == 1) .and. (JBR-JBL == -1)) K_COE = -1
if ((IDL == 1) .and. (IDR == 1) .and. (JBR-JBL == 1)) K_COE = 100
if ((IDL == 2) .and. (IDR == 2) .and. (JBR-JBL == 1)) K_COE = -1
if ((IDL == 2) .and. (IDR == 2) .and. (JBR-JBL == -1)) K_COE = 100
if ((IDL == 2) .and. (IDR == 1)) K_COE = JBR
if ((IDL == 1) .and. (IDR == 2)) K_COE = -JBR-2

return

end function K_COE

subroutine copy_to_drtl()

use gugaci_global, only: ihy, ihyl, ipae, ipael, iy, iyl, jj_sub, jjl_sub, jpad, jpadl, jpae, jpael, jphy, jphyl, no, nohy, norb_inn

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: knode, nod

jpadl = jpad
jpael = jpae
ipael = ipae
ihyl(1:nohy) = ihy(1:nohy)

knode = no(norb_inn+1)
do nod=1,knode
  jjl_sub(1:4,0:knode) = jj_sub(1:4,0:knode)
  iyl(1:4,nod) = iy(1:4,nod)
  jphyl(nod) = jphy(nod)
end do

return

end subroutine copy_to_drtl

!subroutine get_jpadty(jp,ity,jp_ms)
!
!use gugaci_global, only: ns_sm
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: jp
!integer(kind=iwp), intent(out) :: ity, jp_ms
!integer(kind=iwp) :: jpp
!
!if (jp == 1) then
!  ity = 1
!  jp_ms = ns_sm
!  return
!end if
!jpp = jp+15
!if (mod(jpp,8) == 0) then
!  ity = jpp/8-1
!  jp_ms = 8
!else
!  ity = jpp/8
!  jp_ms = mod(jpp,8)
!end if
!
!return
!
!end subroutine get_jpadty

subroutine get_jpty(jpadlr,jptyl,jptyr)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: jpadlr
integer(kind=iwp), intent(out) :: jptyl, jptyr

select case (jpadlr)
  case default ! (1)
    jptyl = 4
    jptyr = 4
  case (2)
    jptyl = 4
    jptyr = 3
  case (3)
    jptyl = 3
    jptyr = 4
  case (4)
    jptyl = 4
    jptyr = 6
  case (5)
    jptyl = 6
    jptyr = 4
  case (6)
    jptyl = 4
    jptyr = 2
  case (7)
    jptyl = 2
    jptyr = 4
  case (8)
    jptyl = 4
    jptyr = 5
  case (9)
    jptyl = 5
    jptyr = 4
  case (10)
    jptyl = 4
    jptyr = 1
  case (11)
    jptyl = 3
    jptyr = 3
  case (12)
    jptyl = 6
    jptyr = 6
  case (13)
    jptyl = 3
    jptyr = 2
  case (14)
    jptyl = 2
    jptyr = 3
  case (15)
    jptyl = 6
    jptyr = 5
  case (16)
    jptyl = 5
    jptyr = 6
  case (17)
    jptyl = 3
    jptyr = 1
  case (18)
    jptyl = 6
    jptyr = 1
  case (19)
    jptyl = 2
    jptyr = 2
  case (20)
    jptyl = 5
    jptyr = 5
  case (21)
    jptyl = 2
    jptyr = 5
  case (22)
    jptyl = 5
    jptyr = 2
  case (23)
    jptyl = 2
    jptyr = 1
  case (24)
    jptyl = 5
    jptyr = 1
  case (25)
    jptyl = 1
    jptyr = 1
end select

return

end subroutine get_jpty
