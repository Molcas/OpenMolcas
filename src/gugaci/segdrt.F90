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

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "stdalloc.fh"
!common/casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
dimension j(4)
allocatable jpihy(:)

jpmax = no(norb_inn+1)
iy = 0
ihy = 0
jj_sub = 0
jphy = 0

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
    if (iy(1,jde) /= 0) goto 304
    do jp0=no(lr-1)+1,no(lr)
      if (jj_sub(1,jp0) == jde) jj_sub(1,jp0) = 0
      if (jj_sub(2,jp0) == jde) jj_sub(2,jp0) = 0
      if (jj_sub(3,jp0) == jde) jj_sub(3,jp0) = 0
      if (jj_sub(4,jp0) == jde) jj_sub(4,jp0) = 0
    end do
    goto 21
304 continue
    if ((j2 == 0) .or. (iy(1,j2) == 0)) goto 31
    iy(2,jde) = iy(1,j1)
31  continue
    if ((j3 == 0) .or. (iy(1,j3) == 0)) goto 32
    iy(3,jde) = iy(1,j1)+iy(1,j2)
32  continue
    if ((j4 == 0) .or. (iy(1,j4) == 0)) goto 303
    iy(4,jde) = iy(1,jde)-iy(1,j4)
303 continue
    if (jde == 1) goto 21
21  continue
  end do
end do
lr = norb_dz+1
jde = jpad
j1 = jj_sub(1,jde)
j2 = jj_sub(2,jde)
j3 = jj_sub(3,jde)
j4 = jj_sub(4,jde)
iy(1,jde) = iy(1,j1)+iy(1,j2)+iy(1,j3)+iy(1,j4)
if (iy(1,jde) /= 0) goto 404
do jp0=no(lr)+1,no(lr+1)
  if (jj_sub(1,jp0) == jde) jj_sub(1,jp0) = 0
  if (jj_sub(2,jp0) == jde) jj_sub(2,jp0) = 0
  if (jj_sub(3,jp0) == jde) jj_sub(3,jp0) = 0
  if (jj_sub(4,jp0) == jde) jj_sub(4,jp0) = 0
end do
goto 40
404 continue
if ((j2 == 0) .or. (iy(1,j2) == 0)) goto 41
iy(2,jde) = iy(1,j1)
41 continue
if ((j3 == 0) .or. (iy(1,j3) == 0)) goto 42
iy(3,jde) = iy(1,j1)+iy(1,j2)
42 continue
if ((j4 == 0) .or. (iy(1,j4) == 0)) goto 403
iy(4,jde) = iy(1,jde)-iy(1,j4)
403 continue
if (jde == 1) goto 40
40 continue

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
  call ajphy(jp,in,jpihy)
  jphy(jp) = ihypos
  ihy(ihypos) = in
  do ij1=1,in
    ihypos = ihypos+1
    ihy(ihypos) = jpihy(ij1)
  end do
  ihypos = ihypos+1
end do
nohy = ihypos-1
!write(6,'(3x,4i10)') jpad,jpae,ndim,nohy
call mma_deallocate(jpihy)

return

end subroutine seg_drt

function iwalk_ad(jdbl,jext,iwa,iwd)

#include "drt_h.fh"
#include "intsort_h.fh"

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

function iwalk_ad_mrso(jdbl,jext,iwa,iwd)

#include "drt_h.fh"
#include "intsort_h.fh"

iwup = jpad_upwei(jdbl)
isegdown = iseg_downwei(jext)
iwalk_ad_mrso = (iwa*iwup+iwd)*isegdown+isub_sta(jdbl,jext)

return

end function iwalk_ad_mrso

subroutine ajphy(jp,in,jpihy)

#include "drt_h.fh"
#include "pl_structure_h.fh"
!common/casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common/sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
dimension iin(0:max_node), jpihy(max_wei)

iin(0) = 0
if (jp == jpad) then
  in = 1
  jpihy(1) = 0
  return
end if
lr = kk(jp)
!write(6,*) '  ajphy,jp,start,end',jp,no(nst-lr)+1,no(nst-lr+1)
jpe = no(lr+1)
iin(jpad:jpe) = 0
iin(jp) = 1
do jpn=no(lr-1),jpad,-1
  iin(jpn) = iin(jj_sub(1,jpn))+iin(jj_sub(2,jpn))+iin(jj_sub(3,jpn))+iin(jj_sub(4,jpn))
end do
in = iin(jpad)
idr = 0
do l=1,in
  jpihy(l) = 0
  jy = l
  nn = jpad
  do i=norb_dz+1,lr-1
    do j=1,4
      if (jj(j,nn) == 0) goto 40
      if (jy > iin(jj_sub(j,nn))) goto 353
      idr = j
      goto 350
353   continue
      jy = jy-iin(jj_sub(j,nn))
40    continue
    end do
350 continue
    if (idr /= 1) jpihy(l) = jpihy(l)+iy(idr,nn)
    nn = jj_sub(idr,nn)
  end do
end do

return

end subroutine ajphy

function K_COE(JBL,JBR,IDDL,IDDR)

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

#include "drt_h.fh"
#include "pl_structure_h.fh"

jpadl = jpad
jpael = jpae
ipael = ipae
ndiml = ndim
nohyl = nohy
ihyl(1:nohyl) = ihy(1:nohy)

knode = no(norb_inn+1)
do nod=1,knode
  jjl_sub(1:4,0:knode) = jj_sub(1:4,0:knode)
  iyl(1:4,nod) = iy(1:4,nod)
  jphyl(nod) = jphy(nod)
end do

return

end subroutine copy_to_drtl

subroutine get_jpadty(jp,ity,jp_ms)

#include "drt_h.fh"

if (jp == 1) then
  ity = 1
  jp_ms = ns_sm
  return
end if
jpp = jp+15
if (mod(jpp,8) == 0) then
  ity = jpp/8-1
  jp_ms = 8
else
  ity = jpp/8
  jp_ms = mod(jpp,8)
end if

return

end subroutine get_jpadty

subroutine get_jpty(jpadlr,jptyl,jptyr)

goto(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25),jpadlr
1 continue
jptyl = 4
jptyr = 4
return
2 continue
jptyl = 4
jptyr = 3
return
3 continue
jptyl = 3
jptyr = 4
return
4 continue
jptyl = 4
jptyr = 6
return
5 continue
jptyl = 6
jptyr = 4
return
6 continue
jptyl = 4
jptyr = 2
return
7 continue
jptyl = 2
jptyr = 4
return
8 continue
jptyl = 4
jptyr = 5
return
9 continue
jptyl = 5
jptyr = 4
return
10 continue
jptyl = 4
jptyr = 1
return
11 continue
jptyl = 3
jptyr = 3
return
12 continue
jptyl = 6
jptyr = 6
return
13 continue
jptyl = 3
jptyr = 2
return
14 continue
jptyl = 2
jptyr = 3
return
15 continue
jptyl = 6
jptyr = 5
return
16 continue
jptyl = 5
jptyr = 6
return
17 continue
jptyl = 3
jptyr = 1
return
18 continue
jptyl = 6
jptyr = 1
return
19 continue
jptyl = 2
jptyr = 2
return
20 continue
jptyl = 5
jptyr = 5
return
21 continue
jptyl = 2
jptyr = 5
return
22 continue
jptyl = 5
jptyr = 2
return
23 continue
jptyl = 2
jptyr = 1
return
24 continue
jptyl = 5
jptyr = 1
return
25 continue
jptyl = 1
jptyr = 1
return

end subroutine get_jpty
