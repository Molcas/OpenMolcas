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

subroutine gugadrt_rcas(id,indd)

use gugadrt_global, only: iseg_downwei, iprint, ja, ja_sys, jb, jb_sys, jc_sys, jd, jj, jm, jroute_sys, js, jt, jv, kk, lsm_inn, &
                          max_innorb, max_node, max_ref, mxnode, n_electron, ng_sm, no, norb_dz, norb_inn, ns_sm, nu_ad
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: id, indd
integer(kind=iwp) :: i, idd, ii, im, imd, ims, imt, inb, it, iysum, j, j1, j2, j3, j4, ja0, jb0, jc0, jde, jds, jk, jp, jp0, jpe, &
                     jps, jq1, jq2, jq3, jq4, k0, l, lr, nc, ndj, nel, nm, node
logical(kind=iwp) :: flag
integer(kind=iwp), allocatable :: ind(:,:), itm(:), iwy(:,:), jc(:), locu(:,:), noh(:)

call mma_allocate(ind,8,max_node,label='ind')
call mma_allocate(itm,[0,max_node],label='itm')
call mma_allocate(iwy,[1,4],[0,max_node],label='iwy')
call mma_allocate(jc,max_node,label='jc')
call mma_allocate(locu,8,max_ref,label='locu')
call mma_allocate(noh,max_innorb,label='noh')
write(u6,*) ' '
write(u6,*) 'now generate distinct row tableau'
noh = 0
itm = 0
iwy = 0
ind = 0
locu = 0
jc = 0

nel = n_electron
nm = ns_sm
no(1:norb_dz-1) = 0
j = 0
ja0 = ja_sys
jb0 = jb_sys
jc0 = jc_sys
! v_node
ja(1) = ja0
jb(1) = jb0
jc(1) = jc0
kk(1) = norb_dz
do i=1,8
  ind(i,1) = 0
end do
jm(1) = nm
! d_node
imd = 0
do node=2,9
  imd = imd+1
  if (nu_ad(node) == 0) cycle
  ja(node) = ja(1)
  jb(node) = jb0+1
  jc(node) = jc0-1
  kk(node) = norb_dz
  do i=1,8
    ind(i,node) = 0
  end do
  jm(node) = imd
end do
! t_node
imt = 0
do node=10,17
  imt = imt+1
  if (nu_ad(node) == 0) cycle
  ja(node) = ja(1)
  jb(node) = jb0+2
  jc(node) = jc0-2
  kk(node) = norb_dz
  do i=1,8
    ind(i,node) = 0
  end do
  jm(node) = imt
end do
! s_node
ims = 0
do node=18,25
  ims = ims+1
  if (nu_ad(node) == 0) cycle
  ja(node) = ja(1)+1
  jb(node) = jb0
  jc(node) = jc0-1
  kk(node) = norb_dz
  do i=1,8
    ind(i,node) = 0
  end do
  jm(node) = ims
end do
if (jroute_sys > 1) then
! d'_node
  imd = 0
  do node=26,33
    imd = imd+1
    if (nu_ad(node) == 0) cycle
    ja(node) = ja(1)+1
    jb(node) = jb0-1
    jc(node) = jc0
    kk(node) = norb_dz
    do i=1,8
      ind(i,node) = 0
    end do
    jm(node) = imd
  end do
end if
if (jroute_sys > 2) then
! t'_node
  imt = 0
  do node=34,41
    imt = imt+1
    if (nu_ad(node) == 0) cycle
    ja(node) = ja(1)+2
    jb(node) = jb0-2
    jc(node) = jc0
    kk(node) = norb_dz
    do i=1,8
      ind(i,node) = 0
    end do
    jm(node) = imt
  end do
end if
no(norb_dz) = mxnode

call gugadrt_ref_gfs(nel,ndj,locu,nm)

j = 0
jk = mxnode
do
  j = j+1
  if (j <= mxnode) then
    if (nu_ad(j) == 0) cycle
  end if
  k0 = kk(j)
  if (jc(j) /= 0) then
    jk = jk+1
    ! ***********
    ! *   d=0   *
    ! ***********
    ja(jk) = ja(j)
    jb(jk) = jb(j)
    jc(jk) = jc(j)-1
    jm(jk) = jm(j)
    do i=1,8
      ind(i,jk) = ind(i,j)
    end do
    inb = 0
    if (kk(j) == norb_inn-1) then
      call gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
    end if
    if (inb > 2) then
      jk = jk-1
      jj(1,j) = 0
    else
      flag = .false.
      outer_1: do ii=no(k0)+1,jk-1
        if (ja(jk) /= ja(ii)) cycle
        if (jb(jk) /= jb(ii)) cycle
        if (jc(jk) /= jc(ii)) cycle
        if (jm(jk) /= jm(ii)) cycle
        if (kk(j) /= norb_inn-1) then
          do i=1,8
            if (ind(i,jk) /= ind(i,ii)) cycle outer_1
          end do
        end if
        jk = jk-1
        jj(1,j) = ii
        flag = .true.
        if (flag) exit
      end do outer_1
      if (.not. flag) then
        jj(1,j) = jk
        kk(jk) = kk(j)+1
        if (kk(jk) /= kk(jk-1)) then
          no(k0) = jk-1
        end if
      end if
    end if
  end if
  ! v=0
  if (jb(j) /= 0) then
    jk = jk+1
    ! ***********
    ! *   d=1   *
    ! ***********
    ja(jk) = ja(j)
    jb(jk) = jb(j)-1
    jc(jk) = jc(j)
    jm(jk) = Mul(lsm_inn(kk(j)+1),jm(j))
    do i=1,8
      ind(i,jk) = ind(i,j)
    end do
    ind(lsm_inn(kk(j)+1),jk) = ind(lsm_inn(kk(j)+1),j)+1
    inb = 0
    if (kk(j) == norb_inn-1) then
      call gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
    end if
    if (inb > 2) then
      jk = jk-1
      jj(2,j) = 0
    else
      flag = .false.
      outer_2: do ii=no(k0)+1,jk-1
        if (ja(jk) /= ja(ii)) cycle
        if (jb(jk) /= jb(ii)) cycle
        if (jc(jk) /= jc(ii)) cycle
        if (jm(jk) /= jm(ii)) cycle
        if (kk(j) /= norb_inn-1) then
          do i=1,8
            if (ind(i,jk) /= ind(i,ii)) cycle outer_2
          end do
        end if
        jk = jk-1
        jj(2,j) = ii
        flag = .true.
        if (flag) exit
      end do outer_2
      if (.not. flag) then
        jj(2,j) = jk
        kk(jk) = kk(j)+1
        if (kk(jk) /= kk(jk-1)) then
          no(k0) = jk-1
        end if
      end if
    end if
  end if
  ! v=0
  if (ja(j) > 0) then
    if (jc(j) /= 0) then
      jk = jk+1
      ! *************
      ! *    d=2    *
      ! *************
      ja(jk) = ja(j)-1
      jb(jk) = jb(j)+1
      jc(jk) = jc(j)-1
      !-----------------------------------------------------------------------
      jm(jk) = Mul(lsm_inn(kk(j)+1),jm(j))
      do i=1,8
        ind(i,jk) = ind(i,j)
      end do
      ind(lsm_inn(kk(j)+1),jk) = ind(lsm_inn(kk(j)+1),j)+1
      inb = 0
      if (kk(j) == norb_inn-1) then
        call gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
      end if
      if (inb > 2) then
        jk = jk-1
        jj(3,j) = 0
      else
        flag = .false.
        outer_3: do ii=no(k0)+1,jk-1
          if (ja(jk) /= ja(ii)) cycle
          if (jb(jk) /= jb(ii)) cycle
          if (jc(jk) /= jc(ii)) cycle
          if (jm(jk) /= jm(ii)) cycle
          if (kk(j) /= norb_inn-1) then
            do i=1,8
              if (ind(i,jk) /= ind(i,ii)) cycle outer_3
            end do
          end if
          jk = jk-1
          jj(3,j) = ii
          flag = .true.
          if (flag) exit
        end do outer_3
        if (.not. flag) then
          jj(3,j) = jk
          kk(jk) = kk(j)+1
          if (kk(jk) /= kk(jk-1)) then
            no(k0) = jk-1
          end if
        end if
      end if
    end if
    ! v=0
    ! ************
    ! *    d=3   *
    ! ************
    jk = jk+1
    ja(jk) = ja(j)-1
    jb(jk) = jb(j)
    jc(jk) = jc(j)
    jm(jk) = jm(j)
    do i=1,8
      ind(i,jk) = ind(i,j)
    end do
    ind(lsm_inn(kk(j)+1),jk) = ind(lsm_inn(kk(j)+1),j)+2
    inb = 0
    if (kk(j) == norb_inn-1) then
      call gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
    end if
    if (inb > 2) then
      jk = jk-1
      jj(4,j) = 0
    else
      flag = .false.
      outer_4: do ii=no(k0)+1,jk-1
        if (ja(jk) /= ja(ii)) cycle
        if (jb(jk) /= jb(ii)) cycle
        if (jc(jk) /= jc(ii)) cycle
        if (jm(jk) /= jm(ii)) cycle
        if (kk(j) /= norb_inn-1) then
          do i=1,8
            if (ind(i,jk) /= ind(i,ii)) cycle outer_4
          end do
        end if
        jk = jk-1
        jj(4,j) = ii
        flag = .true.
        if (flag) exit
      end do outer_4
      if (.not. flag) then
        jj(4,j) = jk
        kk(jk) = kk(j)+1
        if (kk(jk) /= kk(jk-1)) then
          no(k0) = jk-1
        end if
      end if
    end if
    ! v=0
    if (jk > max_node) then
      write(u6,*) '    the number of j exceeds max_node',max_node
      call abend()
    end if
  end if
  if (kk(j) > norb_inn-1) exit
end do
!************** external space  *************

id = no(norb_inn)
write(u6,*)
write(u6,*) '    id=no(norb_inn)',id
do idd=no(norb_inn-1)+1,id
  if ((ja(idd) == 0) .and. (jb(idd) == 0) .and. (jm(idd) == 1)) jv = idd
  if ((ja(idd) == 0) .and. (jb(idd) == 1)) then
    if (jm(idd) == 1) jd(1) = idd
    if (jm(idd) == 2) jd(2) = idd
    if (jm(idd) == 3) jd(3) = idd
    if (jm(idd) == 4) jd(4) = idd
    if (jm(idd) == 5) jd(5) = idd
    if (jm(idd) == 6) jd(6) = idd
    if (jm(idd) == 7) jd(7) = idd
    if (jm(idd) == 8) jd(8) = idd
  end if
  if ((ja(idd) == 0) .and. (jb(idd) == 2)) then
    if (jm(idd) == 1) jt(1) = idd
    if (jm(idd) == 2) jt(2) = idd
    if (jm(idd) == 3) jt(3) = idd
    if (jm(idd) == 4) jt(4) = idd
    if (jm(idd) == 5) jt(5) = idd
    if (jm(idd) == 6) jt(6) = idd
    if (jm(idd) == 7) jt(7) = idd
    if (jm(idd) == 8) jt(8) = idd
  end if
  if ((ja(idd) == 1) .and. (jb(idd) == 0)) then
    if (jm(idd) == 1) js(1) = idd
    if (jm(idd) == 2) js(2) = idd
    if (jm(idd) == 3) js(3) = idd
    if (jm(idd) == 4) js(4) = idd
    if (jm(idd) == 5) js(5) = idd
    if (jm(idd) == 6) js(6) = idd
    if (jm(idd) == 7) js(7) = idd
    if (jm(idd) == 8) js(8) = idd
  end if
end do

iwy(1,jv) = 1
do im=1,ng_sm
  if ((jd(im) /= 0) .and. (iseg_downwei(1+im) /= 0)) iwy(1,jd(im)) = 1
  if ((jt(im) /= 0) .and. (iseg_downwei(9+im) /= 0)) iwy(1,jt(im)) = 1
  if ((js(im) /= 0) .and. (iseg_downwei(17+im) /= 0)) iwy(1,js(im)) = 1
end do

do i=1,4
  iwy(i,0) = 0
end do
do l=norb_inn-1,norb_dz,-1
  jps = no(l-1)+1
  jpe = no(l)
  do jde=jps,jpe
    j1 = jj(1,jde)
    j2 = jj(2,jde)
    j3 = jj(3,jde)
    j4 = jj(4,jde)
    iwy(1,jde) = iwy(1,j1)+iwy(1,j2)+iwy(1,j3)+iwy(1,j4)
    if (iwy(1,jde) == 0) then
      if (l == norb_dz) then
        nc = 0
      else
        nc = no(l-2)
      end if
      do jp0=nc+1,no(l-1)
        if (jj(1,jp0) == jde) jj(1,jp0) = 0
        if (jj(2,jp0) == jde) jj(2,jp0) = 0
        if (jj(3,jp0) == jde) jj(3,jp0) = 0
        if (jj(4,jp0) == jde) jj(4,jp0) = 0
      end do
      cycle
    end if
    if ((j2 /= 0) .and. (iwy(1,j2) /= 0)) iwy(2,jde) = iwy(1,j1)
    if ((j3 /= 0) .and. (iwy(1,j3) /= 0)) iwy(3,jde) = iwy(1,j1)+iwy(1,j2)
    if ((j4 /= 0) .and. (iwy(1,j4) /= 0)) iwy(4,jde) = iwy(1,jde)-iwy(1,j4)
    if (jde == 1) cycle
    do jp=jps,jde-1
      if (iwy(1,jp) /= iwy(1,jde)) cycle
      jq1 = jj(1,jp)
      jq2 = jj(2,jp)
      jq3 = jj(3,jp)
      jq4 = jj(4,jp)
      if ((j1 == jq1) .and. (j2 == jq2) .and. (j3 == jq3) .and. (j4 == jq4)) then
        iwy(1,jde) = 0
        do jp0=no(l-2)+1,no(l-1)
          if (jj(1,jp0) == jde) jj(1,jp0) = jp
          if (jj(2,jp0) == jde) jj(2,jp0) = jp
          if (jj(3,jp0) == jde) jj(3,jp0) = jp
          if (jj(4,jp0) == jde) jj(4,jp0) = jp
        end do
        exit
      end if
    end do
  end do
end do
it = mxnode
itm(1) = 1
noh(norb_dz) = mxnode
do lr=norb_dz,norb_inn-1
  do jp=no(lr)+1,no(lr+1)
    itm(jp) = 0
    if (iwy(1,jp) /= 0) then
      it = it+1
      itm(jp) = it
    end if
  end do
  noh(lr+1) = it
end do

do jpe=mxnode+1,id
  jp = itm(jpe)
  if (jp == jpe) cycle
  if (jp == 0) cycle
  l = kk(jpe)
  jds = noh(l-2)+1
  jde = noh(l-1)
  do jp0=jds,jde
    if (jj(1,jp0) == jpe) jj(1,jp0) = jp
    if (jj(2,jp0) == jpe) jj(2,jp0) = jp
    if (jj(3,jp0) == jpe) jj(3,jp0) = jp
    if (jj(4,jp0) == jpe) jj(4,jp0) = jp
  end do
  ja(jp) = ja(jpe)
  jb(jp) = jb(jpe)
  jc(jp) = jc(jpe)
  jm(jp) = jm(jpe)
  kk(jp) = kk(jpe)
  do i=1,4
    iwy(i,jp) = iwy(i,jpe)
    jj(i,jp) = jj(i,jpe)
    ind(i,jp) = ind(i,jpe)
  end do
  do i=5,8
    ind(i,jp) = ind(i,jpe)
  end do
end do

!open(10,file='rcas.out')
no(norb_dz) = 0
no(norb_dz+1) = mxnode

!write(nf10,'(2x,2i10)')norb_dz+1,no(norb_dz+1)
do lr=norb_dz,norb_inn
  no(lr+1) = noh(lr)
  write(u6,'(2x,2i10)') lr+1,no(lr+1)
end do

itm(0) = 0
jv = itm(jv)
do im=1,8
  jd(im) = itm(jd(im))
  jt(im) = itm(jt(im))
  js(im) = itm(js(im))
end do
id = it
if (it /= no(norb_inn+1)) then
  write(u6,*) '   rcas id is wrong!!   no(norb_inn)=',no(norb_inn),it
  call abend()
end if
iysum = 0
do j=1,mxnode
  iysum = iysum+iwy(1,j)
end do
!write(u6,*) '    end of rcas , node=',id,'  dimension=',iysum
write(u6,*)
indd = no(norb_inn)
!iprint=1
if (iprint == 1) then
  write(u6,*) 'guga drt'
  write(u6,506)
end if

do j=1,id
  kk(j) = kk(j)+1
  if (iprint == 1) then
    write(u6,507) j,kk(j),ja(j),jb(j),jm(j),jj(1,j),jj(2,j),jj(3,j),jj(4,j),iwy(2,j),iwy(3,j),iwy(4,j),iwy(1,j),(ind(i,j),i=1,8)
  end if
end do
write(u6,*) 'end of rcas, drt ..........'
write(u6,*)

!open(21,file='fort.drt',form='unformatted')
!write(21) id
!write(21) ja(1:id),jb(1:id),jm(1:id)
!write(21) jj(1:4,0:id)
!write(21) kk(0:id)
!write(21) no(0:norb_inn+1)
!write(21) jv,jd(1:8),jt(1:8),js(1:8)
!close(21)

call writedrt(id)
call mma_deallocate(ind)
call mma_deallocate(itm)
call mma_deallocate(iwy)
call mma_deallocate(jc)
call mma_deallocate(locu)
call mma_deallocate(noh)

506 format('       j    k   a  b jm    j0   j1   j2   j3         y1        y2        y3         x ind')
507 format(3x,2i5,1x,3i3,1x,4i5,1x,4i10,1x,8i2)

end subroutine gugadrt_rcas
