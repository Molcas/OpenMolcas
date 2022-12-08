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

subroutine gugadrt_rst(id,nndd)
!*******************************************
! 10 may 2007 - revised by wyb

use gugadrt_global, only: iseg_downwei, iprint, iintbit, ja, ja_sys, jb, jb_sys, jc_sys, jd, jj, jm, jroute_sys, js, jt, jv, kk, &
                          lsm_inn, max_innorb, max_node, mxnode, n_ref, n16int, n32int, ng_sm, no, norb_act, norb_all, norb_dz, &
                          norb_inn, ns_sm, nu_ad
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: id, nndd
integer(kind=iwp) :: i, iabcbit, idd, iextbit, ii, iiabkm(1:n16int), iextii(n32int), iextjj(n32int), im, imd, ims, imt, it, &
                     ivalid, iysum, j, j1, j2, j3, j4, ja0, jac, jaj, jajk, jatmp, jb0, jbj, jbjk, jbtmp, jc0, jde, jds, ji, &
                     jjabkm(1:n16int), jk, jkabkm(1:n16int), jmj, jmjk, jmtmp, jp, jp0, jpe, jps, jq1, jq2, jq3, jq4, k0, kj1, &
                     kkj, kkjk, kktmp, kttmp, l, lr, mxtnode, nabcbit, nextbit, nm, node, nrefbit
logical(kind=iwp) :: flag
integer(kind=iwp), allocatable :: jabkm(:,:), ind(:,:), idjj(:,:), iwy(:,:), itm(:), noh(:)

! estimate memory
if (n_ref > 20) then
  if (norb_act >= 10) then
    if (iintbit == 32) mxtnode = 1000000
    if (iintbit == 64) mxtnode = 10000000
  else
    mxtnode = 1000000
  end if
else
  if (norb_act > 10) then
    mxtnode = 1000000
  else
    mxtnode = 500000
  end if
end if
write(u6,*)
write(u6,*) ' now generate distinct row tableau'
call mma_allocate(jabkm,n16int,mxtnode,label='jabkm')
call mma_allocate(ind,[1,n32int],[0,mxtnode],label='ind')
call mma_allocate(idjj,[1,4],[0,mxtnode],label='idjj')
call mma_allocate(iwy,[1,4],[0,mxtnode],label='iwy')
call mma_allocate(itm,[0,mxtnode],label='itm')
call mma_allocate(noh,max_innorb,label='noh')

nm = ns_sm
no(1:norb_dz+1) = 0
jabkm(1:n16int,1:mxtnode) = 0
ind(1:n32int,0:mxtnode) = 0
idjj(1:4,0:mxtnode) = 0
iwy = 0
itm = 0
nrefbit = n32int
jj(1:4,0:max_node) = 0

! 8 bits
iextbit = 2
nextbit = iintbit/iextbit
! 16 bits
iabcbit = 16
nabcbit = iintbit/iabcbit

! v node
j = 0
ja0 = ja_sys
jb0 = jb_sys
jc0 = jc_sys
write(u6,'(4(1x,i4))') ja0,jb0,jc0
! v_node
ja(1) = ja0
jb(1) = jb0
jm(1) = nm
kk(1) = norb_dz
jaj = ja(1)
jbj = jb(1)
jmj = jm(1)
kkj = kk(1)

jjabkm(1:n16int) = 0
call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
jabkm(1:n16int,1) = jjabkm(1:n16int)
ind(1:nrefbit,1) = 0

! d_node
imd = 0
do node=2,9
  imd = imd+1
  if (nu_ad(node) == 0) cycle
  ja(node) = ja(1)
  jb(node) = jb0+1
  jm(node) = imd
  kk(node) = norb_dz
  jaj = ja(node)
  jbj = jb(node)
  jmj = jm(node)
  kkj = kk(node)

  jjabkm(1:n16int) = 0
  call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
  jabkm(1:n16int,node) = jjabkm(1:n16int)
  ind(1:nrefbit,node) = 0
end do

! t_node
imt = 0
do node=10,17
  imt = imt+1
  if (nu_ad(node) == 0) cycle
  ja(node) = ja(1)
  jb(node) = jb0+2
  jm(node) = imt
  kk(node) = norb_dz
  jaj = ja(node)
  jbj = jb(node)
  jmj = jm(node)
  kkj = kk(node)

  jjabkm(1:n16int) = 0
  call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
  jabkm(1:n16int,node) = jjabkm(1:n16int)
  ind(1:nrefbit,node) = 0
end do
! s_node
ims = 0
do node=18,25
  ims = ims+1
  if (nu_ad(node) == 0) cycle
  ja(node) = ja(1)+1
  jb(node) = jb0
  jm(node) = ims
  kk(node) = norb_dz
  jaj = ja(node)
  jbj = jb(node)
  jmj = jm(node)
  kkj = kk(node)
  jjabkm(1:n16int) = 0
  call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
  jabkm(1:n16int,node) = jjabkm(1:n16int)
  ind(1:nrefbit,node) = 0
end do
if (jroute_sys > 1) then
  ! d'_node
  imd = 0
  do node=26,33
    imd = imd+1
    if (nu_ad(node) == 0) cycle
    ja(node) = ja(1)+1
    jb(node) = jb0-1
    jm(node) = imd
    kk(node) = norb_dz
    jaj = ja(node)
    jbj = jb(node)
    jmj = jm(node)
    kkj = kk(node)
    jjabkm(1:n16int) = 0
    call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
    jabkm(1:n16int,node) = jjabkm(1:n16int)
    ind(1:nrefbit,node) = 0
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
    jm(node) = imt
    kk(node) = norb_dz
    jaj = ja(node)
    jbj = jb(node)
    jmj = jm(node)
    kkj = kk(node)
    jjabkm(1:n16int) = 0
    call wrtabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
    jabkm(1:n16int,node) = jjabkm(1:n16int)
    ind(1:nrefbit,node) = 0
  end do
end if
no(norb_dz) = mxnode

j = 0
jk = mxnode
kk(jk) = norb_dz
!call packnod(idkk,jk,norb_dz,nabcbit,iabcbit,mxtnode)

!**********************************************************************

do
  j = j+1
  if (j <= mxnode) then
    if (nu_ad(j) == 0) cycle
  end if
  if (j < mxnode) then
    jaj = ja(j)
    jbj = jb(j)
    jmj = jm(j)
    kkj = kk(j)
    k0 = kkj
  else
    jjabkm(1:n16int) = jabkm(1:n16int,j)
    call redabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
    k0 = kkj
  end if
  jkabkm(1:n16int) = jabkm(1:n16int,jk)
  call redabkm(jkabkm,n16int,nabcbit,iabcbit,jajk,jbjk,jmjk,kkjk)
  if (kkjk /= k0+1) no(k0) = jk
  jk = jk+1

  if (jk > mxtnode) then
    write(u6,*) ' the number of j exceeds mxtnode',mxtnode
    call abend()
    !call errexit(777)
  end if
  !=========================================================
  ! ***********
  ! *   d=0   *
  ! ***********

  jatmp = jaj
  jbtmp = jbj
  jmtmp = jmj
  kktmp = kkj+1

  iextjj(1:nrefbit) = ind(1:nrefbit,j)
  call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,0,kttmp,k0)

  flag = .false.
  if (ivalid == 0) then
    jk = jk-1
    flag = .true.
  else if (k0 == norb_inn-1) then
    ! act||ext
    if (((jatmp == 1) .or. (jbtmp == 2)) .and. (kttmp /= 0)) then
      jk = jk-1
      flag = .true.
    else if ((jbtmp == 1) .and. (kttmp == 2)) then
      jk = jk-1
      flag = .true.
    end if
  end if

  if (.not. flag) then
    call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,jmtmp,kktmp)
    outer_1: do ii=no(k0)+1,jk-1
      iiabkm(1:n16int) = jabkm(1:n16int,ii)
      do i=1,n16int
        if (iiabkm(i) /= jkabkm(i)) cycle outer_1
      end do
      if (k0 /= norb_inn-1) then
        iextii(1:nrefbit) = ind(1:nrefbit,ii)
        do i=1,nrefbit
          if (iextii(i) /= iextjj(i)) cycle outer_1
        end do
      end if
      jk = jk-1
      idjj(1,j) = ii
      flag = .true.
      if (flag) exit
    end do outer_1
  end if

  if (.not. flag) then
    ind(1:nrefbit,jk) = iextjj(1:nrefbit)
    idjj(1,j) = jk
    jabkm(1:n16int,jk) = jkabkm(1:n16int)

    iiabkm(1:n16int) = jabkm(1:n16int,jk-1)
    call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
    if (kktmp /= kj1) no(k0) = jk-1
  end if

  ! =========================================================
  if (jbj /= 0) then
    jk = jk+1
    ! ***********
    ! *   d=1   *
    ! ***********
    jatmp = jaj
    jbtmp = jbj-1
    jmtmp = Mul(lsm_inn(k0+1),jmj)
    kktmp = kkj+1

    iextjj(1:nrefbit) = ind(1:nrefbit,j)
    call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,1,kttmp,k0)

    flag = .false.
    if (ivalid == 0) then
      jk = jk-1
      flag = .true.
    else if (k0 == norb_inn-1) then
      ! act||ext
      if (((jatmp == 1) .or. (jbtmp == 2)) .and. (kttmp /= 0)) then
        jk = jk-1
        flag = .true.
      else if ((jbtmp == 1) .and. (kttmp == 2)) then
        jk = jk-1
        flag = .true.
      end if
    end if

    if (.not. flag) then
      call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,jmtmp,kktmp)
      outer_2: do ii=no(k0)+1,jk-1
        iiabkm(1:n16int) = jabkm(1:n16int,ii)
        do i=1,n16int
          if (iiabkm(i) /= jkabkm(i)) cycle outer_2
        end do
        if (k0 /= norb_inn-1) then
          iextii(1:nrefbit) = ind(1:nrefbit,ii)
          do i=1,nrefbit
            if (iextii(i) /= iextjj(i)) cycle outer_2
          end do
        end if
        jk = jk-1
        idjj(2,j) = ii
        flag = .true.
        if (flag) exit
      end do outer_2
    end if

    if (.not. flag) then
      ind(1:nrefbit,jk) = iextjj(1:nrefbit)
      idjj(2,j) = jk
      jabkm(1:n16int,jk) = jkabkm(1:n16int)

      iiabkm(1:n16int) = jabkm(1:n16int,jk-1)
      call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
      if (kktmp /= kj1) no(k0) = jk-1
    end if
  end if
  ! =========================================================
  jac = norb_all-jaj-jbj
  if (jaj > 0) then
    if (jac /= 0) then
      jk = jk+1
      ! *************
      ! *    d=2    *
      ! *************
      jatmp = jaj-1
      jbtmp = jbj+1
      jmtmp = Mul(lsm_inn(k0+1),jmj)
      kktmp = kkj+1

      iextjj(1:nrefbit) = ind(1:nrefbit,j)
      call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,2,kttmp,k0)

      flag = .false.
      if (ivalid == 0) then
        jk = jk-1
        flag = .true.
      else if (k0 == norb_inn-1) then
        ! act||ext
        if (((jatmp == 1) .or. (jbtmp == 2)) .and. (kttmp /= 0)) then
          jk = jk-1
          flag = .true.
        else if ((jbtmp == 1) .and. (kttmp == 2)) then
          jk = jk-1
          flag = .true.
        end if
      end if

      if (.not. flag) then
        call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,jmtmp,kktmp)
        outer_3: do ii=no(k0)+1,jk-1
          iiabkm(1:n16int) = jabkm(1:n16int,ii)
          do i=1,n16int
            if (iiabkm(i) /= jkabkm(i)) cycle outer_3
          end do
          if (k0 /= norb_inn-1) then
            iextii(1:nrefbit) = ind(1:nrefbit,ii)
            do i=1,nrefbit
              if (iextii(i) /= iextjj(i)) cycle outer_3
            end do
          end if
          jk = jk-1
          idjj(3,j) = ii
          flag = .true.
          if (flag) exit
        end do outer_3
      end if

      if (.not. flag) then
        ind(1:nrefbit,jk) = iextjj(1:nrefbit)
        idjj(3,j) = jk
        jabkm(1:n16int,jk) = jkabkm(1:n16int)

        iiabkm(1:n16int) = jabkm(1:n16int,jk-1)
        call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
        if (kktmp /= kj1) no(k0) = jk-1
      end if
    end if
    ! =========================================================
    ! ************
    ! *    d=3   *
    ! ************

    jk = jk+1
    jatmp = jaj-1
    jbtmp = jbj
    jmtmp = jmj
    kktmp = kkj+1

    iextjj(1:nrefbit) = ind(1:nrefbit,j)
    call gugadrt_njexcit(iextjj,nrefbit,iextbit,nextbit,ivalid,3,kttmp,k0)

    flag = .false.
    if (ivalid == 0) then
      jk = jk-1
      flag = .true.
    else if (k0 == norb_inn-1) then
      ! act||ext
      if (((jatmp == 1) .or. (jbtmp == 2)) .and. (kttmp /= 0)) then
        jk = jk-1
        flag = .true.
      else if ((jbtmp == 1) .and. (kttmp == 2)) then
        jk = jk-1
        flag = .true.
      end if
    end if

    if (.not. flag) then
      call wrtabkm(jkabkm,n16int,nabcbit,iabcbit,jatmp,jbtmp,jmtmp,kktmp)
      outer_4: do ii=no(k0)+1,jk-1
        iiabkm(1:n16int) = jabkm(1:n16int,ii)
        do i=1,n16int
          if (iiabkm(i) /= jkabkm(i)) cycle outer_4
        end do
        if (k0 /= norb_inn-1) then
          iextii(1:nrefbit) = ind(1:nrefbit,ii)
          do i=1,nrefbit
            if (iextii(i) /= iextjj(i)) cycle outer_4
          end do
        end if
        jk = jk-1
        idjj(4,j) = ii
        flag = .true.
        if (flag) exit
      end do outer_4
    end if

    if (.not. flag) then
      ind(1:nrefbit,jk) = iextjj(1:nrefbit)
      idjj(4,j) = jk
      jabkm(1:n16int,jk) = jkabkm(1:n16int)

      iiabkm(1:n16int) = jabkm(1:n16int,jk-1)
      call upacknod(iiabkm,4,kj1,nabcbit,iabcbit,n16int)
      if (kktmp /= kj1) no(k0) = jk-1
    end if
  end if
  !===========================================================
  if (k0 > norb_inn-1) exit
end do

!if(iprint == 1) then
!  open(100,file='tmp.dat')
!  do i=1,jk
!    jkabkm(1:n16int)=jabkm(1:n16int,i)
!    call redabkm(jkabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
!    j1 = idjj(1,i)
!    j2 = idjj(2,i)
!    j3 = idjj(3,i)
!    j4 = idjj(4,i)
!    write(100,'(5(1x,i8))') i,j1,j2,j3,j4
!  end do
!  close(100)
!end if
!write(u6,'(10(1x,i5))') no(1:norb_all)
write(u6,*)
!************** external space  *************
id = no(norb_inn)
!write(u6,508) 'befor,no=',(no(i),i=norb_dz,norb_inn)
do idd=no(norb_inn-1)+1,id
  jjabkm(1:n16int) = jabkm(1:n16int,idd)
  call redabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
  if ((jaj == 0) .and. (jbj == 0) .and. (jmj == 1)) jv = idd
  if ((jaj == 0) .and. (jbj == 1)) then
    do i=1,ng_sm
      if (jmj == i) then
        jd(i) = idd
      end if
    end do
  end if
  if ((jaj == 0) .and. (jbj == 2)) then
    do i=1,ng_sm
      if (jmj == i) then
        jt(i) = idd
      end if
    end do
  end if
  if ((jaj == 1) .and. (jbj == 0)) then
    do i=1,ng_sm
      if (jmj == i) then
        js(i) = idd
      end if
    end do
  end if
end do

iwy(1,jv) = 1
do im=1,ng_sm
  if ((jd(im) /= 0) .and. (iseg_downwei(1+im) /= 0)) iwy(1,jd(im)) = 1
  if ((jt(im) /= 0) .and. (iseg_downwei(9+im) /= 0)) iwy(1,jt(im)) = 1
  if ((js(im) /= 0) .and. (iseg_downwei(17+im) /= 0)) iwy(1,js(im)) = 1
end do

iwy(1:4,0) = 0
do l=norb_inn-1,norb_dz,-1
  jps = no(l-1)+1
  jpe = no(l)
  do jde=jps,jpe
    j1 = idjj(1,jde)
    j2 = idjj(2,jde)
    j3 = idjj(3,jde)
    j4 = idjj(4,jde)
    iwy(1,jde) = iwy(1,j1)+iwy(1,j2)+iwy(1,j3)+iwy(1,j4)
    if (iwy(1,jde) == 0) then
      do jp0=no(l-2)+1,no(l-1)
        if (idjj(1,jp0) == jde) idjj(1,jp0) = 0
        if (idjj(2,jp0) == jde) idjj(2,jp0) = 0
        if (idjj(3,jp0) == jde) idjj(3,jp0) = 0
        if (idjj(4,jp0) == jde) idjj(4,jp0) = 0
      end do
      cycle
    end if
    if ((j2 /= 0) .and. (iwy(1,j2) /= 0)) iwy(2,jde) = iwy(1,j1)
    if ((j3 /= 0) .and. (iwy(1,j3) /= 0)) iwy(3,jde) = iwy(1,j1)+iwy(1,j2)
    if ((j4 /= 0) .and. (iwy(1,j4) /= 0)) iwy(4,jde) = iwy(1,jde)-iwy(1,j4)
    if (jde == 1) cycle
    do jp=jps,jde-1
      if (iwy(1,jp) /= iwy(1,jde)) cycle
      jq1 = idjj(1,jp)
      jq2 = idjj(2,jp)
      jq3 = idjj(3,jp)
      jq4 = idjj(4,jp)
      if ((j1 == jq1) .and. (j2 == jq2) .and. (j3 == jq3) .and. (j4 == jq4)) then
        iwy(1,jde) = 0
        do jp0=no(l-2)+1,no(l-1)
          if (idjj(1,jp0) == jde) idjj(1,jp0) = jp
          if (idjj(2,jp0) == jde) idjj(2,jp0) = jp
          if (idjj(3,jp0) == jde) idjj(3,jp0) = jp
          if (idjj(4,jp0) == jde) idjj(4,jp0) = jp
        end do
        exit
      end if
    end do
  end do
end do

it = mxnode
itm(1) = 1
noh = 0
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

do jpe=1,mxnode
  if (nu_ad(jpe) == 0) cycle
  jj(1,jpe) = idjj(1,jpe)
  jj(2,jpe) = idjj(2,jpe)
  jj(3,jpe) = idjj(3,jpe)
  jj(4,jpe) = idjj(4,jpe)
end do

do jpe=mxnode+1,id
  jp = itm(jpe)
  if (jp == 0) cycle
  j1 = idjj(1,jpe)
  j2 = idjj(2,jpe)
  j3 = idjj(3,jpe)
  j4 = idjj(4,jpe)
  jjabkm(1:n16int) = jabkm(1:n16int,jpe)
  call redabkm(jjabkm,n16int,nabcbit,iabcbit,jaj,jbj,jmj,kkj)
  l = kkj
  jds = noh(l-2)+1
  jde = noh(l-1)
  do jp0=jds,jde
    if (jj(1,jp0) == jpe) jj(1,jp0) = jp
    if (jj(2,jp0) == jpe) jj(2,jp0) = jp
    if (jj(3,jp0) == jpe) jj(3,jp0) = jp
    if (jj(4,jp0) == jpe) jj(4,jp0) = jp
  end do
  ja(jp) = jaj
  jb(jp) = jbj
  jm(jp) = jmj
  kk(jp) = kkj
  do i=1,4
    iwy(i,jp) = iwy(i,jpe)
    ji = idjj(i,jpe)
    if (itm(ji) /= 0) then
      jj(i,jp) = ji
    end if
  end do
end do

!open(10,file='rst.out')
no(norb_dz) = 0
no(norb_dz+1) = mxnode
!write(u6,*) '   end of rst, drt ..........'
!write(u6,'(2x,2i10)')norb_dz+1,no(norb_dz+1)
do lr=norb_dz+1,norb_inn
  no(lr+1) = noh(lr)
end do
jv = itm(jv)
itm(0) = 0
do im=1,8
  jd(im) = itm(jd(im))
  jt(im) = itm(jt(im))
  js(im) = itm(js(im))
end do
iysum = 0
do j=1,mxnode
  iysum = iysum+iwy(1,j)
end do
id = it
if (it /= no(norb_inn+1)) then
  write(u6,*) '   rst id is wrong!!   no(norb_inn)=',no(norb_inn),it
  call abend()
end if

write(u6,*)
!iprint=1
if (iprint == 1) then
  write(u6,*) 'guga drt'
  write(u6,506)
end if
nndd = no(norb_inn)
do j=1,id
  kk(j) = kk(j)+1
  if (iprint == 1) write(u6,510) j,kk(j),ja(j),jb(j),jm(j),jj(1,j),jj(2,j),jj(3,j),jj(4,j)
end do
write(u6,*)
write(u6,*) 'end of rst, drt ..........'

call mma_deallocate(jabkm)
call mma_deallocate(ind)
call mma_deallocate(idjj)
call mma_deallocate(iwy)
call mma_deallocate(itm)
call mma_deallocate(noh)

!open(21,file='fort.drt',form='unformatted')
!write(21) id
!write(21) ja(1:id),jb(1:id),jm(1:id)
!write(21) jj(1:4,0:id)
!write(21) kk(0:id)
!write(21) no(0:norb_inn+1)
!write(21) jv,jd(1:8),jt(1:8),js(1:8)
!close(21)

call writedrt(id)

return

506 format('          j       k       a       b      jm      j0      j1      j2      j3')
!508 format(3x,a10,1x,i5,1x,16i8)
510 format(3x,9(1x,i7))

end subroutine gugadrt_rst
