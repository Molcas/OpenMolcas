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
! Copyright (C) 2007, Bingbing Suo                                     *
!***********************************************************************
! 11 feb 2007 -bsuo- added by suo bing
!                    this file contains subroutines which depend on plat

subroutine version_info()

use Definitions, only: u6

implicit none

write(u6,'(10x,a42)') '*****************************************'
write(u6,'(10x,a42)') '*      Xian-ci mrci program             *'
write(u6,'(10x,a42)') '*     Institute of Modern Physics       *'
write(u6,'(10x,a42)') '*        Northwest University           *'
write(u6,'(10x,a42)') '*        xian, shaanxi, china           *'
write(u6,'(10x,a42)') '*                                       *'
write(u6,'(10x,a42)') '*        report bugs and errors         *'
write(u6,'(10x,a42)') '*           wzy@nwu.edu.cn              *'
write(u6,'(10x,a42)') '*        yubin_wang@hotmail.com         *'
write(u6,'(10x,a42)') '*       bingbing_suo@hotmail.com        *'
write(u6,'(10x,a42)') '*                                       *'
write(u6,'(10x,a42)') '*****************************************'
write(u6,*)
write(u6,*)

return

end subroutine version_info

subroutine allocate_int_memory()

use gugaci_global, only: maxintseg, vint_ci
use stdalloc, only: mma_allocate
use Constants, only: Zero

implicit none

! this subroutine is used to allocate the dynamic memory for integrals
! vint_ci(:) pointer to the base address of the memory of integrals
! maxintseg  maximum length of the integral segment

call mma_allocate(vint_ci,[0,maxintseg+1],label='vint_ci')
vint_ci(:) = Zero

return

end subroutine allocate_int_memory

subroutine deallocate_int_memory()

use gugaci_global, only: vint_ci
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(vint_ci)

return

end subroutine deallocate_int_memory

subroutine read_ml(nf,bv,n,m)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nf, n, m
real(kind=wp), intent(out) :: bv(n)
integer(kind=iwp) :: idisk, irec(64)

idisk = 0
call idafile(nf,2,irec,64,idisk)
idisk = irec(m)
call ddafile(nf,2,bv,n,idisk)

!select case (nf)
!  case (34)
!    call ddafile(nf,2,bv,n,idisk)
!  case (35)
!    call ddafile(nf,2,bv,n,idisk)
!  case default
!end select

return

end subroutine read_ml

subroutine write_ml(nf,bv,n,m)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nf, n, m
real(kind=wp), intent(_IN_) :: bv(n)
integer(kind=iwp) :: idisk, irec(64)

idisk = 0
if (m == 1) then
  irec = 0
  call idafile(nf,1,irec,64,idisk)
  irec(1) = idisk
else
  call idafile(nf,2,irec,64,idisk)
end if

idisk = irec(m)
call ddafile(nf,1,bv,n,idisk)

!select case (nf)
!  case (34)
!    call ddafile(nf,1,bv,n,idisk)
!  case (35)
!    call ddafile(nf,1,bv,n,idisk)
!  case default
!end select
irec(m+1) = idisk
idisk = 0
call idafile(nf,1,irec,64,idisk)
!write(u6,*) 'idisk',m,idisk

return

end subroutine write_ml

subroutine readint(ntyp,vintrd)
! subroutine used to read divided integrals into main memory
! on entry:
! -----------------
! ntyp   - integral type need to be readed
!   1 internal space integrals
!   2 integrals which have one or three indices in external space
!   3 integrals which have two indices in external space
!   4 external space integrals
! on out:
! -------------------
! vintrd(*) - integrals used in hamiltonian matrix calculation

use gugaci_global, only: LuCiInt

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ntyp
real(kind=wp), intent(_OUT_) :: vintrd(*)
integer(kind=iwp) :: idisk, idum(1), idx(4), lenint

idx = 0
idisk = 0
call idafile(luciint,2,idx,4,idisk)
select case (ntyp)
  case (1)
    idisk = idx(1)
    call idafile(luciint,2,idum,1,idisk)
    lenint = idum(1)
    call ddafile(luciint,2,vintrd(2),lenint,idisk)
  case (2)
    ! internal external space 1 and 3 index integrals
    idisk = idx(2)
    call idafile(luciint,2,idum,1,idisk)
    lenint = idum(1)
    call ddafile(luciint,2,vintrd(2),lenint,idisk)
    !write(u6,*) vintrd(1:lenint)
  case (3)
    ! internal external space 2 index integrals
    idisk = idx(3)
    call idafile(luciint,2,idum,1,idisk)
    lenint = idum(1)
    call ddafile(luciint,2,vintrd(2),lenint,idisk)
  case (4)
    ! external space integrals
    idisk = idx(4)
    call idafile(luciint,2,idum,1,idisk)
    lenint = idum(1)
    call ddafile(luciint,2,vintrd(2),lenint,idisk)
end select

!do i=1,lenint+1
!  write(12,'(i8,1x,f15.8)') i,vintrd(i)
!end do

return

end subroutine readint

subroutine gugaciinit()

use gugaci_global, only: FnOneMO, FnTwoMO, LuCiDen, LuCiDia, LuCiInt, LuCiMO, LuCiTv1, LuCiTv2, LuCiVec, LuDrt, LuLoop, LuOneMO, &
                         LuTwoMO
#ifdef MOLPRO
use file_qininit, only:
use molegeom, only:
use groupinfo, only:
use orbinfo, only:
use Definitions, only: iwp
#endif

implicit none
character(len=8) :: fnciden, fncidia, fnciint, fncimo, fncitv1, fncitv2, fncivec, fndrt, fnloop

#ifndef MOLPRO
fnonemo = 'traone'
#endif
fntwomo = 'traint'
fndrt = 'cidrt'
fnciint = 'ciint'
fnloop = 'ciloop'
fncidia = 'cidia'
fncivec = 'civec'
fncitv1 = 'citv1'
fncitv2 = 'citv2'
fnciden = 'ciden'
fncimo = 'cimo'

luonemo = 30
lutwomo = 50
ludrt = 31
luciint = 32
luloop = 33
lucidia = 34
lucivec = 35
lucitv1 = 36
lucitv2 = 37
luciden = 38
lucimo = 39
#ifndef MOLPRO
call daname(ludrt,fndrt)
call daname(luciint,fnciint)
call daname(luloop,fnloop)
call daname(lucidia,fncidia)
call daname(lucivec,fncivec)
call daname(lucitv1,fncitv1)
call daname(lucitv2,fncitv2)
call daname(luciden,fnciden)
call daname(lucimo,fncimo)
#endif
#ifdef _XIANEST_
! open chkfile and drt file
call fileopen(nfchk,anchk,11)
call chkfil_taskctrl(2)  ! read task information

nreps = 0
! read infomation from checkfile and moint file
call chkfil_ciorbinf(2)
nlsm_bas(1:8) = nsbas(1:8)
nlsm_all(1:8) = nsorb(1:8)
ng_sm = nreps

call fileopen(ludrt,fndrt,12)
call fileopen(luciint,fnciint,12)
call fileopen(luloop,fnloop,12)
call fileopen(lucidia,fncidia,12)
call fileopen(lucivec,fncivec,12)
call fileopen(lucitv1,fncitv1,12)
call fileopen(lucitv2,fncitv2,12)
call fileopen(luciden,fnciden,12)
call fileopen(lucimo,fncimo,12)
call fileopen(lutwomo,fntwomo,12)
#endif

return

end subroutine gugaciinit

subroutine gugafinalize()

use gugaci_global, only: LuCiDen, LuCiDia, LuCiInt, LuCiMO, LuCiTv1, LuCiTv2, LuCiVec, LuDrt, LuLoop

implicit none

#ifdef MOLPRO
call fileclos(ludrt,12)
call fileclos(luciint,12)
call fileclos(luloop,12)
call fileclos(lucidia,12)
call fileclos(lucivec,12)
call fileclos(lucitv1,12)
call fileclos(lucitv2,12)
call fileclos(luciden,12)
call fileclos(lucimo,12)
call fileclos(lutwomo,12)
#else
call daclos(ludrt)
call daclos(luciint)
call daclos(luloop)
call daclos(lucidia)
call daclos(lucivec)
call daclos(lucitv1)
call daclos(lucitv2)
call daclos(luciden)
call daclos(lucimo)
#endif

return

end subroutine gugafinalize

subroutine memdengrad_alloc()

use gugaci_global, only: denm1, denm2, norb_all
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nc, ndim

ndim = norb_all*(norb_all+1)/2
call mma_allocate(denm1,ndim,label='denm1')
nc = ndim*(ndim+1)/2
call mma_allocate(denm2,nc,label='denm2')

end subroutine memdengrad_alloc

subroutine memdengrad_free()

use gugaci_global, only: denm1, denm2
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(denm1)
call mma_deallocate(denm2)

end subroutine memdengrad_free

subroutine memcidiag_alloc()

use gugaci_global, only: jeh, jph, jwh, maxpl, th, thh
use stdalloc, only: mma_allocate
use Constants, only: Zero

implicit none

call mma_allocate(jph,maxpl,label='jph')
call mma_allocate(jeh,maxpl,label='jeh')
call mma_allocate(jwh,maxpl,label='jwh')
call mma_allocate(th,maxpl,label='th')
call mma_allocate(thh,maxpl,label='thh')
jph(:) = 0
jeh(:) = 0
jwh(:) = 0
th(:) = Zero
thh(:) = Zero

return

end subroutine memcidiag_alloc

subroutine memcidiag_dealloc()

use gugaci_global, only: jeh, jph, jwh, th, thh
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(jph)
call mma_deallocate(jeh)
call mma_deallocate(jwh)
call mma_deallocate(th)
call mma_deallocate(thh)

return

end subroutine memcidiag_dealloc

subroutine mem_intinnindex_alloc()

use gugaci_global, only: intind_abkk, intind_iabc, intind_iaqq, intind_ijab, intind_ijcc, intind_ijka, intspace_abkk, &
                         intspace_ijab, intspace_ijcc, loij, loij_all, loijk, &
                         loijk_all, nabc, ngw2, ngw3, norb_all, norb_inn, vdint, voint
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: lent

call mma_allocate(loij,50000,label='loij')
call mma_allocate(loijk,1384150,label='loijk')
call mma_allocate(loij_all,50000,label='loij_all')
call mma_allocate(loijk_all,1384150,label='loijk_all')
loij(1:50000) = 0
loijk(1:1384150) = 0
loijk_all(1:50000) = 0
loijk_all(1:1384150) = 0

call mma_allocate(intind_iaqq,50000,label='intind_iaqq')
call mma_allocate(intind_abkk,50000,label='intind_abkk')
lent = norb_inn*nabc+norb_all+ngw2(norb_all)+ngw3(norb_all) !2*norb_ext*norb_ext*norb_ext
call mma_allocate(intind_iabc,lent,label='intind_iabc')
call mma_allocate(intind_ijka,50000,label='intind_ijka')
call mma_allocate(intind_ijcc,50000,label='intind_ijcc')
call mma_allocate(intind_ijab,50000,label='intind_ijab')
intind_iaqq(1:50000) = 0
intind_abkk(1:50000) = 0
intind_iabc(1:lent) = 0
intind_ijka(1:50000) = 0
intind_ijcc(1:50000) = 0
intind_ijab(1:50000) = 0

call mma_allocate(intspace_abkk,50000,label='inspace_abkk')
call mma_allocate(intspace_ijcc,50000,label='inspace_ijcc')
call mma_allocate(intspace_ijab,50000,label='inspace_ijab')
intspace_abkk(1:50000) = 0
intspace_ijcc(1:50000) = 0
intspace_ijab(1:50000) = 0

voint = Zero
vdint = Zero

return

end subroutine mem_intinnindex_alloc

subroutine mem_intinnindex_dealloc()

use gugaci_global, only: intind_abkk, intind_iabc, intind_iaqq, intind_ijab, intind_ijcc, intind_ijka, intspace_abkk, &
                         intspace_ijab, intspace_ijcc, loij, loij_all, loijk, loijk_all
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(loij)
call mma_deallocate(loijk)
call mma_deallocate(loij_all)
call mma_deallocate(loijk_all)

call mma_deallocate(intind_iaqq)
call mma_deallocate(intind_abkk)
call mma_deallocate(intind_iabc)
call mma_deallocate(intind_ijka)
call mma_deallocate(intind_ijcc)
call mma_deallocate(intind_ijab)

call mma_deallocate(intspace_abkk)
call mma_deallocate(intspace_ijcc)
call mma_deallocate(intspace_ijab)

return

end subroutine mem_intinnindex_dealloc

subroutine allocate_casrst()

use gugaci_global, only: ja, jb, jj, jm, kk, max_node
use stdalloc, only: mma_allocate

implicit none

call mma_allocate(ja,max_node,label='ja')
call mma_allocate(jb,max_node,label='jb')
call mma_allocate(jm,[0,max_node],label='jm')
call mma_allocate(jj,[1,4],[0,max_node],label='jj')
call mma_allocate(kk,[0,max_node],label='kk')
ja(:) = 0
jb(:) = 0
jm(:) = 0
jj(:,:) = 0
kk(:) = 0

return

end subroutine allocate_casrst

subroutine deallocate_casrst()

use gugaci_global, only: ja, jb, jj, jm, kk
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(ja)
call mma_deallocate(jb)
call mma_deallocate(jm)
call mma_deallocate(jj)
call mma_deallocate(kk)

return

end subroutine deallocate_casrst

subroutine allocate_subdrt(icase,lent)

use gugaci_global, only: ihy, iy, jj_sub, jphy, max_node, max_wei
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: icase, lent

call mma_allocate(ihy,max_wei,label='ihy')
call mma_allocate(jj_sub,[1,4],[0,max_node],label='jj_sub')
call mma_allocate(iy,[1,4],[0,max_node],label='iy')
if (icase == 1) then
  call mma_allocate(jphy,max_node,label='jphy')
else
  call mma_allocate(jphy,lent,label='jphy')
end if

end subroutine allocate_subdrt

subroutine allocate_subdrtl(icase,lent)

use gugaci_global, only: ihyl, iyl, jjl_sub, jphyl, max_node, max_wei
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: icase, lent

call mma_allocate(ihyl,max_wei,label='ihyl')
call mma_allocate(jjl_sub,[1,4],[0,max_node],label='jjl_sub')
call mma_allocate(iyl,[1,4],[0,max_node],label='iyl')
if (icase == 1) then
  call mma_allocate(jphyl,max_node,label='jphyl')
else
  call mma_allocate(jphyl,lent,label='jphyl')
end if

end subroutine allocate_subdrtl

subroutine deallocate_subdrt()

use gugaci_global, only: ihy, iy, jj_sub, jphy
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(ihy)
call mma_deallocate(jj_sub)
call mma_deallocate(iy)
call mma_deallocate(jphy)

end subroutine deallocate_subdrt

subroutine deallocate_subdrtl()

use gugaci_global, only: ihyl, iyl, jjl_sub, jphyl
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(ihyl)
call mma_deallocate(jjl_sub)
call mma_deallocate(iyl)
call mma_deallocate(jphyl)

end subroutine deallocate_subdrtl

#ifndef MOLPRO

function c_time()

use Definitions, only: wp

implicit none
real(kind=wp) :: c_time
real(kind=wp), external :: seconds

c_time = seconds()

end function c_time

function ipair(i,j)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ipair
integer(kind=iwp), intent(in) :: i, j

if (i >= j) then
  ipair = i*(i-1)/2+j
else
  ipair = j*(j-1)/2+i
end if

end function ipair

!subroutine trimstr(string)
!! delete space character in the head and tail of the string
!
!use Definitions, only: iwp
!
!implicit none
!character(len=*), intent(inout) :: string
!integer(kind=iwp) :: i, j, k
!character(len=128) :: line
!
!k = len_trim(string)
!line(1:k) = string(1:k)
!do i=1,k
!  if (string(i:i) /= ' ') exit
!end do
!string = ' '
!j = k-i+1
!string(1:j) = line(i:k)
!
!return
!
!end subroutine trimstr

#endif
