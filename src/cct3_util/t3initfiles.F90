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

subroutine t3initfiles(length)
! this routine distributes work space WRK for required files
! for fix mediates it defines also mapd and mapi, for help mediates
! it estimates their length and distribute WRK (i.e. def pos0 parameters)
!
! length - overal requirements of work space (O)
!
! !N.B. This routine cannot run with +OP2 level

use CCT3_global, only: dp1, dp2, fk1, fk2, fk3, fk4, fk5, fk6, h1, h2, h3, l1, l2, m1, m2, m3, mchntyp, n, noa, norb, nsym, nvb, &
                       posd0, px, rx1, rx2, rx3, t11, t12, t21, t22, t23, vx, w11, w12, w13, w14, w21, w22, w23, wx
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: length
integer(kind=iwp) :: maxnoa, maxnorb, maxnvb, nhelp1, nhelp2, post, sizeh, sizel, sizem, sizen, sizer, sizew, symp

!1 maps and positions for fix mediated

!1.0 maps for DP - diagonal part
!    N.B. DP has one degree of freedom, while other 1 index has none
!    DP1 - dp(p)a
!    DP2 - dp(p)b

dp1%i(1:nsym,1:nsym,1:nsym) = 0
dp2%i(1:nsym,1:nsym,1:nsym) = 0

post = 1

dp1%pos0 = post
dp1%d(0,1) = 5
dp1%d(0,2) = 0
dp1%d(0,3) = 0
dp1%d(0,4) = 0
dp1%d(0,5) = nsym
dp1%d(0,6) = 0

dp1%d(1:nsym,4:6) = 1
do symp=1,nsym
  dp1%d(symp,1) = post
  dp1%d(symp,2) = norb(symp)
  dp1%d(symp,3) = symp
  dp1%i(symp,1,1) = symp
  post = post+norb(symp)
end do

dp2%pos0 = post
dp2%d(0,1) = 5
dp2%d(0,2) = 0
dp2%d(0,3) = 0
dp2%d(0,4) = 0
dp2%d(0,5) = nsym
dp2%d(0,6) = 0

dp2%d(1:nsym,4:6) = 1
do symp=1,nsym
  dp2%d(symp,1) = post
  dp2%d(symp,2) = norb(symp)
  dp2%d(symp,3) = symp
  dp2%i(symp,1,1) = symp
  post = post+norb(symp)
end do

!1.1 maps for T1
!    T11 - t1oaa(a,i)
!    T12 - t1obb(a,i)

t11%pos0 = post
call cct3_grc0(2,0,3,1,0,0,1,t11,post)
t12%pos0 = post
call cct3_grc0(2,0,4,2,0,0,1,t12,post)

!1.5 maps for FK
!    FK1 - f(a,b)aa
!    FK2 - f(a,b)bb
!    FK3 - f(a,i)aa
!    FK4 - f(a,i)bb
!    FK5 - f(i,j)aa
!    FK6 - f(i,j)bb

fk1%pos0 = post
call cct3_grc0(2,0,3,3,0,0,1,fk1,post)
fk2%pos0 = post
call cct3_grc0(2,0,4,4,0,0,1,fk2,post)
fk3%pos0 = post
call cct3_grc0(2,0,3,1,0,0,1,fk3,post)
fk4%pos0 = post
call cct3_grc0(2,0,4,2,0,0,1,fk4,post)
fk5%pos0 = post
call cct3_grc0(2,0,1,1,0,0,1,fk5,post)
fk6%pos0 = post
call cct3_grc0(2,0,2,2,0,0,1,fk6,post)

!1.6 maps for T2
!    T21 - t2o(ab,ij)aaaa
!    T22 - t2o(ab,ij)bbbb
!    T23 - t2o(a,b,i,j)abab

t21%pos0 = post
call cct3_grc0(4,4,3,3,1,1,1,t21,post)
t22%pos0 = post
call cct3_grc0(4,4,4,4,2,2,1,t22,post)
t23%pos0 = post
call cct3_grc0(4,0,3,4,1,2,1,t23,post)

!1.8 maps for W1
!    W11 - <ie||mn>aaaa
!    W12 - <ie||mn>bbbb
!    W13 - <ie||mn>abab
!    W14 - <ie||mn>baab

w11%pos0 = post
call cct3_grc0(4,3,1,3,1,1,1,w11,post)
w12%pos0 = post
call cct3_grc0(4,3,2,4,2,2,1,w12,post)
w13%pos0 = post
call cct3_grc0(4,0,1,4,1,2,1,w13,post)
w14%pos0 = post
call cct3_grc0(4,0,2,3,1,2,1,w14,post)

!1.9 maps for W2
!    W21 - <ab||ij>aaaa
!    W22 - <ab||ij>bbbb
!    W23 - <a,b|i,j>abab

w21%pos0 = post
call cct3_grc0(4,4,3,3,1,1,1,w21,post)
w22%pos0 = post
call cct3_grc0(4,4,4,4,2,2,1,w22,post)
w23%pos0 = post
call cct3_grc0(4,0,3,4,1,2,1,w23,post)

!2 for help files maps are irrelevant,
!  here only estimation of maximal length is done to
!  define pos0 of help files
!  we have:
!  2  W,V files - of vv2 type
!  2    L files - of vvv (vvo) type
!  3    R files - of vv2+ type
!  3    M files - of vv (vo)  type
!  3    H files - of v (o)  type
!  2  N,P files - of nn   type

wx%pos0 = post

!2.* find maxsize of W,L,M,H
sizew = 0
sizel = 0
sizer = 0
sizem = 0
sizeh = 0

do nhelp1=1,nsym

  ! W,V files
  call cct3_t3grc0(3,2,4,4,4,0,nhelp1,wx,post)
  nhelp2 = post-wx%pos0
  if (nhelp2 > sizew) then
    sizew = nhelp2
  end if

  ! L files
  call cct3_t3grc0(3,0,4,4,4,0,nhelp1,wx,post)
  nhelp2 = post-wx%pos0
  if (nhelp2 > sizel) then
    sizel = nhelp2
  end if
  call cct3_t3grc0(3,0,1,4,4,0,nhelp1,wx,post)
  nhelp2 = post-wx%pos0
  if (nhelp2 > sizel) then
    sizel = nhelp2
  end if

  ! R files
  call cct3_t3grc0(3,8,4,4,4,0,nhelp1,wx,post)
  nhelp2 = post-wx%pos0
  if (nhelp2 > sizer) then
    sizer = nhelp2
  end if

  ! M files
  call cct3_t3grc0(2,0,4,4,0,0,nhelp1,wx,post)
  nhelp2 = post-wx%pos0
  if (nhelp2 > sizem) then
    sizem = nhelp2
  end if
  call cct3_t3grc0(2,0,1,4,0,0,nhelp1,wx,post)
  nhelp2 = post-wx%pos0
  if (nhelp2 > sizem) then
    sizem = nhelp2
  end if

  ! H files
  call cct3_t3grc0(1,0,4,0,0,0,nhelp1,wx,post)
  nhelp2 = post-wx%pos0
  if (nhelp2 > sizeh) then
    sizeh = nhelp2
  end if
  call cct3_t3grc0(1,0,1,0,0,0,nhelp1,wx,post)
  nhelp2 = post-wx%pos0
  if (nhelp2 > sizeh) then
    sizeh = nhelp2
  end if

end do

!2.* def max{noa}, max{norb} ,max{nvb}

maxnoa = noa(1)
maxnvb = nvb(1)
maxnorb = norb(1)
do symp=1,nsym
  if (noa(symp) > maxnoa) maxnoa = noa(symp)
  if (norb(symp) > maxnorb) maxnorb = norb(symp)
  if (nvb(symp) > maxnvb) maxnvb = nvb(symp)
end do

!2.* def lengths of N fils

sizen = 0

do symp=1,nsym
  ! symq is not known for N file
  ! instead of norb(symr) maxnorb will be used so that reallength<=length
  sizen = sizen+norb(symp)*maxnorb
end do

!2.1 WX,VX - files

! wx%pos0 is defined
post = wx%pos0+sizew
vx%pos0 = post
post = post+sizew

!2.2 L - files

l1%pos0 = post
post = post+sizel
l2%pos0 = post
post = post+sizel

!2.3 RX - files

rx1%pos0 = post
post = post+sizer
rx2%pos0 = post
post = post+sizer
rx3%pos0 = post
post = post+sizer

!2.4 M - files

m1%pos0 = post
post = post+sizem
m2%pos0 = post
post = post+sizem
m3%pos0 = post
post = post+sizem

!2.5 H - files

h1%pos0 = post
post = post+sizeh
h2%pos0 = post
post = post+sizeh
h3%pos0 = post
post = post+sizeh

!2.6 N,PX - files

n%pos0 = post
post = post+sizen
px%pos0 = post
post = post+sizen

!2.7 declare space for help matrix D in for matrix multiplication C=AT*B if mchntyp=2

if (mchntyp == 2) then
  posd0 = post
  if (maxnoa <= maxnvb) then
    post = post+maxnoa*maxnoa*maxnvb*maxnvb
  else
    post = post+maxnoa*maxnoa*maxnoa*maxnoa
  end if
end if

length = post-1

return

end subroutine t3initfiles
