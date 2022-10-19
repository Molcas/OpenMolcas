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

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: length
#include "t31.fh"
#include "t32.fh"
integer(kind=iwp) :: maxnoa, maxnorb, maxnvb, nhelp1, nhelp2, post, sizeh, sizel, sizem, sizen, sizer, sizew, symp, symq, symr

!1 maps and positions for fix mediated

!1.0 maps for DP - diagonal part
!    N.B. DP has one degree of freedom, while other 1 index has none
!    DP1 - dp(p)a
!    DP2 - dp(p)b

do symp=1,nsym
  do symq=1,nsym
    do symr=1,nsym
      mapidp1(symp,symq,symr) = 0
      mapidp2(symp,symq,symr) = 0
    end do
  end do
end do

post = 1

posdp10 = post
mapddp1(0,1) = 5
mapddp1(0,2) = 0
mapddp1(0,3) = 0
mapddp1(0,4) = 0
mapddp1(0,5) = nsym
mapddp1(0,6) = 0

do symp=1,nsym
  mapddp1(symp,1) = post
  mapddp1(symp,2) = norb(symp)
  mapddp1(symp,3) = symp
  mapddp1(symp,4) = 1
  mapddp1(symp,5) = 1
  mapddp1(symp,6) = 1
  mapidp1(symp,1,1) = symp
  post = post+norb(symp)
end do

posdp20 = post
mapddp2(0,1) = 5
mapddp2(0,2) = 0
mapddp2(0,3) = 0
mapddp2(0,4) = 0
mapddp2(0,5) = nsym
mapddp2(0,6) = 0

do symp=1,nsym
  mapddp2(symp,1) = post
  mapddp2(symp,2) = norb(symp)
  mapddp2(symp,3) = symp
  mapddp2(symp,4) = 1
  mapddp2(symp,5) = 1
  mapddp2(symp,6) = 1
  mapidp2(symp,1,1) = symp
  post = post+norb(symp)
end do

!1.1 maps for T1
!    T11 - t1oaa(a,i)
!    T12 - t1obb(a,i)

post110 = post
call cct3_grc0(2,0,3,1,0,0,1,post110,post,mapdt11,mapit11)
post120 = post
call cct3_grc0(2,0,4,2,0,0,1,post120,post,mapdt12,mapit12)

!1.5 maps for FK
!    FK1 - f(a,b)aa
!    FK2 - f(a,b)bb
!    FK3 - f(a,i)aa
!    FK4 - f(a,i)bb
!    FK5 - f(i,j)aa
!    FK6 - f(i,j)bb

posfk10 = post
call cct3_grc0(2,0,3,3,0,0,1,posfk10,post,mapdfk1,mapifk1)
posfk20 = post
call cct3_grc0(2,0,4,4,0,0,1,posfk20,post,mapdfk2,mapifk2)
posfk30 = post
call cct3_grc0(2,0,3,1,0,0,1,posfk30,post,mapdfk3,mapifk3)
posfk40 = post
call cct3_grc0(2,0,4,2,0,0,1,posfk40,post,mapdfk4,mapifk4)
posfk50 = post
call cct3_grc0(2,0,1,1,0,0,1,posfk50,post,mapdfk5,mapifk5)
posfk60 = post
call cct3_grc0(2,0,2,2,0,0,1,posfk60,post,mapdfk6,mapifk6)

!1.6 maps for T2
!    T21 - t2o(ab,ij)aaaa
!    T22 - t2o(ab,ij)bbbb
!    T23 - t2o(a,b,i,j)abab

post210 = post
call cct3_grc0(4,4,3,3,1,1,1,post210,post,mapdt21,mapit21)
post220 = post
call cct3_grc0(4,4,4,4,2,2,1,post220,post,mapdt22,mapit22)
post230 = post
call cct3_grc0(4,0,3,4,1,2,1,post230,post,mapdt23,mapit23)

!1.8 maps for W1
!    W11 - <ie||mn>aaaa
!    W12 - <ie||mn>bbbb
!    W13 - <ie||mn>abab
!    W14 - <ie||mn>baab

posw110 = post
call cct3_grc0(4,3,1,3,1,1,1,posw110,post,mapdw11,mapiw11)
posw120 = post
call cct3_grc0(4,3,2,4,2,2,1,posw120,post,mapdw12,mapiw12)
posw130 = post
call cct3_grc0(4,0,1,4,1,2,1,posw130,post,mapdw13,mapiw13)
posw140 = post
call cct3_grc0(4,0,2,3,1,2,1,posw140,post,mapdw14,mapiw14)

!1.9 maps for W2
!    W21 - <ab||ij>aaaa
!    W22 - <ab||ij>bbbb
!    W23 - <a,b|i,j>abab

posw210 = post
call cct3_grc0(4,4,3,3,1,1,1,posw210,post,mapdw21,mapiw21)
posw220 = post
call cct3_grc0(4,4,4,4,2,2,1,posw220,post,mapdw22,mapiw22)
posw230 = post
call cct3_grc0(4,0,3,4,1,2,1,posw230,post,mapdw23,mapiw23)

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

posw0 = post

!2.* find maxsize of W,L,M,H
sizew = 0
sizel = 0
sizer = 0
sizem = 0
sizeh = 0

do nhelp1=1,nsym

  ! W,V files
  call cct3_t3grc0(3,2,4,4,4,0,nhelp1,posw0,post,mapdw,mapiw)
  nhelp2 = post-posw0
  if (nhelp2 > sizew) then
    sizew = nhelp2
  end if

  ! L files
  call cct3_t3grc0(3,0,4,4,4,0,nhelp1,posw0,post,mapdw,mapiw)
  nhelp2 = post-posw0
  if (nhelp2 > sizel) then
    sizel = nhelp2
  end if
  call cct3_t3grc0(3,0,1,4,4,0,nhelp1,posw0,post,mapdw,mapiw)
  nhelp2 = post-posw0
  if (nhelp2 > sizel) then
    sizel = nhelp2
  end if

  ! R files
  call cct3_t3grc0(3,8,4,4,4,0,nhelp1,posw0,post,mapdw,mapiw)
  nhelp2 = post-posw0
  if (nhelp2 > sizer) then
    sizer = nhelp2
  end if

  ! M files
  call cct3_t3grc0(2,0,4,4,0,0,nhelp1,posw0,post,mapdw,mapiw)
  nhelp2 = post-posw0
  if (nhelp2 > sizem) then
    sizem = nhelp2
  end if
  call cct3_t3grc0(2,0,1,4,0,0,nhelp1,posw0,post,mapdw,mapiw)
  nhelp2 = post-posw0
  if (nhelp2 > sizem) then
    sizem = nhelp2
  end if

  ! H files
  call cct3_t3grc0(1,0,4,0,0,0,nhelp1,posw0,post,mapdw,mapiw)
  nhelp2 = post-posw0
  if (nhelp2 > sizeh) then
    sizeh = nhelp2
  end if
  call cct3_t3grc0(1,0,1,0,0,0,nhelp1,posw0,post,mapdw,mapiw)
  nhelp2 = post-posw0
  if (nhelp2 > sizeh) then
    sizeh = nhelp2
  end if

end do

!2.* def max{noa}, max{norb} ,max{nvb}

maxnoa = noa(1)
maxnvb = nvb(1)
maxnorb = norb(1)
do symp=1,nsym
  if (noa(symp) > maxnoa) then
    maxnoa = noa(symp)
  end if
  if (norb(symp) > maxnorb) then
    maxnorb = norb(symp)
  end if
  if (nvb(symp) > maxnvb) then
    maxnvb = nvb(symp)
  end if
end do

!2.* def lengths of N fils

sizen = 0

do symp=1,nsym
  ! symq is not known for N file
  ! instead of norb(symr) maxnorb will be used so that reallength<=length
  sizen = sizen+norb(symp)*maxnorb
end do

!2.1 W,V - files

! posw0 is defined
post = posw0+sizew
posv0 = post
post = post+sizew

!2.2 L - files

posl10 = post
post = post+sizel
posl20 = post
post = post+sizel

!2.3 R - files

posr10 = post
post = post+sizer
posr20 = post
post = post+sizer
posr30 = post
post = post+sizer

!2.4 M - files

posm10 = post
post = post+sizem
posm20 = post
post = post+sizem
posm30 = post
post = post+sizem

!2.5 H - files

posh10 = post
post = post+sizeh
posh20 = post
post = post+sizeh
posh30 = post
post = post+sizeh

!2.6 N,P - files

posn0 = post
post = post+sizen
posp0 = post
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
