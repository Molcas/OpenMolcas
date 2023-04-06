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

subroutine intermezzo(wrk,wrksize,lunw3aaaa,lunw3bbbb,lunw3abba,lunw3baab,lunw3aabb,lunw3bbaa,lunt2o1,lunt2o2,lunt2o3,lunabij1, &
                      lunabij2,lunabij3)
! this routine calculates contributions
! WIII3, WIII4 and T26
!
! assigment of spin combinations:
!
! WIII(m,e,b,j)aaaa  - I
! WIII(m,e,b,j)bbbb  - K
! WIII(m,e,b,j)aabb  - L
! WIII(m,e,b,j)abba  - M
! WIII(m,e,b,j)baab  - N
! WIII(m,e,b,j)bbaa  - J
!
! WIII3
! WIII(m,e,b,j)aaaa <- sum(n-a) [ <mn||je>aaaa . T1o(n,b)aa ]
! WIII(m,e,b,j)bbbb <- sum(n-b) [ <mn||je>bbbb . T1o(n,b)bb ]
! WIII(m,e,b,j)aabb <- sum(n-b) [ <mn||je>abba . T1o(n,b)bb ]
! WIII(m,e,b,j)abba <- sum(n-b) [ <mn||je>abab . T1o(n,b)bb ]
! WIII(m,e,b,j)baab <- -sum(n-a) [ <nm||je>abba . T1o(n,b)aa ]
! WIII(m,e,b,j)bbaa <- -sum(n-a) [ <je||nm>abab . T1o(n,b)aa ]
!
! WIII4
! Q(f,b,j,n)aaaa   <= 0.5 T2o(f,b,j,n)aaaa + T1o(f,j)aa . T1o(b,n)aa
! WIII(m,e,b,j)aaaa <- -sum(n,f-aa)     [ Q(f,b,j,n)aaaa   . <ef||mn>aaaa ]
! <- 0.5 sum(n,f-bb)  [ T2o(b,f,j,n)abab . <ef||mn>abab ]
! WIII(m,e,b,j)aabb <- 0.5 sum(n,f-aa)  [ T2o(f,b,n,j)abab . <ef||mn>aaaa ]
! <- -sum(n,f-bb)     [ Q(f,b,j,n)bbbb   . <ef||mn>abab ]
! Q(f,b,j,n)bbbb   <= 0.5 T2o(f,b,j,n)bbbb + T1o(f,j)nn . T1o(b,n)bb
! WIII(m,e,b,j)bbbb <- -sum(n,f-bb)     [ Q(f,b,j,n)bbbb   . <ef||mn>bbbb ]
! <- 0.5 sum(n,f-aa)  [ T2o(f,b,n,j)abab . <fe||nm>abab ]
! WIII(m,e,b,j)bbaa <- 0.5 sum(n,f-bb)  [ T2o(b,f,j,n)abab . <ef||mn>bbbb ]
! <- - sum(n,f-aa)    [ Q(f,b,j,n)aaaa   . <fe||nm>abab ]
! Q(f,b,j,n)abab   <= 0.5 T2o(f,b,j,n)abab + T1o(f,j)aa . T1o(b,n)bb
! WIII(m,e,b,j)abba <- sum(n,f-ba)      [ Q(f,b,j,n)abab   . <fe||mn>abab ]
! Q(b,f,n,j)abab   <= 0.5 T2o(b,f,n,j)abab + T1o(b,n)aa . T1o(fmj)bb
! WIII(m,e,b,j)baab <- sum(n,f-ab)      [ Q(b,f,n,j)abab   . <ef||nm>abab ]
!
!
! T26
! R1(a,i,b,j)aaaa <= sum(m,e-aa) [ T2o(a,e,i,m)aaaa . WIII(m,e,b,j)aaaa ]
! <- sum(m,e-bb) [ T2o(a,e,i,m)abab . WIII(m,e,b,j)bbaa ]
! T2n(ab,ij)aaaa   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}aaaa
! R1(a,i,b,j)bbbb <= sum(m,e-bb) [ T2o(a,e,i,m)bbbb . WIII(m,e,b,j)bbbb ]
! <- sum(m,e-aa) [ T2o(e,a,m,i)abab . WIII(m,e,b,j)aabb ]
! T2n(ab,ij)bbbb   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}bbbb
! T2n(a,b,i,j)abab <- sum(m,e-aa) [ T2o(a,e,i,m)aaaa . WIII(m,e,b,j)aabb ]
! <- sum(m,e-aa) [ T2o(e,b,m,j)abab . WIII(m,e,a,i)aaaa ]
! <- sum(m,e-bb) [ T2o(a,e,i,m)abab . WIII(m,e,b,j)bbbb ]
! <- sum(m,e-bb) [ T2o(b,e,j,m)bbbb . WIII(m,e,a,i)bbaa ]
! <- sum(m,e-ab) [ T2o(a,e,m,j)abab . WIII(m,e,b,i)abba ]
! <- sum(m,e-ab) [ T2o(e,b,i,m)abab . WIII(m,e,a,j)baab ]
!
! N.B. use and destroy : V1,V2,V3,V4,M1,M2
! N.B. # of read      : 30 + 6
! # of write     : 2

use ccsd_global, only: idaaaa, idaabb, idabba, idbaab, idbbaa, idbbbb, m1, m2, t11, t12, t21, t22, t23, v1, v2, v3, v4, w11, w12, &
                       w13, w14
use Para_Info, only: MyRank
use Constants, only: One, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunw3aaaa, lunw3bbbb, lunw3abba, lunw3baab, lunw3aabb, lunw3bbaa, lunt2o1, lunt2o2, lunt2o3, &
                                   lunabij1, lunabij2, lunabij3
integer(kind=iwp) :: lunqaaaa, lunqbbbb, post, rc, ssc

!A.1 rewind nw3 files
!par
if (myRank == idaaaa) call filemanager(2,lunw3aaaa,rc)

if (myRank == idbaab) call filemanager(2,lunw3baab,rc)

if (myRank == idbbaa) call filemanager(2,lunw3bbaa,rc)

if (myRank == idbbbb) call filemanager(2,lunw3bbbb,rc)

if (myRank == idabba) call filemanager(2,lunw3abba,rc)

if (myRank == idaabb) call filemanager(2,lunw3aabb,rc)

!0.1 map M1(i,a) <- T1o(a,i)aa
call map(wrk,wrksize,2,2,1,0,0,t11,1,m1,post,rc)
!0.2 map M2(i,a) <- T1o(a,i)bb
call map(wrk,wrksize,2,2,1,0,0,t12,1,m2,post,rc)

! part I - W3aaaa
!par all contributions are calculated for idaaaa only,
!    just Qaaaa is calculated both on idaaaa and idbbaa

! par
if (idaaaa == myRank) then

  !I.1 get V1(m,e,a,j) <- W3aaaa(m,e,a,j)
  call getw3(wrk,wrksize,lunw3aaaa,1)

  !I.2 WIII(m,e,b,j)aaaa <- sum(n-a) [ <mn||je>aaaa . T1o(n,b)aa ]

  !I.2.1 expand V2(j,e,m,n) <- <je||mn>aaaa
  call expand(wrk,wrksize,4,3,w11,1,v2,rc)
  !I.2.2 map V3(m,e,j,n) <- V2(j,e,m,n)
  call map(wrk,wrksize,4,3,2,1,4,v2,1,v3,post,rc)
  !I.2.3 mult V2(m,e,j,b) <- V3(m,e,j,n) . M1(n,b)
  call ccmult(wrk,wrksize,4,2,4,1,v3,1,m1,1,v2,ssc,rc)
  !I.2.4 map V3(m,e,b,j) <- V2(m,e,j,b)
  call map(wrk,wrksize,4,1,2,4,3,v2,1,v3,post,rc)
  !I.2.5 add V1(m,e,b,j) <- 1.0 . V3(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,v1,1,rc)
end if
!parend

!I.3 WIII4
!    Q(f,b,j,n)aaaa   <= 0.5 T2o(f,b,j,n)aaaa + T1o(f,j)aa . T1o(b,n)aa
!    WIII(m,e,b,j)aaaa <- -sum(n,f-aa)     [ Q(f,b,j,n)aaaa   . <ef||mn>aaaa ]
!    <- 0.5 sum(n,f-bb)  [ T2o(b,f,j,n)abab . <ef||mn>abab ]

!par
if ((myRank == idaaaa) .or. (myRank == idbbaa)) then
  !I.3.1 get V3(fb,jn) <- T2o(fb,jn)aaaa
  call filemanager(2,lunt2o1,rc)
  call getmediate(wrk,wrksize,lunt2o1,v3,rc)
  !I.3.2 expand V2(f,b,j,n) <- V3(fb,jn)
  call expand(wrk,wrksize,4,4,v3,1,v2,rc)
  !I.3.3 mkQ V2(f,b,j,n) <- 0.5 V2(f,b,j,n) + T1o(f,j)aa . T1o(b,n)aa
  call mkq(wrk,wrksize,v2,t11,t11,Half,rc)
  !I.3.4 map V3(f,n,b,j) <- V2(f,b,j,n)
  call map(wrk,wrksize,4,1,3,4,2,v2,1,v3,post,rc)
end if
!parend

!par
if (myRank == idbbaa) then
  !I.3.5 write V3(f,n,b,j) to lunqaaaa
  call filemanager(1,lunqaaaa,rc)
  call wrtmediate(wrk,wrksize,lunqaaaa,v3,rc)
end if

!par
if (myRank == idaaaa) then
  !I.3.6 get V2(ef,mn) = <ef||mn>aaaa from luna file
  call filemanager(2,lunabij1,rc)
  call getmediate(wrk,wrksize,lunabij1,v2,rc)
  !I.3.7 expand V4(e,f,m,n) <- V2(ef,mn)
  call expand(wrk,wrksize,4,4,v2,1,v4,rc)
  !I.3.8 map V2(m,e,f,n) <- V4(e,f,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v4,1,v2,post,rc)
  !I.3.9 mult V4(m,e,b,j) <- V2(m,e,f,n) . V3(f,n,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v3,1,v4,ssc,rc)
  !I.3.10 add V1(m,e,b,j) <- -1.0 . V4(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,v4,1,v1,1,rc)
  !I.3.11 get V2(e,f,m,n) <- <ef||mn>abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v2,rc)
  !I.3.12 map V3(m,e,f,n) <- V2(e,f,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v2,1,v3,post,rc)
  !I.3.13 get V2(b,f,j,n) <- T2o(b,f,j,n)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v2,rc)
  !I.3.14 map V4(f,n,b,j) <- V2(b,f,j,n)
  call map(wrk,wrksize,4,3,1,4,2,v2,1,v4,post,rc)
  !I.3.15 mult V2(m,e,b,j) <- V3(m,e,f,n) . V4(f,n,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v3,1,v4,1,v2,ssc,rc)
  !I.3.16 add V1(m,e,b,j) <- 0.5 . V2(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,Half,v2,1,v1,1,rc)

  !I.4  R1(a,i,b,j)aaaa <= sum(m,e-aa) [ T2o(a,e,i,m)aaaa . WIII(m,e,b,j)aaaa ]
  !     T2n(ab,ij)aaaa   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}aaaa
  !I.4.1 get V2(ae,im) <- T2o(ae,im)aaaa
  call filemanager(2,lunt2o1,rc)
  call getmediate(wrk,wrksize,lunt2o1,v2,rc)
  !I.4.2 expand V3(a,e,i,m) <- V2(ae,im)
  call expand(wrk,wrksize,4,4,v2,1,v3,rc)
  !I.4.3 map V2(a,i,m,e) <- V3(a,e,i,m)
  call map(wrk,wrksize,4,1,4,2,3,v3,1,v2,post,rc)
  !I.4.4 mult V4(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v4,ssc,rc)
  !I.4.5 map V3(a,b,i,j) <- V4(a,i,b,j)
  call map(wrk,wrksize,4,1,3,2,4,v4,1,v3,post,rc)
  !I.4.6 pack V3(ab,ij) <- V2(ab,i,j) <- V3(a,b,i,j)
  call fack(wrk,wrksize,4,1,v3,1,v2,rc)
  call fack(wrk,wrksize,4,4,v2,1,v3,rc)
  !I.4.7 add T2n(ab,ij)aaaa <- 1.0 V3(ab,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t21,1,rc)

  !I.5 T2n(a,b,i,j)abab <- sum(m,e-aa) [ T2o(e,b,m,j)abab . WIII(m,e,a,i)aaaa ]
  !I.5.1 get V4(e,b,m,j) <- T2o(e,b,m,j)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v4,rc)
  !I.5.2 map V2(b,j,m,e) <- V4(e,b,m,j)
  call map(wrk,wrksize,4,4,1,3,2,v4,1,v2,post,rc)
  !I.5.3 mult V4(b,j,a,i) <- V2(b,j,m,e) . V1(m,e,a,i)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v4,ssc,rc)
  !I.5.4 map V3(a,b,i,j) <- V4(b,j,a,i)
  call map(wrk,wrksize,4,2,4,1,3,v4,1,v3,post,rc)
  !I.5.5 add T2n(a,b,i,j)abab <- 1.0 . V3(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t23,1,rc)
end if
!parend

! J part W3bbaa

! par
if (myRank == idbbaa) then

  !J.1 get V1(m,e,a,j) <- W3bbaa(m,e,a,j)
  call getw3(wrk,wrksize,lunw3bbaa,6)

  !J.2 WIII(m,e,b,j)bbaa <- -sum(n-a) [ <je||nm>abab . T1o(n,b)aa ]
  !J.2.1 map V3(m,e,j,n) <- <j,e||n,m>abab
  call map(wrk,wrksize,4,3,2,4,1,w13,1,v3,post,rc)
  !J.2.2 mult V4(m,e,j,b) <- V3(m,e,j,n) . M1(n,b)
  call ccmult(wrk,wrksize,4,2,4,1,v3,1,m1,1,v4,ssc,rc)
  !J.2.3 map V3(m,e,b,j) <- V4(m,e,j,b)
  call map(wrk,wrksize,4,1,2,4,3,v4,1,v3,post,rc)
  !J.2.4 add V1(m,e,b,j) <- -1.0 V3(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,v3,1,v1,1,rc)

  !J.2 WIII(m,e,b,j)bbaa <- 0.5 sum(n,f-bb)  [ T2o(b,f,j,n)abab . <ef||mn>bbbb ]
  !    <- - sum(n,f-aa)    [ Q(f,b,j,n)aaaa   . <fe||nm>abab ]
  !J.2.1 get V2(f,n,b,j) from lunqaaaa (produced in I step) and close it
  call filemanager(2,lunqaaaa,rc)
  call getmediate(wrk,wrksize,lunqaaaa,v2,rc)
  call filemanager(3,lunqaaaa,rc)
  !J.2.2 get V3(f,e,n,m) <- <fe||nm>abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v3,rc)
  !J.2.3 map V4(m,e,f,n) <- V3(f,e,n,m)
  call map(wrk,wrksize,4,3,2,4,1,v3,1,v4,post,rc)
  !J.2.4 mult V3(m,e,b,j) <- V4(m,e,f,n) . V2(f,n,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v4,1,v2,1,v3,ssc,rc)
  !J.2.5 add V1(m,e,b,j) <- -1.0 . V3(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,v3,1,v1,1,rc)
  !J.2.6 get V2(b,f,j,n) <- T2o(b,f,j,n)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v2,rc)
  !J.2.7 map V3(f,n,b,j) <- V2(b,f,j,n)
  call map(wrk,wrksize,4,3,1,4,2,v2,1,v3,post,rc)
  !J.3.8 get V2(ef,mn) = <ef||mn>bbbb from lunb file
  call filemanager(2,lunabij2,rc)
  call getmediate(wrk,wrksize,lunabij2,v2,rc)
  !J.3.9 exp V4(e,f,m,n) <- V2(ef,mn)
  call expand(wrk,wrksize,4,4,v2,1,v4,rc)
  !J.3.10 map V2(m,e,f,n) <- V4(e,f,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v4,1,v2,post,rc)
  !J.3.11 mult V4(m,e,b,j) <- V2(m,e,f,n) . V3(f,n,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v3,1,v4,ssc,rc)
  !J.3.12 add V1(m,e,b,j) <- 0.5 . V4(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,Half,v4,1,v1,1,rc)

  !J.4 R1(a,i,b,j)aaaa <= sum(m,e-bb) [ T2o(a,e,i,m)abab . WIII(m,e,b,j)bbaa ]
  !    T2n(ab,ij)aaaa   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}aaaa
  !J.4.1 get V3(a,e,i,m) <- T2o(a,e,i,m)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v3,rc)
  !J.4.2 map V2(a,i,m,e) <- V3(a,e,i,m)
  call map(wrk,wrksize,4,1,4,2,3,v3,1,v2,post,rc)
  !J.4.3 mult V4(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v4,ssc,rc)
  !J.4.4 map V3(a,b,i,j) <- V4(a,i,b,j)
  call map(wrk,wrksize,4,1,3,2,4,v4,1,v3,post,rc)
  !J.4.5 pack V3(ab,ij) <- V2(ab,i,j) <- V3(a,b,i,j)
  call fack(wrk,wrksize,4,1,v3,1,v2,rc)
  call fack(wrk,wrksize,4,4,v2,1,v3,rc)
  !J.4.6 add T2n(ab,ij)aaaa <- 1.0 V2(ab,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t21,1,rc)

  !J.5 T2n(a,b,i,j)abab <- sum(m,e-aa) [ T2o(b,e,j,m)bbbb . WIII(m,e,a,i)bbaa ]
  !J.5.1 get V4(be,jm) <- T2o(be,jm)bbbb
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,v4,rc)
  !J.5.2 expand V3(b,e,j,m) <-  V4(be,jm)
  call expand(wrk,wrksize,4,4,v4,1,v3,rc)
  !J.5.2 map V2(b,j,m,e) <- V3(b,e,j,m)
  call map(wrk,wrksize,4,1,4,2,3,v3,1,v2,post,rc)
  !J.5.3 mult V3(b,j,a,i) <- V2(b,j,m,e) . V1(m,e,a,i)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v3,ssc,rc)
  !J.5.4 map V2(a,b,i,j) <- V3(b,j,a,i)
  call map(wrk,wrksize,4,2,4,1,3,v3,1,v2,post,rc)
  !J.5.5 add T2n(a,b,i,j)abab <- 1.0 . V2(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v2,1,t23,1,rc)
end if
!parend

! part K - W3bbbb
!par all contributions are calculated for idbbbb only,
!    just Qbbbb is calculated both on idbbbb and idaabb

!par
if (myRank == idbbbb) then

  !K.1 get V1(m,e,a,j) <- W3bbbb(m,e,a,j)
  call getw3(wrk,wrksize,lunw3bbbb,2)

  !K.2 WIII(m,e,b,j)bbbb <- sum(n-b) [ <mn||je>bbbb . T1o(n,b)bb ]
  !K.2.1 expand V2(j,e,m,n) <- <je||mn>bbbb
  call expand(wrk,wrksize,4,3,w12,1,v2,rc)
  !K.2.2 map V3(m,e,j,n) <- V2(j,e,m,n)
  call map(wrk,wrksize,4,3,2,1,4,v2,1,v3,post,rc)
  !K.2.3 mult V2(m,e,j,b) <- V3(m,e,j,n) . M2(n,b)
  call ccmult(wrk,wrksize,4,2,4,1,v3,1,m2,1,v2,ssc,rc)
  !K.2.4 map V3(m,e,b,j) <- V2(m,e,j,b)
  call map(wrk,wrksize,4,1,2,4,3,v2,1,v3,post,rc)
  !K.2.5 add V1(m,e,b,j) <- 1.0 . V3(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,v1,1,rc)
end if
!parend

!K.3 WIII4
!    Q(f,b,j,n)bbbb   <= 0.5 T2o(f,b,j,n)bbbb + T1o(f,j)bb . T1o(b,n)bb
!    WIII(m,e,b,j)bbbb <- -sum(n,f-bb)     [ Q(f,b,j,n)bbbb   . <ef||mn>bbbb ]
!    <- 0.5 sum(n,f-aa)  [ T2o(f,b,n,j)abab . <fe||nm>abab ]

!par
if ((myRank == idbbbb) .or. (myRank == idaabb)) then
  !K.3.1 get V3(fb,jn) <- T2o(fb,jn)bbbb
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,v3,rc)
  !K.3.2 expand V2(f,b,j,n) <- V3(fb,jn)
  call expand(wrk,wrksize,4,4,v3,1,v2,rc)
  !K.3.3 mkQ V2(f,b,j,n) <- 0.5 V2(f,b,j,n) + T1o(f,j)bb . T1o(b,n)bb
  call mkq(wrk,wrksize,v2,t12,t12,Half,rc)
  !K.3.4 map V3(f,n,b,j) <- V2(f,b,j,n)
  call map(wrk,wrksize,4,1,3,4,2,v2,1,v3,post,rc)
end if
!parend

!par
if (myRank == idaabb) then
  !K.3.5 write V3(f,n,b,j) to lunqbbbb
  call filemanager(1,lunqbbbb,rc)
  call wrtmediate(wrk,wrksize,lunqbbbb,v3,rc)
end if

!par
if (myRank == idbbbb) then
  !K.3.6 get V2(ef,mn) = <ef||mn>bbbb from luna file
  call filemanager(2,lunabij2,rc)
  call getmediate(wrk,wrksize,lunabij2,v2,rc)
  !K.3.7 expand V4(e,f,m,n) <- V2(ef,mn)
  call expand(wrk,wrksize,4,4,v2,1,v4,rc)
  !K.3.8 map V2(m,e,f,n) <- V4(e,f,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v4,1,v2,post,rc)
  !K.3.9 mult V4(m,e,b,j) <- V2(m,e,f,n) . V3(f,n,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v3,1,v4,ssc,rc)
  !K.3.10 add V1(m,e,b,j) <- -1.0 . V4(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,v4,1,v1,1,rc)
  !K.3.11 get V2(f,e,n,m) <- <fe||nm>abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v2,rc)
  !K.3.12 map V3(m,e,f,n) <- V2(f,e,n,m)
  call map(wrk,wrksize,4,3,2,4,1,v2,1,v3,post,rc)
  !K.3.13 get V2(f,b,n,j) <- T2o(f,b,n,j)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v2,rc)
  !K.3.14 map V4(f,n,b,j) <- V2(f,b,n,j)
  call map(wrk,wrksize,4,1,3,2,4,v2,1,v4,post,rc)
  !K.3.15 mult V2(m,e,b,j) <- V3(m,e,f,n) . V4(f,n,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v3,1,v4,1,v2,ssc,rc)
  !K.3.16 add V1(m,e,b,j) <- 0.5 . V2(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,Half,v2,1,v1,1,rc)

  !K.4 R1(a,i,b,j)bbbb <= sum(m,e-bb) [ T2o(a,e,i,m)bbbb . WIII(m,e,b,j)bbbb ]
  !    T2n(ab,ij)bbbb   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}bbbb
  !K.4.1 get V2(ae,im) <- T2o(ae,im)bbbb
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,v2,rc)
  !K.4.2 expand V3(a,e,i,m) <- V2(ae,im)
  call expand(wrk,wrksize,4,4,v2,1,v3,rc)
  !K.4.3 map V2(a,i,m,e) <- V3(a,e,i,m)
  call map(wrk,wrksize,4,1,4,2,3,v3,1,v2,post,rc)
  !K.4.4 mult V4(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v4,ssc,rc)
  !K.4.5 map V3(a,b,i,j) <- V4(a,i,b,j)
  call map(wrk,wrksize,4,1,3,2,4,v4,1,v3,post,rc)
  !K.4.6 pack V3(ab,ij) <- V2(ab,i,j) <- V3(a,b,i,j)
  call fack(wrk,wrksize,4,1,v3,1,v2,rc)
  call fack(wrk,wrksize,4,4,v2,1,v3,rc)
  !K.4.7 add T2n(ab,ij)bbbb <- 1.0 V3(ab,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t22,1,rc)

  !K.5 T2n(a,b,i,j)abab <- sum(m,e-bb) [ T2o(a,e,i,m)abab . WIII(m,e,b,j)bbbb ]
  !K.5.1 get V4(a,e,i,m) <- T2o(a,e,i,m)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v4,rc)
  !K.5.2 map V2(a,i,m,e) <- V4(a,e,i,m)
  call map(wrk,wrksize,4,1,4,2,3,v4,1,v2,post,rc)
  !K.5.3 mult V4(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v4,ssc,rc)
  !K.5.4 map V3(a,b,i,j) <- V4(a,i,b,j)
  call map(wrk,wrksize,4,1,3,2,4,v4,1,v3,post,rc)
  !K.5.5 add T2n(a,b,i,j)abab <- 1.0 . V3(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t23,1,rc)
end if
!parend

! L part W3aabb

!par
if (myRank == idaabb) then

  !L.1 get V1(m,e,a,j) <- W3aabb(m,e,a,j)
  call getw3(wrk,wrksize,lunw3aabb,3)

  !L.2 WIII(m,e,b,j)aabb <- sum(n-b) [ <mn||je>abba . T1o(n,b)bb ]
  !L.2.1 map V3(m,e,j,n) <- <j,e||m,n>baab
  call map(wrk,wrksize,4,3,2,1,4,w14,1,v3,post,rc)
  !L.2.2 mult V4(m,e,j,b) <- V3(m,e,j,n) . M2(n,b)
  call ccmult(wrk,wrksize,4,2,4,1,v3,1,m2,1,v4,ssc,rc)
  !L.2.3 map V3(m,e,b,j) <- V4(m,e,j,b)
  call map(wrk,wrksize,4,1,2,4,3,v4,1,v3,post,rc)
  !L.2.4 add V1(m,e,b,j) <- 1.0 V3(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,v1,1,rc)

  !L.2 WIII(m,e,b,j)aabb <- 0.5 sum(n,f-aa)  [ T2o(f,b,n,j)abab . <ef||mn>aaaa ]
  !    <- - sum(n,f-bb)    [ Q(f,b,j,n)bbbb   . <ef||mn>abab ]
  !L.2.1 get V2(f,n,b,j) from lunqbbbb (produced in K step) and close it
  call filemanager(2,lunqbbbb,rc)
  call getmediate(wrk,wrksize,lunqbbbb,v2,rc)
  call filemanager(3,lunqbbbb,rc)
  !L.2.2 get V3(e,f,m,n) <- <ef||mn>abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v3,rc)
  !L.2.3 map V4(m,e,f,n) <- V3(e,f,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v3,1,v4,post,rc)
  !L.2.4 mult V3(m,e,b,j) <- V4(m,e,f,n) . V2(f,n,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v4,1,v2,1,v3,ssc,rc)
  !L.2.5 add V1(m,e,b,j) <- -1.0 . V3(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,v3,1,v1,1,rc)
  !L.2.6 get V2(f,b,n,j) <- T2o(f,b,n,j)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v2,rc)
  !L.2.7 map V3(f,n,b,j) <- V2(f,b,n,j)
  call map(wrk,wrksize,4,1,3,2,4,v2,1,v3,post,rc)
  !L.3.8 get V2(ef,mn) = <ef||mn>aaaa from lunb file
  call filemanager(2,lunabij1,rc)
  call getmediate(wrk,wrksize,lunabij1,v2,rc)
  !L.3.9 exp V4(e,f,m,n) <- V2(ef,mn)
  call expand(wrk,wrksize,4,4,v2,1,v4,rc)
  !L.3.10 map V2(m,e,f,n) <- V4(e,f,m,n)
  call map(wrk,wrksize,4,2,3,1,4,v4,1,v2,post,rc)
  !L.3.11 mult V4(m,e,b,j) <- V2(m,e,f,n) . V3(f,n,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v3,1,v4,ssc,rc)
  !L.3.12 add V1(m,e,b,j) <- 0.5 . V4(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,Half,v4,1,v1,1,rc)

  !L.4 R1(a,i,b,j)bbbb <- sum(m,e-aa) [ T2o(e,a,m,i)abab . WIII(m,e,b,j)aabb ]
  !    T2n(ab,ij)bbbb   <= {1(a,i,b,j)-R1(b,i,a,j)-R1(a,j,b,i)+R1(b,j,a,i)}bbbb
  !L.4.1 get V3(e,a,m,i) <- T2o(e,a,m,i)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v3,rc)
  !L.4.2 map V2(a,i,m,e) <- V3(e,a,m,i)
  call map(wrk,wrksize,4,4,1,3,2,v3,1,v2,post,rc)
  !L.4.3 mult V4(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v4,ssc,rc)
  !L.4.4 map V3(a,b,i,j) <- V4(a,i,b,j)
  call map(wrk,wrksize,4,1,3,2,4,v4,1,v3,post,rc)
  !L.4.5 pack V3(ab,ij) <- V2(ab,i,j) <- V3(a,b,i,j)
  call fack(wrk,wrksize,4,1,v3,1,v2,rc)
  call fack(wrk,wrksize,4,4,v2,1,v3,rc)
  !L.4.6 add T2n(ab,ij)aaaa <- 1.0 V2(ab,ij)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t22,1,rc)

  !L.5 T2n(a,b,i,j)abab <- sum(m,e-aa) [ T2o(a,e,i,m)aaaa . WIII(m,e,b,j)aabb ]
  !L.5.1 get V4(ae,im) <- T2o(ae,im)aaaa
  call filemanager(2,lunt2o1,rc)
  call getmediate(wrk,wrksize,lunt2o1,v4,rc)
  !L.5.2 expand V3(a,e,i,m) <-  V4(ae,im)
  call expand(wrk,wrksize,4,4,v4,1,v3,rc)
  !L.5.2 map V2(a,i,m,e) <- V3(a,e,i,m)
  call map(wrk,wrksize,4,1,4,2,3,v3,1,v2,post,rc)
  !L.5.3 mult V3(a,i,b,j) <- V2(a,i,m,e) . V1(m,e,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v2,1,v1,1,v3,ssc,rc)
  !L.5.4 map V2(a,b,i,j) <- V3(a,i,b,j)
  call map(wrk,wrksize,4,1,3,2,4,v3,1,v2,post,rc)
  !L.5.5 add T2n(a,b,i,j)abab <- 1.0 . V2(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v2,1,t23,1,rc)
end if
!parend

! M part W3abba

!par
if (myRank == idabba) then

  !M.1 get V1(m,e,a,j) <- W3abba
  call getw3(wrk,wrksize,lunw3abba,4)

  !M.2 WIII(m,e,b,j)abba <- sum(n-b) [ <mn||je>abab . T1o(n,b)bb ]
  !M.2.1 map V3(m,e,j,n) <- <j,e||m,n>abab
  call map(wrk,wrksize,4,3,2,1,4,w13,1,v3,post,rc)
  !M.2.2 mult V4(m,e,j,b) <- V3(m,e,j,n) . M2(n,b)
  call ccmult(wrk,wrksize,4,2,4,1,v3,1,m2,1,v4,ssc,rc)
  !M.2.3 map V3(m,e,b,j) <- V4(m,e,j,b)
  call map(wrk,wrksize,4,1,2,4,3,v4,1,v3,post,rc)
  !M.2.4 add V1(m,e,b,j) <- 1.0 V3(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,v1,1,rc)

  !M.3 Q(f,b,j,n)abab   <= 0.5 T2o(f,b,j,n)abab + T1o(f,j)aa . T1o(b,n)bb
  !    WIII(m,e,b,j)abba <- sum(n,f-ba)      [ Q(f,b,j,n)abab   . <fe||mn>abab ]
  !M.3.1 get V4(f,b,j,n) <- T2o(f,b,j,n)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v4,rc)
  !M.3.2 mkQ V4(f,b,j,n) <- 0.5 V4(f,b,j,n) + T1o(f,j)aa . T1o(b,n)bb
  call mkq(wrk,wrksize,v4,t11,t12,Half,rc)
  !M.3.3 map V3(f,n,b,j) <- V4(f,b,j,n)
  call map(wrk,wrksize,4,1,3,4,2,v4,1,v3,post,rc)
  !M.3.4 get V2(f,e,m,n) <- <fe||mn>abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v2,rc)
  !M.3.5 map V4(m,e,f,n) <- V2(f,e,m,n)
  call map(wrk,wrksize,4,3,2,1,4,v2,1,v4,post,rc)
  !M.3.6 mult V2(m,e,b,j) <- V4(m,e,f,n) . V3(f,n,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v4,1,v3,1,v2,ssc,rc)
  !M.3.7 add V1(m,e,b,j) <- 1.0 V2(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v2,1,v1,1,rc)

  !M.4 T2n(a,b,i,j)abab <- sum(m,e-ab) [ T2o(a,e,m,j)abab . WIII(m,e,b,i)abba ]
  !M.4.1 get V2(a,e,m,j) <- T2o(a,e,m,j)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v2,rc)
  !M.4.2 map V3(a,j,m,e) <- V2(a,e,m,j)
  call map(wrk,wrksize,4,1,4,3,2,v2,1,v3,post,rc)
  !M.4.3 mult V4(a,j,b,i) <- V3(a,j,m,e) . V1(m,e,b,i)
  call ccmult(wrk,wrksize,4,4,4,2,v3,1,v1,1,v4,ssc,rc)
  !M.4.4 map V3(a,b,i,j) <- V4(a,j,b,i)
  call map(wrk,wrksize,4,1,4,2,3,v4,1,v3,post,rc)
  !M.4.5 add T2n(a,b,i,j)abab <- 1.0 V3(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t23,1,rc)
end if
!parend

! N part W3baab

!par
if (myRank == idbaab) then

  !N.1 get V1(m,e,a,j) <- W3baab
  call getw3(wrk,wrksize,lunw3baab,5)

  !N.2 WIII(m,e,b,j)baab <- - sum(n-a) [ <je||nm>baab . T1o(n,b)aa ]
  !N.2.1 map V3(m,e,j,n) <- <j,e||n,m>baab
  call map(wrk,wrksize,4,3,2,4,1,w14,1,v3,post,rc)
  !N.2.2 mult V4(m,e,j,b) <- V3(m,e,j,n) . M1(n,b)
  call ccmult(wrk,wrksize,4,2,4,1,v3,1,m1,1,v4,ssc,rc)
  !N.2.3 map V3(m,e,b,j) <- V4(m,e,j,b)
  call map(wrk,wrksize,4,1,2,4,3,v4,1,v3,post,rc)
  !N.2.4 add V1(m,e,b,j) <- -1.0 V3(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,-One,v3,1,v1,1,rc)

  !N.3  Q(b,f,n,j)abab   <= 0.5 T2o(b,f,n,j)abab + T1o(b,n)aa . T1o(f,j)bb
  !     WIII(m,e,b,j)baab <- sum(n,f-ab)      [ Q(b,f,n,j)abab   . <ef||nm>abab ]
  !N.3.1 get V4(b,f,n,j) <- T2o(b,f,n,j)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v4,rc)
  !N.3.2 mkQ V4(b,f,n,j) <- 0.5 V4(b,f,n,j) + T1o(b,n)aa . T1o(f,j)bb
  call mkq(wrk,wrksize,v4,t11,t12,Half,rc)
  !N.3.3 map V3(f,n,b,j) <- V4(b,f,n,j)
  call map(wrk,wrksize,4,3,1,2,4,v4,1,v3,post,rc)
  !N.3.4 get V2(e,f,n,m) <- <ef||nm>abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,v2,rc)
  !N.3.5 map V4(m,e,f,n) <- V2(e,f,n,m)
  call map(wrk,wrksize,4,2,3,4,1,v2,1,v4,post,rc)
  !N.3.6 mult V2(m,e,b,j) <- V4(m,e,f,n) . V3(f,n,b,j)
  call ccmult(wrk,wrksize,4,4,4,2,v4,1,v3,1,v2,ssc,rc)
  !N.3.7 add V1(m,e,b,j) <- 1.0 V2(m,e,b,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v2,1,v1,1,rc)

  !N.4 T2n(a,b,i,j)abab <- sum(m,e-ab) [ T2o(e,b,i,m)abab . WIII(m,e,a,j)baab ]
  !N.4.1 get V2(e,b,i,m) <- T2o(e,b,i,m)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v2,rc)
  !N.4.2 map V3(b,i,m,e) <- V2(e,b,i,m)
  call map(wrk,wrksize,4,4,1,2,3,v2,1,v3,post,rc)
  !N.4.3 mult V4(b,i,a,j) <- V3(b,i,m,e) . V1(m,e,a,j)
  call ccmult(wrk,wrksize,4,4,4,2,v3,1,v1,1,v4,ssc,rc)
  !N.4.4 map V3(a,b,i,j) <- V4(b,i,a,j)
  call map(wrk,wrksize,4,2,3,1,4,v4,1,v3,post,rc)
  !N.4.5 add T2n(a,b,i,j)abab <- 1.0 V3(a,b,i,j)
  call add(wrk,wrksize,4,4,0,0,0,0,1,1,One,v3,1,t23,1,rc)
end if
!parend

return

end subroutine intermezzo
