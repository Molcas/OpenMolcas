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

subroutine t3sgl(wrk,wrksize,mapdw,mapds1,mapis1,mapds2,mapis2,mapdd1,mapid1,mapdd2,mapid2,typdiv,i,j,k,symi,symj,symk,rc1,mapdm1, &
                 mapim1,posm10,mapdh1,mapih1,posh10,mapdm2,mapim2,posm20,mapdh2,mapih2,posh20,mapdm3,mapim3,posm30,mapdh3,mapih3, &
                 posh30)
! mapdw          - direct map matrix of W (Input)
! mapds1         - direct map matrix of S1 (Input)
! mapis1         - inverse map matrix of S1 (Input)
! mapds2         - direct map matrix of S2 (Input)
! mapis2         - inverse map matrix of S2 (Input)
! mapdd1         - direct map matrix of D1 (Input)
! mapid1         - inverse map matrix of D1 (Input)
! mapdd2         - direct map matrix of D2 (Input)
! mapid2         - inverse map matrix of D2 (Input)
!                  (order is aa>ab>bb, if there is only one spin, use any map's for 2)
! typsgl         - type of operation (see Table) (Input)
! i              - value of occupied index i (Inlut)
! j              - value of occupied index j (Inlut)
! k              - value of occupied index k (Inlut)
! symi           - symmetry of index i (Input)
! symj           - symmetry of index j (Input)
! symk           - symmetry of index k (Input)
! rc1            - return (error) code (Output)
! mapd,mapi,poss - parameters for M1-3,H1-3 files (I)
!
! this routine adds contributions from disconnected
! singles, namely:
!
! W_ijk(abc) <- P [ U1(_i,a) . U2(_jk,bc) ]
!
! for following types of W
!
! typsgl         Operation                   Implemented
! 1     W(abc) <- +s1(_i,a) . d1(_jk,bc)
!                 -s1(_i,b) . d1(_jk,ac)
!                 +s1(_i,c) . d1(_jk,ab)
!
!                 -s1(_j,a) . d1(_ik,bc)
!                 +s1(_j,b) . d1(_ik,ac)
!                 -s1(_j,c) . d1(_ik,ab)
!
!                 +s1(_k,a) . d1(_ij,bc)
!                 -s1(_k,b) . d1(_ij,ac)
!                 +s1(_k,c) . d1(_ij,ab)        Yes
!
! 2     W(ab,c)<- +s1(_i,a) . d2(_j,_k,b,c)
!                 -s1(_i,b) . d2(_j,_k,a,c)
!
!                 -s1(_j,a) . d2(_i,_k,b,c)
!                 +s1(_j,b) . d2(_i,_k,a,c)
!
!                 +s2(_k,c) . d1(_ij,ab)        Yes
!
! 3     W(a,bc)<- +s1(_i,a) . d2(_jk,b,c)
!
!                 +s1(_j,b) . d2(_i,_k,a,c)
!                 -s1(_j,c) . d2(_i,_k,a,b)
!
!                 -s1(_k,b) . d2(_i,_j,a,c)
!                 +s1(_k,c) . d2(_i,_j,a,b)     Yes
!
! N.B. spin combinations aaa,bbb for 1; aab for 2; and abb for 3
! are automatically assumed

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, mapdw(0:512,6), mapds1(0:512,6), mapis1(8,8,8), mapds2(0:512,6), mapis2(8,8,8), mapdd1(0:512,6), &
                     mapid1(8,8,8), mapdd2(0:512,6), mapid2(8,8,8), typdiv, i, j, k, symi, symj, symk, rc1, mapdm1(0:512,6), &
                     mapim1(8,8,8), posm10, mapdh1(0:512,6), mapih1(8,8,8), posh10, mapdm2(0:512,6), mapim2(8,8,8), posm20, &
                     mapdh2(0:512,6), mapih2(8,8,8), posh20, mapdm3(0:512,6), mapim3(8,8,8), posm30, mapdh3(0:512,6), &
                     mapih3(8,8,8), posh30
real(kind=wp) :: wrk(wrksize)
#include "t31.fh"
integer(kind=iwp) :: dima, dimb, dimc, id1, id2, id3, is1, is2, is3, iw, nhelp1, nhelp2, posd1, posd2, posd3, poss1, poss2, poss3, &
                     posw, ssh1, ssh2, ssh3, ssm1, ssm2, ssm3, syma, symab, symac, symb, symbc, symc, symij, symik, symjk

!0.* some tests

if (typdiv == 1) then

  !1 case W(pqr)

  !1.* ext H1(a) <= S1(a,i) for given i
  call ext(wrk,wrksize,2,2,i,0,symi,0,0,mapds1,mapis1,1,posh10,mapdh1,mapih1,ssh1,rc1)

  !1.* ext H2(a) <= S1(a,j) for given j
  call ext(wrk,wrksize,2,2,j,0,symj,0,0,mapds1,mapis1,1,posh20,mapdh2,mapih2,ssh2,rc1)

  !1.* ext H3(a) <= S1(a,k) for given k
  call ext(wrk,wrksize,2,2,k,0,symk,0,0,mapds1,mapis1,1,posh30,mapdh3,mapih3,ssh3,rc1)

  !1.* ext M1(bc) <= D1(bc,jk) for given jk
  call ext(wrk,wrksize,4,7,j,k,symj,symk,0,mapdd1,mapid1,1,posm10,mapdm1,mapim1,ssm1,rc1)

  !1.* ext M2(bc) <= D1(bc,ik) for given ik
  call ext(wrk,wrksize,4,7,i,k,symi,symk,0,mapdd1,mapid1,1,posm20,mapdm2,mapim2,ssm2,rc1)

  !1.* ext M3(bc) <= D1(bc,ij) for given ij
  call ext(wrk,wrksize,4,7,i,j,symi,symj,0,mapdd1,mapid1,1,posm30,mapdm3,mapim3,ssm3,rc1)

  do iw=1,mapdw(0,5)

    !1.* def position of W
    posw = mapdw(iw,1)

    !1.* def symmetry status
    syma = mapdw(iw,3)
    symb = mapdw(iw,4)
    symc = mapdw(iw,5)

    !1.* def dimensions
    dima = dimm(mapdw(0,1),syma)
    dimb = dimm(mapdw(0,2),symb)
    dimc = dimm(mapdw(0,3),symc)
    !1.*
    symij = mmul(symi,symj)
    symik = mmul(symi,symk)
    symjk = mmul(symj,symk)
    symab = mmul(syma,symb)
    symac = mmul(syma,symc)
    symbc = mmul(symb,symc)

    !1.* realize packing

    if (syma == symc) then
      !1.a case syma=symb=symc

      !1.a.1---- 1st. triade   W(abc) <- +s1(_i,a) . d1(_jk,bc)
      !                                  -s1(_i,b) . d1(_jk,ac)
      !                                  +s1(_i,c) . d1(_jk,ab)

      if ((symi == syma) .and. (symjk == symbc)) then
        !1.a.1.* find address for s1,d1
        is1 = mapih1(1,1,1)
        id1 = mapim1(symb,1,1)

        !1.a.1.* def position of s1,d1
        poss1 = mapdh1(is1,1)
        posd1 = mapdm1(id1,1)

        !1.a.1.* def additional dimensions
        nhelp1 = dima*(dima-1)/2
        nhelp2 = dima*(dima-1)*(dima-2)/6

        !1.a.1.* add singly
        call t3sglh11(wrk(posw),dima,nhelp1,nhelp2,wrk(poss1),wrk(posd1),1)
      end if

      !1.a.2---- 2nd. triade   W(abc) <- -s1(_j,a) . d1(_ik,bc)
      !                                  +s1(_j,b) . d1(_ik,ac)
      !                                  -s1(_j,c) . d1(_ik,ab)

      if ((symj == syma) .and. (symik == symbc)) then
        !1.a.2.* find address for s1,d1
        is1 = mapih2(1,1,1)
        id1 = mapim2(symb,1,1)

        !1.a.2.* def position of s1,d1
        poss1 = mapdh2(is1,1)
        posd1 = mapdm2(id1,1)

        !1.a.2.* def additional dimensions
        nhelp1 = dima*(dima-1)/2
        nhelp2 = dima*(dima-1)*(dima-2)/6

        !1.a.2.* add singly
        call t3sglh11(wrk(posw),dima,nhelp1,nhelp2,wrk(poss1),wrk(posd1),-1)
      end if

      !1.a.3---- 3rd. triade   W(abc) <- +s1(_k,a) . d1(_ij,bc)
      !                                  -s1(_k,b) . d1(_ij,ac)
      !                                  +s1(_k,c) . d1(_ij,ab)

      if ((symk == syma) .and. (symij == symbc)) then

        !1.a.3.* find address for s1,d1
        is1 = mapih3(1,1,1)
        id1 = mapim3(symb,1,1)

        !1.a.3.* def position of s1,d1
        poss1 = mapdh3(is1,1)
        posd1 = mapdm3(id1,1)

        !1.a.3.* def additional dimensions
        nhelp1 = dima*(dima-1)/2
        nhelp2 = dima*(dima-1)*(dima-2)/6

        !1.a.3.* add singly
        call t3sglh11(wrk(posw),dima,nhelp1,nhelp2,wrk(poss1),wrk(posd1),1)
      end if

    else if (syma == symb) then
      !1.b case syma=symb/=symc

      !1.b.1---- 1st. triade   W(abc) <- +s1(_i,a) . d1(_jk,bc)
      !                                  -s1(_i,b) . d1(_jk,ac)
      !                                  +s1(_i,c) . d1(_jk,ab)

      if ((symi == syma) .and. (symjk == symbc)) then
        ! if syma=symi, then obviously symc/=symi
        !1.b.1.* find address for s1,d1
        is1 = mapih1(1,1,1)
        id1 = mapim1(syma,1,1)

        !1.b.1.* def position of s1,d1
        poss1 = mapdh1(is1,1)
        posd1 = mapdm1(id1,1)

        !1.b.1.* def additional dimensions
        nhelp1 = dima*(dima-1)/2

        !1.b.1.* add singly
        call t3sglh121(wrk(posw),dima,nhelp1,dimc,wrk(poss1),wrk(posd1),1)
      else if ((symc == symi) .and. (symjk == symab)) then
        !1.b.1.* find address for s3,d3
        is3 = mapih1(1,1,1)
        id3 = mapim1(syma,1,1)

        !1.b.1.* def position of s1,d3
        poss3 = mapdh1(is3,1)
        posd3 = mapdm1(id3,1)

        !1.b.1.* def additional dimensions
        nhelp1 = dima*(dima-1)/2

        !1.b.1.* add singly
        call t3sglh122(wrk(posw),nhelp1,dimc,wrk(poss3),wrk(posd3),1)
      end if

      !1.b.2---- 2nd. triade   W(abc) <- -s1(_j,a) . d1(_ik,bc)
      !                                  +s1(_j,b) . d1(_ik,ac)
      !                                  -s1(_j,c) . d1(_ik,ab)

      if ((symj == syma) .and. (symik == symbc)) then
        ! if syma=symj, then obviously symc/=symj
        !1.b.2.* find address for s1,d1
        is1 = mapih2(1,1,1)
        id1 = mapim2(syma,1,1)

        !1.b.2.* def position of s1,d1
        poss1 = mapdh2(is1,1)
        posd1 = mapdm2(id1,1)

        !1.b.2.* def additional dimensions
        nhelp1 = dima*(dima-1)/2

        !1.b.2.* add singly
        call t3sglh121(wrk(posw),dima,nhelp1,dimc,wrk(poss1),wrk(posd1),-1)
      else if ((symc == symj) .and. (symik == symab)) then
        !1.b.2.* find address for s3,d3
        is3 = mapih2(1,1,1)
        id3 = mapim2(syma,1,1)

        !1.b.2.* def position of s1,d3
        poss3 = mapdh2(is3,1)
        posd3 = mapdm2(id3,1)

        !1.b.2.* def additional dimensions
        nhelp1 = dima*(dima-1)/2

        !1.b.2.* add singly
        call t3sglh122(wrk(posw),nhelp1,dimc,wrk(poss3),wrk(posd3),-1)
      end if

      !1.b.3---- 3rd. triade   W(abc) <- +s1(_k,a) . d1(_ij,bc)
      !                                  -s1(_k,b) . d1(_ij,ac)
      !                                  +s1(_k,c) . d1(_ij,ab)

      if ((symk == syma) .and. (symij == symbc)) then
        ! if syma=symk, then obviously symc/=symk
        !1.b.3.* find address for s1,d1
        is1 = mapih3(1,1,1)
        id1 = mapim3(syma,1,1)

        !1.b.3.* def position of s1,d1
        poss1 = mapdh3(is1,1)
        posd1 = mapdm3(id1,1)

        !1.b.3.* def additional dimensions
        nhelp1 = dima*(dima-1)/2

        !1.b.3.* add singly
        call t3sglh121(wrk(posw),dima,nhelp1,dimc,wrk(poss1),wrk(posd1),1)
      else if ((symc == symk) .and. (symij == symab)) then
        !1.b.3.* find address for s3,d3
        is3 = mapih3(1,1,1)
        id3 = mapim3(syma,1,1)

        !1.b.3.* def position of s1,d3
        poss3 = mapdh3(is3,1)
        posd3 = mapdm3(id3,1)

        !1.b.3.* def additional dimensions
        nhelp1 = dima*(dima-1)/2

        !1.b.3.* add singly
        call t3sglh122(wrk(posw),nhelp1,dimc,wrk(poss3),wrk(posd3),1)
      end if

    else if (symb == symc) then
      !1.c  case syma/=symb=symc

      !1.c.1---- 1st. triade   W(abc) <- +s1(_i,a) . d1(_jk,bc)
      !                                  -s1(_i,b) . d1(_jk,ac)
      !                                  +s1(_i,c) . d1(_jk,ab)

      if ((symi == syma) .and. (symjk == symbc)) then
        ! if syma=symi, then obviously symb(symc)/=symi
        !1.c.1.* find address for s1,d1
        is1 = mapih1(1,1,1)
        id1 = mapim1(symb,1,1)

        !1.c.1.* def position of s1,d1
        poss1 = mapdh1(is1,1)
        posd1 = mapdm1(id1,1)

        !1.c.1.* def additional dimensions
        nhelp1 = dimb*(dimb-1)/2

        !1.c.1.* add singly
        call t3sglh131(wrk(posw),dima,nhelp1,wrk(poss1),wrk(posd1),1)
      else if ((symb == symi) .and. (symjk == symac)) then
        !1.c.1.* find address for s3,d3
        is2 = mapih1(1,1,1)
        id2 = mapim1(syma,1,1)

        !1.c.1.* def position of s1,d3
        poss2 = mapdh1(is2,1)
        posd2 = mapdm1(id2,1)

        !1.c.1.* def additional dimensions
        nhelp1 = dimb*(dimb-1)/2

        !1.c.1.* add singly
        call t3sglh132(wrk(posw),dima,dimb,nhelp1,wrk(poss2),wrk(posd2),1)
      end if

      !1.c.2---- 2nd. triade   W(abc) <- -s1(_j,a) . d1(_ik,bc)
      !                                  +s1(_j,b) . d1(_ik,ac)
      !                                  -s1(_j,c) . d1(_ik,ab)

      if ((symj == syma) .and. (symik == symbc)) then
        ! if syma=symj, then obviously symb(symc)/=symi
        !1.c.2.* find address for s1,d1
        is1 = mapih2(1,1,1)
        id1 = mapim2(symb,1,1)

        !1.c.2.* def position of s1,d1
        poss1 = mapdh2(is1,1)
        posd1 = mapdm2(id1,1)

        !1.c.2.* def additional dimensions
        nhelp1 = dimb*(dimb-1)/2

        !1.c.2.* add singly
        call t3sglh131(wrk(posw),dima,nhelp1,wrk(poss1),wrk(posd1),-1)
      else if ((symb == symj) .and. (symik == symac)) then
        !1.c.2.* find address for s3,d3
        is2 = mapih2(1,1,1)
        id2 = mapim2(syma,1,1)

        !1.c.2.* def position of s1,d3
        poss2 = mapdh2(is2,1)
        posd2 = mapdm2(id2,1)

        !1.c.2.* def additional dimensions
        nhelp1 = dimb*(dimb-1)/2

        !1.c.2.* add singly
        call t3sglh132(wrk(posw),dima,dimb,nhelp1,wrk(poss2),wrk(posd2),-1)
      end if

      !1.c.3---- 3rd. triade   W(abc) <- +s1(_k,a) . d1(_ij,bc)
      !                                  -s1(_k,b) . d1(_ij,ac)
      !                                  +s1(_k,c) . d1(_ij,ab)

      if ((symk == syma) .and. (symij == symbc)) then
        ! if syma=symk, then obviously symb(symc)/=symk
        !1.c.3.* find address for s1,d1
        is1 = mapih3(1,1,1)
        id1 = mapim3(symb,1,1)

        !1.c.3.* def position of s1,d1
        poss1 = mapdh3(is1,1)
        posd1 = mapdm3(id1,1)

        !1.c.3.* def additional dimensions
        nhelp1 = dimb*(dimb-1)/2

        !1.c.3.* add singly
        call t3sglh131(wrk(posw),dima,nhelp1,wrk(poss1),wrk(posd1),1)
      else if ((symb == symk) .and. (symij == symac)) then
        !1.c.3.* find address for s3,d3
        is2 = mapih3(1,1,1)
        id2 = mapim3(syma,1,1)

        !1.c.3.* def position of s1,d3
        poss2 = mapdh3(is2,1)
        posd2 = mapdm3(id2,1)

        !1.c.3.* def additional dimensions
        nhelp1 = dimb*(dimb-1)/2

        !1.c.3.* add singly
        call t3sglh132(wrk(posw),dima,dimb,nhelp1,wrk(poss2),wrk(posd2),1)
      end if

    else
      !1.d case syma/=symb/=symc

      !1.d.1---- 1st. triade   W(abc) <- +s1(_i,a) . d1(_jk,bc)
      !                                  -s1(_i,b) . d1(_jk,ac)
      !                                  +s1(_i,c) . d1(_jk,ab)

      if ((symi == syma) .and. (symjk == symbc)) then
        !1.d.1.* find address for s1,d1
        is1 = mapih1(1,1,1)
        id1 = mapim1(symb,1,1)

        !1.d.1.* def position of s1,d1
        poss1 = mapdh1(is1,1)
        posd1 = mapdm1(id1,1)

        !1.d.1.* add singly
        call t3sglh141(wrk(posw),dima,dimb,dimc,wrk(poss1),wrk(posd1),1)
      else if ((symi == symb) .and. (symjk == symac)) then
        !1.d.1.* find address for s2,d2
        is2 = mapih1(1,1,1)
        id2 = mapim1(syma,1,1)

        !1.d.1.* def position of s1,d1
        poss2 = mapdh1(is2,1)
        posd2 = mapdm1(id2,1)

        !1.d.1.* add singly
        call t3sglh142(wrk(posw),dima,dimb,dimc,wrk(poss2),wrk(posd2),1)
      else if ((symi == symc) .and. (symjk == symab)) then
        !1.d.1.* find address for s1,d1
        is3 = mapih1(1,1,1)
        id3 = mapim1(syma,1,1)

        !1.d.1.* def position of s1,d1
        poss3 = mapdh1(is3,1)
        posd3 = mapdm1(id3,1)

        !1.d.1.* add singly
        call t3sglh143(wrk(posw),dima,dimb,dimc,wrk(poss3),wrk(posd3),1)
      end if

      !1.d.2---- 2nd. triade   W(abc) <- -s1(_j,a) . d1(_ik,bc)
      !                                  +s1(_j,b) . d1(_ik,ac)
      !                                  -s1(_j,c) . d1(_ik,ab)

      if ((symj == syma) .and. (symik == symbc)) then
        !1.d.2.* find address for s1,d1
        is1 = mapih2(1,1,1)
        id1 = mapim2(symb,1,1)

        !1.d.2.* def position of s1,d1
        poss1 = mapdh2(is1,1)
        posd1 = mapdm2(id1,1)

        !1.d.2.* add singly
        call t3sglh141(wrk(posw),dima,dimb,dimc,wrk(poss1),wrk(posd1),-1)
      else if ((symj == symb) .and. (symik == symac)) then
        !1.d.2.* find address for s2,d2
        is2 = mapih2(1,1,1)
        id2 = mapim2(syma,1,1)

        !1.d.2.* def position of s1,d1
        poss2 = mapdh2(is2,1)
        posd2 = mapdm2(id2,1)

        !1.d.2.* add singly
        call t3sglh142(wrk(posw),dima,dimb,dimc,wrk(poss2),wrk(posd2),-1)
      else if ((symj == symc) .and. (symik == symab)) then
        !1.d.2.* find address for s1,d1
        is3 = mapih2(1,1,1)
        id3 = mapim2(syma,1,1)

        !1.d.2.* def position of s1,d1
        poss3 = mapdh2(is3,1)
        posd3 = mapdm2(id3,1)

        !1.d.2.* add singly
        call t3sglh143(wrk(posw),dima,dimb,dimc,wrk(poss3),wrk(posd3),-1)
      end if

      !1.d.3---- 3rd. triade   W(abc) <- +s1(_k,a) . d1(_ij,bc)
      !                                  -s1(_k,b) . d1(_ij,ac)
      !                                  +s1(_k,c) . d1(_ij,ab)

      if ((symk == syma) .and. (symij == symbc)) then
        !1.d.3.* find address for s1,d1
        is1 = mapih3(1,1,1)
        id1 = mapim3(symb,1,1)

        !1.d.3.* def position of s1,d1
        poss1 = mapdh3(is1,1)
        posd1 = mapdm3(id1,1)

        !1.d.3.* add singly
        call t3sglh141(wrk(posw),dima,dimb,dimc,wrk(poss1),wrk(posd1),1)
      else if ((symk == symb) .and. (symij == symac)) then
        !1.d.3.* find address for s2,d2
        is2 = mapih3(1,1,1)
        id2 = mapim3(syma,1,1)

        !1.d.3.* def position of s1,d1
        poss2 = mapdh3(is2,1)
        posd2 = mapdm3(id2,1)

        !1.d.3.* add singly
        call t3sglh142(wrk(posw),dima,dimb,dimc,wrk(poss2),wrk(posd2),1)
      else if ((symk == symc) .and. (symij == symab)) then
        !1.d.3.* find address for s1,d1
        is3 = mapih3(1,1,1)
        id3 = mapim3(syma,1,1)

        !1.d.3.* def position of s1,d1
        poss3 = mapdh3(is3,1)
        posd3 = mapdm3(id3,1)

        !1.d.3.* add singly
        call t3sglh143(wrk(posw),dima,dimb,dimc,wrk(poss3),wrk(posd3),1)
      end if

    end if

  end do

else if (typdiv == 2) then
  !2 case W(pq,r)

  !2.* ext H1(a) <= S1(a,i) for given i
  call ext(wrk,wrksize,2,2,i,0,symi,0,0,mapds1,mapis1,1,posh10,mapdh1,mapih1,ssh1,rc1)

  !2.* ext H2(a) <= S1(a,j) for given j
  call ext(wrk,wrksize,2,2,j,0,symj,0,0,mapds1,mapis1,1,posh20,mapdh2,mapih2,ssh2,rc1)

  !2.* ext H3(a) <= S2(a,k) for given k
  call ext(wrk,wrksize,2,2,k,0,symk,0,0,mapds2,mapis2,1,posh30,mapdh3,mapih3,ssh3,rc1)

  !2.* ext M1(bc) <= D2(bc,jk) for given jk
  call ext(wrk,wrksize,4,7,j,k,symj,symk,0,mapdd2,mapid2,1,posm10,mapdm1,mapim1,ssm1,rc1)

  !2.* ext M2(bc) <= D2(bc,ik) for given ik
  call ext(wrk,wrksize,4,7,i,k,symi,symk,0,mapdd2,mapid2,1,posm20,mapdm2,mapim2,ssm2,rc1)

  !2.* ext M3(bc) <= D1(bc,ij) for given ij
  call ext(wrk,wrksize,4,7,i,j,symi,symj,0,mapdd1,mapid1,1,posm30,mapdm3,mapim3,ssm3,rc1)

  do iw=1,mapdw(0,5)

    !2.* def position of W
    posw = mapdw(iw,1)

    !2.* def symmetry status
    syma = mapdw(iw,3)
    symb = mapdw(iw,4)
    symc = mapdw(iw,5)

    !2.* def dimensions
    dima = dimm(mapdw(0,1),syma)
    dimb = dimm(mapdw(0,2),symb)
    dimc = dimm(mapdw(0,3),symc)
    !2.*
    symij = mmul(symi,symj)
    symik = mmul(symi,symk)
    symjk = mmul(symj,symk)
    symab = mmul(syma,symb)
    symac = mmul(syma,symc)
    symbc = mmul(symb,symc)

    !2.* add singles
    if (syma == symb) then
      !2.a case syma=symb,symc

      !2.a.1 1st. diade  W(ab,c)<- +s1(_i,a) . d2(_j,_k,b,c)
      !                            -s1(_i,b) . d2(_j,_k,a,c)

      if ((syma == symi) .and. (symjk == symbc)) then
        !2.a.1.* find address for s1,d1
        is1 = mapih1(1,1,1)
        id1 = mapim1(syma,1,1)

        !2.a.1.* def position of s1,d1
        poss1 = mapdh1(is1,1)
        posd1 = mapdm1(id1,1)

        !2.a.1.* def additional dimensions
        nhelp1 = dima*(dima-1)/2

        !2.a.1.* add singly
        call t3sglh211(wrk(posw),dima,nhelp1,dimc,wrk(poss1),wrk(posd1),1)
      end if

      !2.a.2 2nd. diade        W(ab,c)<- -s1(_j,a) . d2(_i,_k,b,c)
      !                                  +s1(_j,b) . d2(_i,_k,a,c)

      if ((syma == symj) .and. (symik == symbc)) then
        !2.a.2.* find address for s1,d1
        is1 = mapih2(1,1,1)
        id1 = mapim2(syma,1,1)

        !2.a.2.* def position of s1,d1
        poss1 = mapdh2(is1,1)
        posd1 = mapdm2(id1,1)

        !2.a.1.* def additional dimensions
        nhelp1 = dima*(dima-1)/2

        !2.a.2.* add singly
        call t3sglh211(wrk(posw),dima,nhelp1,dimc,wrk(poss1),wrk(posd1),-1)
      end if

      !2.a.3 3rd. part  W(ab,c)<- +s2(_k,c) . d1(_ij,ab)

      if ((symc == symk) .and. (symij == symab)) then
        !2.a.3.* find address for s1,d1
        is1 = mapih3(1,1,1)
        id1 = mapim3(syma,1,1)

        !2.a.3.* def position of s1,d1
        poss1 = mapdh3(is1,1)
        posd1 = mapdm3(id1,1)

        !2.a.3.* def additional dimensions
        nhelp1 = dima*(dima-1)/2

        !2.a.3.* add singly
        call t3sglh212(wrk(posw),nhelp1,dimc,wrk(poss1),wrk(posd1),1)
      end if

    else
      !2.b case syma>symb,symc

      !2.b.1 1st.diade  W(a>b,c)<- +s1(_i,a) . d2(_j,_k,b,c)
      !                            -s1(_i,b) . d2(_j,_k,a,c)

      if ((syma == symi) .and. (symjk == symbc)) then
        !2.b.1.* find address for s1,d1
        is1 = mapih1(1,1,1)
        id1 = mapim1(symb,1,1)

        !2.b.1.* def position of s1,d1
        poss1 = mapdh1(is1,1)
        posd1 = mapdm1(id1,1)

        !2.b.1.* add singly
        call t3sglh221(wrk(posw),dima,dimb,dimc,wrk(poss1),wrk(posd1),1)
      else if ((symb == symi) .and. (symjk == symac)) then
        !2.b.1.* find address for s2,d2
        is2 = mapih1(1,1,1)
        id2 = mapim1(syma,1,1)

        !2.b.1.* def position of s1,d1
        poss2 = mapdh1(is2,1)
        posd2 = mapdm1(id2,1)

        !2.b.1.* add singly
        call t3sglh222(wrk(posw),dima,dimb,dimc,wrk(poss2),wrk(posd2),1)
      end if

      !2.a.2 2nd. diade        W(ab,c)<- -s1(_j,a) . d2(_i,_k,b,c)
      !                                  +s1(_j,b) . d2(_i,_k,a,c)

      if ((syma == symj) .and. (symik == symbc)) then
        !2.b.2.* find address for s1,d1
        is1 = mapih2(1,1,1)
        id1 = mapim2(symb,1,1)

        !2.b.2.* def position of s1,d1
        poss1 = mapdh2(is1,1)
        posd1 = mapdm2(id1,1)

        !2.b.2.* add singly
        call t3sglh221(wrk(posw),dima,dimb,dimc,wrk(poss1),wrk(posd1),-1)
      else if ((symb == symj) .and. (symik == symac)) then
        !2.b.2.* find address for s2,d2
        is2 = mapih2(1,1,1)
        id2 = mapim2(syma,1,1)

        !2.b.2.* def position of s1,d1
        poss2 = mapdh2(is2,1)
        posd2 = mapdm2(id2,1)

        !2.b.2.* add singly
        call t3sglh222(wrk(posw),dima,dimb,dimc,wrk(poss2),wrk(posd2),-1)
      end if

      !2.a.3 3rd. part  W(ab,c)<- +s2(_k,c) . d1(_ij,ab)

      if ((symc == symk) .and. (symij == symab)) then
        !2.b.3.* find address for s1,d1
        is3 = mapih3(1,1,1)
        id3 = mapim3(syma,1,1)

        !2.b.3.* def position of s1,d1
        poss3 = mapdh3(is3,1)
        posd3 = mapdm3(id3,1)

        !2.b.3.* add singly
        call t3sglh223(wrk(posw),dima,dimb,dimc,wrk(poss3),wrk(posd3),1)
      end if

    end if

  end do

else if (typdiv == 3) then
  !3 case B(p,qr)

  !3.* ext H1(a) <= S1(a,i) for given i
  call ext(wrk,wrksize,2,2,i,0,symi,0,0,mapds1,mapis1,1,posh10,mapdh1,mapih1,ssh1,rc1)

  !3.* ext H2(a) <= S2(a,j) for given j
  call ext(wrk,wrksize,2,2,j,0,symj,0,0,mapds2,mapis2,1,posh20,mapdh2,mapih2,ssh2,rc1)

  !3.* ext H3(a) <= S2(a,k) for given k
  call ext(wrk,wrksize,2,2,k,0,symk,0,0,mapds2,mapis2,1,posh30,mapdh3,mapih3,ssh3,rc1)

  !3.* ext M1(bc) <= D2(bc,jk) for given jk
  call ext(wrk,wrksize,4,7,j,k,symj,symk,0,mapdd2,mapid2,1,posm10,mapdm1,mapim1,ssm1,rc1)

  !3.* ext M2(bc) <= D1(bc,ik) for given ik
  call ext(wrk,wrksize,4,7,i,k,symi,symk,0,mapdd1,mapid1,1,posm20,mapdm2,mapim2,ssm2,rc1)

  !3.* ext M3(bc) <= D1(bc,ij) for given ij
  call ext(wrk,wrksize,4,7,i,j,symi,symj,0,mapdd1,mapid1,1,posm30,mapdm3,mapim3,ssm3,rc1)

  do iw=1,mapdw(0,5)

    !3.* def position of W,V
    posw = mapdw(iw,1)

    !3.* def symmetry status
    syma = mapdw(iw,3)
    symb = mapdw(iw,4)
    symc = mapdw(iw,5)

    !3.* def dimensions
    dima = dimm(mapdw(0,1),syma)
    dimb = dimm(mapdw(0,2),symb)
    dimc = dimm(mapdw(0,3),symc)

    !3.* realize packing

    if (symb == symc) then
      !3.a case syma,symb=symc

      !3.a.1 1st. part  W(a,bc)<- +s1(_i,a) . d2(_jk,b,c)

      if (symi == syma) then
        !3.a.1.* find address for s1,d1
        is1 = mapih1(1,1,1)
        id1 = mapim1(symb,1,1)

        !3.a.1.* def position of s1,d1
        poss1 = mapdh1(is1,1)
        posd1 = mapdm1(id1,1)

        !3.a.1.* def additional dimensions
        nhelp1 = dimb*(dimb-1)/2

        !3.a.1.* add singly
        call t3sglh312(wrk(posw),dima,nhelp1,wrk(poss1),wrk(posd1),1)
      end if

      !3.a.2 2nd. diade        W(a,bc)<- +s1(_j,b) . d2(_i,_k,a,c)
      !                                  -s1(_j,c) . d2(_i,_k,a,b)

      if (symj == symb) then
        !3.a.2.* find address for s1,d1
        is1 = mapih2(1,1,1)
        id1 = mapim2(syma,1,1)

        !3.a.2.* def position of s1,d1
        poss1 = mapdh2(is1,1)
        posd1 = mapdm2(id1,1)

        !3.a.2.* def additional dimensions
        nhelp1 = dimb*(dimb-1)/2

        !3.a.2.* add singly
        call t3sglh311(wrk(posw),dima,dimb,nhelp1,wrk(poss1),wrk(posd1),1)
      end if

      !3.a.3 3rd. diade  W(a,bc)<- -s1(_k,b) . d2(_i,_j,a,c)
      !                            +s1(_k,c) . d2(_i,_j,a,b)

      if (symk == symb) then
        !3.a.3.* find address for s1,d1
        is1 = mapih3(1,1,1)
        id1 = mapim3(syma,1,1)

        !3.a.3.* def position of s1,d1
        poss1 = mapdh3(is1,1)
        posd1 = mapdm3(id1,1)

        !3.a.3.* def additional dimensions
        nhelp1 = dimb*(dimb-1)/2

        !3.a.3.* add singly
        call t3sglh311(wrk(posw),dima,dimb,nhelp1,wrk(poss1),wrk(posd1),-1)
      end if

    else
      !3.b case syma,symb/=symc

      !3.b.1 1st. part  W(a,bc)<- +s1(_i,a) . d2(_jk,b,c)

      if (symi == syma) then
        !3.b.1.* find address for s1,d1
        is1 = mapih1(1,1,1)
        id1 = mapim1(symb,1,1)

        !3.b.1.* def position of s1,d1
        poss1 = mapdh1(is1,1)
        posd1 = mapdm1(id1,1)

        !3.b.1.* add singly
        call t3sglh323(wrk(posw),dima,dimb,dimc,wrk(poss1),wrk(posd1),1)
      end if

      !3.b.2 2nd. diade        W(a,bc)<- +s1(_j,b) . d2(_i,_k,a,c)
      !                                  -s1(_j,c) . d2(_i,_k,a,b)

      if (symj == symb) then
        !3.b.2.* find address for s1,d1
        is1 = mapih2(1,1,1)
        id1 = mapim2(syma,1,1)

        !3.b.2.* def position of s1,d1
        poss1 = mapdh2(is1,1)
        posd1 = mapdm2(id1,1)

        !3.b.2.* add singly
        call t3sglh321(wrk(posw),dima,dimb,dimc,wrk(poss1),wrk(posd1),1)
      else if (symj == symc) then
        !3.b.2.* find address for s1,d1
        is1 = mapih2(1,1,1)
        id1 = mapim2(syma,1,1)

        !3.b.2.* def position of s1,d1
        poss1 = mapdh2(is1,1)
        posd1 = mapdm2(id1,1)

        !3.b.2.* add singly
        call t3sglh322(wrk(posw),dima,dimb,dimc,wrk(poss1),wrk(posd1),1)
      end if

      !3.b.3 3rd. diade  W(a,bc)<- -s1(_k,b) . d2(_i,_j,a,c)
      !                            +s1(_k,c) . d2(_i,_j,a,b)

      if (symk == symb) then
        !3.b.3.* find address for s1,d1
        is1 = mapih3(1,1,1)
        id1 = mapim3(syma,1,1)

        !3.b.3.* def position of s1,d1
        poss1 = mapdh3(is1,1)
        posd1 = mapdm3(id1,1)

        !3.b.3.* add singly
        call t3sglh321(wrk(posw),dima,dimb,dimc,wrk(poss1),wrk(posd1),-1)
      else if (symk == symc) then
        !3.b.3.* find address for s1,d1
        is1 = mapih3(1,1,1)
        id1 = mapim3(syma,1,1)

        !3.b.3.* def position of s1,d1
        poss1 = mapdh3(is1,1)
        posd1 = mapdm3(id1,1)

        !3.b.3.* add singly
        call t3sglh322(wrk(posw),dima,dimb,dimc,wrk(poss1),wrk(posd1),-1)
      end if

    end if

  end do

else
  ! RC=1 , typdiv is not 1,2,3 (NCI)
  return

end if

return

end subroutine t3sgl
