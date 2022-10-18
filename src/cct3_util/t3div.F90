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

subroutine t3div(wrk,wrksize,mapdw,mapdv,ssw,mapdd1,mapid1,mapdd2,mapid2,typdiv,i,j,k,symi,symj,symk,ec,rc)
! mapdw  - direct map matrix of W (Input)
! mapdv  - direct map matrix of V (Input)
! ssw    - overall symmetry state of matrix W (V) (Input)
! mapdd1 - direct map matrix of 1st diagonal el. (Input)
! mapid1 - inverse map matrix of 1st diagonal el. (Input)
! mapdd2 - direct map matrix of 2nd diagonal el. (Input)
! mapid2 - inverse map matrix of 2nd diagonal el. (Input)
!          (if there is only one spin, use any map's for 2)
! typdiv - typ of operation (see Table) (Input)
! i      - value of occupied index i (Inlut)
! j      - value of occupied index j (Inlut)
! k      - value of occupied index k (Inlut)
! symi   - symmetry of index i (Input)
! symj   - symmetry of index j (Input)
! symk   - symmetry of index k (Input)
! ec     - energy contribution from this part (Output)
! rc     - return (error) code (Output)
!
! this routine realizes division by denominators and
! calculate corresponding energy contribution
!
! ec = sum (abc) [ W(abc).V(abc) / Dijkabc ]
!
! for following types of W,V
!
! typdiv         Operation                Implemented
! 1     W(abc)  . V(abc) /Dijkabc           Yes
! 2     W(ab,c) . W(ab,c)/Dijkabc           Yes
! 3     W(a,bc) . V(a,bc)/Dijkabc           Yes
!
! N.B. spin combinations aaa,bbb for 1; aab for 2; and abb for 3
! are automatically assumed

#include "t31.fh"
#include "wrk.fh"
integer ssw, typdiv, i, j, k, symi, symj, symk, rc
integer mapdw(0:512,1:6)
integer mapdv(0:512,1:6)
integer mapdd1(0:512,1:6)
integer mapdd2(0:512,1:6)
!integer mapiw(1:8,1:8,1:8)
!integer mapiv(1:8,1:8,1:8)
integer mapid1(1:8,1:8,1:8)
integer mapid2(1:8,1:8,1:8)
real*8 ec
! help variables
integer iw, possw, possv
integer id1, id2, id3, possd1, possd2, possd3
integer syma, symb, symc, dima, dimb, dimc
integer nhelp1, nhelp2, nhelp3, nhelp4, nhelp5, nhelp6
integer index
real*8 eco, denijk

!0.*  some tests

if (mapdw(0,6) /= mapdv(0,6)) then
  ! RC=1 : typW is not equal to typV (Stup)
  rc = 1
  return
end if

if ((typdiv == 1) .and. (mapdw(0,6) /= 5)) then
  ! RC=2 : typdiv=1, typW is not 5 (Stup)
  rc = 2
  return
end if

if ((typdiv == 2) .and. (mapdw(0,6) /= 1)) then
  ! RC=3 : typdiv=2, typW is not 1 (Stup)
  rc = 3
  return
end if

if ((typdiv == 3) .and. (mapdw(0,6) /= 2)) then
  ! RC=4 : typdiv=3, typW is not 2 (Stup)
  rc = 4
  return
end if

!0.* vanish ec

ec = 0.0d0

!0.* def denijk

if (typdiv == 1) then
  ! cases aaa,bbb

  ! diagonal part i
  id1 = mapid1(symi,1,1)
  possd1 = mapdd1(id1,1)
  index = possd1+i-1
  denijk = wrk(index)

  ! diagonal part j
  id1 = mapid1(symj,1,1)
  possd1 = mapdd1(id1,1)
  index = possd1+j-1
  denijk = denijk+wrk(index)

  ! diagonal part k
  id1 = mapid1(symk,1,1)
  possd1 = mapdd1(id1,1)
  index = possd1+k-1
  denijk = denijk+wrk(index)

else if (typdiv == 2) then
  ! case aab

  ! diagonal part i
  id1 = mapid1(symi,1,1)
  possd1 = mapdd1(id1,1)
  index = possd1+i-1
  denijk = wrk(index)

  ! diagonal part j
  id1 = mapid1(symj,1,1)
  possd1 = mapdd1(id1,1)
  index = possd1+j-1
  denijk = denijk+wrk(index)

  ! diagonal part k
  id2 = mapid2(symk,1,1)
  possd2 = mapdd2(id2,1)
  index = possd2+k-1
  denijk = denijk+wrk(index)

else if (typdiv == 3) then
  ! case abb

  ! diagonal part i
  id1 = mapid1(symi,1,1)
  possd1 = mapdd1(id1,1)
  index = possd1+i-1
  denijk = wrk(index)

  ! diagonal part j
  id2 = mapid2(symj,1,1)
  possd2 = mapdd2(id2,1)
  index = possd2+j-1
  denijk = denijk+wrk(index)

  ! diagonal part k
  id2 = mapid2(symk,1,1)
  possd2 = mapdd2(id2,1)
  index = possd2+k-1
  denijk = denijk+wrk(index)

end if

if (typdiv == 1) then
  !1 case W(pqr). V(pqr)

  do iw=1,mapdw(0,5)

    !1.* def position of W,V
    possw = mapdw(iw,1)
    possv = mapdv(iw,1)

    !1.* def symmetry status
    syma = mapdw(iw,3)
    symb = mapdw(iw,4)
    symc = mapdw(iw,5)

    !1.* def dimensions
    dima = dimm(mapdw(0,1),syma)
    dimb = dimm(mapdw(0,2),symb)
    dimc = dimm(mapdw(0,3),symc)

    !1.* realize packing

    if (syma == symc) then
      !1.a case syma=symb=symc

      !1.a.* find address for d1,2,3 (the same one)
      id1 = mapid1(syma,1,1)

      !1.a.* def position of d1,2,3 (the same one)
      possd1 = mapdd1(id1,1)

      !1.a.* def additional dimensions
      nhelp1 = dima*(dima-1)*(dima-2)/6
      nhelp2 = dimm(5,syma)
      if (mapdw(0,1) == 3) then
        ! alpha cese
        nhelp3 = noa(syma)
      else
        ! beta case
        nhelp3 = nob(syma)
      end if

      !1.a.* do packing
      call t3dhlp4(wrk(possw),wrk(possv),dima,nhelp1,denijk,eco,wrk(possd1),nhelp2,nhelp3)
      ec = ec+eco

    else if (syma == symb) then
      !1.b case syma=symb/=symc

      !1.b.* find address for d1,3
      id1 = mapid1(syma,1,1)
      id3 = mapid1(symc,1,1)

      !1.b.* def position of d1,3
      possd1 = mapdd1(id1,1)
      possd3 = mapdd1(id3,1)

      !1.b.* def additional dimensions
      nhelp1 = dima*(dima-1)/2
      nhelp2 = dimm(5,syma)
      nhelp3 = dimm(5,symc)
      if (mapdw(0,1) == 3) then
        ! alpha cese
        nhelp4 = noa(syma)
        nhelp5 = noa(symc)
      else
        ! beta case
        nhelp4 = nob(syma)
        nhelp5 = nob(symc)
      end if

      !1.b.* do packing
      call t3dhlp2(wrk(possw),wrk(possv),dima,nhelp1,dimc,denijk,eco,wrk(possd1),wrk(possd3),nhelp2,nhelp3,nhelp4,nhelp5)
      ec = ec+eco

    else if (symb == symc) then
      !1.c case syma/=symb=symc

      !1.c.* find address for d1,2
      id1 = mapid1(syma,1,1)
      id2 = mapid1(symb,1,1)

      !1.c.* def position of d1,2
      possd1 = mapdd1(id1,1)
      possd2 = mapdd1(id2,1)

      !1.c.* def additional dimensions
      nhelp1 = dimb*(dimb-1)/2
      nhelp2 = dimm(5,syma)
      nhelp3 = dimm(5,symb)
      if (mapdw(0,1) == 3) then
        ! alpha cese
        nhelp4 = noa(syma)
        nhelp5 = noa(symb)
      else
        ! beta case
        nhelp4 = nob(syma)
        nhelp5 = nob(symb)
      end if

      !1.c.* do packing
      call t3dhlp3(wrk(possw),wrk(possv),dima,dimb,nhelp1,denijk,eco,wrk(possd1),wrk(possd2),nhelp2,nhelp3,nhelp4,nhelp5)
      ec = ec+eco

    else
      !1.d case syma/=symb/=symc

      !1.d.* find address for d1,2,3
      id1 = mapid1(syma,1,1)
      id2 = mapid1(symb,1,1)
      id3 = mapid1(symc,1,1)

      !1.d.* def position of d1,2,3
      possd1 = mapdd1(id1,1)
      possd2 = mapdd1(id2,1)
      possd3 = mapdd1(id3,1)

      !1.d.* def additional dimensions
      nhelp1 = dimm(5,syma)
      nhelp2 = dimm(5,symb)
      nhelp3 = dimm(5,symc)
      if (mapdw(0,1) == 3) then
        ! alpha cese
        nhelp4 = noa(syma)
        nhelp5 = noa(symb)
        nhelp6 = noa(symc)
      else
        ! beta case
        nhelp4 = nob(syma)
        nhelp5 = nob(symb)
        nhelp6 = nob(symc)
      end if

      !1.d.* do packing
      call t3dhlp1(wrk(possw),wrk(possv),dima,dimb,dimc,denijk,eco,wrk(possd1),wrk(possd2),wrk(possd3),nhelp1,nhelp2,nhelp3, &
                   nhelp4,nhelp5,nhelp6)
      ec = ec+eco

    end if

  end do

else if (typdiv == 2) then
  !2 case W(pq,r) . V(pq,r)

  do iw=1,mapdw(0,5)

    !2.* def position of W,W
    possw = mapdw(iw,1)
    possv = mapdv(iw,1)

    !2.* def symmetry status
    syma = mapdw(iw,3)
    symb = mapdw(iw,4)
    symc = mapdw(iw,5)

    !2.* def dimensions
    dima = dimm(mapdw(0,1),syma)
    dimb = dimm(mapdw(0,2),symb)
    dimc = dimm(mapdw(0,3),symc)

    !2.* realize packing

    if (syma == symb) then
      !2.a case syma=symb,symc

      !2.a.* find address for d1,3
      id1 = mapid1(syma,1,1)
      id3 = mapid2(symc,1,1)

      !2.a.* def position of d1,3
      possd1 = mapdd1(id1,1)
      possd3 = mapdd2(id3,1)

      !2.a.* def additional dimensions
      nhelp1 = dima*(dima-1)/2
      nhelp2 = dimm(5,syma)
      nhelp3 = dimm(5,symc)
      nhelp4 = noa(syma)
      nhelp5 = nob(symc)

      !2.a.* do packing
      call t3dhlp2(wrk(possw),wrk(possv),dima,nhelp1,dimc,denijk,eco,wrk(possd1),wrk(possd3),nhelp2,nhelp3,nhelp4,nhelp5)
      ec = ec+eco

    else
      !2.b case syma/=symb,symc

      !2.b.* find address for d1,2,3
      id1 = mapid1(syma,1,1)
      id2 = mapid1(symb,1,1)
      id3 = mapid2(symc,1,1)

      !2.b.* def position of d1,2,3
      possd1 = mapdd1(id1,1)
      possd2 = mapdd1(id2,1)
      possd3 = mapdd2(id3,1)

      !2.b.* def additional dimensions
      nhelp1 = dimm(5,syma)
      nhelp2 = dimm(5,symb)
      nhelp3 = dimm(5,symb)
      nhelp4 = noa(syma)
      nhelp5 = noa(symb)
      nhelp6 = nob(symc)

      !2.b.* do packing
      call t3dhlp1(wrk(possw),wrk(possv),dima,dimb,dimc,denijk,eco,wrk(possd1),wrk(possd2),wrk(possd3),nhelp1,nhelp2,nhelp3, &
                   nhelp4,nhelp5,nhelp6)
      ec = ec+eco

    end if

  end do

else if (typdiv == 3) then
  !3 case B(p,qr) = B(p,qr) + ns*(A(p,r,q)-A(p,q,r))

  do iw=1,mapdw(0,5)

    !3.* def position of W,V
    possw = mapdw(iw,1)
    possv = mapdv(iw,1)

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

      !3.a.* find address for d1,2
      id1 = mapid1(syma,1,1)
      id2 = mapid2(symb,1,1)

      !3.a.* def position of d1,2
      possd1 = mapdd1(id1,1)
      possd2 = mapdd2(id2,1)

      !3.a.* def additional dimensions
      nhelp1 = dimb*(dimb-1)/2
      nhelp2 = dimm(5,syma)
      nhelp3 = dimm(5,symb)
      nhelp4 = noa(syma)
      nhelp5 = nob(symb)

      !3.a.* do packing
      call t3dhlp3(wrk(possw),wrk(possv),dima,dimb,nhelp1,denijk,eco,wrk(possd1),wrk(possd2),nhelp2,nhelp3,nhelp4,nhelp5)
      ec = ec+eco

    else
      !3.b case syma,symb/=symc

      !3.b.* find address for d1,2,3
      id1 = mapid1(syma,1,1)
      id2 = mapid2(symb,1,1)
      id3 = mapid2(symc,1,1)

      !3.b.* def position of d1,2,3
      possd1 = mapdd1(id1,1)
      possd2 = mapdd2(id2,1)
      possd3 = mapdd2(id3,1)

      !3.b.* def additional dimensions
      nhelp1 = dimm(5,syma)
      nhelp2 = dimm(5,symb)
      nhelp3 = dimm(5,symc)
      nhelp4 = noa(syma)
      nhelp5 = nob(symb)
      nhelp6 = nob(symc)

      !3.b.* do packing
      call t3dhlp1(wrk(possw),wrk(possv),dima,dimb,dimc,denijk,eco,wrk(possd1),wrk(possd2),wrk(possd3),nhelp1,nhelp2,nhelp3, &
                   nhelp4,nhelp5,nhelp6)
      ec = ec+eco

    end if

  end do

else
  ! RC=5 , typdiv is not 1,2,3 (NCI)
  rc = 5
  return

end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(ssw)

end subroutine t3div
