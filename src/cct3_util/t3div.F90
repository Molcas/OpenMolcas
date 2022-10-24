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

subroutine t3div(wrk,wrksize,w,v,d1,d2,typdiv,i,j,k,symi,symj,symk,ec,rc)
! w      - W (Input)
! v      - V (Input)
! d1     - 1st diagonal el. (Input)
! d2     - 2nd diagonal el. (Input)
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

use CCT3_global, only: dimm, Map_Type, noa, nob
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, typdiv, i, j, k, symi, symj, symk
real(kind=wp), intent(in) :: wrk(wrksize)
type(Map_Type), intent(in) :: w, v, d1, d2
real(kind=wp), intent(out) :: ec
integer(kind=iwp), intent(inout) :: rc
integer(kind=iwp) :: dima, dimb, dimc, id1, id2, id3, indx, iw, nhelp1, nhelp2, nhelp3, nhelp4, nhelp5, nhelp6, posd1, posd2, &
                     posd3, posv, posw, syma, symb, symc
real(kind=wp) :: denijk, eco

!0.*  some tests

if (w%d(0,6) /= v%d(0,6)) then
  ! RC=1 : typW is not equal to typV (Stup)
  rc = 1
  return
end if

if ((typdiv == 1) .and. (w%d(0,6) /= 5)) then
  ! RC=2 : typdiv=1, typW is not 5 (Stup)
  rc = 2
  return
end if

if ((typdiv == 2) .and. (w%d(0,6) /= 1)) then
  ! RC=3 : typdiv=2, typW is not 1 (Stup)
  rc = 3
  return
end if

if ((typdiv == 3) .and. (w%d(0,6) /= 2)) then
  ! RC=4 : typdiv=3, typW is not 2 (Stup)
  rc = 4
  return
end if

!0.* vanish ec

ec = Zero

!0.* def denijk

if (typdiv == 1) then
  ! cases aaa,bbb

  ! diagonal part i
  id1 = d1%i(symi,1,1)
  posd1 = d1%d(id1,1)
  indx = posd1+i-1
  denijk = wrk(indx)

  ! diagonal part j
  id1 = d1%i(symj,1,1)
  posd1 = d1%d(id1,1)
  indx = posd1+j-1
  denijk = denijk+wrk(indx)

  ! diagonal part k
  id1 = d1%i(symk,1,1)
  posd1 = d1%d(id1,1)
  indx = posd1+k-1
  denijk = denijk+wrk(indx)

else if (typdiv == 2) then
  ! case aab

  ! diagonal part i
  id1 = d1%i(symi,1,1)
  posd1 = d1%d(id1,1)
  indx = posd1+i-1
  denijk = wrk(indx)

  ! diagonal part j
  id1 = d1%i(symj,1,1)
  posd1 = d1%d(id1,1)
  indx = posd1+j-1
  denijk = denijk+wrk(indx)

  ! diagonal part k
  id2 = d2%i(symk,1,1)
  posd2 = d2%d(id2,1)
  indx = posd2+k-1
  denijk = denijk+wrk(indx)

else if (typdiv == 3) then
  ! case abb

  ! diagonal part i
  id1 = d1%i(symi,1,1)
  posd1 = d1%d(id1,1)
  indx = posd1+i-1
  denijk = wrk(indx)

  ! diagonal part j
  id2 = d2%i(symj,1,1)
  posd2 = d2%d(id2,1)
  indx = posd2+j-1
  denijk = denijk+wrk(indx)

  ! diagonal part k
  id2 = d2%i(symk,1,1)
  posd2 = d2%d(id2,1)
  indx = posd2+k-1
  denijk = denijk+wrk(indx)

end if

if (typdiv == 1) then
  !1 case W(pqr). V(pqr)

  do iw=1,w%d(0,5)

    !1.* def position of W,V
    posw = w%d(iw,1)
    posv = v%d(iw,1)

    !1.* def symmetry status
    syma = w%d(iw,3)
    symb = w%d(iw,4)
    symc = w%d(iw,5)

    !1.* def dimensions
    dima = dimm(w%d(0,1),syma)
    dimb = dimm(w%d(0,2),symb)
    dimc = dimm(w%d(0,3),symc)

    !1.* realize packing

    if (syma == symc) then
      !1.a case syma=symb=symc

      !1.a.* find address for d1,2,3 (the same one)
      id1 = d1%i(syma,1,1)

      !1.a.* def position of d1,2,3 (the same one)
      posd1 = d1%d(id1,1)

      !1.a.* def additional dimensions
      nhelp1 = dima*(dima-1)*(dima-2)/6
      nhelp2 = dimm(5,syma)
      if (w%d(0,1) == 3) then
        ! alpha cese
        nhelp3 = noa(syma)
      else
        ! beta case
        nhelp3 = nob(syma)
      end if

      !1.a.* do packing
      call t3dhlp4(wrk(posw),wrk(posv),dima,nhelp1,denijk,eco,wrk(posd1),nhelp2,nhelp3)
      ec = ec+eco

    else if (syma == symb) then
      !1.b case syma=symb/=symc

      !1.b.* find address for d1,3
      id1 = d1%i(syma,1,1)
      id3 = d1%i(symc,1,1)

      !1.b.* def position of d1,3
      posd1 = d1%d(id1,1)
      posd3 = d1%d(id3,1)

      !1.b.* def additional dimensions
      nhelp1 = dima*(dima-1)/2
      nhelp2 = dimm(5,syma)
      nhelp3 = dimm(5,symc)
      if (w%d(0,1) == 3) then
        ! alpha cese
        nhelp4 = noa(syma)
        nhelp5 = noa(symc)
      else
        ! beta case
        nhelp4 = nob(syma)
        nhelp5 = nob(symc)
      end if

      !1.b.* do packing
      call t3dhlp2(wrk(posw),wrk(posv),dima,nhelp1,dimc,denijk,eco,wrk(posd1),wrk(posd3),nhelp2,nhelp3,nhelp4,nhelp5)
      ec = ec+eco

    else if (symb == symc) then
      !1.c case syma/=symb=symc

      !1.c.* find address for d1,2
      id1 = d1%i(syma,1,1)
      id2 = d1%i(symb,1,1)

      !1.c.* def position of d1,2
      posd1 = d1%d(id1,1)
      posd2 = d1%d(id2,1)

      !1.c.* def additional dimensions
      nhelp1 = dimb*(dimb-1)/2
      nhelp2 = dimm(5,syma)
      nhelp3 = dimm(5,symb)
      if (w%d(0,1) == 3) then
        ! alpha cese
        nhelp4 = noa(syma)
        nhelp5 = noa(symb)
      else
        ! beta case
        nhelp4 = nob(syma)
        nhelp5 = nob(symb)
      end if

      !1.c.* do packing
      call t3dhlp3(wrk(posw),wrk(posv),dima,dimb,nhelp1,denijk,eco,wrk(posd1),wrk(posd2),nhelp2,nhelp3,nhelp4,nhelp5)
      ec = ec+eco

    else
      !1.d case syma/=symb/=symc

      !1.d.* find address for d1,2,3
      id1 = d1%i(syma,1,1)
      id2 = d1%i(symb,1,1)
      id3 = d1%i(symc,1,1)

      !1.d.* def position of d1,2,3
      posd1 = d1%d(id1,1)
      posd2 = d1%d(id2,1)
      posd3 = d1%d(id3,1)

      !1.d.* def additional dimensions
      nhelp1 = dimm(5,syma)
      nhelp2 = dimm(5,symb)
      nhelp3 = dimm(5,symc)
      if (w%d(0,1) == 3) then
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
      call t3dhlp1(wrk(posw),wrk(posv),dima,dimb,dimc,denijk,eco,wrk(posd1),wrk(posd2),wrk(posd3),nhelp1,nhelp2,nhelp3,nhelp4, &
                   nhelp5,nhelp6)
      ec = ec+eco

    end if

  end do

else if (typdiv == 2) then
  !2 case W(pq,r) . V(pq,r)

  do iw=1,w%d(0,5)

    !2.* def position of W,W
    posw = w%d(iw,1)
    posv = v%d(iw,1)

    !2.* def symmetry status
    syma = w%d(iw,3)
    symb = w%d(iw,4)
    symc = w%d(iw,5)

    !2.* def dimensions
    dima = dimm(w%d(0,1),syma)
    dimb = dimm(w%d(0,2),symb)
    dimc = dimm(w%d(0,3),symc)

    !2.* realize packing

    if (syma == symb) then
      !2.a case syma=symb,symc

      !2.a.* find address for d1,3
      id1 = d1%i(syma,1,1)
      id3 = d2%i(symc,1,1)

      !2.a.* def position of d1,3
      posd1 = d1%d(id1,1)
      posd3 = d2%d(id3,1)

      !2.a.* def additional dimensions
      nhelp1 = dima*(dima-1)/2
      nhelp2 = dimm(5,syma)
      nhelp3 = dimm(5,symc)
      nhelp4 = noa(syma)
      nhelp5 = nob(symc)

      !2.a.* do packing
      call t3dhlp2(wrk(posw),wrk(posv),dima,nhelp1,dimc,denijk,eco,wrk(posd1),wrk(posd3),nhelp2,nhelp3,nhelp4,nhelp5)
      ec = ec+eco

    else
      !2.b case syma/=symb,symc

      !2.b.* find address for d1,2,3
      id1 = d1%i(syma,1,1)
      id2 = d1%i(symb,1,1)
      id3 = d2%i(symc,1,1)

      !2.b.* def position of d1,2,3
      posd1 = d1%d(id1,1)
      posd2 = d1%d(id2,1)
      posd3 = d2%d(id3,1)

      !2.b.* def additional dimensions
      nhelp1 = dimm(5,syma)
      nhelp2 = dimm(5,symb)
      nhelp3 = dimm(5,symb)
      nhelp4 = noa(syma)
      nhelp5 = noa(symb)
      nhelp6 = nob(symc)

      !2.b.* do packing
      call t3dhlp1(wrk(posw),wrk(posv),dima,dimb,dimc,denijk,eco,wrk(posd1),wrk(posd2),wrk(posd3),nhelp1,nhelp2,nhelp3,nhelp4, &
                   nhelp5,nhelp6)
      ec = ec+eco

    end if

  end do

else if (typdiv == 3) then
  !3 case B(p,qr) = B(p,qr) + ns*(A(p,r,q)-A(p,q,r))

  do iw=1,w%d(0,5)

    !3.* def position of W,V
    posw = w%d(iw,1)
    posv = v%d(iw,1)

    !3.* def symmetry status
    syma = w%d(iw,3)
    symb = w%d(iw,4)
    symc = w%d(iw,5)

    !3.* def dimensions
    dima = dimm(w%d(0,1),syma)
    dimb = dimm(w%d(0,2),symb)
    dimc = dimm(w%d(0,3),symc)

    !3.* realize packing

    if (symb == symc) then
      !3.a case syma,symb=symc

      !3.a.* find address for d1,2
      id1 = d1%i(syma,1,1)
      id2 = d2%i(symb,1,1)

      !3.a.* def position of d1,2
      posd1 = d1%d(id1,1)
      posd2 = d2%d(id2,1)

      !3.a.* def additional dimensions
      nhelp1 = dimb*(dimb-1)/2
      nhelp2 = dimm(5,syma)
      nhelp3 = dimm(5,symb)
      nhelp4 = noa(syma)
      nhelp5 = nob(symb)

      !3.a.* do packing
      call t3dhlp3(wrk(posw),wrk(posv),dima,dimb,nhelp1,denijk,eco,wrk(posd1),wrk(posd2),nhelp2,nhelp3,nhelp4,nhelp5)
      ec = ec+eco

    else
      !3.b case syma,symb/=symc

      !3.b.* find address for d1,2,3
      id1 = d1%i(syma,1,1)
      id2 = d2%i(symb,1,1)
      id3 = d2%i(symc,1,1)

      !3.b.* def position of d1,2,3
      posd1 = d1%d(id1,1)
      posd2 = d2%d(id2,1)
      posd3 = d2%d(id3,1)

      !3.b.* def additional dimensions
      nhelp1 = dimm(5,syma)
      nhelp2 = dimm(5,symb)
      nhelp3 = dimm(5,symc)
      nhelp4 = noa(syma)
      nhelp5 = nob(symb)
      nhelp6 = nob(symc)

      !3.b.* do packing
      call t3dhlp1(wrk(posw),wrk(posv),dima,dimb,dimc,denijk,eco,wrk(posd1),wrk(posd2),wrk(posd3),nhelp1,nhelp2,nhelp3,nhelp4, &
                   nhelp5,nhelp6)
      ec = ec+eco

    end if

  end do

else
  ! RC=5 , typdiv is not 1,2,3 (NCI)
  rc = 5
  return

end if

return

end subroutine t3div
