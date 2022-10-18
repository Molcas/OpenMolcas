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

subroutine t3addpck(wrk,wrksize,nind,typap,mapda,mapia,mapdb,ns,szkey,rc)
! nind  - # of indices in matrices A,B (Input)
! typap - typ of operation (see Table) (Input)
! mapda - direct map matrix corresponding to A (Input)
! mapia - inverse map matrix corresponding to A (Input)
! mapdb - direct map matrix corresponding to B (Input)
! ns    - signum of the operation (+-1) (Input)
! szkey - zet zero key (I)
!         = 0 no vanishing
!         = 1 set B=0 at the beginning
! rc    - return (error) code (Output)
!
! this routine is an addition to CCSD add/pack routine. It does not
! realize operations which are done there (add/pack), but only
! additional ones, required in T3 code, namely:
!
! Operation                                   Nind TypA TypB typap
! B(pqr)  = B + ns*(A(qr,p)-A(pr,q)+A(pq,r))    3    1    5    1
! B(pq,r) = B + ns*(A(q,r,p)-A(p,r,q))          3    1    1    2
! B(p,qr) = B + ns*(     -A(p,r,q)+A(p,q,r))    3    1    2    3
!
! N.B. typab is redundant, it can be determined from mapd's

#include "t31.fh"
#include "wrk.fh"
integer nind, typap, ns, szkey, rc
integer mapda(0:512,1:6)
integer mapdb(0:512,1:6)
integer mapia(1:8,1:8,1:8)
! help variables
integer ib, possb
integer ia1, ia2, ia3, possa1, possa2, possa3
integer syma, symb, symc, dima, dimb, dimc
integer nhelp1, nhelp2

!0 some tests

if (nind /= 3) then
  ! RC=1 : nind is not 3 (NCI)
  rc = 1
  return
end if

if ((typap == 1) .and. ((mapda(0,6) /= 1) .and. (mapdb(0,6) /= 5))) then
  ! RC=2 : typap=1, typA=1, and typB is not 5 (Stup)
  rc = 2
  return
end if

if ((typap == 2) .and. ((mapda(0,6) /= 0) .and. (mapdb(0,6) /= 1))) then
  ! RC=3 : typap=2, typA=0, and typB is not 1 (Stup)
  rc = 3
  return
end if

if ((typap == 3) .and. ((mapda(0,6) /= 0) .and. (mapdb(0,6) /= 2))) then
  ! RC=4 : typap=3, typA=0, and typB is not 2 (Stup)
  rc = 4
  return
end if

if (typap == 1) then
  !1 case B(pqr) = B(pqr) + ns*(A(qr,p)-A(pr,q)+A(pq,r))

  do ib=1,mapdb(0,5)

    !1.* def position of b
    possb = mapdb(ib,1)

    !1.* def symmetry status
    syma = mapdb(ib,3)
    symb = mapdb(ib,4)
    symc = mapdb(ib,5)

    !1.* def dimensions
    dima = dimm(mapdb(0,1),syma)
    dimb = dimm(mapdb(0,2),symb)
    dimc = dimm(mapdb(0,3),symc)

    !1.* realize packing

    if (syma == symc) then
      !1.a case syma=symb=symc

      !1.a.* find address for a1,2,3 (the same one)
      ia1 = mapia(symb,symc,1)

      !1.a.* def position of a1,2,3 (the same one)
      possa1 = mapda(ia1,1)

      !1.a.* def additional dimensions ab,abc
      nhelp1 = dima*(dima-1)/2
      nhelp2 = dima*(dima-1)*(dima-2)/6

      !1.a.* do packing
      call t3aphlp4(wrk(possa1),wrk(possb),dima,nhelp1,nhelp2,ns,szkey)

    else if (syma == symb) then
      !1.b case syma=symb/=symc

      !1.b.* find address for a1,2,3
      ia1 = mapia(symb,symc,1)
      ia2 = mapia(syma,symc,1)
      ia3 = mapia(syma,symb,1)

      !1.b.* def position of a1,2,3
      possa1 = mapda(ia1,1)
      possa2 = mapda(ia2,1)
      possa3 = mapda(ia3,1)

      !1.b.* def additional dimensions ab
      nhelp1 = dima*(dima-1)/2

      !1.b.* do packing
      call t3aphlp2(wrk(possa1),wrk(possa2),wrk(possa3),wrk(possb),dima,dimb,dimc,nhelp1,ns,szkey)

    else if (symb == symc) then
      !1.c case syma/=symb=symc

      !1.c.* find address for a1,2,3
      ia1 = mapia(symb,symc,1)
      ia2 = mapia(syma,symc,1)
      ia3 = mapia(syma,symb,1)

      !1.c.* def position of a1,2,3
      possa1 = mapda(ia1,1)
      possa2 = mapda(ia2,1)
      possa3 = mapda(ia3,1)

      !1.c.* def additional dimensions bc
      nhelp1 = dimb*(dimb-1)/2

      !1.c.* do packing
      call t3aphlp3(wrk(possa1),wrk(possa2),wrk(possa3),wrk(possb),dima,dimb,dimc,nhelp1,ns,szkey)

    else
      !1.d case syma/=symb/=symc

      !1.d.* find address for a1,2,3
      ia1 = mapia(symb,symc,1)
      ia2 = mapia(syma,symc,1)
      ia3 = mapia(syma,symb,1)

      !1.d.* def position of a1,2,3
      possa1 = mapda(ia1,1)
      possa2 = mapda(ia2,1)
      possa3 = mapda(ia3,1)

      !1.d.* do packing
      call t3aphlp1(wrk(possa1),wrk(possa2),wrk(possa3),wrk(possb),dima,dimb,dimc,ns,szkey)

    end if

  end do

else if (typap == 2) then
  !2 case B(pq,r) = B(pq,r) + ns*(A(q,r,p)-A(p,r,q))

  do ib=1,mapdb(0,5)

    !2.* def position of b
    possb = mapdb(ib,1)

    !2.* def symmetry status
    syma = mapdb(ib,3)
    symb = mapdb(ib,4)
    symc = mapdb(ib,5)

    !2.* def dimensions
    dima = dimm(mapdb(0,1),syma)
    dimb = dimm(mapdb(0,2),symb)
    dimc = dimm(mapdb(0,3),symc)

    !2.* realize packing

    if (syma == symb) then
      !2.a case syma=symb,symc

      !2.a.* find address for a1,2
      ia1 = mapia(symb,symc,1)
      ia2 = mapia(syma,symc,1)

      !2.a.* def position of a1,2,3
      possa1 = mapda(ia1,1)
      possa2 = mapda(ia2,1)

      !2.a.* def additional dimensions ab
      nhelp1 = dima*(dima-1)/2

      !2.a.* do packing
      call t3aphlp6(wrk(possa1),wrk(possa2),wrk(possb),dima,dimb,dimc,nhelp1,ns,szkey)

    else
      !2.b case syma/=symb,symc

      !2.b.* find address for a1,2
      ia1 = mapia(symb,symc,1)
      ia2 = mapia(syma,symc,1)

      !2.b.* def position of a1,2,3
      possa1 = mapda(ia1,1)
      possa2 = mapda(ia2,1)

      !2.b.* do packing
      call t3aphlp5(wrk(possa1),wrk(possa2),wrk(possb),dima,dimb,dimc,ns,szkey)

    end if

  end do

else if (typap == 3) then
  !3 case B(p,qr) = B(p,qr) + ns*(-A(p,r,q)+A(p,q,r))

  do ib=1,mapdb(0,5)

    !3.* def position of b
    possb = mapdb(ib,1)

    !3.* def symmetry status
    syma = mapdb(ib,3)
    symb = mapdb(ib,4)
    symc = mapdb(ib,5)

    !3.* def dimensions
    dima = dimm(mapdb(0,1),syma)
    dimb = dimm(mapdb(0,2),symb)
    dimc = dimm(mapdb(0,3),symc)

    !3.* realize packing

    if (symb == symc) then
      !3.a case syma,symb=symc

      !3.a.* find address for a2,3
      ia2 = mapia(syma,symc,1)
      ia3 = mapia(syma,symb,1)

      !3.a.* def position of a2,3
      possa2 = mapda(ia2,1)
      possa3 = mapda(ia3,1)

      !3.a.* def additional dimensions bc
      nhelp1 = dimb*(dimb-1)/2

      !3.a.* do packing
      call t3aphlp8(wrk(possa2),wrk(possb),dima,dimb,nhelp1,ns,szkey)

    else
      !3.b case syma,symb/=symc

      !3.b.* find address for a2,3
      ia2 = mapia(syma,symc,1)
      ia3 = mapia(syma,symb,1)

      !3.b.* def position of a2,3
      possa2 = mapda(ia2,1)
      possa3 = mapda(ia3,1)

      !3.b.* do packing
      call t3aphlp7(wrk(possa2),wrk(possa3),wrk(possb),dima,dimb,dimc,ns,szkey)

    end if

  end do

else
  ! RC=5 , typap is not 1,2,3 (NCI)
  rc = 5
  return

end if

return

end subroutine t3addpck
