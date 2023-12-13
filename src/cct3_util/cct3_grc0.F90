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

subroutine cct3_grc0(nind,typ,typp,typq,typr,typs,stot,med,post)
! this routine defines med%d and med%i for given intermediate
!
! !N.B. (this routine cannot run with +OP2)

use CCT3_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nind, typ, typp, typq, typr, typs, stot
type(Map_Type), intent(inout) :: med
integer(kind=iwp), intent(out) :: post
integer(kind=iwp) :: i, nhelp1, nhelp2, nhelp3, nhelp4, nsymq, nsymr, pos, sp, spq, spqr, sq, sr, ss

! vanishing med%i files

med%i(1:nsym,1:nsym,1:nsym) = 0

select case (nind)

  case (1)

    ! matrix A(p)

    i = 1
    pos = med%pos0
    sp = mmul(stot,1)

    nhelp1 = dimm(typp,sp)

    ! def med%i
    med%i(1,1,1) = i

    ! def position
    med%d(i,1) = pos

    ! def length
    med%d(i,2) = nhelp1

    ! def sym p,q
    med%d(i,3) = sp
    med%d(i,4) = 0
    med%d(i,5) = 0
    med%d(i,6) = 0

    pos = pos+med%d(i,2)
    i = i+1

  case (2)

    ! matrix A(p,q)

    i = 1
    pos = med%pos0

    do sp=1,nsym

      sq = mmul(stot,sp)
      ! Meggie out
      if ((typ == 1) .and. (sp < sq)) cycle

      nhelp1 = dimm(typp,sp)
      nhelp2 = dimm(typq,sq)

      ! def med%i
      med%i(sp,1,1) = i

      ! def position
      med%d(i,1) = pos

      ! def length
      if ((typ == 1) .and. (sp == sq)) then
        med%d(i,2) = nhelp1*(nhelp1-1)/2
      else
        med%d(i,2) = nhelp1*nhelp2
      end if

      ! def sym p,q
      med%d(i,3) = sp
      med%d(i,4) = sq
      med%d(i,5) = 0
      med%d(i,6) = 0

      pos = pos+med%d(i,2)
      i = i+1

    end do

  case (3)

    ! matrix A(p,q,r)

    i = 1
    pos = med%pos0

    do sp=1,nsym
      if (typ == 1) then
        nsymq = sp
      else
        nsymq = nsym
      end if

      do sq=1,nsymq
        spq = mmul(sp,sq)

        sr = mmul(stot,spq)
        ! Meggie out
        if ((typ == 2) .and. (sq < sr)) cycle

        nhelp1 = dimm(typp,sp)
        nhelp2 = dimm(typq,sq)
        nhelp3 = dimm(typr,sr)

        ! def med%i
        med%i(sp,sq,1) = i

        ! def position
        med%d(i,1) = pos

        ! def length
        if ((typ == 1) .and. (sp == sq)) then
          med%d(i,2) = nhelp1*(nhelp1-1)*nhelp3/2
        else if ((typ == 2) .and. (sq == sr)) then
          med%d(i,2) = nhelp1*nhelp2*(nhelp2-1)/2
        else
          med%d(i,2) = nhelp1*nhelp2*nhelp3
        end if

        ! def sym p,q,r
        med%d(i,3) = sp
        med%d(i,4) = sq
        med%d(i,5) = sr
        med%d(i,6) = 0

        pos = pos+med%d(i,2)
        i = i+1

      end do
    end do

  case (4)

    ! matrix A(p,q,r,s)

    i = 1
    pos = med%pos0

    do sp=1,nsym
      if ((typ == 1) .or. (typ == 4)) then
        nsymq = sp
      else
        nsymq = nsym
      end if

      do sq=1,nsymq
        spq = mmul(sp,sq)
        if (typ == 2) then
          nsymr = sq
        else
          nsymr = nsym
        end if

        do sr=1,nsymr
          spqr = mmul(spq,sr)

          ss = mmul(stot,spqr)
          ! Meggie out
          if (((typ == 3) .or. (typ == 4)) .and. (sr < ss)) cycle

          nhelp1 = dimm(typp,sp)
          nhelp2 = dimm(typq,sq)
          nhelp3 = dimm(typr,sr)
          nhelp4 = dimm(typs,ss)

          ! def med%i
          med%i(sp,sq,sr) = i

          ! def position
          med%d(i,1) = pos

          ! def length
          if ((typ == 1) .and. (sp == sq)) then
            med%d(i,2) = nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
          else if ((typ == 2) .and. (sq == sr)) then
            med%d(i,2) = nhelp1*nhelp2*(nhelp3-1)*nhelp4/2
          else if ((typ == 3) .and. (sr == ss)) then
            med%d(i,2) = nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
          else if (typ == 4) then
            if ((sp == sq) .and. (sr == ss)) then
              med%d(i,2) = nhelp1*(nhelp2-1)*nhelp3*(nhelp4-1)/4
            else if (sp == sq) then
              med%d(i,2) = nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
            else if (sr == ss) then
              med%d(i,2) = nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
            else
              med%d(i,2) = nhelp1*nhelp2*nhelp3*nhelp4
            end if
          else
            med%d(i,2) = nhelp1*nhelp2*nhelp3*nhelp4
          end if

          ! def sym p,q,r,s
          med%d(i,3) = sp
          med%d(i,4) = sq
          med%d(i,5) = sr
          med%d(i,6) = ss

          pos = pos+med%d(i,2)
          i = i+1

        end do
      end do
    end do

  case default
    ! To fix some compiler warnings
    pos = 0
    i = 0
end select

post = pos

! definition of other coll

med%d(0,1) = typp
med%d(0,2) = typq
med%d(0,3) = typr
med%d(0,4) = typs
med%d(0,5) = i-1
med%d(0,6) = typ

return

end subroutine cct3_grc0
