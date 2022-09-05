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

subroutine ass1a(D01,D1,PAO,tmp1_,nt,nrys)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nt, nrys
real(kind=wp), intent(in) :: D01(nrys,nt), D1(nrys,nt), PAO(nt)
real(kind=wp), intent(inout) :: tmp1_
integer(kind=iwp) :: iRys, iT
real(kind=wp) :: tmp1

tmp1 = Zero

select case (nRys)

  case (1)
    do iT=1,nt
      tmp1 = tmp1+(D01(1,iT)*PAO(iT))*D1(1,iT)
    end do

  case (2)
    do iT=1,nt
      tmp1 = tmp1+(D01(1,iT)*D1(1,iT)+D01(2,iT)*D1(2,iT))*PAO(iT)
    end do

  case (3)
    do iT=1,nt
      tmp1 = tmp1+(D01(1,iT)*D1(1,iT)+D01(2,iT)*D1(2,iT)+D01(3,iT)*D1(3,iT))*PAO(iT)
    end do

  case (4)
    do iT=1,nt
      tmp1 = tmp1+(D01(1,iT)*D1(1,iT)+D01(2,iT)*D1(2,iT)+D01(3,iT)*D1(3,iT)+D01(4,iT)*D1(4,iT))*PAO(iT)
    end do

  case (5)
    do iT=1,nt
      tmp1 = tmp1+(D01(1,iT)*D1(1,iT)+D01(2,iT)*D1(2,iT)+D01(3,iT)*D1(3,iT)+D01(4,iT)*D1(4,iT)+D01(5,iT)*D1(5,iT))*PAO(iT)
    end do

  case (6)
    do iT=1,nt
      tmp1 = tmp1+(D01(1,iT)*D1(1,iT)+D01(2,iT)*D1(2,iT)+D01(3,iT)*D1(3,iT)+D01(4,iT)*D1(4,iT)+D01(5,iT)*D1(5,iT)+ &
                   D01(6,iT)*D1(6,iT))*PAO(iT)
    end do

  case (7)
    do iT=1,nt
      tmp1 = tmp1+(D01(1,iT)*D1(1,iT)+D01(2,iT)*D1(2,iT)+D01(3,iT)*D1(3,iT)+D01(4,iT)*D1(4,iT)+D01(5,iT)*D1(5,iT)+ &
                   D01(6,iT)*D1(6,iT)+D01(7,iT)*D1(7,iT))*PAO(iT)
    end do

  case (8)
    do iT=1,nt
      tmp1 = tmp1+(D01(1,iT)*D1(1,iT)+D01(2,iT)*D1(2,iT)+D01(3,iT)*D1(3,iT)+D01(4,iT)*D1(4,iT)+D01(5,iT)*D1(5,iT)+ &
                   D01(6,iT)*D1(6,iT)+D01(7,iT)*D1(7,iT)+D01(8,iT)*D1(8,iT))*PAO(iT)
    end do

  case (9)
    do iT=1,nt
      tmp1 = tmp1+(D01(1,iT)*D1(1,iT)+D01(2,iT)*D1(2,iT)+D01(3,iT)*D1(3,iT)+D01(4,iT)*D1(4,iT)+D01(5,iT)*D1(5,iT)+ &
                   D01(6,iT)*D1(6,iT)+D01(7,iT)*D1(7,iT)+D01(8,iT)*D1(8,iT)+D01(9,iT)*D1(9,iT))*PAO(iT)
    end do

  case default
    do iT=1,nt
      do iRys=1,nRys
        tmp1 = tmp1+D01(iRys,iT)*PAO(iT)*D1(iRys,iT)
      end do
    end do

end select

tmp1_ = tmp1_+tmp1

end subroutine ass1a
