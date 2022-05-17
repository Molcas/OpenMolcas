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

subroutine ass1b(D1,PAO,tmp1_,nt,nrys)

implicit real*8(a-h,o-z)
dimension D1(nrys,nt), PAO(nt)

tmp1 = 0.0d0

select case (nRys)

  case (1)
    do it=1,nt
      tmp1 = tmp1+PAO(iT)*D1(1,iT)
    end do

  case (2)
    do it=1,nt
      tmp1 = tmp1+(D1(1,iT)+D1(2,iT))*PAO(It)
    end do

  case (3)
    do it=1,nt
      tmp1 = tmp1+(D1(1,iT)+D1(2,iT)+D1(3,iT))*PAO(It)
    end do

  case (4)
    do it=1,nt
      tmp1 = tmp1+(D1(1,iT)+D1(2,iT)+D1(3,iT)+D1(4,iT))*PAO(It)
    end do

  case (5)
    do it=1,nt
      tmp1 = tmp1+(D1(1,iT)+D1(2,iT)+D1(3,iT)+D1(4,iT)+D1(5,iT))*PAO(It)
    end do

  case (6)
    do it=1,nt
      tmp1 = tmp1+(D1(1,iT)+D1(2,iT)+D1(3,iT)+D1(4,iT)+D1(5,iT)+D1(6,iT))*PAO(It)
    end do

  case (7)
    do it=1,nt
      tmp1 = tmp1+(D1(1,iT)+D1(2,iT)+D1(3,iT)+D1(4,iT)+D1(5,iT)+D1(6,iT)+D1(7,iT))*PAO(It)
    end do

  case (8)
    do it=1,nt
      tmp1 = tmp1+(D1(1,iT)+D1(2,iT)+D1(3,iT)+D1(4,iT)+D1(5,iT)+D1(6,iT)+D1(7,iT)+D1(8,iT))*PAO(It)
    end do

  case (9)
    do it=1,nt
      tmp1 = tmp1+(D1(1,iT)+D1(2,iT)+D1(3,iT)+D1(4,iT)+D1(5,iT)+D1(6,iT)+D1(7,iT)+D1(8,iT)+D1(9,iT))*PAO(It)
    end do

  case default
    do iT=1,nT
      do iRys=1,nRys
        tmp1 = tmp1+PAO(It)*D1(iRys,iT)
      end do
    end do

end select

tmp1_ = tmp1_+tmp1

end subroutine ass1b
