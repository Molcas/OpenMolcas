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

subroutine Sp_Mlt(W_In,ne,W_out,nVec,C,nab)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ne, nVec, nab
real(kind=wp), intent(in) :: W_In(ne,nVec), C(ne,nab)
real(kind=wp), intent(out) :: W_Out(nVec,nab)
integer(kind=iwp), parameter :: meMax = 10
integer(kind=iwp) :: iab, iAux(meMax+1), ie, me

do iab=1,nab
  me = 0
  do ie=1,ne
    if (C(ie,iab) /= Zero) then
      me = me+1
      iAux(me) = ie
      if (me > meMax) exit
    end if
  end do

  select case (me)

    case (1)
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)
    case (2)
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)+C(iAux(2),iab)*W_In(iAux(2),:)
    case (3)
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)+C(iAux(2),iab)*W_In(iAux(2),:)+C(iAux(3),iab)*W_In(iAux(3),:)
    case (4)
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)+C(iAux(2),iab)*W_In(iAux(2),:)+C(iAux(3),iab)*W_In(iAux(3),:)+ &
                     C(iAux(4),iab)*W_In(iAux(4),:)
    case (5)
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)+C(iAux(2),iab)*W_In(iAux(2),:)+C(iAux(3),iab)*W_In(iAux(3),:)+ &
                     C(iAux(4),iab)*W_In(iAux(4),:)+C(iAux(5),iab)*W_In(iAux(5),:)
    case (6)
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)+C(iAux(2),iab)*W_In(iAux(2),:)+C(iAux(3),iab)*W_In(iAux(3),:)+ &
                     C(iAux(4),iab)*W_In(iAux(4),:)+C(iAux(5),iab)*W_In(iAux(5),:)+C(iAux(6),iab)*W_In(iAux(6),:)
    case (7)
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)+C(iAux(2),iab)*W_In(iAux(2),:)+C(iAux(3),iab)*W_In(iAux(3),:)+ &
                     C(iAux(4),iab)*W_In(iAux(4),:)+C(iAux(5),iab)*W_In(iAux(5),:)+C(iAux(6),iab)*W_In(iAux(6),:)+ &
                     C(iAux(7),iab)*W_In(iAux(7),:)
    case (8)
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)+C(iAux(2),iab)*W_In(iAux(2),:)+C(iAux(3),iab)*W_In(iAux(3),:)+ &
                     C(iAux(4),iab)*W_In(iAux(4),:)+C(iAux(5),iab)*W_In(iAux(5),:)+C(iAux(6),iab)*W_In(iAux(6),:)+ &
                     C(iAux(7),iab)*W_In(iAux(7),:)+C(iAux(8),iab)*W_In(iAux(8),:)
    case (9)
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)+C(iAux(2),iab)*W_In(iAux(2),:)+C(iAux(3),iab)*W_In(iAux(3),:)+ &
                     C(iAux(4),iab)*W_In(iAux(4),:)+C(iAux(5),iab)*W_In(iAux(5),:)+C(iAux(6),iab)*W_In(iAux(6),:)+ &
                     C(iAux(7),iab)*W_In(iAux(7),:)+C(iAux(8),iab)*W_In(iAux(8),:)+C(iAux(9),iab)*W_In(iAux(9),:)
    case (10)
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)+C(iAux(2),iab)*W_In(iAux(2),:)+C(iAux(3),iab)*W_In(iAux(3),:)+ &
                     C(iAux(4),iab)*W_In(iAux(4),:)+C(iAux(5),iab)*W_In(iAux(5),:)+C(iAux(6),iab)*W_In(iAux(6),:)+ &
                     C(iAux(7),iab)*W_In(iAux(7),:)+C(iAux(8),iab)*W_In(iAux(8),:)+C(iAux(9),iab)*W_In(iAux(9),:)+ &
                     C(iAux(10),iab)*W_In(iAux(10),:)
    case default
      W_Out(:,iab) = C(iAux(1),iab)*W_In(iAux(1),:)+C(iAux(2),iab)*W_In(iAux(2),:)+C(iAux(3),iab)*W_In(iAux(3),:)+ &
                     C(iAux(4),iab)*W_In(iAux(4),:)+C(iAux(5),iab)*W_In(iAux(5),:)+C(iAux(6),iab)*W_In(iAux(6),:)+ &
                     C(iAux(7),iab)*W_In(iAux(7),:)+C(iAux(8),iab)*W_In(iAux(8),:)+C(iAux(9),iab)*W_In(iAux(9),:)+ &
                     C(iAux(10),iab)*W_In(iAux(10),:)
      do ie=iAux(meMax+1),ne
        if (C(ie,iab) /= Zero) W_Out(:,iab) = W_Out(:,iab)+C(ie,iab)*W_In(ie,:)
      end do
  end select
end do

return

end subroutine Sp_Mlt
