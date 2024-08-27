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

implicit none
integer ne, nVec, nab
real*8 W_In(ne,nVec), W_Out(nVec,nab), C(ne,nab)
integer, parameter :: meMax = 10
integer iAux(meMax+1)
integer iab, me, iVec, ie

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
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)
      end do
    case (2)
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)+C(iAux(2),iab)*W_In(iAux(2),iVec)
      end do
    case (3)
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)+C(iAux(2),iab)*W_In(iAux(2),iVec)+C(iAux(3),iab)*W_In(iAux(3),iVec)
      end do
    case (4)
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)+C(iAux(2),iab)*W_In(iAux(2),iVec)+C(iAux(3),iab)*W_In(iAux(3),iVec)+ &
                          C(iAux(4),iab)*W_In(iAux(4),iVec)
      end do
    case (5)
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)+C(iAux(2),iab)*W_In(iAux(2),iVec)+C(iAux(3),iab)*W_In(iAux(3),iVec)+ &
                          C(iAux(4),iab)*W_In(iAux(4),iVec)+C(iAux(5),iab)*W_In(iAux(5),iVec)
      end do
    case (6)
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)+C(iAux(2),iab)*W_In(iAux(2),iVec)+C(iAux(3),iab)*W_In(iAux(3),iVec)+ &
                          C(iAux(4),iab)*W_In(iAux(4),iVec)+C(iAux(5),iab)*W_In(iAux(5),iVec)+C(iAux(6),iab)*W_In(iAux(6),iVec)
      end do
    case (7)
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)+C(iAux(2),iab)*W_In(iAux(2),iVec)+C(iAux(3),iab)*W_In(iAux(3),iVec)+ &
                          C(iAux(4),iab)*W_In(iAux(4),iVec)+C(iAux(5),iab)*W_In(iAux(5),iVec)+C(iAux(6),iab)*W_In(iAux(6),iVec)+ &
                          C(iAux(7),iab)*W_In(iAux(7),iVec)
      end do
    case (8)
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)+C(iAux(2),iab)*W_In(iAux(2),iVec)+C(iAux(3),iab)*W_In(iAux(3),iVec)+ &
                          C(iAux(4),iab)*W_In(iAux(4),iVec)+C(iAux(5),iab)*W_In(iAux(5),iVec)+C(iAux(6),iab)*W_In(iAux(6),iVec)+ &
                          C(iAux(7),iab)*W_In(iAux(7),iVec)+C(iAux(8),iab)*W_In(iAux(8),iVec)
      end do
    case (9)
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)+C(iAux(2),iab)*W_In(iAux(2),iVec)+C(iAux(3),iab)*W_In(iAux(3),iVec)+ &
                          C(iAux(4),iab)*W_In(iAux(4),iVec)+C(iAux(5),iab)*W_In(iAux(5),iVec)+C(iAux(6),iab)*W_In(iAux(6),iVec)+ &
                          C(iAux(7),iab)*W_In(iAux(7),iVec)+C(iAux(8),iab)*W_In(iAux(8),iVec)+C(iAux(9),iab)*W_In(iAux(9),iVec)
      end do
    case (10)
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)+C(iAux(2),iab)*W_In(iAux(2),iVec)+C(iAux(3),iab)*W_In(iAux(3),iVec)+ &
                          C(iAux(4),iab)*W_In(iAux(4),iVec)+C(iAux(5),iab)*W_In(iAux(5),iVec)+C(iAux(6),iab)*W_In(iAux(6),iVec)+ &
                          C(iAux(7),iab)*W_In(iAux(7),iVec)+C(iAux(8),iab)*W_In(iAux(8),iVec)+C(iAux(9),iab)*W_In(iAux(9),iVec)+ &
                          C(iAux(10),iab)*W_In(iAux(10),iVec)
      end do
    case default
      do iVec=1,nVec
        W_Out(iVec,iab) = C(iAux(1),iab)*W_In(iAux(1),iVec)+C(iAux(2),iab)*W_In(iAux(2),iVec)+C(iAux(3),iab)*W_In(iAux(3),iVec)+ &
                          C(iAux(4),iab)*W_In(iAux(4),iVec)+C(iAux(5),iab)*W_In(iAux(5),iVec)+C(iAux(6),iab)*W_In(iAux(6),iVec)+ &
                          C(iAux(7),iab)*W_In(iAux(7),iVec)+C(iAux(8),iab)*W_In(iAux(8),iVec)+C(iAux(9),iab)*W_In(iAux(9),iVec)+ &
                          C(iAux(10),iab)*W_In(iAux(10),iVec)
      end do
      do ie=iAux(meMax+1),ne
        if (C(ie,iab) /= Zero) call DaXpY_(nVec,C(ie,iab),W_In(ie,1),ne,W_Out(1,iab),1)
      end do
  end select
end do

return

end subroutine Sp_Mlt
