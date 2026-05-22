!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CLagSym(nAshT,DG1,DG2,DF1,DF2, mode)

use Constants, only: Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAshT, mode
real(kind=wp), intent(inout) :: DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT), DF1(nAshT,nAshT), DF2(nAshT,nAshT,nAshT,nAshT)
integer(kind=iwp) :: iI, iJ, iK, iL
real(kind=wp) :: Val

!return
!if (mode == 0) then
do iI=1,nAshT
  do iJ=1,iI-1
    Val = (DG1(iI,iJ)+DG1(iJ,iI))*Half
    DG1(iI,iJ) = Val
    DG1(iJ,iI) = Val
    Val = (DF1(iI,iJ)+DF1(iJ,iI))*Half
    DF1(iI,iJ) = Val
    DF1(iJ,iI) = Val
  end do
end do
!end if

select case (mode)
  case (0)
    !! Follow G2 symmetry
    do iI=1,nAshT
      do iJ=1,nAshT
        do iK=1,nAshT
          do iL=1,nAshT
            Val = (DG2(iI,iJ,iK,iL)+DG2(iJ,iI,iL,iK)+DG2(iK,iL,iI,iJ)+DG2(iL,iK,iJ,iI))*Quart
            DG2(iI,iJ,iK,iL) = Val
            DG2(iJ,iI,iL,iK) = Val
            DG2(iK,iL,iI,iJ) = Val
            DG2(iL,iK,iJ,iI) = Val
            Val = (DF2(iI,iJ,iK,iL)+DF2(iJ,iI,iL,iK)+DF2(iK,iL,iI,iJ)+DF2(iL,iK,iJ,iI))*Quart
            DF2(iI,iJ,iK,iL) = Val
            DF2(iJ,iI,iL,iK) = Val
            DF2(iK,iL,iI,iJ) = Val
            DF2(iL,iK,iJ,iI) = Val
          end do
        end do
      end do
    end do
  case (1)
    !! Follow EtuEyz symmetry
    !do iI=1,nAshT
    !  do iJ=1,nAshT
    !    do iK=1,nAshT
    !      do iL=1,nAshT
    !        Val = (DG2(iI,iJ,iK,iL)+DG2(iL,iK,iJ,iI))*Quart
    !        if ((ii /= il) .and. (ij /= ik)) then
    !          DG2(iI,iJ,iK,iL) = Two*Val
    !          DG2(iL,iK,iJ,iI) = Zero
    !        else
    !          DG2(iI,iJ,iK,iL) = Val
    !          DG2(iL,iK,iJ,iI) = Val
    !        end if
    !        Val = (DF2(iI,iJ,iK,iL)+DF2(iL,iK,iJ,iI))*Quart
    !        if ((ii /= il) .and. (ij /= ik)) then
    !          DF2(iI,iJ,iK,iL) = Two*Val
    !         DF2(iL,iK,iJ,iI) = Zero
    !        else
    !          DF2(iI,iJ,iK,iL) = Val
    !          DF2(iL,iK,iJ,iI) = Val
    !        end if
    !      end do
    !    end do
    !  end do
    !end do
  case (2)
end select

return

end subroutine CLagSym
