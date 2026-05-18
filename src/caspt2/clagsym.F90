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

subroutine CLagSym(nAshT,DG1,DG2,DF1,DF2,mode)

use Constants, only: Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAshT, mode
real(kind=wp), intent(inout) :: DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT), DF1(nAshT,nAshT), DF2(nAshT,nAshT,nAshT,nAshT)
integer(kind=iwp) :: iI, iJ, iK, iL
real(kind=wp) :: Val, Val1, Val2, Val3, Val4

!return
!if (mode == 0) then
do iI=1,nAshT
  do iJ=1,iI-1
    Val1 = DG1(iI,iJ)
    Val2 = DG1(iJ,iI)
    DG1(iI,iJ) = (Val1+Val2)*Half
    DG1(iJ,iI) = (Val1+Val2)*Half
    Val1 = DF1(iI,iJ)
    Val2 = DF1(iJ,iI)
    DF1(iI,iJ) = (Val1+Val2)*Half
    DF1(iJ,iI) = (Val1+Val2)*Half
  end do
end do
!end if

if (mode == 0) then
  !! Follow G2 symmetry
  do iI=1,nAshT
    do iJ=1,nAshT
      do iK=1,nAshT
        do iL=1,nAshT
          Val1 = DG2(iI,iJ,iK,iL)
          Val2 = DG2(iJ,iI,iL,iK)
          Val3 = DG2(iK,iL,iI,iJ)
          Val4 = DG2(iL,iK,iJ,iI)
          Val = (Val1+Val2+Val3+Val4)*Quart
          DG2(iI,iJ,iK,iL) = Val
          DG2(iJ,iI,iL,iK) = Val
          DG2(iK,iL,iI,iJ) = Val
          DG2(iL,iK,iJ,iI) = Val
          Val1 = DF2(iI,iJ,iK,iL)
          Val2 = DF2(iJ,iI,iL,iK)
          Val3 = DF2(iK,iL,iI,iJ)
          Val4 = DF2(iL,iK,iJ,iI)
          Val = (Val1+Val2+Val3+Val4)*Quart
          DF2(iI,iJ,iK,iL) = Val
          DF2(iJ,iI,iL,iK) = Val
          DF2(iK,iL,iI,iJ) = Val
          DF2(iL,iK,iJ,iI) = Val
        end do
      end do
    end do
  end do
else if (mode == 1) then
  !! Follow EtuEyz symmetry
  do iI=1,nAshT
    do iJ=1,nAshT
      do iK=1,nAshT
        do iL=1,nAshT
          Val1 = DG2(iI,iJ,iK,iL)
          Val2 = DG2(iL,iK,iJ,iI)
          Val = (Val1+Val2)*Quart
          ! DG2(iI,iJ,iK,iL) = Val
          ! DG2(iL,iK,iJ,iI) = Val
          !if ((ii /= il) .and. (ij /= ik)) then
          !  DG2(iI,iJ,iK,iL) = Two*val
          !  DG2(iL,iK,iJ,iI) = Zero
          !end if
          Val1 = DF2(iI,iJ,iK,iL)
          Val2 = DF2(iL,iK,iJ,iI)
          Val = (Val1+Val2)*Quart
          ! DF2(iI,iJ,iK,iL) = Val
          ! DF2(iL,iK,iJ,iI) = Val
          !if ((ii /= il) .and. (ij /= ik)) then
          !  DF2(iI,iJ,iK,iL) = Two*val
          !  DF2(iL,iK,iJ,iI) = Zero
          !end if
        end do
      end do
    end do
  end do
end if

return

end subroutine CLagSym
