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

subroutine OrdExp(nExp,rExp,nCntrc,Cff)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nExp, nCntrc
real(kind=wp), intent(inout) :: rExp(nExp), Cff(nExp,nCntrc)
integer(kind=iwp) :: iExp, jExp, kExp
real(kind=wp) :: Exp1, Exp2

! Order exponents
! Make the subsequent change in the contraction matrix

do iExp=1,nExp-1
  Exp1 = rExp(iExp)
  kExp = iExp
  do jExp=iExp+1,nExp
    Exp2 = rExp(jExp)
    if (Exp2 > Exp1) then
      Exp1 = Exp2
      kExp = jExp
    end if
  end do
  if (kExp /= iExp) then
    call DSwap_(1,rExp(iExp),1,rExp(kExp),1)
    call DSwap_(nCntrc,Cff(iExp,1),nExp,Cff(kExp,1),nExp)
  end if
end do

return

end subroutine OrdExp
