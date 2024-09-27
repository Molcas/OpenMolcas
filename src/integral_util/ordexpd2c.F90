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

subroutine OrdExpD2C(nExp,Expn,nCntrc,Cff)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nExp, nCntrc
real(kind=wp), intent(inout) :: Expn(nExp), Cff(nExp,nCntrc)
integer(kind=iwp) :: iExp, jExp, kExp
real(kind=wp) :: Exp1, Exp2
#ifdef _ORDER_BAS_
integer(kind=iwp) :: iCntrc, jCntrc, kCntrc
real(kind=wp) :: Bas1, Bas2
#endif

! Order exponents diffuse to compact
! Make the subsequent change in the contraction matrix

do iExp=1,nExp-1
  Exp1 = Expn(iExp)
  kExp = iExp
  do jExp=iExp+1,nExp
    Exp2 = Expn(jExp)
    if (Exp2 < Exp1) then
      Exp1 = Exp2
      kExp = jExp
    end if
  end do
  if (kExp /= iExp) then
    call DSwap_(1,Expn(iExp),1,Expn(kExp),1)
    call DSwap_(nCntrc,Cff(iExp,1),nExp,Cff(kExp,1),nExp)
  end if
end do

#ifdef _ORDER_BAS_
! Now order the contracted basis functions diffuse to compact

do iCntrc=1,nCntrc-1
  Bas1 = abs(Cff(1,iCntrc))
  kCntrc = iCntrc
  do jCntrc=iCntrc+1,nCntrc
    Bas2 = abs(Cff(1,jCntrc))
    if (Bas2 < Bas1) then
      Bas1 = Bas2
      kCntrc = jCntrc
    end if
  end do
  if (kCntrc /= iCntrc) call DSwap_(nExp,Cff(1,iCntrc),1,Cff(1,kCntrc),1)
end do
#endif

return

end subroutine OrdExpD2C
