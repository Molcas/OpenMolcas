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

subroutine MxLbls(nInter,Grad,Shift,Lbl)

use Slapaf_Info, only: GrdLbl, StpLbl, GrdMax, StpMax
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nInter
real(kind=wp), intent(in) :: Grad(nInter), Shift(nInter)
character(len=8), intent(in) :: Lbl(nInter)
integer(kind=iwp) :: i

#ifdef _DEBUGPRINT_
call RecPrt('MxLbls:Shift',' ',Shift,nInter,1)
call RecPrt('MxLbls:Grad',' ',Grad,nInter,1)
#endif

GrdMax = Zero
StpMax = Zero
do i=1,nInter
  if (abs(Grad(i)) > abs(GrdMax)) then
    GrdMax = Grad(i)
    GrdLbl = Lbl(i)
  end if
  if (abs(Shift(i)) > abs(StpMax)) then
    StpMax = Shift(i)
    StpLbl = Lbl(i)
  end if
end do
#ifdef _DEBUGPRINT_
write(u6,*) ' Tmp output in MxLbls'
write(u6,*) GrdLbl,' ',GrdMax
write(u6,*) StpLbl,' ',StpMax
#endif

return

end subroutine MxLbls
