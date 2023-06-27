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

use Slapaf_Parameters, only: GrdLbl, StpLbl, GrdMax, StpMax

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 Shift(nInter), Grad(nInter)
character Lbl(nInter)*8

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
write(6,*) ' Tmp output in MxLbls'
write(6,*) GrdLbl,' ',GrdMax
write(6,*) StpLbl,' ',StpMax
#endif

return

end subroutine MxLbls
