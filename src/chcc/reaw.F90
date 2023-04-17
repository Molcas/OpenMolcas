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

subroutine ReaW(W,aSGrp,beSGrp,bSGrp,gaSGrp,LunInt)
! simuluje citanie VVVV integralov

use chcc_global, only: DimSGrpa, DimSGrpbe
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: W(1)
integer(kind=iwp) :: aSGrp, beSGrp, bSGrp, gaSGrp, LunInt
integer(kind=iwp) :: dim_

dim_ = DimSGrpa(aSGrp)*DimSGrpbe(beSGrp)*DimSGrpa(bSGrp)*DimSGrpbe(gaSGrp)

call rea1(LunInt,dim_,W(1))

return

end subroutine ReaW
