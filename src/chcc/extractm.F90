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

subroutine ExtractM(M,L2,cGrp,deGrp,cSGrp,deSGrp)
! this routine does:
! extract M(m,c",de") from L2(m,c',de')
!
! parameter description:
! M     - M file (O)
! L     - L2 file (O)
! xGrp  - c,delta Group (I)
! xSGrp - c,delta Group (I)

use chcc_global, only: DimGrpa, DimSGrpa, DimSGrpbe, GrpaLow, GrpbeLow, nc
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: M(*)
real(kind=wp), intent(in) :: L2(*)
integer(kind=iwp), intent(in) :: cGrp, deGrp, cSGrp, deSGrp
integer(kind=iwp) :: depp, i, incL20, k, lenMCpp, posL20, posM0

!1 Initial settings

!1.1 def length of mc" block (=increment for posM0)
lenMCpp = nc*DimSGrpa(cSGrp)

!1.2 def posM0
posM0 = 1

!1.3 def posL20
posL20 = 1
k = 0
if (deSGrp > Grpbelow(deGrp)) then
  do i=Grpbelow(deGrp),deSGrp-1
    k = k+DimSGrpbe(i)
  end do
end if
posL20 = posL20+nc*DimGrpa(cGrp)*k

k = 0
if (cSGrp > Grpalow(cGrp)) then
  do i=Grpalow(cGrp),cSGrp-1
    k = k+DimSGrpa(i)
  end do
end if
posL20 = posL20+nc*k

!1.4 def increment for posL20
incL20 = dimGrpa(cGrp)*nc

!2 cycle over de"
do depp=1,DimSgrpbe(deSGrp)

  !2.1 copy block of #mc" size
  M(posM0:posM0+lenMCpp-1) = L2(PosL20:PosL20+lenMCpp-1)

  !2.2 upgrade posM0 and posL20 for next use
  posM0 = posM0+lenMCpp
  posL20 = posL20+incL20

end do

return

end subroutine ExtractM
