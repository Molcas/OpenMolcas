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

subroutine MakeWwd(Ww,W1,aSGrp,beSGrp,gaSGrp)
! this routine does:
! Make Ww(+)((aa)",(bega)") from W1(a",be",a",ga")
! and index bega" is be">=ga"
! for all cases aSGrp >= bSGrp
! N.B. rutinky MakeWwHlpx niesu prilis vymakane @
!
! parameter description:
! Ww     - array for Ww+(-) (O)
! W1     - array for W1(a",be",a",ga") (I)
! xSGrp  - SubGroup of a",be",ga" (I)

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimSGrpa, DimSGrpbe
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: Ww(*)
real(kind=wp), intent(in) :: W1(*)
integer(kind=iwp), intent(in) :: aSGrp, beSGrp, gaSGrp
integer(kind=iwp) :: dima, dimbe, dimbega, dimga

!1 def dimensions

dimbe = DimSGrpbe(beSGrp)
dimga = DimSGrpbe(gaSGrp)
dima = DimSGrpa(aSGrp)

if (beSGrp == gaSGrp) then
  dimbega = nTri_Elem(dimbe)
else
  dimbega = dimbe*dimga
end if

!2 Make Ww matrix

if (beSGrp == gaSGrp) then
  call MakeWwdHlp1(Ww,W1,dima,dimbe,dimbega)
else
  call MakeWwdHlp2(Ww,W1,dima,dimbe,dimga)
end if

return

end subroutine MakeWwd
