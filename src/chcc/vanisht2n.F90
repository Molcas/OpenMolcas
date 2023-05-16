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

subroutine VanishT2n(T2n1,T2n2,beSGrp,gaSGrp)
! this routine does:
! vanish space for T2n = T2n(-(+)) (i>(>=)j,(be>(>=)ga)")
!
! parameter description:
! T2nx  - arrays for T2+- (O)
! xSGrp - SubGroups of be,ga (I)

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimSGrpbe, no
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: T2n1(*), T2n2(*)
integer(kind=iwp), intent(in) :: beSGrp, gaSGrp
integer(kind=iwp) :: length1, length2

!1 calc lengths
if (beSGrp == gaSGrp) then
  length1 = nTri_Elem(no)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/2
  length2 = nTri_Elem(no-1)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)-1)/2
else
  length1 = nTri_Elem(no)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
  length2 = nTri_Elem(no-1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
end if

!2 vanish
T2n1(1:length1) = Zero
T2n2(1:length2) = Zero

return

end subroutine VanishT2n
