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

subroutine GetT2n(T2n1,T2n2,beSGrp,gaSGrp,LunAux)
! this routine does:
! Read T2n = T2n(-(+)) (i>(>=)j,(be>(>=)ga)")
! from Tmp3Name(be",ga")
!
! parameter description:
! T2nx    - arrays for T2+- (O)
! xSGrp   - SubGroups of be,ga (I)
! LunAux  - Lun (I)

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimSGrpbe, Tmp3Name, no
use Definitions, only: wp, iwp

#include "intent.fh"
implicit none
real(kind=wp), intent(_OUT_) :: T2n1(*), T2n2(*)
integer(kind=iwp), intent(in) :: beSGrp, gaSGrp, LunAux
integer(kind=iwp) :: length1, length2
character(len=6) :: LunName

!1 calc lengths
if (beSGrp == gaSGrp) then
  length1 = nTri_Elem(no)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/2
else
  length1 = nTri_Elem(no)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
end if

if (beSGrp == gaSGrp) then
  length2 = nTri_Elem(no-1)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)-1)/2
else
  length2 = nTri_Elem(no-1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
end if

!2 get T2n1, T2n2
LunName = Tmp3Name(beSGrp,gaSGrp)
call GetX(T2n1,length1,LunAux,LunName,1,0)
call GetX(T2n2,length2,LunAux,LunName,0,1)

return

end subroutine GetT2n
