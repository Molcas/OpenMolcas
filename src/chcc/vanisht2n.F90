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

implicit none
#include "chcc1.fh"
#include "o2v4.fh"
real*8 T2n1(1)
real*8 T2n2(1)
integer beSGrp, gaSGrp
! help variables
integer length1, length2

!1 calc legths
if (beSGrp == gaSGrp) then
  length1 = no*(no+1)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/4
  length2 = no*(no-1)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)-1)/4
else
  length1 = no*(no+1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)/2
  length2 = no*(no-1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)/2
end if

!2 vanish
call mv0zero(length1,length1,T2n1(1))
call mv0zero(length2,length2,T2n2(1))

return

end subroutine VanishT2n
