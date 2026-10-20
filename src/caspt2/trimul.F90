!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine TRIMUL(N,M,ALPHA,ASYM,X,LDX,Y,LDY)
! Multiply symmetric matrix ASYM with matrix X.
! Scale result with ALPHA and add it to matrix Y.

use Index_Functions, only: nTri_Elem
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, M, LDX, LDY
real(kind=wp), intent(in) :: ALPHA, ASYM(nTri_Elem(N)), X(LDX,M)
real(kind=wp), intent(inout) :: Y(LDY,M)
integer(kind=iwp) :: I

do I=1,M
  !call DSLMX(N,ALPHA,ASYM,X(:,I),1,Y(:,I),1)
  call DSPMV_('U',N,ALPHA,ASYM,X(:,I),1,One,Y(:,I),1)
end do

end subroutine TRIMUL
