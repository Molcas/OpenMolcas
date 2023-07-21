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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine CD_Tester_Diag(PDM,Diag,n)

implicit none
integer n
real*8 PDM(n*(n+1)/2), Diag(n)
integer i, ii

do i=1,n
  ii = i*(i-3)/2+2*i
  Diag(i) = PDM(ii)
end do

end subroutine CD_Tester_Diag
