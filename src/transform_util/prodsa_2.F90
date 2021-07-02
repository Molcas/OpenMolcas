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

subroutine ProdsA_2(AB,iA,iB,CMO,nMO,Y)

implicit real*8(a-h,o-z)
implicit integer(i-n)
real*8 AB(iA,iB), CMO(iB,nMO), Y(iA,nMO)

call DGEMM_('N','N',iA,nMO,iB,1.0d0,AB,iA,CMO,iB,0.0d0,Y,iA)

return

end subroutine ProdsA_2
