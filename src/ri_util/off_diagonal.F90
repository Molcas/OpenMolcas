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

subroutine Off_Diagonal(B1,nB,iB1s,iB1e,B2,iB2s,iB2e)

implicit real*8(a-h,o-z)
real*8 B1(nB,iB1s:iB1e), B2(nB,iB2s:iB2e)

do j=iB2s,iB2e
  do i=iB1s,iB1e
    B1(j,i) = B2(i,j)
  end do
end do

return

end subroutine Off_Diagonal
