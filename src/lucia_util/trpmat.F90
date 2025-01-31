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

subroutine TRPMAT(XIN,NROW,NCOL,XOUT)
! XOUT(I,J) = XIN(J,I)

implicit real*8(A-H,O-Z)
dimension XIN(NROW,NCOL), XOUT(NCOL,NROW)

do IROW=1,NROW
  do ICOL=1,NCOL
    XOUT(ICOL,IROW) = XIN(IROW,ICOL)
  end do
end do

end subroutine TRPMAT
