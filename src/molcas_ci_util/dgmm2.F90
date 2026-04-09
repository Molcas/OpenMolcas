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

!#ifdef _DEBUGPRINT_
subroutine DGMM2(A,DIAG,IWAY,NRDIM,NCDIM)
! PRODUCT OF DIAGONAL MATRIX AND MATRIX, IN PLACE :
!
! IWAY = 1 : A(I,J) = DIAG(I)*A(I,J)
! IWAY = 2 : A(I,J) = DIAG(J)*A(I,J)

use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: IWAY, NRDIM, NCDIM
real(kind=wp), intent(inout) :: A(NRDIM,NCDIM)
real(kind=wp), intent(in) :: DIAG(*)
integer(kind=iwp) :: I, J

if (IWAY == 1) then
  do J=1,NCDIM
    do I=1,NRDIM
      A(I,J) = A(I,J)*DIAG(I)
    end do
  end do
end if

if (IWAY == 2) then
  do J=1,NCDIM
    A(:,J) = A(:,J)*DIAG(J)
  end do
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' DIAG A  FROM DGMM2'
if (IWAY == 1) then
  call WRTMAT(DIAG,1,NRDIM,1,NRDIM)
else
  call WRTMAT(DIAG,1,NCDIM,1,NCDIM)
end if
call WRTMAT(A,NRDIM,NCDIM,NRDIM,NCDIM)
#endif

return

end subroutine DGMM2
