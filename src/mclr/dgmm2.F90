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

!#define _DEBUGPRINT_
subroutine DGMM2(AOUT,AIN,DIAG,IWAY,NRDIM,NCDIM)
! PRODUCT OF DIAGONAL MATRIX AND MATRIX :
!
!   IWAY = 1 : AOUT(I,J) = DIAG(I)*AIN(I,J)
!   IWAY = 2 : AOUT(I,J) = DIAG(J)*AIN(I,J)

implicit real*8(A-H,O-Z)
dimension AIN(NRDIM,NCDIM), DIAG(*)
dimension AOUT(NRDIM,NCDIM)

if (IWAY == 1) then
  do J=1,NCDIM
    AOUT(:,J) = DIAG(1:NRDIM)*AIN(:,J)
  end do
else if (IWAY == 2) then
  do J=1,NCDIM
    AOUT(:,J) = DIAG(J)*AIN(:,J)
  end do
end if

#ifdef _DEBUGPRINT_
write(6,*) ' AIN DIAG AOUT  FROM DGMTMT'
call WRTMAT(AIN,NRDIM,NCDIM,NRDIM,NCDIM)
if (IWAY == 1) then
  call WRTMAT(DIAG,1,NRDIM,1,NRDIM)
else
  call WRTMAT(DIAG,1,NCDIM,1,NCDIM)
end if
call WRTMAT(AOUT,NRDIM,NCDIM,NRDIM,NCDIM)
#endif

return

end subroutine DGMM2
