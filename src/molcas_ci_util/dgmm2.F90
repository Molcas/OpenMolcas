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

subroutine DGMM2_MOLCAS(AOUT,AIN,DIAG,IWAY,NRDIM,NCDIM)
! PRODUCT OF DIAGONAL MATRIX AND MATRIX :
!
! IWAY = 1 : AOUT(I,J) = DIAG(I)*AIN(I,J)
! IWAY = 2 : AOUT(I,J) = DIAG(J)*AIN(I,J)

implicit real*8(A-H,O-Z)
dimension AIN(NRDIM,NCDIM), DIAG(*)
dimension AOUT(NRDIM,NCDIM)

if (IWAY == 1) then
  do J=1,NCDIM
    do K=1,NRDIM
      AOUT(K,J) = AIN(K,J)*DIAG(K)
    end do
  end do
end if

if (IWAY == 2) then
  do J=1,NCDIM
    FACTOR = DIAG(J)
    call DCOPY_(NRDIM,AIN(1,J),1,AOUT(1,J),1)
    call DSCAL_(NRDIM,FACTOR,AOUT(1,J),1)
  end do
end if

NTEST = 0
if (NTEST /= 0) then
  write(6,*) ' AIN DIAG AOUT  FROM DGMTMT '
  call WRTMAT(AIN,NRDIM,NCDIM,NRDIM,NCDIM)
  if (IWAY == 1) then
    call WRTMAT(DIAG,1,NRDIM,1,NRDIM)
  else
    call WRTMAT(DIAG,1,NCDIM,1,NCDIM)
  end if
  call WRTMAT(AOUT,NRDIM,NCDIM,NRDIM,NCDIM)
end if

return

end subroutine DGMM2_MOLCAS
