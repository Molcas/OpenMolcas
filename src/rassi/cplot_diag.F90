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

subroutine CPLOT_DIAG(MATR,MATI,DIM,EIGVECR,EIGVECI)

use definitions, only: iwp, wp, u6
use constants, only: Zero

implicit none
integer(kind=iwp), intent(in) :: DIM
real(kind=wp), intent(inout) :: MATR(DIM*(DIM+1)/2), MATI(DIM*(DIM+1)/2)
real(kind=wp), intent(out) :: EIGVECR(DIM,DIM), EIGVECI(DIM,DIM)
real(kind=wp) CEIGVAL(DIM)
complex(kind=wp) MATFULL((DIM*(DIM+1)/2))
complex(kind=wp) CEIGVEC(DIM,DIM)
complex(kind=wp) ZWORK(2*DIM-1)
real(kind=wp) RWORK(3*DIM-2)
integer(kind=iwp) INFO, I, J

do J=1,(DIM*(DIM+1)/2)
  MATFULL(J) = cmplx(MATR(J),MATI(J),kind=wp)
end do

call zhpev_('V','U',DIM,MATFULL,CEIGVAL,CEIGVEC,DIM,ZWORK,RWORK,INFO)

if (INFO /= 0) then
  write(u6,*) 'Error in diagonalization'
  write(u6,*) 'INFO: ',INFO
  call ABEND()
end if

do I=1,DIM
  do J=1,DIM
    EIGVECR(I,J) = real(CEIGVEC(I,J))
    EIGVECI(I,J) = aimag(CEIGVEC(I,J))
  end do
end do

call DCOPY_(DIM*(DIM+1)/2,[Zero],0,MATR,1)
call DCOPY_(DIM*(DIM+1)/2,[Zero],0,MATI,1)

do J=1,DIM
  MATR((J*(J-1)/2)+J) = CEIGVAL(J)
end do

end subroutine CPLOT_DIAG
