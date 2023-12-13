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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Laplace_PRSQ(V,M,N,NDIM)

use ReMez_mod, only: IW
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: M, N, NDIM
real(kind=wp), intent(in) :: V(NDIM,M)
integer(kind=iwp) :: I, IMAX, IMIN, J, NMAX

! ----- PRINT OUT A SQUARE MATRIX -----
! -V- IS -N- ROWS BY -M- COLUMNS, WITH LEADING REAL -NDIM-

NMAX = 10
IMAX = 0
do
  IMIN = IMAX+1
  IMAX = IMAX+NMAX
  if (IMAX > M) IMAX = M
  write(IW,9008)
  write(IW,9028) (I,I=IMIN,IMAX)
  do J=1,N
    write(IW,9048) J,(V(J,I),I=IMIN,IMAX)
  end do
  if (IMAX >= M) exit
end do

return

9008 format(1X)
9028 format(10X,10(4X,I4,4X))
9048 format(I5,1X,10F12.7)

end subroutine Laplace_PRSQ
