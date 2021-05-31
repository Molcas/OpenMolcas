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
! Copyright (C) 2019, Thomas J. Duignan                                *
!               2021, Rulin Feng                                       *
!***********************************************************************

subroutine merge_mag_ints(nb,jz,lt,ut,dotran)
! Splice together square matrices stored artificially as
! lower triangular matrices.  Result is that upper triangular
! portion of lt is set to ut and ut is then transposed so
! the "upper triangular" portion is in the "lower triangular"
! elements to circumvent Molcas' internal integral handling.

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nb, jz
logical(kind=iwp), intent(in) :: dotran
real(kind=wp), intent(inout) :: lt(jz), ut(jz)
integer(kind=iwp) :: im, jm, km, lm

do im=1,nb
  do jm=im,nb
    km = (im-1)*nb+jm
    lt(km) = ut(km)
  end do
end do

if (dotran) then
  do im=1,nb
    do jm=1,nb
      km = (im-1)*nb+jm
      lm = (jm-1)*nb+im
      ut(km) = lt(lm)
    end do
  end do
else
  ut(:) = lt(:)
end if

end subroutine merge_mag_ints
