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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine LoadGtuvx(TUVX,Gtuvx)
! ****************************************************************
! Purpose:                                                       *
! Loading TUVX array to a 4-D tensor.                            *
! Copyied from src/molcas_ci_util/david5.f                       *
! ****************************************************************

use rasscf_global, only: NAC, NACPR2
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: TUVX(NACPR2)
real(kind=wp), intent(out) :: Gtuvx(NAC,NAC,NAC,NAC)
integer(kind=iwp) :: it, ituvx, iu, iv, ix, ixmax

ituvx = 0
do it=1,NAC
  do iu=1,it
    do iv=1,it
      ixmax = iv
      if (it == iv) ixmax = iu
      do ix=1,ixmax
        ituvx = ituvx+1
        Gtuvx(it,iu,iv,ix) = TUVX(ituvx)
        Gtuvx(iu,it,iv,ix) = TUVX(ituvx)
        Gtuvx(it,iu,ix,iv) = TUVX(ituvx)
        Gtuvx(iu,it,ix,iv) = TUVX(ituvx)
        Gtuvx(iv,ix,it,iu) = TUVX(ituvx)
        Gtuvx(ix,iv,it,iu) = TUVX(ituvx)
        Gtuvx(iv,ix,iu,it) = TUVX(ituvx)
        Gtuvx(ix,iv,iu,it) = TUVX(ituvx)
      end do
    end do
  end do
end do

return

end subroutine LoadGtuvx
