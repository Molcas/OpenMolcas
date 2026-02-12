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

use rasscf_global, only: NACPR2, NAC

implicit none
real*8, dimension(NACPR2) :: TUVX
real*8, dimension(NAC,NAC,NAC,NAC) :: Gtuvx
integer it, iu, iv, ix, ituvx, ixmax
#include "warnings.h"

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
