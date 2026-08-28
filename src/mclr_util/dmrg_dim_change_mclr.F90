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

subroutine dmrg_dim_change_mclr(orbspc,ndim,iflag)

use Index_Functions, only: nTri_Elem
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: orbspc(8), iflag
integer(kind=iwp), intent(out) :: ndim

!write(u6,*) '==================================================='
!write(u6,*) ' Currently, only valid for no symmetry calculation'
!write(u6,*) '==================================================='

! I remember it was already valid for all symmetry.
!                  -- yma, need check it again 2015.5.14

!write(u6,*) 'orbspc',orbspc(1:8)

select case (iflag)
  case (0)
    ndim = sum(orbspc(:))
  case (1)
    ndim = orbspc(1)**2
  case (2)
    ndim = orbspc(1)**4
  case (3)
    ndim = nTri_Elem(orbspc(1))
  case (4)
    ndim = nTri_Elem(orbspc(1)**2)
  case default
    ndim = 0
    write(u6,*) 'unknow iflag'
    call Quit_OnUserError()
end select

end subroutine dmrg_dim_change_mclr
