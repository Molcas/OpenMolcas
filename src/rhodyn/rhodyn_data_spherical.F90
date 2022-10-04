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
! Copyright (C) 2022, Vladislav Kochetov                               *
!***********************************************************************

module rhodyn_data_spherical
!***********************************************************************
! Purpose: declaration of variables used between different subroutines *
!          for the propagation of rho in spherical tensors basis       *
!***********************************************************************
! k_max            max rank of spherical tensors used in propagation
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
private

abstract interface
  subroutine equation_func_sph(time,rho_t,res)
    import :: wp
    real(kind=wp), intent(in) :: time
    complex(kind=wp), intent(in) :: rho_t(:,:,:)
    complex(kind=wp), intent(out) :: res(:,:,:)
  end subroutine equation_func_sph
end interface

integer(kind=iwp) :: k_max, len_sph
integer(kind=iwp), allocatable :: k_ranks(:), q_proj(:)
integer(kind=iwp), allocatable :: list_sf_states(:), list_sf_mult(:), list_so_mult(:), list_so_sf(:)
real(kind=wp), allocatable :: list_sf_spin(:), list_so_spin(:), list_so_proj(:)
complex(kind=wp), allocatable :: V_SO_red(:,:,:), rho_sph_t(:,:,:), midk1(:,:,:), midk2(:,:,:), midk3(:,:,:), midk4(:,:,:)

public :: k_max, k_ranks, q_proj, list_sf_states, list_sf_mult, list_so_mult, list_so_sf, list_sf_spin, &
          list_so_spin, list_so_proj, V_SO_red, len_sph, rho_sph_t, midk1, midk2, midk3, midk4, equation_func_sph

end module rhodyn_data_spherical
