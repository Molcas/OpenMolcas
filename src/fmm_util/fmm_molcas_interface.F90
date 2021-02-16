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

!subroutine fmm_call_car_to_sph(CarMoms,SphMoms,ndim,LMAX)
!
!  use fmm_global_paras, only: INTK, REALK
!  use fmm_car_to_sph, only: fmm_init_car_to_sph, &
!                            fmm_free_car_to_sph, &
!                            fmm_transform_car_to_sph
!
!  implicit none
!  integer(INTK), intent(in) :: ndim, LMAX
!  real(REALK), intent(in)   :: CarMoms(ndim,(LMAX+1)*(LMAX+2)/2,0:LMAX)
!  real(REALK), intent(out)  :: SphMoms(ndim,2*LMAX+1,0:LMAX)
!
!  call fmm_init_car_to_sph(LMAX)
!  call fmm_transform_car_to_sph(CarMoms,SphMoms,ndim,LMAX)
!  call fmm_free_car_to_sph()
!
!end subroutine fmm_call_car_to_sph

!------------------------------------------------------------------------------

subroutine fmm_call_get_J_matrix(ndim,nBas,dens,Fock)

#include "macros.fh"

use fmm_global_paras, only: INTK, REALK
!use fmm_interface, only: fmm_get_J_matrix

implicit none
integer(INTK), intent(in)  :: ndim, nBas
real(REALK), intent(in)    :: dens(ndim)
real(REALK), intent(inout) :: Fock(ndim)

unused_var(nBas)
unused_var(dens)
unused_var(Fock)

!print *, 'Calling FMM code...'
!call fmm_get_J_matrix(1,dens(1),Fock(1))
!call fmm_get_J_matrix(nbas,dens,Fock)
!print *, 'Returning from FMM code...'

end subroutine fmm_call_get_J_matrix
