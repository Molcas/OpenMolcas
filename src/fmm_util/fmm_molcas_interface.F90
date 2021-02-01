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
!SUBROUTINE fmm_call_car_to_sph(CarMoms,SphMoms,ndim,LMAX)
!
!   USE fmm_global_paras
!   USE fmm_car_to_sph, ONLY: fmm_init_car_to_sph,      &
!                             fmm_free_car_to_sph,      &
!                             fmm_transform_car_to_sph
!
!   IMPLICIT NONE
!   INTEGER(INTK), INTENT(IN)  :: ndim, LMAX
!   REAL(REALK),   INTENT(IN)  :: CarMoms(ndim, (LMAX+1)*(LMAX+2)/2 , 0:LMAX)
!   REAL(REALK),   INTENT(OUT) :: SphMoms(ndim, 2*LMAX+1 , 0:LMAX)
!
!   CALL fmm_init_car_to_sph(LMAX)
!   CALL fmm_transform_car_to_sph(CarMoms,SphMoms,ndim,LMAX)
!   CALL fmm_free_car_to_sph
!
!END SUBROUTINE fmm_call_car_to_sph
!
!------------------------------------------------------------------------------

SUBROUTINE fmm_call_get_J_matrix(ndim,nBas,dens,Fock)

   USE fmm_global_paras
!   USE fmm_interface, ONLY: fmm_get_J_matrix

   IMPLICIT NONE
   INTEGER(INTK), INTENT(IN)    :: ndim, nBas
   REAL(REALK),   INTENT(IN)    :: dens(ndim)
   REAL(REALK),   INTENT(INOUT) :: Fock(ndim)

!   print *, "Calling FMM code..."
!   CALL fmm_get_J_matrix(1,dens(1),Fock(1))
!   CALL fmm_get_J_matrix(nbas,dens,Fock)
!   print *, "Returning from FMM code..."

! Avoid unused argument warnings
   IF (.FALSE.) THEN
      CALL Unused_integer(nBas)
      CALL Unused_real_array(dens)
      CALL Unused_real_array(Fock)
   END IF
END SUBROUTINE fmm_call_get_J_matrix

!------------------------------------------------------------------------------
