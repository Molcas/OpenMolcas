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
! Copyright (C) 2020, Morgane Vacher                                   *
!***********************************************************************
!
! *********************************************************************
! *                                                                   *
! * Surboutine to project out some nuclear coordinates from the       *
! * velocities and do dynamics in reduced dimensionality.             *
! *                                                                   *
! * 03/04/2020                                                        *
! * Morgane Vacher                                                    *
! *                                                                   *
! *********************************************************************

!   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

      SUBROUTINE project_out_vel(vel,natom)
      IMPLICIT REAL*8 (a-h,o-z)
#include "prgm.fh"
#include "warnings.fh"
#include "Molcas.fh"
      PARAMETER    (ROUTINE='VV_Second')
#include "MD.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
      INTEGER                                   :: i,j,p
      INTEGER                                   :: natom
      REAL*8, ALLOCATABLE                       :: Mass(:)
      REAL*8, DIMENSION(natom*3), INTENT(INOUT) :: vel
      REAL*8, DIMENSION(natom*3)                :: vel_m
      REAL*8, ALLOCATABLE                       :: pcoo(:,:),pcoo_m(:,:)
      REAL*8                                    :: pvel
!
      CALL mma_allocate(pcoo,POUT,natom*3)
      CALL mma_allocate(pcoo_m,POUT,natom*3)
      CALL mma_allocate(Mass,natom)

      CALL Get_dArray('Proj_Coord',pcoo,POUT*natom*3)
      CALL GetMassDx(Mass,natom)

! Mass-weight the velocity vector
      DO i=1, natom
        DO j=1, 3
          vel_m(3*(i-1)+j) = vel(3*(i-1)+j)*sqrt(Mass(i))
        ENDDO
      ENDDO

      DO p = 1,POUT
! Mass-weight the projection vector
        DO i=1, natom
          DO j=1, 3
            pcoo_m(p,3*(i-1)+j) = pcoo(p,3*(i-1)+j)*sqrt(Mass(i))
          ENDDO
        ENDDO
! normalise it (needed or not?)
        pcoo_m(p,:) = pcoo_m(p,:)/                                      &
     &    sqrt(dot_product(pcoo_m(p,:),pcoo_m(p,:)))
! Project out
        pvel = dot_product(pcoo_m(p,:),vel_m)
        IF (pvel.GT.0.000001) THEN
          WRITE(6,'(5X,A,6X,D19.12)') 'Proj comp from velo:',pvel
        ENDIF
        vel_m = vel_m - pvel*pcoo_m(p,:)
      ENDDO

! Un-Mass-weight the velocity vector
      DO i=1, natom
        DO j=1, 3
          vel(3*(i-1)+j) = vel_m(3*(i-1)+j)/sqrt(Mass(i))
        ENDDO
      ENDDO

      CALL mma_deallocate(pcoo)
      CALL mma_deallocate(pcoo_m)
      CALL mma_deallocate(Mass)

      RETURN
      END
