************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2020, Morgane Vacher                                   *
************************************************************************
C
C *********************************************************************
C *                                                                   *
C * Surboutine to keep in only some nuclear coordinates for the       *
C * velocities and do dynamics in reduced dimensionality.             *
C *                                                                   *
C * 27/05/2020                                                        *
C * Morgane Vacher                                                    *
C *                                                                   *
C *********************************************************************

C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

      SUBROUTINE project_in_vel(vel,natom)
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
      REAL*8, DIMENSION(natom*3)                :: vel_m, newvel_m
      REAL*8, ALLOCATABLE                       :: pcoo(:,:),pcoo_m(:,:)
      REAL*8                                    :: pvel
C
      CALL mma_allocate(pcoo,PIN,natom*3)
      CALL mma_allocate(pcoo_m,PIN,natom*3)
      CALL mma_allocate(Mass,natom)

      CALL Get_dArray('Keep_Coord',pcoo,PIN*natom*3)
      CALL GetMassDx(Mass,natom)

      newvel_m = 0.0d0
C Mass-weight the velocity vector
      DO i=1, natom
        DO j=1, 3
          vel_m(3*(i-1)+j) = vel(3*(i-1)+j)*sqrt(Mass(i))
        ENDDO
      ENDDO

      DO p = 1,PIN
C Mass-weight the projection vector
        DO i=1, natom
          DO j=1, 3
            pcoo_m(p,3*(i-1)+j) = pcoo(p,3*(i-1)+j)*sqrt(Mass(i))
          ENDDO
        ENDDO
C normalise it (needed or not?)
        pcoo_m(p,:) = pcoo_m(p,:)/
     &    sqrt(dot_product(pcoo_m(p,:),pcoo_m(p,:)))
C Calculate the projection to keep in
        pvel = dot_product(pcoo_m(p,:),vel_m)
        newvel_m = newvel_m + pvel*pcoo_m(p,:)
      ENDDO

C Un-Mass-weight the velocity vector
      DO i=1, natom
        DO j=1, 3
          vel(3*(i-1)+j) = newvel_m(3*(i-1)+j)/sqrt(Mass(i))
        ENDDO
      ENDDO

      CALL mma_deallocate(pcoo)
      CALL mma_deallocate(pcoo_m)
      CALL mma_deallocate(Mass)

      RETURN
      END
