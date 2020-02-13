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
C * Surboutine to project out some nuclear coordinates from the       *
C * velocities and the forces/mass and do dynamics in reduced         *
C * dimensionality.                                                   *
C *                                                                   *
C * 23/01/2020                                                        *
C * Morgane Vacher                                                    *
C *                                                                   *
C *********************************************************************

C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

      SUBROUTINE project_out(vel,force,natom)
      USE Isotopes
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
      INTEGER                                 :: i,j,p
      INTEGER                                 :: natom
      INTEGER                                 :: Iso
      REAL*8, ALLOCATABLE                     :: Mass(:)
      CHARACTER, ALLOCATABLE                  :: atom(:)*2
      REAL*8, DIMENSION(natom), INTENT(INOUT) :: vel,force
      REAL*8, ALLOCATABLE                     :: pcoo(:,:)
      REAL*8                                  :: pvel,pforce,norm
C
      CALL mma_allocate(atom,natom)
      CALL mma_allocate(pcoo,POUT,natom*3)
      CALL mma_allocate(Mass,natom)

      CALL Get_Name_Full(atom)
      CALL Get_dArray('Proj_Coord',pcoo,POUT*natom*3)
      CALL Get_nAtoms_All(matom)
      CALL Get_Mass_All(Mass,matom)

      DO p = 1,POUT
        norm = 0
        DO i=1, natom
          DO j=1, 3
            norm = norm + pcoo(p,3*(i-1)+j)*pcoo(p,3*(i-1)+j)
          ENDDO
        ENDDO

        pvel = 0
        DO i=1, natom
          IF (i.GT.matom) THEN
            CALL LeftAd(atom(i))
            Iso=0
            CALL Isotope(Iso,atom(i),Mass(i))
          END IF
          DO j=1, 3
            pvel = pvel + pcoo(p,3*(i-1)+j)*vel(3*(i-1)+j)*Mass(i)
          ENDDO
        ENDDO
        pvel = pvel / (norm*sum(Mass))
        DO i=1, natom
          DO j=1, 3
            vel(3*(i-1)+j) = vel(3*(i-1)+j) -
     & pvel*pcoo(p,3*(i-1)+j)
          ENDDO
        ENDDO

        pforce = 0
        DO i=1, natom
          IF (i.GT.matom) THEN
            CALL LeftAd(atom(i))
            Iso=0
            CALL Isotope(Iso,atom(i),Mass(i))
          END IF
          DO j=1, 3
            pforce = pforce + pcoo(p,3*(i-1)+j)*force(3*(i-1)+j)
          ENDDO
        ENDDO
        pforce = pforce / (norm*sum(Mass))
        DO i=1, natom
          IF (i.GT.matom) THEN
            CALL LeftAd(atom(i))
            Iso=0
            CALL Isotope(Iso,atom(i),Mass(i))
          END IF
          DO j=1, 3
            force(3*(i-1)+j) = force(3*(i-1)+j) -
     & pforce*Mass(i)*pcoo(p,3*(i-1)+j)
          ENDDO
        ENDDO
      ENDDO

      CALL mma_deallocate(pcoo)
      CALL mma_deallocate(atom)
      CALL mma_deallocate(Mass)

      RETURN
      END
