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
* Copyright (C) 2006, Igor Schapiro                                    *
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

      SUBROUTINE project_out(vel,force)
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
      EXTERNAL     IsFreeUnit
      INTEGER      i,j,p,natom
      INTEGER      Iso
      LOGICAL      hybrid
      CHARACTER    filname*80
      REAL*8, ALLOCATABLE ::       Mass(:)
      CHARACTER, ALLOCATABLE ::    atom(:)*2
      REAL*8, INTENT(INOUT) ::       vel(:),force(:)
      REAL*8, ALLOCATABLE ::     pcoo(:,:),pvel(:),pforce(:)
C
      CALL mma_allocate(pcoo,POUT,natom*3)
      CALL mma_allocate(pvel,POUT)
      CALL mma_allocate(pforce,POUT)
      CALL mma_allocate(Mass,natom)

      filname = 'comqum.dat'
      CALL F_INQUIRE(filname,hybrid)
      IF (hybrid) THEN
         WRITE(6,'(/,5X,A)') 'Perform QM/MM Molecular Dynamics'
         CALL DxRdNAtomHbrd(natom)
      ELSE
         CALL DxRdNAtomStnd(natom)
      END IF
      CALL Get_dArray('Proj_Coord',pcoo,POUT*natom*3)
      CALL Get_nAtoms_All(matom)
      CALL Get_Mass_All(Mass,matom)

      DO p = 1,POUT
        pvel(p) = dot_product(pcoo(p,:),vel)
     & / dot_product(pcoo(p,:),pcoo(p,:))
        vel(:) = vel(:) - pvel(p)*pcoo(p,:)
        pforce(p) = 0
        DO i=1, natom
          IF (i.GT.matom) THEN
            CALL LeftAd(atom(i))
            Iso=0
            CALL Isotope(Iso,atom(i),Mass(i))
          END IF
          DO j=1, 3
            pforce(p) = pforce(p) + pcoo(p,3*(i-1)+j)*force(3*(i-1)+j)
     & /Mass(i)
          ENDDO
        ENDDO
        pforce(p) = pforce(p) / dot_product(pcoo(p,:),pcoo(p,:))
        DO i=1, natom
          IF (i.GT.matom) THEN
            CALL LeftAd(atom(i))
            Iso=0
            CALL Isotope(Iso,atom(i),Mass(i))
          END IF
          DO j=1, 3
            force(3*(i-1)+j) = force(3*(i-1)+j) -
     & pforce(p)*Mass(i)*pcoo(p,3*(i-1)+j)
          ENDDO
        ENDDO
      ENDDO

      CALL mma_deallocate(pforce)
      CALL mma_deallocate(pvel)
      CALL mma_deallocate(pcoo)

      RETURN
      END
