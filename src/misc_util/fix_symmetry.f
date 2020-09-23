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
* Copyright (C) 2013, Ignacio Fdez. Galvan                             *
************************************************************************
*  Fix_Symmetry
*
*> @brief
*>   Fix the symmetry of a structure by removing off-zero values
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Fix the symmetry of a structure, by making sure that the coordinates that
*> should be zero remain there, removing numerical inaccuracies.
*>
*> @param[in,out] Coord Cartesian coordinates of the structure to fix (unique atoms)
*> @param[in]     nAt   Number of atoms in the structure
*> @param[in]     Stab  Stabilizers for each atom
************************************************************************
      SUBROUTINE Fix_Symmetry(Coord,nAt,Stab)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER nAt,Stab(nAt)
      REAL*8 Coord(3,nAt),thr
      PARAMETER ( thr = 1.0D-12 )
#include "real.fh"
#include "WrkSpc.fh"
*
      DO iAt=1,nAt
        DO j=0,2
          IF (IAND(Stab(iAt),2**j).GT.0) THEN
            IF (ABS(Coord(j+1,iAt)).GT.thr) THEN
              CALL WarningMessage(1,
     &             'Significant deviation from symmetry axis.')
            END IF
            Coord(j+1,iAt)=Zero
          END IF
        END DO
      END DO
*
      END
