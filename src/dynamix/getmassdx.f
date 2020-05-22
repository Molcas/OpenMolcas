************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE GetMassDx(Mass,natom)
      USE Isotopes
      IMPLICIT NONE
#include "stdalloc.fh"
      REAL*8, INTENT(INOUT) :: Mass(*)
      INTEGER, INTENT(IN) :: natom
      INTEGER :: matom,i,Iso
      CHARACTER, ALLOCATABLE :: atom(:)*2

      CALL mma_allocate(atom,natom)
      CALL Get_Name_Full(atom)
      CALL Get_nAtoms_All(matom)
      CALL Get_Mass_All(Mass,matom)
      DO i=1,natom
        IF (i.GT.matom) THEN
          CALL LeftAd(atom(i))
          Iso=0
          CALL Isotope(Iso,atom(i),Mass(i))
        END IF
      END DO
      CALL mma_deallocate(atom)

      END
