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
* Copyright (C) 1998, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE GASSPC
*
*
* Divide orbital spaces into
*
*  Inactive spaces : Orbitals that are doubly occupied in all CI spaces
*  Active orbitals : Orbitals that have variable occ in atleast some spaces.
*  Secondary spaces: Orbitals that are unoccupied in all spaces
*
* I_IAD : Division based upon occupation in Compound CI spaces IGSOCC
* I_IADX: Division based upon occupation in First CI space
*
* Jeppe Olsen, Summer of 98 ( not much of an summer !)
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.fh"
#include "cgas.fh"
#include "strinp.fh"
#include "orbinp.fh"
*
* Some dummy initializtions
      NEL_MAX = 0 ! jwk-cleanup
*
      NEL_REF = NELEC(1) + NELEC(2)
*
* For compound space
*
      DO IGAS = 1, NGAS
*
       IF(IGAS.EQ.1) THEN
         NEL_MAX = 2*NGSOBT(IGAS)
       ELSE
         NEL_MAX = NEL_MAX + 2*NGSOBT(IGAS)
       END IF
*
       IF(IGSOCC(IGAS,1) .EQ. NEL_MAX  .AND.
     &    IGSOCC(IGAS,2) .EQ. NEL_MAX       ) THEN
*. Inactive  space
          I_IAD(IGAS) = 1
       ELSE IF(IGAS.GT.1.AND.IGSOCC(IGAS-1,1) .EQ. NEL_REF ) THEN
*. Delete space
          I_IAD(IGAS) = 3
       ELSE
*. Active space
          I_IAD(IGAS) = 2
       END IF
*
      END DO
*
* For First CI space
*
      DO IGAS = 1, NGAS
*
       IF(IGAS.EQ.1) THEN
         NEL_MAX = 2*NGSOBT(IGAS)
       ELSE
         NEL_MAX = NEL_MAX + 2*NGSOBT(IGAS)
       END IF
*
       IF(IGSOCCX(IGAS,1,1) .EQ. NEL_MAX  .AND.
     &    IGSOCCX(IGAS,2,1) .EQ. NEL_MAX       ) THEN
*. Inactive  space
          I_IADX(IGAS) = 1
       ELSE IF(IGAS.GT.1.AND.IGSOCCX(IGAS-1,1,1) .EQ. NEL_REF ) THEN
*. Delete space
          I_IADX(IGAS) = 3
       ELSE
*. Active space
          I_IADX(IGAS) = 2
       END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
     &  ' Division of orbitals according to compound CI space'
        WRITE(6,*)
     &  ' ================================================== '
        WRITE(6,*)
        WRITE(6,*) ' Inactive = 1, Active = 2, Delete = 3 '
        WRITE(6,*)
        CALL IWRTMA(I_IAD,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*)
     &  ' Division of orbitals according to first CI space'
        WRITE(6,*)
     &  ' ================================================== '
        WRITE(6,*)
        WRITE(6,*) ' Inactive = 1, Active = 2, Delete = 3 '
        WRITE(6,*)
        CALL IWRTMA(I_IADX,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
