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
* Copyright (C) 1992,1994, Per Ake Malmqvist                           *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE POLY1(CI,NCI)
      use gugx, only: SGS
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nG1, iAdr10, cLab10
      use definitions, only: iwp, wp
      IMPLICIT NONE
* PER-AAKE MALMQUIST, 92-12-07
* THIS PROGRAM CALCULATES THE 1-EL DENSITY
* MATRIX FOR A CASSCF WAVE FUNCTION.

      INTEGER(kind=iwp), INTENT(IN) :: NCI
      REAL(kind=wp), INTENT(IN) :: CI(NCI)

      Integer(kind=iwp) :: nLev
      Real(kind=wp), Allocatable:: SGM1(:), G1TMP(:)

      nLev = SGS%nLev

      IF(NLEV>0) THEN
        CALL mma_allocate(SGM1 ,NCI,LABEL='SGM1')
        CALL mma_allocate(G1TMP,NG1,LABEL='G1TMP')
        CALL DENS1_RPT2(CI,nCI,SGM1,nCI,G1TMP,nLev)
      END IF

* REINITIALIZE USE OF DMAT.
* The fields IADR10 and CLAB10 are kept in caspt2_module.F90
      IADR10(:,1)=-1
      IADR10(:,2)=0
      CLAB10(:)='   EMPTY'
      IADR10(1,1)=0
* HENCEFORTH, THE CALL PUT(NSIZE,LABEL,ARRAY) WILL ENTER AN
* ARRAY ON LUDMAT AND UPDATE THE TOC.
      IF(NLEV>0) THEN
        CALL PT2_PUT(NG1,' GAMMA1',G1TMP)

        CALL mma_deallocate(SGM1)
        CALL mma_deallocate(G1TMP)
      END IF

      END SUBROUTINE POLY1

