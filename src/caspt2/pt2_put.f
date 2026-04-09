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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE PT2_PUT(NSIZE,LAB,VEC)
      use definitions, only: iwp, wp, u6
      use caspt2_global, only: LUDMAT
      use pt2_guga, only: IADR10, cLab10
      IMPLICIT None
      integer(kind=iwp), intent(in):: NSIZE
      CHARACTER(len=*), intent(in):: LAB
      real(kind=wp):: VEC(*)

      CHARACTER(len=8) LAB1
      integer(kind=iwp) I, IAD

      I=9-LEN(LAB)
      IF(I>=1) THEN
        LAB1='        '
        LAB1(I:8)=LAB
      ELSE
        LAB1=LAB(1:8)
      END IF

C FIND DISK ADDRESS:
      DO I=1,64
        IF(CLAB10(I)=='   EMPTY') THEN
          CLAB10(I)=LAB1
          IAD=IADR10(I,1)
          IADR10(I,2)=NSIZE
          CALL DDAFILE(LUDMAT,1,VEC,NSIZE,IAD)
          IF(I<64) IADR10(I+1,1)=IAD
          RETURN
        ELSE IF (CLAB10(I)==LAB1) THEN
          IF(NSIZE.GT.IADR10(I,2)) THEN
             WRITE(u6,*)' ATTEMPT TO INCREASE SIZE OF A FIELD.'
             WRITE(u6,*)' SUBROUTINE PUT FAILS.'
             CALL ABEND()
          End If
          IAD=IADR10(I,1)
          IADR10(I,2)=NSIZE
          CALL DDAFILE(LUDMAT,1,VEC,NSIZE,IAD)
          RETURN
        END IF
      END DO

      WRITE(u6,*)' NO MORE AVAILABLE FIELDS ON FILE DENS.'
      WRITE(u6,*)' SUBROUTINE PUT FAILS.'
      CALL ABEND()


      END SUBROUTINE PT2_PUT
