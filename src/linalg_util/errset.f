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
      SUBROUTINE ERRSET(ierno, inoal, inomes, itrace, iusadr, irange)
C
C     DUMMY SUBROUTINE TO REPLACE THE ERRSET ROUTINE
C     AVAILABLE IN VM FORTRAN.
C
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_integer(ierno)
         CALL Unused_integer(inoal)
         CALL Unused_integer(inomes)
         CALL Unused_integer(itrace)
         CALL Unused_integer(iusadr)
         CALL Unused_integer(irange)
      END IF
      END
