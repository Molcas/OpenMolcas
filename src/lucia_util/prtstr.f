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
      SUBROUTINE PRTSTR(ISTR,NEL,NSTR)
*
* Print NSTR strings each containing NEL electrons
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ISTR(NEL,NSTR)
*
      DO JSTR = 1, NSTR
        WRITE(6,'(1H0,A,I6,A,4X,10(2X,I4),/,(1H ,19X,10(2X,I4)))' )
     &   ' String ',JSTR,' : ',(ISTR(IEL,JSTR),IEL=1,NEL)
      END DO
*
      RETURN
      END
