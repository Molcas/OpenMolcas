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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE IFOCK(FC,NI,NJ,NK,FINI,II)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FC(*)
      IF(NI.GT.0)RETURN
      IF(NJ.LE.0.OR.NK.LE.0)RETURN
      JKPOS=NJ*(NJ-1)/2+NK
      IF(NK.GT.NJ)JKPOS=NK*(NK-1)/2+NJ
      IF(II.EQ.0)GO TO 10
      FC(JKPOS)=FC(JKPOS)+FINI+FINI
      RETURN
10    FC(JKPOS)=FC(JKPOS)-FINI
      RETURN
      END
