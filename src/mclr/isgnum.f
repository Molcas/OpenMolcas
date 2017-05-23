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
      FUNCTION ISGNUM2
     &         (NLEV,NVERT,MIDLEV,MIDV1,MIDV2,NMIDV,MXUP,MXDWN,
     &          IDOWN,IUP,IDAW,IRAW,IUSGN,ILSGN,IWALK)
C
C     PURPOSE: FOR ANY GIVEN WALK (STEP VECTOR) COMPUTE THE
C              LEXICAL NUMBER IN THE SPLIT GUGA REPRESENTATION
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION IDOWN(NVERT,0:3),IUP(NVERT,0:3)
      DIMENSION IDAW(NVERT,0:4),IRAW(NVERT,0:4)
      DIMENSION IUSGN(MXUP,NMIDV),ILSGN(MXDWN,NMIDV)
      DIMENSION IWALK(NLEV)
C
C
C     FIND THE MIDVERTEX AND THE COMBINED WALK SYMMETRY
C
      MIDV=1
      DO 100 LEV=NLEV,(MIDLEV+1),-1
        ICASE=IWALK(LEV)
        MIDV=IDOWN(MIDV,ICASE)
100   CONTINUE
      MIDV=MIDV-MIDV1+1
C
C     FIND REVERSE ARC WEIGHT FOR THE UPPER WALK
C
      IRAWSUM=1
      LV=1
      DO 110 LEV=NLEV,(MIDLEV+1),-1
        IC=IWALK(LEV)
        LV=IDOWN(LV,IC)
        IRAWSUM=IRAWSUM+IRAW(LV,IC)
110   CONTINUE
      IUW=IUSGN(IRAWSUM,MIDV)
C
C     FIND DIRECT ARC WEIGHT FOR THE LOWER WALK
C
      IDAWSUM=1
      LV=NVERT
      DO 120 LEV=1,MIDLEV
        IC=IWALK(LEV)
        LV=IUP(LV,IC)
        IDAWSUM=IDAWSUM+IDAW(LV,IC)
120   CONTINUE
      ILW=ILSGN(IDAWSUM,MIDV)
C
C     COMPUTE LEXICAL ORDERING NUMBER
C
      ISGNUM2=IUW+ILW
C
C
C     EXIT
C
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(MIDV2)
      END
