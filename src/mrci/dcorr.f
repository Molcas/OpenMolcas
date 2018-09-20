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
      SUBROUTINE DCORR(JREFX,AREF,ICSPCK,INTSYM,INDX,DMO)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION AREF(*),JREFX(*),ICSPCK(*),
     *          INTSYM(*),INDX(*),DMO(*)
*
      JO(L)=ICUNP(ICSPCK,L)
*
C CORRECTION TO DENSITY MATRIX IN ACPF CASE.
      IF(IPRINT.GE.7) WRITE(6,*)' ENP IN DENS =',ENP
      FAC=1.0D00-(1.0D00/ENP)
      IAD27=0
      CALL dDAFILE(Lu_27,2,AREF,NREF,IAD27)
      IK=0
      DO 40 INDA=1,IRC(1)
        II1=(INDA-1)*LN
        IF(JREFX(INDA).NE.0) THEN
          IK=IK+1
          TSUM=AREF(IK)*AREF(IK)*FAC
          IJ=0
          DO 110 I=1,LN
            IOC=(1+JO(II1+I))/2
            IJ=IJ+I
            DMO(IJ)=DMO(IJ)+IOC*TSUM
110       CONTINUE
        END IF
40    CONTINUE
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer_array(INDX)
      END
