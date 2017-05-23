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
      FUNCTION GTH1ES(IREOTS,IPNT,H,IBSO,MXPNGAS,
     &           IBTSOB,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,IJSM)
*
* one electron integral between orbitals (iorb,itp,ism,jorb,jsm,jtp)
* correct combination of row and column symmetry is assumed
* IJSM = 1 => Lower triangular packed
*      else=> No triangular packing
*
* Last Revision January 98 (IJSM added )
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER IREOTS(*),IPNT(*),IBTSOB(MXPNGAS,*),IBSO(*)
      INTEGER NACOBS(*)
      DIMENSION H(*)
*
      IABS = IORB+IBTSOB(ITP,ISM)-1
      IREO = IREOTS(IABS)
      JABS = JORB+IBTSOB(JTP,JSM)-1
      JREO = IREOTS(JABS)
C?    write(6,*) ' GTH1ES : IREO JREO ',IREO,JREO
*
C?    write(6,*) ' GTH1ES : IBSO ', IBSO(ISM)
      IJ=-2**30
      IF(IJSM.EQ.1) THEN
        IF(ISM.GT.JSM) THEN
          NI = NACOBS(ISM)
          IJ = IPNT(ISM)-1+(JREO-IBSO(JSM))*NI+IREO-IBSO(ISM)+1
        ELSE IF(ISM.EQ.JSM) THEN
          IJMAX = MAX(IREO-IBSO(ISM)+1,JREO-IBSO(JSM)+1)
          IJMIN = MIN(IREO-IBSO(ISM)+1,JREO-IBSO(JSM)+1)
          IJ = IPNT(ISM)-1+IJMAX*(IJMAX-1)/2+IJMIN
        ELSE IF (ISM.LT.JSM) THEN
          NJ = NACOBS(JSM)
          IJ = IPNT(JSM)-1+(IREO-IBSO(ISM))*NJ+JREO-IBSO(JSM)+1
        END IF
      ELSE
         NI = NACOBS(ISM)
         IJ = IPNT(ISM)-1+(JREO-IBSO(JSM))*NI+IREO-IBSO(ISM)+1
      END IF
*
      GTH1ES = H(IJ)
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' one electron integral '
        WRITE(6,*) ' IORB ITP ISM ',IORB,ITP,ISM
        WRITE(6,*) ' JORB JTP JSM ',JORB,JTP,JSM
        WRITE(6,*) ' IJ and H(IJ) ', IJ,H(IJ)
      END IF
*
      RETURN
      END
