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
      FUNCTION GTH1ES_MCLR(IREOTS,IPNT,H,IBSO,
     &           IBTSOB,IORB,ITP,ISM,JORB,JTP,JSM)
* ireots
* ipnt
* H
* ibso
* ibtsob
*
* This is a fantastic solution of the problem hehehehe!!
*  EAW
*
* one electron integral between orbitals (iorb,itp,ism,jorb,jsm,jtp)
*
* correct combination of row and column symmetry is assumed
      IMPLICIT REAL*8(A-H,O-Z)

#include "Input.fh"
*.Input
      INTEGER IREOTS(*),IPNT(*),IBTSOB(3,*),IBSO(*)
      DIMENSION H(*)
*
      IABS = IORB+IBTSOB(ITP,ISM)-1
      IREO = IREOTS(IABS)
      JABS = JORB+IBTSOB(JTP,JSM)-1
      JREO = IREOTS(JABS)
      I1=IREO-IBSO(ISM)+1+nISH(ISM)
      J1=JREO-IBSO(JSM)+1+nISH(jSM)
      IJ=IPNT(ISM)-1+(J1-1)*NORB(ISM)+I1
      GTH1ES_MCLR = H(IJ)
      RETURN
      END
