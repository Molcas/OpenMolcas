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
      SUBROUTINE SYMINF_MCLR(NIRREP,IPRNT)
*
* Information about number of symmetries
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* NSMSX : number of symmetries of single excitations
* NSMDX : Number of symmetries of double excitations
* NSMST : Number of symmetries of strings
* NSMCI : NUmber of symmetries of CI spaces
* ITSSX : Total symmetrix single excitation
* ITSDX : Total symmetrix double excitation
*
#include "detdim.fh"
#include "csm.fh"
#include "csmprd.fh"
*
* ADASX : symmetry of orbs i and i => symmetry of a+iaj
* ASXAD : symmetry of orb j and excit a+iaj => symmetry of i
* ADSXA : symmetry of orb i and excit a+iaj => symmetry of j
*
* SXSXDX : Symmetry of two single excitations
*          => symmetry of double  excitation
* SXDXSX : Symmetry of single excitation and double excitation
*          => symmetry of single  excitation
*.
      NSMSX = NIRREP
      NSMDX = NIRREP
      NSMST = NIRREP
      NSMCI = NIRREP
      NSMXT = NIRREP
      ITSSX = 1
      ITSDX = 1
      ITSXT = 1
*
      Call iCopy(MXPOBS*MXPOBS,[0],0,ADASX,1)
      Call iCopy(MXPOBS*2*MXPOBS,[0],0,ADSXA,1)
      Call iCopy(MXPOBS*2*MXPOBS,[0],0,ASXAD,1)
      Call iCopy(2*MXPOBS*2*MXPOBS,[0],0,SXSXDX,1)
      Call iCopy(2*MXPOBS*4*MXPOBS,[0],0,SXDXSX,1)
      DO 10 ISYM=1,8
         DO 20 JSYM=1,8
            IJSYM=1+IEOR(ISYM-1,JSYM-1)
            ADASX(ISYM,JSYM)=IJSYM
            ADSXA(ISYM,JSYM)=IJSYM
            ASXAD(ISYM,JSYM)=IJSYM
            SXSXDX(ISYM,JSYM)=IJSYM
            SXDXSX(ISYM,JSYM)=IJSYM
20       CONTINUE
10    CONTINUE
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(IPRNT)
      END
