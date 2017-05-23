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
      FUNCTION IOCTP2_MCLR(STRING,NEL,ITYP)
*
* Obtain occupation type for STRING .
* For forbidden strings a zero is returned
*
* New version allowing general set of strings
*
#include "detdim.fh"
*. Specific input
      INTEGER  STRING(*)
      Logical Reduce_Prt
      External Reduce_Prt
*. General input
#include "strinp_mclr.fh"
#include "orbinp_mclr.fh"
*
      IF(ITYP.LE.0) THEN
        Write (6,*) 'IOCTP2: ITYP.LE.0'
        Write (6,*) 'ITYP=',ITYP
************************************************************************
*     The argument iPL was missing so I inserted this piece inside the
*     stars to calculate it the same way as in the start of mclr.f
*     //Jonas B
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=iPL-1
*                                                                      *
************************************************************************
*                                                                      *
        Call PrInp_MCLR(iPL)
        Call QTrace
        Call Abend()
      END IF
*. Number of electrons in RAS1 and RAS 3
      IEL1 = 0
      IEL3 = 0
      DO 20 IEL = 1,NEL
        IF(STRING(IEL) .LE. NORB1) IEL1 = IEL1 +1
        IF(NORB1+NORB2+1 .LE. STRING(IEL)) IEL3 = IEL3 + 1
20    CONTINUE
*. Type
      IF((IEL1.GE.MNRS1(ITYP).AND.IEL1.LE.MXRS1(ITYP)).AND.
     *   (IEL3.GE.MNRS3(ITYP).AND.IEL3.LE.MXRS3(ITYP))) THEN
          ITYP2 = (MXRS1(ITYP)-IEL1)
     *         * (MXRS3(ITYP)-MNRS3(ITYP)+1 )
     *         + IEL3-MNRS3(ITYP)+1
      ELSE
          ITYP2 = 0
      END IF
      IOCTP2_MCLR = ITYP2
*
      RETURN
      END
