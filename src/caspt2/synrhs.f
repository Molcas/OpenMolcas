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
      SUBROUTINE PCOLLVEC(IVEC,iTYPE)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"


***************************************************************
      DO ICASE=1,NCASES
       DO ISYM=1,NSYM
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 100
        IF (ITYPE.EQ.0) THEN
          NAS=NINDEP(ISYM,ICASE)
        ELSE
          NAS=NASUP(ISYM,ICASE)
        END IF
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NW.EQ.0) GOTO 100
        CALL DRA2SOLV (NAS,NIS,iCASE,iSYM,iVEC)
 100    CONTINUE
       END DO
      END DO


      RETURN
      END

#if 0
      SUBROUTINE PDISTVEC(IVEC,iTYPE)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"


***************************************************************
      DO ICASE=1,NCASES
       DO ISYM=1,NSYM
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 100
        IF (ITYPE.EQ.0) THEN
          NAS=NINDEP(ISYM,ICASE)
        ELSE
          NAS=NASUP(ISYM,ICASE)
        END IF
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NW.EQ.0) GOTO 100
        CALL SOLV2DRA (NAS,NIS,iCASE,iSYM,iVEC)
 100    CONTINUE
       END DO
      END DO


      RETURN
      END
#endif
