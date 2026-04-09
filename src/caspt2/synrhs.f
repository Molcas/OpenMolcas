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
      use definitions, only: iwp
      use caspt2_module, only: nCases, nSym, nInDep, nASup, nISup
      IMPLICIT None
      integer(kind=iwp), intent(in):: iVec, iType

      integer(kind=iwp) iCase, iSym, NAS, NIS, NW

***************************************************************
      DO ICASE=1,NCASES
       DO ISYM=1,NSYM
        IF(NINDEP(ISYM,ICASE).EQ.0) Cycle
        IF (ITYPE.EQ.0) THEN
          NAS=NINDEP(ISYM,ICASE)
        ELSE
          NAS=NASUP(ISYM,ICASE)
        END IF
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NW.EQ.0) Cycle
        CALL DRA2SOLV (NAS,NIS,iCASE,iSYM,iVEC)
       END DO
      END DO

      END SUBROUTINE PCOLLVEC

#if 0
      SUBROUTINE PDISTVEC(IVEC,iTYPE)
      use definitions, only: iwp
      use caspt2_module, only: nCases, nSym, nInDep, nASup, nISup
      IMPLICIT None
      integer(kind=iwp), intent(in):: iVec, iType

      integer(kind=iwp) iCase, iSym, NAS, NIS, NW

***************************************************************
      DO ICASE=1,NCASES
       DO ISYM=1,NSYM
        IF(NINDEP(ISYM,ICASE).EQ.0) Cycle
        IF (ITYPE.EQ.0) THEN
          NAS=NINDEP(ISYM,ICASE)
        ELSE
          NAS=NASUP(ISYM,ICASE)
        END IF
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NW.EQ.0) Cycle
        CALL SOLV2DRA (NAS,NIS,iCASE,iSYM,iVEC)
       END DO
      END DO

      END SUBROUTINE PDISTVEC
#endif
