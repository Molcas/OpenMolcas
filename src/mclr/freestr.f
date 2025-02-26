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
* Copyright (C) 1990,1994,1995, Jeppe Olsen                            *
************************************************************************
      SUBROUTINE FREESTR()
      Use Str_Info, only: NSTTYP,STR,ITYP_DUMMY,IUNIQMP,INDMAP,INUMAP,
     &                    ISTAC,IUNIQTP
      use stdalloc, only: mma_deallocate
*
* Free pointers for saving information about strings and
* their mappings
*
*========
* Input :
*========
* Number and types of strings defined by /STRINP/
* Symmetry information stored in         /CSM/
* String information stored in           /STINF/
*=========
* Output
*=========
* Pointers stored in common block /STRBAS/
*
* Jeppe Olsen, Winter of 1990
*
* Last Revision, Dec 24 1990 , Almaden
*
* Updated with iuniqtp, dec 11, 1994
* Modified for deallocation, Sept. 25, 2005.
      IMPLICIT None
*
      INTEGER ITYP,IITYP,IIIITEST,IMNEW,JJTYP,IANEQ,ICREQ
*. Start of string information
* =====================================================================
*
* 1 : String information
*
* =====================================================================
*
      DO ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
*.  Offsets for occupation of strings and reordering array
          Call mma_deallocate(Str(ITYP)%OCSTR_Hidden)
          Call mma_deallocate(Str(ITYP)%STREO_Hidden)
*. Symmetry and class of each string
          Call mma_deallocate(Str(ITYP)%STSM_Hidden)
          Call mma_deallocate(Str(ITYP)%STCL_Hidden)
        END IF
        nullify(Str(ITYP)%OCSTR,Str(ITYP)%STREO,Str(ITYP)%STSM,
     &          Str(ITYP)%STCL)
      END DO

*. Number of strings per symmetry and occupation
      DO ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
          Call mma_deallocate(Str(ITYP)%NSTSO_Hidden)
*. Offset of strings per symmetry and occupation
          Call mma_deallocate(Str(ITYP)%ISTSO_Hidden)
*. Number of electrons in RAS1 and RAS3 per sub type, is sub-type active
          Call mma_deallocate(Str(ITYP)%EL1_Hidden)
          Call mma_deallocate(Str(ITYP)%EL3_Hidden)
          Call mma_deallocate(Str(ITYP)%ACTP_Hidden)
CMS: New array introduced according to Jeppes new strinfo representation
          Call mma_deallocate(Str(ITYP)%EL123_Hidden)
**. Lexical addressing of arrays: NB! Not allocated here in Jeppes new version!
          Call mma_deallocate(Str(ITYP)%Z_Hidden)
        ELSE
*. redirect
          IITYP = - IUNIQTP(ITYP)
        END IF
        nullify(Str(ITYP)%NSTSO,Str(ITYP)%ISTSO,Str(ITYP)%EL1,
     &          Str(ITYP)%EL3,Str(ITYP)%ACTP,Str(ITYP)%EL123,
     &          Str(ITYP)%Z)
      END DO

*. Mappings between different string types
      DO ITYP = 1, NSTTYP
          IF(ISTAC(ITYP,2).NE.0.AND.ISTAC(ITYP,1).NE.0) THEN
*.creation on string allowed , use full orbital notation
             Call mma_deallocate(Str(ITYP)%STSTMI)
             Call mma_deallocate(Str(ITYP)%STSTMN)
          ELSE IF(ISTAC(ITYP,1).NE.0.AND.ISTAC(ITYP,2).EQ.0) THEN

*. only annihilation allowed, use compact scheme
             Call mma_deallocate(Str(ITYP)%STSTMI)
             Call mma_deallocate(Str(ITYP)%STSTMN)
CMS: New else block
          ELSE IF (ISTAC(ITYP,1).EQ.0.AND.ISTAC(ITYP,2).NE.0) THEN
*. Only creation allowed, use compact scheme with offsets
*
*. Explicit offsets and lengths
             Call mma_deallocate(Str(ITYP)%STSTMI)
             Call mma_deallocate(Str(ITYP)%STSTMN)
          END IF
*. has this map been constructed before ?
          IIIITEST = 0
          IF(IUNIQTP(ITYP).EQ.ITYP.OR.IIIITEST.EQ.1) THEN
            IMNEW = 1
            IUNIQMP(ITYP) = ITYP
          ELSE
*. check type of previous map
            DO JJTYP = 1, ITYP-1
            IITYP = -IUNIQTP(ITYP)
            IF(ABS(IUNIQTP(JJTYP)).EQ.IITYP.AND.
     &          IUNIQMP(JJTYP).EQ.JJTYP) THEN
            IF((ISTAC(ITYP,1).EQ.0.AND.ISTAC(JJTYP,1).EQ.0).OR.
     &         (ISTAC(ITYP,1).NE.0.AND.ISTAC(JJTYP,1).NE.0.AND.
     &          ABS(IUNIQTP(ISTAC(ITYP,1))).EQ.
     &          ABS(IUNIQTP(ISTAC(JJTYP,1)))) ) THEN
                IANEQ = 1
            ELSE
                IANEQ = 0
            END IF
            IF((ISTAC(ITYP,2).EQ.0.AND.ISTAC(JJTYP,2).EQ.0).OR.
     &         (ISTAC(ITYP,2).NE.0.AND.ISTAC(JJTYP,2).NE.0.AND.
     &          ABS(IUNIQTP(ISTAC(ITYP,2))).EQ.
     &          ABS(IUNIQTP(ISTAC(JJTYP,2)))) ) THEN
                ICREQ = 1
            ELSE
                ICREQ = 0
            END IF
            IF(IANEQ.EQ.1.AND.ICREQ.EQ.1) THEN
              IMNEW = 0
              IUNIQMP(ITYP) = -JJTYP
              GOTO 1211
            END IF
*
            END IF
            END DO
*. Normal exit from DO loop only if no identical map was found
            IMNEW = 1
            IUNIQMP(ITYP) = ITYP
 1211       CONTINUE
          END IF
          IF(IMNEW.EQ.1) Call mma_deallocate(Str(ITYP)%STSTM_Hidden)
          nullify(Str(ITYP)%STSTM)
      END DO
*. Symmetry of conjugated orbitals and orbital excitations
*     COBSM,NIFSJ,IFSJ,IFSJO
!     Call mma_deallocate(COBSM)
!     Call mma_deallocate(NIFSJ)
!     Call mma_deallocate(IFSJ)
!     Call mma_deallocate(IFSJO)
*. Symmetry of excitation connecting  strings of given symmetry
!     Call mma_deallocate(STSTX)
*
**. Up and down mappings of strings containing the same number of electrons
*
      DO 70 ITYP = 1, NSTTYP
       IF(INUMAP(ITYP).NE.0) Call mma_deallocate(Str(ITYP)%NUMAP)
       IF(INDMAP(ITYP).NE.0) Call mma_deallocate(Str(ITYP)%NDMAP)
   70 CONTINUE
*
*
*     Some dummy dallocations
*
      ITYP=ITYP_Dummy
      Call mma_deallocate(Str(ITYP)%NSTSO_Hidden)
      Call mma_deallocate(Str(ITYP)%EL1_Hidden)
      Call mma_deallocate(Str(ITYP)%EL3_Hidden)
      nullify(Str(ITYP)%NSTSO,Str(ITYP)%EL1,Str(ITYP)%EL3)

      END SUBROUTINE FREESTR
