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
* Copyright (C) 1990,1994, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE MEMSTR()
      Use Str_Info
*
*
* Construct pointers for saving information about strings and
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
* Jeppe Olsen , Winter of 1990
*
* Last Revision , Dec 24 1990 , Almaden
*
* Updated with iuniqtp, dec 11, 1994
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "detdim.fh"
#include "stdalloc.fh"
#include "orbinp_mclr.fh"
#include "csm.fh"
      LENGTH=100
*. Start of string information
*
* =====================================================================
*
* 1 : String information
*
* =====================================================================
*
*     Some calls are done with points which are out of bounds. To make
*     this to be strictly secure we add a dummy layer and point to that
*     in the case that we are out of bounds. The code seems, however,
*     not to touch these arrays.
*
      ITYP_Dummy=NSTTYP+1
      Allocate(Str(1:NSTTYP+1))

      IIITEST = 1
      DO ITYP = 1, NSTTYP
        IF (IUNIQTP(ITYP).EQ.ITYP) THEN
        NSTRIN = NUMST3(NELEC(ITYP),NORB1,MNRS1(ITYP),MXRS1(ITYP),
     &                  NORB2,NORB3,MNRS3(ITYP),MXRS3(ITYP) )
        LSTRIN = NSTRIN * NELEC(ITYP)
*.  Offsets for occupation of strings and reordering array
          Call mma_allocate(Str(ITYP)%OCSTR_Hidden,LSTRIN,Label='OCSTR')
          Str(ITYP)%OCSTR => Str(ITYP)%OCSTR_Hidden
          Call mma_allocate(Str(ITYP)%STREO_Hidden,NSTRIN,Label='STREO')
          Str(ITYP)%STREO => Str(ITYP)%STREO_Hidden

*. Symmetry and class of each string
          Call mma_allocate(Str(ITYP)%STSM_Hidden,NSTRIN,Label='STSM')
          Str(ITYP)%STSM  => Str(ITYP)%STSM_Hidden
          Call mma_allocate(Str(ITYP)%STCL_Hidden,NSTRIN,Label='STCL')
          Str(ITYP)%STCL  => Str(ITYP)%STCL_Hidden
        ELSE
          IITYP = - IUNIQTP(ITYP)
          Str(ITYP)%OCSTR => Str(IITYP)%OCSTR_Hidden
          Str(ITYP)%STREO => Str(IITYP)%STREO_Hidden
          Str(ITYP)%STSM  => Str(IITYP)%STSM_Hidden
          Str(ITYP)%STCL  => Str(IITYP)%STCL_Hidden
        END IF
      END DO

*. Number of strings per symmetry and occupation
      DO ITYP = 1, NSTTYP
        IF (IUNIQTP(ITYP).EQ.ITYP) THEN
           Call mma_allocate(Str(ITYP)%NSTSO_Hidden,NOCTYP(ITYP)*NSMST,
     &                       Label='NSTSO')
           Str(ITYP)%NSTSO => Str(ITYP)%NSTSO_Hidden
*. Offset of strings per symmetry and occupation
           Call mma_allocate(Str(ITYP)%ISTSO_Hidden,NOCTYP(ITYP)*NSMST,
     &                       Label='ISTSO')
           Str(ITYP)%ISTSO => Str(ITYP)%ISTSO_Hidden
*. Number of electrons in RAS1 and RAS3 per sub type, is sub-type active
           Call mma_allocate(Str(ITYP)%EL1_Hidden,NOCTYP(ITYP),
     &                       Label='EL1')
           Str(ITYP)%EL1  => Str(ITYP)%EL1_Hidden
           Call mma_allocate(Str(ITYP)%EL3_Hidden,NOCTYP(ITYP),
     &                       Label='EL3')
           Str(ITYP)%EL3  => Str(ITYP)%EL3_Hidden
           Call mma_allocate(Str(ITYP)%ACTP_Hidden,NOCTYP(ITYP),
     &                       Label='ACTP')
           Str(ITYP)%ACTP => Str(ITYP)%ACTP_Hidden
CMS: New array introduced according to Jeppes new strinfo representation
           Call mma_allocate(Str(ITYP)%EL123_Hidden,3*NOCTYP(ITYP),
     &                       Label='EL123')
           Str(ITYP)%EL123=> Str(ITYP)%EL123_Hidden
**. Lexical adressing of arrays: NB! Not allocated here in Jeppes new version!
           Call mma_allocate(Str(ITYP)%Z_Hidden,NACOB*NELEC(ITYP),
     &                       Label='Z')
           Str(ITYP)%Z=> Str(ITYP)%Z_Hidden
        ELSE
*. redirect
          IITYP = - IUNIQTP(ITYP)
          Str(ITYP)%NSTSO => Str(IITYP)%NSTSO_Hidden
          Str(ITYP)%ISTSO => Str(IITYP)%ISTSO_Hidden
          Str(ITYP)%EL1   => Str(IITYP)%EL1_Hidden
          Str(ITYP)%EL3   => Str(IITYP)%EL3_Hidden
          Str(ITYP)%ACTP  => Str(IITYP)%ACTP_Hidden
          Str(ITYP)%EL123 => Str(IITYP)%EL123_Hidden
          Str(ITYP)%Z     => Str(IITYP)%Z_Hidden
        END IF
      END DO

CMS: Introduced according to Jeppes new concept.
CMS: NB! Str(ITYP)%EL123 added to IEL13 parameter list!
CMS: Be aware that IEL13 is also called in STRINF
      DO  ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
        CALL IEL13(MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),
     &             NELEC(ITYP),NOCTYP(ITYP),Str(ITYP)%EL1,
     &             Str(ITYP)%EL3,Str(ITYP)%EL123,
     &             Str(ITYP)%ACTP)
        END IF
      END DO

*. Mappings between different string types
      DO ITYP = 1, NSTTYP
c        write(6,*) nelec(ityp),nstrin
         NSTRIN = NSTFTP(ITYP)

         IF (ISTAC(ITYP,2).NE.0.AND.ISTAC(ITYP,1).NE.0) THEN
*.creation on string allowed , use full orbital notation
            LENGTH = NACOB*NSTRIN
            Call mma_allocate(Str(ITYP)%STSTMI,1,Label='STSTMI')
            Call mma_allocate(Str(ITYP)%STSTMN,1,Label='STSTMN')
         ELSE IF(ISTAC(ITYP,1).NE.0.AND.ISTAC(ITYP,2).EQ.0) THEN

*. only annihilation allowed, use compact scheme
            LENGTH = NELEC(ITYP)*NSTRIN
            Call mma_allocate(Str(ITYP)%STSTMI,1,Label='STSTMI')
            Call mma_allocate(Str(ITYP)%STSTMN,1,Label='STSTMN')
CMS: New else block
          ELSE IF (ISTAC(ITYP,1).EQ.0.AND.ISTAC(ITYP,2).NE.0) THEN
*. Only creation allowed, use compact scheme with offsets
*
             CALL NUMST4_MCLR(NELEC(ITYP),NORB1,MNRS1(ITYP),MXRS1(ITYP),
     &                     NORB2,NORB3,MNRS3(ITYP),MXRS3(ITYP),
     &                     Str(ITYP)%NSTSO)
            LENGTH = NCASTR_MCLR(2,Str(ITYP)%NSTSO,NOCTYP(ITYP),
     &                      ITYP,NOBPT,3,Str(ITYP)%EL123)
*. Explicit offsets and lengths
            Call mma_allocate(Str(ITYP)%STSTMI,NSTRIN,Label='STSTMI')
            Call mma_allocate(Str(ITYP)%STSTMN,NSTRIN,Label='STSTMN')
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

          IF(IMNEW.EQ.1) THEN
            Call mma_allocate(Str(ITYP)%STSTM_Hidden,LENGTH,2,
     &                        Label='STSTM')
            Str(ITYP)%STSTM => Str(ITYP)%STSTM_Hidden
          ELSE
            IITYP = -IUNIQMP(ITYP)
            Str(ITYP)%STSTM => Str(IITYP)%STSTM_Hidden
          END IF

      END DO

*. Symmetry of conjugated orbitals and orbital excitations
*     COBSM,NIFSJ,IFSJ,IFSJO
!     Call mma_allocate(COBSM,NACOB,Label='COBSM')
!     Call mma_allocate(NIFSJ,NACOB*NSMSX,Label='NIFSJ')
!     Call mma_allocate(IFSJ,NACOB**2,Label='IFSJ')
!     Call mma_allocate(IFSJO,NACOB*NSMSX,Label='IFSJO')
*. Symmetry of excitation connecting  strings of given symmetry
!     Call mma_allocate(STSTX,NSMST*NSMST,Label='STSTX')
*
**. Up and down mappings of strings containing the same number of electrons
*
      DO ITYP = 1, NSTTYP
         IF(INUMAP(ITYP).NE.0)
     &     Call mma_allocate(Str(ITYP)%NUMAP,NSTFTP(ITYP),Label='NUMAP')
         IF(INDMAP(ITYP).NE.0)
     &     Call mma_allocate(Str(ITYP)%NDMAP,NSTFTP(ITYP),Label='NDMAP')
      END DO
*
*     Some dummy allocations
*
      ITYP=ITYP_Dummy
      Call mma_allocate(Str(ITYP)%NSTSO_Hidden,1,Label='NSTSO')
      Str(ITYP)%NSTSO => Str(ITYP)%NSTSO_Hidden
      Call mma_allocate(Str(ITYP)%EL1_Hidden,1,Label='EL1')
      Str(ITYP)%EL1  => Str(ITYP)%EL1_Hidden
      Call mma_allocate(Str(ITYP)%EL3_Hidden,1,Label='EL3')
      Str(ITYP)%EL3  => Str(ITYP)%EL3_Hidden

*. Last word of string information
      RETURN
      END
