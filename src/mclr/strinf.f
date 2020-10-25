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
      SUBROUTINE STRINF(IPRNT)
      Use Str_Info
*
* Strings for internal space.
* Information is stored in
* Largest allowed length is MSTINF
*
*.Input
* /LUCINP/,/ORBINP/,/CSM/
*.Output
* /STRINP/,/STINF/,/STRBAS/ and string information in STIN
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
*     (and /LUCINP/ not occuring here )
#include "detdim.fh"
#include "orbinp_mclr.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
*
* ======
* Output
* ======
*
#include "strinp_mclr.fh"
#include "stinf_mclr.fh"
#include "strbas_mclr.fh"
#include "csm.fh"
*
      Integer ISGSTI(1),ISGSTO(1)
      Integer, Allocatable:: KFREEL(:)
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
*
*
**.2 : Number of classes per string type and mappings between
**.    string types (/STINF/)
*
      CALL ZSTINF_MCLR(IPRNT)
*
*
**.3 : Static memory for string information
*
      CALL MEMSTR
*
**.4 :Reverse lexical adresing schemes for each type of string
*
*.First free address
      Call mma_MaxINT(imax)
      CALL mma_allocate(KFREEL,imax,Label='KFREEL')
      DO 20 ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
        CALL WEIGHT_mclr(Str(ITYP)%Z,NELEC(ITYP),NORB1,NORB2,NORB3,
     &                   MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),
     &                   MXRS3(ITYP),KFREEL )
        END IF
   20 CONTINUE
*
**.5 : Number of electrons in RAS1 and RAS3 per string sub type
*
      DO 25 ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
        CALL IEL13(MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),
     &             NELEC(ITYP),NOCTYP(ITYP),Str(ITYP)%EL1,
     &             Str(ITYP)%EL3,Str(ITYP)%EL123,
     &             Str(ITYP)%ACTP)
        END IF
   25 CONTINUE
*
**.6 : Number of strings per type and symmetry for a given string type
*
      DO 30 ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
        CALL NSTRSO_MCLR(NELEC(ITYP),NORB1,NORB2,NORB3,
     &                   MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),
     &                   MXRS3(ITYP),KFREEL,NACOB,Str(ITYP)%NSTSO,
     &                   NOCTYP(ITYP),NSMST,ITYP,IPRNT)
*. Corresponding offset array
        CALL ZBASE(Str(ITYP)%NSTSO,Str(ITYP)%ISTSO,
     &             NSMST*NOCTYP(ITYP) )
*. Symmetry and class index for each string
         CALL ZSMCL(NSMST,NOCTYP(ITYP),Str(ITYP)%NSTSO,
     &              Str(ITYP)%STSM,Str(ITYP)%STCL )
        END IF
   30 CONTINUE
*
**.7 Construct strings, ordered according to symmetry and class
*
      DO 40 ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
        CALL GENSTR_MCLR(NELEC(ITYP),MNRS1(ITYP),MXRS1(ITYP),
     &                   MNRS3(ITYP),MXRS3(ITYP),Str(ITYP)%ISTSO,
     &                   NOCTYP(ITYP),NSMST,Str(ITYP)%Z,KFREEL,
     &                   Str(ITYP)%STREO,Str(ITYP)%OCSTR,
     &                   KFREEL(1+NOCTYP(ITYP)*NSMST),ITYP,IPRNT)
        END IF
   40 CONTINUE
*
*
**.8 Internal annihilation arrays between types of strings
*
      DO 50 ITYP = 1, NSTTYP
        IF(IUNIQMP(ITYP).EQ.ITYP) THEN
        IF(ISTAC(ITYP,1).NE.0) THEN
          JTYP = ISTAC(ITYP,1)
          IF(IPRNT.GE.2) THEN
          WRITE(6,*) ' Annihilator arrays between types ',ITYP,JTYP
          WRITE(6,*) ' ==========================================='
          END IF
          IF(ISTAC(ITYP,2).EQ.0) THEN
            LROW = NELEC(ITYP)
          ELSE
            LROW = NACOB
          END IF
          Call iCOPY(NSTFTP(ITYP)*LROW,[0],0,iWork(KSTSTM(ITYP,1)),1)
          Call iCOPY(NSTFTP(ITYP)*LROW,[0],0,iWork(KSTSTM(ITYP,2)),1)
          CALL ANNSTR(Str(ITYP)%OCSTR,NSTFTP(ITYP),NSTFTP(JTYP),
     &                NELEC(ITYP),NACOB,Str(JTYP)%Z,
     &                Str(JTYP)%STREO,LROW,0,ISGSTI,ISGSTO,
     &                iwork(KSTSTM(ITYP,1)),iwork(KSTSTM(ITYP,2)),
     &                JTYP,IPRNT)
        END IF
        END IF
   50 CONTINUE
*
** 6 : Creation arrays
*
      DO 60 ITYP = 1, NSTTYP
        IF(IUNIQMP(ITYP).EQ.ITYP) THEN
        IF(ISTAC(ITYP,2).NE.0) THEN
*. Type of creation map
          IF (ISTAC(ITYP,1).EQ.0) THEN
*. Only creation map, compact scheme with offsets
            LROW = -1
          ELSE IF (ISTAC(ITYP,1).NE.0) THEN
*. Both annihilation and creation, use full form
            LROW = NACOB
          END IF
          JTYP = ISTAC(ITYP,2)
          IF(IPRNT.GE.2) THEN
          WRITE(6,*) ' Creator  arrays between types ',ITYP,JTYP
          WRITE(6,*) ' ==========================================='
          END IF
          CALL CRESTR(Str(ITYP)%OCSTR,NSTFTP(ITYP),NSTFTP(JTYP),
     &                NELEC(ITYP),NACOB,Str(JTYP)%Z,
     &                Str(JTYP)%STREO,0,ISGSTI,ISGSTO,
     &                iwork(KSTSTM(ITYP,1)),iwork(KSTSTM(ITYP,2)),
     &                Str(ITYP)%STSTMN,Str(ITYP)%STSTMI,
     &                LROW,JTYP,IPRNT)
        END IF
        END IF
   60 CONTINUE
      CALL mma_deallocate(KFREEL)
      RETURN
      END
