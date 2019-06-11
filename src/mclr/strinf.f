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
      DIMENSION ISGSTI(1),ISGSTO(1)
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
      CALL GetMem('MAXMEM ','MAX    ','INTEGER',KFREEL,imax)
      CALL GetMem('MAXMEM ','ALLO','INTEGER',KFREEL,imax)
      DO 20 ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
        CALL WEIGHT_mclr(iWork(KZ(ITYP)),NELEC(ITYP),NORB1,NORB2,NORB3,
     &              MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),
     &              iWork(KFREEL) )
        END IF
   20 CONTINUE
*
**.5 : Number of electrons in RAS1 and RAS3 per string sub type
*
      DO 25 ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
        CALL IEL13(MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),
     &             NELEC(ITYP),NOCTYP(ITYP),iwork(KEL1(ITYP)),
     &             iWork(KEL3(ITYP)),iwork(KEL123(ITYP)),
     &             iWork(KACTP(ITYP)) )
        END IF
   25 CONTINUE
*
**.6 : Number of strings per type and symmetry for a given string type
*
      DO 30 ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
        CALL NSTRSO_MCLR(NELEC(ITYP),NORB1,NORB2,NORB3,
     &              MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),
     &              iWork(KFREEL),NACOB,iwork(KNSTSO(ITYP)),
     &              NOCTYP(ITYP),NSMST,ITYP,IPRNT)
*. Corresponding offset array
        CALL ZBASE(iWork(KNSTSO(ITYP)),iwork(KISTSO(ITYP)),
     &             NSMST*NOCTYP(ITYP) )
*. Symmetry and class index for each string
         CALL ZSMCL(NSMST,NOCTYP(ITYP),iWork(KNSTSO(ITYP)),
     &        iwork(KSTSM(ITYP)),iWork(KSTCL(ITYP)) )
        END IF
   30 CONTINUE
*
**.7 Construct strings, ordered according to symmetry and class
*
      DO 40 ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
        CALL GENSTR_MCLR(NELEC(ITYP),MNRS1(ITYP),MXRS1(ITYP),
     &              MNRS3(ITYP),MXRS3(ITYP),iWork(KISTSO(ITYP)),
     &              NOCTYP(ITYP),NSMST,iwork(KZ(ITYP)),iwork(KFREEL),
     &              iWork(KSTREO(ITYP)),iWork(KOCSTR(ITYP)),
     &              iwork(KFREEL+NOCTYP(ITYP)*NSMST),ITYP,IPRNT)
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
          CALL ANNSTR(iWork(KOCSTR(ITYP)),NSTFTP(ITYP),NSTFTP(JTYP),
     &                NELEC(ITYP),NACOB,iWork(KZ(JTYP)),
     &                iwork(KSTREO(JTYP)),LROW,0,ISGSTI,ISGSTO,
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
          CALL CRESTR(iwork(KOCSTR(ITYP)),NSTFTP(ITYP),NSTFTP(JTYP),
     &                NELEC(ITYP),NACOB,iwork(KZ(JTYP)),
     &                iwork(KSTREO(JTYP)),0,ISGSTI,ISGSTO,
     &                iwork(KSTSTM(ITYP,1)),iwork(KSTSTM(ITYP,2)),
     &                iWORK(KSTSTMN(ITYP)),iWORK(KSTSTMI(ITYP)),
     &                LROW,JTYP,IPRNT)
        END IF
        END IF
   60 CONTINUE
      CALL GetMem('MAXMEM ','FREE','INTEGER',KFREEL,imax)
      RETURN
      END
