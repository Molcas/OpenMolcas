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
* Copyright (C) Jeppe Olsen                                            *
************************************************************************
      SUBROUTINE STRTYP(MS2,NACTEL,MNRS10,MXRS30,IPRNT)
*
* construct input common blocks /STRINP/
* from /LUCINP/ and /ORBINP/
*
      IMPLICIT REAL*8 (A-H,O-Z)
      Logical Reduce_Prt
      External Reduce_Prt
#include "detdim.fh"
#include "orbinp_mclr.fh"
*========
* Output
*========
*   Common block /STRINP/
*
* Jeppe Olsen ,  Dec.24 ,Almaden
*                Last Revision March 31
*
#include "strinp_mclr.fh"
* Where INTXC is internal excitation level 1 => no int exc
*                                          2 => single int exc
*                                          3 => double int exc
* DELTA : number of electrons in string - reference + 5
* ISTTP = 0 => zero  order space
* ISTTP = 1 => reference  space ,no internal excitations
* ISTTP = 2 => reference  space ,single internal excitations
* ISTTP = 3 => reference  space ,single internal excitations
*
*
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
*. Number of alpha and beta electrons
      NAEL = (MS2 + NACTEL ) / 2
      NBEL = (NACTEL - MS2 ) / 2
      IF (NAEL + NBEL .NE. NACTEL ) THEN
         Write (6,*) 'STRTYP: NAEL + NBEL .NE. NACTEL'
         Write (6,*) 'NAEL,NBEL,NACTEL=',NAEL,NBEL,NACTEL
************************************************************************
*     The argument iPL was missing so I inserted this piece inside the
*     stars to calculate it the same way as in the start of mclr.f
*     //Jonas B
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=iPL-1
*                                                                      *
************************************************************************
        Call PrInp_MCLR(iPL)
        Call Abend()
      END IF
* Default on alternatice RAS limits
      MXRS10 = MAX(NACTEL,2*NORB1)
      MNRS30 = 0
* =============================
*. Strings in zero order space
* =============================
*. Type : alpha-strings
      ITYPE = 1
      IAZTP = ITYPE
      NELEC(ITYPE) = NAEL
      MNRS1(ITYPE) = MAX(0,MNRS10-MIN(NBEL,NORB1))
      MXRS1(ITYPE) = MIN(NAEL,NORB1,MXRS10)
      MNRS3(ITYPE) = MAX(0,MNRS30-MIN(NBEL,NORB3))
      MXRS3(ITYPE) = MIN(NAEL,NORB3,MXRS30)

      IZORR(ITYPE) = 1
      ISTTP(ITYPE) = 0
*. Type : single annihilated alphastrings
      IF(NAEL.GE.1) THEN
        ITYPE = ITYPE + 1
        NELEC(ITYPE) = NAEL -1
        MNRS1(ITYPE) = MAX(0,MNRS1(1)-1)
        MXRS1(ITYPE) = MIN(NAEL-1,MXRS1(1))
        MNRS3(ITYPE) = MAX(0,MNRS3(1)-1)
        MXRS3(ITYPE) = MIN(NAEL-1,MXRS3(1))
        IZORR(ITYPE) = 1
        ISTTP(ITYPE) = 0
        IATPM1=ITYPE
      END IF
*. Type : double annihilated alphastrings
      IF(NAEL.GE.2) THEN
        ITYPE = ITYPE + 1
        NELEC(ITYPE) = NAEL -2
        MNRS1(ITYPE) = MAX(0,MNRS1(1)-2)
        MXRS1(ITYPE) = MIN(NAEL-2,MXRS1(1))
        MNRS3(ITYPE) = MAX(0,MNRS3(1)-2)
        MXRS3(ITYPE) = MIN(NAEL-2,MXRS3(1))
        IZORR(ITYPE) = 1
        ISTTP(ITYPE) = 0
        IATPM2=ITYPE
      END IF
*. Type : beta strings
      IF(NAEL.EQ.NBEL) THEN
        IBZTP = IAZTP
        IBTPM1=IATPM1
        IBTPM2=IATPM2
      ELSE
        ITYPE = ITYPE + 1
        IBZTP = ITYPE
        NELEC(ITYPE) = NBEL
        MNRS1(ITYPE) = MAX(0,MNRS10-MIN(NAEL,NORB1))
        MXRS1(ITYPE) = MIN(NBEL,NORB1,MXRS10)
        MNRS3(ITYPE) = MAX(0,MNRS30-MIN(NAEL,NORB3))
        MXRS3(ITYPE) = MIN(NBEL,NORB3,MXRS30)
        IZORR(ITYPE) = 1
        ISTTP(ITYPE) = 0
*. Type : single annihilated betastrings
        IF(NBEL.GE.1) THEN
          ITYPE = ITYPE + 1
          NELEC(ITYPE) = NBEL -1
          MNRS1(ITYPE) = MAX(0,MNRS1(IBZTP)-1)
          MXRS1(ITYPE) = MIN(NBEL-1,MXRS1(IBZTP))
          MNRS3(ITYPE) = MAX(0,MNRS3(IBZTP)-1)
          MXRS3(ITYPE) = MIN(NBEL-1,MXRS3(IBZTP))
          IZORR(ITYPE) = 1
          ISTTP(ITYPE) = 0
          IBTPM1=ITYPE
        END IF
*. Type : double annihilated alphastrings
        IF(NBEL.GE.2) THEN
          ITYPE = ITYPE + 1
          NELEC(ITYPE) = NBEL -2
          MNRS1(ITYPE) = MAX(0,MNRS1(IBZTP)-2)
          MXRS1(ITYPE) = MIN(NBEL-2,MXRS1(IBZTP))
          MNRS3(ITYPE) = MAX(0,MNRS3(IBZTP)-2)
          MXRS3(ITYPE) = MIN(NBEL-2,MXRS3(IBZTP))
          IZORR(ITYPE) = 1
          ISTTP(ITYPE) = 0
          IBTPM2=ITYPE
        END IF
      END IF
*
      NSTTYP = ITYPE
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' Information about string types generated '
        WRITE(6,*) ' ========================================='
        WRITE(6,*)
        WRITE(6,'(A,I3)') ' Number of types generated ', NSTTYP
        WRITE(6,*)
          WRITE(6,'(A)' )
     *    ' ============================================'
        WRITE(6,'(A)' )
     * '  Type  NELEC MNRS1 MXRS1 MNRS3 MXRS3 ISTTP '
          WRITE(6,'(A)' )
     *    ' ============================================'
        DO 100 ITYP = 1, NSTTYP
          WRITE(6,'(7I6)') ITYP,NELEC(ITYP),MNRS1(ITYP),MXRS1(ITYP),
     *                     MNRS3(ITYP),MXRS3(ITYP),ISTTP(ITYP)
100     CONTINUE
        IF(NTEST.GE.2) THEN
          WRITE(6,*) ' IARTP IBRTP '
          CALL IWRTMA(IARTP,3,7,3,10)
          CALL IWRTMA(IBRTP,3,7,3,10)
        END IF
      END IF
*
*EAW
         DO ITYP = 1, NSTTYP
           IUNIQMP(ITYP) = ITYP
           IUNIQTP(ITYP) = ITYP
         END DO
*EAW
*
      RETURN
      END
