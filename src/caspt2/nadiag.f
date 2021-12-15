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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE NADIAG
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

#include "SysDef.fh"

C Set up non-active diagonal elements of H0.


      DO 3000 ICASE=1,13
        DO 3001 ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) GOTO 3001
          NIS=NISUP(ISYM,ICASE)
          NAS=NASUP(ISYM,ICASE)
          IF(ICASE.GT.11)
     &        CALL GETMEM('LBD','ALLO','REAL',LBD,NAS)
          CALL GETMEM('LID','ALLO','REAL',LID,NIS)
          GOTO (100,200,300,400,500,600,700,800,900,1000,
     &          1100,1200,1300) ICASE
          WRITE(6,*)'NADIAG: Fall through computed GOTO.'
          CALL ABEND()

C VJTU CASE:
 100  CONTINUE
      DO 110 IIS=1,NIS
        IIQ=IIS+NIES(ISYM)
        WORK(LID-1+IIS)= -EPSI(IIQ)
 110  CONTINUE
      GOTO 2000

C VJTIP CASE:
 200  CONTINUE
      DO 210 IIS=1,NIS
        IIQ=IIS+NIGEJES(ISYM)
        IIABS=MIGEJ(1,IIQ)
        IJABS=MIGEJ(2,IIQ)
        WORK(LID-1+IIS)= -EPSI(IIABS)-EPSI(IJABS)
 210  CONTINUE
      GOTO 2000

C VJTIM CASE:
 300  CONTINUE
      DO 310 IIS=1,NIS
        IIQ=IIS+NIGTJES(ISYM)
        IIABS=MIGTJ(1,IIQ)
        IJABS=MIGTJ(2,IIQ)
        WORK(LID-1+IIS)= -EPSI(IIABS)-EPSI(IJABS)
 310  CONTINUE
      GOTO 2000

C ATVX  CASE:
 400  CONTINUE
      DO 410 IIS=1,NIS
        IIQ=IIS+NSES(ISYM)
        WORK(LID-1+IIS)= +EPSE(IIQ)
 410  CONTINUE
      GOTO 2000

C AIVX  CASE:
 500  CONTINUE
      IIS=0
      DO 510 ISYMA=1,NSYM
        ISYMI=MUL(ISYMA,ISYM)
        DO 511 IA=1,NSSH(ISYMA)
          IAABS=IA+NSES(ISYMA)
          DO 512 II=1,NISH(ISYMI)
            IIABS=II+NIES(ISYMI)
            IIS=IIS+1
            EDIAG= -EPSI(IIABS)+EPSE(IAABS)
            WORK(LID-1+IIS)= EDIAG
 512      CONTINUE
 511    CONTINUE
 510  CONTINUE
      GOTO 2000

C VJAIP CASE:
 600  CONTINUE
      IIS=0
      DO 610 ISYMA=1,NSYM
        ISYMIJ=MUL(ISYMA,ISYM)
        DO 611 I2=1,NIGEJ(ISYMIJ)
          I2ABS=I2+NIGEJES(ISYMIJ)
          IIABS=MIGEJ(1,I2ABS)
          IJABS=MIGEJ(2,I2ABS)
          DO 612 IA=1,NSSH(ISYMA)
            IAABS=IA+NSES(ISYMA)
            IIS=IIS+1
            EDIAG= -EPSI(IIABS)-EPSI(IJABS)+EPSE(IAABS)
            WORK(LID-1+IIS)= EDIAG
 612      CONTINUE
 611    CONTINUE
 610  CONTINUE
      GOTO 2000

C VJAIM CASE:
 700  CONTINUE
      IIS=0
      DO 710 ISYMA=1,NSYM
        ISYMIJ=MUL(ISYMA,ISYM)
        DO 711 I2=1,NIGTJ(ISYMIJ)
          I2ABS=I2+NIGTJES(ISYMIJ)
          IIABS=MIGTJ(1,I2ABS)
          IJABS=MIGTJ(2,I2ABS)
          DO 712 IA=1,NSSH(ISYMA)
            IAABS=IA+NSES(ISYMA)
            IIS=IIS+1
            EDIAG= -EPSI(IIABS)-EPSI(IJABS)+EPSE(IAABS)
            WORK(LID-1+IIS)= EDIAG
 712      CONTINUE
 711    CONTINUE
 710  CONTINUE
      GOTO 2000

C BVATP CASE:
 800  CONTINUE
      DO 810 IIS=1,NIS
        IIQ=IIS+NAGEBES(ISYM)
        IAABS=MAGEB(1,IIQ)
        IBABS=MAGEB(2,IIQ)
        WORK(LID-1+IIS)= +EPSE(IAABS)+EPSE(IBABS)
 810  CONTINUE
      GOTO 2000

C BVATM CASE:
 900  CONTINUE
      DO 910 IIS=1,NIS
        IIQ=IIS+NAGTBES(ISYM)
        IAABS=MAGTB(1,IIQ)
        IBABS=MAGTB(2,IIQ)
        WORK(LID-1+IIS)= +EPSE(IAABS)+EPSE(IBABS)
 910  CONTINUE
      GOTO 2000

C BJATP CASE:
 1000 CONTINUE
      IIS=0
      DO 1010 ISYMI=1,NSYM
        ISYMAB=MUL(ISYMI,ISYM)
        DO 1011 I2=1,NAGEB(ISYMAB)
          I2ABS=I2+NAGEBES(ISYMAB)
          IAABS=MAGEB(1,I2ABS)
          IBABS=MAGEB(2,I2ABS)
          DO 1012 II=1,NISH(ISYMI)
            IIABS=II+NIES(ISYMI)
            IIS=IIS+1
            EDIAG= -EPSI(IIABS)+EPSE(IAABS)+EPSE(IBABS)
            WORK(LID-1+IIS)= EDIAG
 1012     CONTINUE
 1011   CONTINUE
 1010 CONTINUE
      GOTO 2000

C BJATM CASE:
 1100 CONTINUE
      IIS=0
      DO 1110 ISYMI=1,NSYM
        ISYMAB=MUL(ISYMI,ISYM)
        DO 1111 I2=1,NAGTB(ISYMAB)
          I2ABS=I2+NAGTBES(ISYMAB)
          IAABS=MAGTB(1,I2ABS)
          IBABS=MAGTB(2,I2ABS)
          DO 1112 II=1,NISH(ISYMI)
            IIABS=II+NIES(ISYMI)
            IIS=IIS+1
            EDIAG= -EPSI(IIABS)+EPSE(IAABS)+EPSE(IBABS)
            WORK(LID-1+IIS)= EDIAG
 1112     CONTINUE
 1111   CONTINUE
 1110 CONTINUE
      GOTO 2000

C BJAIP CASE:
 1200 CONTINUE
      DO 1210 IAB=1,NAGEB(ISYM)
        IABQ=IAB+NAGEBES(ISYM)
        IAABS=MAGEB(1,IABQ)
        IBABS=MAGEB(2,IABQ)
        IIS=IIS+1
        EDIAG= EPSE(IAABS)+EPSE(IBABS)
        WORK(LBD-1+IAB)= EDIAG
 1210 CONTINUE
      DO 1220 IIJ=1,NIGEJ(ISYM)
        IIJQ=IIJ+NIGEJES(ISYM)
        IIABS=MIGEJ(1,IIJQ)
        IJABS=MIGEJ(2,IIJQ)
        EDIAG= -EPSI(IIABS)-EPSI(IJABS)
        WORK(LID-1+IIJ)= EDIAG
 1220 CONTINUE
      GOTO 2000

C BJAIM CASE:
 1300 CONTINUE
      DO 1310 IAB=1,NAGTB(ISYM)
        IABQ=IAB+NAGTBES(ISYM)
        IAABS=MAGTB(1,IABQ)
        IBABS=MAGTB(2,IABQ)
        IIS=IIS+1
        EDIAG= EPSE(IAABS)+EPSE(IBABS)
        WORK(LBD-1+IAB)= EDIAG
 1310 CONTINUE
      DO 1320 IIJ=1,NIGTJ(ISYM)
        IIJQ=IIJ+NIGTJES(ISYM)
        IIABS=MIGTJ(1,IIJQ)
        IJABS=MIGTJ(2,IIJQ)
        EDIAG= -EPSI(IIABS)-EPSI(IJABS)
        WORK(LID-1+IIJ)= EDIAG
 1320 CONTINUE
      GOTO 2000

 2000 CONTINUE
C NOTE: BDIAG elements used in cases 12 & 13.
      IDID=IDBMAT(ISYM,ICASE)
      IF(ICASE.GT.11) THEN
        CALL DDAFILE(LUSBT,1,WORK(LBD),NAS,IDID)
        CALL GETMEM('LBD','FREE','REAL',LBD,NAS)
      ELSE
C Dummy read the BDIAG elements. NOTE: NAS, not NIN.
        CALL DDAFILE(LUSBT,0,WORK(1),NAS,IDID)
      END IF
      CALL DDAFILE(LUSBT,1,WORK(LID),NIS,IDID)
      CALL GETMEM('LID','FREE','REAL',LID,NIS)

 3001   CONTINUE
 3000 CONTINUE


      RETURN
      END
