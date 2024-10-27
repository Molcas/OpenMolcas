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
      SUBROUTINE NADIAG()
      USE SUPERINDEX
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)

#include "caspt2.fh"

      Real*8 Dummy(1)
      Real*8, ALLOCATABLE:: BD(:), ID(:)

C Set up non-active diagonal elements of H0.


      DO 3000 ICASE=1,13
        DO 3001 ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) CYCLE
          NIS=NISUP(ISYM,ICASE)
          NAS=NASUP(ISYM,ICASE)
          IF(ICASE.GT.11)CALL mma_allocate(BD,NAS,LABEL='BD')
          CALL mma_allocate(ID,NIS,LABEL='ID')


          SELECT CASE (ICASE)
          CASE (1)
C VJTU CASE:
             DO 110 IIS=1,NIS
               IIQ=IIS+NIES(ISYM)
               ID(IIS)= -EPSI(IIQ)
 110         CONTINUE

          CASE (2)
C VJTIP CASE:
             DO 210 IIS=1,NIS
               IIQ=IIS+NIGEJES(ISYM)
               IIABS=MIGEJ(1,IIQ)
               IJABS=MIGEJ(2,IIQ)
               ID(IIS)= -EPSI(IIABS)-EPSI(IJABS)
 210         CONTINUE

          CASE (3)
C VJTIM CASE:
             DO 310 IIS=1,NIS
               IIQ=IIS+NIGTJES(ISYM)
               IIABS=MIGTJ(1,IIQ)
               IJABS=MIGTJ(2,IIQ)
               ID(IIS)= -EPSI(IIABS)-EPSI(IJABS)
 310         CONTINUE

          CASE (4)
C ATVX  CASE:
             DO 410 IIS=1,NIS
                IIQ=IIS+NSES(ISYM)
                ID(IIS)= +EPSE(IIQ)
 410          CONTINUE

          CASE (5)
C AIVX  CASE:
             IIS=0
             DO 510 ISYMA=1,NSYM
               ISYMI=MUL(ISYMA,ISYM)
               DO 511 IA=1,NSSH(ISYMA)
                 IAABS=IA+NSES(ISYMA)
                 DO 512 II=1,NISH(ISYMI)
                   IIABS=II+NIES(ISYMI)
                   IIS=IIS+1
                   EDIAG= -EPSI(IIABS)+EPSE(IAABS)
                   ID(IIS)= EDIAG
 512             CONTINUE
 511           CONTINUE
 510         CONTINUE

          CASE (6)
C VJAIP CASE:
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
                   ID(IIS)= EDIAG
 612             CONTINUE
 611           CONTINUE
 610         CONTINUE

          CASE (7)
C VJAIM CASE:
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
                   ID(IIS)= EDIAG
 712             CONTINUE
 711           CONTINUE
 710         CONTINUE

          CASE (8)
C BVATP CASE:
             DO 810 IIS=1,NIS
               IIQ=IIS+NAGEBES(ISYM)
               IAABS=MAGEB(1,IIQ)
               IBABS=MAGEB(2,IIQ)
               ID(IIS)= +EPSE(IAABS)+EPSE(IBABS)
 810         CONTINUE

          CASE (9)
C BVATM CASE:
             DO 910 IIS=1,NIS
               IIQ=IIS+NAGTBES(ISYM)
               IAABS=MAGTB(1,IIQ)
               IBABS=MAGTB(2,IIQ)
               ID(IIS)= +EPSE(IAABS)+EPSE(IBABS)
 910         CONTINUE

          CASE (10)
C BJATP CASE:
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
                   ID(IIS)= EDIAG
 1012            CONTINUE
 1011          CONTINUE
 1010        CONTINUE

          CASE (11)
C BJATM CASE:
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
                   ID(IIS)= EDIAG
 1112            CONTINUE
 1111          CONTINUE
 1110        CONTINUE

          CASE (12)
C BJAIP CASE:
             DO 1210 IAB=1,NAGEB(ISYM)
               IABQ=IAB+NAGEBES(ISYM)
               IAABS=MAGEB(1,IABQ)
               IBABS=MAGEB(2,IABQ)
               IIS=IIS+1
               EDIAG= EPSE(IAABS)+EPSE(IBABS)
               BD(IAB)= EDIAG
 1210        CONTINUE
             DO 1220 IIJ=1,NIGEJ(ISYM)
               IIJQ=IIJ+NIGEJES(ISYM)
               IIABS=MIGEJ(1,IIJQ)
               IJABS=MIGEJ(2,IIJQ)
               EDIAG= -EPSI(IIABS)-EPSI(IJABS)
               ID(IIJ)= EDIAG
 1220        CONTINUE

          CASE (13)
C BJAIM CASE:
             DO 1310 IAB=1,NAGTB(ISYM)
               IABQ=IAB+NAGTBES(ISYM)
               IAABS=MAGTB(1,IABQ)
               IBABS=MAGTB(2,IABQ)
               IIS=IIS+1
               EDIAG= EPSE(IAABS)+EPSE(IBABS)
               BD(IAB)= EDIAG
 1310        CONTINUE
             DO 1320 IIJ=1,NIGTJ(ISYM)
               IIJQ=IIJ+NIGTJES(ISYM)
               IIABS=MIGTJ(1,IIJQ)
               IJABS=MIGTJ(2,IIJQ)
               EDIAG= -EPSI(IIABS)-EPSI(IJABS)
               ID(IIJ)= EDIAG
 1320        CONTINUE

          CASE DEFAULT
             WRITE(6,*)'NADIAG: illegal case number'
             CALL ABEND()
          END SELECT

C NOTE: BDIAG elements used in cases 12 & 13.
      IDID=IDBMAT(ISYM,ICASE)
      IF(ICASE.GT.11) THEN
        CALL DDAFILE(LUSBT,1,BD,NAS,IDID)
        CALL mma_deallocate(BD)
      ELSE
C Dummy read the BDIAG elements. NOTE: NAS, not NIN.
        CALL DDAFILE(LUSBT,0,Dummy,NAS,IDID)
      END IF
      CALL DDAFILE(LUSBT,1,ID,NIS,IDID)
      CALL mma_deallocate(ID)

 3001   CONTINUE
 3000 CONTINUE

      END SUBROUTINE NADIAG
