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
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp, u6
      USE SUPERINDEX, only: MIGEJ, MIGTJ, MAGEB, MAGTB
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDBMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM, NINDEP, NISUP, NASUP, NIES, EPSI,
     &                         NIGEJES, NSES, NIGTJES, EPSE, NSSH,
     &                         NISH, NIGEJ, NAGEBES, NAGEBES, NAGEB,
     &                         NIGTJ, NAGTBES, NAGTB
      IMPLICIT NONE

      real(kind=wp) Dummy(1), EDIAG
      real(kind=wp), ALLOCATABLE:: BD(:), ID(:)
      integer(kind=iwp) ICASE, ISYM, I2, I2ABS, IA, IAABS, IAB, IABQ,
     &                  IBABS, IDID, II, IIABS, IIJ, IIJQ, IIQ, IIS,
     &                  IJABS, ISYMA, ISYMAB, ISYMI, ISYMIJ, NAS, NIN,
     &                  NIS

C Set up non-active diagonal elements of H0.


      DO ICASE=1,13
        DO  ISYM=1,NSYM

          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) CYCLE
          NIS=NISUP(ISYM,ICASE)
          NAS=NASUP(ISYM,ICASE)
          IF(ICASE.GT.11)CALL mma_allocate(BD,NAS,LABEL='BD')
          CALL mma_allocate(ID,NIS,LABEL='ID')


          SELECT CASE (ICASE)
          CASE (1)
C VJTU CASE:
             DO IIS=1,NIS
               IIQ=IIS+NIES(ISYM)
               ID(IIS)= -EPSI(IIQ)
             END DO

          CASE (2)
C VJTIP CASE:
             DO IIS=1,NIS
               IIQ=IIS+NIGEJES(ISYM)
               IIABS=MIGEJ(1,IIQ)
               IJABS=MIGEJ(2,IIQ)
               ID(IIS)= -EPSI(IIABS)-EPSI(IJABS)
             END DO

          CASE (3)
C VJTIM CASE:
             DO IIS=1,NIS
               IIQ=IIS+NIGTJES(ISYM)
               IIABS=MIGTJ(1,IIQ)
               IJABS=MIGTJ(2,IIQ)
               ID(IIS)= -EPSI(IIABS)-EPSI(IJABS)
             END DO

          CASE (4)
C ATVX  CASE:
             DO IIS=1,NIS
                IIQ=IIS+NSES(ISYM)
                ID(IIS)= +EPSE(IIQ)
             END DO

          CASE (5)
C AIVX  CASE:
             IIS=0
             DO ISYMA=1,NSYM
               ISYMI=Mul(ISYMA,ISYM)
               DO IA=1,NSSH(ISYMA)
                 IAABS=IA+NSES(ISYMA)
                 DO II=1,NISH(ISYMI)
                   IIABS=II+NIES(ISYMI)
                   IIS=IIS+1
                   EDIAG= -EPSI(IIABS)+EPSE(IAABS)
                   ID(IIS)= EDIAG
                 END DO
               END DO
             END DO

          CASE (6)
C VJAIP CASE:
             IIS=0
             DO ISYMA=1,NSYM
               ISYMIJ=Mul(ISYMA,ISYM)
               DO I2=1,NIGEJ(ISYMIJ)
                 I2ABS=I2+NIGEJES(ISYMIJ)
                 IIABS=MIGEJ(1,I2ABS)
                 IJABS=MIGEJ(2,I2ABS)
                 DO IA=1,NSSH(ISYMA)
                   IAABS=IA+NSES(ISYMA)
                   IIS=IIS+1
                   EDIAG= -EPSI(IIABS)-EPSI(IJABS)+EPSE(IAABS)
                   ID(IIS)= EDIAG
                 END DO
               END DO
             END DO

          CASE (7)
C VJAIM CASE:
             IIS=0
             DO ISYMA=1,NSYM
               ISYMIJ=Mul(ISYMA,ISYM)
               DO I2=1,NIGTJ(ISYMIJ)
                 I2ABS=I2+NIGTJES(ISYMIJ)
                 IIABS=MIGTJ(1,I2ABS)
                 IJABS=MIGTJ(2,I2ABS)
                 DO IA=1,NSSH(ISYMA)
                   IAABS=IA+NSES(ISYMA)
                   IIS=IIS+1
                   EDIAG= -EPSI(IIABS)-EPSI(IJABS)+EPSE(IAABS)
                   ID(IIS)= EDIAG
                 END DO
               END DO
             END DO

          CASE (8)
C BVATP CASE:
             DO IIS=1,NIS
               IIQ=IIS+NAGEBES(ISYM)
               IAABS=MAGEB(1,IIQ)
               IBABS=MAGEB(2,IIQ)
               ID(IIS)= +EPSE(IAABS)+EPSE(IBABS)
             END DO

          CASE (9)
C BVATM CASE:
             DO IIS=1,NIS
               IIQ=IIS+NAGTBES(ISYM)
               IAABS=MAGTB(1,IIQ)
               IBABS=MAGTB(2,IIQ)
               ID(IIS)= +EPSE(IAABS)+EPSE(IBABS)
             END DO

          CASE (10)
C BJATP CASE:
             IIS=0
             DO ISYMI=1,NSYM
               ISYMAB=Mul(ISYMI,ISYM)
               DO I2=1,NAGEB(ISYMAB)
                 I2ABS=I2+NAGEBES(ISYMAB)
                 IAABS=MAGEB(1,I2ABS)
                 IBABS=MAGEB(2,I2ABS)
                 DO II=1,NISH(ISYMI)
                   IIABS=II+NIES(ISYMI)
                   IIS=IIS+1
                   EDIAG= -EPSI(IIABS)+EPSE(IAABS)+EPSE(IBABS)
                   ID(IIS)= EDIAG
                 END DO
               END DO
             END DO

          CASE (11)
C BJATM CASE:
             IIS=0
             DO ISYMI=1,NSYM
               ISYMAB=Mul(ISYMI,ISYM)
               DO I2=1,NAGTB(ISYMAB)
                 I2ABS=I2+NAGTBES(ISYMAB)
                 IAABS=MAGTB(1,I2ABS)
                 IBABS=MAGTB(2,I2ABS)
                 DO II=1,NISH(ISYMI)
                   IIABS=II+NIES(ISYMI)
                   IIS=IIS+1
                   EDIAG= -EPSI(IIABS)+EPSE(IAABS)+EPSE(IBABS)
                   ID(IIS)= EDIAG
                 END DO
               END DO
             END DO

          CASE (12)
C BJAIP CASE:
             DO IAB=1,NAGEB(ISYM)
               IABQ=IAB+NAGEBES(ISYM)
               IAABS=MAGEB(1,IABQ)
               IBABS=MAGEB(2,IABQ)
               IIS=IIS+1
               EDIAG= EPSE(IAABS)+EPSE(IBABS)
               BD(IAB)= EDIAG
             END DO
             DO IIJ=1,NIGEJ(ISYM)
               IIJQ=IIJ+NIGEJES(ISYM)
               IIABS=MIGEJ(1,IIJQ)
               IJABS=MIGEJ(2,IIJQ)
               EDIAG= -EPSI(IIABS)-EPSI(IJABS)
               ID(IIJ)= EDIAG
             END DO

          CASE (13)
C BJAIM CASE:
             DO IAB=1,NAGTB(ISYM)
               IABQ=IAB+NAGTBES(ISYM)
               IAABS=MAGTB(1,IABQ)
               IBABS=MAGTB(2,IABQ)
               IIS=IIS+1
               EDIAG= EPSE(IAABS)+EPSE(IBABS)
               BD(IAB)= EDIAG
             END DO
             DO IIJ=1,NIGTJ(ISYM)
               IIJQ=IIJ+NIGTJES(ISYM)
               IIABS=MIGTJ(1,IIJQ)
               IJABS=MIGTJ(2,IIJQ)
               EDIAG= -EPSI(IIABS)-EPSI(IJABS)
               ID(IIJ)= EDIAG
             END DO

          CASE DEFAULT
             WRITE(u6,*)'NADIAG: illegal case number'
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

        END DO
      END DO

      END SUBROUTINE NADIAG
