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
      SUBROUTINE MKLIST(LIST)
      USE SUPERINDEX
      IMPLICIT NONE
C Subroutine for setting up the 17 lists of coupling
C  coefficients -- See sgm.f and sgm.ol for usage.
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"

      INTEGER LIST(*)

      INTEGER IA,IB,II,IJ,IT,IU
      INTEGER IAB,IIJ,ITU,ITU1,ITU2,IUT1,IUT2,IUV,IUU2
      INTEGER ITUV,IUTV,IUVT,IVTU,IVUT
      INTEGER IAQ,IBQ,IIQ,IJQ,ITQ,IUQ,IVQ,IUVQ
      INTEGER ILIST,ISL1,ISL2,ISL3
      INTEGER LADR,LADR1,LADR2,LADR3,LADR4,LADR5,LADR6,LADR7,LADR8,
     &        LADR9,LADR10,LADR11,LADR12,LADR13,LADR14,LADR15,
     &        LADR16,LADR17
      INTEGER NOFF


      LADR=1
      DO ILIST=1,17
       DO ISL1=1,NSYM
        DO ISL3=1,NSYM
         LLIST(ISL1,ISL3,ILIST)=LADR
         LADR=LADR+4*NLIST(ISL1,ISL3,ILIST)
        END DO
       END DO
      END DO

C Lists 1 and 2. TUV/TU*
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        LADR1 =LLIST(ISL1,ISL3,1)
        LADR2 =LLIST(ISL1,ISL3,2)
        ISL2=MUL(ISL1,ISL3)
        DO IT=1,NASH(ISL2)
          ITQ=IT+NAES(ISL2)
          DO IUV=1,NTU(ISL3)
            IUVQ=IUV+NTUES(ISL3)
            IUQ=MTU(1,IUVQ)
            IVQ=MTU(2,IUVQ)
            ITUV=KTUV(ITQ,IUQ,IVQ)-NTUVES(ISL1)
            IUTV=KTUV(IUQ,ITQ,IVQ)-NTUVES(ISL1)
C Add to list 1: ITUV,IT,IUV,V=1
C Add to list 1: IUTV,IT,IUV+NTU(ISL3),V=1
            LIST(LADR1)=ITUV
            LIST(LADR1+1)=IT
            LIST(LADR1+2)=IUV
            LIST(LADR1+3)=1
            LADR1=LADR1+4
            LIST(LADR1)=IUTV
            LIST(LADR1+1)=IT
            LIST(LADR1+2)=IUV+NTU(ISL3)
            LIST(LADR1+3)=1
            LADR1=LADR1+4
            IVUT=KTUV(IVQ,IUQ,ITQ)-NTUVES(ISL1)
C Add to list 2: ITUV,IT,IUV,V=1
C Add to list 2: IVUT,IT,IUV+NTU(ISL3),V=-1
            LIST(LADR2)=ITUV
            LIST(LADR2+1)=IT
            LIST(LADR2+2)=IUV
            LIST(LADR2+3)=1
            LADR2=LADR2+4
            LIST(LADR2)=IVUT
            LIST(LADR2+1)=IT
            LIST(LADR2+2)=IUV+NTU(ISL3)
            LIST(LADR2+3)=2
            LADR2=LADR2+4
          END DO
         END DO
        END DO
       END DO

C Lists 3 and 5. TUV/TU+
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        LADR3 =LLIST(ISL1,ISL3,3)
        LADR5 =LLIST(ISL1,ISL3,5)
        ISL2=MUL(ISL1,ISL3)
        DO IT=1,NASH(ISL2)
          ITQ=IT+NAES(ISL2)
          DO IUV=1,NTGEU(ISL3)
            IUVQ=IUV+NTGEUES(ISL3)
            IUQ=MTGEU(1,IUVQ)
            IVQ=MTGEU(2,IUVQ)
            IUVT=KTUV(IUQ,IVQ,ITQ)-NTUVES(ISL1)
            IUTV=KTUV(IUQ,ITQ,IVQ)-NTUVES(ISL1)
            IF(IUQ.EQ.IVQ) THEN
C Add to list 3: IUVT,IT,IUV,V=2
C Add to list 5: IUTV,IT,IUV,V=2
              LIST(LADR3)=IUVT
              LIST(LADR3+1)=IT
              LIST(LADR3+2)=IUV
              LIST(LADR3+3)=2
              LADR3=LADR3+4
              LIST(LADR5)=IUTV
              LIST(LADR5+1)=IT
              LIST(LADR5+2)=IUV
              LIST(LADR5+3)=2
              LADR5=LADR5+4
            ELSE
              IVUT=KTUV(IVQ,IUQ,ITQ)-NTUVES(ISL1)
              IVTU=KTUV(IVQ,ITQ,IUQ)-NTUVES(ISL1)
C Add to list 3: IUVT,IT,IUV,V=1
C Add to list 3: IVUT,IT,IUV,V=1
C Add to list 5: IUTV,IT,IUV,V=1
C Add to list 5: IVTU,IT,IUV,V=1
              LIST(LADR3)=IUVT
              LIST(LADR3+1)=IT
              LIST(LADR3+2)=IUV
              LIST(LADR3+3)=1
              LADR3=LADR3+4
              LIST(LADR3)=IVUT
              LIST(LADR3+1)=IT
              LIST(LADR3+2)=IUV
              LIST(LADR3+3)=1
              LADR3=LADR3+4
              LIST(LADR5)=IUTV
              LIST(LADR5+1)=IT
              LIST(LADR5+2)=IUV
              LIST(LADR5+3)=1
              LADR5=LADR5+4
              LIST(LADR5)=IVTU
              LIST(LADR5+1)=IT
              LIST(LADR5+2)=IUV
              LIST(LADR5+3)=1
              LADR5=LADR5+4
            END IF
          END DO
         END DO
        END DO
       END DO

C Lists 4 and 6. TUV/TU-
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        LADR4 =LLIST(ISL1,ISL3,4)
        LADR6 =LLIST(ISL1,ISL3,6)
        ISL2=MUL(ISL1,ISL3)
        DO IT=1,NASH(ISL2)
          ITQ=IT+NAES(ISL2)
          DO IUV=1,NTGTU(ISL3)
            IUVQ=IUV+NTGTUES(ISL3)
            IUQ=MTGTU(1,IUVQ)
            IVQ=MTGTU(2,IUVQ)
            IUVT=KTUV(IUQ,IVQ,ITQ)-NTUVES(ISL1)
            IUTV=KTUV(IUQ,ITQ,IVQ)-NTUVES(ISL1)
            IVUT=KTUV(IVQ,IUQ,ITQ)-NTUVES(ISL1)
            IVTU=KTUV(IVQ,ITQ,IUQ)-NTUVES(ISL1)
C Add to list 4: IUVT,IT,IUV,V=1
C Add to list 4: IVUT,IT,IUV,V=-1
C Add to list 6: IUTV,IT,IUV,V=1
C Add to list 6: IVTU,IT,IUV,V=-1
            LIST(LADR4)=IUVT
            LIST(LADR4+1)=IT
            LIST(LADR4+2)=IUV
            LIST(LADR4+3)=1
            LADR4=LADR4+4
            LIST(LADR4)=IVUT
            LIST(LADR4+1)=IT
            LIST(LADR4+2)=IUV
            LIST(LADR4+3)=2
            LADR4=LADR4+4
            LIST(LADR6)=IUTV
            LIST(LADR6+1)=IT
            LIST(LADR6+2)=IUV
            LIST(LADR6+3)=1
            LADR6=LADR6+4
            LIST(LADR6)=IVTU
            LIST(LADR6+1)=IT
            LIST(LADR6+2)=IUV
            LIST(LADR6+3)=2
            LADR6=LADR6+4
          END DO
         END DO
        END DO
       END DO

C Lists 7 and 8. TU*/T
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        LADR7 =LLIST(ISL1,ISL3,7)
        LADR8 =LLIST(ISL1,ISL3,8)
        ISL2=MUL(ISL1,ISL3)
        NOFF=NTU(ISL1)
        DO IT=1,NASH(ISL2)
          ITQ=IT+NAES(ISL2)
          DO IU=1,NASH(ISL3)
            IUQ=IU+NAES(ISL3)
            IUT1=KTU(IUQ,ITQ)-NTUES(ISL1)
            IUT2=IUT1+NOFF
C Add to list 7: IUT1,IT,IU,V=
C Add to list 7: IUT2,IT,IU,V=
            LIST(LADR7)=IUT1
            LIST(LADR7+1)=IT
            LIST(LADR7+2)=IU
            LIST(LADR7+3)=1
            LADR7=LADR7+4
            LIST(LADR7)=IUT2
            LIST(LADR7+1)=IT
            LIST(LADR7+2)=IU
            LIST(LADR7+3)=2
            LADR7=LADR7+4
C Add to list 8: ITU1,IT,IU,V=
C Add to list 8: ITU2,IT,IU,V=
            ITU1=KTU(ITQ,IUQ)-NTUES(ISL1)
            ITU2=ITU1+NOFF
            LIST(LADR8)=ITU1
            LIST(LADR8+1)=IT
            LIST(LADR8+2)=IU
            LIST(LADR8+3)=1
            LADR8=LADR8+4
            LIST(LADR8)=ITU2
            LIST(LADR8+1)=IT
            LIST(LADR8+2)=IU
            LIST(LADR8+3)=2
            LADR8=LADR8+4
          END DO
         END DO
        END DO
       END DO

C Lists  9 and 10. TU+-/T
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        LADR9 =LLIST(ISL1,ISL3,9)
        LADR10=LLIST(ISL1,ISL3,10)
        ISL2=MUL(ISL1,ISL3)
        DO IT=1,NASH(ISL2)
          ITQ=IT+NAES(ISL2)
          DO IU=1,NASH(ISL3)
            IUQ=IU+NAES(ISL3)
            IF(ITQ.GT.IUQ) THEN
              ITU=KTGEU(ITQ,IUQ)-NTGEUES(ISL1)
C Add to list  9: ITU,IT,IU,V= 1/Sqr(2)
              LIST(LADR9)=ITU
              LIST(LADR9+1)=IT
              LIST(LADR9+2)=IU
              LIST(LADR9+3)=1
              LADR9=LADR9+4
              ITU=KTGTU(ITQ,IUQ)-NTGTUES(ISL1)
C Add to list 10: ITU,IT,IU,V= 1/Sqr(6)
              LIST(LADR10)=ITU
              LIST(LADR10+1)=IT
              LIST(LADR10+2)=IU
              LIST(LADR10+3)=1
              LADR10=LADR10+4
            ELSE IF(ITQ.LT.IUQ) THEN
              ITU=KTGEU(IUQ,ITQ)-NTGEUES(ISL1)
C Add to list  9: ITU,IT,IU,V= 1/Sqr(2)
              LIST(LADR9)=ITU
              LIST(LADR9+1)=IT
              LIST(LADR9+2)=IU
              LIST(LADR9+3)=1
              LADR9=LADR9+4
              ITU=KTGTU(IUQ,ITQ)-NTGTUES(ISL1)
C Add to list 10: ITU,IT,IU,V=-1/Sqr(6)
              LIST(LADR10)=ITU
              LIST(LADR10+1)=IT
              LIST(LADR10+2)=IU
              LIST(LADR10+3)=2
              LADR10=LADR10+4
            ELSE
              ITU=KTGEU(ITQ,IUQ)-NTGEUES(ISL1)
C Add to list  9: ITU,IT,IU,V= 1/Sqr(2)
              LIST(LADR9)=ITU
              LIST(LADR9+1)=IT
              LIST(LADR9+2)=IU
              LIST(LADR9+3)=2
              LADR9=LADR9+4
            END IF
          END DO
         END DO
        END DO
       END DO

C Lists 12 and 13. T/TU+-
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        LADR12=LLIST(ISL1,ISL3,12)
        LADR13=LLIST(ISL1,ISL3,13)
        ISL2=MUL(ISL1,ISL3)
        DO IT=1,NASH(ISL1)
          ITQ=IT+NAES(ISL1)
          DO IU=1,NASH(ISL2)
            IUQ=IU+NAES(ISL2)
            IF(ITQ.GT.IUQ) THEN
              ITU=KTGEU(ITQ,IUQ)-NTGEUES(ISL3)
C Add to list 12: IT,IU,ITU,V= 1
              LIST(LADR12)=IT
              LIST(LADR12+1)=IU
              LIST(LADR12+2)=ITU
              LIST(LADR12+3)=1
              LADR12=LADR12+4
              ITU=KTGTU(ITQ,IUQ)-NTGTUES(ISL3)
C Add to list 13: IT,IU,ITU,V= 1
              LIST(LADR13)=IT
              LIST(LADR13+1)=IU
              LIST(LADR13+2)=ITU
              LIST(LADR13+3)=1
              LADR13=LADR13+4
            ELSE IF(ITQ.LT.IUQ) THEN
              ITU=KTGEU(IUQ,ITQ)-NTGEUES(ISL3)
C Add to list 12: IT,IU,ITU,V= 1
              LIST(LADR12)=IT
              LIST(LADR12+1)=IU
              LIST(LADR12+2)=ITU
              LIST(LADR12+3)=1
              LADR12=LADR12+4
              ITU=KTGTU(IUQ,ITQ)-NTGTUES(ISL3)
C Add to list 13: IT,IU,ITU,V=-1
              LIST(LADR13)=IT
              LIST(LADR13+1)=IU
              LIST(LADR13+2)=ITU
              LIST(LADR13+3)=2
              LADR13=LADR13+4
            ELSE
              ITU=KTGEU(ITQ,IUQ)-NTGEUES(ISL3)
C Add to list 12: IT,IU,ITU,V= 2
              LIST(LADR12)=IT
              LIST(LADR12+1)=IU
              LIST(LADR12+2)=ITU
              LIST(LADR12+3)=2
              LADR12=LADR12+4
            END IF
          END DO
         END DO
        END DO
       END DO

C List 11. T/TU*
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        LADR11=LLIST(ISL1,ISL3,11)
        ISL2=MUL(ISL1,ISL3)
        NOFF=NTU(ISL3)
        DO IT=1,NASH(ISL1)
          ITQ=IT+NAES(ISL1)
          DO IU=1,NASH(ISL2)
            IUQ=IU+NAES(ISL2)
            IUT2=NOFF+KTU(IUQ,ITQ)-NTUES(ISL3)
C Add to list 11: IT,IU,IUT2,V= 2
            LIST(LADR11)=IT
            LIST(LADR11+1)=IU
            LIST(LADR11+2)=IUT2
            LIST(LADR11+3)=1
            LADR11=LADR11+4
          END DO
          IF(ISL3.EQ.1) THEN
           DO IUQ=1,NASHT
            IUU2=NOFF+KTU(IUQ,IUQ)-NTUES(ISL3)
C Add to list 11: IT,IT,IUU2,V= 1
            LIST(LADR11)=IT
            LIST(LADR11+1)=IT
            LIST(LADR11+2)=IUU2
            LIST(LADR11+3)=2
            LADR11=LADR11+4
           END DO
          END IF
         END DO
        END DO
       END DO

C Lists 14 and 15. I/IJ+-
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        LADR14=LLIST(ISL1,ISL3,14)
        LADR15=LLIST(ISL1,ISL3,15)
        ISL2=MUL(ISL1,ISL3)
        DO II=1,NISH(ISL1)
          IIQ=II+NIES(ISL1)
          DO IJ=1,NISH(ISL2)
            IJQ=IJ+NIES(ISL2)
            IF(IIQ.GT.IJQ) THEN
              IIJ=KIGEJ(IIQ,IJQ)-NIGEJES(ISL3)
C Add to list 14: II,IJ,IIJ,V= 1
              LIST(LADR14)=II
              LIST(LADR14+1)=IJ
              LIST(LADR14+2)=IIJ
              LIST(LADR14+3)=1
              LADR14=LADR14+4
              IIJ=KIGTJ(IIQ,IJQ)-NIGTJES(ISL3)
C Add to list 15: II,IJ,IIJ,V= 1
              LIST(LADR15)=II
              LIST(LADR15+1)=IJ
              LIST(LADR15+2)=IIJ
              LIST(LADR15+3)=1
              LADR15=LADR15+4
            ELSE IF(IIQ.LT.IJQ) THEN
              IIJ=KIGEJ(IJQ,IIQ)-NIGEJES(ISL3)
C Add to list 14: II,IJ,IIJ,V= 1
              LIST(LADR14)=II
              LIST(LADR14+1)=IJ
              LIST(LADR14+2)=IIJ
              LIST(LADR14+3)=1
              LADR14=LADR14+4
              IIJ=KIGTJ(IJQ,IIQ)-NIGTJES(ISL3)
C Add to list 15: II,IJ,IIJ,V=-1
              LIST(LADR15)=II
              LIST(LADR15+1)=IJ
              LIST(LADR15+2)=IIJ
              LIST(LADR15+3)=2
              LADR15=LADR15+4
            ELSE
              IIJ=KIGEJ(IIQ,IJQ)-NIGEJES(ISL3)
C Add to list 14: II,IJ,IIJ,V= Sqr(2)
              LIST(LADR14)=II
              LIST(LADR14+1)=IJ
              LIST(LADR14+2)=IIJ
              LIST(LADR14+3)=2
              LADR14=LADR14+4
            END IF
          END DO
         END DO
        END DO
       END DO

C Lists 16 and 17. A/AB+-
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        LADR16=LLIST(ISL1,ISL3,16)
        LADR17=LLIST(ISL1,ISL3,17)
        ISL2=MUL(ISL1,ISL3)
        DO IA=1,NSSH(ISL1)
          IAQ=IA+NSES(ISL1)
          DO IB=1,NSSH(ISL2)
            IBQ=IB+NSES(ISL2)
            IF(IAQ.GT.IBQ) THEN
              IAB=KAGEB(IAQ,IBQ)-NAGEBES(ISL3)
C Add to list 16: IA,IB,IAB,V= 1
              LIST(LADR16)=IA
              LIST(LADR16+1)=IB
              LIST(LADR16+2)=IAB
              LIST(LADR16+3)=1
              LADR16=LADR16+4
              IAB=KAGTB(IAQ,IBQ)-NAGTBES(ISL3)
C Add to list 17: IA,IB,IAB,V= 1
              LIST(LADR17)=IA
              LIST(LADR17+1)=IB
              LIST(LADR17+2)=IAB
              LIST(LADR17+3)=1
              LADR17=LADR17+4
            ELSE IF(IAQ.LT.IBQ) THEN
              IAB=KAGEB(IBQ,IAQ)-NAGEBES(ISL3)
C Add to list 16: IA,IB,IAB,V= 1
              LIST(LADR16)=IA
              LIST(LADR16+1)=IB
              LIST(LADR16+2)=IAB
              LIST(LADR16+3)=1
              LADR16=LADR16+4
              IAB=KAGTB(IBQ,IAQ)-NAGTBES(ISL3)
C Add to list 17: IA,IB,IAB,V=-1
              LIST(LADR17)=IA
              LIST(LADR17+1)=IB
              LIST(LADR17+2)=IAB
              LIST(LADR17+3)=2
              LADR17=LADR17+4
            ELSE
              IAB=KAGEB(IAQ,IBQ)-NAGEBES(ISL3)
C Add to list 16: IA,IB,IAB,V= Sqr(2)
              LIST(LADR16)=IA
              LIST(LADR16+1)=IB
              LIST(LADR16+2)=IAB
              LIST(LADR16+3)=2
              LADR16=LADR16+4
            END IF
          END DO
         END DO
        END DO
       END DO


      RETURN
      END
