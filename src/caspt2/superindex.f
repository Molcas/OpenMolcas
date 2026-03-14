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
      MODULE SUPERINDEX
      use definitions, only: iwp
      IMPLICIT NONE
      INTEGER(kind=iwp), ALLOCATABLE, SAVE ::
     &  KTU(:,:),    MTU(:,:),
     &  KTUV(:,:,:), MTUV(:,:),
     &  KTGEU(:,:),  MTGEU(:,:),
     &  KTGTU(:,:),  MTGTU(:,:),
     &  KAGEB(:,:),  MAGEB(:,:),
     &  KAGTB(:,:),  MAGTB(:,:),
     &  KIGEJ(:,:),  MIGEJ(:,:),
     &  KIGTJ(:,:),  MIGTJ(:,:),
     &  KIA(:,:),    MIA(:,:),
     &  MIREL(:,:), MTREL(:,:), MAREL(:,:)

      CONTAINS
      SUBROUTINE SUPINI
      use stdalloc, only: mma_allocate
      use Symmetry_Info, only: Mul
      use caspt2_module, only: nInDep, MxCase, nAshT, nCases, nIshT,
     &                         nAshT, nSym, nTUVEs, nAsh, nSshT,
     &                         nAes, nTUV, nTUES,
     &                         nTGEUES, nTGTUES, nTU, nTGEU,
     &                         nTGTU, nIGEJES, nIGTJES, nIsh,
     &                         nIES, nIGEJ, nIGTJ, nAGEBES,
     &                         nAGTBES, nSsh, nSES, nAGEB, nIAES,
     &                         nSES, nSsh, nAGTB, Cases, nASUP, nISUP
      IMPLICIT NONE
      CHARACTER(LEN=8), Parameter :: CSNAME(MXCASE)=[
     &              'VJTU    ','VJTIP   ','VJTIM   ',
     &     'ATVX    ','AIVX    ','VJAIP   ','VJAIM   ','BVATP   ',
     &     'BVATM   ','BJATP   ','BJATM   ','BJAIP   ','BJAIM   ']
      INTEGER(kind=iwp) ISYM,ICASE
      INTEGER(kind=iwp) NSUM,NCOUNT
      INTEGER(kind=iwp) IA,IB,II,IJ,IT,IU,IV
      INTEGER(kind=iwp) IAQ,IBQ,IIQ,IJQ,ITQ,IUQ,IVQ
      INTEGER(kind=iwp) ISA,ISB,ISI,ISJ,IST,ISU,ISV
      INTEGER(kind=iwp) IAGEB,IAGTB,NMAGEB,NMAGTB
      INTEGER(kind=iwp) IIGEJ,IIGTJ,NMIGEJ,NMIGTJ
      INTEGER(kind=iwp) ITGEU,ITGTU,NMTGEU,NMTGTU,ITUV
      INTEGER(kind=iwp) IIA,ITU,NMIA
      INTEGER(kind=iwp) IS1,IS2,ISYA,ISYI,ISUV
      INTEGER(kind=iwp) JC0,JCM,JCP
      INTEGER(kind=iwp) NC0,NCM,NCP
      INTEGER(kind=iwp) JCOUNT
      INTEGER(kind=iwp) N,N5,N6,N7,N10,N11
      INTEGER(kind=iwp) NAT,NAU,NAV,NII,NIJ,NSA,NSB


      CALL MMA_ALLOCATE(KTUV,NASHT,NASHT,NASHT,Label='KTUV')
      CALL MMA_ALLOCATE(MTUV,3,NASHT**3,Label='MTUV')
      ITUV=0
      DO ISYM=1,NSYM
        NCOUNT=0
        NTUVES(ISYM)=ITUV
        DO ISV=1,NSYM
          NAV=NASH(ISV)
          DO ISU=1,NSYM
            NAU=NASH(ISU)
            ISUV=Mul(ISU,ISV)
            IST=Mul(ISUV,ISYM)
            NAT=NASH(IST)
            JCOUNT=NAV*NAU*NAT
            IF (JCOUNT.EQ.0) CYCLE
            NCOUNT=NCOUNT+JCOUNT
            DO IV=1,NAV
              IVQ=NAES(ISV)+IV
              DO IU=1,NAU
                IUQ=NAES(ISU)+IU
                DO IT=1,NAT
                  ITQ=NAES(IST)+IT
                  ITUV=ITUV+1
                  KTUV(ITQ,IUQ,IVQ)=ITUV
                  MTUV(1,ITUV)=ITQ
                  MTUV(2,ITUV)=IUQ
                  MTUV(3,ITUV)=IVQ
                END DO
              END DO
            END DO
          END DO
        END DO
        NTUV(ISYM)=NCOUNT
      END DO

      CALL MMA_ALLOCATE(KTU,NASHT,NASHT,Label='KTU')
      CALL MMA_ALLOCATE(MTU,3,NASHT**2,Label='MTU')

      CALL MMA_ALLOCATE(KTGEU,NASHT,NASHT,Label='KTGEU')
      CALL MMA_ALLOCATE(KTGTU,NASHT,NASHT,Label='KTGTU')
      NMTGEU=(NASHT*(NASHT+1))/2
      NMTGTU=(NASHT*(NASHT-1))/2
      CALL MMA_ALLOCATE(MTGEU,2,NMTGEU,Label='MTGEU')
      CALL MMA_ALLOCATE(MTGTU,2,NMTGTU,Label='MTGTU')

      ITU=0
      ITGEU=0
      ITGTU=0
      DO ISYM=1,NSYM
        NC0=0
        NCP=0
        NCM=0
        NTUES(ISYM)=ITU
        NTGEUES(ISYM)=ITGEU
        NTGTUES(ISYM)=ITGTU
        DO ISU=1,NSYM
          NAU=NASH(ISU)
          IST=Mul(ISU,ISYM)
          NAT=NASH(IST)
          JC0=0
          JCP=0
          JCM=0
          DO IU=1,NAU
            IUQ=IU+NAES(ISU)
            DO IT=1,NAT
              ITQ=IT+NAES(IST)
              JC0=JC0+1
              ITU=ITU+1
              KTU(ITQ,IUQ)=ITU
              MTU(1,ITU)=ITQ
              MTU(2,ITU)=IUQ
              IF(ITQ.LT.IUQ) CYCLE
              JCP=JCP+1
              ITGEU=ITGEU+1
              KTGEU(ITQ,IUQ)=ITGEU
              MTGEU(1,ITGEU)=ITQ
              MTGEU(2,ITGEU)=IUQ
              IF(ITQ.LE.IUQ) CYCLE
              JCM=JCM+1
              ITGTU=ITGTU+1
              KTGTU(ITQ,IUQ)=ITGTU
              MTGTU(1,ITGTU)=ITQ
              MTGTU(2,ITGTU)=IUQ
            END DO
          END DO
          NC0=NC0+JC0
          NCP=NCP+JCP
          NCM=NCM+JCM
        END DO
        NTU(ISYM)=NC0
        NTGEU(ISYM)=NCP
        NTGTU(ISYM)=NCM
      END DO

CPAM99 Use allocated workspace instead of MAGEB, MAGTB:
      NMAGEB=(NSSHT*(NSSHT+1))/2
      NMAGTB=(NSSHT*(NSSHT-1))/2
      CALL MMA_ALLOCATE(MAGEB,2,NMAGEB,Label='MAGEB')
      CALL MMA_ALLOCATE(MAGTB,2,NMAGTB,Label='MAGTB')
      NMIGEJ=(NISHT*(NISHT+1))/2
      NMIGTJ=(NISHT*(NISHT-1))/2
      CALL MMA_ALLOCATE(MIGEJ,2,NMIGEJ,Label='MIGEJ')
      CALL MMA_ALLOCATE(MIGTJ,2,NMIGTJ,Label='MIGTJ')

      CALL MMA_ALLOCATE(KIGEJ,NISHT,NISHT,Label='KIGEJ')
      CALL MMA_ALLOCATE(KIGTJ,NISHT,NISHT,Label='KIGTJ')
      CALL MMA_ALLOCATE(KAGEB,NSSHT,NSSHT,Label='KAGEB')
      CALL MMA_ALLOCATE(KAGTB,NSSHT,NSSHT,Label='KAGTB')
      KIGEJ(:,:) = 0
      KIGTJ(:,:) = 0

      CALL MMA_ALLOCATE(KIA,NISHT,NSSHT,Label='KIA')
      NMIA=NISHT*NSSHT
      CALL MMA_ALLOCATE(MIA,2,NMIA,Label='MIA')

C Construct tables for inactive and secondary pair indices:
      IIGEJ=0
      IIGTJ=0
      IAGEB=0
      IAGTB=0
      IIA=0
      DO ISYM=1,NSYM
C Inactive pair indices:
        NIGEJES(ISYM)=IIGEJ
        NIGTJES(ISYM)=IIGTJ
        NCM=0
        NCP=0
        DO ISI=1,NSYM
          ISJ=Mul(ISI,ISYM)
          IF(ISI.LT.ISJ) CYCLE
          NII=NISH(ISI)
          NIJ=NISH(ISJ)
          DO II=1,NII
            IIQ=NIES(ISI)+II
            DO IJ=1,NIJ
              IJQ=NIES(ISJ)+IJ
              NCP=NCP+1
              IIGEJ=IIGEJ+1
              KIGEJ(IIQ,IJQ)=IIGEJ
              MIGEJ(1,IIGEJ)=IIQ
              MIGEJ(2,IIGEJ)=IJQ
              IF(IIQ.LE.IJQ) EXIT
              NCM=NCM+1
              IIGTJ=IIGTJ+1
              KIGTJ(IIQ,IJQ)=IIGTJ
              MIGTJ(1,IIGTJ)=IIQ
              MIGTJ(2,IIGTJ)=IJQ
            END DO
          END DO
        END DO

        NIGEJ(ISYM)=NCP
        NIGTJ(ISYM)=NCM

C Secondary pair indices:
        NAGEBES(ISYM)=IAGEB
        NAGTBES(ISYM)=IAGTB
        NCM=0
        NCP=0
        DO ISA=1,NSYM
          ISB=Mul(ISA,ISYM)
          IF(ISA.LT.ISB) CYCLE
          NSA=NSSH(ISA)
          NSB=NSSH(ISB)
          DO IA=1,NSA
            IAQ=NSES(ISA)+IA
            DO IB=1,NSB
              IBQ=NSES(ISB)+IB
              NCP=NCP+1
              IAGEB=IAGEB+1
              KAGEB(IAQ,IBQ)=IAGEB
              MAGEB(1,IAGEB)=IAQ
              MAGEB(2,IAGEB)=IBQ
              IF(IAQ.LE.IBQ) EXIT
              NCM=NCM+1
              IAGTB=IAGTB+1
              KAGTB(IAQ,IBQ)=IAGTB
              MAGTB(1,IAGTB)=IAQ
              MAGTB(2,IAGTB)=IBQ
            END DO
          END DO
        END DO
        NAGEB(ISYM)=NCP
        NAGTB(ISYM)=NCM

C Inactive-Secondary pair indices:
        NIAES(ISYM)=IIA
        DO ISYA=1,NSYM
          ISYI=Mul(ISYA,ISYM)
          DO IA=1,NSSH(ISYA)
            IAQ=IA+NSES(ISYA)
            DO II=1,NISH(ISYI)
              IIQ=II+NIES(ISYI)
              IIA=IIA+1
              KIA(IIQ,IAQ)=IIA
              MIA(1,IIA)=IIQ
              MIA(2,IIA)=IAQ
            END DO
          END DO
        END DO

      END DO
* End of loop over symmetries.


      DO ICASE=1,NCASES
        CASES(ICASE)=CSNAME(ICASE)
      END DO

      DO ISYM=1,NSYM
        NASUP(ISYM,1 )= NTUV(ISYM)
        NASUP(ISYM,2 )= NTGEU(ISYM)
        NASUP(ISYM,3 )= NTGTU(ISYM)
        NASUP(ISYM,4 )= NTUV(ISYM)
        NASUP(ISYM,5 )= 2*NTU(ISYM)
        NASUP(ISYM,6 )= NASH(ISYM)
        NASUP(ISYM,7 )= NASH(ISYM)
        NASUP(ISYM,8 )= NTGEU(ISYM)
        NASUP(ISYM,9 )= NTGTU(ISYM)
        NASUP(ISYM,10)= NASH(ISYM)
        NASUP(ISYM,11)= NASH(ISYM)
        NASUP(ISYM,12)= NAGEB(ISYM)
        NASUP(ISYM,13)= NAGTB(ISYM)
        N5 =0
        N6 =0
        N7 =0
        N10=0
        N11=0
        DO IS1=1,NSYM
          IS2=Mul(IS1,ISYM)
          N5 =N5 +NSSH(IS1)*NISH(IS2)
          N6 =N6 +NSSH(IS1)*NIGEJ(IS2)
          N7 =N7 +NSSH(IS1)*NIGTJ(IS2)
          N10=N10+NAGEB(IS1)*NISH(IS2)
          N11=N11+NAGTB(IS1)*NISH(IS2)
        END DO
        NISUP(ISYM,1 )= NISH(ISYM)
        NISUP(ISYM,2 )= NIGEJ(ISYM)
        NISUP(ISYM,3 )= NIGTJ(ISYM)
        NISUP(ISYM,4 )= NSSH(ISYM)
        NISUP(ISYM,5 )= N5
        NISUP(ISYM,6 )= N6
        NISUP(ISYM,7 )= N7
        NISUP(ISYM,8 )= NAGEB(ISYM)
        NISUP(ISYM,9 )= NAGTB(ISYM)
        NISUP(ISYM,10)= N10
        NISUP(ISYM,11)= N11
        NISUP(ISYM,12)= NIGEJ(ISYM)
        NISUP(ISYM,13)= NIGTJ(ISYM)
      END DO
      DO ICASE=1,NCASES
        NSUM=0
        DO ISYM=1,NSYM
          N=NASUP(ISYM,ICASE)*NISUP(ISYM,ICASE)
C Preliminary value for NINDEP: Nr of independent active params:
          NINDEP(ISYM,ICASE)=NASUP(ISYM,ICASE)
          IF(N.EQ.0) NINDEP(ISYM,ICASE)=0
          NSUM=NSUM+N
        END DO
      END DO

CSVC: prepare tables to translate from absolute indices to
C(index,symmetry) pairs.

      CALL MMA_ALLOCATE(MIREL,2,NISHT,Label='MIREL')
      CALL MMA_ALLOCATE(MTREL,2,NASHT,Label='MTREL')
      CALL MMA_ALLOCATE(MAREL,2,NSSHT,Label='MAREL')

      DO ISYM=1,NSYM
        DO II=1,NISH(ISYM)
          IIQ=II+NIES(ISYM)
          MIREL(1,IIQ)=II
          MIREL(2,IIQ)=ISYM
        END DO
        DO IT=1,NASH(ISYM)
          ITQ=IT+NAES(ISYM)
          MTREL(1,ITQ)=IT
          MTREL(2,ITQ)=ISYM
        END DO
        DO IA=1,NSSH(ISYM)
          IAQ=IA+NSES(ISYM)
          MAREL(1,IAQ)=IA
          MAREL(2,IAQ)=ISYM
        END DO
      END DO


      END SUBROUTINE SUPINI

      SUBROUTINE SUPFREE()
      use stdalloc, only: mma_deallocate
      IMPLICIT NONE
      ! deallocate the superindex tables
      CALL MMA_DEALLOCATE(KIGEJ)
      CALL MMA_DEALLOCATE(KIGTJ)
      CALL MMA_DEALLOCATE(MIGEJ)
      CALL MMA_DEALLOCATE(MIGTJ)
      CALL MMA_DEALLOCATE(MAGEB)
      CALL MMA_DEALLOCATE(MAGTB)
      CALL MMA_DEALLOCATE(KAGEB)
      CALL MMA_DEALLOCATE(KAGTB)
      CALL MMA_DEALLOCATE(MTGEU)
      CALL MMA_DEALLOCATE(MTGTU)
      CALL MMA_DEALLOCATE(KTGEU)
      CALL MMA_DEALLOCATE(KTGTU)
      CALL MMA_DEALLOCATE(KTU)
      CALL MMA_DEALLOCATE(MTU)
      CALL MMA_DEALLOCATE(KTUV)
      CALL MMA_DEALLOCATE(MTUV)
      CALL MMA_DEALLOCATE(KIA)
      CALL MMA_DEALLOCATE(MIA)
      CALL MMA_DEALLOCATE(MIREL)
      CALL MMA_DEALLOCATE(MTREL)
      CALL MMA_DEALLOCATE(MAREL)
      END SUBROUTINE SUPFREE

      END MODULE SUPERINDEX
