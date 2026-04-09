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
      SUBROUTINE NSIND(INS,ISYM,ICASE,IP,IQ,IR)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, u6
      USE SUPERINDEX, only: MIGEJ,MIGTJ,MAGTB,MAGEB
      use caspt2_module, only: NIES,IINAIS,NIGEJES,NSES,NSYM,NIGTJES,
     &                         IEXTIS,NISH,NSSH,NIGEJ,NIGTJ,NAGEBES,
     &                         NAGTBES,NAGEB,NAGTB
      IMPLICIT None
      integer(kind=iwp), intent(in):: INS,ISYM,ICASE
      integer(kind=iwp), intent(out):: IP,IQ,IR

      integer(kind=iwp) IA, IAABS, IAB, IABABS, IAIJ, IBABS, II, IIA,
     &                  IIAB, IIABS, IIJ, IIJABS, IJABS, ISYMA, ISYMAB,
     &                  ISYMI, ISYMIJ, NA, NAB, NAIJ, NI, NIA, NIAB, NIJ

      Select case(ICASE)
      Case(1)
      IIABS=INS+NIES(ISYM)
      IP=IINAIS(IIABS)
      IQ=0
      IR=0
      Case(2,12)
      IIJABS=INS+NIGEJES(ISYM)
      IIABS=MIGEJ(1,IIJABS)
      IJABS=MIGEJ(2,IIJABS)
      IP=IINAIS(IIABS)
      IQ=IINAIS(IJABS)
      IR=0
      Case(3,13)
      IIJABS=INS+NIGTJES(ISYM)
      IIABS=MIGTJ(1,IIJABS)
      IJABS=MIGTJ(2,IIJABS)
      IP=IINAIS(IIABS)
      IQ=IINAIS(IJABS)
      IR=0
      Case(4)
      IAABS=INS+NSES(ISYM)
      IP=IEXTIS(IAABS)
      IQ=0
      IR=0
      RETURN
      Case(5)
      IIA=INS
      DO ISYMA=1,NSYM
        ISYMI=Mul(ISYMA,ISYM)
        NIA=NISH(ISYMI)*NSSH(ISYMA)
        IF(IIA.LE.NIA) THEN
          IA=1+(IIA-1)/NISH(ISYMI)
          II=IIA-NISH(ISYMI)*(IA-1)
          IAABS=IA+NSES(ISYMA)
          IIABS=II+NIES(ISYMI)
          IP=IINAIS(IIABS)
          IQ=IEXTIS(IAABS)
          IR=0
          RETURN
        END IF
        IIA=IIA-NIA
      END DO
      WRITE(u6,*)'NSIND AIVX: Impossible situation.'
      CALL ABEND()
      Case(6,7)
      IAIJ=INS
      NIJ = 0 ! dummy initialize
      DO ISYMA=1,NSYM
        ISYMIJ=Mul(ISYMA,ISYM)
        IF(ICASE.EQ.6) NIJ=NIGEJ(ISYMIJ)
        IF(ICASE.EQ.7) NIJ=NIGTJ(ISYMIJ)
        NA=NSSH(ISYMA)
        NAIJ=NA*NIJ
        IF(IAIJ.LE.NAIJ) THEN
          IIJ=1+(IAIJ-1)/NA
          IA=IAIJ-NA*(IIJ-1)
          IAABS=IA+NSES(ISYMA)
          IF(ICASE.EQ.6) THEN
            IIJABS=IIJ+NIGEJES(ISYMIJ)
            IIABS=MIGEJ(1,IIJABS)
            IJABS=MIGEJ(2,IIJABS)
          ELSE
            IIJABS=IIJ+NIGTJES(ISYMIJ)
            IIABS=MIGTJ(1,IIJABS)
            IJABS=MIGTJ(2,IIJABS)
          END IF
          IP=IEXTIS(IAABS)
          IQ=IINAIS(IIABS)
          IR=IINAIS(IJABS)
          RETURN
        END IF
        IAIJ=IAIJ-NAIJ
      END DO
      WRITE(u6,*)'NSIND VJAI: Impossible situation.'
      CALL ABEND()
      CASE(8)
      IABABS=INS+NAGEBES(ISYM)
      IAABS=MAGEB(1,IABABS)
      IBABS=MAGEB(2,IABABS)
      IP=IEXTIS(IAABS)
      IQ=IEXTIS(IBABS)
      IR=0
      CASE(9)
      IABABS=INS+NAGTBES(ISYM)
      IAABS=MAGTB(1,IABABS)
      IBABS=MAGTB(2,IABABS)
      IP=IEXTIS(IAABS)
      IQ=IEXTIS(IBABS)
      IR=0
      CASE(10,11)
      IIAB=INS
      NAB = 0 ! dummy initialize
      DO ISYMI=1,NSYM
        ISYMAB=Mul(ISYMI,ISYM)
        IF(ICASE.EQ.10) NAB=NAGEB(ISYMAB)
        IF(ICASE.EQ.11) NAB=NAGTB(ISYMAB)
        NI=NISH(ISYMI)
        NIAB=NI*NAB
        IF(IIAB.LE.NIAB) THEN
          IAB=1+(IIAB-1)/NI
          II=IIAB-NI*(IAB-1)
          IIABS=II+NIES(ISYMI)
          IF(ICASE.EQ.10) THEN
            IABABS=IAB+NAGEBES(ISYMAB)
            IAABS=MAGEB(1,IABABS)
            IBABS=MAGEB(2,IABABS)
          ELSE
            IABABS=IAB+NAGTBES(ISYMAB)
            IAABS=MAGTB(1,IABABS)
            IBABS=MAGTB(2,IABABS)
          END IF
          IP=IINAIS(IIABS)
          IQ=IEXTIS(IAABS)
          IR=IEXTIS(IBABS)
          RETURN
        END IF
        IIAB=IIAB-NIAB
      END DO
      CASE DEFAULT
      WRITE(u6,*)'NSIND BJAT: Impossible situation.'
      CALL ABEND()
      End Select

      END SUBROUTINE NSIND
