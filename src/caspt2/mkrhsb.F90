!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------
! 1998  PER-AAKE MALMQUIST
! DEPARTMENT OF THEORETICAL CHEMISTRY
! UNIVERSITY OF LUND
! SWEDEN
!--------------------------------------------
      SUBROUTINE MKRHSB(IVEC,ERI,nERI,SCR,nSCR)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: Quart, Half, Two
      USE SUPERINDEX, only: KTGEU,KTGTU,KIGEJ,KIGTJ
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,                  &
     &                            Deallocate_GA_Array
      use caspt2_module, only: NSYM,NINDEP,NTGEU,NIGEJ,NTGTU,NIGTJ,     &
     &                         NASH,NISH,NAES,NTGEUES,NTGTUES,NIGEJES,  &
     &                         NIES,NORB,NIGTJES
      IMPLICIT None
      integer(kind=iwp), intent(in):: IVEC, nERI, nSCR
      real(kind=wp), Intent(inout):: ERI(nERI), SCR(nSCR)

      real(kind=wp), parameter:: SQ2=SQRT(Two)
      integer(kind=iwp) ISYM,NINP,NINM,NASP,NISP,NVP,NASM,NISM,NVM,LWP, &
     &                  LWM,ISYMT,ISYMU,ISYMI,ISYMJ,IT,ITABS,ITTOT,IU,  &
     &                  IUABS,IUTOT,ITUP,ITUM,II,IIABS,IJ,IJABS,IBUF,   &
     &                  IIJP,JWP,IIJM,IWM,ICASE
      real(kind=wp) VALUE

! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV for cases 2 and 3 (VJTI).

! VJTI CASE:
      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,2)
        NINM=NINDEP(ISYM,3)
        IF(NINP+NINM.EQ.0) CYCLE
          NASP=NTGEU(ISYM)
          NISP=NIGEJ(ISYM)
          NVP=NASP*NISP
          IF(NVP.EQ.0) CYCLE
          NASM=NTGTU(ISYM)
          NISM=NIGTJ(ISYM)
          NVM=NASM*NISM
!   Allocate WP,WM
          LWP=Allocate_GA_Array(NVP,'WBP')
          LWM=Allocate_GA_Array(NVM,'WBM')
!   Let  W(tu,i,j)=(it,ju):
!   WP(tu,ij)=(W(tu,i,j)+W(tu,j,i))*(1-Kron(t,u)/2) /2
! With new normalisation, replace /2 with /(2*SQRT(1+Kron(ij))
!   WM(tu,ij)=(W(tu,i,j)-W(tu,j,i))*(1-Kron(t,u)/2) /2
          DO ISYMT=1,NSYM
            ISYMU=Mul(ISYMT,ISYM)
            IF(ISYMT.LT.ISYMU) CYCLE
            IF(NASH(ISYMT)*NASH(ISYMU).EQ.0) CYCLE
            DO ISYMI=1,NSYM
              ISYMJ=Mul(ISYMI,ISYM)
              IF(NISH(ISYMI)*NISH(ISYMJ).EQ.0) CYCLE
              DO IT=1,NASH(ISYMT)
                ITABS=IT+NAES(ISYMT)
                ITTOT=IT+NISH(ISYMT)
                DO IU=1,NASH(ISYMU)
                  IUABS=IU+NAES(ISYMU)
                  IUTOT=IU+NISH(ISYMU)
                  IF(ITABS.LT.IUABS) EXIT
                  ITUP=KTGEU(ITABS,IUABS)-NTGEUES(ISYM)
                  ITUM=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
                  CALL EXCH(ISYMI,ISYMT,ISYMJ,ISYMU,                    &
     &                      ITTOT,IUTOT,ERI,SCR)
                  IF(ITABS.NE.IUABS) THEN
                   DO II=1,NISH(ISYMI)
                    IIABS=II+NIES(ISYMI)
                    DO IJ=1,NISH(ISYMJ)
                      IJABS=IJ+NIES(ISYMJ)
                      IBUF=II+NORB(ISYMI)*(IJ-1)
                      VALUE=Half*ERI(IBUF)
                      IF(IIABS.GE.IJABS) THEN
                        IIJP=KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
                        JWP=ITUP+NASP*(IIJP-1)
                        IF(IIABS.GT.IJABS) THEN
                          GA_Arrays(LWP)%A(JWP)=                        &
     &                       GA_Arrays(LWP)%A(JWP)+VALUE
                          IIJM=KIGTJ(IIABS,IJABS)-NIGTJES(ISYM)
                          IWM=ITUM+NASM*(IIJM-1)
                          GA_Arrays(LWM)%A(IWM)=                        &
     &                       GA_Arrays(LWM)%A(IWM)+VALUE
                        ELSE
                          GA_Arrays(LWP)%A(JWP)=                        &
     &                       GA_Arrays(LWP)%A(JWP)+SQ2*VALUE
                        END IF
                      ELSE
                        IIJP=KIGEJ(IJABS,IIABS)-NIGEJES(ISYM)
                        JWP=ITUP+NASP*(IIJP-1)
                        GA_Arrays(LWP)%A(JWP)=                          &
     &                     GA_Arrays(LWP)%A(JWP)+VALUE
                        IIJM=KIGTJ(IJABS,IIABS)-NIGTJES(ISYM)
                        IWM=ITUM+NASM*(IIJM-1)
                        GA_Arrays(LWM)%A(IWM)=                          &
     &                     GA_Arrays(LWM)%A(IWM)-VALUE
                      END IF
                    END DO
                   END DO
                  ELSE
                   DO II=1,NISH(ISYMI)
                    IIABS=II+NIES(ISYMI)
                    DO IJ=1,NISH(ISYMJ)
                      IJABS=IJ+NIES(ISYMJ)
                      IBUF=II+NORB(ISYMI)*(IJ-1)
                      VALUE=Quart*ERI(IBUF)
                      IF(IIABS.GE.IJABS) THEN
                        IIJP=KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
                        JWP=ITUP+NASP*(IIJP-1)
                        IF(IIABS.GT.IJABS) THEN
                          GA_Arrays(LWP)%A(JWP)=                        &
     &                       GA_Arrays(LWP)%A(JWP)+VALUE
                        ELSE
                          GA_Arrays(LWP)%A(JWP)=                        &
     &                       GA_Arrays(LWP)%A(JWP)+SQ2*VALUE
                        END IF
                      ELSE
                        IIJP=KIGEJ(IJABS,IIABS)-NIGEJES(ISYM)
                        JWP=ITUP+NASP*(IIJP-1)
                        GA_Arrays(LWP)%A(JWP)=                          &
     &                     GA_Arrays(LWP)%A(JWP)+VALUE
                      END IF
                    END DO
                   END DO
                  END IF
                END DO
              END DO
            END DO
          END DO
!   Put WP on disk
          ICASE=2
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
!  Put WM on disk
          IF(NINM.GT.0) THEN
            ICASE=3
            CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
          END IF
          Call Deallocate_GA_Array(LWM)
          Call Deallocate_GA_Array(LWP)
      END DO

      END SUBROUTINE MKRHSB
