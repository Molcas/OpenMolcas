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
      SUBROUTINE MKRHSF(IVEC,ERI1,nERI1,ERI2,nERI2,SCR,nSCR)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only:  half, One, two
      USE SUPERINDEX, only: KTGEU,KAGEB,KTGTU,KAGTB
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,                  &
     &                            Deallocate_GA_Array
      use caspt2_module, only: NSYM,NINDEP,NASUP,NISUP,NASH,NAES,       &
     &                         NISH,NSSH,NSES,NORB,NTGEUES,NAGEBES,     &
     &                         NTGTUES,NAGTBES
      IMPLICIT NONE
      integer(kind=iwp), intent(in):: IVEC, nERI1, nERI2, nSCR
      real(kind=wp), Intent(inout):: ERI1(nERI1),ERI2(nERI2), SCR(nSCR)

      real(kind=wp), parameter:: SQ2=SQRT(Two), SQI2=One/SQ2
      integer(kind=iwp) ISYM,NINP,NINM,NASP,NISP,NASM,NISM,NVP,NVM,LWP, &
     &                  ISYMA,ISYMB,ISYMT,ISYMU,IT,ITABS,ITTOT,IU,IUABS,&
     &                  IUTOT,IA,IAABS,IATOT,IB,IBABS,IBTOT,IBUF,IWAP,  &
     &                  IWIP,JWP,IWAM,IWIM,IWM,LWM,ICASE
      real(kind=wp) A, B
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV, for cases 8 and 9 (BVAT).


      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,8)
        NINM=NINDEP(ISYM,9)
        IF(NINP+NINM.EQ.0) CYCLE
          NASP=NASUP(ISYM,8)
          NISP=NISUP(ISYM,8)
          NASM=NASUP(ISYM,9)
          NISM=NISUP(ISYM,9)
          NVP=NASP*NISP
          IF(NVP.EQ.0) CYCLE
          NVM=NASM*NISM
          LWP=Allocate_GA_Array(NVP,'WFP')
          IF(NVM.GT.0) LWM=Allocate_GA_Array(NVM,'WFM')
!   Let W(t,u,ab)=(aubt)
!   WP(tu,ab)=(W(t,u,ab)+W(u,t,ab))*(1-Kron(t,u)/2) /2
! With new normalisation, replace /2 with /(2*SQRT(1+Kron(ab))
!   WM(tu,ab)=(W(t,u,ab)-W(u,t,ab))*(1-Kron(t,u)/2) /2
          DO ISYMA=1,NSYM
            ISYMB=Mul(ISYMA,ISYM)
            IF(ISYMA.LT.ISYMB) CYCLE
            DO ISYMT=1,NSYM
              ISYMU=Mul(ISYMT,ISYM)
              IF(ISYMT.LT.ISYMU) CYCLE
              DO IT=1,NASH(ISYMT)
                ITABS=IT+NAES(ISYMT)
                ITTOT=IT+NISH(ISYMT)
                DO IU=1,NASH(ISYMU)
                  IUABS=IU+NAES(ISYMU)
                  IUTOT=IU+NISH(ISYMU)
                  IF(ITABS.LT.IUABS) EXIT
                  CALL EXCH(ISYMA,ISYMU,ISYMB,ISYMT,                    &
     &                      IUTOT,ITTOT,ERI1,SCR)
                  CALL EXCH(ISYMA,ISYMT,ISYMB,ISYMU,                    &
     &                      ITTOT,IUTOT,ERI2,SCR)
                  DO IA=1,NSSH(ISYMA)
                    IAABS=IA+NSES(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    DO IB=1,NSSH(ISYMB)
                      IBABS=IB+NSES(ISYMB)
                      IBTOT=IB+NISH(ISYMB)+NASH(ISYMB)
                      IF(IAABS.LT.IBABS) EXIT
                      IBUF=IATOT+NORB(ISYMA)*(IBTOT-1)
                      A=Half*(ERI1(IBUF)+ERI2(IBUF))
                      IF(ITABS.EQ.IUABS) A=Half*A
                      IWAP=KTGEU(ITABS,IUABS)-NTGEUES(ISYM)
                      IWIP=KAGEB(IAABS,IBABS)-NAGEBES(ISYM)
                      JWP=IWAP+NASP*(IWIP-1)
                      IF(IAABS.NE.IBABS) THEN
                        GA_Arrays(LWP)%A(JWP)=A
                        IF(ITABS.NE.IUABS) THEN
                          B=Half*(ERI1(IBUF)-ERI2(IBUF))
                          IWAM=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
                          IWIM=KAGTB(IAABS,IBABS)-NAGTBES(ISYM)
                          IWM=IWAM+NASM*(IWIM-1)
                          GA_Arrays(LWM)%A(IWM)=B
                        END IF
                      ELSE
                        GA_Arrays(LWP)%A(JWP)=SQI2*A
                      END IF
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
!   Put WP on disk
          ICASE=8
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
          Call Deallocate_GA_Array(LWP)
          IF(NINM.GT.0) THEN
!   Put WM on disk
            ICASE=9
            CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
          END IF
          IF(NVM.GT.0) Call Deallocate_GA_Array(LWM)
      END DO

      END SUBROUTINE MKRHSF
