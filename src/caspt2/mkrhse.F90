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
      SUBROUTINE MKRHSE(IVEC,ERI1,nERI1,ERI2,nERI2,SCR,nSCR)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: half, One, two, three
      USE SUPERINDEX, only: KIGEJ, KIGTJ
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,                  &
     &                            Deallocate_GA_Array
      use caspt2_module, only: NSYM,NINDEP,NISUP,NASH,NISH,NSSH,        &
     &                         NORB,NIGEJ,NIES,NIGEJES,NIGTJES,NIGTJ
      IMPLICIT NONE
      integer(kind=iwp), intent(in):: IVEC, nERI1, nERI2, nSCR
      real(kind=wp), Intent(inout):: ERI1(nERI1),ERI2(nERI2), SCR(nSCR)

      integer(kind=iwp) IOFF1(8),IOFF2(8)
      real(kind=wp), parameter:: SQ2=SQRT(Two), SQI2=One/SQ2,           &
     &                           SQ3=SQRT(Three), SQ32=SQ3*SQI2
      integer(kind=iwp) ISYM,IO1,IO2,NAS,NISP,NISM,NVP,NVM,LWP,LWM,     &
     &                  ISYMA,ISYMI,IT,ITTOT,II,IA,IGEJ,IGTJ,IIABS,     &
     &                  IATOT,IBUF,IWA,IWIP,JWP,IJ,IJABS,ISYMIJ,ISYMJ,  &
     &                  IWIM,IWM,ICASE
      real(kind=wp) A, B
!#define _KIGEJ_
!#define _KIGTJ_
!#include "mig_kig.fh"

! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV, for cases 6 and 7 (VJAI).


      DO ISYM=1,NSYM
        IF(NINDEP(ISYM,6)+NINDEP(ISYM,7).EQ.0) CYCLE
! Set up offset table:
          IO1=0
          IO2=0
          DO ISYMA=1,NSYM
            IOFF1(ISYMA)=IO1
            IOFF2(ISYMA)=IO2
            ISYMIJ=Mul(ISYMA,ISYM)
            IO1=IO1+NSSH(ISYMA)*NIGEJ(ISYMIJ)
            IO2=IO2+NSSH(ISYMA)*NIGTJ(ISYMIJ)
          END DO
!   Allocate W with parts WP,WM
          NAS=NASH(ISYM)
          NISP=NISUP(ISYM,6)
          NISM=NISUP(ISYM,7)
          NVP=NAS*NISP
          IF(NVP.EQ.0) CYCLE
          NVM=NAS*NISM
          LWP=Allocate_GA_Array(NVP,'WEP')
          LWM=Allocate_GA_Array(NVM,'WEM')
!  Let W(t,i,j,a)=(aitj)
!   WP(t,ij,a)=  (W(t,i,j,a)+W(t,j,i,a))
! With new normalisation, divide by /SQRT(2+2*Kron(ij))
!   WM(t,ij,a)=3*(W(t,i,j,a)-W(t,j,i,a))
! With new normalisation, divide by /SQRT(6)
          DO ISYMA=1,NSYM
            ISYMIJ=Mul(ISYMA,ISYM)
            DO ISYMI=1,NSYM
              ISYMJ=Mul(ISYMI,ISYMIJ)
              IF(ISYMI.LT.ISYMJ) CYCLE
              DO II=1,NISH(ISYMI)
                IIABS=II+NIES(ISYMI)
                DO IJ=1,NISH(ISYMJ)
                  IJABS=IJ+NIES(ISYMJ)
                  IF(IIABS.LT.IJABS) EXIT
                  CALL EXCH(ISYMA,ISYMI,ISYM,ISYMJ,II,IJ,ERI1,SCR)
                  CALL EXCH(ISYMA,ISYMJ,ISYM,ISYMI,IJ,II,ERI2,SCR)
                  IGEJ=KIGEJ(IIABS,IJABS)-NIGEJES(ISYMIJ)
                  IGTJ=KIGTJ(IIABS,IJABS)-NIGTJES(ISYMIJ)
                  DO IA=1,NSSH(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    DO IT=1,NASH(ISYM)
                      ITTOT=IT+NISH(ISYM)
                      IBUF=IATOT+NORB(ISYMA)*(ITTOT-1)
                      A=ERI1(IBUF)+ERI2(IBUF)
                      IWA=IT
                      IWIP=IA+NSSH(ISYMA)*(IGEJ-1)+IOFF1(ISYMA)
                      JWP=IWA+NAS*(IWIP-1)
                      IF(IIABS.GT.IJABS) THEN
                        GA_Arrays(LWP)%A(JWP)=SQI2*A
                        B=ERI1(IBUF)-ERI2(IBUF)
                        IWIM=IA+NSSH(ISYMA)*(IGTJ-1)+IOFF2(ISYMA)
                        IWM=IWA+NAS*(IWIM-1)
                        GA_Arrays(LWM)%A(IWM)=SQ32*B
                      ELSE
                        GA_Arrays(LWP)%A(JWP)=half*A
                      END IF
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
!   Put WP and WM on disk.
          ICASE=6
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
          IF(NVM.GT.0) THEN
            ICASE=7
            CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
          END IF
          Call Deallocate_GA_Array(LWP)
          Call Deallocate_GA_Array(LWM)
      END DO

      END SUBROUTINE MKRHSE
