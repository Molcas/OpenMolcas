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
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------
* 1998  PER-AAKE MALMQUIST
* DEPARTMENT OF THEORETICAL CHEMISTRY
* UNIVERSITY OF LUND
* SWEDEN
*--------------------------------------------
* NOTE: This new MKRHS code produces ONLY the
* contravariant components. 980928, P-A Malmqvist
*--------------------------------------------
      SUBROUTINE MKRHS(IVEC)
      use definitions, only: iwp, wp, u6
      use caspt2_global, only:iPrGlb
      use caspt2_global, only: FIMO
      use PrintLevel, only: verbose
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT,NOMX
      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) NERI, NFIMO
      real(kind=wp), ALLOCATABLE, TARGET:: ERI(:)
      real(kind=wp), POINTER:: ERI0(:), ERI1(:), ERI2(:), SCR(:)

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV. The coupling matrix elements from the
C root state to the 1st order interacting space are computed, as
C combinations of MO integrals.
C This is the RHS vector in contravariant representation.


      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(u6,'(1X,A)') ' Using conventional MKRHS algorithm'
      END IF

C INTEGRAL BUFFERS:
      NERI=NOMX**2
      CALL mma_allocate(ERI,3*NERI,Label='ERI')
      ERI0(1:2*NERI)=>ERI(1:2*NERI)
      ERI1(1:NERI)=>ERI(1:NERI)
      ERI2(1:NERI)=>ERI(NERI+1:2*NERI)
      SCR(1:NERI)=>ERI(2*NERI+1:3*NERI)

      IF(NASHT.GT.0) THEN
        NFIMO=SIZE(FIMO)
        CALL MKRHSA(IVEC,FIMO,NFIMO,ERI0,SCR)
        CALL MKRHSB(IVEC,ERI0,SCR)
        CALL MKRHSC(IVEC,FIMO,NFIMO,ERI0,SCR)
        CALL MKRHSD(IVEC,FIMO,NFIMO,ERI1,ERI2,SCR)
        CALL MKRHSE(IVEC,ERI1,ERI2,SCR)
        CALL MKRHSF(IVEC,ERI1,ERI2,SCR)
        CALL MKRHSG(IVEC,ERI1,ERI2,SCR)
      END IF
      CALL MKRHSH(IVEC,ERI1,ERI2,SCR)

      ERI0=>Null()
      ERI1=>Null()
      ERI2=>Null()
      SCR=>Null()
      CALL mma_deallocate(ERI)

      END SUBROUTINE MKRHS

      SUBROUTINE MKRHSA(IVEC,FIMO,NFIMO,ERI,SCR)
      use definitions, only: iwp, wp
      USE SUPERINDEX
      use EQSOLV
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,
     &                            Deallocate_GA_Array
      use caspt2_module

      IMPLICIT REAL*8 (A-H,O-Z)

      integer(kind=iwp), intent(in)::IVEC, NFIMO
      real(kind=wp), intent(inout):: FIMO(NFIMO), ERI(*), SCR(*)

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for case 1 (VJTU).

      NFNXT=0
      DO ISYM=1,NSYM
        NFIMOES=NFNXT
        NFNXT=NFNXT+(NORB(ISYM)*(NORB(ISYM)+1))/2
        IF(NINDEP(ISYM,1).EQ.0) CYCLE
          NAS=NTUV(ISYM)
          NIS=NISH(ISYM)
          NV=NAS*NIS
          IF(NV.EQ.0) CYCLE
C Set up a matrix FWI(w,i)=FIMO(wi)
          NI=NISH(ISYM)

C Compute W(tuv,i)=(ti,uv) + FIMO(t,i)*delta(u,v)/NACTEL
          LW=Allocate_GA_Array(NV,'WA')
          DO ISYMT=1,NSYM
            ISYMUV=MUL(ISYMT,ISYM)
            DO ISYMU=1,NSYM
              ISYMV=MUL(ISYMU,ISYMUV)
              DO IT=1,NASH(ISYMT)
                ITTOT=IT+NISH(ISYMT)
                ITABS=IT+NAES(ISYMT)
                DO II=1,NI
                  CALL COUL(ISYMU,ISYMV,ISYMT,ISYM,ITTOT,II,ERI,SCR)
                  ONEADD=0.0D0
                  IF(ISYMT.EQ.ISYM) THEN
                    FTI=FIMO(NFIMOES+(ITTOT*(ITTOT-1))/2+II)
                    ONEADD=FTI/DBLE(MAX(1,NACTEL))
                  END IF
                  DO IU=1,NASH(ISYMU)
                    IUTOT=IU+NISH(ISYMU)
                    IUABS=IU+NAES(ISYMU)
                    DO IV=1,NASH(ISYMV)
                      IVTOT=IV+NISH(ISYMV)
                      IVABS=IV+NAES(ISYMV)
                      IW1=KTUV(ITABS,IUABS,IVABS)-NTUVES(ISYM)
                      IW2=II
                      IW=IW1+NAS*(IW2-1)
                      IBUF=IUTOT+NORB(ISYMU)*(IVTOT-1)
                      WTUVI=ERI(IBUF)
                      IF(IVABS.EQ.IUABS) WTUVI=WTUVI+ONEADD
                      GA_Arrays(LW)%A(IW)=WTUVI
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C Put W on disk:
          ICASE=1
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
          Call Deallocate_GA_Array(LW)
      END DO

      END SUBROUTINE MKRHSA

      SUBROUTINE MKRHSB(IVEC,ERI,SCR)
      use definitions, only: iwp, wp
      USE SUPERINDEX
      use EQSOLV
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,
     &                            Deallocate_GA_Array
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)
      integer(kind=iwp), intent(in):: IVEC
      real(kind=wp), Intent(inout):: ERI(*), SCR(*)
*#define _KIGEJ_
*#define _KIGTJ_
*#include "mig_kig.fh"

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV for cases 2 and 3 (VJTI).


      SQ2=SQRT(2.0D00)
C VJTI CASE:
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
C   Allocate WP,WM
          LWP=Allocate_GA_Array(NVP,'WBP')
          LWM=Allocate_GA_Array(NVM,'WBM')
C   Let  W(tu,i,j)=(it,ju):
C   WP(tu,ij)=(W(tu,i,j)+W(tu,j,i))*(1-Kron(t,u)/2) /2
C With new normalisation, replace /2 with /(2*SQRT(1+Kron(ij))
C   WM(tu,ij)=(W(tu,i,j)-W(tu,j,i))*(1-Kron(t,u)/2) /2
          DO ISYMT=1,NSYM
            ISYMU=MUL(ISYMT,ISYM)
            IF(ISYMT.LT.ISYMU) CYCLE
            IF(NASH(ISYMT)*NASH(ISYMU).EQ.0) CYCLE
            DO ISYMI=1,NSYM
              ISYMJ=MUL(ISYMI,ISYM)
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
                  CALL EXCH(ISYMI,ISYMT,ISYMJ,ISYMU,
     &                      ITTOT,IUTOT,ERI,SCR)
                  IF(ITABS.NE.IUABS) THEN
                   DO II=1,NISH(ISYMI)
                    IIABS=II+NIES(ISYMI)
                    DO IJ=1,NISH(ISYMJ)
                      IJABS=IJ+NIES(ISYMJ)
                      IBUF=II+NORB(ISYMI)*(IJ-1)
                      VALUE=0.5D0*ERI(IBUF)
                      IF(IIABS.GE.IJABS) THEN
                        IIJP=KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
                        JWP=ITUP+NASP*(IIJP-1)
                        IF(IIABS.GT.IJABS) THEN
                          GA_Arrays(LWP)%A(JWP)=
     &                       GA_Arrays(LWP)%A(JWP)+VALUE
                          IIJM=KIGTJ(IIABS,IJABS)-NIGTJES(ISYM)
                          IWM=ITUM+NASM*(IIJM-1)
                          GA_Arrays(LWM)%A(IWM)=
     &                       GA_Arrays(LWM)%A(IWM)+VALUE
                        ELSE
                          GA_Arrays(LWP)%A(JWP)=
     &                       GA_Arrays(LWP)%A(JWP)+SQ2*VALUE
                        END IF
                      ELSE
                        IIJP=KIGEJ(IJABS,IIABS)-NIGEJES(ISYM)
                        JWP=ITUP+NASP*(IIJP-1)
                        GA_Arrays(LWP)%A(JWP)=
     &                     GA_Arrays(LWP)%A(JWP)+VALUE
                        IIJM=KIGTJ(IJABS,IIABS)-NIGTJES(ISYM)
                        IWM=ITUM+NASM*(IIJM-1)
                        GA_Arrays(LWM)%A(IWM)=
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
                      VALUE=0.25D0*ERI(IBUF)
                      IF(IIABS.GE.IJABS) THEN
                        IIJP=KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
                        JWP=ITUP+NASP*(IIJP-1)
                        IF(IIABS.GT.IJABS) THEN
                          GA_Arrays(LWP)%A(JWP)=
     &                       GA_Arrays(LWP)%A(JWP)+VALUE
                        ELSE
                          GA_Arrays(LWP)%A(JWP)=
     &                       GA_Arrays(LWP)%A(JWP)+SQ2*VALUE
                        END IF
                      ELSE
                        IIJP=KIGEJ(IJABS,IIABS)-NIGEJES(ISYM)
                        JWP=ITUP+NASP*(IIJP-1)
                        GA_Arrays(LWP)%A(JWP)=
     &                     GA_Arrays(LWP)%A(JWP)+VALUE
                      END IF
                    END DO
                   END DO
                  END IF
                END DO
              END DO
            END DO
          END DO
C   Put WP on disk
          ICASE=2
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
C  Put WM on disk
          IF(NINM.GT.0) THEN
            ICASE=3
            CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
          END IF
          Call Deallocate_GA_Array(LWM)
          Call Deallocate_GA_Array(LWP)
      END DO

      END SUBROUTINE MKRHSB

      SUBROUTINE MKRHSC(IVEC,FIMO,NFIMO,ERI,SCR)
      use definitions, only: iwp, wp
      USE SUPERINDEX, only: KTUV
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,
     &                            Deallocate_GA_Array
      use caspt2_module, only: NSYM,NORB,NINDEP,NTUV,NSSH,MUL,NASH,NISH,
     &                         NAES,NSSH,NTUVES,NASHT,NACTEL
      IMPLICIT REAL*8 (A-H,O-Z)
      integer(kind=iwp), intent(in):: IVEC, NFIMO
      real(kind=wp), intent(inout):: FIMO(NFIMO),ERI(*), SCR(*)

      integer(kind=iwp) NFNXT,ISYM,NFIMOES,NAS,NIS,NV,LW,ISYMT,
     &                  ISYMUV,ISYMU,ISYMV,IU,IUTOT,IUABS,IV,IVTOT,
     &                  IVABS,IA,IATOT,IT,ITTOT,ITABS,IW1,IW2,IW,IBUF,
     &                  IFIMO,IYABS,IYYW,IYYW
      real(kind=wp) SUM,ONEADD
C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV for case 4 (ATVX).

      NFNXT=0
      DO ISYM=1,NSYM
        NFIMOES=NFNXT
        NFNXT=NFNXT+(NORB(ISYM)*(NORB(ISYM)+1))/2
        IF(NINDEP(ISYM,4).EQ.0) CYCLE
          NAS=NTUV(ISYM)
          NIS=NSSH(ISYM)
          NV=NAS*NIS
          IF(NV.EQ.0) CYCLE

C   Allocate W. Put in W(tuv,a)=(at,uv) +
C             (FIMO(a,t)-sum(y)(ay,yt))*delta(u,v)/NACTEL.
C First, just the two-electron integrals. Later, add correction.

          LW=Allocate_GA_Array(NV,'WC')
          DO ISYMT=1,NSYM
            ISYMUV=MUL(ISYMT,ISYM)
            DO ISYMU=1,NSYM
              ISYMV=MUL(ISYMU,ISYMUV)
              DO IU=1,NASH(ISYMU)
                IUTOT=IU+NISH(ISYMU)
                IUABS=IU+NAES(ISYMU)
                DO IV=1,NASH(ISYMV)
                  IVTOT=IV+NISH(ISYMV)
                  IVABS=IV+NAES(ISYMV)
                  CALL COUL(ISYM,ISYMT,ISYMU,ISYMV,
     &                      IUTOT,IVTOT,ERI,SCR)
                  DO IA=1,NSSH(ISYM)
                    IATOT=IA+NISH(ISYM)+NASH(ISYM)
                    DO IT=1,NASH(ISYMT)
                      ITTOT=IT+NISH(ISYMT)
                      ITABS=IT+NAES(ISYMT)
                      IW1=KTUV(ITABS,IUABS,IVABS)-NTUVES(ISYM)
                      IW2=IA
                      IW=IW1+NAS*(IW2-1)
                      IBUF=IATOT+NORB(ISYM)*(ITTOT-1)
                      GA_Arrays(LW)%A(IW)=ERI(IBUF)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO

          DO IT=1,NASH(ISYM)
            ITTOT=IT+NISH(ISYM)
            ITABS=IT+NAES(ISYM)
            DO IA=1,NSSH(ISYM)
              IATOT=IA+NISH(ISYM)+NASH(ISYM)
              IFIMO=NFIMOES+(IATOT*(IATOT-1))/2+ITTOT
              SUM=FIMO(IFIMO)
              DO IYABS=1,NASHT
                IYYW=KTUV(IYABS,IYABS,ITABS)-NTUVES(ISYM)
                IYYWA=IYYW+NAS*(IA-1)
                SUM=SUM-GA_Arrays(LW)%A(IYYWA)
              END DO
              ONEADD=SUM/DBLE(MAX(1,NACTEL))
              DO ISYMU=1,NSYM
                DO IU=1,NASH(ISYMU)
                  IUABS=IU+NAES(ISYMU)
                  IW1=KTUV(ITABS,IUABS,IUABS)-NTUVES(ISYM)
                  IW2=IA
                  IW=IW1+NAS*(IW2-1)
                  GA_Arrays(LW)%A(IW)=GA_Arrays(LW)%A(IW)+ONEADD
                END DO
              END DO
            END DO
          END DO

C   Put W on disk
          ICASE=4
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LW)

          Call Deallocate_GA_Array(LW)
      END DO

      END SUBROUTINE MKRHSC

      SUBROUTINE MKRHSD(IVEC,FIMO,NFIMO,ERI1,ERI2,SCR)
      use definitions, only: iwp, wp
      use constants, only: Zero
      USE SUPERINDEX, only: KTU
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,
     &                            Deallocate_GA_Array
      use caspt2_module, only: NSYM,NINDEP,MUL,NSSH,NISH,NTU,NISUP,
     &                         NACTEL,NORB,NAES,NASH,NTUES
      IMPLICIT None
      integer(kind=iwp), intent(in):: IVEC, NFIMO
      real(kind=wp), intent(inout):: FIMO(NFIMO)
      real(kind=wp), intent(inout):: ERI1(*),ERI2(*), SCR(*)

      integer(kind=iwp) IOFF(8)
      integer(kind=iwp) ISYM,IO,ISYMI,NAS1,NAS,NIS,NV,LW,NFSUM,
     &                  NFIMOES,ISYMA,ISYMU,ISYMT,II,IU,IUABS,
     &                  IUTOT,IA,IATOT,IT,ITABS,ITTOT,IWA,IWI,IW1,IW2,
     &                  IBUF1,IBUF2,ICASE
      real(kind=wp) ONEADD,FAI,WAITU
C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for case 5, AIVX.

      DO ISYM=1,NSYM
        IF(NINDEP(ISYM,5).EQ.0) CYCLE
C Set up offset table:
          IO=0
          DO ISYMA=1,NSYM
            IOFF(ISYMA)=IO
            ISYMI=MUL(ISYMA,ISYM)
            IO=IO+NSSH(ISYMA)*NISH(ISYMI)
          END DO
C   Allocate W; W subdivided into W1,W2.
          NAS1=NTU(ISYM)
          NAS=2*NAS1
          NIS=NISUP(ISYM,5)
          NV=NAS*NIS
          IF(NV.EQ.0) CYCLE
C Compute W1(tu,ai)=(ai,tu) + FIMO(a,i)*delta(t,u)/NACTEL
C Compute W2(tu,ai)=(ti,au)
          LW=Allocate_GA_Array(NV,'WD')
          NFSUM=0
          DO ISYMI=1,NSYM
            NFIMOES=NFSUM
            NFSUM=NFSUM+(NORB(ISYMI)*(NORB(ISYMI)+1))/2
            ISYMA=MUL(ISYMI,ISYM)
            DO ISYMU=1,NSYM
              ISYMT=MUL(ISYMU,ISYM)
              DO II=1,NISH(ISYMI)
                DO IU=1,NASH(ISYMU)
                  IUABS=IU+NAES(ISYMU)
                  IUTOT=IU+NISH(ISYMU)
                  CALL EXCH(ISYMA,ISYMI,ISYMT,ISYMU,
     &                      II,IUTOT,ERI1,SCR)
                  CALL EXCH(ISYMT,ISYMI,ISYMA,ISYMU,
     &                      II,IUTOT,ERI2,SCR)
                  DO IA=1,NSSH(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    ONEADD=Zero
                    IF(ISYM.EQ.1) THEN
                      FAI=FIMO(NFIMOES+(IATOT*(IATOT-1))/2+II)
                      ONEADD=FAI/DBLE(MAX(1,NACTEL))
                    END IF
                    DO IT=1,NASH(ISYMT)
                      ITABS=IT+NAES(ISYMT)
                      ITTOT=IT+NISH(ISYMT)
                      IWA=KTU(ITABS,IUABS)-NTUES(ISYM)
                      IWI=II+NISH(ISYMI)*(IA-1)+IOFF(ISYMA)
                      IW1=IWA+NAS*(IWI-1)
                      IW2=IW1+NAS1
                      IBUF1=IATOT+NORB(ISYMA)*(ITTOT-1)
                      IBUF2=ITTOT+NORB(ISYMT)*(IATOT-1)
                      WAITU=ERI1(IBUF1)
                      IF(ITABS.EQ.IUABS) WAITU=WAITU+ONEADD
                      GA_Arrays(LW)%A(IW1)=WAITU
                      GA_Arrays(LW)%A(IW2)=ERI2(IBUF2)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C   Put W on disk.
          ICASE=5
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
          Call Deallocate_GA_Array(LW)
      END DO

      END SUBROUTINE MKRHSD

      SUBROUTINE MKRHSE(IVEC,ERI1,ERI2,SCR)
      use definitions, only: iwp, wp
      use constants, only: half, One, two, three
      USE SUPERINDEX, only: KIGEJ, KIGTJ
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,
     &                            Deallocate_GA_Array
      use caspt2_module, only: NSYM,NINDEP,NISUP,MUL,NASH,NISH,NSSH,
     &                         NORB,NIGEJ,NIES,NIGEJES,NIGTJES,NIGTJ
      IMPLICIT NONE
      integer(kind=iwp), intent(in):: IVEC
      real(kind=wp), Intent(inout):: ERI1(*),ERI2(*), SCR(*)

      integer(kind=iwp) IOFF1(8),IOFF2(8)
      real(kind=wp), parameter:: SQ2=SQRT(Two), SQI2=One/SQ2,
     &                           SQ3=SQRT(Three), SQ32=SQ3*SQI2
      integer(kind=iwp) ISYM,IO1,IO2,NAS,NISP,NISM,NVP,NVM,LWP,LWM,
     &                  ISYMA,ISYMI,IT,ITTOT,II,IA,IGEJ,IGTJ,IIABS,
     &                  IATOT,IBUF,IWA,IWIP,JWP,IJ,IJABS,ISYMIJ,ISYMJ,
     &                  IWIM,IWM,ICASE
      real(kind=wp) A, B
*#define _KIGEJ_
*#define _KIGTJ_
*#include "mig_kig.fh"

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for cases 6 and 7 (VJAI).


      DO ISYM=1,NSYM
        IF(NINDEP(ISYM,6)+NINDEP(ISYM,7).EQ.0) CYCLE
C Set up offset table:
          IO1=0
          IO2=0
          DO ISYMA=1,NSYM
            IOFF1(ISYMA)=IO1
            IOFF2(ISYMA)=IO2
            ISYMIJ=MUL(ISYMA,ISYM)
            IO1=IO1+NSSH(ISYMA)*NIGEJ(ISYMIJ)
            IO2=IO2+NSSH(ISYMA)*NIGTJ(ISYMIJ)
          END DO
C   Allocate W with parts WP,WM
          NAS=NASH(ISYM)
          NISP=NISUP(ISYM,6)
          NISM=NISUP(ISYM,7)
          NVP=NAS*NISP
          IF(NVP.EQ.0) CYCLE
          NVM=NAS*NISM
          LWP=Allocate_GA_Array(NVP,'WEP')
          LWM=Allocate_GA_Array(NVM,'WEM')
C  Let W(t,i,j,a)=(aitj)
C   WP(t,ij,a)=  (W(t,i,j,a)+W(t,j,i,a))
C With new normalisation, divide by /SQRT(2+2*Kron(ij))
C   WM(t,ij,a)=3*(W(t,i,j,a)-W(t,j,i,a))
C With new normalisation, divide by /SQRT(6)
          DO ISYMA=1,NSYM
            ISYMIJ=MUL(ISYMA,ISYM)
            DO ISYMI=1,NSYM
              ISYMJ=MUL(ISYMI,ISYMIJ)
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
C   Put WP and WM on disk.
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

      SUBROUTINE MKRHSF(IVEC,ERI1,ERI2,SCR)
      use definitions, only: iwp, wp
      use constants, only:  half, One, two
      USE SUPERINDEX, only: KTGEU,KAGEB,KTGTU,KAGTB
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,
     &                            Deallocate_GA_Array
      use caspt2_module, only: NSYM,NINDEP,NASUP,NISUP,MUL,NASH,NAES,
     &                         NISH,NSSH,NSES,NORB,NTGEUES,NAGEBES,
     &                         NTGTUES,NAGTBES
      IMPLICIT NONE
      integer(kind=iwp), intent(in):: IVEC
      real(kind=wp), Intent(inout):: ERI1(*),ERI2(*), SCR(*)

      real(kind=wp), parameter:: SQ2=SQRT(Two), SQI2=One/SQ2
      integer(kind=iwp) ISYM,NINP,NINM,NASP,NISP,NASM,NISM,NVP,NVM,LWP,
     &                  ISYMA,ISYMB,ISYMT,ISYMU,IT,ITABS,ITTOT,IU,IUABS,
     &                  IUTOT,IA,IAABS,IATOT,IB,IBABS,IBTOT,IBUF,IWAP,
     &                  IWIP,JWP,IWAM,IWIM,IWM,LWM,ICASE
      real(kind=wp) A, B
C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for cases 8 and 9 (BVAT).


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
C   Let W(t,u,ab)=(aubt)
C   WP(tu,ab)=(W(t,u,ab)+W(u,t,ab))*(1-Kron(t,u)/2) /2
C With new normalisation, replace /2 with /(2*SQRT(1+Kron(ab))
C   WM(tu,ab)=(W(t,u,ab)-W(u,t,ab))*(1-Kron(t,u)/2) /2
          DO ISYMA=1,NSYM
            ISYMB=MUL(ISYMA,ISYM)
            IF(ISYMA.LT.ISYMB) CYCLE
            DO ISYMT=1,NSYM
              ISYMU=MUL(ISYMT,ISYM)
              IF(ISYMT.LT.ISYMU) CYCLE
              DO IT=1,NASH(ISYMT)
                ITABS=IT+NAES(ISYMT)
                ITTOT=IT+NISH(ISYMT)
                DO IU=1,NASH(ISYMU)
                  IUABS=IU+NAES(ISYMU)
                  IUTOT=IU+NISH(ISYMU)
                  IF(ITABS.LT.IUABS) EXIT
                  CALL EXCH(ISYMA,ISYMU,ISYMB,ISYMT,
     &                      IUTOT,ITTOT,ERI1,SCR)
                  CALL EXCH(ISYMA,ISYMT,ISYMB,ISYMU,
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
C   Put WP on disk
          ICASE=8
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
          Call Deallocate_GA_Array(LWP)
          IF(NINM.GT.0) THEN
C   Put WM on disk
            ICASE=9
            CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
          END IF
          IF(NVM.GT.0) Call Deallocate_GA_Array(LWM)
      END DO

      END SUBROUTINE MKRHSF

      SUBROUTINE MKRHSG(IVEC,ERI1,ERI2,SCR)
      use definitions, only: iwp, wp
      use constants, only: half, One, two, three
      USE SUPERINDEX, only: KAGEB,KAGTB
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,
     &                            Deallocate_GA_Array
      use caspt2_module, only: NSYM,NINDEP,MUL,NISH,NAGEB,NAGTB,NASH,
     &                         NISUP,NISUP,NSSH,NSES,NORB,NAGEBES,
     &                         NAGTBES
      IMPLICIT NONE
      integer(kind=iwp), intent(in):: IVEC
      real(kind=wp), Intent(inout):: ERI1(*),ERI2(*), SCR(*)

      integer(kind=iwp) IOFF1(8),IOFF2(8)
      real(kind=wp), parameter:: SQ2=SQRT(Two), SQI2=One/SQ2,
     &                           SQ3=SQRT(Three), SQ32=SQ3*SQI2
      integer(kind=iwp) ISYM,IO1,IO2,NAS,NISP,NISM,NVP,NVM,LWP,LWM,
     &                  ISYMA,ISYMB,ISYMAB,ISYMI,IT,ITTOT,II,IA,IAABS,
     &                  IATOT,IB,IBABS,IBTOT,IBUF,IWA,IAGEB,IWIP,JWP,
     &                  IAGTB,IWIM,IWM,ICASE
      real(kind=wp) A,B

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for cases 10 and 11 (BJAT).


      DO ISYM=1,NSYM
        IF(NINDEP(ISYM,10)+NINDEP(ISYM,11).EQ.0) CYCLE
C Set up offset table:
          IO1=0
          IO2=0
          DO ISYMI=1,NSYM
            IOFF1(ISYMI)=IO1
            IOFF2(ISYMI)=IO2
            ISYMAB=MUL(ISYMI,ISYM)
            IO1=IO1+NISH(ISYMI)*NAGEB(ISYMAB)
            IO2=IO2+NISH(ISYMI)*NAGTB(ISYMAB)
          END DO
C   Allocate W with parts WP,WM
          NAS=NASH(ISYM)
          NISP=NISUP(ISYM,10)
          NISM=NISUP(ISYM,11)
          NVP=NAS*NISP
          IF(NVP.EQ.0) CYCLE
          NVM=NAS*NISM
          LWP=Allocate_GA_Array(NVP,'WGP')
          LWM=Allocate_GA_Array(NVM,'WGM')
C   Let  W(t,i,a,b)=(atbi)
C   WP(t,i,ab)=  (W(t,i,a,b)+W(t,i,b,a))
C With new normalisation, divide by /SQRT(2+2*Kron(ab))
C   WM(t,i,ab)=3*(W(t,i,a,b)-W(t,i,b,a))
C With new normalisation, divide by /SQRT(6)
          DO ISYMA=1,NSYM
            DO ISYMB=1,ISYMA
              ISYMAB=MUL(ISYMA,ISYMB)
              ISYMI=MUL(ISYMAB,ISYM)
              DO IT=1,NASH(ISYM)
                ITTOT=IT+NISH(ISYM)
                DO II=1,NISH(ISYMI)
                  CALL EXCH(ISYMA,ISYM ,ISYMB,ISYMI,
     &                      ITTOT,II,ERI1,SCR)
                  CALL EXCH(ISYMA,ISYMI,ISYMB,ISYM ,
     &                      II,ITTOT,ERI2,SCR)
                  DO IA=1,NSSH(ISYMA)
                    IAABS=IA+NSES(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    DO IB=1,NSSH(ISYMB)
                      IBABS=IB+NSES(ISYMB)
                      IBTOT=IB+NISH(ISYMB)+NASH(ISYMB)
                      IF(IAABS.LT.IBABS) EXIT
                      IBUF=IATOT+NORB(ISYMA)*(IBTOT-1)
                      IWA=IT
                      IAGEB=KAGEB(IAABS,IBABS)-NAGEBES(ISYMAB)
                      IWIP=II+NISH(ISYMI)*(IAGEB-1)+IOFF1(ISYMI)
                      JWP=IWA+NAS*(IWIP-1)
                      A=ERI1(IBUF)+ERI2(IBUF)
                      IF(IAABS.NE.IBABS) THEN
                        GA_Arrays(LWP)%A(JWP)=SQI2*A
                        IAGTB=KAGTB(IAABS,IBABS)-NAGTBES(ISYMAB)
                        IWIM=II+NISH(ISYMI)*(IAGTB-1)+IOFF2(ISYMI)
                        IWM=IWA+NAS*(IWIM-1)
                        B=ERI1(IBUF)-ERI2(IBUF)
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
C   Put WP and WM on disk.
          ICASE=10
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
          IF(NVM.GT.0) THEN
           ICASE=11
           CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
          END IF
          Call Deallocate_GA_Array(LWP)
          Call Deallocate_GA_Array(LWM)
      END DO

      END SUBROUTINE MKRHSG

      SUBROUTINE MKRHSH(IVEC,ERI1,ERI2,SCR)
      use definitions, only: iwp, wp
      use constants, only: half, One, two, three
      USE SUPERINDEX, only: KAGEB,KIGEJ,KAGTB,KIGTJ
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,
     &                            Deallocate_GA_Array
      use caspt2_module, only: NSYM,NAGEB,NIGEJ,NAGTB,NIGTJ,MUL,NISH,
     &                         NIES,NSES,NSSH,NORB,NASH,NAGEBES,NIGEJES,
     &                         NAGTBES,NIGTJES
      IMPLICIT None
      integer(kind=iwp), intent(in):: IVEC
      real(kind=wp), Intent(inout):: ERI1(*),ERI2(*), SCR(*)

      real(kind=wp), parameter:: SQ2=SQRT(Two), SQI2=One/SQ2,
     &                           SQ3=SQRT(Three)
      integer(kind=iwp) ISYM,NASP,NISP,NVP,NASM,NISM,NVM,LVP,ISYMA,
     &                  ISYMB,ISYMI,ISYMJ,II,IIABS,IJ,IJABS,
     &                  IA,IAABS,IATOT,IB,IBABS,IBTOT,IBUF,IVAP,IVIP,
     &                  IVP,IVAM,IVIM,IVM,LVM,ICASE
      real(kind=wp) A,B

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for cases 12 and 13 (BJAI).

      DO ISYM=1,NSYM
          NASP=NAGEB(ISYM)
          NISP=NIGEJ(ISYM)
          NVP=NASP*NISP
          IF(NVP.EQ.0) CYCLE
          NASM=NAGTB(ISYM)
          NISM=NIGTJ(ISYM)
          NVM=NASM*NISM
          LVP=Allocate_GA_Array(NVP,'WHP')
          IF(NVM.GT.0) LVM=Allocate_GA_Array(NVM,'WHM')
C   VP(ij,ab)=2*((aibj)+(ajbi))
C With new norm., divide by /SQRT(4*(1+Kron(ij))*(1+Kron(ab))
C   VM(ij,ab)=6*((aibj)-(ajbi))
C With new norm., divide by /SQRT(12)
          DO ISYMA=1,NSYM
            ISYMB=MUL(ISYMA,ISYM)
            IF(ISYMA.LT.ISYMB) CYCLE
            DO ISYMI=1,NSYM
              ISYMJ=MUL(ISYMI,ISYM)
              IF(ISYMI.LT.ISYMJ) CYCLE
              DO II=1,NISH(ISYMI)
                IIABS=II+NIES(ISYMI)
                DO IJ=1,NISH(ISYMJ)
                  IJABS=IJ+NIES(ISYMJ)
                  IF(IIABS.LT.IJABS) EXIT
                  CALL EXCH(ISYMA,ISYMI,ISYMB,ISYMJ,II,IJ,ERI1,SCR)
                  CALL EXCH(ISYMA,ISYMJ,ISYMB,ISYMI,IJ,II,ERI2,SCR)
                  DO IA=1,NSSH(ISYMA)
                    IAABS=IA+NSES(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    DO IB=1,NSSH(ISYMB)
                      IBABS=IB+NSES(ISYMB)
                      IF(IAABS.LT.IBABS) EXIT
                      IBTOT=IB+NISH(ISYMB)+NASH(ISYMB)
                      IBUF=IATOT+NORB(ISYMA)*(IBTOT-1)
                      IVAP=KAGEB(IAABS,IBABS)-NAGEBES(ISYM)
                      IVIP=KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
                      IVP=IVAP+NAGEB(ISYM)*(IVIP-1)
                      A=ERI1(IBUF)+ERI2(IBUF)
                      IF(IIABS.NE.IJABS) THEN
                        IF(IAABS.NE.IBABS) THEN
                          GA_Arrays(LVP)%A(IVP)=A
                          IVAM=KAGTB(IAABS,IBABS)-NAGTBES(ISYM)
                          IVIM=KIGTJ(IIABS,IJABS)-NIGTJES(ISYM)
                          IVM=IVAM+NAGTB(ISYM)*(IVIM-1)
                          B=ERI1(IBUF)-ERI2(IBUF)
                          GA_Arrays(LVM)%A(IVM)=SQ3*B
                        ELSE
                          GA_Arrays(LVP)%A(IVP)=SQI2*A
                        END IF
                      ELSE
                        IF(IAABS.NE.IBABS) THEN
                          GA_Arrays(LVP)%A(IVP)=SQI2*A
                        ELSE
                          GA_Arrays(LVP)%A(IVP)=Half*A
                        END IF
                      END IF
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO

          ICASE=12
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LVP)
          Call Deallocate_GA_Array(LVP)
          IF(NVM.GT.0) THEN
           ICASE=13
           CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LVM)
           Call Deallocate_GA_Array(LVM)
          END IF
      END DO

      END SUBROUTINE MKRHSH

      SUBROUTINE MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
CSVC: special routine to save the RHS array. MKRHS works in serial, so
C in case of a true parallel run we need to put the local array in a
C global array and then save that to disk in a distributed fashion.
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module
      IMPLICIT None

      integer(kind=iwp), intent(in):: ICASE,ISYM,IVEC,LW

      integer(kind=iwp) NAS, NIS, lg_w

      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)

#ifdef _MOLCAS_MPP_
      IF (IS_REAL_PAR()) THEN
        CALL RHS_ALLO(NAS,NIS,lg_W)
        CALL RHS_PUT(NAS,NIS,lg_W,GA_Arrays(LW)%A)
      ELSE
#endif
        lg_W=LW
#ifdef _MOLCAS_MPP_
      END IF
#endif

      CALL RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)

#ifdef _MOLCAS_MPP_
      IF (IS_REAL_PAR()) THEN
        CALL RHS_FREE(lg_W)
      END IF
#endif
      END SUBROUTINE MKRHS_SAVE

