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
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV. The coupling matrix elements from the
C root state to the 1st order interacting space are computed, as
C combinations of MO integrals.
C This is the RHS vector in contravariant representation.

      CALL QENTER('MKRHS')

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,'(1X,A)') ' Using conventional MKRHS algorithm'
      END IF

C INTEGRAL BUFFERS:
      NERI=NOMX**2
      CALL GETMEM('ERI','ALLO','REAL',LERI,3*NERI)
      LERI1=LERI
      LERI2=LERI+NERI
      LSCR =LERI+NERI*2

      IF(NASHT.GT.0) THEN
        CALL MKRHSA(IVEC,WORK(LFIMO),WORK(LERI),WORK(LSCR))
        CALL MKRHSB(IVEC,WORK(LERI),WORK(LSCR))
        CALL MKRHSC(IVEC,WORK(LFIMO),WORK(LERI),WORK(LSCR))
        CALL MKRHSD(IVEC,WORK(LFIMO),WORK(LERI1),WORK(LERI2),WORK(LSCR))
        CALL MKRHSE(IVEC,WORK(LERI1),WORK(LERI2),WORK(LSCR))
        CALL MKRHSF(IVEC,WORK(LERI1),WORK(LERI2),WORK(LSCR))
        CALL MKRHSG(IVEC,WORK(LERI1),WORK(LERI2),WORK(LSCR))
      END IF
      CALL MKRHSH(IVEC,WORK(LERI1),WORK(LERI2),WORK(LSCR))

      CALL GETMEM('ERI','FREE','REAL',LERI,2*NERI)

      CALL QEXIT('MKRHS')

      RETURN
      END

      SUBROUTINE MKRHSA(IVEC,FIMO,ERI,SCR)
      USE SUPERINDEX

      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
      DIMENSION FIMO(NFIMO), ERI(*), SCR(*)

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for case 1 (VJTU).

      CALL QENTER('MKRHSA')

      NFNXT=0
      DO 190 ISYM=1,NSYM
        NFIMOES=NFNXT
        NFNXT=NFNXT+(NORB(ISYM)*(NORB(ISYM)+1))/2
        IF(NINDEP(ISYM,1).EQ.0) GOTO 190
          NAS=NTUV(ISYM)
          NIS=NISH(ISYM)
          NV=NAS*NIS
          IF(NV.EQ.0) GOTO 190
C Set up a matrix FWI(w,i)=FIMO(wi)
          NI=NISH(ISYM)

C Compute W(tuv,i)=(ti,uv) + FIMO(t,i)*delta(u,v)/NACTEL
          CALL GETMEM('WA','ALLO','REAL',LW,NV)
          DO 130 ISYMT=1,NSYM
            ISYMUV=MUL(ISYMT,ISYM)
            DO 130 ISYMU=1,NSYM
              ISYMV=MUL(ISYMU,ISYMUV)
              DO 130 IT=1,NASH(ISYMT)
                ITTOT=IT+NISH(ISYMT)
                ITABS=IT+NAES(ISYMT)
                DO 130 II=1,NI
                  CALL COUL(ISYMU,ISYMV,ISYMT,ISYM,ITTOT,II,ERI,SCR)
                  ONEADD=0.0D0
                  IF(ISYMT.EQ.ISYM) THEN
                    FTI=FIMO(NFIMOES+(ITTOT*(ITTOT-1))/2+II)
                    ONEADD=FTI/DBLE(MAX(1,NACTEL))
                  END IF
                  DO 130 IU=1,NASH(ISYMU)
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
                      WORK(LW-1+IW)=WTUVI
                    END DO
 130      CONTINUE
C Put W on disk:
          ICASE=1
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
          CALL GETMEM('WA','FREE','REAL',LW,NV)
 190    CONTINUE

      CALL QEXIT('MKRHSA')

      RETURN
      END

      SUBROUTINE MKRHSB(IVEC,ERI,SCR)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
      DIMENSION ERI(*), SCR(*)
*#define _KIGEJ_
*#define _KIGTJ_
*#include <mig_kig.fh>

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV for cases 2 and 3 (VJTI).

      CALL QENTER('MKRHSB')

      SQ2=SQRT(2.0D00)
C VJTI CASE:
      DO 290 ISYM=1,NSYM
        NINP=NINDEP(ISYM,2)
        NINM=NINDEP(ISYM,3)
        IF(NINP+NINM.EQ.0) GOTO 290
          NASP=NTGEU(ISYM)
          NISP=NIGEJ(ISYM)
          NVP=NASP*NISP
          IF(NVP.EQ.0) GOTO 290
          NSBP=(NASP*(NASP+1))/2
          NASM=NTGTU(ISYM)
          NISM=NIGTJ(ISYM)
          NVM=NASM*NISM
          NSBM=(NASM*(NASM+1))/2
C   Allocate WP,WM
          NV=NVP+NVM
          CALL GETMEM('WB','ALLO','REAL',LW,NV)
          CALL DCOPY_(NV,[0.0D0],0,WORK(LW),1)
          LWP=LW
          LWM=LW+NVP
C   Let  W(tu,i,j)=(it,ju):
C   WP(tu,ij)=(W(tu,i,j)+W(tu,j,i))*(1-Kron(t,u)/2) /2
C With new normalisation, replace /2 with /(2*SQRT(1+Kron(ij))
C   WM(tu,ij)=(W(tu,i,j)-W(tu,j,i))*(1-Kron(t,u)/2) /2
          DO 240 ISYMT=1,NSYM
            ISYMU=MUL(ISYMT,ISYM)
            IF(ISYMT.LT.ISYMU) GOTO 240
            IF(NASH(ISYMT)*NASH(ISYMU).EQ.0) GOTO 240
            DO 230 ISYMI=1,NSYM
              ISYMJ=MUL(ISYMI,ISYM)
              IF(NISH(ISYMI)*NISH(ISYMJ).EQ.0) GOTO 230
              DO 220 IT=1,NASH(ISYMT)
                ITABS=IT+NAES(ISYMT)
                ITTOT=IT+NISH(ISYMT)
                DO 210 IU=1,NASH(ISYMU)
                  IUABS=IU+NAES(ISYMU)
                  IUTOT=IU+NISH(ISYMU)
                  IF(ITABS.LT.IUABS) GOTO 220
                  ITUP=KTGEU(ITABS,IUABS)-NTGEUES(ISYM)
                  ITUM=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
                  CALL EXCH(ISYMI,ISYMT,ISYMJ,ISYMU,
     &                      ITTOT,IUTOT,ERI,SCR)
                  IF(ITABS.NE.IUABS) THEN
                   DO 205 II=1,NISH(ISYMI)
                    IIABS=II+NIES(ISYMI)
                    DO 205 IJ=1,NISH(ISYMJ)
                      IJABS=IJ+NIES(ISYMJ)
                      IBUF=II+NORB(ISYMI)*(IJ-1)
                      VALUE=0.5D0*ERI(IBUF)
                      IF(IIABS.GE.IJABS) THEN
                        IIJP=KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
                        IWP=ITUP+NASP*(IIJP-1)
                        IF(IIABS.GT.IJABS) THEN
                          WORK(LWP-1+IWP)=WORK(LWP-1+IWP)+VALUE
                          IIJM=KIGTJ(IIABS,IJABS)-NIGTJES(ISYM)
                          IWM=ITUM+NASM*(IIJM-1)
                          WORK(LWM-1+IWM)=WORK(LWM-1+IWM)+VALUE
                        ELSE
                          WORK(LWP-1+IWP)=WORK(LWP-1+IWP)+SQ2*VALUE
                        END IF
                      ELSE
                        IIJP=KIGEJ(IJABS,IIABS)-NIGEJES(ISYM)
                        IWP=ITUP+NASP*(IIJP-1)
                        WORK(LWP-1+IWP)=WORK(LWP-1+IWP)+VALUE
                        IIJM=KIGTJ(IJABS,IIABS)-NIGTJES(ISYM)
                        IWM=ITUM+NASM*(IIJM-1)
                        WORK(LWM-1+IWM)=WORK(LWM-1+IWM)-VALUE
                      END IF
 205               CONTINUE
                  ELSE
                   DO 215 II=1,NISH(ISYMI)
                    IIABS=II+NIES(ISYMI)
                    DO 215 IJ=1,NISH(ISYMJ)
                      IJABS=IJ+NIES(ISYMJ)
                      IBUF=II+NORB(ISYMI)*(IJ-1)
                      VALUE=0.25D0*ERI(IBUF)
                      IF(IIABS.GE.IJABS) THEN
                        IIJP=KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
                        IWP=ITUP+NASP*(IIJP-1)
                        IF(IIABS.GT.IJABS) THEN
                          WORK(LWP-1+IWP)=WORK(LWP-1+IWP)+VALUE
                        ELSE
                          WORK(LWP-1+IWP)=WORK(LWP-1+IWP)+SQ2*VALUE
                        END IF
                      ELSE
                        IIJP=KIGEJ(IJABS,IIABS)-NIGEJES(ISYM)
                        IWP=ITUP+NASP*(IIJP-1)
                        WORK(LWP-1+IWP)=WORK(LWP-1+IWP)+VALUE
                      END IF
 215               CONTINUE
                  END IF
 210            CONTINUE
 220          CONTINUE
 230        CONTINUE
 240      CONTINUE
C   Put WP on disk
          ICASE=2
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
C  Put WM on disk
          IF(NINM.GT.0) THEN
            ICASE=3
            CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
          END IF

          CALL GETMEM('WB','FREE','REAL',LW,NV)
 290    CONTINUE

      CALL QEXIT('MKRHSB')

      RETURN
      END

      SUBROUTINE MKRHSC(IVEC,FIMO,ERI,SCR)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
      DIMENSION FIMO(NFIMO),ERI(*), SCR(*)

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV for case 4 (ATVX).

      CALL QENTER('MKRHSC')

      NFNXT=0
      DO 390 ISYM=1,NSYM
        NFIMOES=NFNXT
        NFNXT=NFNXT+(NORB(ISYM)*(NORB(ISYM)+1))/2
        IF(NINDEP(ISYM,4).EQ.0) GOTO 390
          NAS=NTUV(ISYM)
          NIS=NSSH(ISYM)
          NV=NAS*NIS
          IF(NV.EQ.0) GOTO 390

C   Allocate W. Put in W(tuv,a)=(at,uv) +
C             (FIMO(a,t)-sum(y)(ay,yt))*delta(u,v)/NACTEL.
C First, just the two-electron integrals. Later, add correction.

          CALL GETMEM('WC','ALLO','REAL',LW,NV)
          DO 310 ISYMT=1,NSYM
            ISYMUV=MUL(ISYMT,ISYM)
            DO 310 ISYMU=1,NSYM
              ISYMV=MUL(ISYMU,ISYMUV)
              DO 310 IU=1,NASH(ISYMU)
                IUTOT=IU+NISH(ISYMU)
                IUABS=IU+NAES(ISYMU)
                DO 310 IV=1,NASH(ISYMV)
                  IVTOT=IV+NISH(ISYMV)
                  IVABS=IV+NAES(ISYMV)
                  CALL COUL(ISYM,ISYMT,ISYMU,ISYMV,
     &                      IUTOT,IVTOT,ERI,SCR)
                  DO 310 IA=1,NSSH(ISYM)
                    IATOT=IA+NISH(ISYM)+NASH(ISYM)
                    DO 310 IT=1,NASH(ISYMT)
                      ITTOT=IT+NISH(ISYMT)
                      ITABS=IT+NAES(ISYMT)
                      IW1=KTUV(ITABS,IUABS,IVABS)-NTUVES(ISYM)
                      IW2=IA
                      IW=IW1+NAS*(IW2-1)
                      IBUF=IATOT+NORB(ISYM)*(ITTOT-1)
                      WORK(LW-1+IW)=ERI(IBUF)
 310      CONTINUE

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
                SUM=SUM-WORK(LW-1+IYYWA)
              END DO
              ONEADD=SUM/DBLE(MAX(1,NACTEL))
              DO ISYMU=1,NSYM
                DO IU=1,NASH(ISYMU)
                  IUABS=IU+NAES(ISYMU)
                  IW1=KTUV(ITABS,IUABS,IUABS)-NTUVES(ISYM)
                  IW2=IA
                  IW=IW1+NAS*(IW2-1)
                  WORK(LW-1+IW)=WORK(LW-1+IW)+ONEADD
                END DO
              END DO
            END DO
          END DO

C   Put W on disk
          ICASE=4
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LW)

          CALL GETMEM('WC','FREE','REAL',LW,NV)
 390    CONTINUE

      CALL QEXIT('MKRHSC')

      RETURN
      END

      SUBROUTINE MKRHSD(IVEC,FIMO,ERI1,ERI2,SCR)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
      DIMENSION IOFF(8),FIMO(NFIMO)
      DIMENSION ERI1(*),ERI2(*), SCR(*)

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for case 5, AIVX.

      CALL QENTER('MKRHSD')

      DO 490 ISYM=1,NSYM
        IF(NINDEP(ISYM,5).EQ.0) GOTO 490
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
          IF(NV.EQ.0) GOTO 490
          NSD=(NAS*(NAS+1))/2
C Compute W1(tu,ai)=(ai,tu) + FIMO(a,i)*delta(t,u)/NACTEL
C Compute W2(tu,ai)=(ti,au)
          CALL GETMEM('WD','ALLO','REAL',LW,NV)
          NFSUM=0
          DO 410 ISYMI=1,NSYM
            NFIMOES=NFSUM
            NFSUM=NFSUM+(NORB(ISYMI)*(NORB(ISYMI)+1))/2
            ISYMA=MUL(ISYMI,ISYM)
            DO 410 ISYMU=1,NSYM
              ISYMT=MUL(ISYMU,ISYM)
              DO 410 II=1,NISH(ISYMI)
                DO 410 IU=1,NASH(ISYMU)
                  IUABS=IU+NAES(ISYMU)
                  IUTOT=IU+NISH(ISYMU)
                  CALL EXCH(ISYMA,ISYMI,ISYMT,ISYMU,
     &                      II,IUTOT,ERI1,SCR)
                  CALL EXCH(ISYMT,ISYMI,ISYMA,ISYMU,
     &                      II,IUTOT,ERI2,SCR)
                  DO 410 IA=1,NSSH(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    ONEADD=0.0D0
                    IF(ISYM.EQ.1) THEN
                      FAI=FIMO(NFIMOES+(IATOT*(IATOT-1))/2+II)
                      ONEADD=FAI/DBLE(MAX(1,NACTEL))
                    END IF
                    DO 410 IT=1,NASH(ISYMT)
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
                      WORK(LW-1+IW1)=WAITU
                      WORK(LW-1+IW2)=ERI2(IBUF2)
 410      CONTINUE
C   Put W on disk.
          ICASE=5
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
          CALL GETMEM('WD','FREE','REAL',LW,NV)
 490    CONTINUE

      CALL QEXIT('MKRHSD')

      RETURN
      END

      SUBROUTINE MKRHSE(IVEC,ERI1,ERI2,SCR)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
      DIMENSION IOFF1(8),IOFF2(8)
      DIMENSION ERI1(*),ERI2(*), SCR(*)
*#define _KIGEJ_
*#define _KIGTJ_
*#include <mig_kig.fh>

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for cases 6 and 7 (VJAI).

      CALL QENTER('MKRHSE')

      SQ2=SQRT(2.0D00)
      SQI2=1.0D0/SQ2
      SQ3=SQRT(3.0D00)
      SQ32=SQ3*SQI2
      DO 590 ISYM=1,NSYM
        IF(NINDEP(ISYM,6)+NINDEP(ISYM,7).EQ.0) GOTO 590
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
          NIS=NISP+NISM
          NVP=NAS*NISP
          IF(NVP.EQ.0) GOTO 590
          NVM=NAS*NISM
          NSE=(NAS*(NAS+1))/2
          NV=NVP+NVM
          CALL GETMEM('WE','ALLO','REAL',LW,NV)
          LWP=LW
          LWM=LW+NVP
C  Let W(t,i,j,a)=(aitj)
C   WP(t,ij,a)=  (W(t,i,j,a)+W(t,j,i,a))
C With new normalisation, divide by /SQRT(2+2*Kron(ij))
C   WM(t,ij,a)=3*(W(t,i,j,a)-W(t,j,i,a))
C With new normalisation, divide by /SQRT(6)
          DO 540 ISYMA=1,NSYM
            ISYMIJ=MUL(ISYMA,ISYM)
            DO 530 ISYMI=1,NSYM
              ISYMJ=MUL(ISYMI,ISYMIJ)
              IF(ISYMI.LT.ISYMJ) GOTO 530
              DO 520 II=1,NISH(ISYMI)
                IIABS=II+NIES(ISYMI)
                DO 510 IJ=1,NISH(ISYMJ)
                  IJABS=IJ+NIES(ISYMJ)
                  IF(IIABS.LT.IJABS) GOTO 520
                  CALL EXCH(ISYMA,ISYMI,ISYM,ISYMJ,II,IJ,ERI1,SCR)
                  CALL EXCH(ISYMA,ISYMJ,ISYM,ISYMI,IJ,II,ERI2,SCR)
                  IGEJ=KIGEJ(IIABS,IJABS)-NIGEJES(ISYMIJ)
                  IGTJ=KIGTJ(IIABS,IJABS)-NIGTJES(ISYMIJ)
                  DO 510 IA=1,NSSH(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    DO 510 IT=1,NASH(ISYM)
                      ITTOT=IT+NISH(ISYM)
                      IBUF=IATOT+NORB(ISYMA)*(ITTOT-1)
                      A=ERI1(IBUF)+ERI2(IBUF)
                      IWA=IT
                      IWIP=IA+NSSH(ISYMA)*(IGEJ-1)+IOFF1(ISYMA)
                      IWP=IWA+NAS*(IWIP-1)
                      IF(IIABS.GT.IJABS) THEN
                        WORK(LWP-1+IWP)=SQI2*A
                        B=ERI1(IBUF)-ERI2(IBUF)
                        IWIM=IA+NSSH(ISYMA)*(IGTJ-1)+IOFF2(ISYMA)
                        IWM=IWA+NAS*(IWIM-1)
                        WORK(LWM-1+IWM)=SQ32*B
                      ELSE
                        WORK(LWP-1+IWP)=0.5D0*A
                      END IF
 510            CONTINUE
 520          CONTINUE
 530        CONTINUE
 540      CONTINUE
C   Put WP and WM on disk.
          ICASE=6
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
          IF(NVM.GT.0) THEN
            ICASE=7
            CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
          END IF
          CALL GETMEM('WE','FREE','REAL',LW,NV)
 590    CONTINUE

      CALL QEXIT('MKRHSE')

      RETURN
      END

      SUBROUTINE MKRHSF(IVEC,ERI1,ERI2,SCR)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
      DIMENSION ERI1(*),ERI2(*), SCR(*)

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for cases 8 and 9 (BVAT).

      CALL QENTER('MKRHSF')

      SQ2=SQRT(2.0D00)
      SQI2=1.0D0/SQ2

      DO 690 ISYM=1,NSYM
        NINP=NINDEP(ISYM,8)
        NINM=NINDEP(ISYM,9)
        IF(NINP+NINM.EQ.0) GOTO 690
          NASP=NASUP(ISYM,8)
          NISP=NISUP(ISYM,8)
          NASM=NASUP(ISYM,9)
          NISM=NISUP(ISYM,9)
          NVP=NASP*NISP
          IF(NVP.EQ.0)GOTO 690
          NVM=NASM*NISM
          NSFP=(NASP*(NASP+1))/2
          NSFM=(NASM*(NASM+1))/2
          CALL GETMEM('WFP','ALLO','REAL',LWP,NVP)
          IF(NVM.GT.0) CALL GETMEM('WFM','ALLO','REAL',LWM,NVM)
C   Let W(t,u,ab)=(aubt)
C   WP(tu,ab)=(W(t,u,ab)+W(u,t,ab))*(1-Kron(t,u)/2) /2
C With new normalisation, replace /2 with /(2*SQRT(1+Kron(ab))
C   WM(tu,ab)=(W(t,u,ab)-W(u,t,ab))*(1-Kron(t,u)/2) /2
          DO 640 ISYMA=1,NSYM
            ISYMB=MUL(ISYMA,ISYM)
            IF(ISYMA.LT.ISYMB) GOTO 640
            DO 630 ISYMT=1,NSYM
              ISYMU=MUL(ISYMT,ISYM)
              IF(ISYMT.LT.ISYMU) GOTO 630
              DO 620 IT=1,NASH(ISYMT)
                ITABS=IT+NAES(ISYMT)
                ITTOT=IT+NISH(ISYMT)
                DO 610 IU=1,NASH(ISYMU)
                  IUABS=IU+NAES(ISYMU)
                  IUTOT=IU+NISH(ISYMU)
                  IF(ITABS.LT.IUABS) GOTO 620
                  CALL EXCH(ISYMA,ISYMU,ISYMB,ISYMT,
     &                      IUTOT,ITTOT,ERI1,SCR)
                  CALL EXCH(ISYMA,ISYMT,ISYMB,ISYMU,
     &                      ITTOT,IUTOT,ERI2,SCR)
                  DO 610 IA=1,NSSH(ISYMA)
                    IAABS=IA+NSES(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    DO 600 IB=1,NSSH(ISYMB)
                      IBABS=IB+NSES(ISYMB)
                      IBTOT=IB+NISH(ISYMB)+NASH(ISYMB)
                      IF(IAABS.LT.IBABS) GOTO 610
                      IBUF=IATOT+NORB(ISYMA)*(IBTOT-1)
                      A=0.5D0*(ERI1(IBUF)+ERI2(IBUF))
                      IF(ITABS.EQ.IUABS) A=0.5D0*A
                      IWAP=KTGEU(ITABS,IUABS)-NTGEUES(ISYM)
                      IWIP=KAGEB(IAABS,IBABS)-NAGEBES(ISYM)
                      IWP=IWAP+NASP*(IWIP-1)
                      IF(IAABS.NE.IBABS) THEN
                        WORK(LWP-1+IWP)=A
                        IF(ITABS.NE.IUABS) THEN
                          B=0.5D0*(ERI1(IBUF)-ERI2(IBUF))
                          IWAM=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
                          IWIM=KAGTB(IAABS,IBABS)-NAGTBES(ISYM)
                          IWM=IWAM+NASM*(IWIM-1)
                          WORK(LWM-1+IWM)=B
                        END IF
                      ELSE
                        WORK(LWP-1+IWP)=SQI2*A
                      END IF
 600                CONTINUE
 610            CONTINUE
 620          CONTINUE
 630        CONTINUE
 640      CONTINUE
C   Put WP on disk
          ICASE=8
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
          CALL GETMEM('WFP','FREE','REAL',LWP,NVP)
          IF(NINM.GT.0) THEN
C   Put WM on disk
            ICASE=9
            CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
          END IF
          IF(NVM.GT.0) CALL GETMEM('WFM','FREE','REAL',LWM,NVM)
 690    CONTINUE

      CALL QEXIT('MKRHSF')

      RETURN
      END

      SUBROUTINE MKRHSG(IVEC,ERI1,ERI2,SCR)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
      DIMENSION IOFF1(8),IOFF2(8)
      DIMENSION ERI1(*),ERI2(*), SCR(*)

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for cases 10 and 11 (BJAT).

      CALL QENTER('MKRHSG')

      SQ2=SQRT(2.0D00)
      SQI2=1.0D0/SQ2
      SQ3=SQRT(3.0D00)
      SQ32=SQ3*SQI2
      DO 790 ISYM=1,NSYM
        IF(NINDEP(ISYM,10)+NINDEP(ISYM,11).EQ.0) GOTO 790
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
          NIS=NISP+NISM
          NVP=NAS*NISP
          IF(NVP.EQ.0) GOTO 790
          NVM=NAS*NISM
          NSG=(NAS*(NAS+1))/2
          NV=NVP+NVM
          CALL GETMEM('WG','ALLO','REAL',LW,NV)
          CALL DCOPY_(NV,[0.0D0],0,WORK(LW),1)
          LWP=LW
          LWM=LW+NVP
C   Let  W(t,i,a,b)=(atbi)
C   WP(t,i,ab)=  (W(t,i,a,b)+W(t,i,b,a))
C With new normalisation, divide by /SQRT(2+2*Kron(ab))
C   WM(t,i,ab)=3*(W(t,i,a,b)-W(t,i,b,a))
C With new normalisation, divide by /SQRT(6)
          DO 730 ISYMA=1,NSYM
            DO 730 ISYMB=1,ISYMA
              ISYMAB=MUL(ISYMA,ISYMB)
              ISYMI=MUL(ISYMAB,ISYM)
              DO 730 IT=1,NASH(ISYM)
                ITTOT=IT+NISH(ISYM)
                DO 730 II=1,NISH(ISYMI)
                  CALL EXCH(ISYMA,ISYM ,ISYMB,ISYMI,
     &                      ITTOT,II,ERI1,SCR)
                  CALL EXCH(ISYMA,ISYMI,ISYMB,ISYM ,
     &                      II,ITTOT,ERI2,SCR)
                  DO 720 IA=1,NSSH(ISYMA)
                    IAABS=IA+NSES(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    DO 710 IB=1,NSSH(ISYMB)
                      IBABS=IB+NSES(ISYMB)
                      IBTOT=IB+NISH(ISYMB)+NASH(ISYMB)
                      IF(IAABS.LT.IBABS) GOTO 720
                      IBUF=IATOT+NORB(ISYMA)*(IBTOT-1)
                      IWA=IT
                      IAGEB=KAGEB(IAABS,IBABS)-NAGEBES(ISYMAB)
                      IWIP=II+NISH(ISYMI)*(IAGEB-1)+IOFF1(ISYMI)
                      IWP=IWA+NAS*(IWIP-1)
                      A=ERI1(IBUF)+ERI2(IBUF)
                      IF(IAABS.NE.IBABS) THEN
                        WORK(LWP-1+IWP)=SQI2*A
                        IAGTB=KAGTB(IAABS,IBABS)-NAGTBES(ISYMAB)
                        IWIM=II+NISH(ISYMI)*(IAGTB-1)+IOFF2(ISYMI)
                        IWM=IWA+NAS*(IWIM-1)
                        B=ERI1(IBUF)-ERI2(IBUF)
                        WORK(LWM-1+IWM)=SQ32*B
                      ELSE
                        WORK(LWP-1+IWP)=0.5D0*A
                      END IF
 710                CONTINUE
 720              CONTINUE
 730      CONTINUE
C   Put WP and WM on disk.
          ICASE=10
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWP)
          IF(NVM.GT.0) THEN
           ICASE=11
           CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LWM)
          END IF
          CALL GETMEM('WG','FREE','REAL',LW,NV)
 790    CONTINUE

      CALL QEXIT('MKRHSG')

      RETURN
      END

      SUBROUTINE MKRHSH(IVEC,ERI1,ERI2,SCR)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
      DIMENSION ERI1(*),ERI2(*), SCR(*)
*#define _KIGEJ_
*#define _KIGTJ_
*#include <mig_kig.fh>

C Set up RHS vector of PT2 Linear Equation System, in vector
C number IVEC of LUSOLV, for cases 12 and 13 (BJAI).

      CALL QENTER('MKRHSH')

      SQ2=SQRT(2.0D00)
      SQI2=1.0D0/SQ2
      SQ3=SQRT(3.0D00)

      DO 890 ISYM=1,NSYM
          NASP=NAGEB(ISYM)
          NISP=NIGEJ(ISYM)
          NVP=NASP*NISP
          IF(NVP.EQ.0) GOTO 890
          NASM=NAGTB(ISYM)
          NISM=NIGTJ(ISYM)
          NVM=NASM*NISM
          CALL GETMEM('VP','ALLO','REAL',LVP,NVP)
          IF(NVM.GT.0) CALL GETMEM('VM','ALLO','REAL',LVM,NVM)
C   VP(ij,ab)=2*((aibj)+(ajbi))
C With new norm., divide by /SQRT(4*(1+Kron(ij))*(1+Kron(ab))
C   VM(ij,ab)=6*((aibj)-(ajbi))
C With new norm., divide by /SQRT(12)
          DO 840 ISYMA=1,NSYM
            ISYMB=MUL(ISYMA,ISYM)
            IF(ISYMA.LT.ISYMB) GOTO 840
            DO 830 ISYMI=1,NSYM
              ISYMJ=MUL(ISYMI,ISYM)
              IF(ISYMI.LT.ISYMJ) GOTO 830
              DO 820 II=1,NISH(ISYMI)
                IIABS=II+NIES(ISYMI)
                DO 810 IJ=1,NISH(ISYMJ)
                  IJABS=IJ+NIES(ISYMJ)
                  IF(IIABS.LT.IJABS) GOTO 820
                  CALL EXCH(ISYMA,ISYMI,ISYMB,ISYMJ,II,IJ,ERI1,SCR)
                  CALL EXCH(ISYMA,ISYMJ,ISYMB,ISYMI,IJ,II,ERI2,SCR)
                  DO 810 IA=1,NSSH(ISYMA)
                    IAABS=IA+NSES(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    DO 800 IB=1,NSSH(ISYMB)
                      IBABS=IB+NSES(ISYMB)
                      IF(IAABS.LT.IBABS) GOTO 810
                      IBTOT=IB+NISH(ISYMB)+NASH(ISYMB)
                      IBUF=IATOT+NORB(ISYMA)*(IBTOT-1)
                      IVAP=KAGEB(IAABS,IBABS)-NAGEBES(ISYM)
                      IVIP=KIGEJ(IIABS,IJABS)-NIGEJES(ISYM)
                      IVP=IVAP+NAGEB(ISYM)*(IVIP-1)
                      A=ERI1(IBUF)+ERI2(IBUF)
                      IF(IIABS.NE.IJABS) THEN
                        IF(IAABS.NE.IBABS) THEN
                          WORK(LVP-1+IVP)=A
                          IVAM=KAGTB(IAABS,IBABS)-NAGTBES(ISYM)
                          IVIM=KIGTJ(IIABS,IJABS)-NIGTJES(ISYM)
                          IVM=IVAM+NAGTB(ISYM)*(IVIM-1)
                          B=ERI1(IBUF)-ERI2(IBUF)
                          WORK(LVM-1+IVM)=SQ3*B
                        ELSE
                          WORK(LVP-1+IVP)=SQI2*A
                        END IF
                      ELSE
                        IF(IAABS.NE.IBABS) THEN
                          WORK(LVP-1+IVP)=SQI2*A
                        ELSE
                          WORK(LVP-1+IVP)=0.5D0*A
                        END IF
                      END IF
 800                CONTINUE
 810            CONTINUE
 820          CONTINUE
 830        CONTINUE
 840      CONTINUE

          ICASE=12
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LVP)
          CALL GETMEM('VP','FREE','REAL',LVP,NVP)
          IF(NVM.GT.0) THEN
           ICASE=13
           CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LVM)
           CALL GETMEM('VM','FREE','REAL',LVM,NVM)
          END IF
 890    CONTINUE

      CALL QEXIT('MKRHSH')

      RETURN
      END

      SUBROUTINE MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
CSVC: special routine to save the RHS array. MKRHS works in serial, so
C in case of a true parallel run we need to put the local array in a
C global array and then save that to disk in a distributed fashion.
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "WrkSpc.fh"
#include "caspt2.fh"
#include "para_info.fh"

      NAS=NASUP(ISYM,ICASE)
      NIS=NISUP(ISYM,ICASE)

#ifdef _MOLCAS_MPP_
      IF (IS_REAL_PAR()) THEN
        CALL RHS_ALLO(NAS,NIS,lg_W)
        CALL RHS_PUT(NAS,NIS,lg_W,WORK(LW))
      ELSE
        lg_W=LW
      END IF
#else
      lg_W=LW
#endif

      CALL RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)

#ifdef _MOLCAS_MPP_
      IF (IS_REAL_PAR()) THEN
        CALL RHS_FREE(NAS,NIS,lg_W)
      END IF
#endif
      END
