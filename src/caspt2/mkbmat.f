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
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MKBMAT()
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug, verbose
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_global, only: DREF, PREF
      use caspt2_global, only: LUSOLV, LUSBT
      use EQSOLV
      IMPLICIT REAL*8 (A-H,O-Z)
C Set up B matrices for cases 1..13.

#include "caspt2.fh"
#include "pt2_guga.fh"
      REAL*8 DUM(1)
      INTEGER*1, ALLOCATABLE :: idxG3(:,:)
      REAL*8, ALLOCATABLE:: F1(:), F2(:), F3(:), FD(:), FP(:)


      IF(IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,*)' Construct B matrices'
      END IF

      IF(NASHT.EQ.0) GOTO 100

      CALL mma_allocate(F1,NG1,Label='F1')
      NFD=SIZE(PREF)
      CALL mma_allocate(FD,NFD,Label='FD')
      CALL PT2_GET(NG1,'DELTA1',F1)
      CALL MKDREF_RPT2(NASHT,F1,FD)
      CALL mma_deallocate(F1)
      CALL mma_allocate(F2,NG2,Label='F2')
      CALL PT2_GET(NG2,'DELTA2',F2)
      NFP=SIZE(PREF)
      CALL mma_allocate(FP,NFP,Label='FP')
      CALL MKPREF_RPT2(NASHT,F2,FP)
      CALL mma_deallocate(F2)
      CALL mma_allocate(F3,NG3,Label='F3')
      CALL PT2_GET(NG3,'DELTA3',F3)

      IF(IPRGLB.GE.DEBUG) THEN
        WRITE(6,'("DEBUG> ",A)') 'CASE SYM B-MATRIX NORM'
        WRITE(6,'("DEBUG> ",A)') '==== === ============='
      END IF

      CALL mma_allocate(idxG3,6,NG3,label='idxG3')
      iLUID=0
      CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

      CALL MKBA(DREF,SIZE(DREF),PREF,SIZE(PREF),FD,FP,NG3,F3,idxG3)
      CALL MKBC(DREF,SIZE(DREF),PREF,SIZE(PREF),FD,FP,NG3,F3,idxG3)

      CALL mma_deallocate(F3)
      CALL mma_deallocate(idxG3)

      CALL MKBB(DREF,SIZE(DREF),PREF,SIZE(PREF),FD,FP)
      CALL MKBD(DREF,SIZE(DREF),PREF,SIZE(PREF),FD,FP)
      CALL MKBE(DREF,SIZE(DREF),FD)
      CALL MKBF(DREF,SIZE(DREF),PREF,SIZE(PREF),FP)
      CALL MKBG(DREF,SIZE(DREF),FD)
      CALL mma_deallocate(FP)
      CALL mma_deallocate(FD)

 100  CONTINUE

C For completeness, even case H has formally S and B
C matrices. This costs nothing, and saves conditional
C looping, etc in the rest  of the routines.
      DUM(1)=0.0D0
      DO ISYM=1,NSYM
        DO ICASE=12,13
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.GT.0) THEN
            IDISK=IDBMAT(ISYM,ICASE)
            CALL DDAFILE(LUSBT,1,DUM,1,IDISK)
          END IF
        END DO
      END DO

      END SUBROUTINE MKBMAT

********************************************************************************
* Case A (ICASE=1)
********************************************************************************
      SUBROUTINE MKBA(DREF,NDREF,PREF,NPREF,FD,FP,NG3,F3,idxG3)
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug
      USE SUPERINDEX
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      use EQSOLV
      IMPLICIT REAL*8 (A-H,O-Z)
#include "caspt2.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      INTEGER NDREF,NPREF, NG3
      Real*8 DREF(NDREF),PREF(NPREF),F3(NG3)
      Real*8 FD(NDREF),FP(NPREF)
      INTEGER*1 idxG3(6,NG3)
#ifdef _MOLCAS_MPP_
      Real*8 Dummy(1)
#endif

      ICASE=1
C LONG loop over superindex symmetry.
      DO ISYM=1,NSYM
        NIN=NINDEP(ISYM,ICASE)
        IF(NIN.EQ.0) CYCLE
        NAS=NTUV(ISYM)
        NBA=(NAS*(NAS+1))/2
        IF(NBA.LE.0) CYCLE

C Set up the matrix BA(tuv,xyz) defined by the expression
C <ituv|H0-E0|kxyz> = dik ( alpha(a) SA(tuv,xyz) + BA(tuv,xyz) )
C Formula used:
C BA(tuv,xyz) = (Ey+Eu+Ex+Et-EASUM)*SA(tuv,xyz) - Fvuxtyz
C - dyu ( Fvzxt - Eu Gvzxt ) - dyt ( Fvuxz - Et Gvuxz )
C - dxu ( Fvtyz - Eu Gvtyz ) - dxu dyt ( Fvz - (Et+Eu) Gvz )
C + 2dxt ( Fvuyz - Et Gvuyz ) + 2dxt dyu ( Fvz - (Et+Eu) Gvz )

C where dyu = Kronecker(y,u) etc. Gvutxyz=<Evutxyz>, etc.
C Similarly, Fvutxyz= Sum(w)(EPSA(w)<Evutxyzww>, etc.

        CALL PSBMAT_GETMEM('BA',lg_BA,NAS)
        CALL PSBMAT_READ('S',iCase,iSym,lg_BA,NAS)

        ! fill in the 3-el parts
#ifdef _MOLCAS_MPP_
        IF (IS_REAL_PAR()) THEN
          MYRANK = GA_NODEID()
          CALL GA_DISTRIBUTION (LG_BA,MYRANK,ILO,IHI,JLO,JHI)
          IF (JLO.NE.0 .AND. (JHI-JLO+1).NE.NAS) THEN
            WRITE(6,*) 'MKBA: MISMATCH IN RANGE OF THE SUPERINDICES'
            CALL ABEND()
          END IF
          IF (ILO.GT.0 .AND. JLO.GT.0) THEN
            CALL GA_ACCESS (LG_BA,ILO,IHI,JLO,JHI,MA,LDA)
            CALL MKBA_DP(DREF,NDREF,PREF,NPREF,FD,FP,iSYM,
     &                   DBL_MB(MA),ILO,IHI,JLO,JHI,LDA)
            CALL MKBA_F3_MPP(ISYM,DBL_MB(MA),ILO,IHI,JLO,JHI,LDA,
     &                       NG3,F3,IDXG3)
            CALL GA_RELEASE_UPDATE (LG_BA,ILO,IHI,JLO,JHI)
          ELSE
            CALL MKBA_F3_MPP(ISYM,DUMMY,ILO,IHI,JLO,JHI,LDA,
     &                       NG3,F3,IDXG3)
          END IF
        ELSE
#endif
          CALL MKBA_DP(DREF,NDREF,PREF,NPREF,FD,FP,
     &                 ISYM,GA_Arrays(lg_BA)%A(:),1,NAS,1,NAS,0)
          CALL MKBA_F3(ISYM,GA_Arrays(lg_BA)%A(:),NG3,F3,IDXG3)
#ifdef _MOLCAS_MPP_
        END IF
#endif

        CALL PSBMAT_WRITE('B',iCase,iSYM,lg_BA,NAS)

        IF(IPRGLB.GE.DEBUG) THEN
          DBA=PSBMAT_FPRINT(lg_BA,NAS)
          WRITE(6,'("DEBUG> ",A4,1X,I3,1X,ES21.14)') 'A', ISYM, DBA
        END IF

        CALL PSBMAT_FREEMEM(lg_BA)
      END DO

      END SUBROUTINE MKBA

      SUBROUTINE MKBA_DP (DREF,NDREF,PREF,NPREF,FD,FP,iSYM,
     &                    BA,iLo,iHi,jLo,jHi,LDA)
      USE SUPERINDEX
      use caspt2_global, only:ipea_shift
      use EQSOLV
      IMPLICIT REAL*8 (A-H,O-Z)
#include "caspt2.fh"
      INTEGER NDREF, NPREF, iSYM, iLo, iHi, jLo, jHi, LDA
      DIMENSION DREF(NDREF),PREF(NPREF)
      DIMENSION FD(NDREF),FP(NPREF)
      DIMENSION BA(*)

CSV.20100831: fill in the F2 and F1 corrections for this BA block
C on entry, BA should contain SA!!
      DO IXYZ=jLo,jHi
        IXYZABS=IXYZ+NTUVES(ISYM)
        IXABS=MTUV(1,IXYZABS)
        IYABS=MTUV(2,IXYZABS)
        IZABS=MTUV(3,IXYZABS)
        EX=EPSA(IXABS)
        EY=EPSA(IYABS)
        DO ITUV=iLo,iHi
          ITUVABS=ITUV+NTUVES(ISYM)
          ITABS=MTUV(1,ITUVABS)
          IUABS=MTUV(2,ITUVABS)
          IVABS=MTUV(3,ITUVABS)
          ET=EPSA(ITABS)
          EU=EPSA(IUABS)
          ETU=ET+EU
          FACT=EY+EU+EX+ET-EASUM
          ISADR=1+iTUV-iLo+LDA*(iXYZ-jLo)
          IF (LDA.EQ.0) THEN
            IF (iXYZ.LE.iTUV) THEN
              ISADR=(ITUV*(ITUV-1))/2+IXYZ
            ELSE
              CYCLE
            END IF
          END IF
          VALUE=FACT*BA(ISADR)
          IF(IYABS.EQ.IUABS) THEN
            IVZ=IVABS+NASHT*(IZABS-1)
            IXT=IXABS+NASHT*(ITABS-1)
            IP1=MAX(IVZ,IXT)
            IP2=MIN(IVZ,IXT)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE-2.0D0*(FP(IP)-EU*PREF(IP))
            IF(IXABS.EQ.ITABS) THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE+2.0D0*(FD(ID)-ETU*DREF(ID))
            END IF
          END IF
C Add  dyt ( -Fvuxz + Et*Gvuxz +dxu (-Fvz+(Et+Eu)*Gvz))
          IF(IYABS.EQ.ITABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            IXZ=IXABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,IXZ)
            IP2=MIN(IVU,IXZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE-2.0D0*(FP(IP)-ET*PREF(IP))
            IF(IXABS.EQ.IUABS) THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE - (FD(ID)-ETU*DREF(ID))
            END IF
          END IF
C Add  dxu ( -Fvtyz + Eu*Gvtyz )
          IF(IXABS.EQ.IUABS) THEN
            IVT=IVABS+NASHT*(ITABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVT,IYZ)
            IP2=MIN(IVT,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE-2.0D0*(FP(IP)-EU*PREF(IP))
          END IF
C Add  2dtx ( Fvuyz-Et*Gvuyz )
          IF(ITABS.EQ.IXABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,IYZ)
            IP2=MIN(IVU,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+4.0D0*(FP(IP)-ET*PREF(IP))
          END IF
CGG.Nov03
          IF (ITUV.eq.IXYZ) THEN
            IDT=(ITABS*(ITABS+1))/2
            IDU=(IUABS*(IUABS+1))/2
            IDV=(IVABS*(IVABS+1))/2
            VALUE=VALUE+ipea_shift*0.5d0*BA(ISADR)*
     &                  (2.0d0-DREF(IDV)+DREF(IDT)+DREF(IDU))
          ENDIF
CGG End
          BA(ISADR)=VALUE
        END DO
      END DO
      END SUBROUTINE MKBA_DP

      SUBROUTINE MKBA_F3(ISYM,BA,NG3,F3,idxG3)
      USE SUPERINDEX
      use EQSOLV
      IMPLICIT REAL*8 (A-H,O-Z)
#include "caspt2.fh"

      REAL*8 BA(*)
      REAL*8 F3(NG3)
      INTEGER*1 idxG3(6,NG3)

C-SVC20100831: determine indices in SA where a certain F3 value will end up
      DO iG3=1,NG3
        iT=idxG3(1,iG3)
        iU=idxG3(2,iG3)
        iV=idxG3(3,iG3)
        iX=idxG3(4,iG3)
        iY=idxG3(5,iG3)
        iZ=idxG3(6,iG3)
        iST=IASYM(iT)
        iSU=IASYM(iU)
        iSV=IASYM(iV)
        iSX=IASYM(iX)
        iSY=IASYM(iY)
        iSZ=IASYM(iZ)
        ituvs=MUL(IST,MUL(ISU,ISV))
        ixyzs=MUL(ISX,MUL(ISY,ISZ))
        if(ituvs.ne.ixyzs) goto 500
        iTU=iT+NASHT*(iU-1)
        iVX=iV+NASHT*(iX-1)
        iYZ=iY+NASHT*(iZ-1)
        F3VAL=-F3(iG3)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - F(tuvxyz) -> BA(xut,vyz)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iV,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
C  - F(vxtuyz) -> BA(uxv,tyz)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iT,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(yzvxtu) -> BA(xzy,vtu)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iV,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(tuyzvx) -> BA(zut,yvx)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iY,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
 200   CONTINUE
C  - F(yztuvx) -> BA(uzy,tvx)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iT,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(vxyztu) -> BA(zxv,ytu)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iY,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
 300   CONTINUE
        if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 500
        if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 500
        if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 500
        if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 500
C  - F(utxvzy) -> BA(vtu,xzy)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iX,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 500
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
C  - F(xvutzy) -> BA(tvx,uzy)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iU,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(zyxvut) -> BA(vyz,xut)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iX,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(utzyxv) -> BA(ytu,zxv)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
 400   CONTINUE
C  - F(zyutxv) -> BA(tyz,uxv)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iU,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(xvzyut) -> BA(yvx,zut)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BA(ISADR)=BA(ISADR)+F3VAL
          END IF
        ENDIF
 500   CONTINUE
      END DO

      END SUBROUTINE MKBA_F3

#ifdef _MOLCAS_MPP_
      SUBROUTINE MKBA_F3_MPP(ISYM,BA,iLo,iHi,jLo,jHi,LDA,
     &                       NG3,F3,idxG3)
      USE MPI
      USE SUPERINDEX
      use stdalloc, only: mma_MaxDBLE
      use EQSOLV
      use definitions, only: MPIInt,RtoB,wp
      IMPLICIT REAL*8 (A-H,O-Z)
#include "caspt2.fh"

#include "global.fh"
#include "mafdecls.fh"

      DIMENSION BA(LDA,*)
      DIMENSION F3(NG3)
      INTEGER*1 idxG3(6,NG3)

      integer(kind=MPIInt), ALLOCATABLE :: SCOUNTS(:), RCOUNTS(:)
      integer(kind=MPIInt), ALLOCATABLE :: SCOUNTS2(:), RCOUNTS2(:)
      integer(kind=MPIInt), ALLOCATABLE :: SDISPLS(:), RDISPLS(:)
      integer(kind=MPIInt), ALLOCATABLE :: SDISPLS2(:), RDISPLS2(:)

      integer(kind=MPIInt), ALLOCATABLE :: SENDIDX(:), RECVIDX(:)
      REAL*8,    ALLOCATABLE :: SENDVAL(:), RECVVAL(:)

      integer(kind=MPIInt), PARAMETER :: ONE4=1, TWO4=2
      integer(kind=MPIInt) :: IERROR4
      INTEGER, PARAMETER :: I4=KIND(ONE4)

      INTEGER, ALLOCATABLE :: IBUF(:)

#include "mpi_interfaces.fh"

      ! Since we are stuck with collective calls to MPI_Alltoallv in
      ! order to gather the elements, each process needs to loop over
      ! the same number of blocks.
      NG3MAX=NG3
      CALL GAIGOP_SCAL(NG3MAX,'max')
      IF (NG3MAX.EQ.0) RETURN

      ! basic information
      NPROCS=GA_NNODES()
      MYRANK=GA_NODEID()

      ALLOCATE(SCOUNTS(NPROCS))
      ALLOCATE(RCOUNTS(NPROCS))
      ALLOCATE(SCOUNTS2(NPROCS))
      ALLOCATE(RCOUNTS2(NPROCS))
      ALLOCATE(SDISPLS(NPROCS))
      ALLOCATE(RDISPLS(NPROCS))
      ALLOCATE(SDISPLS2(NPROCS))
      ALLOCATE(RDISPLS2(NPROCS))

      ALLOCATE(IBUF(NPROCS))

      ! The global SA matrix has already been allocated, so we need to
      ! find out how much memory is left for buffering (4 equally sized
      ! buffers for sending and receiving values and indices)
      CALL mma_MaxDBLE(MAXMEM)
      iscal = (MPIInt*4 + wp*2)/RtoB
      MAXBUF=MIN(NINT(0.95D0*MAXMEM)/iscal,2000000000/8)

      ! Loop over blocks NG3B of NG3, so that 12*NG3B < MAXBUF/NPROCS.
      ! This guarantees that e.g. if all processes send all their data
      ! to one other, that process receives NPROCS*NG3B*12 elements
      ! in the receive buffer.
      NG3B=MAXBUF/(NPROCS*12)
      NG3B=MIN(NG3B,NG3MAX)
      CALL GAIGOP_SCAL(NG3B,'min')
      NBUF=12*NG3B

      ALLOCATE(SENDVAL(NBUF))
      ALLOCATE(SENDIDX(2*NBUF))

      ! Finally, we need some info on the layout of the global array in
      ! order to compute the process row of the row index.
      NAS=jHi-jLo+1
      NQOT=NAS/NPROCS
      NREM=NAS-NPROCS*NQOT

      NBLOCKS=(NG3MAX-1)/NG3B+1
      DO IBLOCK=1,NBLOCKS
        IG3STA=1+(IBLOCK-1)*NG3B
        IG3END=MIN(IG3STA+NG3B-1,NG3)

        SCOUNTS=0
        ! First pass to determine how many values will need to be sent
        ! to other processes. This is necessary to be able to allocate
        ! the buffer size and offsets.
        DO iG3=IG3STA,IG3END
          iT=idxG3(1,iG3)
          iU=idxG3(2,iG3)
          iV=idxG3(3,iG3)
          iX=idxG3(4,iG3)
          iY=idxG3(5,iG3)
          iZ=idxG3(6,iG3)
          iST=IASYM(iT)
          iSU=IASYM(iU)
          iSV=IASYM(iV)
          iSX=IASYM(iX)
          iSY=IASYM(iY)
          iSZ=IASYM(iZ)
          ituvs=MUL(IST,MUL(ISU,ISV))
          ixyzs=MUL(ISX,MUL(ISY,ISZ))
          if(ituvs.ne.ixyzs) CYCLE
          iTU=iT+NASHT*(iU-1)
          iVX=iV+NASHT*(iX-1)
          iYZ=iY+NASHT*(iZ-1)
          ! There are 12 equivalent cases, of which the second half
          ! reflects the S(tuv,xyz)=S(xyz,tuv) symmetry.

          ! - F(tuvxyz) -> BA(xut,vyz)
          jSYM=MUL(iSX,MUL(iSU,iST))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
          ! - F(vxtuyz) -> BA(uxv,tyz)
          jSYM=MUL(iSU,MUL(iSX,iSV))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - F(yzvxtu) -> BA(xzy,vtu)
          jSYM=MUL(iSX,MUL(iSZ,iSY))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - F(tuyzvx) -> BA(zut,yvx)
          jSYM=MUL(iSZ,MUL(iSU,iST))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 200      CONTINUE
          ! - F(yztuvx) -> BA(uzy,tvx)
          jSYM=MUL(iSU,MUL(iSZ,iSY))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - F(vxyztu) -> BA(zxv,ytu)
          jSYM=MUL(iSZ,MUL(iSX,iSV))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 300      CONTINUE
          if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) CYCLE
          if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) CYCLE
          if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) CYCLE
          if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) CYCLE
          ! - F(utxvzy) -> BA(vtu,xzy)
          jSYM=MUL(iSV,MUL(iST,iSU))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) CYCLE
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
          ! - F(xvutzy) -> BA(tvx,uzy)
          jSYM=MUL(iST,MUL(iSV,iSX))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - F(zyxvut) -> BA(vyz,xut)
          jSYM=MUL(iSV,MUL(iSY,iSZ))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - F(utzyxv) -> BA(ytu,zxv)
          jSYM=MUL(iSY,MUL(iST,iSU))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 400      CONTINUE
          ! - F(zyutxv) -> BA(tyz,uxv)
          jSYM=MUL(iST,MUL(iSY,iSZ))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - F(xvzyut) -> BA(yvx,zut)
          jSYM=MUL(iSY,MUL(iSV,iSX))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
        END DO

        ! At this point, SCOUNTS contains the number of values generated
        ! for each process. Use them to determine the send offsets.
        IOFFSET=0
        DO I=1,NPROCS
          SDISPLS(I)=INT(IOFFSET,I4)
          IBUF(I)=IOFFSET
          IOFFSET=IOFFSET+SCOUNTS(I)
        END DO

        ! Second pass fills the buffers with values and indices
        DO iG3=IG3STA,IG3END
          iT=idxG3(1,iG3)
          iU=idxG3(2,iG3)
          iV=idxG3(3,iG3)
          iX=idxG3(4,iG3)
          iY=idxG3(5,iG3)
          iZ=idxG3(6,iG3)
          iST=IASYM(iT)
          iSU=IASYM(iU)
          iSV=IASYM(iV)
          iSX=IASYM(iX)
          iSY=IASYM(iY)
          iSZ=IASYM(iZ)
          ituvs=MUL(IST,MUL(ISU,ISV))
          ixyzs=MUL(ISX,MUL(ISY,ISZ))
          if(ituvs.ne.ixyzs) CYCLE
          iTU=iT+NASHT*(iU-1)
          iVX=iV+NASHT*(iX-1)
          iYZ=iY+NASHT*(iZ-1)
          F3VAL=-F3(iG3)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - F(tuvxyz) -> BA(xut,vyz)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 301
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 201
C  - F(vxtuyz) -> BA(uxv,tyz)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(yzvxtu) -> BA(xzy,vtu)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iV,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(tuyzvx) -> BA(zut,yvx)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iY,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
 201      CONTINUE
C  - F(yztuvx) -> BA(uzy,tvx)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iT,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(vxyztu) -> BA(zxv,ytu)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iY,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
 301      CONTINUE
          if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) CYCLE
          if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) CYCLE
          if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) CYCLE
          if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) CYCLE
C  - F(utxvzy) -> BA(vtu,xzy)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) CYCLE
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 401
C  - F(xvutzy) -> BA(tvx,uzy)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(zyxvut) -> BA(vyz,xut)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iX,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(utzyxv) -> BA(ytu,zxv)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
 401      CONTINUE
C  - F(zyutxv) -> BA(tyz,uxv)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iU,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(xvzyut) -> BA(yvx,zut)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
        END DO

        ! Now we need to determine the receive counts.
        CALL MPI_ALLTOALL(SCOUNTS, ONE4, MPI_INTEGER,
     &                    RCOUNTS, ONE4, MPI_INTEGER,
     &                    MPI_COMM_WORLD, IERROR4)

        IOFFSET=0
        DO I=1,NPROCS
          RDISPLS(I)=INT(IOFFSET,I4)
          IOFFSET=IOFFSET+RCOUNTS(I)
          SCOUNTS2(I)=TWO4*SCOUNTS(I)
          RCOUNTS2(I)=TWO4*RCOUNTS(I)
          SDISPLS2(I)=TWO4*SDISPLS(I)
          RDISPLS2(I)=TWO4*RDISPLS(I)
        END DO
        NRECV=IOFFSET

        ALLOCATE(RECVVAL(NRECV))
        ALLOCATE(RECVIDX(2*NRECV))

        ! Now, it is time to collect the appropriate values and indices
        ! in their respective receive buffers.
        CALL MPI_ALLTOALLV(SENDVAL, SCOUNTS, SDISPLS, MPI_REAL8,
     &                     RECVVAL, RCOUNTS, RDISPLS, MPI_REAL8,
     &                     MPI_COMM_WORLD, IERROR4)
        CALL MPI_ALLTOALLV(SENDIDX, SCOUNTS2, SDISPLS2, MPI_INTEGER,
     &                     RECVIDX, RCOUNTS2, RDISPLS2, MPI_INTEGER,
     &                     MPI_COMM_WORLD, IERROR4)

        ! Finally, fill the local chunk of the SA matrix (block of rows)
        ! with the received values at their appropriate place.
        DO I=1,NRECV
          ISUP=RECVIDX(2*I-1)
          JSUP=RECVIDX(2*I)
          BA(ISUP-ILO+1,JSUP-JLO+1)=
     &      BA(ISUP-ILO+1,JSUP-JLO+1)+RECVVAL(I)
        END DO

        DEALLOCATE(RECVVAL)
        DEALLOCATE(RECVIDX)

      END DO ! end loop over blocks of F3 values

      DEALLOCATE(SENDVAL)
      DEALLOCATE(SENDIDX)

      DEALLOCATE(SCOUNTS)
      DEALLOCATE(RCOUNTS)
      DEALLOCATE(SCOUNTS2)
      DEALLOCATE(RCOUNTS2)
      DEALLOCATE(SDISPLS)
      DEALLOCATE(RDISPLS)
      DEALLOCATE(SDISPLS2)
      DEALLOCATE(RDISPLS2)

      DEALLOCATE(IBUF)

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL UNUSED_INTEGER(iHi)

      CONTAINS

      PURE INTEGER FUNCTION IPROW(IROW,NQOT,NREM)
      INTEGER, INTENT(IN) :: IROW, NQOT, NREM
      INTEGER :: TMP
      TMP=IROW-NREM*(NQOT+1)
      IF (TMP.GT.0) THEN
        IPROW=(TMP-1)/NQOT+NREM+1
      ELSE
        IPROW=(IROW-1)/(NQOT+1)+1
      END IF
      END FUNCTION

      END SUBROUTINE MKBA_F3_MPP
#endif

********************************************************************************
* Case C (ICASE=4)
********************************************************************************
      SUBROUTINE MKBC(DREF,NDREF,PREF,NPREF,FD,FP,NG3,F3,idxG3)
      USE SUPERINDEX
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      use EQSOLV
      IMPLICIT REAL*8 (A-H,O-Z)
#include "caspt2.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      INTEGER NDREF, NPREF, NG3
      Real*8 DREF(NDREF),PREF(NPREF),F3(NG3)
      Real*8 FD(NDREF),FP(NPREF)
      INTEGER*1 idxG3(6,NG3)
#ifdef _MOLCAS_MPP_
      Real*8 Dummy(1)
#endif

      ICASE=4
C LONG loop over superindex symmetry.
      DO ISYM=1,NSYM
        NIN=NINDEP(ISYM,ICASE)
        IF(NIN.EQ.0) CYCLE
        NAS=NTUV(ISYM)
        NBC=(NAS*(NAS+1))/2
        IF(NBC.LE.0) CYCLE

C Set up the matrix BC(tuv,xyz) defined by the expression
C <atuv|H0-E0|cxyz> = dac ( alpha(a) SC(tuv,xyz) + BC(tuv,xyz) )
C Formula used:
C    BC(tuv,xyz)
C    = Fvutxyz +dyu Fvztx + dyx Fvutz + dtu Fvxyz + dtu dyx Fvz
C    +(Ey+Eu-EASUM)*SC(tuv,xyz)
C    -Eu*( dyu Gvztx + dtu Gvxyz )
C    -Ey dyx Gvutz
C    -(Eu+Ey)*( dtu dyx Gvz )

C where dyu = Kronecker(y,u) etc. Gvutxyz=<Evutxyz>, etc.
C Similarly, Fvutxyz= Sum(w)(EPSA(w)<Evutxyzww>, etc.

        CALL PSBMAT_GETMEM('BC',lg_BC,NAS)
        CALL PSBMAT_READ('S',iCase,iSym,lg_BC,NAS)

        ! fill in the 3-el parts
#ifdef _MOLCAS_MPP_
        IF (IS_REAL_PAR()) THEN
          MYRANK = GA_NODEID()
          CALL GA_DISTRIBUTION (LG_BC,MYRANK,ILO,IHI,JLO,JHI)
          IF (JLO.NE.0 .AND. (JHI-JLO+1).NE.NAS) THEN
            WRITE(6,*) 'MKBC: MISMATCH IN RANGE OF THE SUPERINDICES'
            CALL ABEND()
          END IF
          IF (ILO.GT.0 .AND. JLO.GT.0) THEN
            CALL GA_ACCESS (LG_BC,ILO,IHI,JLO,JHI,MA,LDA)
            CALL MKBC_DP(DREF,NDREF,PREF,NPREF,FD,FP,iSYM,
     &                   DBL_MB(MA),ILO,IHI,JLO,JHI,LDA)
            CALL MKBC_F3_MPP(ISYM,DBL_MB(MA),ILO,IHI,JLO,JHI,LDA,
     &                       NG3,F3,IDXG3)
            CALL GA_RELEASE_UPDATE (LG_BC,ILO,IHI,JLO,JHI)
          ELSE
            CALL MKBC_F3_MPP(ISYM,DUMMY,ILO,IHI,JLO,JHI,LDA,
     &                       NG3,F3,IDXG3)
          END IF
        ELSE
#endif
          CALL MKBC_DP(DREF,NDREF,PREF,NPREF,FD,FP,
     &                 ISYM,GA_Arrays(lg_BC)%A(:),1,NAS,1,NAS,0)
          CALL MKBC_F3(ISYM,GA_Arrays(lg_BC)%A(:),NG3,F3,IDXG3)

#ifdef _MOLCAS_MPP_
        END IF
#endif

        CALL PSBMAT_WRITE('B',iCase,iSYM,lg_BC,NAS)

        IF(IPRGLB.GE.DEBUG) THEN
          DBC=PSBMAT_FPRINT(lg_BC,NAS)
          WRITE(6,'("DEBUG> ",A4,1X,I3,1X,ES21.14)') 'C', ISYM, DBC
        END IF

        CALL PSBMAT_FREEMEM(lg_BC)
      END DO

      END SUBROUTINE MKBC

      SUBROUTINE MKBC_DP (DREF,NDREF,PREF,NPREF,FD,FP,iSYM,
     &                    BC,iLo,iHi,jLo,jHi,LDC)
      USE SUPERINDEX
      use caspt2_global, only:ipea_shift
      use EQSOLV
      IMPLICIT REAL*8 (A-H,O-Z)
#include "caspt2.fh"
      INTEGER NDREF,NPREF, iSYM,iLo,iHi,jLo,jHi,LDC
      REAL*8 DREF(NDREF),PREF(NPREF)
      REAL*8 FD(NDREF),FP(NPREF)
      REAL*8 BC(*)

      DO IXYZ=jLo,jHi
        IXYZABS=IXYZ+NTUVES(ISYM)
        IXABS=MTUV(1,IXYZABS)
        IYABS=MTUV(2,IXYZABS)
        IZABS=MTUV(3,IXYZABS)
        EY=EPSA(IYABS)
        DO ITUV=iLo,iHi
          ITUVABS=ITUV+NTUVES(ISYM)
          ITABS=MTUV(1,ITUVABS)
          IUABS=MTUV(2,ITUVABS)
          IVABS=MTUV(3,ITUVABS)
          EU=EPSA(IUABS)
          EYU=EY + EU
          FACT=EYU-EASUM
          ISADR=1+iTUV-iLo+LDC*(iXYZ-jLo)
          IF (LDC.EQ.0) THEN
            IF (iXYZ.LE.iTUV) THEN
              ISADR=(ITUV*(ITUV-1))/2+IXYZ
            ELSE
              CYCLE
            END IF
          END IF
          VALUE=FACT*BC(ISADR)
C VALUE= Fvutxyz + (EPSA(y)+EPSA(u))*SC(tuv,xyz)
C Add  dyu ( Fvztx - EPSA(u)*Gvztx )
          IF(IYABS.EQ.IUABS) THEN
            IVZ=IVABS+NASHT*(IZABS-1)
            ITX=ITABS+NASHT*(IXABS-1)
            IP1=MAX(IVZ,ITX)
            IP2=MIN(IVZ,ITX)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+2.0D0*(FP(IP)-EU*PREF(IP))
          END IF
C Add  dyx ( Fvutz - EPSA(y)*Gvutz )
          IF(IYABS.EQ.IXABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            ITZ=ITABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,ITZ)
            IP2=MIN(IVU,ITZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+2.0D0*(FP(IP)-EY*PREF(IP))
          END IF
C Add  dtu ( Fvxyz - EPSA(u)*Gvxyz + dyx Fvz -
C             (EPSA(u)+EPSA(y)*dyz Gvz)
          IF(ITABS.EQ.IUABS) THEN
            IVX=IVABS+NASHT*(IXABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVX,IYZ)
            IP2=MIN(IVX,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+2.0D0*(FP(IP)-EU*PREF(IP))
            IF(IYABS.EQ.IXABS) THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE+FD(ID)-EYU*DREF(ID)
            END IF
          END IF
CGG.Nov03
          IF (ITUV.eq.IXYZ) THEN
            IDT=(ITABS*(ITABS+1))/2
            IDU=(IUABS*(IUABS+1))/2
            IDV=(IVABS*(IVABS+1))/2
            VALUE=VALUE+ipea_shift*0.5d0*BC(ISADR)*
     &                    (4.0d0-DREF(IDT)-DREF(IDV)+DREF(IDU))
          ENDIF
CGG End
          BC(ISADR)=VALUE
        END DO
      END DO
      END SUBROUTINE MKBC_DP

      SUBROUTINE MKBC_F3(ISYM,BC,NG3,F3,idxG3)
      USE SUPERINDEX
      use EQSOLV
      IMPLICIT REAL*8 (A-H,O-Z)
#include "caspt2.fh"

      DIMENSION BC(*)
      DIMENSION F3(NG3)
      INTEGER*1 idxG3(6,NG3)

C-SVC20100831: determine indices in BC where a certain G3 value will end up
      DO iG3=1,NG3
        iT=idxG3(1,iG3)
        iU=idxG3(2,iG3)
        iV=idxG3(3,iG3)
        iX=idxG3(4,iG3)
        iY=idxG3(5,iG3)
        iZ=idxG3(6,iG3)
        iST=IASYM(iT)
        iSU=IASYM(iU)
        iSV=IASYM(iV)
        iSX=IASYM(iX)
        iSY=IASYM(iY)
        iSZ=IASYM(iZ)
        ituvs=MUL(IST,MUL(ISU,ISV))
        ixyzs=MUL(ISX,MUL(ISY,ISZ))
        if(ituvs.ne.ixyzs) goto 500
        iTU=iT+NASHT*(iU-1)
        iVX=iV+NASHT*(iX-1)
        iYZ=iY+NASHT*(iZ-1)
        F3VAL=F3(iG3)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - F(tuvxyz) -> BC(vut,xyz)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iX,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
C  - F(vxtuyz) -> BC(txv,uyz)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iU,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(yzvxtu) -> BC(vzy,xtu)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iX,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(tuyzvx) -> BC(yut,zvx)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
 200   CONTINUE
C  - F(yztuvx) -> BC(tzy,uvx)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iU,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(vxyztu) -> BC(yxv,ztu)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
 300   CONTINUE
        if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 500
        if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 500
        if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 500
        if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 500
C  - F(utxvzy) -> BC(xtu,vzy)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iV,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 500
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
C  - F(xvutzy) -> BC(uvx,tzy)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iT,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(zyxvut) -> BC(xyz,vut)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iV,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(utzyxv) -> BC(ztu,yxv)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iY,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
 400   CONTINUE
C  - F(zyutxv) -> BC(uyz,txv)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iT,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
C  - F(xvzyut) -> BC(zvx,yut)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iY,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            BC(ISADR)=BC(ISADR)+F3VAL
          END IF
        ENDIF
 500   CONTINUE
      END DO

      END SUBROUTINE MKBC_F3

#ifdef _MOLCAS_MPP_
      SUBROUTINE MKBC_F3_MPP(ISYM,BC,iLo,iHi,jLo,jHi,LDC,
     &                       NG3,F3,idxG3)
      USE MPI
      USE SUPERINDEX
      use stdalloc, only: mma_MaxDBLE
      use EQSOLV
      use definitions, only: MPIInt,RtoB,wp
      IMPLICIT REAL*8 (A-H,O-Z)
#include "caspt2.fh"

#include "global.fh"
#include "mafdecls.fh"

      DIMENSION BC(LDC,*)
      DIMENSION F3(NG3)
      INTEGER*1 idxG3(6,NG3)

      integer(kind=MPIInt), ALLOCATABLE :: SCOUNTS(:), RCOUNTS(:)
      integer(kind=MPIInt), ALLOCATABLE :: SCOUNTS2(:), RCOUNTS2(:)
      integer(kind=MPIInt), ALLOCATABLE :: SDISPLS(:), RDISPLS(:)
      integer(kind=MPIInt), ALLOCATABLE :: SDISPLS2(:), RDISPLS2(:)

      integer(kind=MPIInt), ALLOCATABLE :: SENDIDX(:), RECVIDX(:)
      REAL*8,    ALLOCATABLE :: SENDVAL(:), RECVVAL(:)

      integer(kind=MPIInt), PARAMETER :: ONE4=1, TWO4=2
      integer(kind=MPIInt) :: IERROR4
      INTEGER, PARAMETER :: I4=KIND(ONE4)

      INTEGER, ALLOCATABLE :: IBUF(:)

#include "mpi_interfaces.fh"

      ! Since we are stuck with collective calls to MPI_Alltoallv in
      ! order to gather the elements, each process needs to loop over
      ! the same number of blocks.
      NG3MAX=NG3
      CALL GAIGOP_SCAL(NG3MAX,'max')
      IF (NG3MAX.EQ.0) RETURN

      ! basic information
      NPROCS=GA_NNODES()
      MYRANK=GA_NODEID()

      ALLOCATE(SCOUNTS(NPROCS))
      ALLOCATE(RCOUNTS(NPROCS))
      ALLOCATE(SCOUNTS2(NPROCS))
      ALLOCATE(RCOUNTS2(NPROCS))
      ALLOCATE(SDISPLS(NPROCS))
      ALLOCATE(RDISPLS(NPROCS))
      ALLOCATE(SDISPLS2(NPROCS))
      ALLOCATE(RDISPLS2(NPROCS))

      ALLOCATE(IBUF(NPROCS))

      ! The global SC matrix has already been allocated, so we need to
      ! find out how much memory is left for buffering (4 equally sized
      ! buffers for sending and receiving values and indices)
      CALL mma_MaxDBLE(MAXMEM)
      iscal = (MPIInt*4 + wp*2)/RtoB
      MAXBUF=MIN(NINT(0.95D0*MAXMEM)/iscal,2000000000/8)

      ! Loop over blocks NG3B of NG3, so that 12*NG3B < MAXBUF/NPROCS.
      ! This guarantees that e.g. if all processes send all their data
      ! to one other, that process receives NPROCS*NG3B*12 elements
      ! in the receive buffer.
      NG3B=MAXBUF/(NPROCS*12)
      NG3B=MIN(NG3B,NG3MAX)
      CALL GAIGOP_SCAL(NG3B,'min')
      NBUF=12*NG3B

      ALLOCATE(SENDVAL(NBUF))
      ALLOCATE(SENDIDX(2*NBUF))

      ! Finally, we need some info on the layout of the global array in
      ! order to compute the process row of the row index.
      NAS=jHi-jLo+1
      NQOT=NAS/NPROCS
      NREM=NAS-NPROCS*NQOT

      NBLOCKS=(NG3MAX-1)/NG3B+1
      DO IBLOCK=1,NBLOCKS
        IG3STA=1+(IBLOCK-1)*NG3B
        IG3END=MIN(IG3STA+NG3B-1,NG3)

        SCOUNTS=0
        ! First pass to determine how many values will need to be sent
        ! to other processes. This is necessary to be able to allocate
        ! the buffer size and offsets.
        DO iG3=IG3STA,IG3END
          iT=idxG3(1,iG3)
          iU=idxG3(2,iG3)
          iV=idxG3(3,iG3)
          iX=idxG3(4,iG3)
          iY=idxG3(5,iG3)
          iZ=idxG3(6,iG3)
          iST=IASYM(iT)
          iSU=IASYM(iU)
          iSV=IASYM(iV)
          iSX=IASYM(iX)
          iSY=IASYM(iY)
          iSZ=IASYM(iZ)
          ituvs=MUL(IST,MUL(ISU,ISV))
          ixyzs=MUL(ISX,MUL(ISY,ISZ))
          if(ituvs.ne.ixyzs) goto 500
          iTU=iT+NASHT*(iU-1)
          iVX=iV+NASHT*(iX-1)
          iYZ=iY+NASHT*(iZ-1)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the B(tuv,xyz)=B(xyz,tuv) symmetry:
C  - F(tuvxyz) -> BC(vut,xyz)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
C  - F(vxtuyz) -> BC(txv,uyz)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - F(yzvxtu) -> BC(vzy,xtu)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - F(tuyzvx) -> BC(yut,zvx)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 200   CONTINUE
C  - F(yztuvx) -> BC(tzy,uvx)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - F(vxyztu) -> BC(yxv,ztu)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 300   CONTINUE
          if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 500
          if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 500
          if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 500
          if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 500
C  - F(utxvzy) -> BC(xtu,vzy)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 500
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
C  - F(xvutzy) -> BC(uvx,tzy)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - F(zyxvut) -> BC(xyz,vut)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - F(utzyxv) -> BC(ztu,yxv)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 400   CONTINUE
C  - F(zyutxv) -> BC(uyz,txv)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - F(xvzyut) -> BC(zvx,yut)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 500   CONTINUE
        END DO

        ! At this point, SCOUNTS contains the number of values generated
        ! for each process. Use them to determine the send offsets.
        IOFFSET=0
        DO I=1,NPROCS
          SDISPLS(I)=INT(IOFFSET,I4)
          IBUF(I)=IOFFSET
          IOFFSET=IOFFSET+SCOUNTS(I)
        END DO

        ! Second pass fills the buffers with values and indices
        DO iG3=IG3STA,IG3END
          iT=idxG3(1,iG3)
          iU=idxG3(2,iG3)
          iV=idxG3(3,iG3)
          iX=idxG3(4,iG3)
          iY=idxG3(5,iG3)
          iZ=idxG3(6,iG3)
          iST=IASYM(iT)
          iSU=IASYM(iU)
          iSV=IASYM(iV)
          iSX=IASYM(iX)
          iSY=IASYM(iY)
          iSZ=IASYM(iZ)
          ituvs=MUL(IST,MUL(ISU,ISV))
          ixyzs=MUL(ISX,MUL(ISY,ISZ))
          if(ituvs.ne.ixyzs) goto 501
          iTU=iT+NASHT*(iU-1)
          iVX=iV+NASHT*(iX-1)
          iYZ=iY+NASHT*(iZ-1)
          F3VAL=F3(iG3)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - F(tuvxyz) -> BC(vut,xyz)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iX,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 301
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 201
C  - F(vxtuyz) -> BC(txv,uyz)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iU,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(yzvxtu) -> BC(vzy,xtu)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iX,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(tuyzvx) -> BC(yut,zvx)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
 201   CONTINUE
C  - F(yztuvx) -> BC(tzy,uvx)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iU,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(vxyztu) -> BC(yxv,ztu)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
 301   CONTINUE
          if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 501
          if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 501
          if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 501
          if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 501
C  - F(utxvzy) -> BC(xtu,vzy)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iV,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 501
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 401
C  - F(xvutzy) -> BC(uvx,tzy)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iT,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(zyxvut) -> BC(xyz,vut)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iV,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(utzyxv) -> BC(ztu,yxv)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iY,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
 401   CONTINUE
C  - F(zyutxv) -> BC(uyz,txv)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iT,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - F(xvzyut) -> BC(zvx,yut)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iY,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=F3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
 501   CONTINUE
        END DO

        ! Now we need to determine the receive counts.
        CALL MPI_ALLTOALL(SCOUNTS, ONE4, MPI_INTEGER,
     &                    RCOUNTS, ONE4, MPI_INTEGER,
     &                    MPI_COMM_WORLD, IERROR4)

        IOFFSET=0
        DO I=1,NPROCS
          RDISPLS(I)=INT(IOFFSET,I4)
          IOFFSET=IOFFSET+RCOUNTS(I)
          SCOUNTS2(I)=TWO4*SCOUNTS(I)
          RCOUNTS2(I)=TWO4*RCOUNTS(I)
          SDISPLS2(I)=TWO4*SDISPLS(I)
          RDISPLS2(I)=TWO4*RDISPLS(I)
        END DO
        NRECV=IOFFSET

        ALLOCATE(RECVVAL(NRECV))
        ALLOCATE(RECVIDX(2*NRECV))

        ! Now, it is time to collect the appropriate values and indices
        ! in their respective receive buffers.
        CALL MPI_ALLTOALLV(SENDVAL, SCOUNTS, SDISPLS, MPI_REAL8,
     &                     RECVVAL, RCOUNTS, RDISPLS, MPI_REAL8,
     &                     MPI_COMM_WORLD, IERROR4)
        CALL MPI_ALLTOALLV(SENDIDX, SCOUNTS2, SDISPLS2, MPI_INTEGER,
     &                     RECVIDX, RCOUNTS2, RDISPLS2, MPI_INTEGER,
     &                     MPI_COMM_WORLD, IERROR4)

        ! Finally, fill the local chunk of the SC matrix (block of rows)
        ! with the received values at their appropriate place.
        DO I=1,NRECV
          ISUP=RECVIDX(2*I-1)
          JSUP=RECVIDX(2*I)
          BC(ISUP-ILO+1,JSUP-JLO+1)=
     &      BC(ISUP-ILO+1,JSUP-JLO+1)+RECVVAL(I)
        END DO

        DEALLOCATE(RECVVAL)
        DEALLOCATE(RECVIDX)

      END DO ! end loop over blocks of G3 values

      DEALLOCATE(SENDVAL)
      DEALLOCATE(SENDIDX)

      DEALLOCATE(SCOUNTS)
      DEALLOCATE(RCOUNTS)
      DEALLOCATE(SCOUNTS2)
      DEALLOCATE(RCOUNTS2)
      DEALLOCATE(SDISPLS)
      DEALLOCATE(RDISPLS)
      DEALLOCATE(SDISPLS2)
      DEALLOCATE(RDISPLS2)

      DEALLOCATE(IBUF)
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL UNUSED_INTEGER(iHi)

      CONTAINS

      PURE INTEGER FUNCTION IPROW(IROW,NQOT,NREM)
      INTEGER, INTENT(IN) :: IROW, NQOT, NREM
      INTEGER :: TMP
      TMP=IROW-NREM*(NQOT+1)
      IF (TMP.GT.0) THEN
        IPROW=(TMP-1)/NQOT+NREM+1
      ELSE
        IPROW=(IROW-1)/(NQOT+1)+1
      END IF
      END FUNCTION

      END SUBROUTINE MKBC_F3_MPP
#endif

      SUBROUTINE MKBB(DREF,NDREF,PREF,NPREF,FD,FP)
      USE SUPERINDEX
      use caspt2_global, only:ipea_shift
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)

#include "caspt2.fh"

      INTEGER NDREF,NPREF
      REAL*8 DREF(NDREF),PREF(NPREF)
      REAL*8 FD(NDREF),FP(NPREF)

      REAL*8, ALLOCATABLE:: BB(:), BBP(:), SP(:), SDP(:),
     &                             BBM(:), SM(:), SDM(:)

C Set up the matrices BBP(tu,xy) and BBM(tu,xy)
C Formulae used:
C    BB(tu,xy)= 2*( Fyuxt - (A-Et-Eu-Ex-Ey)*Gyuxt )
C      + 4*dxt ( (A-Et-Ey-Eu)*Dyu - Fyu)
C      + 4*dyu ( (A-Et-Ey-Ex)*Dxt - Fxt)
C      - 2*dyt ( (A-Et-Eu-Ex)*Dxu - Fxu)
C      - 2*dxu ( (A-Et-Eu-Ey)*Dyt - Fyt)
C      + 8*dxt*dyu (Et+Ey)
C      - 4*dxu*dyt (Et+Ex)
C where A= EASUM= sum over active w of (Ew*Dww).
C    BBP(tu,xy)=BB(tu,xy)+BB(tu,yx)
C    BBM(tu,xy)=BB(tu,xy)-BB(tu,yx)


C Loop over superindex symmetry.
      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,2)
        IF(NINP.EQ.0) GOTO 1000
        NAS=NTU(ISYM)
        NBB=(NAS*(NAS+1))/2
        IF(NBB.GT.0) THEN
          CALL mma_allocate(BB,NBB,Label='BB')
        END IF
        DO ITU=1,NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          ET=EPSA(ITABS)
          EU=EPSA(IUABS)
          DO IXY=1,ITU
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            EX=EPSA(IXABS)
            EY=EPSA(IYABS)
            IBADR=(ITU*(ITU-1))/2+IXY
            IXT=IXABS+NASHT*(ITABS-1)
            IYU=IYABS+NASHT*(IUABS-1)
            IP1=MAX(IXT,IYU)
            IP2=MIN(IXT,IYU)
            IP=(IP1*(IP1-1))/2+IP2
            ATUXY=EASUM-ET-EU-EX-EY
            VALUE=4.0D0*(FP(IP)-ATUXY*PREF(IP))
C Add  + 4*dxt ( (A-Et-Ey-Eu)*Dyu - Fyu)
            IF(IXABS.EQ.ITABS) THEN
              ID1=MAX(IYABS,IUABS)
              ID2=MIN(IYABS,IUABS)
              ID=(ID1*(ID1-1))/2+ID2
              ATYU=EASUM-ET-EY-EU
              VALUE=VALUE+4.0D0*(ATYU*DREF(ID)-FD(ID))
C Add  + 8*dxt*dyu (Et+Ey)
              IF(IYABS.EQ.IUABS) THEN
                VALUE=VALUE+8.0D0*(ET+EY)
              END IF
            END IF
C Add  + 4*dyu ( (A-Et-Ey-Ex)*Dxt - Fxt)
            IF(IYABS.EQ.IUABS) THEN
              ID1=MAX(IXABS,ITABS)
              ID2=MIN(IXABS,ITABS)
              ID=(ID1*(ID1-1))/2+ID2
              ATYX=EASUM-ET-EY-EX
              VALUE=VALUE+4.0D0*(ATYX*DREF(ID)-FD(ID))
            END IF
C Add  - 2*dyt ( (A-Et-Eu-Ex)*Dxu - Fxu)
            IF(IYABS.EQ.ITABS) THEN
              ID1=MAX(IXABS,IUABS)
              ID2=MIN(IXABS,IUABS)
              ID=(ID1*(ID1-1))/2+ID2
              ATUX=EASUM-ET-EU-EX
              VALUE=VALUE-2.0D0*(ATUX*DREF(ID)-FD(ID))
C Add  - 4*dxu*dyt (Et+Ex)
              IF(IXABS.EQ.IUABS) THEN
                VALUE=VALUE-4.0D0*(ET+EX)
              END IF
            END IF
C Add  - 2*dxu ( (A-Et-Eu-Ey)*Dyt - Fyt)
            IF(IXABS.EQ.IUABS) THEN
              ID1=MAX(IYABS,ITABS)
              ID2=MIN(IYABS,ITABS)
              ID=(ID1*(ID1-1))/2+ID2
              ATUY=EASUM-ET-EU-EY
              VALUE=VALUE-2.0D0*(ATUY*DREF(ID)-FD(ID))
            END IF
            BB(IBADR)=VALUE
          END DO
        END DO
        NASP=NTGEU(ISYM)
        NBBP=(NASP*(NASP+1))/2
        IF(NBBP.GT.0) THEN
          CALL mma_allocate(BBP,NBBP,Label='BBP')
CGG.Nov03  Load in SDP the diagonal elements of SBP matrix:
          NSP=(NASP*(NASP+1))/2
          CALL mma_allocate(SP,NSP,Label='SP')
          CALL mma_allocate(SDP,NASP,Label='SDP')
          IDSP=IDSMAT(ISYM,2)
          CALL DDAFILE(LUSBT,2,SP,NSP,IDSP)
          IDIAG=0
          DO I=1,NASP
            IDIAG=IDIAG+I
            SDP(I)=SP(IDIAG)
          END DO
          CALL mma_deallocate(SP)
CGG End
        END IF
        NASM=NTGTU(ISYM)
        NBBM=(NASM*(NASM+1))/2
        IF(NBBM.GT.0) THEN
          CALL mma_allocate(BBM,NBBM,Label='BBM')
CGG.Nov03  Load in SDM the diagonal elements of SBM matrix:
          NSM=(NASM*(NASM+1))/2
          CALL mma_allocate(SM,NSM,Label='SM')
          CALL mma_allocate(SDM,NASM,Label='SDM')
          IDSM=IDSMAT(ISYM,3)
          CALL DDAFILE(LUSBT,2,SM,NSM,IDSM)
          IDIAG=0
          DO I=1,NASM
            IDIAG=IDIAG+I
            SDM(I)=SM(IDIAG)
          END DO
          CALL mma_deallocate(SM)
CGG End
        END IF
        INSM=1
        DO ITGEU=1,NASP
          ITGEUABS=ITGEU+NTGEUES(ISYM)
          ITABS=MTGEU(1,ITGEUABS)
          IUABS=MTGEU(2,ITGEUABS)
          ITU=KTU(ITABS,IUABS)-NTUES(ISYM)
          DO IXGEY=1,ITGEU
            IXGEYABS=IXGEY+NTGEUES(ISYM)
            IXABS=MTGEU(1,IXGEYABS)
            IYABS=MTGEU(2,IXGEYABS)
            IXY=KTU(IXABS,IYABS)-NTUES(ISYM)
            IYX=KTU(IYABS,IXABS)-NTUES(ISYM)
            IF(ITU.GE.IXY) THEN
              IBADR=(ITU*(ITU-1))/2+IXY
            ELSE
              IBADR=(IXY*(IXY-1))/2+ITU
            END IF
            BTUXY=BB(IBADR)
            IF(ITU.GE.IYX) THEN
              IBADR=(ITU*(ITU-1))/2+IYX
            ELSE
              IBADR=(IYX*(IYX-1))/2+ITU
            END IF
            BTUYX=BB(IBADR)
            IBPADR=(ITGEU*(ITGEU-1))/2+IXGEY
            BBP(IBPADR)=BTUXY+BTUYX
CGG.Nov03
            IF (ITGEU.eq.IXGEY) THEN
              IDT=(ITABS*(ITABS+1))/2
              IDU=(IUABS*(IUABS+1))/2
              BBP(IBPADR)=BBP(IBPADR)+ipea_shift*0.5d0*
     &                          (DREF(IDT)+DREF(IDU))*SDP(ITGEU)
            ENDIF
CGG End
            IF(ITABS.EQ.IUABS) GOTO 200
            IF(IXABS.EQ.IYABS) GOTO 200
            ITGTU=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
            IXGTY=KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
            IBMADR=(ITGTU*(ITGTU-1))/2+IXGTY
            BBM(IBMADR)=BTUXY-BTUYX
CGG.Nov03
            IF (ITGEU.eq.IXGEY) THEN
              IDT=(ITABS*(ITABS+1))/2
              IDU=(IUABS*(IUABS+1))/2
              BBM(IBMADR)=BBM(IBMADR)+ipea_shift*0.5d0*
     &                           (DREF(IDT)+DREF(IDU))*SDM(INSM)
              INSM=INSM+1
            ENDIF
CGG.End
 200        CONTINUE
          END DO
        END DO
        IF(NBB.GT.0) CALL mma_deallocate(BB)

C Write to disk, and save size and address.
        IF(NBBP.GT.0.and.NINDEP(ISYM,2).GT.0) THEN
          IDISK=IDBMAT(ISYM,2)
          CALL DDAFILE(LUSBT,1,BBP,NBBP,IDISK)
          CALL mma_deallocate(BBP)
CGG.Nov03 DisAlloc SDP
          CALL mma_deallocate(SDP)
CGG End
        END IF
        IF(NBBM.GT.0) THEN
          IF(NINDEP(ISYM,3).GT.0) THEN
            IDISK=IDBMAT(ISYM,3)
            CALL DDAFILE(LUSBT,1,BBM,NBBM,IDISK)
          END IF
          CALL mma_deallocate(BBM)
CGG.Nov03 DisAlloc SDM
          CALL mma_deallocate(SDM)
CGG End
        END IF
 1000 CONTINUE
      END DO

      END SUBROUTINE MKBB

      SUBROUTINE MKBD(DREF,NDREF,PREF,NPREF,FD,FP)
      USE SUPERINDEX
      use caspt2_global, only:ipea_shift
      use caspt2_global, only:LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)

#include "caspt2.fh"

      INTEGER NDREF,NPREF
      REAL*8 DREF(NDREF),PREF(NPREF)
      REAL*8 FD(NDREF),FP(NPREF)

      REAL*8, ALLOCATABLE:: BD(:), S(:), SD(:)

C Set up the matrix BD(tuP,xyQ),P and Q are 1 or 2,
C Formulae used:
C    BD(tu1,xy1)=
C      2*Futxy + 2*(Ex+Et-A)*Gutxy + 2*dxt (Fuy + (Et-A)*Duy)
C    BD(tu2,xy1)= -BD(tu1,xy1)/2
C    BD(tu2,xy2)=
C       -Fxtuy - (Ex+Et-A)*Gxtuy + 2*dxt (Fuy + (Ex-A)*Duy)
C where A=EASUM=Sum(w) of (Ew*Dww)


C Loop over superindex symmetry.
      DO ISYM=1,NSYM
        NIN=NINDEP(ISYM,5)
        IF(NIN.EQ.0) GOTO 1000
        NAS=NTU(ISYM)
        NBD=(2*NAS*(2*NAS+1))/2
        IF(NBD.GT.0) THEN
          CALL mma_allocate(BD,NBD,Label='BD')
CGG.Nov03  Load in SD the diagonal elements of SD matrix:
          NS2=(2*NAS*(2*NAS+1))/2
          NAS2=2*NAS
          CALL mma_allocate(S,NS2,Label='S')
          CALL mma_allocate(SD,NAS2,Label='SD')
          IDS=IDSMAT(ISYM,5)
          CALL DDAFILE(LUSBT,2,S,NS2,IDS)
          IDIAG=0
          DO I=1,NAS2
            IDIAG=IDIAG+I
            SD(I)=S(IDIAG)
          END DO
          CALL mma_deallocate(S)
CGG End
        END IF
        DO ITU=1,NAS
          ITU2=ITU+NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          ET=EPSA(ITABS)
          DO IXY=1,ITU
            IXY2=IXY+NAS
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            EX=EPSA(IXABS)
            IB11=(ITU*(ITU-1))/2+IXY
            IB21=(ITU2*(ITU2-1))/2+IXY
            IB12=(IXY2*(IXY2-1))/2+ITU
            IB22=(ITU2*(ITU2-1))/2+IXY2
            IUTP=IUABS+NASHT*(ITABS-1)
            IXYP=IXABS+NASHT*(IYABS-1)
            IP1=MAX(IUTP,IXYP)
            IP2=MIN(IUTP,IXYP)
            IP=(IP1*(IP1-1))/2+IP2
            ETX=ET+EX
            B11=4.0D0*(FP(IP)+(ETX-EASUM)*PREF(IP))
            IXTP=IXABS+NASHT*(ITABS-1)
            IUYP=IUABS+NASHT*(IYABS-1)
            IP1=MAX(IXTP,IUYP)
            IP2=MIN(IXTP,IUYP)
            IP=(IP1*(IP1-1))/2+IP2
            B22=-2.0D0*(FP(IP)+(ETX-EASUM)*PREF(IP))
            IF(IXABS.EQ.ITABS) THEN
              ID1=MAX(IUABS,IYABS)
              ID2=MIN(IUABS,IYABS)
              ID=(ID1*(ID1-1))/2+ID2
              FUY=FD(ID)
              DUY=DREF(ID)
              B11=B11+2.0D0*(FUY+(ET-EASUM)*DUY)
              B22=B22+2.0D0*(FUY+(EX-EASUM)*DUY)
            END IF
            BD(IB11)= B11
            BD(IB21)=-0.5D0*B11
            BD(IB12)=-0.5D0*B11
            BD(IB22)= B22
CGG.Nov03
            IF (ITU.eq.IXY) THEN
              IDT=(ITABS*(ITABS+1))/2
              IDU=(IUABS*(IUABS+1))/2
              BD(IB11)=BD(IB11)+ipea_shift*0.5d0*
     &                       (2.0d0-DREF(IDU)+DREF(IDT))*SD(ITU)
              BD(IB22)=BD(IB22)+ipea_shift*0.5d0*
     &                   (2.0d0-DREF(IDU)+DREF(IDT))*SD(ITU+NAS)
            ENDIF
CGG End
          END DO
        END DO

C Write to disk
        IF(NBD.GT.0.and.NINDEP(ISYM,5).GT.0) THEN
          IDISK=IDBMAT(ISYM,5)
          CALL DDAFILE(LUSBT,1,BD,NBD,IDISK)
          CALL mma_deallocate(BD)
CGG.Nov03 DisAlloc SD
          CALL mma_deallocate(SD)
CGG End
        END IF
 1000 CONTINUE
      END DO


      END SUBROUTINE MKBD

      SUBROUTINE MKBE(DREF,NDREF,FD)
      use caspt2_global, only:ipea_shift
      use caspt2_global, only:LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)

#include "caspt2.fh"

      INTEGER NDREF
      REAL*8 DREF(NDREF),FD(NDREF)

      Real*8, ALLOCATABLE:: BE(:), S(:), SD(:)

C Set up the matrix BE(t,x)
C Formula used:
C    BE(t,x)=-Ftx + (EASUM-Ex-Et)*Dtx
C            + 2dtx Ex


      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,6)
        IF(NINP.EQ.0) GOTO 1000
        NINM=NINDEP(ISYM,7)
        NAS=NASH(ISYM)
        NBE=(NAS*(NAS+1))/2
        IF(NBE.GT.0) THEN
          CALL mma_allocate(BE,NBE,LABEL='BE')
CGG.Nov03  Load in SD the diagonal elements of SE matrix:
          NS=(NAS*(NAS+1))/2
          CALL mma_allocate(S,NS,Label='S')
          CALL mma_allocate(SD,NAS,Label='SD')
          IDS=IDSMAT(ISYM,6)
          CALL DDAFILE(LUSBT,2,S,NS,IDS)
          IDIAG=0
          DO I=1,NAS
            IDIAG=IDIAG+I
            SD(I)=S(IDIAG)
          END DO
          CALL mma_deallocate(S)
        ENDIF
CGG End
        DO IT=1,NAS
          ITABS=IT+NAES(ISYM)
          ET=EPSA(ITABS)
          DO IX=1,IT
            IXABS=IX+NAES(ISYM)
            EX=EPSA(IXABS)
            IBE=(IT*(IT-1))/2+IX
            ID=(ITABS*(ITABS-1))/2+IXABS
            VALUE=-FD(ID)+(EASUM-EX-ET)*DREF(ID)
            IF(ITABS.EQ.IXABS) THEN
              VALUE=VALUE+2.0D0*EX
            END IF
CGG.Nov03
            IF (IT.eq.IX) THEN
              IDT=(ITABS*(ITABS+1))/2
              VALUE=VALUE+ipea_shift*0.5d0*DREF(IDT)*SD(IT)
            ENDIF
CGG End
            BE(IBE)=VALUE
          END DO
        END DO

C Write to disk
        IF(NBE.GT.0.and.NINDEP(ISYM,6).GT.0) THEN
          IDISK=IDBMAT(ISYM,6)
          CALL DDAFILE(LUSBT,1,BE,NBE,IDISK)
          IF(NINM.GT.0.and.NINDEP(ISYM,7).GT.0) THEN
            IDISK=IDBMAT(ISYM,7)
            CALL DDAFILE(LUSBT,1,BE,NBE,IDISK)
          END IF
          CALL mma_deallocate(BE)
CGG.Nov03 DisAlloc SD
          CALL mma_deallocate(SD)
CGG End
        END IF
 1000 CONTINUE
      END DO

      END SUBROUTINE MKBE

      SUBROUTINE MKBF(DREF,NDREF,PREF,NPREF,FP)
      USE SUPERINDEX
      use caspt2_global, only:ipea_shift
      use caspt2_global, only:LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)

#include "caspt2.fh"

      INTEGER NDREF,NPREF
      REAL*8 PREF(NPREF),FP(NPREF),DREF(NDREF)

      REAL*8, ALLOCATABLE:: BF(:), BFP(:), SDP(:),
     &                      SP(:), BFM(:), SDM(:),
     &                      SM(:)

C Set up the matrices BFP(tu,xy) and BFM(tu,xy)
C Formulae used:
C    BF(tu,xy)= 2*(Ftxuy - EASUM*Gtxuy)
C    BFP(tu,xy)=BF(tu,xy)+BF(tu,yx)
C    BFM(tu,xy)=BF(tu,xy)-BF(tu,yx)


C Loop over superindex symmetry.
      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,8)
        IF(NINP.EQ.0) GOTO 1000
        NAS=NTU(ISYM)
        NBF=(NAS*(NAS+1))/2
        IF(NBF.GT.0) THEN
          CALL mma_allocate(BF,NBF,LABEL='BF')
        END IF
        DO ITU=1,NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          DO IXY=1,ITU
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            IBADR=(ITU*(ITU-1))/2+IXY
            ITX=ITABS+NASHT*(IXABS-1)
            IUY=IUABS+NASHT*(IYABS-1)
            IP1=MAX(ITX,IUY)
            IP2=MIN(ITX,IUY)
            IP=(IP1*(IP1-1))/2+IP2
            BF(IBADR)=4.0D0*(FP(IP)-EASUM*PREF(IP))
          END DO
        END DO
        NASP=NTGEU(ISYM)
        NBFP=(NASP*(NASP+1))/2
        IF(NBFP.GT.0) THEN
          CALL mma_allocate(BFP,NBFP,Label='BFP')
CGG.Nov03  Load in SDP the diagonal elements of SFP matrix:
          NSP=(NASP*(NASP+1))/2
          CALL mma_allocate(SP,NSP,Label='SP')
          CALL mma_allocate(SDP,NASP,Label='SDP')
          IDSP=IDSMAT(ISYM,8)
          CALL DDAFILE(LUSBT,2,SP,NSP,IDSP)
          IDIAG=0
          DO I=1,NASP
            IDIAG=IDIAG+I
            SDP(I)=SP(IDIAG)
          END DO
          CALL mma_deallocate(SP)
CGG End
        END IF
        NASM=NTGTU(ISYM)
        NBFM=(NASM*(NASM+1))/2
        IF(NBFM.GT.0) THEN
          CALL mma_allocate(BFM,NBFM,Label='BFM')
CGG.Nov03  Load in SDM the diagonal elements of SFM matrix:
          NSM=(NASM*(NASM+1))/2
          CALL mma_allocate(SM,NSM,Label='SM')
          CALL mma_allocate(SDM,NASM,Label='SDM')
          IDSM=IDSMAT(ISYM,9)
          CALL DDAFILE(LUSBT,2,SM,NSM,IDSM)
          IDIAG=0
          DO I=1,NASM
            IDIAG=IDIAG+I
            SDM(I)=SM(IDIAG)
          END DO
          CALL mma_deallocate(SM)
CGG End
        END IF
        INSM=1
        DO ITGEU=1,NASP
          ITGEUABS=ITGEU+NTGEUES(ISYM)
          ITABS=MTGEU(1,ITGEUABS)
          IUABS=MTGEU(2,ITGEUABS)
          ITU=KTU(ITABS,IUABS)-NTUES(ISYM)
          DO IXGEY=1,ITGEU
            IXGEYABS=IXGEY+NTGEUES(ISYM)
            IXABS=MTGEU(1,IXGEYABS)
            IYABS=MTGEU(2,IXGEYABS)
            IXY=KTU(IXABS,IYABS)-NTUES(ISYM)
            IYX=KTU(IYABS,IXABS)-NTUES(ISYM)
            IF(ITU.GE.IXY) THEN
              IBADR=(ITU*(ITU-1))/2+IXY
            ELSE
              IBADR=(IXY*(IXY-1))/2+ITU
            END IF
            BTUXY=BF(IBADR)
            IF(ITU.GE.IYX) THEN
              IBADR=(ITU*(ITU-1))/2+IYX
            ELSE
              IBADR=(IYX*(IYX-1))/2+ITU
            END IF
            BTUYX=BF(IBADR)
            IBPADR=(ITGEU*(ITGEU-1))/2+IXGEY
            BFP(IBPADR)=BTUXY+BTUYX
CGG.Nov03
            IF (ITGEU.eq.IXGEY) THEN
              IDT=(ITABS*(ITABS+1))/2
              IDU=(IUABS*(IUABS+1))/2
              BFP(IBPADR)=BFP(IBPADR)+ipea_shift*0.5d0*
     &                    (4.0d0-DREF(IDT)-DREF(IDU))*SDP(ITGEU)
            ENDIF
CGG End
            IF(ITABS.EQ.IUABS) GOTO 200
            IF(IXABS.EQ.IYABS) GOTO 200
            ITGTU=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
            IXGTY=KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
            IBMADR=(ITGTU*(ITGTU-1))/2+IXGTY
            BFM(IBMADR)=BTUXY-BTUYX
CGG.Nov03
            IF (ITGEU.eq.IXGEY) THEN
              IDT=(ITABS*(ITABS+1))/2
              IDU=(IUABS*(IUABS+1))/2
              BFM(IBMADR)=BFM(IBMADR)+ipea_shift*0.5d0*
     &                    (4.0d0-DREF(IDT)-DREF(IDU))*SDM(INSM)
              INSM=INSM+1
            ENDIF

 200        CONTINUE
          END DO
        END DO
        IF(NBF.GT.0) CALL mma_deallocate(BF)

C Write to disk
        IF(NBFP.GT.0.and.NINDEP(ISYM,8).GT.0) THEN
          IDISK=IDBMAT(ISYM,8)
          CALL DDAFILE(LUSBT,1,BFP,NBFP,IDISK)
          CALL mma_deallocate(BFP)
CGG.Nov03 DisAlloc SDP
          CALL mma_deallocate(SDP)
CGG End
        END IF
        IF(NBFM.GT.0) THEN
         IF(NINDEP(ISYM,9).GT.0) THEN
          IDISK=IDBMAT(ISYM,9)
          CALL DDAFILE(LUSBT,1,BFM,NBFM,IDISK)
         END IF
         CALL mma_deallocate(BFM)
CGG.Nov03 DisAlloc SDM
         CALL mma_deallocate(SDM)
CGG End
        END IF
 1000 CONTINUE
      END DO

      END SUBROUTINE MKBF

      SUBROUTINE MKBG(DREF,NDREF,FD)
      use caspt2_global, only:ipea_shift
      use caspt2_global, only:LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)

#include "caspt2.fh"

      INTEGER NDREF
      REAL*8 DREF(NDREF),FD(NDREF)

      REAL*8, ALLOCATABLE:: BG(:), S(:), SD(:)

C     Set up the matrix BG(t,x)
C     Formula used:
C     BG(t,x)= Ftx -EASUM*Dtx


      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,10)
        IF(NINP.EQ.0) GOTO 1000
        NINM=NINDEP(ISYM,11)
        NAS=NASH(ISYM)
        NBG=(NAS*(NAS+1))/2
        IF(NBG.GT.0) THEN
          CALL mma_Allocate(BG,NBG,LABEL='BG')
CGG.Nov03  Load in SD the diagonal elements of SG matrix:
          NS=(NAS*(NAS+1))/2
          CALL mma_allocate(S,NS,Label='S')
          CALL mma_allocate(SD,NAS,Label='SD')
          IDS=IDSMAT(ISYM,10)
          CALL DDAFILE(LUSBT,2,S,NS,IDS)
          IDIAG=0
          DO I=1,NAS
            IDIAG=IDIAG+I
            SD(I)=S(IDIAG)
          END DO
          CALL mma_deallocate(S)
        ENDIF
CGG End
        DO IT=1,NAS
          ITABS=IT+NAES(ISYM)
          DO IX=1,IT
            IXABS=IX+NAES(ISYM)
            IBG=(IT*(IT-1))/2+IX
            ID=(ITABS*(ITABS-1))/2+IXABS
CGG.Nov03
            VALUE=FD(ID)-EASUM*DREF(ID)
            IF (IT.eq.IX) THEN
              IDT=(ITABS*(ITABS+1))/2
              VALUE = VALUE +
     &                ipea_shift*0.5d0*(2.0d0-DREF(IDT))*SD(IT)
            ENDIF
            BG(IBG)=VALUE
C           BG(BG)=FD(ID)-EASUM*DREF(ID)
CGG End
          END DO
        END DO

C Write to disk, and save size and address.
        IF(NBG.GT.0) THEN
         IF(NINDEP(ISYM,10).GT.0) THEN
          IDISK=IDBMAT(ISYM,10)
          CALL DDAFILE(LUSBT,1,BG,NBG,IDISK)
         END IF
         IF(NINM.GT.0.and.NINDEP(ISYM,11).GT.0) THEN
           IDISK=IDBMAT(ISYM,11)
           CALL DDAFILE(LUSBT,1,BG,NBG,IDISK)
         END IF
CGG.Nov03 DisAlloc SD
         CALL mma_deallocate(SD)
CGG End
         CALL mma_deallocate(BG)
        END IF
 1000 CONTINUE
      END DO

      END SUBROUTINE MKBG
