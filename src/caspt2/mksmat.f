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
      SUBROUTINE MKSMAT()
      use definitions, only: iwp, wp, u6, byte
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug, verbose
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_global, only: DREF, PREF
      use caspt2_global, only: LUSOLV, LUSBT
      use EQSOLV
      use caspt2_module
      use pt2_guga
      IMPLICIT REAL*8 (A-H,O-Z)
C     Set up S matrices for cases 1..13.
      real(kind=wp) DUM(1)
      INTEGER(kind=byte), ALLOCATABLE :: idxG3(:,:)

      real(kind=wp), ALLOCATABLE:: G3(:)

      IF(IPRGLB.GE.VERBOSE) THEN
        WRITE(u6,*)
        WRITE(u6,*)' Construct S matrices'
      END IF

      IF(NASHT.GT.0) THEN
CSVC: print header for debug info
        IF(IPRGLB.GE.DEBUG) THEN
          WRITE(u6,'("DEBUG> ",A)') 'CASE SYM S-MATRIX NORM'
          WRITE(u6,'("DEBUG> ",A)') '==== === ============='
        END IF
C For the cases A and C, begin by reading in the local storage
C  part of the three-electron density matrix G3:
        CALL mma_allocate(G3,NG3,Label='G3')
        CALL PT2_GET(NG3,'GAMMA3',G3)
        CALL mma_allocate(idxG3,6,NG3,label='idxG3')
        iLUID=0
        CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

        CALL MKSA(DREF,SIZE(DREF),PREF,SIZE(PREF),NG3,G3,idxG3)
        CALL MKSC(DREF,SIZE(DREF),PREF,SIZE(PREF),NG3,G3,idxG3)

        CALL mma_deallocate(G3)
        CALL mma_deallocate(idxG3)

C-SVC20100902: For the remaining cases that do not need G3, use replicate arrays
        CALL MKSB(DREF,SIZE(DREF),PREF,SIZE(PREF))
        CALL MKSD(DREF,SIZE(DREF),PREF,SIZE(PREF))
        CALL MKSE(DREF,SIZE(DREF))
        CALL MKSF(PREF,SIZE(PREF))
        CALL MKSG(DREF,SIZE(DREF))
      END IF

C For completeness, even case H has formally S and B
C matrices. This costs nothing, and saves conditional
C looping, etc in the rest  of the routines.
      DUM(1)=1.0D00
      DO ISYM=1,NSYM
        DO ICASE=12,13
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.GT.0) THEN
            IDISK=IDSMAT(ISYM,ICASE)
            CALL DDAFILE(LUSBT,1,DUM,1,IDISK)
          END IF
        END DO
      END DO


      RETURN
      END SUBROUTINE MKSMAT

********************************************************************************
* Case A (ICASE=1)
********************************************************************************
      SUBROUTINE MKSA(DREF,NDREF,PREF,NPREF,NG3,G3,idxG3)
      use definitions, only: iwp, wp, u6, Byte
      USE SUPERINDEX
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use EQSOLV
      use fake_GA, only: GA_Arrays
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp) NDREF,NPREF, NG3
      real(kind=wp) DREF(NDREF),PREF(NPREF),G3(NG3)
      INTEGER(kind=Byte) idxG3(6,NG3)
#ifdef _MOLCAS_MPP_
      real(kind=wp) Dummy(1)
#endif

      ICASE=1
C LONG loop over superindex symmetry.
      DO ISYM=1,NSYM
        NIN=NINDEP(ISYM,ICASE)
        IF(NIN.EQ.0) CYCLE
        NAS=NTUV(ISYM)
        NSA=(NAS*(NAS+1))/2
        IF(NSA.LE.0) CYCLE
C Set up the matrix SA(tuv,xyz) defined by the expression
C <ituv|kxyz> = dik SA(tuv,xyz)
C Formula used:
C    SA(tuv,xyz) =  -Gvuxtyz -dyu Gvzxt - dyt Gvuxz -
C         - dxu Gvtyz - dxu dyt Gvz +2 dtx Gvuyz + 2 dtx dyu Gvz

        CALL PSBMAT_GETMEM('SA',lg_SA,NAS)

        ! fill in the 3-el parts
#ifdef _MOLCAS_MPP_
        IF (IS_REAL_PAR()) THEN
          MYRANK = GA_NODEID()
          CALL GA_DISTRIBUTION (LG_SA,MYRANK,ILO,IHI,JLO,JHI)
          IF (JLO.NE.0 .AND. (JHI-JLO+1).NE.NAS) THEN
            WRITE(u6,*) 'MKSA: MISMATCH IN RANGE OF THE SUPERINDICES'
            CALL ABEND()
          END IF
          IF (ILO.GT.0 .AND. JLO.GT.0) THEN
            CALL GA_ACCESS (LG_SA,ILO,IHI,JLO,JHI,MA,LDA)
            CALL MKSA_G3_MPP(ISYM,DBL_MB(MA),ILO,IHI,JLO,JHI,LDA,
     &                       NG3,G3,IDXG3)
            CALL MKSA_DP(DREF,NDREF,PREF,NPREF,
     &                   ISYM,DBL_MB(MA),ILO,IHI,JLO,JHI,LDA)
            CALL GA_RELEASE_UPDATE (LG_SA,ILO,IHI,JLO,JHI)
          ELSE
            CALL MKSA_G3_MPP(ISYM,DUMMY,ILO,IHI,JLO,JHI,LDA,
     &                       NG3,G3,IDXG3)
          END IF
        ELSE
#endif
          CALL MKSA_G3(ISYM,GA_Arrays(lg_SA)%A(:),NG3,G3,IDXG3)
          CALL MKSA_DP(DREF,NDREF,PREF,NPREF,
     &                 ISYM,GA_Arrays(lg_SA)%A(:),1,NAS,1,NAS,0)
#ifdef _MOLCAS_MPP_
        END IF
#endif

        CALL PSBMAT_WRITE('S',iCase,iSYM,lg_SA,NAS)

        IF(IPRGLB.GE.DEBUG) THEN
          DSA=PSBMAT_FPRINT(lg_SA,NAS)
          WRITE(u6,'("DEBUG> ",A4,1X,I3,1X,ES21.14)') 'A', ISYM, DSA
        END IF

        CALL PSBMAT_FREEMEM(lg_SA)
      END DO

      END SUBROUTINE MKSA

      SUBROUTINE MKSA_G3(ISYM,SA,NG3,G3,idxG3)
      use definitions, only: iwp, wp, Byte
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)

      real(kind=wp) SA(*)
      real(kind=wp) G3(NG3)
      INTEGER(kind=Byte) idxG3(6,NG3)

      integer(kind=iwp) iG3

C-SVC20100831: determine indices in SA where a certain G3 value will end up
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
        G3VAL=-G3(iG3)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - G(tuvxyz) -> SA(xut,vyz)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iV,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
C  - G(vxtuyz) -> SA(uxv,tyz)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iT,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(yzvxtu) -> SA(xzy,vtu)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iV,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(tuyzvx) -> SA(zut,yvx)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iY,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
 200   CONTINUE
C  - G(yztuvx) -> SA(uzy,tvx)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iT,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(vxyztu) -> SA(zxv,ytu)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iY,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
 300   CONTINUE
        if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 500
        if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 500
        if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 500
        if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 500
C  - G(utxvzy) -> SA(vtu,xzy)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iX,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 500
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
C  - G(xvutzy) -> SA(tvx,uzy)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iU,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(zyxvut) -> SA(vyz,xut)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iX,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(utzyxv) -> SA(ytu,zxv)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
 400   CONTINUE
C  - G(zyutxv) -> SA(tyz,uxv)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iU,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(xvzyut) -> SA(yvx,zut)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SA(ISADR)=G3VAL
          END IF
        ENDIF
 500   CONTINUE
      END DO

      END SUBROUTINE MKSA_G3

#ifdef _MOLCAS_MPP_
      SUBROUTINE MKSA_G3_MPP(ISYM,SA,iLo,iHi,jLo,jHi,LDA,
     &                       NG3,G3,idxG3)
      use definitions, only: iwp, wp, MPIInt, RtoB, Byte
      USE MPI
      USE SUPERINDEX
      use stdalloc, only: mma_MaxDBLE
      use EQSOLV
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)

#include "global.fh"
#include "mafdecls.fh"

      real(kind=wp) SA(LDA,*)
      real(kind=wp) G3(NG3)
      INTEGER(kind=Byte) idxG3(6,NG3)

      integer(kind=MPIInt), ALLOCATABLE :: SCOUNTS(:), RCOUNTS(:)
      integer(kind=MPIInt), ALLOCATABLE :: SCOUNTS2(:), RCOUNTS2(:)
      integer(kind=MPIInt), ALLOCATABLE :: SDISPLS(:), RDISPLS(:)
      integer(kind=MPIInt), ALLOCATABLE :: SDISPLS2(:), RDISPLS2(:)

      integer(kind=MPIInt), ALLOCATABLE :: SENDIDX(:), RECVIDX(:)
      real(kind=wp),    ALLOCATABLE :: SENDVAL(:), RECVVAL(:)

      integer(kind=MPIInt), PARAMETER :: ONE4=1, TWO4=2
      integer(kind=MPIInt) :: IERROR4

      integer(kind=iwp), ALLOCATABLE :: IBUF(:)

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
      ! we need two real and four integer values per element
      iscal = (MPIInt*4 + wp*2)/RtoB ! RtoB from module Definitions
      !MAXBUF=MIN(NINT(0.95D0*MAXMEM)/4,2000000000/8)
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

          ! - G(tuvxyz) -> SA(xut,vyz)
          jSYM=MUL(iSX,MUL(iSU,iST))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
          ! - G(vxtuyz) -> SA(uxv,tyz)
          jSYM=MUL(iSU,MUL(iSX,iSV))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(yzvxtu) -> SA(xzy,vtu)
          jSYM=MUL(iSX,MUL(iSZ,iSY))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(tuyzvx) -> SA(zut,yvx)
          jSYM=MUL(iSZ,MUL(iSU,iST))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 200      CONTINUE
          ! - G(yztuvx) -> SA(uzy,tvx)
          jSYM=MUL(iSU,MUL(iSZ,iSY))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(vxyztu) -> SA(zxv,ytu)
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
          ! - G(utxvzy) -> SA(vtu,xzy)
          jSYM=MUL(iSV,MUL(iST,iSU))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) CYCLE
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
          ! - G(xvutzy) -> SA(tvx,uzy)
          jSYM=MUL(iST,MUL(iSV,iSX))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(zyxvut) -> SA(vyz,xut)
          jSYM=MUL(iSV,MUL(iSY,iSZ))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(utzyxv) -> SA(ytu,zxv)
          jSYM=MUL(iSY,MUL(iST,iSU))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 400      CONTINUE
          ! - G(zyutxv) -> SA(tyz,uxv)
          jSYM=MUL(iST,MUL(iSY,iSZ))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(xvzyut) -> SA(yvx,zut)
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
          SDISPLS(I)=INT(IOFFSET,kind=MPIInt)
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
          G3VAL=-G3(iG3)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - G(tuvxyz) -> SA(xut,vyz)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 301
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 201
C  - G(vxtuyz) -> SA(uxv,tyz)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
C  - G(yzvxtu) -> SA(xzy,vtu)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iV,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
C  - G(tuyzvx) -> SA(zut,yvx)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iY,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
 201      CONTINUE
C  - G(yztuvx) -> SA(uzy,tvx)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iT,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
C  - G(vxyztu) -> SA(zxv,ytu)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iY,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
 301      CONTINUE
          if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) CYCLE
          if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) CYCLE
          if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) CYCLE
          if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) CYCLE
C  - G(utxvzy) -> SA(vtu,xzy)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) CYCLE
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 401
C  - G(xvutzy) -> SA(tvx,uzy)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
C  - G(zyxvut) -> SA(vyz,xut)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iX,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
C  - G(utzyxv) -> SA(ytu,zxv)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
 401      CONTINUE
C  - G(zyutxv) -> SA(tyz,uxv)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iU,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
C  - G(xvzyut) -> SA(yvx,zut)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
        END DO

        ! Now we need to determine the receive counts.
        CALL MPI_ALLTOALL(SCOUNTS, ONE4, MPI_INTEGER,
     &                    RCOUNTS, ONE4, MPI_INTEGER,
     &                    MPI_COMM_WORLD, IERROR4)

        IOFFSET=0
        DO I=1,NPROCS
          RDISPLS(I)=INT(IOFFSET,kind=MPIInt)
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
          SA(ISUP-ILO+1,JSUP-JLO+1)=RECVVAL(I)
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

      PURE FUNCTION IPROW(IROW,NQOT,NREM)
      use definitions, only: iwp
      implicit None
      integer(kind=iwp) IPROW
      integer(kind=iwp), INTENT(IN) :: IROW, NQOT, NREM
      integer(kind=iwp) :: TMP
      TMP=IROW-NREM*(NQOT+1)
      IF (TMP.GT.0) THEN
        IPROW=(TMP-1)/NQOT+NREM+1
      ELSE
        IPROW=(IROW-1)/(NQOT+1)+1
      END IF
      END FUNCTION IPROW

      END SUBROUTINE MKSA_G3_MPP
#endif

      SUBROUTINE MKSA_DP (DREF,NDREF,PREF,NPREF,
     &                    iSYM,SA,iLo,iHi,jLo,jHi,LDA)
C In parallel, this subroutine is called on a local chunk of memory
C and LDA is set. In serial, the whole array is passed but then the
C storage uses a triangular scheme, and the LDA passed is zero.
      use definitions, only: iwp, wp
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)
      integer(kind=iwp) NDREF,NPREF,iSYM,iLo,iHi,jLo,jHi,LDA
      real(kind=wp) DREF(NDREF),PREF(NPREF)
      real(kind=wp) SA(*)

      ISADR=0
C-SVC20100831: fill in the G2 and G1 corrections for SA
      DO 100 IXYZ=jLo,jHi
        IXYZABS=IXYZ+NTUVES(ISYM)
        IXABS=MTUV(1,IXYZABS)
        IYABS=MTUV(2,IXYZABS)
        IZABS=MTUV(3,IXYZABS)
        DO 101 ITUV=iLo,iHi
          ITUVABS=ITUV+NTUVES(ISYM)
          ITABS=MTUV(1,ITUVABS)
          IUABS=MTUV(2,ITUVABS)
          IVABS=MTUV(3,ITUVABS)
C Add  2 dtx Gvuyz + 2 dtx dyu Gvz
          IF (LDA.NE.0) THEN
            VALUE=SA(1+(iTUV-iLo)+LDA*(iXYZ-jLo))
          ELSE
            IF (IXYZ.LE.ITUV) THEN
              ISADR=(ITUV*(ITUV-1))/2+IXYZ
              VALUE=SA(ISADR)
            ELSE
              GOTO 101
            ENDIF
          END IF
          IF(ITABS.EQ.IXABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,IYZ)
            IP2=MIN(IVU,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+4.0D0*PREF(IP)
            IF(IYABS.EQ.IUABS)THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE+2.0D0*DREF(ID)
            END IF
          END IF
C Add  -dxu Gvtyz -dxu dyt Gvz
          IF(IXABS.EQ.IUABS) THEN
            IVT=IVABS+NASHT*(ITABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVT,IYZ)
            IP2=MIN(IVT,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE - 2.0D0*PREF(IP)
            IF(IYABS.EQ.ITABS)THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE - DREF(ID)
            END IF
          END IF
C Add  -dyt Gvuxz
          IF(IYABS.EQ.ITABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            IXZ=IXABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,IXZ)
            IP2=MIN(IVU,IXZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE - 2.0D0*PREF(IP)
          END IF
C Add -dyu Gvzxt
          IF(IYABS.EQ.IUABS) THEN
            IVZ=IVABS+NASHT*(IZABS-1)
            IXT=IXABS+NASHT*(ITABS-1)
            IP1=MAX(IVZ,IXT)
            IP2=MIN(IVZ,IXT)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE - 2.0D0*PREF(IP)
          END IF
          IF (LDA.NE.0) THEN
            SA(1+(iTUV-iLo)+LDA*(iXYZ-jLo))=VALUE
          ELSE
            SA(ISADR)=VALUE
          END IF
 101    CONTINUE
 100  CONTINUE

      END SUBROUTINE MKSA_DP

********************************************************************************
* Case C (ICASE=4)
********************************************************************************
      SUBROUTINE MKSC(DREF,NDREF,PREF,NPREF,NG3,G3,idxG3)
      use definitions, only: iwp, wp, u6, Byte
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug
      USE SUPERINDEX
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use EQSOLV
      use fake_GA, only: GA_Arrays
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp) NDREF,NPREF, NG3
      real(kind=wp) DREF(NDREF),PREF(NPREF),G3(NG3)
      INTEGER(kind=Byte) idxG3(6,NG3)
#ifdef _MOLCAS_MPP_
      real(kind=wp) Dummy(1)
#endif

      ICASE=4
C LONG loop over superindex symmetry.
      DO ISYM=1,NSYM
        NIN=NINDEP(ISYM,ICASE)
        IF(NIN.EQ.0) CYCLE
        NAS=NTUV(ISYM)
        NSC=(NAS*(NAS+1))/2
        IF(NSC.LE.0) CYCLE

C Set up the matrix SC(tuv,xyz) defined by the expression
C <atuv|cxyz> = dac SC(tuv,xyz)
C Formula used:
C    SC(tuv,xyz)
C    = Gvutxyz +dyu Gvztx + dyx Gvutz + dtu Gvxyz + dtu dyx Gvz

        CALL PSBMAT_GETMEM('SC',lg_SC,NAS)

        ! fill in the 3-el parts
#ifdef _MOLCAS_MPP_
        IF (IS_REAL_PAR()) THEN
          MYRANK = GA_NODEID()
          CALL GA_DISTRIBUTION (LG_SC,MYRANK,ILO,IHI,JLO,JHI)
          IF (JLO.NE.0 .AND. (JHI-JLO+1).NE.NAS) THEN
            WRITE(u6,*) 'MKSC: MISMATCH IN RANGE OF THE SUPERINDICES'
            CALL ABEND()
          END IF
          IF (ILO.GT.0 .AND. JLO.GT.0) THEN
            CALL GA_ACCESS (LG_SC,ILO,IHI,JLO,JHI,MC,LDC)
            CALL MKSC_G3_MPP(ISYM,DBL_MB(MC),ILO,IHI,JLO,JHI,LDC,
     &                       NG3,G3,IDXG3)
            CALL MKSC_DP(DREF,NDREF,PREF,NPREF,
     &                   ISYM,DBL_MB(MC),ILO,IHI,JLO,JHI,LDC)
            CALL GA_RELEASE_UPDATE (LG_SC,ILO,IHI,JLO,JHI)
          ELSE
            CALL MKSC_G3_MPP(ISYM,DUMMY,ILO,IHI,JLO,JHI,LDC,
     &                       NG3,G3,IDXG3)
          END IF
        ELSE
#endif
          CALL MKSC_G3(ISYM,GA_Arrays(lg_SC)%A(:),NG3,G3,IDXG3)
          CALL MKSC_DP(DREF,NDREF,PREF,NPREF,
     &                 ISYM,GA_Arrays(lg_SC)%A(:),1,NAS,1,NAS,0)
#ifdef _MOLCAS_MPP_
        END IF
#endif

        CALL PSBMAT_WRITE('S',iCase,iSYM,lg_SC,NAS)

        IF(IPRGLB.GE.DEBUG) THEN
          DSC=PSBMAT_FPRINT(lg_SC,NAS)
          WRITE(u6,'("DEBUG> ",A4,1X,I3,1X,ES21.14)') 'C', ISYM, DSC
        END IF

        CALL PSBMAT_FREEMEM(lg_SC)
      END DO

      END SUBROUTINE MKSC

      SUBROUTINE MKSC_G3(ISYM,SC,NG3,G3,idxG3)
      use definitions, only: iwp, wp, Byte
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)

      integer(kind=iwp) ISYM,NG3
      real(kind=wp) SC(*)
      real(kind=wp) G3(NG3)
      INTEGER(kind=Byte) idxG3(6,NG3)

C-SVC20100831: determine indices in SC where a certain G3 value will end up
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
        G3VAL=G3(iG3)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - G(tuvxyz) -> SC(vut,xyz)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iX,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
C  - G(vxtuyz) -> SC(txv,uyz)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iU,iY,iZ)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(yzvxtu) -> SC(vzy,xtu)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iX,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(tuyzvx) -> SC(yut,zvx)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
 200   CONTINUE
C  - G(yztuvx) -> SC(tzy,uvx)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iU,iV,iX)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(vxyztu) -> SC(yxv,ztu)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iT,iU)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
 300   CONTINUE
        if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 500
        if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 500
        if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 500
        if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 500
C  - G(utxvzy) -> SC(xtu,vzy)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iV,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 500
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
C  - G(xvutzy) -> SC(uvx,tzy)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iT,iZ,iY)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(zyxvut) -> SC(xyz,vut)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iV,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(utzyxv) -> SC(ztu,yxv)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iY,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
 400   CONTINUE
C  - G(zyutxv) -> SC(uyz,txv)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iT,iX,iV)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
C  - G(xvzyut) -> SC(zvx,yut)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iY,iU,iT)-nTUVES(jSYM)
          IF (JSUP.LE.ISUP) THEN
            ISADR=(ISUP*(ISUP-1))/2+JSUP
            SC(ISADR)=G3VAL
          END IF
        ENDIF
 500   CONTINUE
      END DO

      END SUBROUTINE MKSC_G3

#ifdef _MOLCAS_MPP_
      SUBROUTINE MKSC_G3_MPP(ISYM,SC,iLo,iHi,jLo,jHi,LDC,
     &                       NG3,G3,idxG3)
      use definitions, only: iwp, wp, Byte
      USE MPI
      USE SUPERINDEX
      use stdalloc, only: mma_MaxDBLE
      use EQSOLV
      use definitions, only: MPIInt,RtoB,wp
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)

#include "global.fh"
#include "mafdecls.fh"

      real(kind=wp) SC(LDC,*)
      real(kind=wp) G3(NG3)
      INTEGER(kind=Byte) idxG3(6,NG3)

      integer(kind=MPIInt), ALLOCATABLE :: SCOUNTS(:), RCOUNTS(:)
      integer(kind=MPIInt), ALLOCATABLE :: SCOUNTS2(:), RCOUNTS2(:)
      integer(kind=MPIInt), ALLOCATABLE :: SDISPLS(:), RDISPLS(:)
      integer(kind=MPIInt), ALLOCATABLE :: SDISPLS2(:), RDISPLS2(:)

      integer(kind=MPIInt), ALLOCATABLE :: SENDIDX(:), RECVIDX(:)
      real(kind=wp),    ALLOCATABLE :: SENDVAL(:), RECVVAL(:)

      integer(kind=MPIInt), PARAMETER :: ONE4=1, TWO4=2
      integer(kind=MPIInt) :: IERROR4
      INTEGER, PARAMETER :: I4=KIND(ONE4)

      INTEGER(kind=iwp), ALLOCATABLE :: IBUF(:)

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
      ! 12 corresponds to the number of if (jSym.Eq.iSym) branch
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
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - G(tuvxyz) -> SC(vut,xyz)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
C  - G(vxtuyz) -> SC(txv,uyz)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - G(yzvxtu) -> SC(vzy,xtu)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - G(tuyzvx) -> SC(yut,zvx)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 200   CONTINUE
C  - G(yztuvx) -> SC(tzy,uvx)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - G(vxyztu) -> SC(yxv,ztu)
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
C  - G(utxvzy) -> SC(xtu,vzy)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 500
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
C  - G(xvutzy) -> SC(uvx,tzy)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - G(zyxvut) -> SC(xyz,vut)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - G(utzyxv) -> SC(ztu,yxv)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
 400   CONTINUE
C  - G(zyutxv) -> SC(uyz,txv)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
C  - G(xvzyut) -> SC(zvx,yut)
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
          G3VAL=G3(iG3)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - G(tuvxyz) -> SC(vut,xyz)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iX,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 301
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 201
C  - G(vxtuyz) -> SC(txv,uyz)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iU,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - G(yzvxtu) -> SC(vzy,xtu)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iX,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - G(tuyzvx) -> SC(yut,zvx)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
 201   CONTINUE
C  - G(yztuvx) -> SC(tzy,uvx)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iU,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - G(vxyztu) -> SC(yxv,ztu)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
 301   CONTINUE
          if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 501
          if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 501
          if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 501
          if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 501
C  - G(utxvzy) -> SC(xtu,vzy)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iV,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 501
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 401
C  - G(xvutzy) -> SC(uvx,tzy)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iT,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - G(zyxvut) -> SC(xyz,vut)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iV,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - G(utzyxv) -> SC(ztu,yxv)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iY,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
 401   CONTINUE
C  - G(zyutxv) -> SC(uyz,txv)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iT,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,I4)
            SENDIDX(2*IBUF(IP))=INT(ICOL,I4)
          ENDIF
C  - G(xvzyut) -> SC(zvx,yut)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iY,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
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
          SC(ISUP-ILO+1,JSUP-JLO+1)=RECVVAL(I)
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

      PURE FUNCTION IPROW(IROW,NQOT,NREM)
      use definitions, only: iwp
      INTEGER(kind=iwp) :: IPROW
      INTEGER(kind=iwp), INTENT(IN) :: IROW, NQOT, NREM
      INTEGER(kind=iwp) :: TMP
      TMP=IROW-NREM*(NQOT+1)
      IF (TMP.GT.0) THEN
        IPROW=(TMP-1)/NQOT+NREM+1
      ELSE
        IPROW=(IROW-1)/(NQOT+1)+1
      END IF
      END FUNCTION IPROW

      END SUBROUTINE MKSC_G3_MPP
#endif

      SUBROUTINE MKSC_DP (DREF,NDREF,PREF,NPREF,
     &                    iSYM,SC,iLo,iHi,jLo,jHi,LDC)
C In parallel, this subroutine is called on a local chunk of memory
C and LDC is set. In serial, the whole array is passed but then the
C storage uses a triangular scheme, and the LDC passed is zero.
      use definitions, only: iwp, wp
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)
      integer(kind=iwp) :: NDREF,NPREF,iSYM,iLo,iHi,jLo,jHi,LDC
      real(kind=wp) DREF(NDREF),PREF(NPREF)
      real(kind=wp) SC(*)

      ISADR=0
C-SVC20100831: fill in the G2 and G1 corrections for this SC block
      DO 100 IXYZ=jLo,jHi
        IXYZABS=IXYZ+NTUVES(ISYM)
        IXABS=MTUV(1,IXYZABS)
        IYABS=MTUV(2,IXYZABS)
        IZABS=MTUV(3,IXYZABS)
        DO 101 ITUV=iLo,iHi
          ITUVABS=ITUV+NTUVES(ISYM)
          ITABS=MTUV(1,ITUVABS)
          IUABS=MTUV(2,ITUVABS)
          IVABS=MTUV(3,ITUVABS)
          IF (LDC.NE.0) THEN
            VALUE=SC(1+iTUV-iLo+LDC*(iXYZ-jLo))
          ELSE
            IF (IXYZ.LE.ITUV) THEN
              ISADR=(ITUV*(ITUV-1))/2+IXYZ
              VALUE=SC(ISADR)
            ELSE
              GOTO 101
            ENDIF
          END IF
C Add  dyu Gvztx
          IF(IYABS.EQ.IUABS) THEN
            IVZ=IVABS+NASHT*(IZABS-1)
            ITX=ITABS+NASHT*(IXABS-1)
            IP1=MAX(IVZ,ITX)
            IP2=MIN(IVZ,ITX)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+2.0D0*PREF(IP)
          END IF
C Add  dyx Gvutz
          IF(IYABS.EQ.IXABS) THEN
            IVU=IVABS+NASHT*(IUABS-1)
            ITZ=ITABS+NASHT*(IZABS-1)
            IP1=MAX(IVU,ITZ)
            IP2=MIN(IVU,ITZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+2.0D0*PREF(IP)
          END IF
C Add  dtu Gvxyz + dtu dyx Gvz
          IF(ITABS.EQ.IUABS) THEN
            IVX=IVABS+NASHT*(IXABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IP1=MAX(IVX,IYZ)
            IP2=MIN(IVX,IYZ)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=VALUE+2.0D0*PREF(IP)
            IF(IYABS.EQ.IXABS) THEN
              ID1=MAX(IVABS,IZABS)
              ID2=MIN(IVABS,IZABS)
              VALUE=VALUE+DREF((ID1*(ID1-1))/2+ID2)
            END IF
          END IF
          IF (LDC.NE.0) THEN
            SC(1+iTUV-iLo+LDC*(iXYZ-jLo))=VALUE
          ELSE
            SC(ISADR)=VALUE
          END IF
 101    CONTINUE
 100  CONTINUE
      END SUBROUTINE MKSC_DP

********************************************************************************
* Case B (ICASE=2,3)
********************************************************************************
      SUBROUTINE MKSB(DREF,NDREF,PREF,NPREF)
      use definitions, only: iwp, wp
      USE SUPERINDEX
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER(kind=iwp) NDREF,NPREF
      REAL(kind=wp) DREF(NDREF),PREF(NPREF)

      REAL(kind=wp), ALLOCATABLE:: SB(:), SBP(:), SBM(:)
C Set up the matrices SBP(tu,xy) and SBM(tu,xy)
C Formulae used:
C    SB(tu,xy)=
C    = 4 Pxtyu -4dxt Dyu -4dyu Dxt +2dyt Dxu + 8 dxt dyu
C      -4dxu dyt + 2dxu Dyt
C    SBP(tu,xy)=SB(tu,xy)+SB(tu,yx)
C    SBM(tu,xy)=SB(tu,xy)-SB(tu,yx)


C Loop over superindex symmetry.
      DO 1000 ISYM=1,NSYM
        NINP=NINDEP(ISYM,2)
        IF(NINP.EQ.0) GOTO 1000
        NAS=NTU(ISYM)
        NSB=(NAS*(NAS+1))/2
        IF(NSB.GT.0) CALL mma_allocate(SB,NSB,Label='SB')
        DO 100 ITU=1,NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          DO 101 IXY=1,ITU
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            ISADR=(ITU*(ITU-1))/2+IXY
            IXT=IXABS+NASHT*(ITABS-1)
            IYU=IYABS+NASHT*(IUABS-1)
            IP1=MAX(IXT,IYU)
            IP2=MIN(IXT,IYU)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=4.0D0*PREF(IP)
C Add  -4 dxt Dyu + 8dxt dyu
            IF(IXABS.EQ.ITABS) THEN
              ID1=MAX(IYABS,IUABS)
              ID2=MIN(IYABS,IUABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE-4.0D0*DREF(ID)
              IF(IYABS.EQ.IUABS) VALUE=VALUE+8.0D00
            END IF
C Add  -4 dyu Dxt
            IF(IYABS.EQ.IUABS) THEN
              ID1=MAX(IXABS,ITABS)
              ID2=MIN(IXABS,ITABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE-4.0D0*DREF(ID)
            END IF
C Add  +2 dyt Dxu
            IF(IYABS.EQ.ITABS) THEN
              ID1=MAX(IXABS,IUABS)
              ID2=MIN(IXABS,IUABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE+2.0D0*DREF(ID)
            END IF
C Add  -4dxu dyt + 2dxu Dyt
            IF(IXABS.EQ.IUABS) THEN
              ID1=MAX(IYABS,ITABS)
              ID2=MIN(IYABS,ITABS)
              ID=(ID1*(ID1-1))/2+ID2
              VALUE=VALUE+2.0D0*DREF(ID)
              IF(IYABS.EQ.ITABS) VALUE=VALUE-4.0D00
            END IF
            ISADR=(ITU*(ITU-1))/2+IXY
            SB(ISADR)=VALUE
 101      CONTINUE
 100    CONTINUE
        NASP=NTGEU(ISYM)
        NSBP=(NASP*(NASP+1))/2
        IF(NSBP.GT.0) CALL mma_allocate(SBP,NSBP,Label='SBP')
        NASM=NTGTU(ISYM)
        NSBM=(NASM*(NASM+1))/2
        IF(NSBM.GT.0) CALL mma_allocate(SBM,NSBM,Label='SBM')
        DO 200 ITGEU=1,NASP
          ITGEUABS=ITGEU+NTGEUES(ISYM)
          ITABS=MTGEU(1,ITGEUABS)
          IUABS=MTGEU(2,ITGEUABS)
          ITU=KTU(ITABS,IUABS)-NTUES(ISYM)
          DO 201 IXGEY=1,ITGEU
            IXGEYABS=IXGEY+NTGEUES(ISYM)
            IXABS=MTGEU(1,IXGEYABS)
            IYABS=MTGEU(2,IXGEYABS)
            IXY=KTU(IXABS,IYABS)-NTUES(ISYM)
            IYX=KTU(IYABS,IXABS)-NTUES(ISYM)
            IF(ITU.GE.IXY) THEN
              ISADR=(ITU*(ITU-1))/2+IXY
            ELSE
              ISADR=(IXY*(IXY-1))/2+ITU
            END IF
            STUXY=SB(ISADR)
            IF(ITU.GE.IYX) THEN
              ISADR=(ITU*(ITU-1))/2+IYX
            ELSE
              ISADR=(IYX*(IYX-1))/2+ITU
            END IF
            STUYX=SB(ISADR)
            ISPADR=(ITGEU*(ITGEU-1))/2+IXGEY
            SBP(ISPADR)=STUXY+STUYX
            IF(ITABS.EQ.IUABS) GOTO 201
            IF(IXABS.EQ.IYABS) GOTO 201
            ITGTU=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
            IXGTY=KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
            ISMADR=(ITGTU*(ITGTU-1))/2+IXGTY
            SBM(ISMADR)=STUXY-STUYX
 201      CONTINUE
 200    CONTINUE
        IF(NSB.GT.0) CALL mma_deallocate(SB)

C Write to disk, and save size and address.
        IF(NSBP.GT.0) THEN
          IDISK=IDSMAT(ISYM,2)
          CALL DDAFILE(LUSBT,1,SBP,NSBP,IDISK)
          CALL mma_deallocate(SBP)
        END IF
        IF(NSBM.GT.0) THEN
          IF(NINDEP(ISYM,3).GT.0) THEN
            IDISK=IDSMAT(ISYM,3)
            CALL DDAFILE(LUSBT,1,SBM,NSBM,IDISK)
          END IF
          CALL mma_deallocate(SBM)
        END IF
 1000 CONTINUE


      RETURN
      END SUBROUTINE MKSB

      SUBROUTINE MKSD(DREF,NDREF,PREF,NPREF)
      use definitions, only: iwp, wp
      USE SUPERINDEX
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER(kind=iwp) NDREF,NPREF
      REAL(kind=wp) DREF(NDREF),PREF(NPREF)

      REAL(kind=wp), ALLOCATABLE:: SD(:)
C Set up the matrix SD(tuP,xyQ),P and Q are 1 or 2,
C Formulae used:
C    SD(tu1,xy1)=2*(Gutxy + dxt Duy)
C    SD(tu2,xy1)= -(Gutxy + dxt Duy)
C    SD(tu2,xy2)= -Gxtuy +2*dxt Duy


C Loop over superindex symmetry.
      DO 1000 ISYM=1,NSYM
        NIN=NINDEP(ISYM,5)
        IF(NIN.EQ.0) GOTO 1000
        NAS=NTU(ISYM)
        NSD=(2*NAS*(2*NAS+1))/2
        IF(NSD.GT.0) Call mma_allocate(SD,NSD,LABEL='SD')
        DO 100 ITU=1,NAS
        ITU2=ITU+NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          DO 101 IXY=1,ITU
            IXY2=IXY+NAS
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            IS11=(ITU*(ITU-1))/2+IXY
            IS21=(ITU2*(ITU2-1))/2+IXY
            IS12=(IXY2*(IXY2-1))/2+ITU
            IS22=(ITU2*(ITU2-1))/2+IXY2
            IUTP=IUABS+NASHT*(ITABS-1)
            IXYP=IXABS+NASHT*(IYABS-1)
            IP1=MAX(IUTP,IXYP)
            IP2=MIN(IUTP,IXYP)
            IP=(IP1*(IP1-1))/2+IP2
            GUTXY=2.0D0*PREF(IP)
            IXTP=IXABS+NASHT*(ITABS-1)
            IUYP=IUABS+NASHT*(IYABS-1)
            IP1=MAX(IXTP,IUYP)
            IP2=MIN(IXTP,IUYP)
            IP=(IP1*(IP1-1))/2+IP2
            GXTUY=2.0D0*PREF(IP)
            S11=2.0D0*GUTXY
            S22=-GXTUY
            IF(IXABS.EQ.ITABS) THEN
              ID1=MAX(IUABS,IYABS)
              ID2=MIN(IUABS,IYABS)
              ID=(ID1*(ID1-1))/2+ID2
              DUY=DREF(ID)
              S11=S11+2.0D0*DUY
              S22=S22+2.0D0*DUY
            END IF
C    SD(tu1,xy1)=2*(Gutxy + dtx Duy)
            SD(IS11)= S11
C    SD(tu2,xy1)= -(Gutxy + dtx Duy)
            SD(IS21)=-0.5D0*S11
            SD(IS12)=-0.5D0*S11
C    SD(tu2,xy2)= -Gxtuy +2*dtx Duy
            SD(IS22)= S22
 101      CONTINUE
 100    CONTINUE

C Write to disk
        IF(NSD.GT.0) THEN
         IF(NINDEP(ISYM,5).GT.0) THEN
          IDISK=IDSMAT(ISYM,5)
          CALL DDAFILE(LUSBT,1,SD,NSD,IDISK)
         END IF
         CALL mma_deallocate(SD)
        END IF
 1000 CONTINUE


      RETURN
      END SUBROUTINE MKSD

      SUBROUTINE MKSE(DREF,NDREF)
      use definitions, only: iwp, wp
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER(kind=iwp) NDREF
      REAL(kind=wp) DREF(NDREF)

      REAL(kind=wp), ALLOCATABLE:: SE(:)
C Set up the matrix SE(t,x)
C Formula used:
C    SE(t,x)=2*dtx - Dtx



      DO 1000 ISYM=1,NSYM
        NINP=NINDEP(ISYM,6)
        IF(NINP.EQ.0) GOTO 1000
        NINM=NINDEP(ISYM,7)
        NAS=NASH(ISYM)
        NSE=(NAS*(NAS+1))/2
        IF(NSE.GT.0) CALL mma_allocate(SE,NSE,Label='SE')
        DO 100 IT=1,NAS
          ITABS=IT+NAES(ISYM)
          DO 101 IX=1,IT
            IXABS=IX+NAES(ISYM)
            ISE=(IT*(IT-1))/2+IX
            ID=(ITABS*(ITABS-1))/2+IXABS
            IF(ITABS.EQ.IXABS) THEN
              SE(ISE)=2.0D00-DREF(ID)
            ELSE
              SE(ISE)=-DREF(ID)
            END IF
 101      CONTINUE
 100    CONTINUE

C Write to disk
        IF(NSE.GT.0.and.NINDEP(ISYM,6).GT.0) THEN
          IDISK=IDSMAT(ISYM,6)
          CALL DDAFILE(LUSBT,1,SE,NSE,IDISK)
          IF(NINM.GT.0.and.NINDEP(ISYM,7).GT.0) THEN
            IDISK=IDSMAT(ISYM,7)
            CALL DDAFILE(LUSBT,1,SE,NSE,IDISK)
          END IF
          CALL mma_deallocate(SE)
        END IF
 1000 CONTINUE

      END SUBROUTINE MKSE

      SUBROUTINE MKSF(PREF,NPREF)
      use definitions, only: iwp, wp
      USE SUPERINDEX
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER(kind=iwp) NPREF
      REAL(kind=wp) PREF(NPREF)

      REAL(kind=wp), ALLOCATABLE:: SF(:), SFP(:), SFM(:)
C Set up the matrices SFP(tu,xy) and SFM(tu,xy)
C Formulae used:
C    SF(tu,xy)= 4 Ptxuy
C    SFP(tu,xy)=SF(tu,xy)+SF(tu,yx)
C    SFM(tu,xy)=SF(tu,xy)-SF(tu,yx)



C Loop over superindex symmetry.
      DO 1000 ISYM=1,NSYM
        NINP=NINDEP(ISYM,8)
        IF(NINP.EQ.0) GOTO 1000
        NAS=NTU(ISYM)
        NSF=(NAS*(NAS+1))/2
        IF(NSF.GT.0) CALL mma_allocate(SF,NSF,Label='SF')
        DO 100 ITU=1,NAS
          ITUABS=ITU+NTUES(ISYM)
          ITABS=MTU(1,ITUABS)
          IUABS=MTU(2,ITUABS)
          DO 101 IXY=1,ITU
            IXYABS=IXY+NTUES(ISYM)
            IXABS=MTU(1,IXYABS)
            IYABS=MTU(2,IXYABS)
            ISADR=(ITU*(ITU-1))/2+IXY
            ITX=ITABS+NASHT*(IXABS-1)
            IUY=IUABS+NASHT*(IYABS-1)
            IP1=MAX(ITX,IUY)
            IP2=MIN(ITX,IUY)
            IP=(IP1*(IP1-1))/2+IP2
            VALUE=4.0D0*PREF(IP)
            SF(ISADR)=VALUE
 101      CONTINUE
 100    CONTINUE
        NASP=NTGEU(ISYM)
        NSFP=(NASP*(NASP+1))/2
        IF(NSFP.GT.0) CALL mma_allocate(SFP,NSFP,Label='SFP')
        NASM=NTGTU(ISYM)
        NSFM=(NASM*(NASM+1))/2
        IF(NSFM.GT.0) THEN
          CALL mma_allocate(SFM,NSFM,Label='SFM')
        END IF
        DO 200 ITGEU=1,NASP
          ITGEUABS=ITGEU+NTGEUES(ISYM)
          ITABS=MTGEU(1,ITGEUABS)
          IUABS=MTGEU(2,ITGEUABS)
          ITU=KTU(ITABS,IUABS)-NTUES(ISYM)
          DO 201 IXGEY=1,ITGEU
            IXGEYABS=IXGEY+NTGEUES(ISYM)
            IXABS=MTGEU(1,IXGEYABS)
            IYABS=MTGEU(2,IXGEYABS)
            IXY=KTU(IXABS,IYABS)-NTUES(ISYM)
            IYX=KTU(IYABS,IXABS)-NTUES(ISYM)
            IF(ITU.GE.IXY) THEN
              ISADR=(ITU*(ITU-1))/2+IXY
            ELSE
              ISADR=(IXY*(IXY-1))/2+ITU
            END IF
            STUXY=SF(ISADR)
            IF(ITU.GE.IYX) THEN
              ISADR=(ITU*(ITU-1))/2+IYX
            ELSE
              ISADR=(IYX*(IYX-1))/2+ITU
            END IF
            STUYX=SF(ISADR)
            ISPADR=(ITGEU*(ITGEU-1))/2+IXGEY
            SFP(ISPADR)=STUXY+STUYX
            IF(ITABS.EQ.IUABS) GOTO 201
            IF(IXABS.EQ.IYABS) GOTO 201
            ITGTU=KTGTU(ITABS,IUABS)-NTGTUES(ISYM)
            IXGTY=KTGTU(IXABS,IYABS)-NTGTUES(ISYM)
            ISMADR=(ITGTU*(ITGTU-1))/2+IXGTY
            SFM(ISMADR)=STUXY-STUYX
 201      CONTINUE
 200    CONTINUE
        IF(NSF.GT.0) CALL mma_deallocate(SF)

C Write to disk
        IF(NSFP.GT.0.and.NINDEP(ISYM,8).GT.0) THEN
          IDISK=IDSMAT(ISYM,8)
          CALL DDAFILE(LUSBT,1,SFP,NSFP,IDISK)
          CALL mma_deallocate(SFP)
        END IF
        IF(NSFM.GT.0) THEN
          IF(NINDEP(ISYM,9).GT.0) THEN
           IDISK=IDSMAT(ISYM,9)
           CALL DDAFILE(LUSBT,1,SFM,NSFM,IDISK)
          END IF
          CALL mma_deallocate(SFM)
        END IF
 1000 CONTINUE

      END SUBROUTINE MKSF

      SUBROUTINE MKSG(DREF,NDREF)
      use definitions, only: iwp, wp
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER(kind=iwp) NDREF
      REAL(kind=wp) DREF(NDREF)

      REAL(kind=wp), ALLOCATABLE:: SG(:)
C Set up the matrix SG(t,x)
C Formula used:
C    SG(t,x)= Dtx


      DO 1000 ISYM=1,NSYM
        NINP=NINDEP(ISYM,10)
        IF(NINP.EQ.0) GOTO 1000
        NINM=NINDEP(ISYM,11)
        NAS=NASH(ISYM)
        NSG=(NAS*(NAS+1))/2
        IF(NSG.GT.0) CALL mma_allocate(SG,NSG,Label='SG')
        DO 100 IT=1,NAS
          ITABS=IT+NAES(ISYM)
          DO 101 IX=1,IT
            IXABS=IX+NAES(ISYM)
            ISG=(IT*(IT-1))/2+IX
            ID=(ITABS*(ITABS-1))/2+IXABS
            SG(ISG)= DREF(ID)
 101      CONTINUE
 100    CONTINUE

C Write to disk
        IF(NSG.GT.0.and.NINDEP(ISYM,10).GT.0) THEN
          IDISK=IDSMAT(ISYM,10)
          CALL DDAFILE(LUSBT,1,SG,NSG,IDISK)
          IF(NINM.GT.0.and.NINDEP(ISYM,11).GT.0) THEN
            IDISK=IDSMAT(ISYM,11)
            CALL DDAFILE(LUSBT,1,SG,NSG,IDISK)
          END IF
          CALL mma_deallocate(SG)
        END IF
 1000 CONTINUE

      END SUBROUTINE MKSG
