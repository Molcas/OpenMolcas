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
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
#ifdef _MOLCAS_MPP_
      SUBROUTINE MKSA_G3_MPP(ISYM,SA,iLo,iHi,NAS,LDA,                   &
     &                       NG3,G3,idxG3)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp, MPIInt, Byte
      USE MPI
      USE SUPERINDEX, only: KTUV
      use stdalloc, only: mma_allocate,mma_deallocate,mma_MaxDBLE
      use caspt2_module, only: IASYM,NASHT,nTUVES
      IMPLICIT None

#include "global.fh"
#include "mafdecls.fh"

      integer(kind=iwp), intent(in):: ISYM,iLo,iHi,NAS,LDA,NG3
      real(kind=wp), intent(out):: SA(LDA,NAS)
      real(kind=wp), intent(in):: G3(NG3)
      INTEGER(kind=Byte), intent(in):: idxG3(6,NG3)

      integer(kind=MPIInt), ALLOCATABLE :: SCOUNTS(:), RCOUNTS(:)
      integer(kind=MPIInt), ALLOCATABLE :: SCOUNTS2(:), RCOUNTS2(:)
      integer(kind=MPIInt), ALLOCATABLE :: SDISPLS(:), RDISPLS(:)
      integer(kind=MPIInt), ALLOCATABLE :: SDISPLS2(:), RDISPLS2(:)

      integer(kind=MPIInt), ALLOCATABLE :: SENDIDX(:), RECVIDX(:)
      real(kind=wp),    ALLOCATABLE :: SENDVAL(:), RECVVAL(:)

      integer(kind=MPIInt), PARAMETER :: ONE4=1, TWO4=2
      integer(kind=MPIInt) :: IERROR4

      integer(kind=iwp), ALLOCATABLE :: IBUF(:)
      integer(kind=iwp) iG3,iT,iU,iV,iX,iY,iZ,iST,iSU,iSV,iSX,iSY,iSZ,  &
     &                  ituvs,ixyzs,iTU,iVX,iYZ,JSYM,ISUP,JSUP
      integer(kind=iwp) NG3MAX,NPROCS,                                  &
     &                  MAXMEM,iscal,MAXBUF,NG3B,NBUF,NQOT,NREM,        &
     &                  NBLOCKS,IBLOCK,IG3STA,IG3END,IROW,IP,           &
     &                  IOFFSET,I,ICOL,NRECV
      real(kind=wp)     G3VAL

#include "mpi_interfaces.fh"

      ! Since we are stuck with collective calls to MPI_Alltoallv in
      ! order to gather the elements, each process needs to loop over
      ! the same number of blocks.
      NG3MAX=NG3
      CALL GAIGOP_SCAL(NG3MAX,'max')
      IF (NG3MAX.EQ.0) RETURN

      ! basic information
      NPROCS=GA_NNODES()

      call MMA_ALLOCATE(SCOUNTS,NPROCS,Label='SCOUNTS')
      call MMA_ALLOCATE(RCOUNTS,NPROCS,Label='RCOUNTS')
      call MMA_ALLOCATE(SCOUNTS2,NPROCS,Label='SCOUNTS2')
      call MMA_ALLOCATE(RCOUNTS2,NPROCS,Label='RCOUNTS2')
      call MMA_ALLOCATE(SDISPLS,NPROCS,Label='SDISPLS')
      call MMA_ALLOCATE(RDISPLS,NPROCS,Label='RDISPLS')
      call MMA_ALLOCATE(SDISPLS2,NPROCS,Label='SDISPLS2')
      call MMA_ALLOCATE(RDISPLS2,NPROCS,Label='RDISPLS2')

      call MMA_ALLOCATE(IBUF,NPROCS,Label='IBUF')

      ! The global SA matrix has already been allocated, so we need to
      ! find out how much memory is left for buffering (4 equally sized
      ! buffers for sending and receiving values and indices)
      CALL mma_MaxDBLE(MAXMEM)
      ! we need two real and four integer values per element
      iscal = (storage_size(SENDVAL)+2*storage_size(SENDIDX)+           &
     &         storage_size(RECVVAL)+2*storage_size(RECVIDX))/          &
     &        storage_size(1.0_wp)
      !MAXBUF=MIN(NINT(0.95D0*MAXMEM)/4,2000000000/8)
      MAXBUF=MIN(NINT(0.95E0_wp*MAXMEM)/iscal,2000000000/8)

      ! Loop over blocks NG3B of NG3, so that 12*NG3B < MAXBUF/NPROCS.
      ! This guarantees that e.g. if all processes send all their data
      ! to one other, that process receives NPROCS*NG3B*12 elements
      ! in the receive buffer.
      NG3B=MAXBUF/(NPROCS*12)
      NG3B=MIN(NG3B,NG3MAX)
      CALL GAIGOP_SCAL(NG3B,'min')
      NBUF=12*NG3B

      call MMA_ALLOCATE(SENDVAL,NBUF,Label='SENDVAL')
      call MMA_ALLOCATE(SENDIDX,2*NBUF,Label='SENDIDX')

      ! Finally, we need some info on the layout of the global array in
      ! order to compute the process row of the row index.
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
          ituvs=Mul(IST,Mul(ISU,ISV))
          ixyzs=Mul(ISX,Mul(ISY,ISZ))
          if(ituvs.ne.ixyzs) CYCLE
          iTU=iT+NASHT*(iU-1)
          iVX=iV+NASHT*(iX-1)
          iYZ=iY+NASHT*(iZ-1)
          ! There are 12 equivalent cases, of which the second half
          ! reflects the S(tuv,xyz)=S(xyz,tuv) symmetry.

          ! - G(tuvxyz) -> SA(xut,vyz)
          jSYM=Mul(iSX,Mul(iSU,iST))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          if (.NOT.(iTU.eq.iVX.and.iVX.eq.iYZ)) THEN

          if (.NOT.(iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ)) THEN
          ! - G(vxtuyz) -> SA(uxv,tyz)
          jSYM=Mul(iSU,Mul(iSX,iSV))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(yzvxtu) -> SA(xzy,vtu)
          jSYM=Mul(iSX,Mul(iSZ,iSY))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(tuyzvx) -> SA(zut,yvx)
          jSYM=Mul(iSZ,Mul(iSU,iST))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ENDIF

          ! - G(yztuvx) -> SA(uzy,tvx)
          jSYM=Mul(iSU,Mul(iSZ,iSY))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(vxyztu) -> SA(zxv,ytu)
          jSYM=Mul(iSZ,Mul(iSX,iSV))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF

          ENDIF

          if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) CYCLE
          if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) CYCLE
          if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) CYCLE
          if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) CYCLE
          ! - G(utxvzy) -> SA(vtu,xzy)
          jSYM=Mul(iSV,Mul(iST,iSU))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) CYCLE
          if (.NOT.(iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ)) THEN
          ! - G(xvutzy) -> SA(tvx,uzy)
          jSYM=Mul(iST,Mul(iSV,iSX))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(zyxvut) -> SA(vyz,xut)
          jSYM=Mul(iSV,Mul(iSY,iSZ))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(utzyxv) -> SA(ytu,zxv)
          jSYM=Mul(iSY,Mul(iST,iSU))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ENDIF
          ! - G(zyutxv) -> SA(tyz,uxv)
          jSYM=Mul(iST,Mul(iSY,iSZ))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            SCOUNTS(IP)=SCOUNTS(IP)+ONE4
          ENDIF
          ! - G(xvzyut) -> SA(yvx,zut)
          jSYM=Mul(iSY,Mul(iSV,iSX))
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
          ituvs=Mul(IST,Mul(ISU,ISV))
          ixyzs=Mul(ISX,Mul(ISY,ISZ))
          if(ituvs.ne.ixyzs) CYCLE
          iTU=iT+NASHT*(iU-1)
          iVX=iV+NASHT*(iX-1)
          iYZ=iY+NASHT*(iZ-1)
          G3VAL=-G3(iG3)
!-SVC20100829: 12 equivalent cases, of which the second
!  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
!  - G(tuvxyz) -> SA(xut,vyz)
          jSYM=Mul(IASYM(iX),Mul(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
          if (.NOT.(iTU.eq.iVX.and.iVX.eq.iYZ)) THEN

          if (.NOT.(iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ)) THEN
!  - G(vxtuyz) -> SA(uxv,tyz)
          jSYM=Mul(IASYM(iU),Mul(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
!  - G(yzvxtu) -> SA(xzy,vtu)
          jSYM=Mul(IASYM(iX),Mul(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iV,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
!  - G(tuyzvx) -> SA(zut,yvx)
          jSYM=Mul(IASYM(iZ),Mul(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iY,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
          ENDIF

!  - G(yztuvx) -> SA(uzy,tvx)
          jSYM=Mul(IASYM(iU),Mul(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iT,iV,iX)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
!  - G(vxyztu) -> SA(zxv,ytu)
          jSYM=Mul(IASYM(iZ),Mul(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iY,iT,iU)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF

          ENDIF

          if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) CYCLE
          if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) CYCLE
          if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) CYCLE
          if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) CYCLE
!  - G(utxvzy) -> SA(vtu,xzy)
          jSYM=Mul(IASYM(iV),Mul(IASYM(iT),IASYM(iU)))
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

          if (.NOT.(iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ)) THEN
!  - G(xvutzy) -> SA(tvx,uzy)
          jSYM=Mul(IASYM(iT),Mul(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
!  - G(zyxvut) -> SA(vyz,xut)
          jSYM=Mul(IASYM(iV),Mul(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iX,iU,iT)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
!  - G(utzyxv) -> SA(ytu,zxv)
          jSYM=Mul(IASYM(iY),Mul(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF

          ENDIF

!  - G(zyutxv) -> SA(tyz,uxv)
          jSYM=Mul(IASYM(iT),Mul(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iU,iX,iV)-nTUVES(jSYM)
            IP=IPROW(IROW,NQOT,NREM)
            IBUF(IP)=IBUF(IP)+1
            SENDVAL(IBUF(IP))=G3VAL
            SENDIDX(2*IBUF(IP)-1)=INT(IROW,kind=MPIInt)
            SENDIDX(2*IBUF(IP))=INT(ICOL,kind=MPIInt)
          ENDIF
!  - G(xvzyut) -> SA(yvx,zut)
          jSYM=Mul(IASYM(iY),Mul(IASYM(iV),IASYM(iX)))
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
        CALL MPI_ALLTOALL(SCOUNTS, ONE4, MPI_INTEGER,                   &
     &                    RCOUNTS, ONE4, MPI_INTEGER,                   &
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

        call MMA_ALLOCATE(RECVVAL,NRECV,Label='RECVVAL')
        call MMA_ALLOCATE(RECVIDX,2*NRECV,Label='RECVIDX')

        ! Now, it is time to collect the appropriate values and indices
        ! in their respective receive buffers.
        CALL MPI_ALLTOALLV(SENDVAL, SCOUNTS, SDISPLS, MPI_REAL8,        &
     &                     RECVVAL, RCOUNTS, RDISPLS, MPI_REAL8,        &
     &                     MPI_COMM_WORLD, IERROR4)
        CALL MPI_ALLTOALLV(SENDIDX, SCOUNTS2, SDISPLS2, MPI_INTEGER,    &
     &                     RECVIDX, RCOUNTS2, RDISPLS2, MPI_INTEGER,    &
     &                     MPI_COMM_WORLD, IERROR4)

        ! Finally, fill the local chunk of the SA matrix (block of rows)
        ! with the received values at their appropriate place.
        DO I=1,NRECV
          ISUP=RECVIDX(2*I-1)
          JSUP=RECVIDX(2*I)
          SA(ISUP-ILO+1,JSUP)=RECVVAL(I)
        END DO

        call MMA_DEALLOCATE(RECVVAL)
        call MMA_DEALLOCATE(RECVIDX)

      END DO ! end loop over blocks of G3 values

      call MMA_DEALLOCATE(SENDVAL)
      call MMA_DEALLOCATE(SENDIDX)

      call MMA_DEALLOCATE(SCOUNTS)
      call MMA_DEALLOCATE(RCOUNTS)
      call MMA_DEALLOCATE(SCOUNTS2)
      call MMA_DEALLOCATE(RCOUNTS2)
      call MMA_DEALLOCATE(SDISPLS)
      call MMA_DEALLOCATE(RDISPLS)
      call MMA_DEALLOCATE(SDISPLS2)
      call MMA_DEALLOCATE(RDISPLS2)

      call MMA_DEALLOCATE(IBUF)

! Avoid unused argument warnings
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

#elif defined (NAGFOR)
      ! Some compilers do not like empty files
      subroutine empty_mksa_g3_mpp()
      end subroutine empty_mksa_g3_mpp
#endif
