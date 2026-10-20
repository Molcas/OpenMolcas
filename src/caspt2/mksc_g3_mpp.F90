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

#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

subroutine MKSC_G3_MPP(ISYM,SC,iLo,NAS,LDC,NG3,G3,idxG3)

use Symmetry_Info, only: Mul
use MPI_Wrapper, only: MPI_AllToAll, MPI_AllToAllV, MPI_COMM_WORLD, MPI_INTEGER, MPI_REAL8
use GA_Wrapper, only: GA_NNodes
use SUPERINDEX, only: KTUV
use caspt2_module, only: IASYM, NASHT, nTUVES
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use Definitions, only: wp, iwp, byte, MPIInt

implicit none
integer(kind=iwp), intent(in) :: ISYM, iLo, NAS, LDC, NG3
real(kind=wp), intent(out) :: SC(LDC,NAS)
real(kind=wp), intent(in) :: G3(NG3)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) :: I, IBLOCK, ICOL, iG3, IG3END, IG3STA, IOFFSET, IP, IROW, ISCAL, iST, iSU, ISUP, iSV, iSX, iSY, iSZ, iT, iTU, &
                     ituvs, iU, iV, iVX, iX, ixyzs, iY, iYZ, iZ, JSUP, JSYM, MAXBUF, MAXMEM, NBLOCKS, NBUF, NG3B, NG3MAX, NPROCS, &
                     NQOT, NRECV, NREM
integer(kind=MPIInt) :: IERROR4
real(kind=wp) :: G3VAL
integer(kind=iwp), allocatable :: IBUF(:)
integer(kind=MPIInt), allocatable :: RCOUNTS(:), RCOUNTS2(:), RDISPLS(:), RDISPLS2(:), RECVIDX(:), SCOUNTS(:), SCOUNTS2(:), &
                                     SDISPLS(:), SDISPLS2(:), SENDIDX(:)
real(kind=wp), allocatable :: RECVVAL(:), SENDVAL(:)
integer(kind=iwp), external :: IPROW

! Since we are stuck with collective calls to MPI_Alltoallv in
! order to gather the elements, each process needs to loop over
! the same number of blocks.
NG3MAX = NG3
call GAIGOP_SCAL(NG3MAX,'max')
if (NG3MAX == 0) return

! basic information
NPROCS = GA_NNODES()

call MMA_ALLOCATE(SCOUNTS,NPROCS,Label='SCOUNTS')
call MMA_ALLOCATE(RCOUNTS,NPROCS,Label='RCOUNTS')
call MMA_ALLOCATE(SCOUNTS2,NPROCS,Label='SCOUNTS2')
call MMA_ALLOCATE(RCOUNTS2,NPROCS,Label='RCOUNTS2')
call MMA_ALLOCATE(SDISPLS,NPROCS,Label='SDISPLS')
call MMA_ALLOCATE(RDISPLS,NPROCS,Label='RDISPLS')
call MMA_ALLOCATE(SDISPLS2,NPROCS,Label='SDISPLS2')
call MMA_ALLOCATE(RDISPLS2,NPROCS,Label='RDISPLS2')

call MMA_ALLOCATE(IBUF,NPROCS,Label='IBUF')

! The global SC matrix has already been allocated, so we need to
! find out how much memory is left for buffering (4 equally sized
! buffers for sending and receiving values and indices)
call mma_MaxDBLE(MAXMEM)
iscal = (storage_size(SENDVAL)+2*storage_size(SENDIDX)+storage_size(RECVVAL)+2*storage_size(RECVIDX))/storage_size(1.0_wp)
MAXBUF = min(nint(0.95_wp*MAXMEM)/iscal,2000000000/8)

! Loop over blocks NG3B of NG3, so that 12*NG3B < MAXBUF/NPROCS.
! This guarantees that e.g. if all processes send all their data
! to one other, that process receives NPROCS*NG3B*12 elements
! in the receive buffer.
NG3B = MAXBUF/(NPROCS*12)
NG3B = min(NG3B,NG3MAX)
call GAIGOP_SCAL(NG3B,'min')
! 12 corresponds to the number of if (jSym == iSym) branch
NBUF = 12*NG3B

call MMA_ALLOCATE(SENDVAL,NBUF,Label='SENDVAL')
call MMA_ALLOCATE(SENDIDX,2*NBUF,Label='SENDIDX')

! Finally, we need some info on the layout of the global array in
! order to compute the process row of the row index.
NQOT = NAS/NPROCS
NREM = NAS-NPROCS*NQOT

NBLOCKS = (NG3MAX-1)/NG3B+1
do IBLOCK=1,NBLOCKS
  IG3STA = 1+(IBLOCK-1)*NG3B
  IG3END = min(IG3STA+NG3B-1,NG3)

  SCOUNTS = 0
  ! First pass to determine how many values will need to be sent
  ! to other processes. This is necessary to be able to allocate
  ! the buffer size and offsets.
  do iG3=IG3STA,IG3END
    iT = idxG3(1,iG3)
    iU = idxG3(2,iG3)
    iV = idxG3(3,iG3)
    iX = idxG3(4,iG3)
    iY = idxG3(5,iG3)
    iZ = idxG3(6,iG3)
    iST = IASYM(iT)
    iSU = IASYM(iU)
    iSV = IASYM(iV)
    iSX = IASYM(iX)
    iSY = IASYM(iY)
    iSZ = IASYM(iZ)
    ituvs = Mul(IST,Mul(ISU,ISV))
    ixyzs = Mul(ISX,Mul(ISY,ISZ))
    if (ituvs /= ixyzs) cycle
    iTU = iT+NASHT*(iU-1)
    iVX = iV+NASHT*(iX-1)
    iYZ = iY+NASHT*(iZ-1)
    !-SVC20100829: 12 equivalent cases, of which the second
    !  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
    !  - G(tuvxyz) -> SC(vut,xyz)
    jSYM = Mul(IASYM(iV),Mul(IASYM(iU),IASYM(iT)))
    if (jSYM == iSYM) then
      IROW = KTUV(iV,iU,iT)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
    end if
    if (.not. ((iTU == iVX) .and. (iVX == iYZ))) then
      if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
        !  - G(vxtuyz) -> SC(txv,uyz)
        jSYM = Mul(IASYM(iT),Mul(IASYM(iX),IASYM(iV)))
        if (jSYM == iSYM) then
          IROW = KTUV(iT,iX,iV)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
        end if
        !  - G(yzvxtu) -> SC(vzy,xtu)
        jSYM = Mul(IASYM(iV),Mul(IASYM(iZ),IASYM(iY)))
        if (jSYM == iSYM) then
          IROW = KTUV(iV,iZ,iY)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
        end if
        !  - G(tuyzvx) -> SC(yut,zvx)
        jSYM = Mul(IASYM(iY),Mul(IASYM(iU),IASYM(iT)))
        if (jSYM == iSYM) then
          IROW = KTUV(iY,iU,iT)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
        end if
      end if
      !  - G(yztuvx) -> SC(tzy,uvx)
      jSYM = Mul(IASYM(iT),Mul(IASYM(iZ),IASYM(iY)))
      if (jSYM == iSYM) then
        IROW = KTUV(iT,iZ,iY)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
      end if
      !  - G(vxyztu) -> SC(yxv,ztu)
      jSYM = Mul(IASYM(iY),Mul(IASYM(iX),IASYM(iV)))
      if (jSYM == iSYM) then
        IROW = KTUV(iY,iX,iV)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
      end if
    end if
    if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
    if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
    if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
    if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle
    !  - G(utxvzy) -> SC(xtu,vzy)
    jSYM = Mul(IASYM(iX),Mul(IASYM(iT),IASYM(iU)))
    if (jSYM == iSYM) then
      IROW = KTUV(iX,iT,iU)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
    end if
    if ((iTU == iVX) .and. (iVX == iYZ)) cycle
    if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
      !  - G(xvutzy) -> SC(uvx,tzy)
      jSYM = Mul(IASYM(iU),Mul(IASYM(iV),IASYM(iX)))
      if (jSYM == iSYM) then
        IROW = KTUV(iU,iV,iX)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
      end if
      !  - G(zyxvut) -> SC(xyz,vut)
      jSYM = Mul(IASYM(iX),Mul(IASYM(iY),IASYM(iZ)))
      if (jSYM == iSYM) then
        IROW = KTUV(iX,iY,iZ)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
      end if
      !  - G(utzyxv) -> SC(ztu,yxv)
      jSYM = Mul(IASYM(iZ),Mul(IASYM(iT),IASYM(iU)))
      if (jSYM == iSYM) then
        IROW = KTUV(iZ,iT,iU)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
      end if
    end if
    !  - G(zyutxv) -> SC(uyz,txv)
    jSYM = Mul(IASYM(iU),Mul(IASYM(iY),IASYM(iZ)))
    if (jSYM == iSYM) then
      IROW = KTUV(iU,iY,iZ)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
    end if
    !  - G(xvzyut) -> SC(zvx,yut)
    jSYM = Mul(IASYM(iZ),Mul(IASYM(iV),IASYM(iX)))
    if (jSYM == iSYM) then
      IROW = KTUV(iZ,iV,iX)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
    end if
  end do

  ! At this point, SCOUNTS contains the number of values generated
  ! for each process. Use them to determine the send offsets.
  IOFFSET = 0
  do I=1,NPROCS
    SDISPLS(I) = int(IOFFSET,kind=MPIInt)
    IBUF(I) = IOFFSET
    IOFFSET = IOFFSET+SCOUNTS(I)
  end do

  ! Second pass fills the buffers with values and indices
  do iG3=IG3STA,IG3END
    iT = idxG3(1,iG3)
    iU = idxG3(2,iG3)
    iV = idxG3(3,iG3)
    iX = idxG3(4,iG3)
    iY = idxG3(5,iG3)
    iZ = idxG3(6,iG3)
    iST = IASYM(iT)
    iSU = IASYM(iU)
    iSV = IASYM(iV)
    iSX = IASYM(iX)
    iSY = IASYM(iY)
    iSZ = IASYM(iZ)
    ituvs = Mul(IST,Mul(ISU,ISV))
    ixyzs = Mul(ISX,Mul(ISY,ISZ))
    if (ituvs /= ixyzs) cycle
    iTU = iT+NASHT*(iU-1)
    iVX = iV+NASHT*(iX-1)
    iYZ = iY+NASHT*(iZ-1)
    G3VAL = G3(iG3)
    !-SVC20100829: 12 equivalent cases, of which the second
    !  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
    !  - G(tuvxyz) -> SC(vut,xyz)
    jSYM = Mul(IASYM(iV),Mul(IASYM(iU),IASYM(iT)))
    if (jSYM == iSYM) then
      IROW = KTUV(iV,iU,iT)-nTUVES(jSYM)
      ICOL = KTUV(iX,iY,iZ)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = G3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
    if (.not. ((iTU == iVX) .and. (iVX == iYZ))) then
      if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
        !  - G(vxtuyz) -> SC(txv,uyz)
        jSYM = Mul(IASYM(iT),Mul(IASYM(iX),IASYM(iV)))
        if (jSYM == iSYM) then
          IROW = KTUV(iT,iX,iV)-nTUVES(jSYM)
          ICOL = KTUV(iU,iY,iZ)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          IBUF(IP) = IBUF(IP)+1
          SENDVAL(IBUF(IP)) = G3VAL
          SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
          SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
        end if
        !  - G(yzvxtu) -> SC(vzy,xtu)
        jSYM = Mul(IASYM(iV),Mul(IASYM(iZ),IASYM(iY)))
        if (jSYM == iSYM) then
          IROW = KTUV(iV,iZ,iY)-nTUVES(jSYM)
          ICOL = KTUV(iX,iT,iU)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          IBUF(IP) = IBUF(IP)+1
          SENDVAL(IBUF(IP)) = G3VAL
          SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
          SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
        end if
        !  - G(tuyzvx) -> SC(yut,zvx)
        jSYM = Mul(IASYM(iY),Mul(IASYM(iU),IASYM(iT)))
        if (jSYM == iSYM) then
          IROW = KTUV(iY,iU,iT)-nTUVES(jSYM)
          ICOL = KTUV(iZ,iV,iX)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          IBUF(IP) = IBUF(IP)+1
          SENDVAL(IBUF(IP)) = G3VAL
          SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
          SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
        end if
      end if
      !  - G(yztuvx) -> SC(tzy,uvx)
      jSYM = Mul(IASYM(iT),Mul(IASYM(iZ),IASYM(iY)))
      if (jSYM == iSYM) then
        IROW = KTUV(iT,iZ,iY)-nTUVES(jSYM)
        ICOL = KTUV(iU,iV,iX)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = G3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
      !  - G(vxyztu) -> SC(yxv,ztu)
      jSYM = Mul(IASYM(iY),Mul(IASYM(iX),IASYM(iV)))
      if (jSYM == iSYM) then
        IROW = KTUV(iY,iX,iV)-nTUVES(jSYM)
        ICOL = KTUV(iZ,iT,iU)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = G3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
    end if
    if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
    if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
    if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
    if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle
    !  - G(utxvzy) -> SC(xtu,vzy)
    jSYM = Mul(IASYM(iX),Mul(IASYM(iT),IASYM(iU)))
    if (jSYM == iSYM) then
      IROW = KTUV(iX,iT,iU)-nTUVES(jSYM)
      ICOL = KTUV(iV,iZ,iY)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = G3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
    if ((iTU == iVX) .and. (iVX == iYZ)) cycle
    if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
      !  - G(xvutzy) -> SC(uvx,tzy)
      jSYM = Mul(IASYM(iU),Mul(IASYM(iV),IASYM(iX)))
      if (jSYM == iSYM) then
        IROW = KTUV(iU,iV,iX)-nTUVES(jSYM)
        ICOL = KTUV(iT,iZ,iY)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = G3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
      !  - G(zyxvut) -> SC(xyz,vut)
      jSYM = Mul(IASYM(iX),Mul(IASYM(iY),IASYM(iZ)))
      if (jSYM == iSYM) then
        IROW = KTUV(iX,iY,iZ)-nTUVES(jSYM)
        ICOL = KTUV(iV,iU,iT)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = G3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
      !  - G(utzyxv) -> SC(ztu,yxv)
      jSYM = Mul(IASYM(iZ),Mul(IASYM(iT),IASYM(iU)))
      if (jSYM == iSYM) then
        IROW = KTUV(iZ,iT,iU)-nTUVES(jSYM)
        ICOL = KTUV(iY,iX,iV)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = G3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
    end if
    !  - G(zyutxv) -> SC(uyz,txv)
    jSYM = Mul(IASYM(iU),Mul(IASYM(iY),IASYM(iZ)))
    if (jSYM == iSYM) then
      IROW = KTUV(iU,iY,iZ)-nTUVES(jSYM)
      ICOL = KTUV(iT,iX,iV)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = G3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
    !  - G(xvzyut) -> SC(zvx,yut)
    jSYM = Mul(IASYM(iZ),Mul(IASYM(iV),IASYM(iX)))
    if (jSYM == iSYM) then
      IROW = KTUV(iZ,iV,iX)-nTUVES(jSYM)
      ICOL = KTUV(iY,iU,iT)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = G3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
  end do

  ! Now we need to determine the receive counts.
  call MPI_ALLTOALL(SCOUNTS,1_MPIInt,MPI_INTEGER,RCOUNTS,1_MPIInt,MPI_INTEGER,MPI_COMM_WORLD,IERROR4)

  IOFFSET = 0
  do I=1,NPROCS
    RDISPLS(I) = int(IOFFSET,kind=MPIInt)
    IOFFSET = IOFFSET+RCOUNTS(I)
  end do
  SCOUNTS2(:) = 2_MPIInt*SCOUNTS(:)
  RCOUNTS2(:) = 2_MPIInt*RCOUNTS(:)
  SDISPLS2(:) = 2_MPIInt*SDISPLS(:)
  RDISPLS2(:) = 2_MPIInt*RDISPLS(:)
  NRECV = IOFFSET

  call MMA_ALLOCATE(RECVVAL,NRECV,Label='RECVVAL')
  call MMA_ALLOCATE(RECVIDX,2*NRECV,Label='RECVIDX')

  ! Now, it is time to collect the appropriate values and indices
  ! in their respective receive buffers.
  call MPI_ALLTOALLV(SENDVAL,SCOUNTS,SDISPLS,MPI_REAL8,RECVVAL,RCOUNTS,RDISPLS,MPI_REAL8,MPI_COMM_WORLD,IERROR4)
  call MPI_ALLTOALLV(SENDIDX,SCOUNTS2,SDISPLS2,MPI_INTEGER,RECVIDX,RCOUNTS2,RDISPLS2,MPI_INTEGER,MPI_COMM_WORLD,IERROR4)

  ! Finally, fill the local chunk of the SC matrix (block of rows)
  ! with the received values at their appropriate place.
  do I=1,NRECV
    ISUP = RECVIDX(2*I-1)
    JSUP = RECVIDX(2*I)
    SC(ISUP-ILO+1,JSUP) = RECVVAL(I)
  end do

  call MMA_DEALLOCATE(RECVVAL)
  call MMA_DEALLOCATE(RECVIDX)

end do ! end loop over blocks of G3 values

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

end subroutine MKSC_G3_MPP

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(MKSC_G3_MPP)

#endif
