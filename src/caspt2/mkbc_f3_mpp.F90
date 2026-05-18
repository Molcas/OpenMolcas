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
subroutine MKBC_F3_MPP(ISYM,BC,iLo,iHi,NAS,LDC,NG3,F3,idxG3)

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp, Byte
use MPI
use SUPERINDEX, only: KTUV
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use definitions, only: MPIInt, wp
use caspt2_module, only: IASYM, NASHT, NTUVES

implicit none
#include "global.fh"
#include "mafdecls.fh"
integer(kind=iwp), intent(in) :: ISYM, iLo, iHi, NAS, LDC, NG3
real(kind=wp), intent(inout) :: BC(LDC,NAS)
real(kind=wp), intent(in) :: F3(NG3)
integer(kind=Byte), intent(in) :: idxG3(6,NG3)
integer(kind=MPIInt), allocatable :: SCOUNTS(:), RCOUNTS(:)
integer(kind=MPIInt), allocatable :: SCOUNTS2(:), RCOUNTS2(:)
integer(kind=MPIInt), allocatable :: SDISPLS(:), RDISPLS(:)
integer(kind=MPIInt), allocatable :: SDISPLS2(:), RDISPLS2(:)
integer(kind=MPIInt), allocatable :: SENDIDX(:), RECVIDX(:)
real(kind=wp), allocatable :: SENDVAL(:), RECVVAL(:)
integer(kind=MPIInt), parameter :: ONE4 = 1, TWO4 = 2
integer(kind=MPIInt) :: IERROR4
integer(kind=iwp), allocatable :: IBUF(:)
integer(kind=iwp) NG3MAX, NPROCS, iscal, MAXBUF, NG3B, NBUF, NQOT, NREM, NBLOCKS, IBLOCK, IG3STA, IG3END, MAXMEM, IT, IG3, IU, IV, &
                  IX, IY, IZ, IST, ISU, ISV, ISX, ISY, ISZ, ITUVS, IXYZS, ITU, IVX, IYZ, JSYM, IROW, IP, IOFFSET, I, ICOL, NRECV, &
                  ISUP, JSUP
real(kind=wp) F3VAL
integer(kind=iwp), external :: IPROW
#include "mpi_interfaces.fh"

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
MAXBUF = min(nint(0.95d0*MAXMEM)/iscal,2000000000/8)

! Loop over blocks NG3B of NG3, so that 12*NG3B < MAXBUF/NPROCS.
! This guarantees that e.g. if all processes send all their data
! to one other, that process receives NPROCS*NG3B*12 elements
! in the receive buffer.
NG3B = MAXBUF/(NPROCS*12)
NG3B = min(NG3B,NG3MAX)
call GAIGOP_SCAL(NG3B,'min')
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
    !  half reflects the B(tuv,xyz)=B(xyz,tuv) symmetry:
    !  - F(tuvxyz) -> BC(vut,xyz)
    jSYM = Mul(IASYM(iV),Mul(IASYM(iU),IASYM(iT)))
    if (jSYM == iSYM) then
      IROW = KTUV(iV,iU,iT)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+ONE4
    end if
    if (.not. ((iTU == iVX) .and. (iVX == iYZ))) then
      if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
        !  - F(vxtuyz) -> BC(txv,uyz)
        jSYM = Mul(IASYM(iT),Mul(IASYM(iX),IASYM(iV)))
        if (jSYM == iSYM) then
          IROW = KTUV(iT,iX,iV)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          SCOUNTS(IP) = SCOUNTS(IP)+ONE4
        end if
        !  - F(yzvxtu) -> BC(vzy,xtu)
        jSYM = Mul(IASYM(iV),Mul(IASYM(iZ),IASYM(iY)))
        if (jSYM == iSYM) then
          IROW = KTUV(iV,iZ,iY)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          SCOUNTS(IP) = SCOUNTS(IP)+ONE4
        end if
        !  - F(tuyzvx) -> BC(yut,zvx)
        jSYM = Mul(IASYM(iY),Mul(IASYM(iU),IASYM(iT)))
        if (jSYM == iSYM) then
          IROW = KTUV(iY,iU,iT)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          SCOUNTS(IP) = SCOUNTS(IP)+ONE4
        end if
      end if
      !  - F(yztuvx) -> BC(tzy,uvx)
      jSYM = Mul(IASYM(iT),Mul(IASYM(iZ),IASYM(iY)))
      if (jSYM == iSYM) then
        IROW = KTUV(iT,iZ,iY)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+ONE4
      end if
      !  - F(vxyztu) -> BC(yxv,ztu)
      jSYM = Mul(IASYM(iY),Mul(IASYM(iX),IASYM(iV)))
      if (jSYM == iSYM) then
        IROW = KTUV(iY,iX,iV)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+ONE4
      end if
    end if
    if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
    if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
    if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
    if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle
    !  - F(utxvzy) -> BC(xtu,vzy)
    jSYM = Mul(IASYM(iX),Mul(IASYM(iT),IASYM(iU)))
    if (jSYM == iSYM) then
      IROW = KTUV(iX,iT,iU)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+ONE4
    end if
    if ((iTU == iVX) .and. (iVX == iYZ)) cycle
    if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
      !  - F(xvutzy) -> BC(uvx,tzy)
      jSYM = Mul(IASYM(iU),Mul(IASYM(iV),IASYM(iX)))
      if (jSYM == iSYM) then
        IROW = KTUV(iU,iV,iX)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+ONE4
      end if
      !  - F(zyxvut) -> BC(xyz,vut)
      jSYM = Mul(IASYM(iX),Mul(IASYM(iY),IASYM(iZ)))
      if (jSYM == iSYM) then
        IROW = KTUV(iX,iY,iZ)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+ONE4
      end if
      !  - F(utzyxv) -> BC(ztu,yxv)
      jSYM = Mul(IASYM(iZ),Mul(IASYM(iT),IASYM(iU)))
      if (jSYM == iSYM) then
        IROW = KTUV(iZ,iT,iU)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+ONE4
      end if
    end if
    !  - F(zyutxv) -> BC(uyz,txv)
    jSYM = Mul(IASYM(iU),Mul(IASYM(iY),IASYM(iZ)))
    if (jSYM == iSYM) then
      IROW = KTUV(iU,iY,iZ)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+ONE4
    end if
    !  - F(xvzyut) -> BC(zvx,yut)
    jSYM = Mul(IASYM(iZ),Mul(IASYM(iV),IASYM(iX)))
    if (jSYM == iSYM) then
      IROW = KTUV(iZ,iV,iX)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+ONE4
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
    F3VAL = F3(iG3)
    !-SVC20100829: 12 equivalent cases, of which the second
    !  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
    !  - F(tuvxyz) -> BC(vut,xyz)
    jSYM = Mul(IASYM(iV),Mul(IASYM(iU),IASYM(iT)))
    if (jSYM == iSYM) then
      IROW = KTUV(iV,iU,iT)-nTUVES(jSYM)
      ICOL = KTUV(iX,iY,iZ)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = F3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
    if (.not. ((iTU == iVX) .and. (iVX == iYZ))) then
      if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
        !  - F(vxtuyz) -> BC(txv,uyz)
        jSYM = Mul(IASYM(iT),Mul(IASYM(iX),IASYM(iV)))
        if (jSYM == iSYM) then
          IROW = KTUV(iT,iX,iV)-nTUVES(jSYM)
          ICOL = KTUV(iU,iY,iZ)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          IBUF(IP) = IBUF(IP)+1
          SENDVAL(IBUF(IP)) = F3VAL
          SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
          SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
        end if
        !  - F(yzvxtu) -> BC(vzy,xtu)
        jSYM = Mul(IASYM(iV),Mul(IASYM(iZ),IASYM(iY)))
        if (jSYM == iSYM) then
          IROW = KTUV(iV,iZ,iY)-nTUVES(jSYM)
          ICOL = KTUV(iX,iT,iU)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          IBUF(IP) = IBUF(IP)+1
          SENDVAL(IBUF(IP)) = F3VAL
          SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
          SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
        end if
        !  - F(tuyzvx) -> BC(yut,zvx)
        jSYM = Mul(IASYM(iY),Mul(IASYM(iU),IASYM(iT)))
        if (jSYM == iSYM) then
          IROW = KTUV(iY,iU,iT)-nTUVES(jSYM)
          ICOL = KTUV(iZ,iV,iX)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          IBUF(IP) = IBUF(IP)+1
          SENDVAL(IBUF(IP)) = F3VAL
          SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
          SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
        end if
      end if
      !  - F(yztuvx) -> BC(tzy,uvx)
      jSYM = Mul(IASYM(iT),Mul(IASYM(iZ),IASYM(iY)))
      if (jSYM == iSYM) then
        IROW = KTUV(iT,iZ,iY)-nTUVES(jSYM)
        ICOL = KTUV(iU,iV,iX)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = F3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
      !  - F(vxyztu) -> BC(yxv,ztu)
      jSYM = Mul(IASYM(iY),Mul(IASYM(iX),IASYM(iV)))
      if (jSYM == iSYM) then
        IROW = KTUV(iY,iX,iV)-nTUVES(jSYM)
        ICOL = KTUV(iZ,iT,iU)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = F3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
    end if
    if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
    if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
    if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
    if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle
    !  - F(utxvzy) -> BC(xtu,vzy)
    jSYM = Mul(IASYM(iX),Mul(IASYM(iT),IASYM(iU)))
    if (jSYM == iSYM) then
      IROW = KTUV(iX,iT,iU)-nTUVES(jSYM)
      ICOL = KTUV(iV,iZ,iY)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = F3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
    if ((iTU == iVX) .and. (iVX == iYZ)) cycle
    if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
      !  - F(xvutzy) -> BC(uvx,tzy)
      jSYM = Mul(IASYM(iU),Mul(IASYM(iV),IASYM(iX)))
      if (jSYM == iSYM) then
        IROW = KTUV(iU,iV,iX)-nTUVES(jSYM)
        ICOL = KTUV(iT,iZ,iY)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = F3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
      !  - F(zyxvut) -> BC(xyz,vut)
      jSYM = Mul(IASYM(iX),Mul(IASYM(iY),IASYM(iZ)))
      if (jSYM == iSYM) then
        IROW = KTUV(iX,iY,iZ)-nTUVES(jSYM)
        ICOL = KTUV(iV,iU,iT)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = F3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
      !  - F(utzyxv) -> BC(ztu,yxv)
      jSYM = Mul(IASYM(iZ),Mul(IASYM(iT),IASYM(iU)))
      if (jSYM == iSYM) then
        IROW = KTUV(iZ,iT,iU)-nTUVES(jSYM)
        ICOL = KTUV(iY,iX,iV)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = F3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
    end if
    !  - F(zyutxv) -> BC(uyz,txv)
    jSYM = Mul(IASYM(iU),Mul(IASYM(iY),IASYM(iZ)))
    if (jSYM == iSYM) then
      IROW = KTUV(iU,iY,iZ)-nTUVES(jSYM)
      ICOL = KTUV(iT,iX,iV)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = F3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
    !  - F(xvzyut) -> BC(zvx,yut)
    jSYM = Mul(IASYM(iZ),Mul(IASYM(iV),IASYM(iX)))
    if (jSYM == iSYM) then
      IROW = KTUV(iZ,iV,iX)-nTUVES(jSYM)
      ICOL = KTUV(iY,iU,iT)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = F3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
  end do

  ! Now we need to determine the receive counts.
  call MPI_ALLTOALL(SCOUNTS,ONE4,MPI_INTEGER,RCOUNTS,ONE4,MPI_INTEGER,MPI_COMM_WORLD,IERROR4)

  IOFFSET = 0
  do I=1,NPROCS
    RDISPLS(I) = int(IOFFSET,kind=MPIInt)
    IOFFSET = IOFFSET+RCOUNTS(I)
    SCOUNTS2(I) = TWO4*SCOUNTS(I)
    RCOUNTS2(I) = TWO4*RCOUNTS(I)
    SDISPLS2(I) = TWO4*SDISPLS(I)
    RDISPLS2(I) = TWO4*RDISPLS(I)
  end do
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
    BC(ISUP-ILO+1,JSUP) = BC(ISUP-ILO+1,JSUP)+RECVVAL(I)
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

return
! Avoid unused argument warnings
if (.false.) call UNUSED_INTEGER(iHi)

end subroutine MKBC_F3_MPP

#elif defined (NAGFOR)

! Some compilers do not like empty files
subroutine empty_mkbc_f3_mpp()
end subroutine empty_mkbc_f3_mpp

#endif
