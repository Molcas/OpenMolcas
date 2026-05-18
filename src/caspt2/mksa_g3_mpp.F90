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
subroutine MKSA_G3_MPP(ISYM,SA,iLo,iHi,NAS,LDA,NG3,G3,idxG3)

use Symmetry_Info, only: Mul
use definitions, only: iwp, wp, MPIInt, Byte
use MPI
use SUPERINDEX, only: KTUV
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use caspt2_module, only: IASYM, NASHT, nTUVES

implicit none
#include "global.fh"
#include "mafdecls.fh"
integer(kind=iwp), intent(in) :: ISYM, iLo, iHi, NAS, LDA, NG3
real(kind=wp), intent(out) :: SA(LDA,NAS)
real(kind=wp), intent(in) :: G3(NG3)
integer(kind=Byte), intent(in) :: idxG3(6,NG3)
integer(kind=MPIInt), allocatable :: SCOUNTS(:), RCOUNTS(:)
integer(kind=MPIInt), allocatable :: SCOUNTS2(:), RCOUNTS2(:)
integer(kind=MPIInt), allocatable :: SDISPLS(:), RDISPLS(:)
integer(kind=MPIInt), allocatable :: SDISPLS2(:), RDISPLS2(:)
integer(kind=MPIInt), allocatable :: SENDIDX(:), RECVIDX(:)
real(kind=wp), allocatable :: SENDVAL(:), RECVVAL(:)
integer(kind=MPIInt) :: IERROR4
integer(kind=iwp), allocatable :: IBUF(:)
integer(kind=iwp) iG3, iT, iU, iV, iX, iY, iZ, iST, iSU, iSV, iSX, iSY, iSZ, ituvs, ixyzs, iTU, iVX, iYZ, JSYM, ISUP, JSUP
integer(kind=iwp) NG3MAX, NPROCS, MAXMEM, iscal, MAXBUF, NG3B, NBUF, NQOT, NREM, NBLOCKS, IBLOCK, IG3STA, IG3END, IROW, IP, &
                  IOFFSET, I, ICOL, NRECV
real(kind=wp) G3VAL
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

! The global SA matrix has already been allocated, so we need to
! find out how much memory is left for buffering (4 equally sized
! buffers for sending and receiving values and indices)
call mma_MaxDBLE(MAXMEM)
! we need two real and four integer values per element
iscal = (storage_size(SENDVAL)+2*storage_size(SENDIDX)+storage_size(RECVVAL)+2*storage_size(RECVIDX))/storage_size(1.0_wp)
!MAXBUF=MIN(NINT(0.95_wp*MAXMEM)/4,2000000000/8)
MAXBUF = min(nint(0.95_wp*MAXMEM)/iscal,2000000000/8)

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
    ! There are 12 equivalent cases, of which the second half
    ! reflects the S(tuv,xyz)=S(xyz,tuv) symmetry.

    ! - G(tuvxyz) -> SA(xut,vyz)
    jSYM = Mul(iSX,Mul(iSU,iST))
    if (jSYM == iSYM) then
      IROW = KTUV(iX,iU,iT)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
    end if
    if (.not. ((iTU == iVX) .and. (iVX == iYZ))) then

      if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
        ! - G(vxtuyz) -> SA(uxv,tyz)
        jSYM = Mul(iSU,Mul(iSX,iSV))
        if (jSYM == iSYM) then
          IROW = KTUV(iU,iX,iV)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
        end if
        ! - G(yzvxtu) -> SA(xzy,vtu)
        jSYM = Mul(iSX,Mul(iSZ,iSY))
        if (jSYM == iSYM) then
          IROW = KTUV(iX,iZ,iY)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
        end if
        ! - G(tuyzvx) -> SA(zut,yvx)
        jSYM = Mul(iSZ,Mul(iSU,iST))
        if (jSYM == iSYM) then
          IROW = KTUV(iZ,iU,iT)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
        end if
      end if

      ! - G(yztuvx) -> SA(uzy,tvx)
      jSYM = Mul(iSU,Mul(iSZ,iSY))
      if (jSYM == iSYM) then
        IROW = KTUV(iU,iZ,iY)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
      end if
      ! - G(vxyztu) -> SA(zxv,ytu)
      jSYM = Mul(iSZ,Mul(iSX,iSV))
      if (jSYM == iSYM) then
        IROW = KTUV(iZ,iX,iV)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
      end if

    end if

    if ((iT == iU) .and. (iV == iX) .and. (iY == iZ)) cycle
    if ((iT == iU) .and. (iV == iZ) .and. (iX == iY)) cycle
    if ((iX == iV) .and. (iT == iZ) .and. (iU == iY)) cycle
    if ((iZ == iY) .and. (iV == iU) .and. (iX == iT)) cycle
    ! - G(utxvzy) -> SA(vtu,xzy)
    jSYM = Mul(iSV,Mul(iST,iSU))
    if (jSYM == iSYM) then
      IROW = KTUV(iV,iT,iU)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
    end if
    if ((iTU == iVX) .and. (iVX == iYZ)) cycle
    if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
      ! - G(xvutzy) -> SA(tvx,uzy)
      jSYM = Mul(iST,Mul(iSV,iSX))
      if (jSYM == iSYM) then
        IROW = KTUV(iT,iV,iX)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
      end if
      ! - G(zyxvut) -> SA(vyz,xut)
      jSYM = Mul(iSV,Mul(iSY,iSZ))
      if (jSYM == iSYM) then
        IROW = KTUV(iV,iY,iZ)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
      end if
      ! - G(utzyxv) -> SA(ytu,zxv)
      jSYM = Mul(iSY,Mul(iST,iSU))
      if (jSYM == iSYM) then
        IROW = KTUV(iY,iT,iU)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
      end if
    end if
    ! - G(zyutxv) -> SA(tyz,uxv)
    jSYM = Mul(iST,Mul(iSY,iSZ))
    if (jSYM == iSYM) then
      IROW = KTUV(iT,iY,iZ)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      SCOUNTS(IP) = SCOUNTS(IP)+1_MPIInt
    end if
    ! - G(xvzyut) -> SA(yvx,zut)
    jSYM = Mul(iSY,Mul(iSV,iSX))
    if (jSYM == iSYM) then
      IROW = KTUV(iY,iV,iX)-nTUVES(jSYM)
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
    G3VAL = -G3(iG3)
    !-SVC20100829: 12 equivalent cases, of which the second
    !  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
    !  - G(tuvxyz) -> SA(xut,vyz)
    jSYM = Mul(IASYM(iX),Mul(IASYM(iU),IASYM(iT)))
    if (jSYM == iSYM) then
      IROW = KTUV(iX,iU,iT)-nTUVES(jSYM)
      ICOL = KTUV(iV,iY,iZ)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = G3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
    if (.not. ((iTU == iVX) .and. (iVX == iYZ))) then

      if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
        !  - G(vxtuyz) -> SA(uxv,tyz)
        jSYM = Mul(IASYM(iU),Mul(IASYM(iX),IASYM(iV)))
        if (jSYM == iSYM) then
          IROW = KTUV(iU,iX,iV)-nTUVES(jSYM)
          ICOL = KTUV(iT,iY,iZ)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          IBUF(IP) = IBUF(IP)+1
          SENDVAL(IBUF(IP)) = G3VAL
          SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
          SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
        end if
        !  - G(yzvxtu) -> SA(xzy,vtu)
        jSYM = Mul(IASYM(iX),Mul(IASYM(iZ),IASYM(iY)))
        if (jSYM == iSYM) then
          IROW = KTUV(iX,iZ,iY)-nTUVES(jSYM)
          ICOL = KTUV(iV,iT,iU)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          IBUF(IP) = IBUF(IP)+1
          SENDVAL(IBUF(IP)) = G3VAL
          SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
          SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
        end if
        !  - G(tuyzvx) -> SA(zut,yvx)
        jSYM = Mul(IASYM(iZ),Mul(IASYM(iU),IASYM(iT)))
        if (jSYM == iSYM) then
          IROW = KTUV(iZ,iU,iT)-nTUVES(jSYM)
          ICOL = KTUV(iY,iV,iX)-nTUVES(jSYM)
          IP = IPROW(IROW,NQOT,NREM)
          IBUF(IP) = IBUF(IP)+1
          SENDVAL(IBUF(IP)) = G3VAL
          SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
          SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
        end if
      end if

      !  - G(yztuvx) -> SA(uzy,tvx)
      jSYM = Mul(IASYM(iU),Mul(IASYM(iZ),IASYM(iY)))
      if (jSYM == iSYM) then
        IROW = KTUV(iU,iZ,iY)-nTUVES(jSYM)
        ICOL = KTUV(iT,iV,iX)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = G3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
      !  - G(vxyztu) -> SA(zxv,ytu)
      jSYM = Mul(IASYM(iZ),Mul(IASYM(iX),IASYM(iV)))
      if (jSYM == iSYM) then
        IROW = KTUV(iZ,iX,iV)-nTUVES(jSYM)
        ICOL = KTUV(iY,iT,iU)-nTUVES(jSYM)
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
    !  - G(utxvzy) -> SA(vtu,xzy)
    jSYM = Mul(IASYM(iV),Mul(IASYM(iT),IASYM(iU)))
    if (jSYM == iSYM) then
      IROW = KTUV(iV,iT,iU)-nTUVES(jSYM)
      ICOL = KTUV(iX,iZ,iY)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = G3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
    if ((iTU == iVX) .and. (iVX == iYZ)) cycle

    if (.not. ((iTU == iVX) .or. (iTU == iYZ) .or. (iVX == iYZ))) then
      !  - G(xvutzy) -> SA(tvx,uzy)
      jSYM = Mul(IASYM(iT),Mul(IASYM(iV),IASYM(iX)))
      if (jSYM == iSYM) then
        IROW = KTUV(iT,iV,iX)-nTUVES(jSYM)
        ICOL = KTUV(iU,iZ,iY)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = G3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
      !  - G(zyxvut) -> SA(vyz,xut)
      jSYM = Mul(IASYM(iV),Mul(IASYM(iY),IASYM(iZ)))
      if (jSYM == iSYM) then
        IROW = KTUV(iV,iY,iZ)-nTUVES(jSYM)
        ICOL = KTUV(iX,iU,iT)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = G3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if
      !  - G(utzyxv) -> SA(ytu,zxv)
      jSYM = Mul(IASYM(iY),Mul(IASYM(iT),IASYM(iU)))
      if (jSYM == iSYM) then
        IROW = KTUV(iY,iT,iU)-nTUVES(jSYM)
        ICOL = KTUV(iZ,iX,iV)-nTUVES(jSYM)
        IP = IPROW(IROW,NQOT,NREM)
        IBUF(IP) = IBUF(IP)+1
        SENDVAL(IBUF(IP)) = G3VAL
        SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
        SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
      end if

    end if

    !  - G(zyutxv) -> SA(tyz,uxv)
    jSYM = Mul(IASYM(iT),Mul(IASYM(iY),IASYM(iZ)))
    if (jSYM == iSYM) then
      IROW = KTUV(iT,iY,iZ)-nTUVES(jSYM)
      ICOL = KTUV(iU,iX,iV)-nTUVES(jSYM)
      IP = IPROW(IROW,NQOT,NREM)
      IBUF(IP) = IBUF(IP)+1
      SENDVAL(IBUF(IP)) = G3VAL
      SENDIDX(2*IBUF(IP)-1) = int(IROW,kind=MPIInt)
      SENDIDX(2*IBUF(IP)) = int(ICOL,kind=MPIInt)
    end if
    !  - G(xvzyut) -> SA(yvx,zut)
    jSYM = Mul(IASYM(iY),Mul(IASYM(iV),IASYM(iX)))
    if (jSYM == iSYM) then
      IROW = KTUV(iY,iV,iX)-nTUVES(jSYM)
      ICOL = KTUV(iZ,iU,iT)-nTUVES(jSYM)
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
    SCOUNTS2(I) = 2_MPIInt*SCOUNTS(I)
    RCOUNTS2(I) = 2_MPIInt*RCOUNTS(I)
    SDISPLS2(I) = 2_MPIInt*SDISPLS(I)
    RDISPLS2(I) = 2_MPIInt*RDISPLS(I)
  end do
  NRECV = IOFFSET

  call MMA_ALLOCATE(RECVVAL,NRECV,Label='RECVVAL')
  call MMA_ALLOCATE(RECVIDX,2*NRECV,Label='RECVIDX')

  ! Now, it is time to collect the appropriate values and indices
  ! in their respective receive buffers.
  call MPI_ALLTOALLV(SENDVAL,SCOUNTS,SDISPLS,MPI_REAL8,RECVVAL,RCOUNTS,RDISPLS,MPI_REAL8,MPI_COMM_WORLD,IERROR4)
  call MPI_ALLTOALLV(SENDIDX,SCOUNTS2,SDISPLS2,MPI_INTEGER,RECVIDX,RCOUNTS2,RDISPLS2,MPI_INTEGER,MPI_COMM_WORLD,IERROR4)

  ! Finally, fill the local chunk of the SA matrix (block of rows)
  ! with the received values at their appropriate place.
  do I=1,NRECV
    ISUP = RECVIDX(2*I-1)
    JSUP = RECVIDX(2*I)
    SA(ISUP-ILO+1,JSUP) = RECVVAL(I)
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

! Avoid unused argument warnings
if (.false.) call UNUSED_INTEGER(iHi)

end subroutine MKSA_G3_MPP

#elif defined (NAGFOR)

! Some compilers do not like empty files
subroutine empty_mksa_g3_mpp()
end subroutine empty_mksa_g3_mpp

#endif
