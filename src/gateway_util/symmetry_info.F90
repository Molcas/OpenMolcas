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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

module Symmetry_Info

implicit none
private

public :: nIrrep, iOper, iChTbl, iChCar, Mul, iChBas, lIrrep, lBsFnc, SymLab, iSkip, Symmetry_Info_Set, Symmetry_Info_Dmp, &
          Symmetry_Info_Get, Symmetry_Info_Back, Symmetry_Info_Free, Symmetry_Info_Setup, VarR, VarT

#include "stdalloc.fh"
integer :: nIrrep = 1
integer :: iOper(0:7) = [0,0,0,0,0,0,0,0]
integer :: iChTbl(0:7,0:7) = reshape([0,0,0,0,0,0,0,0, &
                                      0,0,0,0,0,0,0,0, &
                                      0,0,0,0,0,0,0,0, &
                                      0,0,0,0,0,0,0,0, &
                                      0,0,0,0,0,0,0,0, &
                                      0,0,0,0,0,0,0,0, &
                                      0,0,0,0,0,0,0,0, &
                                      0,0,0,0,0,0,0,0],[8,8])
integer :: iChCar(3) = [0,0,0]
integer :: MxFnc
integer, parameter :: Mul(8,8) = reshape([1,2,3,4,5,6,7,8, &
                                          2,1,4,3,6,5,8,7, &
                                          3,4,1,2,7,8,5,6, &
                                          4,3,2,1,8,7,6,5, &
                                          5,6,7,8,1,2,3,4, &
                                          6,5,8,7,2,1,4,3, &
                                          7,8,5,6,3,4,1,2, &
                                          8,7,6,5,4,3,2,1],[8,8])
integer, allocatable :: iChBas(:)
character(LEN=3) :: lIrrep(0:7) = ['','','','','','','','']
character(LEN=80) :: lBsFnc(0:7) = ['','','','','','','','']
character(LEN=3) SymLab
integer :: iSkip(0:7) = [0,0,0,0,0,0,0,0]
logical :: VarR = .false., VarT = .false.

interface
  subroutine Abend()
  end subroutine Abend
  subroutine Put_iArray(Label,data,nData)
    character*(*) Label
    integer nData
    integer data(nData)
  end subroutine Put_iArray
  subroutine Get_iArray(Label,data,nData)
    character*(*) Label
    integer nData
    integer data(nData)
  end subroutine Get_iArray
  subroutine Qpg_iArray(Label,Found,nData)
    character*(*) Label
    logical Found
    integer nData
  end subroutine Qpg_iArray
end interface

!***********************************************************************
!***********************************************************************

contains

!***********************************************************************
!***********************************************************************

! temporary routine!
subroutine Symmetry_Info_Back(mIrrep)

  integer :: mIrrep

  mIrrep = nIrrep
# ifdef _DEBUGPRINT_
  write(6,*) 'Call Symmetry_Info_Back'
  write(6,'(A,I4)') 'nIrrep=',nIrrep
# endif

end subroutine Symmetry_Info_Back

!***********************************************************************
!***********************************************************************

subroutine Symmetry_Info_Dmp()

  integer i, j, k, liDmp, lcDmp
  integer, allocatable :: iDmp(:)
  character(LEN=1), allocatable :: cDmp(:)

  liDmp = 1+8+8*8+3+MxFnc+8+2
  call mma_allocate(iDmp,liDmp,Label='iDmp')

  i = 0
  iDmp(i+1) = nIrrep
  i = i+1
  iDmp(i+1:i+8) = iOper(:)
  i = i+8
  iDmp(i+1:i+8) = iChTbl(:,0)
  i = i+8
  iDmp(i+1:i+8) = iChTbl(:,1)
  i = i+8
  iDmp(i+1:i+8) = iChTbl(:,2)
  i = i+8
  iDmp(i+1:i+8) = iChTbl(:,3)
  i = i+8
  iDmp(i+1:i+8) = iChTbl(:,4)
  i = i+8
  iDmp(i+1:i+8) = iChTbl(:,5)
  i = i+8
  iDmp(i+1:i+8) = iChTbl(:,6)
  i = i+8
  iDmp(i+1:i+8) = iChTbl(:,7)
  i = i+8
  iDmp(i+1:i+3) = iChCar(1:3)
  i = i+3
  iDmp(i+1:i+MxFnc) = iChBas(1:MxFnc)
  i = i+MxFnc
  iDmp(i+1:i+8) = iSkip(0:7)
  i = i+8
  iDmp(i+1:i+2) = merge([1,1],[0,0],[VarR,VarT])
  i = i+2

# ifdef _DEBUGPRINT_
  write(6,*) 'Symmetry_Info_Dmp'
  write(6,*) 'liDmp=',liDmp
  write(6,*) 'MxFnc=',MxFnc
  write(6,*) 'nIrrep=',nIrrep
  write(6,*) 'iOper:'
  write(6,'(8I4)') (iOper(i),i=0,nIrrep-1)
  write(6,*)
  write(6,'(9I4)') (iDmp(i),i=1,liDmp)
  write(6,*)
  write(6,*) 'lIrrep:'
  do i=0,nIrrep-1
    write(6,'(A)') lIrrep(i)
  end do
  write(6,*) 'lBsFnc:'
  do i=0,nIrrep-1
    write(6,'(A)') lBsFnc(i)
  end do
  write(6,'(2A)') 'SymLab:',SymLab
# endif

  call Put_iArray('Symmetry Info',iDmp,liDmp)
  call mma_deallocate(iDmp)

  lcDmp = 3*8+80*8+3
  call mma_allocate(cDmp,lcDmp,Label='cDmp')
  k = 0
  do i=0,7
    do j=1,3
      cDmp(j+k) = lIrrep(i)(j:j)
    end do
    k = k+3
  end do
  do i=0,7
    do j=1,80
      cDmp(j+k) = lBsFnc(i)(j:j)
    end do
    k = k+80
  end do
  do i=1,3
    cDmp(i+k) = SymLab(i:i)
  end do
  k = k+3
  call put_cArray('SymmetryCInfo',cDmp(1),lcDmp)
  call mma_deallocate(cDmp)

end subroutine Symmetry_Info_Dmp

!***********************************************************************
!***********************************************************************

subroutine Symmetry_Info_Get()

  integer i, j, k, liDmp, lcDmp
  integer, allocatable :: iDmp(:)
  logical Found
  character(LEN=1), allocatable :: cDmp(:)

  if (allocated(iChBas)) return
  call Qpg_iArray('Symmetry Info',Found,liDmp)
  call mma_allocate(iDmp,liDmp,Label='iDmp')
  call Get_iArray('Symmetry Info',iDmp,liDmp)

  MxFnc = liDmp-(1+8+8*8+3+8+2)
  call mma_allocate(iChBas,MxFnc,Label='iChBas')

  i = 0
  nIrrep = iDmp(i+1)
  i = i+1
  iOper(:) = iDmp(i+1:i+8)
  i = i+8
  iChTbl(:,0) = iDmp(i+1:i+8)
  i = i+8
  iChTbl(:,1) = iDmp(i+1:i+8)
  i = i+8
  iChTbl(:,2) = iDmp(i+1:i+8)
  i = i+8
  iChTbl(:,3) = iDmp(i+1:i+8)
  i = i+8
  iChTbl(:,4) = iDmp(i+1:i+8)
  i = i+8
  iChTbl(:,5) = iDmp(i+1:i+8)
  i = i+8
  iChTbl(:,6) = iDmp(i+1:i+8)
  i = i+8
  iChTbl(:,7) = iDmp(i+1:i+8)
  i = i+8
  iChCar(1:3) = iDmp(i+1:i+3)
  i = i+3
  iChBas(1:MxFnc) = iDmp(i+1:i+MxFnc)
  i = i+MxFnc
  iSKip(0:7) = iDmp(i+1:i+8)
  i = i+8
  VarR = iDmp(i+1) /= 0
  VarT = iDmp(i+2) /= 0
  i = i+2
  call mma_deallocate(iDmp)

  lcDmp = 3*8+80*8+3
  call mma_allocate(cDmp,lcDmp,Label='cDmp')
  call get_carray('SymmetryCInfo',cDmp(1),lcDmp)
  k = 0
  do i=0,7
    do j=1,3
      lIrrep(i)(j:j) = cDmp(j+k)
    end do
    k = k+3
  end do
  do i=0,7
    do j=1,80
      lBsFnc(i)(j:j) = cDmp(j+k)
    end do
    k = k+80
  end do
  do i=1,3
    SymLab(i:i) = cDmp(i+k)
  end do
  call mma_deallocate(cDmp)
# ifdef _DEBUGPRINT_
  write(6,*) 'Symmetry_Info_Get'
  write(6,*) 'liDmp=',liDmp
  write(6,*) 'MxFnc=',MxFnc
  write(6,*) 'nIrrep=',nIrrep
  write(6,*)
  write(6,'(2A)') 'SymLab:',SymLab
  write(6,*) 'iOper:'
  write(6,'(8I4)') (iOper(i),i=0,nIrrep-1)
  write(6,*)
  write(6,*) 'lIrrep:'
  do i=0,nIrrep-1
    write(6,'(A)') lIrrep(i)
  end do
  do i=0,nIrrep-1
    write(6,'(A)') lBsFnc(i)
  end do
# endif

end subroutine Symmetry_Info_Get

!***********************************************************************
!***********************************************************************

subroutine Symmetry_Info_Free()

  if (.not. allocated(iChBas)) return
  call mma_deallocate(iChBas)
  MxFnc = 0

end subroutine Symmetry_Info_Free

!***********************************************************************
!***********************************************************************

subroutine Symmetry_Info_Setup(nOper,Oper,iAng)

  implicit none
  integer :: nOper, iAng, i, j
  character(LEN=3) :: Oper(3)

  if (allocated(iChBas)) return   ! Return if already initiated.

  nIrrep = 2**nOper
  call Put_iScalar('NSYM',nIrrep)

  iOper(0) = 0
  do i=1,nOper
    iOper(i) = 0
    do j=1,3
      if (Oper(i)(j:j) == 'X') iOper(i) = iOper(i)+1
      if (Oper(i)(j:j) == 'Y') iOper(i) = iOper(i)+2
      if (Oper(i)(j:j) == 'Z') iOper(i) = iOper(i)+4
    end do
    if (iOper(i) == 0) then
      call WarningMessage(2,'RdCtl: Illegal symmetry operator!')
      write(6,*) 'Oper=',Oper(i)
      write(6,*)
      call Abend()
    end if
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Generate all operations of the group

  if (nOper >= 2) then
    iOper(4) = iOper(3)
    iOper(3) = ieor(iOper(1),iOper(2))
  end if
  if (nOper == 3) then
    iOper(5) = ieor(iOper(1),iOper(4))
    iOper(6) = ieor(iOper(2),iOper(4))
    iOper(7) = ieor(iOper(1),ieor(iOper(2),iOper(4)))
  end if

  call Put_iArray('Symmetry operations',iOper,nIrrep)

  ! Generate the Character table for all Irreps, iChTbl

  ! All Irreps are one dimensional, i.e. the Character for the
  ! unit operator is 1 in all irreps.
  ! The totally symmetric representation will have the character
  ! of 1 for any given operation
  ! Now, the Irreps are due to classes of operations and will
  ! present the character of this class. In case of Abelian groups
  ! or other one dimensional groups the classes will have one
  ! and only one operation. Hence, the operations themselves can
  ! be used to present the character of the Irreps.

  ! Generate iChTbl, lIrrep, lBsFnc (iSigma)
  call ChTab(iOper,nIrrep,iChTbl)

  ! Generate iChCar, iChBas, MxFnc
  call Symmetry_Info_Set(iAng)

end subroutine Symmetry_Info_Setup

!***********************************************************************
!***********************************************************************

subroutine Symmetry_Info_Set(iAng)

  integer :: iIrrep, jIrrep
  integer :: iSymX, iSymY, iSymZ, i
  integer :: iAng, lxyz, ixyz, ix, jx, iyMax, iy, jy, iz, jz, jxyz

  if (allocated(iChBas)) return

  ! Setup characteristics for cartesian basis functions.
  ! Observe that this is affected by the defined generators.
  ! In the array we will set the bit corresponding to a symop
  ! if that symop will alter the sign of the basis function.

  iSymX = 0
  iSymY = 0
  iSymZ = 0
  do i=0,nIrrep-1
    if (iand(iOper(i),1) /= 0) iSymX = 1
    if (iand(iOper(i),2) /= 0) iSymY = 2
    if (iand(iOper(i),4) /= 0) iSymZ = 4
  end do
  iChCar(1) = iSymX
  iChCar(2) = iSymY
  iChCar(3) = iSymZ

  MxFnc = (iAng+1)*(iAng+2)*(iAng+3)/6
  call mma_allocate(iChBas,MxFnc,Label='iChBas')
# ifdef _DEBUGPRINT_
  write(6,*) 'Symmetry_Info_Set:'
  write(6,*) 'iAng,MxFnc=',iAng,MxFnc
# endif

  lxyz = 0
  do ixyz=0,iAng
    do ix=ixyz,0,-1
      jx = mod(ix,2)
      iyMax = ixyz-ix
      do iy=iyMax,0,-1
        jy = mod(iy,2)
        lxyz = lxyz+1
        iz = ixyz-ix-iy
        jz = mod(iz,2)
        jxyz = jx*iSymX+jy*iSymY+jz*iSymZ
        iChBas(lxyz) = jxyz
      end do
    end do
  end do

  do iIrrep=0,nIrrep-2
    do jIrrep=iIrrep+1,nIrrep-1
      if (iOper(iIrrep) == iOper(jIrrep)) then
        call WarningMessage(2,' The generators of the point group are over defined, correct input!;'// &
                            'Abend: correct symmetry specifications!')
        call Quit_OnUserError()
      end if
    end do
  end do
# ifdef _DEBUGPRINT_
  write(6,*) 'Symmetry_Info_Set:'
  write(6,*) 'MxFnc=',MxFnc
  write(6,*) 'nIrrep=',nIrrep
  write(6,'(A,8I4)') 'iOper:',iOper(0:nIrrep-1)
  write(6,*) 'iChTbl:'
  do i=0,nIrrep-1
    write(6,'(8I4)') iChTbl(0:nIrrep-1,i)
  end do
# endif

end subroutine Symmetry_Info_Set

!***********************************************************************
!***********************************************************************

subroutine ChTab(iOper,nIrrep,iChTbl)
  !*********************************************************************
  !                                                                    *
  ! Object: to generate the character table of a point group within    *
  !         D2h.                                                       *
  !                                                                    *
  !     Author: Roland Lindh, Dept. of Theoretical Chemistry,          *
  !             University of Lund, SWEDEN                             *
  !             September '91                                          *
  !*********************************************************************

  implicit none
  integer nIrrep
  integer iOper(nIrrep), iChTbl(1:8,1:8) ! ugly dimensions change to 0:7!
  integer iTest(8)
  integer :: iSigma = 1
  character(Len=80) Tmp
  logical Inv, Rot
  character(LEN=6) :: xyz(0:7) = ['      ','x     ','y     ','xy, Rz','z     ','xz, Ry','yz, Rx','I     ']
  integer i, i1, i2, ia, ib, iCh, iFnc, iIrrep, iRot, iSub, iSymX, iSymY, iSymZ, ix, iy, iz, jIrrep
  integer j, jx, jy, jz, Lenlbs, LenlIrr, LenTmp
  integer iclast
  external iclast

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (nIrrep == 1) then
    SymLab = 'C1 '
    iSigma = 1
  else if (nIrrep == 2) then
    if (iOper(2) == 7) then
      SymLab = 'Ci '
      iSigma = 1
    else if ((iOper(2) == 1) .or. (iOper(2) == 2) .or. (iOper(2) == 4)) then
      SymLab = 'Cs'
      iSigma = 1
    else
      SymLab = 'C2'
      iSigma = 2
    end if
  else if (nIrrep == 4) then
    if ((iOper(2) == 7) .or. (iOper(3) == 7) .or. (iOper(4) == 7)) then
      SymLab = 'C2h'
      iSigma = 2
    else
      Rot = .true.
      do i=1,nIrrep
        if ((iOper(i) == 1) .or. (iOper(i) == 2) .or. (iOper(i) == 4)) Rot = .false.
      end do
      if (Rot) then
        SymLab = 'D2 '
        iSigma = 2
      else
        SymLab = 'C2v'
        iSigma = 2
      end if
    end if
  else if (nIrrep == 8) then
    SymLab = 'D2h'
    iSigma = 2
  else
    call WarningMessage(2,'ChTab: Illegal value of nIrrep')
    write(6,*) 'nIrrep=',nIrrep
    call Abend()
  end if
  ichTbl(:,:) = 0

  ! Go through the functions x, y, and z, and the dyadic functions.

  iSymX = 0
  iSymY = 0
  iSymZ = 0
  do i=1,nIrrep
    if (iand(iOper(i),1) /= 0) iSymX = 1
    if (iand(iOper(i),2) /= 0) iSymY = 2
    if (iand(iOper(i),4) /= 0) iSymZ = 4
  end do

  ! Loop over basis functions (a' la Malmqvist)

  lBsFnc(0:nIrrep-1) = '' ! For this to work we need a clean slate.
  do iFnc=0,7
    Tmp = xyz(iFnc)

    ! Generate a row in the character table of this function

    ix = iand(iFnc,iSymX)
    iy = iand(iFnc,iSymY)/2
    iz = iand(iFnc,iSymZ)/4
    ! Loop over all operators
    do i=1,nIrrep
      jx = iand(iOper(i),iSymX)
      jy = iand(iOper(i),iSymY)/2
      jz = iand(iOper(i),iSymZ)/4
      iCh = 1
      if ((ix /= 0) .and. (jx /= 0)) iCh = -iCh
      if ((iy /= 0) .and. (jy /= 0)) iCh = -iCh
      if ((iz /= 0) .and. (jz /= 0)) iCh = -iCh
      iTest(i) = iCh
    end do

    ! Compute place of Irrep

    if (nIrrep == 1) then
      jIrrep = 1
    else if (nIrrep == 2) then
      jIrrep = 1+(1-iTest(2))/2
    else if (nIrrep == 4) then
      jIrrep = 1+((1-iTest(2))+2*(1-iTest(3)))/2
    else if (nIrrep == 8) then
      jIrrep = 1+((1-iTest(2))+2*(1-iTest(3))+4*(1-iTest(5)))/2
    else
      jIrrep = -1
      call WarningMessage(2,'ChTab: Illegal nIrrep value!')
      write(6,*) 'nIrrep=',nIrrep
      call Abend()
    end if
    if (lBsFnc(jIrrep-1)(1:1) == ' ') then
      lBsFnc(jIrrep-1) = Tmp
      call ICopy(nIrrep,iTest,1,iChTbl(jIrrep,1),8)
    else
      LenlBs = len(lBsFnc(jIrrep-1))
      LenTmp = len(Tmp)
      i1 = iCLast(lBsFnc(jIrrep-1),LenlBs)
      i2 = iCLast(Tmp,LenTmp)
      lBsFnc(jIrrep-1) = lBsFnc(jIrrep-1)(1:i1)//', '//Tmp(1:i2)
    end if
  end do

  ! Set up some Mulliken symbols for the irreps

  do iIrrep=1,nIrrep
    lIrrep(iIrrep-1) = 'a'
    do i=1,nIrrep

      ! If the character of an rotation in an irreps is -1 then
      ! the irreps is assigned the character B, otherwise A.

      if (((iOper(i) == 3) .or. (iOper(i) == 5) .or. (iOper(i) == 6)) .and. (iChTbl(iIrrep,i) == -1)) lIrrep(iIrrep-1) = 'b'

    end do
  end do
  iSub = 0

  ! Subscript according to C2 operations

  Rot = .false.
  do i=1,nIrrep
    if ((iOper(i) == 3) .or. (iOper(i) == 5) .or. (iOper(i) == 6)) Rot = .true.
  end do
  if (Rot .and. (SymLab /= 'C2v')) then
    iSub = iSub+1

    ! Find the number of A's and B's

    ia = 0
    ib = 0
    do i=1,nIrrep
      if (lIrrep(i-1)(1:1) == 'a') ia = ia+1
      if (lIrrep(i-1)(1:1) == 'b') ib = ib+1
    end do
    if (nIrrep == 8) then
      ia = ia/2
      ib = ib/2
    end if
    if (SymLab == 'C2h') then
      ia = ia/2
      ib = ib/2
    end if

    ! Find the rotations

    iRot = 0
    do i=1,nIrrep
      if ((iOper(i) == 3) .or. (iOper(i) == 5) .or. (iOper(i) == 6)) then
        iRot = iRot+1
        write(Tmp,'(I1)') iRot
        if (ia > 1) then
          do j=1,nIrrep
            if ((lIrrep(j-1)(1:1) == 'a') .and. (iChTbl(j,i) == 1)) lIrrep(j-1) = lIrrep(j-1)(1:1)//Tmp(1:1)
          end do
        end if
        if (ib > 1) then
          do j=1,nIrrep
            if ((lIrrep(j-1)(1:1) == 'b') .and. (iChTbl(j,i) == 1)) lIrrep(j-1) = lIrrep(j-1)(1:1)//Tmp(1:1)
          end do
        end if
      end if
    end do
  else if (Rot .and. (SymLab == 'C2v')) then

    ! Find the Rotation

    iRot = -1
    do i=1,nIrrep
      if ((iOper(i) == 3) .or. (iOper(i) == 5) .or. (iOper(i) == 6)) iRot = iOper(i)
    end do

    ! Find the first vertical mirror plane to this axis

    do i=1,nIrrep
      if ((iOper(i) /= 3) .and. (iOper(i) /= 5) .and. (iOper(i) /= 6) .and. (iOper(i) /= 7) .and. (iand(iOper(i),iRot) /= 1)) &
        iRot = i
    end do
    do i=1,nIrrep
      if (iChTbl(i,iRot) == 1) then
        j = 1
      else
        j = 2
      end if
      write(Tmp,'(I1)') j
      lIrrep(i-1) = lIrrep(i-1)(1:1)//Tmp(1:1)
    end do
  end if

  ! Subscript according to inversion if present

  Inv = .false.
  do i=1,nIrrep
    Inv = (iOper(i) == 7) .or. Inv
  end do
  if (Inv) then
    iSub = iSub+1

    ! Loop over each Irrep

    do iIrrep=1,nIrrep
      LenlIrr = len(lIrrep(iIrrep-1))
      i1 = 1+iCLast(lIrrep(iIrrep-1),LenlIrr)

      ! Loop over operators

      do i=1,nIrrep
        if (iOper(i) == 7) then
#         ifdef _WARNING_WORKAROUND_
          ! see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=101827
          i1 = min(i1,len(lIrrep))
#         endif
          if (iChTbl(iIrrep,i) == 1) then
            lIrrep(iIrrep-1)(i1:i1) = 'g'
          else if (iChTbl(iIrrep,i) == -1) then
            lIrrep(iIrrep-1)(i1:i1) = 'u'
          end if
        end if
      end do
    end do
  end if

  ! Fix labels for Cs

  if (SymLab(1:2) == 'Cs') then
    lIrrep(0) = "a'"
    lIrrep(1) = 'a"'
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call Put_iScalar('Rotational Symmetry Number',iSigma)
  call Put_cArray('Irreps',lIrrep(0),24)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  return

end subroutine ChTab

!***********************************************************************
!***********************************************************************

end module Symmetry_Info
