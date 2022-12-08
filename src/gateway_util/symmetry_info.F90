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

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: iChCar(3) = 0, iChTbl(0:7,0:7) = 0, iOper(0:7) = 0, iSkip(0:7) = 0, MxFnc, nIrrep = 1
logical(kind=iwp) :: VarR = .false., VarT = .false.
character(len=80) :: lBsFnc(0:7) = ''
character(len=3) :: lIrrep(0:7) = '', SymLab
integer(kind=iwp), allocatable :: iChBas(:)
integer(kind=iwp), parameter :: Mul(8,8) = reshape([1,2,3,4,5,6,7,8, &
                                                    2,1,4,3,6,5,8,7, &
                                                    3,4,1,2,7,8,5,6, &
                                                    4,3,2,1,8,7,6,5, &
                                                    5,6,7,8,1,2,3,4, &
                                                    6,5,8,7,2,1,4,3, &
                                                    7,8,5,6,3,4,1,2, &
                                                    8,7,6,5,4,3,2,1],[8,8]), &
                                Prmt(0:7,0:7) = reshape([1, 1, 1, 1, 1, 1, 1, 1, &
                                                         1,-1, 1,-1, 1,-1, 1,-1, &
                                                         1, 1,-1,-1, 1, 1,-1,-1, &
                                                         1,-1,-1, 1, 1,-1,-1, 1, &
                                                         1, 1, 1, 1,-1,-1,-1,-1, &
                                                         1,-1, 1,-1,-1, 1,-1, 1, &
                                                         1, 1,-1,-1,-1,-1, 1, 1, &
                                                         1,-1,-1, 1,-1, 1, 1,-1],[8,8])

public :: iChBas, iChCar, iChTbl, iOper, iSkip, lBsFnc, lIrrep, Mul, nIrrep, Prmt, SymLab, Symmetry_Info_Dmp, Symmetry_Info_Free, &
          Symmetry_Info_Get, Symmetry_Info_Set, Symmetry_Info_Setup, VarR, VarT

!***********************************************************************
!***********************************************************************

contains

!***********************************************************************
!***********************************************************************

subroutine Symmetry_Info_Dmp()

  use fortran_strings, only: char_array
  use stdalloc, only: mma_allocate, mma_deallocate
# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  integer(kind=iwp) :: i, j, lcDmp, liDmp
  integer(kind=iwp), allocatable :: iDmp(:)
  character, allocatable :: cDmp(:)

  if (.not. allocated(iChBas)) then
    call WarningMessage(2,'Symmetry_Info_Dmp: iChBas is not allocated!')
    call Abend()
  end if

  liDmp = 1+8+8*8+3+MxFnc+8+2
  call mma_allocate(iDmp,liDmp,Label='iDmp')

  i = 0
  iDmp(i+1) = nIrrep
  i = i+1
  iDmp(i+1:i+8) = iOper(:)
  i = i+8
  iDmp(i+1:i+8*8) = reshape(iChTbl(:,:),[8*8])
  i = i+8*8
  iDmp(i+1:i+3) = iChCar(1:3)
  i = i+3
  iDmp(i+1:i+MxFnc) = iChBas(1:MxFnc)
  i = i+MxFnc
  iDmp(i+1:i+8) = iSkip(0:7)
  i = i+8
  iDmp(i+1) = merge(1,0,VarR)
  i = i+1
  iDmp(i+1) = merge(1,0,VarT)
  i = i+1

# ifdef _DEBUGPRINT_
  write(u6,*) 'Symmetry_Info_Dmp'
  write(u6,*) 'liDmp=',liDmp
  write(u6,*) 'MxFnc=',MxFnc
  write(u6,*) 'nIrrep=',nIrrep
  write(u6,*) 'iOper:'
  write(u6,'(8I4)') (iOper(i),i=0,nIrrep-1)
  write(u6,*)
  write(u6,'(9I4)') (iDmp(i),i=1,liDmp)
  write(u6,*)
  write(u6,*) 'lIrrep:'
  do i=0,nIrrep-1
    write(u6,'(A)') lIrrep(i)
  end do
  write(u6,*) 'lBsFnc:'
  do i=0,nIrrep-1
    write(u6,'(A)') lBsFnc(i)
  end do
  write(u6,'(2A)') 'SymLab:',SymLab
# endif

  call Put_iArray('Symmetry Info',iDmp,liDmp)
  call mma_deallocate(iDmp)

  lcDmp = 3*8+80*8+3
  call mma_allocate(cDmp,lcDmp,Label='cDmp')
  j = 0
  do i=0,7
    cDmp(j+1:j+3) = char_array(lIrrep(i))
    j = j+3
  end do
  do i=0,7
    cDmp(j+1:j+80) = char_array(lBsFnc(i))
    j = j+80
  end do
  cDmp(j+1:j+3) = char_array(SymLab)
  j = j+3
  call put_cArray('SymmetryCInfo',cDmp(1),lcDmp)
  call mma_deallocate(cDmp)

end subroutine Symmetry_Info_Dmp

!***********************************************************************
!***********************************************************************

subroutine Symmetry_Info_Get()

  use fortran_strings, only: str
  use stdalloc, only: mma_allocate, mma_deallocate
# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  integer(kind=iwp) :: i, j, liDmp, lcDmp
  logical(kind=iwp) :: Found
  integer(kind=iwp), allocatable :: iDmp(:)
  character, allocatable :: cDmp(:)

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
  iChTbl(:,:) = reshape(iDmp(i+1:i+8*8),[8,8])
  i = i+8*8
  iChCar(1:3) = iDmp(i+1:i+3)
  i = i+3
  iChBas(1:MxFnc) = iDmp(i+1:i+MxFnc)
  i = i+MxFnc
  iSKip(0:7) = iDmp(i+1:i+8)
  i = i+8
  VarR = iDmp(i+1) > 0
  i = i+1
  VarT = iDmp(i+1) > 0
  i = i+1
  call mma_deallocate(iDmp)

  lcDmp = 3*8+80*8+3
  call mma_allocate(cDmp,lcDmp,Label='cDmp')
  call get_carray('SymmetryCInfo',cDmp(1),lcDmp)
  j = 0
  do i=0,7
    lIrrep(i) = str(cDmp(j+1:j+3))
    j = j+3
  end do
  do i=0,7
    lBsFnc(i) = str(cDmp(j+1:j+80))
    j = j+80
  end do
  SymLab = str(cDmp(j+1:j+3))
  j = j+3
  call mma_deallocate(cDmp)
# ifdef _DEBUGPRINT_
  write(u6,*) 'Symmetry_Info_Get'
  write(u6,*) 'liDmp=',liDmp
  write(u6,*) 'MxFnc=',MxFnc
  write(u6,*) 'nIrrep=',nIrrep
  write(u6,*)
  write(u6,'(2A)') 'SymLab:',SymLab
  write(u6,*) 'iOper:'
  write(u6,'(8I4)') (iOper(i),i=0,nIrrep-1)
  write(u6,*)
  write(u6,*) 'lIrrep:'
  do i=0,nIrrep-1
    write(u6,'(A)') lIrrep(i)
  end do
  do i=0,nIrrep-1
    write(u6,'(A)') lBsFnc(i)
  end do
# endif

end subroutine Symmetry_Info_Get

!***********************************************************************
!***********************************************************************

subroutine Symmetry_Info_Free()

  use stdalloc, only: mma_deallocate

  if (.not. allocated(iChBas)) return
  call mma_deallocate(iChBas)
  MxFnc = 0

end subroutine Symmetry_Info_Free

!***********************************************************************
!***********************************************************************

subroutine Symmetry_Info_Setup(nOper,Oper,iAng)

  use Definitions, only: u6

  integer(kind=iwp), intent(in) :: nOper, iAng
  character(len=3), intent(in) :: Oper(3)
  integer(kind=iwp) :: i, j

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
      write(u6,*) 'Oper=',Oper(i)
      write(u6,*)
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

  use stdalloc, only: mma_allocate
# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  integer(kind=iwp), intent(in) :: iAng
  integer(kind=iwp) :: i, iIrrep, iSymX, iSymY, iSymZ, ix, ixyz, iy, iyMax, iz, jIrrep, jx, jxyz, jy, jz, lxyz

  if (allocated(iChBas)) return

  ! Setup characteristics for cartesian basis functions.
  ! Observe that this is affected by the defined generators.
  ! In the array we will set the bit corresponding to a symop
  ! if that symop will alter the sign of the basis function.

  iSymX = 0
  iSymY = 0
  iSymZ = 0
  do i=0,nIrrep-1
    if (btest(iOper(i),0)) iSymX = 2**0
    if (btest(iOper(i),1)) iSymY = 2**1
    if (btest(iOper(i),2)) iSymZ = 2**2
  end do
  iChCar(1) = iSymX
  iChCar(2) = iSymY
  iChCar(3) = iSymZ

  MxFnc = (iAng+1)*(iAng+2)*(iAng+3)/6
  call mma_allocate(iChBas,MxFnc,Label='iChBas')
# ifdef _DEBUGPRINT_
  write(u6,*) 'Symmetry_Info_Set:'
  write(u6,*) 'iAng,MxFnc=',iAng,MxFnc
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
  write(u6,*) 'Symmetry_Info_Set:'
  write(u6,*) 'MxFnc=',MxFnc
  write(u6,*) 'nIrrep=',nIrrep
  write(u6,'(A,8I4)') 'iOper:',iOper(0:nIrrep-1)
  write(u6,*) 'iChTbl:'
  do i=0,nIrrep-1
    write(u6,'(8I4)') iChTbl(0:nIrrep-1,i)
  end do
# endif

end subroutine Symmetry_Info_Set

!***********************************************************************
!***********************************************************************

subroutine ChTab(iOper,nIrrep,outChTbl)
  !*********************************************************************
  !                                                                    *
  ! Object: to generate the character table of a point group within    *
  !         D2h.                                                       *
  !                                                                    *
  !     Author: Roland Lindh, Dept. of Theoretical Chemistry,          *
  !             University of Lund, SWEDEN                             *
  !             September '91                                          *
  !*********************************************************************

  use Definitions, only: u6

  integer(kind=iwp), intent(in) :: nIrrep, iOper(nIrrep)
  integer(kind=iwp), intent(out) :: outChTbl(1:8,1:8) ! ugly dimensions change to 0:7!
  integer(kind=iwp) :: i, i1, ia, ib, iCh, iFnc, iIrrep, iRot, iSigma = 1, iSub, iTest(8), ix, iy, iz, j, jIrrep, jx, jy, jz
  logical(kind=iwp) :: Inv, Rot, SymX, SymY, SymZ
  character(len=80) :: Tmp
  character(len=*), parameter :: xyz(0:7) = ['      ', &
                                             'x     ', &
                                             'y     ', &
                                             'xy, Rz', &
                                             'z     ', &
                                             'xz, Ry', &
                                             'yz, Rx', &
                                             'I     ']

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
    write(u6,*) 'nIrrep=',nIrrep
    call Abend()
  end if
  outChTbl(:,:) = 0

  ! Go through the functions x, y, and z, and the dyadic functions.

  SymX = .false.
  SymY = .false.
  SymZ = .false.
  do i=1,nIrrep
    SymX = SymX .or. btest(iOper(i),0)
    SymY = SymY .or. btest(iOper(i),1)
    SymZ = SymZ .or. btest(iOper(i),2)
  end do

  ! Loop over basis functions (a' la Malmqvist)

  lBsFnc(0:nIrrep-1) = '' ! For this to work we need a clean slate.
  do iFnc=0,7
    Tmp = xyz(iFnc)

    ! Generate a row in the character table of this function

    ix = merge(1,0,SymX .and. btest(iFnc,0))
    iy = merge(1,0,SymY .and. btest(iFnc,1))
    iz = merge(1,0,SymZ .and. btest(iFnc,2))
    ! Loop over all operators
    do i=1,nIrrep
      jx = merge(1,0,SymX .and. btest(iOper(i),0))
      jy = merge(1,0,SymY .and. btest(iOper(i),1))
      jz = merge(1,0,SymZ .and. btest(iOper(i),2))
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
      write(u6,*) 'nIrrep=',nIrrep
      call Abend()
    end if
    if (lBsFnc(jIrrep-1)(1:1) == ' ') then
      lBsFnc(jIrrep-1) = Tmp
      outChTbl(jIrrep,1:nIrrep) = iTest(1:nIrrep)
    else
      lBsFnc(jIrrep-1) = trim(lBsFnc(jIrrep-1))//', '//trim(Tmp)
    end if
  end do

  ! Set up some Mulliken symbols for the irreps

  do iIrrep=1,nIrrep
    lIrrep(iIrrep-1) = 'a'
    do i=1,nIrrep

      ! If the character of an rotation in an irreps is -1 then
      ! the irreps is assigned the character B, otherwise A.

      if (((iOper(i) == 3) .or. (iOper(i) == 5) .or. (iOper(i) == 6)) .and. (outChTbl(iIrrep,i) == -1)) lIrrep(iIrrep-1) = 'b'

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
            if ((lIrrep(j-1)(1:1) == 'a') .and. (outChTbl(j,i) == 1)) lIrrep(j-1) = lIrrep(j-1)(1:1)//Tmp(1:1)
          end do
        end if
        if (ib > 1) then
          do j=1,nIrrep
            if ((lIrrep(j-1)(1:1) == 'b') .and. (outChTbl(j,i) == 1)) lIrrep(j-1) = lIrrep(j-1)(1:1)//Tmp(1:1)
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
      if (outChTbl(i,iRot) == 1) then
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
      i1 = 1+len_trim(lIrrep(iIrrep-1))

      ! Loop over operators

      do i=1,nIrrep
        if (iOper(i) == 7) then
#         ifdef _WARNING_WORKAROUND_
          ! see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=101827
          i1 = min(i1,len(lIrrep))
#         endif
          if (outChTbl(iIrrep,i) == 1) then
            lIrrep(iIrrep-1)(i1:i1) = 'g'
          else if (outChTbl(iIrrep,i) == -1) then
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
