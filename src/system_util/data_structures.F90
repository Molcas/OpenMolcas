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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************

module Data_Structures

use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
private

public :: Alloc1DiArray_Type, Alloc1DArray_Type, Alloc2DArray_Type, Alloc4DArray_Type, Allocate_DT, Deallocate_DT, DSBA_Type, &
          G2_Type, Integer_Pointer, NDSBA_Type, SBA_Type, twxy_Type, V1, V2

type Integer_Pointer
  integer(kind=iwp), contiguous, pointer :: I1(:) => null()
end type Integer_Pointer

type SB_Type
  real(kind=wp), contiguous, pointer :: A3(:,:,:) => null()
  real(kind=wp), contiguous, pointer :: A2(:,:) => null()
  real(kind=wp), contiguous, pointer :: A1(:) => null()
end type SB_Type

type DSB_Type
  real(kind=wp), contiguous, pointer :: A2(:,:) => null()
  real(kind=wp), contiguous, pointer :: A1(:) => null()
end type DSB_Type

type V1
  real(kind=wp), contiguous, pointer :: A(:) => null()
end type V1

type V2
  real(kind=wp), contiguous, pointer :: A(:,:) => null()
end type V2

type G2_pointers
  real(kind=wp), contiguous, pointer :: A4(:,:,:,:) => null()
  real(kind=wp), contiguous, pointer :: A2(:,:) => null()
end type G2_pointers

type NDSBA_Type
  integer(kind=iwp) :: iCase = 0
  integer(kind=iwp) :: nSym = 0
  real(kind=wp), allocatable :: A0(:)
  type(DSB_Type) :: SB(8,8)
end type NDSBA_Type

type DSBA_Type
  integer(kind=iwp) :: iCase = 0
  integer(kind=iwp) :: nSym = 0
  logical(kind=iwp) :: Fake = .false.
  logical(kind=iwp) :: Active = .false.
  real(kind=wp), allocatable :: A00(:)
  real(kind=wp), contiguous, pointer :: A0(:) => null()
  type(DSB_Type) :: SB(8)
end type DSBA_Type

type SBA_Type
  integer(kind=iwp) :: iCase = 0
  integer(kind=iwp) :: iSym = 0
  integer(kind=iwp) :: nSym = 0
  integer(kind=iwp) :: ipOff(8) = 0
  real(kind=wp), allocatable :: A0(:)
  type(SB_Type) :: SB(8)
end type SBA_Type

type twxy_type
  integer(kind=iwp) :: iCase = 0
  integer(kind=iwp) :: JSYM = 0
  integer(kind=iwp) :: nSym = 0
  real(kind=wp), allocatable :: twxy_full(:)
  type(V2) :: SB(8,8)
end type twxy_type

type G2_type
  integer(kind=iwp) :: iCase = 0
  integer(kind=iwp) :: nSym = 0
  real(kind=wp), allocatable :: A0(:)
  type(G2_Pointers) :: SB(8,8,8)
end type G2_type

type Alloc1DArray_Type
  real(kind=wp), allocatable :: A(:)
end type Alloc1DArray_Type

type Alloc2DArray_Type
  real(kind=wp), allocatable :: A(:,:)
end type Alloc2DArray_Type

type Alloc4DArray_Type
  real(kind=wp), allocatable :: A(:,:,:,:)
end type Alloc4DArray_Type

type Alloc1DiArray_Type
  integer(kind=iwp), allocatable :: A(:)
end type Alloc1DiArray_Type

! Allocate/deallocate data types
interface Allocate_DT
  module procedure :: Allocate_DSBA, Allocate_SBA, Allocate_twxy, Allocate_NDSBA, Allocate_G2, Alloc_Alloc_DSBA, &
                      Alloc_Alloc1DArray, Alloc2D_Alloc1DArray, Alloc_Alloc2DArray, Alloc_Alloc1DiArray
end interface Allocate_DT
interface Deallocate_DT
  module procedure :: Deallocate_DSBA, Deallocate_SBA, Deallocate_twxy, Deallocate_NDSBA, Deallocate_G2, Free_Alloc_DSBA, &
                      Free_Alloc1DArray, Free2D_Alloc1DArray, Free_Alloc2DArray, Free_Alloc1DiArray
end interface Deallocate_DT

! Private extensions to mma interfaces

interface mma_allocate
  module procedure :: dsba_mma_allo_1D, dsba_mma_allo_1D_lim, a1da_mma_allo_1D, a1da_mma_allo_1D_lim, a1da_mma_allo_2D, &
                      a1da_mma_allo_2D_lim, a2da_mma_allo_1D, a2da_mma_allo_1D_lim, a1ia_mma_allo_1D, a1ia_mma_allo_1D_lim
end interface
interface mma_deallocate
  module procedure :: dsba_mma_free_1D, a1da_mma_free_1D, a1da_mma_free_2D, a2da_mma_free_1D, a1ia_mma_free_1D
end interface

#include "compiler_features.h"

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                N D S B A - T Y P E   S E C T I O N                   !
!                                                                      !
!                Non-Diagonal Symmetry Blocked Arrays                  !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Allocate_NDSBA(Adam,n,m,nSym,Label)

  type(NDSBA_Type), target, intent(out) :: Adam
  integer(kind=iwp), intent(in) :: nSym, n(nSym), m(nSym)
  character(len=*), intent(in), optional :: Label
  integer(kind=iwp) :: iE, iS, iSym, jSym, MemTot

  Adam%iCase = 1
  Adam%nSym = nSym

  MemTot = 0
  do jSym=1,nSym
    do iSym=jSym,nSym
      MemTot = MemTot+n(iSym)*m(jSym)
    end do
  end do
  if (present(Label)) then
    call mma_allocate(Adam%A0,MemTot,Label=Label)
  else
    call mma_allocate(Adam%A0,MemTot,Label='%A0')
  end if

  iE = 0
  do jSym=1,nSym
    do iSym=jSym,nSym
      iS = iE+1
      iE = iE+n(iSym)*m(jSym)

      Adam%SB(jSym,iSym)%A2(1:n(iSym),1:m(jSym)) => Adam%A0(iS:iE)
      Adam%SB(jSym,iSym)%A1(1:n(iSym)*m(jSym)) => Adam%A0(iS:iE)
      Adam%SB(iSym,jSym)%A2(1:n(iSym),1:m(jSym)) => Adam%A0(iS:iE)
      Adam%SB(iSym,jSym)%A1(1:n(iSym)*m(jSym)) => Adam%A0(iS:iE)
    end do
  end do

end subroutine Allocate_NDSBA

subroutine Deallocate_NDSBA(Adam)

  type(NDSBA_Type), intent(inout) :: Adam
  integer(kind=iwp) :: iSym, jSym

  do iSym=1,Adam%nSym
    do jSym=1,Adam%nSym
      nullify(Adam%SB(iSym,jSym)%A2,Adam%SB(iSym,jSym)%A1)
    end do
  end do
  call mma_deallocate(Adam%A0)
  Adam%nSym = 0
  Adam%iCase = 0

end subroutine Deallocate_NDSBA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                D S B A - T Y P E   S E C T I O N                     !
!                                                                      !
!                Diagonal Symmetry Blocked Arrays                      !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Allocate_DSBA(Adam,n,m,nSym,aCase,Ref,Label)

  type(DSBA_Type), target, intent(out) :: Adam
  integer(kind=iwp), intent(in) :: nSym, n(nSym), m(nSym)
  character(len=3), intent(in), optional :: aCase
  real(kind=wp), target, intent(in), optional :: Ref(*)
  character(len=*), intent(in), optional :: Label
  integer(kind=iwp) :: iE, iS, iSym, MemTot, iCase = 0
  character(len=3) :: aCase_

  if (Adam%Active) then
    write(u6,*) 'DSBA-Type double allocate'
    call abend()
  end if

  if (present(aCase)) then
    aCase_ = aCase
  else
    aCase_ = 'REC'
  end if
  select case (aCase_)
    case ('TRI')
      iCase = 2
      do iSym=1,nSym
        if (n(iSym) /= m(iSym)) then
          write(u6,*) 'Allocate_DSBA: n(iSym)/=m(iSym), illegal if aCase="TRI".'
          call Abend()
        end if
      end do
    case ('REC')
      iCase = 1
    case ('ONE')
      iCase = 0
    case default
      write(u6,*) 'Allocate_DSBA: Illegal aCase parameter, aCase=',aCase_
      write(u6,*) 'Allowed value are "TRI", "REC", and "ONE".'
      call Abend()
  end select
  Adam%iCase = iCase
  Adam%nSym = nSym

  MemTot = 0
  select case (iCase)
    case (0)
      do iSym=1,nSym
        MemTot = MemTot+n(iSym)
      end do
    case (1)
      do iSym=1,nSym
        MemTot = MemTot+n(iSym)*m(iSym)
      end do
    case (2)
      do iSym=1,nSym
        MemTot = MemTot+n(iSym)*(n(iSym)+1)/2
      end do
  end select
  if (present(Ref)) then
    Adam%Fake = .true.
    Adam%A0(1:MemTot) => Ref(1:MemTot)
  else
    nullify(Adam%A0)
    if (present(Label)) then
      call mma_allocate(Adam%A00,MemTot,Label=Label)
    else
      call mma_allocate(Adam%A00,MemTot,Label='%A00')
    end if
    Adam%A0(1:MemTot) => Adam%A00(1:MemTot)
  end if

  Adam%Active = .true.
  iE = 0
  select case (iCase)
    case (0)
      do iSym=1,nSym
        iS = iE+1
        iE = iE+n(iSym)

        Adam%SB(iSym)%A1(1:n(iSym)) => Adam%A0(iS:iE)
      end do
    case (1)
      do iSym=1,nSym
        iS = iE+1
        iE = iE+n(iSym)*m(iSym)

        Adam%SB(iSym)%A2(1:n(iSym),1:m(iSym)) => Adam%A0(iS:iE)
        Adam%SB(iSym)%A1(1:n(iSym)*m(iSym)) => Adam%A0(iS:iE)
      end do
    case (2)
      do iSym=1,nSym
        iS = iE+1
        iE = iE+n(iSym)*(n(iSym)+1)/2

        Adam%SB(iSym)%A1(1:n(iSym)*(n(iSym)+1)/2) => Adam%A0(iS:iE)
      end do
  end select

end subroutine Allocate_DSBA

subroutine Deallocate_DSBA(Adam)

  type(DSBA_Type), intent(inout) :: Adam
  integer(kind=iwp) :: iSym

  if (.not. Adam%Active) return

  Adam%Active = .false.
  select case (Adam%iCase)
    case (0,2)
      do iSym=1,Adam%nSym
        nullify(Adam%SB(iSym)%A1)
      end do
    case (1)
      do iSym=1,Adam%nSym
        nullify(Adam%SB(iSym)%A2,Adam%SB(iSym)%A1)
      end do
  end select
  nullify(Adam%A0)
  if (Adam%Fake) then
    Adam%Fake = .false.
  else
    call mma_deallocate(Adam%A00)
  end if
  Adam%nSym = 0
  Adam%iCase = 0

end subroutine Deallocate_DSBA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                  S B A - T Y P E   S E C T I O N                     !
!                                                                      !
!                    Symmetry Blocked Arrays                           !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Allocate_SBA(Adam,n,m,NUMV,iSym,nSym,iCase,Memory,Label)

  type(SBA_Type), target, intent(out) :: Adam
  integer(kind=iwp), intent(in) :: nSym, n(nSym), m(nSym), NUMV, iSym, iCase
  integer(kind=iwp), intent(out), optional :: Memory
  character(len=*), intent(in), optional :: Label
  integer(kind=iwp) :: iE, iS, iSyma, iSymb, MemTot, n2Dim, n3Dim

  MemTot = 0

  select case (iCase)
    case (0)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        MemTot = MemTot+n(iSyma)*m(iSymb)*NUMV
      end do
    case (1)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        MemTot = MemTot+m(iSyma)*n(iSymb)*NUMV
      end do
    case (2)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        MemTot = MemTot+n(iSyma)*NUMV*m(iSymb)
      end do
    case (3)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        MemTot = MemTot+m(iSyma)*NUMV*n(iSymb)
      end do
    case (4)
      do iSyma=1,nSym
        if (n(iSyma) /= m(iSyma)) then
          write(u6,*) 'Allocate_SBA: iCase=4 only valid if n(:)=m(:).'
          call abend()
        end if
        iSymb = Mul(iSym,iSyma)
        if (iSyma == iSymb) then
          n2Dim = n(iSyma)*(n(iSyma)+1)/2
        else
          n2Dim = n(iSyma)*n(iSymb)
        end if
        MemTot = MemTot+n2dim*NUMV
      end do
    case (5)
      do iSyma=1,nSym
        if (n(iSyma) /= m(iSyma)) then
          write(u6,*) 'Allocate_SBA: iCase=5 only valid if n(:)=m(:).'
          call abend()
        end if
        iSymb = Mul(iSym,iSyma)
        if (iSyma == iSymb) then
          n2Dim = n(iSyma)*(n(iSyma)+1)/2
        else if (iSymb > iSyma) then
          n2Dim = n(iSyma)*n(iSymb)
        else
          n2Dim = 0
        end if
        MemTot = MemTot+n2dim*NUMV
      end do
    case (6)
      do iSyma=1,nSym
        if (n(iSyma) /= m(iSyma)) then
          write(u6,*) 'Allocate_SBA: iCase=6 only valid if n(:)=m(:).'
          call abend()
        end if
        iSymb = Mul(iSym,iSyma)
        if (iSyma >= iSymb) then
          n2Dim = n(iSyma)*n(iSymb)
        else
          n2Dim = 0
        end if
        MemTot = MemTot+n2dim*NUMV
      end do
    case default
      write(u6,*) 'Allocate_SBA: Illegal case.'
      call Abend()
  end select

  if (present(Memory)) then
    Memory = MemTot
    return
  end if

  Adam%iSym = iSym
  Adam%nSym = nSym
  Adam%iCase = iCase

  if (present(Label)) then
    call mma_allocate(Adam%A0,MemTot,Label=Label)
  else
    call mma_allocate(Adam%A0,MemTot,Label='%A0')
  end if

  iE = 0

  select case (iCase)
    case (0)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        iS = iE+1
        iE = iE+n(iSyma)*m(iSymb)*NUMV
        n2Dim = n(iSyma)*m(iSymb)
        n3Dim = n2Dim*NUMV
        Adam%SB(iSyma)%A1(1:n3Dim) => Adam%A0(iS:iE)
        Adam%SB(iSyma)%A2(1:n2Dim,1:NUMV) => Adam%A0(iS:iE)
        Adam%SB(iSyma)%A3(1:n(iSyma),1:m(iSymb),1:NUMV) => Adam%A0(iS:iE)
      end do
    case (1)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        iS = iE+1
        iE = iE+m(iSyma)*n(iSymb)*NUMV
        n2Dim = m(iSyma)*n(iSymb)
        n3Dim = n2Dim*NUMV
        Adam%SB(iSyma)%A1(1:n3Dim) => Adam%A0(iS:iE)
        Adam%SB(iSyma)%A2(1:n2Dim,1:NUMV) => Adam%A0(iS:iE)
        Adam%SB(iSyma)%A3(1:m(iSyma),1:n(iSymb),1:NUMV) => Adam%A0(iS:iE)
      end do
    case (2)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        iS = iE+1
        iE = iE+n(iSyma)*NUMV*m(iSymb)
        Adam%SB(iSyma)%A3(1:n(iSyma),1:NUMV,1:m(iSymb)) => Adam%A0(iS:iE)
      end do
    case (3)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        iS = iE+1
        iE = iE+m(iSyma)*NUMV*n(iSymb)
        Adam%SB(iSyma)%A3(1:m(iSyma),1:NUMV,1:n(iSymb)) => Adam%A0(iS:iE)
      end do
    case (4)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        iS = iE+1
        if (iSyma == iSymb) then
          n2Dim = n(iSyma)*(n(iSyma)+1)/2
        else
          n2Dim = n(iSyma)*n(iSymb)
        end if
        iE = iE+n2Dim*NUMV
        Adam%SB(iSymb)%A2(1:n2Dim,1:NUMV) => Adam%A0(iS:iE)
      end do
    case (5)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        if (iSymb > iSyma) cycle
        iS = iE+1
        if (iSyma == iSymb) then
          n2Dim = n(iSyma)*(n(iSyma)+1)/2
        else
          n2Dim = n(iSyma)*n(iSymb)
        end if
        iE = iE+n2Dim*NUMV
        Adam%SB(iSymb)%A2(1:n2Dim,1:NUMV) => Adam%A0(iS:iE)
      end do
    case (6)
      do iSyma=1,nSym
        iSymb = Mul(iSym,iSyma)
        if (iSymb > iSyma) cycle
        iS = iE+1
        n2Dim = n(iSyma)*n(iSymb)
        iE = iE+n2Dim*NUMV
        Adam%SB(iSymb)%A2(1:n2Dim,1:NUMV) => Adam%A0(iS:iE)
      end do
    case default
      write(u6,*) 'Allocate_SBA: Illegal case.'
      call Abend()
  end select

end subroutine Allocate_SBA

subroutine Deallocate_SBA(Adam)

  type(SBA_Type), intent(inout) :: Adam
  integer(kind=iwp) :: iSym

  do iSym=1,Adam%nSym
    nullify(Adam%SB(iSym)%A1,Adam%SB(iSym)%A2,Adam%SB(iSym)%A3)
  end do
  call mma_deallocate(Adam%A0)
  Adam%iCase = 0
  Adam%iSym = 0
  Adam%nSym = 0

end subroutine Deallocate_SBA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                  T W X Y - T Y P E   S E C T I O N                   !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Allocate_twxy(twxy,n,m,JSYM,nSym,iCase)

  type(twxy_type), target, intent(out) :: twxy
  integer(kind=iwp), intent(in) :: nSym, n(nSym), m(nSym), JSYM, iCase
  integer(kind=iwp) :: iSymx, iSymy, iSymt, iSymw, iSyma, mtwxy, iS, iE, n1, n2

  twxy%iCase = iCase
  twxy%JSYM = JSYM
  twxy%nSym = nSym
  ! *** memory for the (tw|xy) integrals --- temporary array
  mtwxy = 0
  select case (iCase)
    case (0)  ! twxy
      do iSymy=1,nSym
        if (n(iSymy) /= m(iSymy)) then
          write(u6,*) 'Allocate_twxy: iCase=0 only valid if n(:)=m(:).'
          call abend()
        end if
        iSymx = Mul(iSymy,JSYM)
        n2 = n(iSymx)*n(iSymy)
        do iSymw=iSymy,nSym    ! iSymw >= iSymy (particle symmetry)
          iSymt = Mul(isymw,JSYM)
          n1 = n(iSymt)*n(iSymw)
          mtwxy = mtwxy+n1*n2
        end do
      end do
    case (1) ! waxy
      do iSymy=1,nSym
        iSymx = Mul(iSymy,JSYM)
        n2 = n(iSymx)*n(iSymy)
        if (iSymx == iSymy) n2 = n(iSymx)*(n(iSymx)+1)/2
        if (iSymx <= iSymy) then
          do iSyma=1,nSym
            iSymw = Mul(iSyma,JSYM)
            n1 = n(iSymw)*m(iSyma)
            mtwxy = mtwxy+n1*n2
          end do
        end if
      end do
    case (2)  ! twxy
      do iSymy=1,nSym
        if (n(iSymy) /= m(iSymy)) then
          write(u6,*) 'Allocate_twxy: iCase=2 only valid if n(:)=m(:).'
          call abend()
        end if
        iSymx = Mul(iSymy,JSYM)
        n2 = n(iSymx)*n(iSymy)
        if (iSymx == iSymy) n2 = n(iSymx)*(n(iSymx)+1)/2
        if (iSymx >= iSymy) then
          do iSymw=iSymy,nSym
            iSymt = Mul(iSymw,JSYM)
            if (iSymt >= iSymw) then
              n1 = n(iSymt)*n(iSymw)
              if (iSymt == iSymw) n1 = n(iSymt)*(n(iSymt)+1)/2
              mtwxy = mtwxy+n1*n2
            end if
          end do
        end if
      end do
    case default
      write(u6,*) 'Allocate_twxy: Illegal case.'
      call Abend()
  end select

  call mma_allocate(twxy%twxy_full,mtwxy,Label='twxy')
  twxy%twxy_full(:) = Zero

  ! *** setup pointers to the symmetry blocks of (tw|xy)

  iE = 0
  select case (iCase)
    case (0)
      do iSymy=1,nSym
        iSymx = Mul(iSymy,JSYM)
        n2 = n(iSymx)*n(iSymy)
        do iSymw=iSymy,nSym   ! iSymw >= iSymy (particle symmetry)
          iSymt = Mul(isymw,JSYM)
          n1 = n(iSymt)*n(iSymw)
          iS = iE+1
          iE = iE+n1*n2
          twxy%SB(iSymw,iSymy)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE)
        end do
      end do
    case (1)
      do iSymy=1,nSym
        iSymx = Mul(iSymy,JSYM)
        n2 = n(iSymx)*n(iSymy)
        if (iSymx == iSymy) n2 = n(iSymx)*(n(iSymx)+1)/2
        if (iSymx <= iSymy) then
          do iSyma=1,nSym
            iSymw = Mul(iSyma,JSYM)
            n1 = n(iSymw)*m(iSyma)
            iS = iE+1
            iE = iE+n1*n2
            twxy%SB(iSymw,iSymx)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE)
          end do
        end if
      end do
    case (2) ! twxy
      do iSymy=1,nSym
        iSymx = Mul(iSymy,JSYM)
        n2 = n(iSymx)*n(iSymy)
        if (iSymx == iSymy) n2 = n(iSymx)*(n(iSymx)+1)/2
        if (iSymx >= iSymy) then
          do iSymw=iSymy,nSym ! iSymw >= iSymy
            iSymt = Mul(iSymw,JSYM)
            if (iSymt >= iSymw) then
              n1 = n(iSymt)*n(iSymw)
              if (iSymt == iSymw) n1 = n(iSymt)*(n(iSymt)+1)/2
              iS = iE+1
              iE = iE+n1*n2
              twxy%SB(iSymw,iSymy)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE)
              twxy%SB(iSymy,iSymw)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE) ! symmetrization
            end if
          end do
        end if
      end do
  end select

end subroutine Allocate_twxy

subroutine Deallocate_twxy(twxy)

  type(twxy_type), intent(inout) :: twxy
  integer(kind=iwp) :: iSymy, iSymw

  call mma_deallocate(twxy%twxy_full)

  ! *** setup pointers to the symmetry blocks of (tw|xy)

  do iSymy=1,8
    do iSymw=1,8
      nullify(twxy%SB(iSymw,iSymy)%A)
    end do
  end do

end subroutine Deallocate_twxy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                      G 2 - T Y P E   S E C T I O N                   !
!                                                                      !
!              Symmetry block 2-particle-like arrays                   !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Allocate_G2(Adam,n,nSym,iCase)

  type(G2_Type), target, intent(out) :: Adam
  integer(kind=iwp), intent(in) :: nSym, n(nSym), iCase
  integer(kind=iwp) :: MemTot, ijSym, iSym, jSym, kSym, lSym, iE, iS, n1, n2, n3, n4, n12, n34

  Adam%nSym = nSym
  Adam%iCase = iCase

  MemTot = 0
  select case (iCase)

    case (1)

      do ijsym=1,nsym
        do isym=1,nsym
          jsym = Mul(isym,ijsym)
          n12 = n(iSym)*n(jSym)
          do kSym=1,nSym
            lSym = Mul(kSym,ijSym)
            n34 = n(kSym)*n(lSym)
            MemTot = MemTot+n12*n34
          end do
        end do
      end do

    case default

      write(u6,*) 'Allocate_G2: illegal case valeu=',iCase
      call Abend()

  end select

  call mma_allocate(Adam%A0,MemTot,Label='G2%A0')

  iE = 0
  select case (iCase)

    case (1)

      do ijsym=1,nsym
        do isym=1,nsym
          jsym = Mul(isym,ijsym)
          n1 = n(iSym)
          n2 = n(jSym)
          n12 = n1*n2
          do kSym=1,nSym
            lSym = Mul(kSym,ijSym)
            n3 = n(kSym)
            n4 = n(lSym)
            n34 = n3*n4
            iS = iE+1
            iE = iE+n12*n34
            Adam%SB(iSym,jSym,kSym)%A4(1:n1,1:n2,1:n3,1:n4) => Adam%A0(iS:iE)
            Adam%SB(iSym,jSym,kSym)%A2(1:n12,1:n34) => Adam%A0(iS:iE)
          end do
        end do
      end do

    case default

      write(u6,*) 'What?'
      call Abend()

  end select

end subroutine Allocate_G2

subroutine Deallocate_G2(Adam)

  type(G2_Type), intent(inout) :: Adam
  integer(kind=iwp) :: iSym, jSym, kSym

  Adam%iCase = 0

  call mma_deallocate(Adam%A0)

  do iSym=1,Adam%nSym
    do jSym=1,Adam%nSym
      do kSym=1,Adam%nSym
        nullify(Adam%SB(iSym,jSym,kSym)%A4,Adam%SB(iSym,jSym,kSym)%A2)
      end do
    end do
  end do
  Adam%nSym = 0

end subroutine Deallocate_G2

subroutine Alloc_Alloc_DSBA(Array,n_Array,n,m,nSym,aCase,Label)

  type(DSBA_Type), allocatable, intent(out) :: Array(:)
  integer(kind=iwp), intent(in) :: n_Array, nSym, n(nSym), m(nSym)
  character(len=3), intent(in), optional :: aCase
  character(len=*), intent(in), optional :: Label
  integer(kind=iwp) :: i

  if (present(Label)) then
    call mma_allocate(Array,n_Array,label=Label)
  else
    call mma_allocate(Array,n_Array,label='DSBA(:)')
  end if

  if (present(aCase)) then
    do i=1,n_Array
      call Allocate_DT(Array(i),n,m,nSym,aCase)
    end do
  else
    do i=1,n_Array
      call Allocate_DT(Array(i),n,m,nSym)
    end do
  end if

# include "macros.fh"
  unused_proc(mma_allocate(Array,[0,0]))

end subroutine Alloc_Alloc_DSBA

subroutine Free_Alloc_DSBA(Array)

  type(DSBA_Type), allocatable, intent(inout) :: Array(:)
  integer(kind=iwp) :: i

  do i=lbound(Array,1),ubound(Array,1)
    call Deallocate_DT(Array(i))
  end do
  call mma_deallocate(Array)

end subroutine Free_Alloc_DSBA

subroutine Alloc_Alloc1DArray(Array,N,Label)

  type(Alloc1DArray_Type), allocatable, intent(inout) :: Array(:)
  integer(kind=iwp), intent(in) :: N(2)
  character(len=*), intent(in) :: Label
# ifdef _GARBLE_
  interface
    subroutine c_null_alloc(A) _BIND_C_
      import :: wp
      real(kind=wp), allocatable :: A(:)
    end subroutine c_null_alloc
  end interface
  integer(kind=iwp) :: i
# endif

  call mma_allocate(Array,N,label=Label)
# ifdef _GARBLE_
  ! Garbling corrupts the allocation status of allocatable components, use a hack to reset it
  do i=N(1),N(2)
    call c_null_alloc(Array(i)%A)
  end do
# endif

# include "macros.fh"
  unused_proc(mma_allocate(Array,0))

end subroutine Alloc_Alloc1DArray

subroutine Alloc2D_Alloc1DArray(Array,N1,N2,Label)

  type(Alloc1DArray_Type), allocatable, intent(inout) :: Array(:,:)
  integer(kind=iwp), intent(in) :: N1(2), N2(2)
  character(len=*), intent(in) :: Label
# ifdef _GARBLE_
  interface
    subroutine c_null_alloc(A) _BIND_C_
      import :: wp
      real(kind=wp), allocatable :: A(:)
    end subroutine c_null_alloc
  end interface
  integer(kind=iwp) :: i, j
# endif

  call mma_allocate(Array,N1,N2,label=Label)
# ifdef _GARBLE_
  ! Garbling corrupts the allocation status of allocatable components, use a hack to reset it
  do j=N2(1),N2(2)
    do i=N1(1),N1(2)
      call c_null_alloc(Array(i,j)%A)
    end do
  end do
# endif

# include "macros.fh"
  unused_proc(mma_allocate(Array,0,0))

end subroutine Alloc2D_Alloc1DArray

subroutine Alloc_Alloc2DArray(Array,N,Label)

  type(Alloc2DArray_Type), allocatable, intent(inout) :: Array(:)
  integer(kind=iwp), intent(in) :: N(2)
  character(len=*), intent(in) :: Label
# ifdef _GARBLE_
  interface
    subroutine c_null_alloc2(A) _BIND_C_
      import :: wp
      real(kind=wp), allocatable :: A(:,:)
    end subroutine c_null_alloc2
  end interface
  integer(kind=iwp) :: i
# endif

  call mma_allocate(Array,N,label=Label)
# ifdef _GARBLE_
  ! Garbling corrupts the allocation status of allocatable components, use a hack to reset it
  ! (note c_null_alloc2 is just the same as c_null_alloc)
  do i=N(1),N(2)
    call c_null_alloc2(Array(i)%A)
  end do
# endif

# include "macros.fh"
  unused_proc(mma_allocate(Array,0))

end subroutine Alloc_Alloc2DArray

subroutine Alloc_Alloc1DiArray(Array,N,Label)

  type(Alloc1DiArray_Type), allocatable, intent(inout) :: Array(:)
  integer(kind=iwp), intent(in) :: N(2)
  character(len=*), intent(in) :: Label
# ifdef _GARBLE_
  interface
    subroutine c_null_alloc3(A) _BIND_C_
      import :: iwp
      integer(kind=iwp), allocatable :: A(:)
    end subroutine c_null_alloc3
  end interface
  integer(kind=iwp) :: i
# endif

  call mma_allocate(Array,N,label=Label)
# ifdef _GARBLE_
  ! Garbling corrupts the allocation status of allocatable components, use a hack to reset it
  ! (note c_null_alloc3 is just the same as c_null_alloc)
  do i=N(1),N(2)
    call c_null_alloc3(Array(i)%A)
  end do
# endif

# include "macros.fh"
  unused_proc(mma_allocate(Array,0))

end subroutine Alloc_Alloc1DiArray

subroutine Free_Alloc1DArray(Array)

  type(Alloc1DArray_Type), allocatable, intent(inout) :: Array(:)
  integer(kind=iwp) :: i

  do i=lbound(Array,1),ubound(Array,1)
    call mma_deallocate(Array(i)%A,safe='*')
  end do
  call mma_deallocate(Array)

end subroutine Free_Alloc1DArray

subroutine Free2D_Alloc1DArray(Array)

  type(Alloc1DArray_Type), allocatable, intent(inout) :: Array(:,:)
  integer(kind=iwp) :: i, j

  do j=lbound(Array,2),ubound(Array,2)
    do i=lbound(Array,1),ubound(Array,1)
      call mma_deallocate(Array(i,j)%A,safe='*')
    end do
  end do
  call mma_deallocate(Array)

end subroutine Free2D_Alloc1DArray

subroutine Free_Alloc2DArray(Array)

  type(Alloc2DArray_Type), allocatable, intent(inout) :: Array(:)
  integer(kind=iwp) :: i

  do i=lbound(Array,1),ubound(Array,1)
    call mma_deallocate(Array(i)%A,safe='*')
  end do
  call mma_deallocate(Array)

end subroutine Free_Alloc2DArray

subroutine Free_Alloc1DiArray(Array)

  type(Alloc1DiArray_Type), allocatable, intent(inout) :: Array(:)
  integer(kind=iwp) :: i

  do i=lbound(Array,1),ubound(Array,1)
    call mma_deallocate(Array(i)%A,safe='*')
  end do
  call mma_deallocate(Array)

end subroutine Free_Alloc1DiArray

! Define dsba_mma_allo_1D, dsba_mma_allo_1D_lim, dsba_mma_free_1D
!        a1da_mma_allo_1D, a1da_mma_allo_1D_lim, a1da_mma_free_1D,
!        a1da_mma_allo_2D, a1da_mma_allo_2D_lim, a1da_mma_free_2D
!        a2da_mma_allo_1D, a2da_mma_allo_1D_lim, a2da_mma_free_1D
!        a1ia_mma_allo_1D, a1ia_mma_allo_1D_lim, a1ia_mma_free_1D,

! (using _NO_GARBLE_ because all members are initialized)
#define _TYPE_ type(DSBA_Type)
#  define _NO_GARBLE_
#  define _SUBR_NAME_ dsba_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'dsba_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

#define _TYPE_ type(Alloc1DArray_Type)
#  define _SUBR_NAME_ a1da_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'a1da_mma'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  define _DIMENSIONS_ 2
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

#define _TYPE_ type(Alloc2DArray_Type)
#  define _SUBR_NAME_ a2da_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'a2da_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

#define _TYPE_ type(Alloc1DiArray_Type)
#  define _SUBR_NAME_ a1ia_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'a1ia_mma'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module Data_Structures
