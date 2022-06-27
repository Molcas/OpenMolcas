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
! Copyright (C) 2005, Giovanni Ghigo                                   *
!***********************************************************************

module MkSubs
! Module containing MkCouSB* and MkExSB* subroutines, to avoid explicit interfaces

use Cho_Tra, only: nAsh, nIsh, nSsh, TCVX
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
private

public :: MkCouSB11, MkCouSB21, MkCouSB22, MkCouSB31, MkCouSB32, MkCouSB33, &
          MkExSB11, MkExSB12, MkExSB13, MkExSB21, MkExSB22, MkExSB23, MkExSB31, MkExSB32, MkExSB33

contains

subroutine MkCouSB11(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           July 2005                                                  *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(1,1) (p,q inactive) of the      *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: LenAB, LenSB
  real(kind=wp), allocatable :: Lij(:)

  ! SubBlock 1 1
  LenSB = nIsh(iSymA)*nIsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  ! Define Lab
  LenAB = LenSB
  !---------------------------------------------------------------------
  !if (IfTest) then
  !  write(u6,*) '     MkCouSB11: TCVA(',iSymA,',',iSymB,')'
  !  write(u6,'(8F10.6)') TCVX(1,iSymA,iSymB)%A(:,:)
  !  call XFlush(u6)
  !end if
  !---------------------------------------------------------------------

  ! Build Lij
  call mma_allocate(Lij,NumV,Label='Lij')
  call MkLij(iSymI,iSymJ,iI,iJ,numV,Lij)

  ! Generate the SubBlock
  call DGEMM_('N','N',LenAB,1,numV,One,TCVX(1,iSymA,iSymB)%A,LenAB,Lij,numV,Zero,AddSB,LenSB)

  call mma_deallocate(Lij)

  return

end subroutine MkCouSB11

subroutine MkCouSB21(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           July 2005                                                  *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(2,1) (p active, q inactive) of  *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: LenAB, LenSB
  real(kind=wp), allocatable :: AddSBt(:), Lij(:)

  ! SubBlock 2 1
  LenSB = nAsh(iSymA)*nIsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  call mma_allocate(AddSBt,LenSB,Label='AddSBt')

  ! Define Lab
  LenAB = LenSB
  !---------------------------------------------------------------------
  !if (IfTest) then
  !  write(u6,*) '     MkCouSB21: TCVB(',iSymA,',',iSymB,')'
  !  write(u6,'(8F10.6)') TCVX(2,iSymA,iSymB)%A(:,:)
  !  call XFlush(u6)
  !end if
  !---------------------------------------------------------------------

  ! Build Lij
  call mma_allocate(Lij,NumV,Label='Lij')
  call MkLij(iSymI,iSymJ,iI,iJ,numV,Lij)

  ! Generate the SubBlock
  call DGEMM_('N','N',LenAB,1,numV,One,TCVX(2,iSymA,iSymB)%A,LenAB,Lij,NumV,Zero,AddSBt,LenSB)
  call Trnsps(nAsh(iSymA),nIsh(iSymB),AddSBt,AddSB)

  call mma_deallocate(Lij)
  call mma_deallocate(AddSBt)

  return

end subroutine MkCouSB21

subroutine MkCouSB22(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           July 2005                                                  *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(2,2) (p,q active) of the        *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: LenAB, LenSB
  real(kind=wp), allocatable :: Lij(:)

  ! SubBlock 2 2
  LenSB = nAsh(iSymA)*nAsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  ! Define Lab
  LenAB = LenSB
  !---------------------------------------------------------------------
  !if (IfTest) then
  !  write(u6,*) '     MkCouSB22: TCVD(',iSymA,',',iSymB,')'
  !  write(u6,'(8F10.6)') TCVX(4,iSymA,iSymB)%A(:,:)
  !  call XFlush(u6)
  !end if
  !---------------------------------------------------------------------

  ! Build Lij
  call mma_allocate(Lij,NumV,Label='Lij')
  call MkLij(iSymI,iSymJ,iI,iJ,numV,Lij)

  ! Generate the SubBlock
  call DGEMM_('N','N',LenAB,1,numV,One,TCVX(4,iSymA,iSymB)%A,LenAB,Lij,NumV,Zero,AddSB,LenSB)

  call mma_deallocate(Lij)

  return

end subroutine MkCouSB22

subroutine MkCouSB31(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           July 2005                                                  *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(3,1) (p secondary, q inactive)  *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: LenAB, LenSB
  real(kind=wp), allocatable :: AddSBt(:), Lij(:)

  ! SubBlock 3 1
  LenSB = nSsh(iSymA)*nIsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  call mma_allocate(AddSBt,LenSB,Label='AddSBt')

  ! Define Lab
  LenAB = LenSB
  !---------------------------------------------------------------------
  !if (IfTest) then
  !  write(u6,*) '     MkCouSB31: TCVB(',iSymA,',',iSymB,')'
  !  write(u6,'(8F10.6)') TCVX(3,iSymA,iSymB)%A(:,:)
  !  call XFlush(u6)
  !end if
  !---------------------------------------------------------------------

  ! Build Lij
  call mma_allocate(Lij,NumV,Label='Lij')
  call MkLij(iSymI,iSymJ,iI,iJ,numV,Lij)

  ! Generate the SubBlock
  call DGEMM_('N','N',LenAB,1,numV,One,TCVX(3,iSymA,iSymB)%A,LenAB,Lij,NumV,Zero,AddSBt,LenSB)
  call Trnsps(nSsh(iSymA),nIsh(iSymB),AddSBt,AddSB)

  call mma_deallocate(Lij)
  call mma_deallocate(AddSBt)

  return

end subroutine MkCouSB31

subroutine MkCouSB32(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           July 2005                                                  *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(3,2) (p secondary, q active) of *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: LenAB, LenSB
  real(kind=wp), allocatable :: AddSBt(:), Lij(:)

  ! SubBlock 3 2
  LenSB = nSsh(iSymA)*nAsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  call mma_allocate(AddSBt,LenSB,Label='AddSBt')

  ! Define Lab
  LenAB = LenSB
  !---------------------------------------------------------------------
  !if (IfTest) then
  !  write(u6,*) '     MkCouSB32: TCVE(',iSymA,',',iSymB,')'
  !  write(u6,'(8F10.6)') TCVX(5,iSymA,iSymB)%A(:,:)
  !  call XFlush(u6)
  !end if
  !---------------------------------------------------------------------

  ! Build Lij
  call mma_allocate(Lij,NumV,Label='Lij')
  call MkLij(iSymI,iSymJ,iI,iJ,numV,Lij)

  ! Generate the SubBlock
  call DGEMM_('N','N',LenAB,1,numV,One,TCVX(5,iSymA,iSymB)%A,LenAB,Lij,NumV,Zero,AddSBt,LenSB)
  call Trnsps(nSsh(iSymA),nAsh(iSymB),AddSBt,AddSB)

  call mma_deallocate(Lij)
  call mma_deallocate(AddSBt)

  return

end subroutine MkCouSB32

subroutine MkCouSB33(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           July 2005                                                  *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(3,3) (p,q secondary) of the     *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: LenAB, LenSB
  real(kind=wp), allocatable :: Lij(:)

  ! SubBlock 3 3
  LenSB = nSsh(iSymA)*nSsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  ! Define Lab
  LenAB = LenSB
  !---------------------------------------------------------------------
  !if (IfTest) then
  !  write(u6,*) '     MkCouSB33: TCVF(',iSymA,',',iSymB,')'
  !  write(u6,'(8F10.6)') TCVX(6,iSymA,iSymB)%A(:,:)
  !  call XFlush(u6)
  !end if
  !-------------------------------------------------------------------

  ! Build Lij
  call mma_allocate(Lij,NumV,Label='Lij')
  call MkLij(iSymI,iSymJ,iI,iJ,numV,Lij)

  ! Generate the SubBlock
  call DGEMM_('N','N',LenAB,1,numV,One,TCVX(6,iSymA,iSymB)%A,LenAB,Lij,NumV,Zero,AddSB,LenSB)

  call mma_deallocate(Lij)

  return

end subroutine MkCouSB33

subroutine MkExSB11(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(1,1) (p,q inactive) of the      *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: iIx, LenSB, LxType
  logical(kind=iwp) :: SameLx
  real(kind=wp), allocatable :: Lx0(:), Ly0(:)

  ! SubBlock 1 1
  LenSB = nIsh(iSymA)*nIsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='LenSB')

  ! Build Lx
  call mma_allocate(Lx0,nIsh(iSymA)*numV,Label='Lx0')
  LxType = 0
  iIx = 0
  SameLx = .false.
  call MkL1(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

  ! Build Ly
  call mma_allocate(Ly0,nIsh(iSymB)*numV,Label='Ly0')
  if (iSymA == iSymB) SameLx = .true.
  call MkL1(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

  ! Generate the SubBlock (Ly*Lx)
  if (.not. SameLx) then
    call DGEMM_('N','T',nIsh(iSymB),nIsh(iSymA),numV,One,Ly0,nIsh(iSymB),Lx0,nIsh(iSymA),Zero,AddSB,nIsh(iSymB))
  else
    call DGEMM_('N','T',nIsh(iSymA),nIsh(iSymA),numV,One,Lx0,nIsh(iSymA),Lx0,nIsh(iSymA),Zero,AddSB,nIsh(iSymA))
  end if

  call mma_deallocate(Ly0)
  call mma_deallocate(Lx0)

  return

end subroutine MkExSB11

subroutine MkExSB12(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(1,2) (p,q inactive,active) of   *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: iIx, LenSB, LxType
  logical(kind=iwp) :: SameLx
  real(kind=wp), allocatable :: Lx0(:), Ly0(:)

  ! SubBlock 1 2
  LenSB = nIsh(iSymA)*nAsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  ! Build Lx
  call mma_allocate(Lx0,nIsh(iSymA)*numV,Label='Lx0')
  LxType = 0
  iIx = 0
  SameLx = .false.
  call MkL1(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

  ! Build Ly
  call mma_allocate(Ly0,nAsh(iSymB)*numV,Label='Ly0')
  call MkL2(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

  ! Generate the SubBlock
  call DGEMM_('N','T',nAsh(iSymB),nIsh(iSymA),numV,One,Ly0,nAsh(iSymB),Lx0,nIsh(iSymA),Zero,AddSB,nAsh(iSymB))

  call mma_deallocate(Ly0)
  call mma_deallocate(Lx0)

  return

end subroutine MkExSB12

subroutine MkExSB13(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(1,3) (p,q inactive,secondary)   *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: iIx, LenSB, LxType
  logical(kind=iwp) :: SameLx
  real(kind=wp), allocatable :: Lx0(:), Ly0(:)

  ! SubBlock 1 3
  LenSB = nIsh(iSymA)*nSsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  ! Build Lx
  call mma_allocate(Lx0,nIsh(iSymA)*numV,Label='Lx0')
  LxType = 0
  iIx = 0
  SameLx = .false.
  call MkL1(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

  ! Build Ly
  call mma_allocate(Ly0,nSsh(iSymB)*numV,Label='Ly0')
  call MkL3(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

  ! Generate the SubBlock
  call DGEMM_('N','T',nSsh(iSymB),nIsh(iSymA),numV,One,Ly0,nSsh(iSymB),Lx0,nIsh(iSymA),Zero,AddSB,nSsh(iSymB))

  call mma_deallocate(Ly0)
  call mma_deallocate(Lx0)

  return

end subroutine MkExSB13

subroutine MkExSB21(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddSBt)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(2,1) (p,q active,inactive) of   *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  real(kind=wp), intent(in) :: AddSBt(*)
  integer(kind=iwp) :: iIx, LenSB, LxType
  logical(kind=iwp) :: SameLx
  real(kind=wp), allocatable :: Lx0(:), Ly0(:)

  ! SubBlock 2 1
  LenSB = nAsh(iSymA)*nIsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')
  if ((iSymA == iSymB) .and. (iSymI == iSymJ) .and. (iI == iJ)) then
    ! SB 2,1 = (SB 1,2)+
    call Trnsps(nAsh(iSymB),nIsh(iSymA),AddSBt,AddSB)
    return
  end if

  ! Build Lx
  call mma_allocate(Lx0,nAsh(iSymA)*numV,Label='Lx0')
  LxType = 0
  iIx = 0
  SameLx = .false.
  call MkL2(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

  ! Build Ly
  call mma_allocate(Ly0,nIsh(iSymB)*numV,Label='Ly0')
  call MkL1(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

  ! Generate the SubBlock
  call DGEMM_('N','T',nIsh(iSymB),nAsh(iSymA),numV,One,Ly0,nIsh(iSymB),Lx0,nAsh(iSymA),Zero,AddSB,nIsh(iSymB))

  call mma_deallocate(Ly0)
  call mma_deallocate(Lx0)

  return

end subroutine MkExSB21

subroutine MkExSB22(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(2,2) (p,q active) of the        *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: iIx, LenSB, LxType
  logical(kind=iwp) :: SameLx
  real(kind=wp), allocatable :: Lx0(:), Ly0(:)

  ! SubBlock 2 2
  LenSB = nAsh(iSymA)*nAsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  ! Build Lx
  call mma_allocate(Lx0,nAsh(iSymA)*numV,Label='Lx0')
  LxType = 0
  iIx = 0
  SameLx = .false.
  call MkL2(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

  ! Build Ly
  call mma_allocate(Ly0,nAsh(iSymB)*numV,Label='Ly0')
  if (iSymA == iSymB) SameLx = .true.
  call MkL2(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

  ! Generate the SubBlock
  if (.not. SameLx) then
    call DGEMM_('N','T',nAsh(iSymB),nAsh(iSymA),numV,One,Ly0,nAsh(iSymB),Lx0,nAsh(iSymA),Zero,AddSB,nAsh(iSymB))
  else
    call DGEMM_('N','T',nAsh(iSymA),nAsh(iSymA),numV,One,Lx0,nAsh(iSymA),Lx0,nAsh(iSymA),Zero,AddSB,nAsh(iSymA))
  end if

  call mma_deallocate(Ly0)
  call mma_deallocate(Lx0)

  return

end subroutine MkExSB22

subroutine MkExSB23(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(2,3) (p,q active,secondary) of  *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: iIx, LenSB, LxType
  logical(kind=iwp) :: SameLx
  real(kind=wp), allocatable :: Lx0(:), Ly0(:)

  ! SubBlock 2 3
  LenSB = nAsh(iSymA)*nSsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  ! Build Lx
  call mma_allocate(Lx0,nAsh(iSymA)*numV,Label='Lx0')
  LxType = 0
  iIx = 0
  SameLx = .false.
  call MkL2(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

  ! Build Ly
  call mma_allocate(Ly0,nSsh(iSymB)*numV,Label='Ly0')
  call MkL3(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

  ! Generate the SubBlock
  call DGEMM_('N','T',nSsh(iSymB),nAsh(iSymA),numV,One,Ly0,nSsh(iSymB),Lx0,nAsh(iSymA),Zero,AddSB,nSsh(iSymB))

  call mma_deallocate(Ly0)
  call mma_deallocate(Lx0)

  return

end subroutine MkExSB23

subroutine MkExSB31(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddSBt)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(3,1) (p,q secondary,inactive)   *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  real(kind=wp), intent(in) :: AddSBt(*)
  integer(kind=iwp) :: iIx, LenSB, LxType
  logical(kind=iwp) :: SameLx
  real(kind=wp), allocatable :: Lx0(:), Ly0(:)

  ! SubBlock 3 1
  LenSB = nSsh(iSymA)*nIsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')
  if ((iSymA == iSymB) .and. (iSymI == iSymJ) .and. (iI == iJ)) then
    ! SB 3,1 = (SB 1,3)+
    call Trnsps(nSsh(iSymB),nIsh(iSymA),AddSBt,AddSB)
    return
  end if

  ! Build Lx
  call mma_allocate(Lx0,nSsh(iSymA)*numV,Label='Lx0')
  LxType = 0
  iIx = 0
  SameLx = .false.
  call MkL3(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

  ! Build Ly
  call mma_allocate(Ly0,nIsh(iSymB)*numV,Label='Ly0')
  call MkL1(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

  ! Generate the SubBlock
  call DGEMM_('N','T',nIsh(iSymB),nSsh(iSymA),numV,One,Ly0,nIsh(iSymB),Lx0,nSsh(iSymA),Zero,AddSB,nIsh(iSymB))

  call mma_deallocate(Ly0)
  call mma_deallocate(Lx0)

  return

end subroutine MkExSB31

subroutine MkExSB32(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddSBt)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           February 2005                                              *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(3,2) (p,q secondary,active) of  *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  real(kind=wp), intent(in) :: AddSBt(*)
  integer(kind=iwp) :: iIx, LenSB, LxType
  logical(kind=iwp) :: SameLx
  real(kind=wp), allocatable :: Lx0(:), Ly0(:)

  ! SubBlock 3 2
  LenSB = nSsh(iSymA)*nAsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')
  if ((iSymA == iSymB) .and. (iSymI == iSymJ) .and. (iI == iJ)) then
    ! SB 3,2 = (SB 2,3)+
    call Trnsps(nSsh(iSymB),nAsh(iSymA),AddSBt,AddSB)
    return
  end if

  ! Build Lx
  call mma_allocate(Lx0,nSsh(iSymA)*numV,Label='Lx0')
  LxType = 0
  iIx = 0
  SameLx = .false.
  call MkL3(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

  ! Build Ly
  call mma_allocate(Ly0,nAsh(iSymB)*numV,Label='Ly0')
  call MkL2(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

  ! Generate the SubBlock
  call DGEMM_('N','T',nAsh(iSymB),nSsh(iSymA),numV,One,Ly0,nAsh(iSymB),Lx0,nSsh(iSymA),Zero,AddSB,nAsh(iSymB))

  call mma_deallocate(Ly0)
  call mma_deallocate(Lx0)

  return

end subroutine MkExSB32

subroutine MkExSB33(AddSB,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           January-February 2005                                      *
!----------------------------------------------------------------------*
! Purpuse:  Generation of the SubBlock(3,3) (p,q secondary) of the     *
!           two-electron integral matrix for each i,j occupied couple. *
!***********************************************************************

  real(kind=wp), allocatable, intent(out) :: AddSB(:)
  integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV
  integer(kind=iwp) :: iIx, LenSB, LxType
  logical(kind=iwp) :: SameLx
  real(kind=wp), allocatable :: Lx0(:), Ly0(:)

  ! SubBlock 3 3
  LenSB = nSsh(iSymA)*nSsh(iSymB)
  call mma_allocate(AddSB,LenSB,Label='AddSB')

  ! Build Lx
  call mma_allocate(Lx0,nSsh(iSymA)*numV,Label='Lx0')
  LxType = 0
  iIx = 0
  SameLx = .false.
  call MkL3(iSymA,iSymI,iI,numV,LxType,iIx,Lx0,SameLx)

  ! Build Ly
  call mma_allocate(Ly0,nSsh(iSymB)*numV,Label='Ly0')
  if (iSymA == iSymB) SameLx = .true.
  call MkL3(iSymB,iSymJ,iJ,numV,LxType,iIx,Ly0,SameLx)

  ! Generate the SubBlock
  if (.not. SameLx) then
    call DGEMM_('N','T',nSsh(iSymB),nSsh(iSymA),numV,One,Ly0,nSsh(iSymB),Lx0,nSsh(iSymA),Zero,AddSB,nSsh(iSymB))
  else
    call DGEMM_('N','T',nSsh(iSymA),nSsh(iSymA),numV,One,Lx0,nSsh(iSymA),Lx0,nSsh(iSymA),Zero,AddSB,nSsh(iSymA))
  end if

  call mma_deallocate(Ly0)
  call mma_deallocate(Lx0)

  return

end subroutine MkExSB33

end module MkSubs
