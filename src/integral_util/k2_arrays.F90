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

module k2_arrays

! DeDe is an array with desymmetrized 1-particle densities used
! in the integral direct construction of the Fock matrix. In this
! the array contain the subblocks with their associatated base
! pointer.
! ipDeDe: the set of densities to be used sorted according to the
!  integral shell pairs they contribute. Pointers to the individual
!  matrices are stored in ipOffD.
! ipD00: some of the pointers in ipOffD point to an "empty" slot
!  and this pointer is to this part of DeDe.
! ipDijS: points to an auxiliary pice of memory which is used in
!  case that used a subset of the elements of a matrix is used. In
!  this case picky_ will extract those elements and put them into
!  this part of DeDe on the fly.

implicit none

integer, allocatable :: ipOffD(:,:)
real*8, allocatable :: FT(:), DeDe(:)
integer ipDeDe, ipD00, ipDijS
integer nDeDe, nDeDe_DFT, MaxDe, MxDij, MxFT, nFT
logical :: DoGrad_ = .false., DoHess_ = .false.
real*8, target, allocatable :: Fq(:), Dq(:)
real*8, pointer :: pFq(:) => null(), pDq(:) => null()
real*8, allocatable :: Aux(:)
integer, allocatable :: iSOSym(:,:)
logical :: XMem = .false.
real*8, allocatable, target :: Sew_Scr(:)

type BraKet_Type
  real*8, pointer :: Zeta(:)
  real*8, pointer :: ZInv(:)
  real*8, pointer :: KappaAB(:)
  real*8, pointer :: P(:,:)
  real*8, pointer :: xA(:)
  real*8, pointer :: xB(:)
  real*8, pointer :: Eta(:)
  real*8, pointer :: EInv(:)
  real*8, pointer :: KappaCD(:)
  real*8, pointer :: Q(:,:)
  real*8, pointer :: xG(:)
  real*8, pointer :: xD(:)
  real*8, pointer :: xPre(:)
  integer, pointer :: IndZet(:)
  integer, pointer :: IndEta(:)
end type BraKet_Type

type(BraKet_Type) BraKet
real*8, allocatable, target :: BraKet_Base_R(:)
integer, allocatable, target :: BraKet_Base_I(:)

contains

subroutine Create_BraKet_Base(nZeta)

  use stdalloc, only: mma_allocate

  integer nZeta
  integer Mem

  Mem = nZeta*16
  if (DoHess_) Mem = Mem+nZeta**2
  call mma_allocate(BraKet_Base_R,Mem,Label='Base_R')
  Mem = (nZeta+1)*2
  call mma_allocate(BraKet_Base_I,Mem,Label='Base_I')

end subroutine Create_BraKet_Base

subroutine Destroy_BraKet_Base()

  use stdalloc, only: mma_deallocate

  if (allocated(BraKet_Base_R)) call mma_deallocate(BraKet_Base_R)
  if (allocated(BraKet_Base_I)) call mma_deallocate(BraKet_Base_I)

end subroutine Destroy_BraKet_Base

subroutine Create_BraKet(nZeta,nEta)

  integer nZeta, nEta
  integer iS, iE

  if ((.not. allocated(BraKet_base_R)) .or. (.not. allocated(BraKet_base_I))) then
    write(6,*) 'Braket_Base not allocated!'
    call Abend()
  end if

  if (nZeta*nEta == 0) return

  iE = 0

  if (nZeta /= 0) then
    iS = iE+1
    iE = iE+nZeta
    Braket%Zeta(1:nZeta) => BraKet_Base_R(iS:iE)
    iS = iE+1
    iE = iE+nZeta
    Braket%ZInv(1:nZeta) => BraKet_Base_R(iS:iE)
    iS = iE+1
    iE = iE+nZeta
    Braket%KappaAB(1:nZeta) => BraKet_Base_R(iS:iE)
    iS = iE+1
    iE = iE+3*nZeta
    Braket%P(1:nZeta,1:3) => BraKet_Base_R(iS:iE)
    iS = iE+1
    iE = iE+nZeta
    Braket%xA(1:nZeta) => BraKet_Base_R(iS:iE)
    iS = iE+1
    iE = iE+nZeta
    Braket%xB(1:nZeta) => BraKet_Base_R(iS:iE)
  end if

  if (nEta /= 0) then
    iS = iE+1
    iE = iE+nEta
    Braket%Eta(1:nEta) => BraKet_Base_R(iS:iE)
    iS = iE+1
    iE = iE+nEta
    Braket%EInv(1:nEta) => BraKet_Base_R(iS:iE)
    iS = iE+1
    iE = iE+nEta
    Braket%KappaCD(1:nEta) => BraKet_Base_R(iS:iE)
    iS = iE+1
    iE = iE+3*nEta
    Braket%Q(1:nEta,1:3) => BraKet_Base_R(iS:iE)
    iS = iE+1
    iE = iE+nEta
    Braket%xG(1:nEta) => BraKet_Base_R(iS:iE)
    iS = iE+1
    iE = iE+nEta
    Braket%xD(1:nEta) => BraKet_Base_R(iS:iE)
  end if

  if ((nZeta*nEta /= 0) .and. DoHess_) then
    iS = iE+1
    iE = iE+nZeta*nEta
    Braket%xPre(1:nEta) => BraKet_Base_R(iS:iE)
  end if

  iE = 0

  if (nZeta /= 0) then
    iS = iE+1
    iE = iE+nZeta+1
    Braket%IndZet(1:nZeta+1) => BraKet_Base_I(iS:iE)
  end if

  if (nEta /= 0) then
    iS = iE+1
    iE = iE+nEta+1
    Braket%IndEta(1:nEta+1) => BraKet_Base_I(iS:iE)
  end if

end subroutine Create_BraKet

subroutine Destroy_BraKet()

  nullify(Braket%Zeta,Braket%ZInv,Braket%KappaAB,Braket%P,Braket%xA,Braket%xB,Braket%Eta,Braket%EInv,Braket%KappaCD,Braket%Q, &
          Braket%xG,Braket%xD,Braket%IndZet,Braket%IndEta,BraKet%xPre)

end subroutine Destroy_BraKet

end module k2_arrays
