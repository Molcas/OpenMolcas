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
!  this case picky will extract those elements and put them into
!  this part of DeDe on the fly.

use Definitions, only: wp, iwp

implicit none
private

type BraKet_Type
  real(kind=wp), pointer :: Zeta(:) => null(), ZInv(:) => null(), KappaAB(:) => null(), P(:,:) => null(), xA(:) => null(), &
                            xB(:) => null(), Eta(:) => null(), EInv(:) => null(), KappaCD(:) => null(), Q(:,:) => null(), &
                            xG(:) => null(), xD(:) => null(), xPre(:) => null()
  integer(kind=iwp), pointer :: IndZet(:) => null(), IndEta(:) => null()
end type BraKet_Type

integer(kind=iwp) :: ipD00, ipDeDe, ipDijS, ipDijS2, MaxDe, MxDij, MxFT, nDeDe, nDeDe_DFT, nFT
real(kind=wp), pointer :: pDq(:) => null(), pFq(:) => null()
logical(kind=iwp) :: DoGrad_ = .false., DoHess_ = .false., XMem = .false.
type(BraKet_Type) :: BraKet
integer(kind=iwp), allocatable :: ipOffD(:,:), ipOffDA(:,:), iSOSym(:,:)
integer(kind=iwp), allocatable, target :: BraKet_Base_I(:)
real(kind=wp), allocatable, target :: Aux(:), BraKet_Base_R(:), DeDe(:), DeDe2(:), Dq(:), Fq(:), FT(:), Sew_Scr(:)

public :: Aux, BraKet, Create_BraKet, Create_BraKet_Base, DeDe, DeDe2, Destroy_BraKet, Destroy_BraKet_Base, DoGrad_, DoHess_, Dq, &
          Fq, FT, ipD00, ipDeDe, ipDijS, ipDijS2, ipOffD, ipOffDA, iSOSym, MaxDe, MxDij, MxFT, nDeDe, nDeDe_DFT, nFT, pDq, pFq, &
          Sew_Scr, XMem

contains

subroutine Create_BraKet_Base(nZeta)

  use stdalloc, only: mma_allocate

  integer(kind=iwp), intent(in) :: nZeta
  integer(kind=iwp) :: Mem

  Mem = nZeta*16
  if (DoHess_) Mem = Mem+nZeta**2
  call mma_allocate(BraKet_Base_R,Mem,Label='Base_R')
  Mem = (nZeta+1)*2
  call mma_allocate(BraKet_Base_I,Mem,Label='Base_I')

end subroutine Create_BraKet_Base

subroutine Destroy_BraKet_Base()

  use stdalloc, only: mma_deallocate

  call mma_deallocate(BraKet_Base_R,safe='*')
  call mma_deallocate(BraKet_Base_I,safe='*')

end subroutine Destroy_BraKet_Base

subroutine Create_BraKet(nZeta,nEta)

  use Definitions, only: u6

  integer(kind=iwp), intent(in) :: nZeta, nEta
  integer(kind=iwp) :: iE, iS

  if (.not. (allocated(BraKet_base_R) .and. allocated(BraKet_base_I))) then
    write(u6,*) 'Braket_Base not allocated!'
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
