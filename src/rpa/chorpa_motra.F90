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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoRPA_MOTra(includeFrozen,includeDeleted)

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Transform Cholesky vectors to MO basis.
!
! TODO/FIXME:
! 1. This routine computes all MO blocks (ij, ai, ab), even though
!    we may only need some of them. A more flexible interface would
!    be nice to have. Presumably not a performance issue in RPA,
!    though (remains to be verified).
! 2. For unrestricted calculations, the alpha and beta
!    transformations are done separately, which means that the AO
!    vectors are read twice. Simultaneous transformation would be
!    desirable!

use RPA_globals, only: CMO, nBas, nDel, nFro, nOcc, nOrb, nSym, nVir
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: includeFrozen, includeDeleted
integer(kind=iwp) :: nSpin, iSym, l_lCMO, lnFro, lnOcc, Zeros, lnVir, lnDel, iSpin
character(len=6) BName(2)
integer(kind=iwp), allocatable :: loc(:,:)
real(kind=wp), allocatable :: lCMO(:)
integer(kind=iwp), external :: RPA_iUHF
character(len=*), parameter :: SecNam = 'ChoRPA_MOTra'

nSpin = RPA_iUHF()
if (nSpin == 1) then
  BName(1) = 'MOVECS'
  BName(2) = 'unused'
else if (nSpin == 2) then
  BName(1) = 'MOVECa'
  BName(2) = 'MOVECb'
else
  call RPA_Warn(3,SecNam//': illegal nSpin')
  BName(1) = 'unused'
  BName(2) = 'unused'
end if
l_lCMO = nBas(1)**2
do iSym=2,nSym
  l_lCMO = l_lCMO+nBas(iSym)**2
end do
call mma_allocate(lCMO,l_lCMO,label='locCMO')
call mma_allocate(loc,nSym,5,label='local')
lnFro = 1
lnOcc = 2
Zeros = 3
lnVir = 4
lnDel = 5
loc(:,Zeros) = 0
if (includeFrozen) lnFro = Zeros
if (includeDeleted) lnDel = Zeros

do iSpin=1,nSpin

  ! Set orbital blocks
  if (includeFrozen) then
    do iSym=1,nSym
      loc(iSym,lnOcc) = nFro(iSym,iSpin)+nOcc(iSym,iSpin)
    end do
  else
    loc(:,lnFro) = nFro(:,iSpin)
    loc(:,lnOcc) = nOcc(:,iSpin)
  end if
  if (includeDeleted) then
    do iSym=1,nSym
      loc(:,lnVir-1+iSym) = nVir(iSym,iSpin)+nDel(iSym,iSpin)
    end do
  else
    loc(:,lnVir) = nVir(:,iSpin)
    loc(:,lnDel) = nDel(:,iSpin)
  end if
  ! Reorder CMO array
  call ChoRPA_MOTra_ReorderCMO(nSym,nBas,nOrb,CMO(:,iSpin),lCMO)
  ! Set base name for MO files
  ! Transform Cholesky vectors
  call Cho_MOTra_Inner(lCMO,l_lCMO,nSym,nBas,loc(:,lnFro),loc(:,lnOcc),loc(:,Zeros),loc(:,lnVir),loc(:,lnDel),BName(iSpin), &
                       .false.,0,.false.)

end do

call mma_deallocate(lCMO)
call mma_deallocate(loc)

end subroutine ChoRPA_MOTra
