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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine DFroz(Dlt,nDlt,CMO,nCMO,OccNo)
!***********************************************************************
!                                                                      *
!     purpose: Compute contribution to the density matrix from frozen  *
!              orbitals                                                *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use InfSCF, only: nBas, nFro, nnB, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDlt, nCMO
real(kind=wp), intent(out) :: Dlt(nDlt)
real(kind=wp), intent(in) :: CMO(nCMO)
integer(kind=iwp), intent(in) :: OccNo(*)
integer(kind=iwp) :: i, iCol, ij, iOrb, ipCff, ipDlt, ipOcc, iRow, iStrtN, iStrtO, iSym, lth, nBs, nFr, nOF, nOr
real(kind=wp) :: rSum
real(kind=wp), allocatable :: NewOcc(:)
real(kind=wp), parameter :: Scal = Two, SScale = One

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Allocate memory for new occupation numbers
call mma_allocate(NewOcc,nnB,Label='NewOcc')

! Set proper occupation numbers
iStrtN = 0
iStrtO = 0
do iSym=1,nSym
  do iOrb=1,nOrb(iSym)
    NewOcc(iStrtN+iOrb) = Zero
    if (iOrb <= nFro(iSym)) NewOcc(iStrtN+iOrb) = OccNo(iStrtO+iOrb)
  end do
  iStrtN = iStrtN+nOrb(iSym)
  iStrtO = iStrtO+nOrb(iSym)
end do

! Compute density contribution

ipCff = 0
ipDlt = 0
ipOcc = 0
do iSym=1,nSym

  nBs = nBas(iSym)
  nOr = nOrb(iSym)
  nFr = nFro(iSym)

  nOF = nOr-nFr
  lth = nTri_Elem(nBs)

  ipCff = ipCff+nBs*nFr
  do iRow=1,nBs
    rSum = Zero
    ij = -1
    do i=nFr+1,nOr
      ij = ij+1
      !if (NewOcc(ipOcc+i) == Zero) exit
      rSum = rSum+NewOcc(ipOcc+i)*CMO(ipCff+iRow+ij*nBs)*CMO(ipCff+iRow+ij*nBs)
    end do
    Dlt(ipDlt+iTri(iRow,iRow)) = rSum*SScale

    do iCol=1,iRow-1
      rSum = Zero
      ij = -1
      do i=nFr+1,nOr
        ij = ij+1
        !if (NewOcc(ipOcc+i) == Zero) exit
        rSum = rSum+NewOcc(ipOcc+i)*CMO(ipCff+iRow+ij*nBs)*CMO(ipCff+iCol+ij*nBs)
      end do
      Dlt(ipDlt+iTri(iRow,iCol)) = Scal*rSum
    end do
  end do
# ifdef _DEBUGPRINT_
  call NrmClc(Dlt(ipDlt),nBs,'DFroz','Dlt(ipDlt)')
# endif

  ipCff = ipCff+nBs*nOF
  ipDlt = ipDlt+lth
  ipOcc = ipOcc+nOr
end do

! Deallocate memory for new occupation numbers
call mma_deallocate(NewOcc)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine DFroz
