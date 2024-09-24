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

use InfSCF, only: nBas, nnB, nSym, nFro, nOrb
use Constants, only: Zero, One, Two
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer nDlt, nCMO
real*8 Dlt(nDlt), CMO(nCMO)
integer OccNo(*)
integer iOrb, iStrtN, iStrtO, iSym
real*8, dimension(:), allocatable :: NewOcc

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
call DOne_SCF_froz(CMO,nCMO,NewOcc,nnB)

! Deallocate memory for new occupation numbers
call mma_deallocate(NewOcc)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

contains

subroutine DOne_SCF_froz(Cff,nCff,Occ,nnB)
  !*********************************************************************
  !                                                                    *
  !   purpose: Compute density matrix in AO basis                      *
  !                                                                    *
  !   input:                                                           *
  !     Cff     : molecular orbitals                                   *
  !     Occ     : occupation numbers                                   *
  !                                                                    *
  !   output:                                                          *
  !     Dlt     : density matrix in triangular storrage                *
  !                                                                    *
  !   called from: DFroz                                               *
  !                                                                    *
  !*********************************************************************

  implicit none
  integer nCff, nnB
  real*8 Cff(nCff), Occ(nnB)
  integer i, j, iCol, ij, ipCff, ipDlt, ipOcc, iRow, lth, nBs, nFr, nOr, nOF, iSym
  integer Ind
  real*8 Scale, SScale, Sum
  ! Statement function for triangular storrage
  Ind(i,j) = i*(i-1)/2+j

  !--------------------------------------------------------------------*
  !     Start                                                          *
  !--------------------------------------------------------------------*

  Scale = Two
  SScale = One
  ipCff = 0
  ipDlt = 0
  ipOcc = 0
  do iSym=1,nSym

    nBs = nBas(iSym)
    nOr = nOrb(iSym)
    nFr = nFro(iSym)

    nOF = nOr-nFr
    lth = nBs*(nBs+1)/2

    ipCff = ipCff+nBs*nFr
    do iRow=1,nBs
      Sum = Zero
      ij = -1
      do i=nFr+1,nOr
        ij = ij+1
        !if (Occ(ipOcc+i) == Zero) Go To 100
        Sum = Sum+Occ(ipOcc+i)*Cff(ipCff+iRow+ij*nBs)*Cff(ipCff+iRow+ij*nBs)
      end do
      !100 continue
      Dlt(ipDlt+Ind(iRow,iRow)) = Sum*SScale

      do iCol=1,iRow-1
        Sum = Zero
        ij = -1
        do i=nFr+1,nOr
          ij = ij+1
          !if (Occ(ipOcc+i) == Zero) Go To 200
          Sum = Sum+Occ(ipOcc+i)*Cff(ipCff+iRow+ij*nBs)*Cff(ipCff+iCol+ij*nBs)
        end do
        !200 continue
        Dlt(ipDlt+Ind(iRow,iCol)) = Scale*Sum
      end do
    end do
#   ifdef _DEBUGPRINT_
    call NrmClc(Dlt(ipDlt),nBs,'DOne_SCF_froz','Dlt(ipDlt)')
#   endif

    ipCff = ipCff+nBs*nOF
    ipDlt = ipDlt+lth
    ipOcc = ipOcc+nOr
  end do

  !--------------------------------------------------------------------*
  !     Exit                                                           *
  !--------------------------------------------------------------------*

  return

end subroutine DOne_SCF_froz

end subroutine DFroz
