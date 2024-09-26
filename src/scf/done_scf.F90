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
!               2011, Francesco Aquilante                              *
!***********************************************************************

#include "compiler_features.h"
#ifdef _IN_MODULE_

!#define _DEBUGPRINT_
subroutine DOne_SCF(nSym,nBas,nOrb,nFro,CMO,nCMO,Occ,Dlt,alpha_density)
!***********************************************************************
!                                                                      *
!     purpose: Compute density matrix in AO basis                      *
!                                                                      *
!     input:                                                           *
!       nSym    : number of symmetries                                 *
!       nBas(i) : number of basis functions (i = 1, nSym)              *
!       nOrb(i) : number of orbitals (i = 1, nSym)                     *
!       nFro(i) : number of frozen orbitals (i = 1, nSym)              *
!       CMO     : molecular orbitals                                   *
!       Occ     : occupation numbers                                   *
!                                                                      *
!       alpha_density: .true. iff alpha MOs were sent in               *
!                                                                      *
!     output:                                                          *
!       Dlt     : density matrix in triangular storrage                *
!***********************************************************************

#include "compiler_features.h"

#ifndef POINTER_REMAP
use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
#endif
use SpinAV, only: Do_SpinAV, DSc
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nSym, nBas(nSym), nOrb(nSym), nFro(nSym), nCMO
real(kind=wp), target :: CMO(nCMO), Occ(*), Dlt(*)
logical(kind=iwp) :: alpha_density
integer(kind=iwp) :: i, iCol, iDSc, iOff, iOffD, ipCMO, ipDlt, ipDScc, ipOcc, iROw, iSym, j, ji, jj, lOff, lth, nBs, nFr, nOr
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nOcc
#endif
real(kind=wp) :: rSum, xsign
real(kind=wp), pointer :: pCMO(:,:), pOcc(:), pDlt(:)
! Statement function for triangular storage
integer(kind=iwp) :: Ind
Ind(i,j) = i*(i-1)/2+j

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

#ifdef _DEBUGPRINT_
call NrmClc(CMO,nCMO,'DOne_SCF','CMO')
nOcc = 0
do iSym=1,nSym
  nOcc = nOcc+nOrb(iSym)
end do
call NrmClc(Occ,nOcc,'DOne_SCF','Occ')
#endif

ipCMO = 1
ipDlt = 1
ipOcc = 1
do iSym=1,nSym

  nBs = nBas(iSym)
  nOr = nOrb(iSym)
  nFr = nFro(iSym)

  lth = nBs*(nBs+1)/2
  call FZero(Dlt(ipDlt),lth)
  if (nOr == 0) Go To 100

# ifdef POINTER_REMAP
  pCMO(1:nBs,1:nBs) => CMO(ipCMO:ipCMO+nBs**2-1)
# else
  call c_f_pointer(c_loc(CMO(ipCMO)),pCMO,[nBs,nBs])
# endif
  pOcc => Occ(ipOcc:ipOcc+nOr-1)
  pDlt => Dlt(ipDlt:ipDlt+lth-1)

  do iRow=1,nBs
    rSum = Zero
    do i=nFr+1,nOr
      rSum = rSum+pOcc(i)*pCMO(iRow,i)**2
      !rSum = rSum+pOcc(i)*pCMO(iRow,i)*pCMO(iRow,i)
    end do
    pDlt(Ind(iRow,iRow)) = rSum

    do iCol=1,iRow-1
      rSum = Zero
      do i=nFr+1,nOr
        rSum = rSum+pOcc(i)*pCMO(iRow,i)*pCMO(iCol,i)
      end do
      pDlt(Ind(iRow,iCol)) = Two*rSum
    end do
  end do
# ifdef _DEBUGPRINT_
  call NrmClc(pDlt,nBs*(nBs+1)/2,'DOne_SCF','Dlt')
  !call RecPrt('CMO',' ',pCMO,nBs,nBs)
  !call RecPrt('Occ',' ',pOcc,1,nOr)
  !call TriPrt('Dlt',' ',pDlt,nBs)
# endif

100 continue

  nullify(pCMO,pOcc,pDlt)
  ipCMO = ipCMO+nBs**2
  ipDlt = ipDlt+lth
  ipOcc = ipOcc+nOr
end do

xsign = One
if (Do_SpinAV) then

  if (alpha_density) xsign = -One
  iOffD = 1
  lOff = 0
  do iSym=1,nSym
    lth = nBas(iSym)*(nBas(iSym)+1)/2

    pDlt => Dlt(iOffD:iOffD+lth-1)

    ipDScc = 1+lOff
    do j=1,nBas(iSym)
      do i=1,j-1
        ji = j*(j-1)/2+i
        iDSc = ipDScc-1+nBas(iSym)*(j-1)+i
        pDlt(ji) = pDlt(ji)+xsign*Two*DSc(iDSc)
      end do
      jj = j*(j+1)/2
      iDSc = ipDScc-1+nBas(iSym)*(j-1)+j
      pDlt(jj) = pDlt(jj)+xsign*DSc(iDSc)
    end do
    iOff = iOff+lth
    lOff = lOff+nBas(iSym)**2
  end do

end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine DOne_SCF

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(DOne_SCF)

#endif
