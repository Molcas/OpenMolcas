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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

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
!       Dlt     : density matrix in triangular storage                 *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use SpinAV, only: Do_SpinAV, DSc
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym), nFro(nSym), nCMO
real(kind=wp), target, intent(in) :: CMO(nCMO), Occ(*)
real(kind=wp), target, intent(_OUT_) :: Dlt(*)
logical(kind=iwp), intent(in) :: alpha_density
integer(kind=iwp) :: i, iCol, iDSc, iOff, iOffD, ipCMO, ipDlt, ipDScc, ipOcc, iROw, iSym, j, ji, jj, lOff, lth, nBs, nFr, nOr
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nOcc
#endif
real(kind=wp) :: rSum, xsign
real(kind=wp), pointer :: pCMO(:,:), pOcc(:), pDlt(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

#ifdef _DEBUGPRINT_
call NrmClc(CMO,nCMO,'DOne_SCF','CMO')
nOcc = sum(nOrb(:))
call NrmClc(Occ,nOcc,'DOne_SCF','Occ')
#endif

ipCMO = 1
ipDlt = 1
ipOcc = 1
do iSym=1,nSym

  nBs = nBas(iSym)
  nOr = nOrb(iSym)
  nFr = nFro(iSym)

  lth = nTri_Elem(nBs)
  Dlt(ipDlt:ipDlt+lth-1) = Zero

  if (nOr /= 0) then
    pCMO(1:nBs,1:nBs) => CMO(ipCMO:ipCMO+nBs**2-1)
    pOcc => Occ(ipOcc:ipOcc+nOr-1)
    pDlt => Dlt(ipDlt:ipDlt+lth-1)

    do iRow=1,nBs
      do iCol=1,iRow
        rSum = sum(pOcc(nFr+1:nOr)*pCMO(iRow,nFr+1:nOr)*pCMO(iCol,nFr+1:nOr))
        if (iCol == iRow) then
          pDlt(iTri(iRow,iRow)) = rSum
        else
          pDlt(iTri(iRow,iCol)) = Two*rSum
        end if
      end do
    end do
#   ifdef _DEBUGPRINT_
    call NrmClc(pDlt,nTri_Elem(nBs),'DOne_SCF','Dlt')
    !call RecPrt('CMO',' ',pCMO,nBs,nBs)
    !call RecPrt('Occ',' ',pOcc,1,nOr)
    !call TriPrt('Dlt',' ',pDlt,nBs)
#   endif
    nullify(pCMO,pOcc,pDlt)
  end if

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
    lth = nTri_Elem(nBas(iSym))

    pDlt => Dlt(iOffD:iOffD+lth-1)

    ipDScc = 1+lOff
    do j=1,nBas(iSym)
      do i=1,j-1
        ji = iTri(j,i)
        iDSc = ipDScc-1+nBas(iSym)*(j-1)+i
        pDlt(ji) = pDlt(ji)+xsign*Two*DSc(iDSc)
      end do
      jj = nTri_Elem(j)
      iDSc = ipDScc-1+nBas(iSym)*(j-1)+j
      pDlt(jj) = pDlt(jj)+xsign*DSc(iDSc)
    end do
    iOff = iOff+lth
    lOff = lOff+nBas(iSym)**2
  end do

  nullify(pDlt)

end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

#undef _DEBUGPRINT_
end subroutine DOne_SCF
