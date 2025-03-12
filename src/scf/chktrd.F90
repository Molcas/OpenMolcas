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
! Copyright (C) 1996, Martin Schuetz                                   *
!***********************************************************************

subroutine ChkTrD(nSym,nBas,nOrb,Occ,nOcc,Dlt,nDlt)
!***********************************************************************
!                                                                      *
!     purpose: Compute trace of density matrix and compare with sum    *
!              over occupation numbers...                              *
!                                                                      *
!     input:                                                           *
!       nSym    : number of symmetries                                 *
!       nBas(i) : number of basis functions (i = 1, nSym)              *
!       Occ     : occupation numbers                                   *
!       Dlt     : density matrix in triangular storage                 *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use InfSCF, only: Ovrlp
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym), nOcc, nDlt
real(kind=wp), intent(in) :: Occ(nOcc), Dlt(nDlt)
integer(kind=iwp) :: ipDlt, ipOcc, ipOvl, iSym, lth, nBs, nOr
real(kind=wp) :: Scal, SumOcc, TrDns
real(kind=wp), parameter :: ThrDif = 1.0e-7_wp
real(kind=wp), external :: DDot_

ipDlt = 1
ipOvl = 1
ipOcc = 0
Scal = One
do iSym=1,nSym
  nBs = nBas(iSym)
  if (nBs < 1) cycle
  nOr = nOrb(iSym)
  lth = nTri_Elem(nBs)
  ! count occupation number...
  SumOcc = sum(Occ(ipOcc+1:ipOcc+nOr))*Scal
  ! do trace of PS for symmetry block...
  TrDns = DDOT_(lth,Dlt(ipDlt),1,Ovrlp(ipOvl),1)
  ipDlt = ipDlt+lth
  ipOvl = ipOvl+lth
  ipOcc = ipOcc+nOr
  if (abs(SumOcc-TrDns) > ThrDif) then
    write(u6,*) abs(SumOcc-TrDns)
    ! print Warning...
    call WarningMessage(1,'WARNING: trace of density is inconsistent with occupation !')
    write(u6,'(A,I1,A,3F12.7)') 'SymBlock: ',iSym,' deviation: ',SumOcc-TrDns,SumOcc,TrDns
  end if
end do

end subroutine ChkTrD
