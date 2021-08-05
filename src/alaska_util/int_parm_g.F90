!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Int_Parm_g(iSD4,nSD,iAnga,iCmpa,iShlla,iShela,iPrimi,jPrimj,kPrimk,lPriml,k2ij,nDCRR,k2kl,nDCRS,mdci,mdcj,mdck,mdcl, &
                      AeqB,CeqD,nZeta,nEta,ipZeta,ipZI,ipP,ipEta,ipEI,ipQ,ipiZet,ipiEta,ipxA,ipxB,ipxG,ipxD,l2DI,nab,nHmab,ncd, &
                      nHmcd,nIrrep)

use k2_setup, only: Indk2
use Basis_Info, only: Shells
use Index_Functions, only: nTri_Elem1, nTri3_Elem1
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSD, iSD4(0:nSD,4), ipZeta, ipiZet, nIrrep
integer(kind=iwp), intent(out) :: iAnga(4), iCmpa(4), iShlla(4), iShela(4), iPrimi, jPrimj, kPrimk, lPriml, k2ij, nDCRR, k2kl, &
                                  nDCRS, mdci, mdcj, mdck, mdcl, nZeta, nEta, ipZI, ipP, ipEta, ipEI, ipQ, ipiEta, ipxA, ipxB, &
                                  ipxG, ipxD, nab, nHmab, ncd, nHmcd
logical(kind=iwp), intent(out) :: AeqB, CeqD
logical(kind=iwp), intent(in) :: l2DI
integer(kind=iwp) :: iAng, iCmp, ijShll, iShell, jAng, jCmp, jShell, kAng, kCmp, klShll, kShell, lAng, lCmp, lShell

iAnga(:) = iSD4(1,:)
iCmpa(:) = iSD4(2,:)
iShlla(:) = iSD4(0,:)
iShela(:) = iSD4(11,:)
iPrimi = Shells(iSD4(0,1))%nExp
jPrimj = Shells(iSD4(0,2))%nExp
kPrimk = Shells(iSD4(0,3))%nExp
lPriml = Shells(iSD4(0,4))%nExp
iShell = iSD4(11,1)
jShell = iSD4(11,2)
kShell = iSD4(11,3)
lShell = iSD4(11,4)
if (iShell >= jShell) then
  ijShll = iShell*(iShell-1)/2+jShell
else
  ijShll = jShell*(jShell-1)/2+iShell
end if
if (kShell >= lShell) then
  klShll = kShell*(kShell-1)/2+lShell
else
  klShll = lShell*(lShell-1)/2+kShell
end if
iAng = iSD4(1,1)
jAng = iSD4(1,2)
kAng = iSD4(1,3)
lAng = iSD4(1,4)
iCmp = iSD4(2,1)
jCmp = iSD4(2,2)
kCmp = iSD4(2,3)
lCmp = iSD4(2,4)

nab = nTri_Elem1(iAng)*nTri_Elem1(jAng)
ncd = nTri_Elem1(kAng)*nTri_Elem1(lAng)
nHmab = iCmp*jCmp*(nTri3_Elem1(iAng+jAng)-nTri3_Elem1(max(iAng,jAng)-1))
nHmab = nHmab*nIrrep
nHmcd = kCmp*lCmp*(nTri3_Elem1(kAng+lAng)-nTri3_Elem1(max(kAng,lAng)-1))
nHmcd = nHmcd*nIrrep
if (.not. l2DI) then
  nab = 0
  ncd = 0
end if
k2ij = Indk2(1,ijShll)
nDCRR = Indk2(2,ijShll)
k2kl = Indk2(1,klShll)
nDCRS = Indk2(2,klShll)
mdci = iSD4(10,1)
mdcj = iSD4(10,2)
mdck = iSD4(10,3)
mdcl = iSD4(10,4)
AeqB = (iSD4(13,1) == iSD4(13,2)) .and. (mdci == mdcj)
CeqD = (iSD4(13,3) == iSD4(13,4)) .and. (mdck == mdcl)
nZeta = iPrimi*jPrimj
nEta = kPrimk*lPriml
ipZI = ipZeta+nZeta
ipP = ipZI+nZeta
ipxA = ipP+nZeta*3
ipxB = ipxA+nZeta
ipEta = ipxB+nZeta
ipEI = ipEta+nEta
ipQ = ipEI+nEta
ipxG = ipQ+nEta*3
ipxD = ipxG+nEta

ipiEta = ipiZet+nZeta+1

return

end subroutine Int_Parm_g
