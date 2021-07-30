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

subroutine Int_Parm_g(iSD4,nSD,iAnga,iCmpa,iShlla,iShela,iPrimi,jPrimj,kPrimk,lPriml,indij,k2ij,nDCRR,k2kl,nDCRS,mdci,mdcj,mdck, &
                      mdcl,AeqB,CeqD,nZeta,nEta,ipZeta,ipZI,ipP,ipEta,ipEI,ipQ,ipiZet,ipiEta,ipxA,ipxB,ipxG,ipxD,l2DI,nab,nHmab, &
                      ncd,nHmcd,nIrrep)

use k2_setup, only: Indk2
use Basis_Info, only: Shells
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nSD, iSD4(0:nSD,4), iAnga(4), iCmpa(4), iShlla(4), iShela(4), iPrimi, jPrimj, kPrimk, lPriml, indij, k2ij, &
                     nDCRR, k2kl, nDCRS, mdci, mdcj, mdck, mdcl, nZeta, nEta, ipZeta, ipZI, ipP, ipEta, ipEI, ipQ, ipiZet, ipiEta, &
                     ipxA, ipxB, ipxG, ipxD, nab, nHmab, ncd, nHmcd, nIrrep
logical(kind=iwp) :: AeqB, CeqD, l2DI
integer(kind=iwp) :: iAng, iCmp, ijShll, iShell, jAng, jCmp, jShell, kAng, kCmp, klShll, kShell, lAng, lCmp, lShell
! Statement functions
integer(kind=iwp) :: nElem, ixyz, nabSz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

call ICopy(4,iSD4(1,1),nSD+1,iAnga,1)
call ICopy(4,iSD4(2,1),nSD+1,iCmpa,1)
call ICopy(4,iSD4(0,1),nSD+1,iShlla,1)
call ICopy(4,iSD4(11,1),nSD+1,iShela,1)
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

nab = nElem(iAng)*nElem(jAng)
ncd = nElem(kAng)*nElem(lAng)
nHmab = iCmp*jCmp*(nabSz(iAng+jAng)-nabSz(max(iAng,jAng)-1))
nHmab = nHmab*nIrrep
nHmcd = kCmp*lCmp*(nabSz(kAng+lAng)-nabSz(max(kAng,lAng)-1))
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
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(indij)
end if

end subroutine Int_Parm_g
