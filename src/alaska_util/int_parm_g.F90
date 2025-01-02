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

subroutine Int_Parm_g(iSD4,nSD,iPrimi,jPrimj,kPrimk,lPriml, &
                      nZeta,nEta,l2DI,nab,nHmab,ncd,nHmcd,nIrrep)

use k2_arrays, only: Create_BraKet
use Basis_Info, only: Shells
use Index_Functions, only: nTri_Elem1, nTri3_Elem1
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSD, iSD4(0:nSD,4), nIrrep
integer(kind=iwp), intent(out) :: iPrimi, jPrimj, kPrimk, lPriml, &
                                  nZeta, nEta, nab, nHmab, ncd, nHmcd
logical(kind=iwp), intent(in) :: l2DI
integer(kind=iwp) :: iAng, iCmp, jAng, jCmp, kAng, kCmp, lAng, lCmp

iPrimi = Shells(iSD4(0,1))%nExp
jPrimj = Shells(iSD4(0,2))%nExp
kPrimk = Shells(iSD4(0,3))%nExp
lPriml = Shells(iSD4(0,4))%nExp

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

nZeta = iPrimi*jPrimj
nEta = kPrimk*lPriml
call Create_BraKet(nZeta,nEta)

end subroutine Int_Parm_g
