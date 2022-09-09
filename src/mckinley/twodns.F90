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
! Copyright (C) 1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine TwoDns(ianga,iCmp,shijij,ishll,ishell,iAO,nop,iBasi,jBasj,kBask,lBasl,Aux,nAux,Work2,nWork2,Work3,nWork3,Work4,nWork4, &
                  PSO,nPSO,Fact)

use Index_Functions, only: nTri_Elem1
use Basis_Info, only: Shells
use Real_Spherical, only: ipSph, RSph
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iAnga(4), iCmp(4), iShll(4), iShell(4), iAO(4), nOp(4), iBasi, jBasj, kBask, lBasl, nAux, nWork2, &
                                 nWork3, nWork4, nPSO
logical(kind=iwp), intent(in) :: Shijij
real(kind=wp), intent(out) :: Aux(nAux), Work2(nWork2), Work3(nWork3), Work4(nWork4)
real(kind=wp), intent(in) :: PSO(nPSO), Fact
integer(kind=iwp) :: iCmpa, ijklab, iShlla, jCmpb, jShllb, kCmpc, kShllc, la, lb, lc, lCmpd, ld, lShlld, mab, mcd, n, nijkl

nijkl = iBasi*jBasj*kBask*lBasl
iShlla = iShll(1)
jShllb = iShll(2)
kShllc = iShll(3)
lShlld = iShll(4)
la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)
iCmpa = iCmp(1)
jCmpb = iCmp(2)
kCmpc = iCmp(3)
lCmpd = iCmp(4)
mab = nTri_Elem1(la)*nTri_Elem1(lb)
mcd = nTri_Elem1(lc)*nTri_Elem1(ld)

!----------------------------------------------------------------------*
!
! Fix the second order density matrix
!
!----------------------------------------------------------------------*

! Desymmetrize the second order density matrix

! (faA fbR(B) | fcT(C) fdTS(D))ijkl
! PSO->Work2

call DesymP(iAnga,iCmp(1),iCmp(2),iCmp(3),iCmp(4),Shijij,iShll,iShell,iAO,nOp,nijkl,Aux,nAux,Work2,PSO,nPSO)

if (Fact /= One) then
  n = nijkl*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)
  Work2(1:n) = Fact*Work2(1:n)
end if

! Backtransform 2nd order density matrix from spherical
! harmonic gaussians to cartesian gaussians.

ijklab = nijkl*iCmp(1)*iCmp(2)
! Work2->Work2  (Work3:Scratch)
call SphCr1(Work2,ijklab,Work3,nWork3,RSph(ipSph(lc)),nTri_Elem1(lc),kCmpc,Shells(kShllc)%Transf,Shells(kShllc)%Prjct, &
            RSph(ipSph(ld)),nTri_Elem1(ld),lCmpd,Shells(lShlld)%Transf,Shells(lShlld)%Prjct,Work2,mcd)
! Work2->Work4  (Work3:Scratch)
call SphCr2(Work2,nijkl,mcd,Work3,nWork3,RSph(ipSph(la)),nTri_Elem1(la),iCmpa,Shells(iShlla)%Transf,Shells(iShlla)%Prjct, &
            RSph(ipSph(lb)),nTri_Elem1(lb),jCmpb,Shells(jShllb)%Transf,Shells(jShllb)%Prjct,Work4,mab)

!----------------------------------------------------------------------*

! P is now in cartesian AO base

return

end subroutine TwoDns
