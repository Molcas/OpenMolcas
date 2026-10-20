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
! Copyright (C) 2026, Lila Zapp                                        *
!***********************************************************************

subroutine GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,Gradient)
! Purpose: compute the gradient of the Pipek-Mezey functional w.r.t. elements of the kappa matrix

use Index_Functions, only: nTri_Elem
use Localisation_globals, only: Debug
use Constants, only: Zero, Four, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: GradNorm, Gradient(nTri_Elem(nOrb2Loc-1))
integer(kind=iwp) :: iAtom, k, kl, l
real(kind=wp) :: Q_kk, Q_kl, Q_ll

Q_ll = Zero
Q_kk = Zero
Q_kl = Zero

! gradient calculation according to DOI: 10.1002/jcc.23281 equation (15)
! the Gradient matrix is antisymmetric
kl = 0
Gradient(:) = Zero
do k=1,nOrb2Loc-1
  do l=k+1,nOrb2Loc
    kl = kl+1
    do iAtom=1,nAtoms
      Q_kk = PA(k,k,iAtom)
      Q_ll = PA(l,l,iAtom)
      Q_kl = PA(k,l,iAtom)
      Gradient(kl) = Gradient(kl)+(Q_kk-Q_ll)*Q_kl
      !write(u6,*) 'k,l,kl, Gradient(kl) = ',k,l,kl,Gradient(kl)
    end do
  end do
end do
Gradient(:) = Four*Gradient(:)

! If a gradient element is truly small this could be an indication of that the
! coordinate is symmetry breaking. Set it to zero to avoid symmetry breaking due to
! numerical noise.
kl = 0
do k=1,nOrb2Loc-1
  do l=k+1,nOrb2Loc
    kl = kl+1
    if (abs(Gradient(kl)) < 1.0e-12_wp) Gradient(kl) = Zero
  end do
end do

! GradientNorm
GradNorm = sqrt(dot_product(Gradient,Gradient))

if (Debug) then
  write(u6,*)
  write(u6,*) 'In GetGrad_PM'
  write(u6,*) '-------------'
  call RecPrt('Gradient',' ',Gradient(:)*Quart,size(Gradient),1)
  write(u6,*) 'gradnorm =',GradNorm
  write(u6,*)
end if

end subroutine GetGrad_PM
