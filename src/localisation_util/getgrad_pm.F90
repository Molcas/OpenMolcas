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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,Gradient,H_diag)
! LilaZapp, December 2025/26.
!
! Purpose: compute the gradient of the Pipek-Mezey functional.

use Constants, only: Zero, Four, Eight
use Definitions, only: wp, iwp, u6
use Localisation_globals, only: Debug

implicit none

integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: GradNorm,Gradient(nOrb2Loc*(nOrb2Loc-1)/2),H_diag(nOrb2Loc*(nOrb2Loc-1)/2)
integer(kind=iwp) :: iAtom, k,l,kl
real(kind=wp) :: Q_ll, Q_kk, Q_kl

!gradient and Hessian - needed only for new optimizer
Q_ll = Zero
Q_kk = Zero
Q_kl = Zero

!gradient calculation according to DOI: 10.1002/jcc.23281 equation (15)
!the Gradient matrix is antisymmetric
kl = 0
Gradient(:) = Zero
do k=1,nOrb2Loc-1
    do l=k+1,nOrb2Loc
        kl = kl + 1
        do iAtom=1,nAtoms
            Q_kk=PA(k,k,iAtom)
            Q_ll=PA(l,l,iAtom)
            Q_kl=PA(k,l,iAtom)
            Gradient(kl)=Gradient(kl)+(Q_kk-Q_ll)*Q_kl
            !write(u6,*) "k,l,kl, Gradient(kl) = ",k,l,kl,Gradient(kl)
        end do
    end do
end do
Gradient(:)=Four*Gradient(:)

!Second derivative for GEK optimization: Later put this into an "if GEK=true" environment
!Hessian diagonal according to DOI: 10.1002/jcc.23281 equation (17)
H_diag(:) = Zero
kl = 0
do k=1,nOrb2Loc-1
   do l=k+1,nOrb2Loc
      kl = kl + 1 !listindex
      do iAtom=1,nAtoms
          Q_kk=PA(k,k,iAtom)
          Q_ll=PA(l,l,iAtom)
          Q_kl=PA(k,l,iAtom)
          H_diag(kl)=H_diag(kl) + Four*Q_ll*(Q_kk-Q_ll) + Four*Q_kk*(Q_ll-Q_kk) + Eight*Q_kl**2
      end do
      !write(u6,"(A,I5,I5,I5,3X,A,F18.8)") "k,l,kl = ",k,l,kl,"H_diag(kl)",H_diag(kl)
!     Make sure that element has a negative value -- we are maximizing the target function
!     Make sure that the element is not too small, this would yield a too large displacement.
      If (H_diag(kl)>0.0) Then
!        Write (*,*) 'H_diag(k,l)=',H_diag(k,l)
         H_diag(kl)=-H_diag(kl)
      End If
      If (Abs(H_diag(kl))<1.0e-2_wp) Then
!        Write (*,*) 'H_diag(k,l)=',H_diag(k,l)
         H_diag(kl)=-1.0e-2_wp
      End If

   end do
end do

!GradientNorm - needed for all optimization schemes,
! here calculated as the vector norm
kl=0
GradNorm = Zero
do k = 1, nOrb2Loc-1
    do l = k+1, nOrb2Loc
        kl = kl+1
        GradNorm = GradNorm + Gradient(kl)**2
    end do
end do
GradNorm = sqrt(GradNorm)

if (Debug) then
    write(u6,*) ' '
    write(u6,*) 'In GetGrad_PM'
    write(u6,*) '-------------'
    call RecPrt('Gradient',' ',Gradient(:), nOrb2Loc*(nOrb2Loc-1)/2,1)  !this is also printed in the Gradientlist
    call RecPrt('H_diag',' ',H_diag(:), nOrb2Loc*(nOrb2Loc-1)/2,1)
    write(u6,*) "gradnorm =",GradNorm
    write(u6,*) ' '
end if

end subroutine GetGrad_PM
