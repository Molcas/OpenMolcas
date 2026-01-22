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

subroutine GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,Rmat,Gradient,H_diag)
! Thomas Bondo Pedersen, December 2005.
!
! Purpose: compute the gradient of the Pipek-Mezey functional.

use Constants, only: Zero, Four, Eight
use Definitions, only: wp, iwp, u6
use Localisation_globals, only: Debug

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: GradNorm, Rmat(nOrb2Loc,nOrb2Loc), Gradient(nOrb2Loc, nOrb2Loc), &
                              H_diag(nOrb2Loc, nOrb2Loc)
integer(kind=iwp) :: iAtom, i,j,k,l
real(kind=wp) :: Fun, Rjj, Q_ll, Q_kk, Q_kl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gradient and Hessian - needed only for new optimizer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Q_ll = Zero
Q_kk = Zero
Q_kl = Zero

!New gradient calculation according to DOI: 10.1002/jcc.23281 equation (15)
!the Gradient matrix is antisymmetric
Gradient(:,:) = Zero
do iAtom=1,nAtoms
    do k=1,nOrb2Loc
        Q_kk=PA(k,k,iAtom)
        do l=1,nOrb2Loc
            Q_ll=PA(l,l,iAtom)
            Q_kl=PA(k,l,iAtom)
            Gradient(k,l)=Gradient(k,l)+(Q_kk-Q_ll)*Q_kl
        end do
    end do
end do
Gradient(:,:)=Four*Gradient(:,:)

!Second derivative for GEK optimization: Later put this into an "if GEK=true" environment
!Hessian diagonal according to DOI: 10.1002/jcc.23281 equation (17)
H_diag(:,:) = Zero
do k=1,nOrb2Loc
   do l=1,nOrb2Loc
      do iAtom=1,nAtoms
          Q_kk=PA(k,k,iAtom)
          Q_ll=PA(l,l,iAtom)
          Q_kl=PA(k,l,iAtom)
          H_diag(k,l)=H_diag(k,l) + Four*Q_ll*(Q_kk-Q_ll) + Four*Q_kk*(Q_ll-Q_kk) + Eight*Q_kl**2
      end do
!     Make sure that element has a negative value -- we are maximizing the target function
!     Make sure that the element is not too small, this would yield a too large displacement.
      If (H_diag(k,l)>0.0) Then
!        Write (*,*) 'H_diag(k,l)=',H_diag(k,l)
         H_diag(k,l)=-H_diag(k,l)
      End If
      If (Abs(H_diag(k,l))<1.0e-2_wp) Then
!        Write (*,*) 'H_diag(k,l)=',H_diag(k,l)
         H_diag(k,l)=-1.0e-2_wp
      End If
   end do
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GradientNorm - needed for all optimization schemes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate the R matrix used to calculate the gradient norm
RMat(:,:) = Zero
do iAtom=1,nAtoms
  do j=1,nOrb2Loc
    Rjj = PA(j,j,iAtom)
    do i=1,nOrb2Loc
      Rmat(i,j) = Rmat(i,j)+PA(i,j,iAtom)*Rjj
    end do
  end do
end do

!The Gradient Norm (always positive!) is used as threshold criterium
GradNorm = Zero
do i=1,nOrb2Loc-1
  do j=i+1,nOrb2Loc
    GradNorm = GradNorm + (Rmat(i,j) - Rmat(j,i))**2 !square to avoid negative norm
  end do
end do
GradNorm = Four*sqrt(GradNorm) !sqrt to complement previous **2



if (Debug) then
    write(u6,*) ' '
    write(u6,*) 'In GetGrad_PM'
    write(u6,*) '-------------'
    call RecPrt('Gradient',' ',Gradient(:,:), nOrb2Loc, nOrb2Loc)  !this is also printed in the Gradientlist
    call RecPrt('H_diag',' ',H_diag(:,:), nOrb2Loc, nOrb2Loc)
    write(u6,*) ' '
    write(u6,*) 'GradNorm = ',gradnorm
    ! Trace of Rmat is the value of the PM functional
    Fun=Zero
    do i=1,nOrb2Loc
        Fun = Fun+Rmat(i,i)
    end do
    !call RecPrt('RMat',' ',RMat(:,:), nOrb2Loc, nOrb2Loc)
    !write(u6,*) 'PM_functional = Tr(R) = ',Fun
end if



end subroutine GetGrad_PM
