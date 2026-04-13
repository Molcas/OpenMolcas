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

subroutine GetHdiag_PM(nAtoms,nOrb2Loc,PA,H_diag)
!
! Purpose: compute the Hessian diagonal elements of the Pipek-Mezey functional w.r.t. elements of the kappa matrix

use Constants, only: Zero, Four, Eight
use Definitions, only: wp, iwp, u6
use Localisation_globals, only: Debug

implicit none

integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: H_diag(nOrb2Loc*(nOrb2Loc-1)/2)
integer(kind=iwp) :: iAtom, k,l,kl
real(kind=wp) :: Q_ll, Q_kk, Q_kl, mat(nOrb2Loc,nOrb2Loc)

Q_ll = Zero
Q_kk = Zero
Q_kl = Zero

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
        !Write (u6,*) 'H_diag(k,l)=',H_diag(kl)
         H_diag(kl)=-1.0e-2_wp
      End If

   end do
end do

if (Debug) then
    write(u6,*) ' '
    write(u6,*) 'In GetHdiag_PM'
    write(u6,*) '-------------'
    call vec2upper_triag(mat,norb2loc,h_diag,(nOrb2Loc*(nOrb2Loc-1)/2),.false.)
    call RecPrt('H_diag',' ',mat(:,:), nOrb2Loc,nOrb2Loc)
    write(u6,*) ' '
end if

end subroutine GetHdiag_PM
