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
subroutine GetHdiag_PM(nAtoms,nOrb2Loc,PA,H_diag,npos,gradnorm,modify)
!
! Purpose: compute the Hessian diagonal elements of the Pipek-Mezey functional w.r.t. elements of the kappa matrix

use Constants, only: Zero, Four
use Definitions, only: wp, iwp, u6
use Localisation_globals, only: Debug

implicit none

integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms),gradnorm
real(kind=wp), intent(out) :: H_diag(nOrb2Loc*(nOrb2Loc-1)/2)
logical(kind=iwp), intent(in) :: modify
integer(kind=iwp),intent(out) :: npos
integer(kind=iwp) :: iAtom, k,l,kl
real(kind=wp) :: Q_ll, Q_kk, Q_kl, Thr
logical(kind=iwp) :: SORange,prnt=.false.,prnt2=.false.

#ifdef _NOTUSED_
real(kind=wp) :: maxel
#endif

npos= 0

! set to false, if positive diagonal elements
SORange = .true.

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
!         H_diag(kl)=H_diag(kl) + Four*Q_ll*(Q_kk-Q_ll) + Four*Q_kk*(Q_ll-Q_kk) + Four*Four*Q_kl!**2
!          H_diag(kl)=H_diag(kl) + Four*(Q_ll*(Q_kk-Q_ll) + Q_kk*(Q_ll-Q_kk) + Four*Q_kl**2)
          H_diag(kl)=H_diag(kl) + Four*(-(Q_kk-Q_ll)**2 + Four*Q_kl**2)
      end do

      if (modify) then
          !write(u6,"(A,I5,I5,I5,3X,A,F18.8)") "k,l,kl = ",k,l,kl,"H_diag(kl)",H_diag(kl)
    !     Make sure that element has a negative value -- we are maximizing the target function
    !     Make sure that the element is not too small, this would yield a too large displacement.
          If (H_diag(kl)>Zero) Then
            npos = npos + 1
            if (prnt2) write(u6,*) "flip sign at",kl,'H_diag(kl)=',H_diag(kl)
            H_diag(kl)=-H_diag(kl)
            SORange = .false.
          End If
      end if

   end do
end do

if (modify) then
    if (SORange .or. gradnorm < 1.0e-2_wp) then
      ! higher trust in the hessian now and allow faster convergence
      thr = gradnorm
      if (prnt) write(u6,*) "in SORange: no positive diagonal elements"
    else
    ! outside of quadratic region: hessian not so accurate because of flipped signs
    ! higher threshold to avoid unjustified large steps
      thr = 4.0e-2_wp
      if (prnt) write(u6,*) "outside of SORange: multiple positive EVs"
    end if

    kl = 0
    do k=1,nOrb2Loc-1
       do l=k+1,nOrb2Loc
          kl = kl + 1 !listindex
              If (Abs(H_diag(kl))<Thr) Then
                 if (prnt2) write(u6,*) "lower limit  ",kl, 'H_diag(k,l)=',H_diag(kl)
                 H_diag(kl) = -Thr
              End If

       end do
    end do
end if !modify



#ifdef _NOTUSED_
maxel = maxval(H_diag)
if (maxel > Zero) Then
    !write(u6,*) "max(H_diag)",maxel
    !call RecPrt("H_diag before","",H_diag,nOrb2Loc*(nOrb2Loc-1)/2,1)
    H_diag(:) = H_diag(:)-maxel-1.0e-2_wp
    !call RecPrt("H_diag after","",H_diag,nOrb2Loc*(nOrb2Loc-1)/2,1)
end if
#endif

if (Debug) then
    write(u6,*) ' '
    write(u6,*) 'In GetHdiag_PM'
    write(u6,*) '-------------'
    call RecPrt("Hdiag","",H_diag,nOrb2Loc*(nOrb2Loc-1)/2,1)
    write(u6,*) ' '
end if

end subroutine GetHdiag_PM
