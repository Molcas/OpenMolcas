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

subroutine mkangLmin(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,angintSO,angintOO,Lfirst,Llast,Lblocks,ncont1,ncont2,ncont3,ncont4,caseaSO, &
                     caseb1SO,caseb2SO,casecSO,caseaOO,caseb1OO,caseb2OO,casecOO,preroots,clebsch,dummy,bonn,breit,sameorb)
!bs subroutine for combining radial intgrls with angular
!bs factors for the block with l1,l2,l3,l4,m1,m2,m3m,m4
!bs this routine mkangLmin = make angular factors for the L- -part
!bs includes both, spin-same and spin-other-orbit parts.

use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lmax, l1, l2, l3, l4, m1, m2, m3, m4, Lfirst(4), Llast(4), Lblocks(4), ncont1, ncont2, ncont3, &
                                 ncont4
real(kind=wp), intent(out) :: angintSO(ncont1,ncont2,ncont3,ncont4), angintOO(ncont1,ncont2,ncont3,ncont4), dummy(0:2*Lmax+1)
real(kind=wp), intent(in) :: caseaSO(ncont1*ncont2*ncont3*ncont4,*), caseb1SO(ncont1*ncont2*ncont3*ncont4,*), &
                             caseb2SO(ncont1*ncont2*ncont3*ncont4,*), casecSO(ncont1*ncont2*ncont3*ncont4,*), &
                             caseaOO(ncont1*ncont2*ncont3*ncont4,*), caseb1OO(ncont1*ncont2*ncont3*ncont4,*), &
                             caseb2OO(ncont1*ncont2*ncont3*ncont4,*), casecOO(ncont1*ncont2*ncont3*ncont4,*), preroots(2,0:Lmax), &
                             clebsch(3,2,-Lmax:Lmax,0:Lmax)
logical(kind=iwp), intent(in) :: bonn, breit, sameorb
integer(kind=iwp) :: Kfirst, Klast, L, Lrun, M, ncontall
real(kind=wp) :: cheater, factor
real(kind=wp), parameter :: root2 = sqrt(Two), root2inv = One/root2
real(kind=wp), external :: LMdepang
!bs all the arrays with the radial intgrls for
!bs this combination of l-values
! caseaSO:  (2,0)   intgrls with alpha1*alpha3
! caseb1SO: (0,0)   intgrls with alpha1
! caseb2SO: (0,0)   intgrls with alpha3
! casecSO:  (-2,0)  intgrls with factor 1
! caseaOO:  (2,0)   intgrls with alpha1*alpha3
! caseb1OO: (0,0)   intgrls with alpha1
! caseb2OO: (0,0)   intgrls with alpha3
! casecOO:  (-2,0)  intgrls with factor 1
! preroots: some prefactors: sqrt( (l(+1))/(2l+1))
! clebsch:  some clebsch gordans, that appear regulary

!write(u6,*) 'begin mkangL- ',l1,l2,l3,l4,m1,m2,m3,m4

ncontall = ncont1*ncont2*ncont3*ncont4
!bs cheater introduced to correct signs, because they were different from HERMIT
if (mod(l1+l2+l3+l4,4) == 2) then
  cheater = One
else
  cheater = -One
end if
!bs cleaning up
if (bonn .or. breit .or. sameorb) then
  angintSO(:,:,:,:) = Zero
else
  angintSO(:,:,:,:) = Zero
  angintOO(:,:,:,:) = Zero
end if
!bs starting with the same-orbit-contributions
!bs first term: ########################################################
factor = -root2inv*preroots(2,l1)*preroots(2,l3)*clebsch(3,2,m1,l1)*clebsch(2,2,m3,l3)
if (factor /= Zero) then
  dummy(:) = Zero
  !bs get the L,M dependent coefficients
  if (Lblocks(1) > 0) then
    M = m2-m4
    Lrun = 1
    do L=Lfirst(1),Llast(1),2
      dummy(L) = LMdepang(L,M,l1+1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,Four*factor*dummy(L),caseaOO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
end if
!bs second term: #######################################################
factor = -root2inv*preroots(1,l1)*preroots(2,l3)*clebsch(3,1,m1,l1)*clebsch(2,2,m3,l3)
if (factor /= Zero) then
  dummy(:) = Zero
  Klast = 0
  Kfirst = 2*Lmax+1 ! just to be sure ..
  !bs get the L,M dependent coefficients
  if (Lblocks(1) > 0) then
    M = m2-m4
    Kfirst = Lfirst(1)
    Klast = Llast(1)
    Lrun = 1
    do L=Lfirst(1),Llast(1),2
      dummy(L) = LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,Four*factor*dummy(L),caseaOO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
  if (Lblocks(3) > 0) then
    M = m2-m4
    if (Lfirst(3) < Kfirst) then
      do L=Lfirst(3),Kfirst,2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
      end do
      Kfirst = Lfirst(3)
    end if
    if (Llast(3) > Klast) then
      do L=Klast,Llast(3),2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
      end do
      Klast = Llast(3)
    end if
    Lrun = 1
    do L=Lfirst(3),Llast(3),2
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2SO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2SO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2OO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
end if
!bs third term: ########################################################
factor = -root2inv*preroots(2,l1)*preroots(1,l3)*clebsch(3,2,m1,l1)*clebsch(2,1,m3,l3)
if (factor /= Zero) then
  dummy(:) = Zero
  Klast = 0
  Kfirst = 2*Lmax+1 ! just to be sure ..
  !bs get the L,M dependent coefficients
  if (Lblocks(1) > 0) then
    M = m2-m4
    Kfirst = Lfirst(1)
    Klast = Llast(1)
    Lrun = 1
    do L=Lfirst(1),Llast(1),2
      dummy(L) = LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,Four*factor*dummy(L),caseaOO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
  if (Lblocks(2) > 0) then
    M = m2-m4
    if (Lfirst(2) < Kfirst) then
      do L=Lfirst(2),Kfirst,2
        dummy(L) = LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3,m4,Cheater)
      end do
      Kfirst = Lfirst(2)
    end if
    if (Llast(2) > Klast) then
      do L=Klast,Llast(2),2
        dummy(L) = LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      end do
      Klast = Llast(2)
    end if
    Lrun = 1
    do L=Lfirst(2),Llast(2),2
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1SO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1SO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1OO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
end if
!bs fourth term: #######################################################
factor = -root2inv*preroots(1,l1)*preroots(1,l3)*clebsch(3,1,m1,l1)*clebsch(2,1,m3,l3)
if (factor /= Zero) then
  dummy(:) = Zero
  Klast = 0
  Kfirst = 2*Lmax+1 ! just to be sure ..
  !bs get the L,M dependent coefficients
  if (Lblocks(1) > 0) then
    M = m2-m4
    Kfirst = Lfirst(1)
    Klast = Llast(1)
    Lrun = 1
    do L=Lfirst(1),Llast(1),2
      dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,Four*factor*dummy(L),caseaOO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
  if (Lblocks(2) > 0) then
    M = m2-m4
    if (Lfirst(2) < Kfirst) then
      do L=Lfirst(2),Kfirst,2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      end do
      Kfirst = Lfirst(2)
    end if
    if (Llast(2) > Klast) then
      do L=Klast,Llast(2),2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      end do
      Klast = Llast(2)
    end if
    Lrun = 1
    do L=Lfirst(2),Llast(2),2
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1SO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1SO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1OO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
  if (Lblocks(3) > 0) then
    M = m2-m4
    if (Lfirst(3) < Kfirst) then
      do L=Lfirst(3),Kfirst,2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      end do
      Kfirst = Lfirst(3)
    end if
    if (Llast(3) > Klast) then
      do L=Klast,Llast(3),2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      end do
      Klast = Llast(3)
    end if
    Lrun = 1
    do L=Lfirst(3),Llast(3),2
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2SO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2SO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2OO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
  if (Lblocks(4) > 0) then
    M = m2-m4
    if (Lfirst(4) < Kfirst) then
      do L=Lfirst(4),Kfirst,2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      end do
      Kfirst = Lfirst(4)
    end if
    if (Llast(4) > Klast) then
      do L=Klast,Llast(4),2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      end do
      Klast = Llast(4)
    end if
    Lrun = 1
    do L=Lfirst(4),Llast(4),2
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,real(4*l1*l3+2*l1+2*l3+1,kind=wp)*factor*dummy(L),casecSO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,real(4*l1*l3+2*l1+2*l3+1,kind=wp)*factor*dummy(L),casecSO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,real(4*l1*l3+2*l1+2*l3+1,kind=wp)*factor*dummy(L),casecOO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
end if
!bs fifth term: ########################################################
factor = -root2inv*preroots(2,l1)*preroots(2,l3)*clebsch(2,2,m1,l1)*clebsch(1,2,m3,l3)
if (factor /= Zero) then
  dummy(:) = Zero
  !bs get the L,M dependent coefficients
  if (Lblocks(1) > 0) then
    M = m2-m4
    Lrun = 1
    do L=Lfirst(1),Llast(1),2
      dummy(L) = LMdepang(L,M,l1+1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,Four*factor*dummy(L),caseaOO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
end if
!bs sixth term: ########################################################
factor = -root2inv*preroots(1,l1)*preroots(2,l3)*clebsch(2,1,m1,l1)*clebsch(1,2,m3,l3)
if (factor /= Zero) then
  dummy(:) = Zero
  Klast = 0
  Kfirst = 2*Lmax+1 ! just to be sure ..
  !bs get the L,M dependent coefficients
  if (Lblocks(1) > 0) then
    M = m2-m4
    Kfirst = Lfirst(1)
    Klast = Llast(1)
    Lrun = 1
    do L=Lfirst(1),Llast(1),2
      dummy(L) = LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,Four*factor*dummy(L),caseaOO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
  if (Lblocks(3) > 0) then
    M = m2-m4
    if (Lfirst(3) < Kfirst) then
      do L=Lfirst(3),Kfirst,2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      end do
      Kfirst = Lfirst(3)
    end if
    if (Llast(3) > Klast) then
      do L=Klast,Llast(3),2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      end do
      Klast = Llast(3)
    end if
    Lrun = 1
    do L=Lfirst(3),Llast(3),2
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2SO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2SO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2OO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
end if
!bs seventh term: ######################################################
factor = -root2inv*preroots(2,l1)*preroots(1,l3)*clebsch(2,2,m1,l1)*clebsch(1,1,m3,l3)
if (factor /= Zero) then
  dummy(:) = Zero
  Klast = 0
  Kfirst = 2*Lmax+1 ! just to be sure ..
  !bs get the L,M dependent coefficients
  if (Lblocks(1) > 0) then
    M = m2-m4
    Kfirst = Lfirst(1)
    Klast = Llast(1)
    Lrun = 1
    do L=Lfirst(1),Llast(1),2
      dummy(L) = LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,Four*factor*dummy(L),caseaOO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
  if (Lblocks(2) > 0) then
    M = m2-m4
    if (Lfirst(2) < Kfirst) then
      do L=Lfirst(2),Kfirst,2
        dummy(L) = LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      end do
      Kfirst = Lfirst(2)
    end if
    if (Llast(2) > Klast) then
      do L=Klast,Llast(2),2
        dummy(L) = LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      end do
      Klast = Llast(2)
    end if
    Lrun = 1
    do L=Lfirst(2),Llast(2),2
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1SO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1SO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1OO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
end if
!bs eighth term: #######################################################
factor = -root2inv*preroots(1,l1)*preroots(1,l3)*clebsch(2,1,m1,l1)*clebsch(1,1,m3,l3)
if (factor /= Zero) then
  dummy(:) = Zero
  Klast = 0
  Kfirst = 2*Lmax+1 ! just to be sure ..
  !bs get the L,M dependent coefficients
  if (Lblocks(1) > 0) then
    M = m2-m4
    Kfirst = Lfirst(1)
    Klast = Llast(1)
    Lrun = 1
    do L=Lfirst(1),Llast(1),2
      dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,Four*factor*dummy(L),caseaSO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,Four*factor*dummy(L),caseaOO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
  if (Lblocks(2) > 0) then
    M = m2-m4
    if (Lfirst(2) < Kfirst) then
      do L=Lfirst(2),Kfirst,2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      end do
      Kfirst = Lfirst(2)
    end if
    if (Llast(2) > Klast) then
      do L=Klast,Llast(2),2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      end do
      Klast = Llast(2)
    end if
    Lrun = 1
    do L=Lfirst(2),Llast(2),2
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1SO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1SO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,-real(2+4*l3,kind=wp)*factor*dummy(L),caseb1OO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
  if (Lblocks(3) > 0) then
    M = m2-m4
    if (Lfirst(3) < Kfirst) then
      do L=Lfirst(3),Kfirst,2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      end do
      Kfirst = Lfirst(3)
    end if
    if (Llast(3) > Klast) then
      do L=Klast,Llast(3),2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      end do
      Klast = Llast(3)
    end if
    Lrun = 1
    do L=Lfirst(3),Llast(3),2
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2SO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2SO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,-real(2+4*l1,kind=wp)*factor*dummy(L),caseb2OO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
  if (Lblocks(4) > 0) then
    M = m2-m4
    if (Lfirst(4) < Kfirst) then
      do L=Lfirst(4),Kfirst,2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      end do
      Kfirst = Lfirst(4)
    end if
    if (Llast(4) > Klast) then
      do L=Klast,Llast(4),2
        dummy(L) = LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      end do
      Klast = Llast(4)
    end if
    Lrun = 1
    do L=Lfirst(4),Llast(4),2
      if (dummy(L) /= Zero) then
        if (bonn .or. breit .or. sameorb) then
          call daxpy_(ncontall,real(4*l1*l3+2*l1+2*l3+1,kind=wp)*factor*dummy(L),casecSO(:,Lrun),1,angintSO,1)
        else
          call daxpy_(ncontall,real(4*l1*l3+2*l1+2*l3+1,kind=wp)*factor*dummy(L),casecSO(:,Lrun),1,angintSO,1)
          call daxpy_(ncontall,real(4*l1*l3+2*l1+2*l3+1,kind=wp)*factor*dummy(L),casecOO(:,Lrun),1,angintOO,1)
        end if
      end if
      Lrun = Lrun+1
    end do
  end if
end if

return

end subroutine mkangLmin
