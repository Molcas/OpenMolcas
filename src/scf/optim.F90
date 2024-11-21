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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Optim(E_Pred,G,H,C,n)
!***********************************************************************
!                                                                      *
!    purpose: Solve a set of non-linear equations to get interpolation *
!             coefficients for damping                                 *
!                                                                      *
!     input:                                                           *
!       G       : linear terms                                         *
!       H       : bilinear terms                                       *
!       n       : size of matrices                                     *
!                                                                      *
!     output:                                                          *
!       C       : interpolation coefficients                           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! The KISS method is employed here. Do a scan of the energy surface    *
! along lines parallel to the edges. Reduce the steplength gradually.  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark                                                     *
!     University of Lund, Sweden, 2003                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(out) :: E_Pred, C(n)
real(kind=wp), intent(in) :: G(n), H(n,n)
integer(kind=iwp) :: i, Iter, j
real(kind=wp) :: Ci, Cj, Eminus, Eplus, Eref, rSum, Step, Step_m, Step_p
logical(kind=iwp) :: Debug, Debug2, DidChange
real(kind=wp), external :: Optim_E

!----------------------------------------------------------------------*
! Initialize                                                           *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
Debug = .true.
Debug2 = .true.
#else
Debug = .false.
Debug2 = .false.
#endif
!----------------------------------------------------------------------*
! Make first guess.                                                    *
!----------------------------------------------------------------------*
C(:) = Zero
j = 1
do i=1,n
  if (G(i)+Half*H(i,i) < G(j)+half*H(j,j)) j = i
end do
C(j) = 0.9_wp
do i=1,n
  if (i /= j) C(i) = (One-C(j))/(n-1)
end do
if (Debug) then
  write(u6,'(a)') 'Start C:'
  write(u6,'(6F15.6)') (C(i),i=1,n)
  write(u6,'(a)') 'Start G:'
  write(u6,'(6F26.16)') (G(i),i=1,n)
  call RecPrt('H',' ',H,n,n)
end if
!----------------------------------------------------------------------*
! Compute start energy.                                                *
!----------------------------------------------------------------------*
Eref = sum(C(:)*G(:))
do i=1,n
  Eref = Eref+Half*C(i)*sum(C(:)*H(i,:))
end do
if (Debug) write(u6,'(a,F24.16)') 'Eref = ',Eref
!----------------------------------------------------------------------*
! Do a scan in all pairwise directions                                 *
!----------------------------------------------------------------------*
Iter = 0
DidChange = .true.
Step = 0.1_wp
Eminus = One
Eplus = One

do while ((Iter < 500) .and. DidChange)

  Iter = Iter+1
  DidChange = .false.
  do i=1,n-1
    do j=i+1,n

      Ci = C(i)
      Cj = C(j)
      Step_p = min(Step,One-C(i),C(j))
      C(i) = C(i)+Step_p
      C(j) = C(j)-Step_p

      Eplus = optim_E(C,G,H,n)

      C(i) = Ci
      C(j) = Cj

      Step_m = min(Step,C(i),One-C(j))

      C(i) = C(i)-Step_m
      C(j) = C(j)+Step_m

      EMinus = optim_E(C,G,H,n)

      C(i) = Ci
      C(j) = Cj

      if (Debug2) then
        write(u6,'(a,2I3,3F26.16)') 'i,j,Eref,Eplus,Eminus',i,j,Eref,Eplus,Eminus
        write(u6,'(a,2F26.16)') 'Step_p, Step_m=',Step_p,Step_m
        write(u6,*) 'Eplus < EMinus=',Eplus < EMinus
        write(u6,*) 'Eplus < Eref=',Eplus < Eref
        write(u6,*) 'Eminus < Eref=',Eminus < Eref
      end if
      if (abs(Eplus-EMinus) > 1.0e-12_wp) then
        if (Eplus < Eminus) then
          if (Eplus < Eref) then
            C(i) = C(i)+Step_p
            C(j) = C(j)-Step_p
            Eref = Eplus
            DidChange = .true.
          end if
        else
          if (Eminus < Eref) then
            C(i) = C(i)-Step_m
            C(j) = C(j)+Step_m
            Eref = Eminus
            DidChange = .true.
          end if
        end if
      end if

    end do ! j
  end do ! i

  if (.not. DidChange) then
    if (Step > 0.9e-4_wp) then
      Step = 0.1_wp*Step
      DidChange = .true.
      if (Debug2) write(u6,*) 'Step is',Step
    end if
  end if

  ! Check that the constraint is fulfilled.

  rSum = Zero
  do i=1,n
    if (C(i) < Zero) C(i) = Zero
    if (C(i) > One) C(i) = One
    rSum = rSum+C(i)
  end do
# ifdef _DEBUGPRINT_
  write(u6,*) 'optim: rSum-1',rSum-One
# endif
  C(:) = C(:)/rSum
end do
#ifdef _DEBUGPRINT_
write(u6,*) 'ERef=',ERef
#endif
E_Pred = ERef

end subroutine Optim
