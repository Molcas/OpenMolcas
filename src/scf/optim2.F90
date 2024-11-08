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
!               2017, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Optim2(E_Pred,G,H,D,C,n,n0,n1,r2)
!***********************************************************************
!                                                                      *
!    purpose: Solve a set of non-linear equations to get interpolation *
!             coefficients for damping                                 *
!                                                                      *
!     input:                                                           *
!       G       : linear energy terms                                  *
!       H       : bilinear energy terms                                *
!       D       : bilinear density terms                               *
!       n       : size of matrices                                     *
!       n0      : index number of the lowest energy                    *
!       n1      : index number of the second lowest energy             *
!       r2      : square of density change constraint                  *
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
!     hacked for optimization on a sphere by:                          *
!     Roland Lindh                                                     *
!     Uppsala University, Uppsala, Sweden, 2017.                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, n0, n1
real(kind=wp), intent(out) :: E_Pred, C(n)
real(kind=wp), intent(in) :: G(n), H(n,n), D(n,n), r2
integer(kind=iwp) :: i, Iter, j, k, l, nLow
real(kind=wp) :: A, AI, AJ1, AJ2, AJK, AK1, AK2, B, BI, BJ1, BJ2, Eref, Es(4), Step, Step_mi, Step_mj, Step_mk, Step_pi, Step_pj, &
                 Step_pk, Steps(3,4)
#ifdef _DEBUGPRINT_
real(kind=wp) :: rSum
#endif
logical(kind=iwp) :: DidChange
real(kind=wp), external :: Optim_E

!----------------------------------------------------------------------*
! Initialize                                                           *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
! Make first guess.                                                    *
!----------------------------------------------------------------------*
!
! The initial guessed values of the Cs are set by finding the
! coeffcients for the densities with the two lowest energies,
! n0 and n1. All other coefficients are assumed to be zero.
! The coefficient are found by solving the two simultaneous
! equations,
!
! Sum_i C(i) = 0
!
! Sum_i C(i) C(j) D(i,j) = r2

C(:) = Zero
A = Two*(D(n0,n1)-D(n1,n1))/(D(n0,n0)-Two*D(n0,n1)+D(n1,n1))
B = (D(n1,n1)-r2)/(D(n0,n0)-Two*D(n0,n1)+D(n1,n1))
C(n0) = -(A/Two)+sqrt((A/Two)**2-B)
if (C(n0) > One) C(n0) = -(A/Two)-sqrt((A/Two)**2-B)
if (C(n0) < Zero) then
  write(u6,*) 'C(n0) < Zero'
  call Abend()
end if
C(n1) = One-C(n0)
if (C(n1) > One) then
  write(u6,*) 'C(n1) > One'
  call Abend()
end if
if (C(n1) < Zero) then
  write(u6,*) 'C(n1) < Zero'
  call Abend()
end if

! At this point we have a set of Cs which fulfil the two contraints.
#ifdef _DEBUGPRINT_
write(u6,'(a)') 'Start C:'
write(u6,'(6F15.6)') (C(i),i=1,n)
write(u6,'(a)') 'Start G:'
write(u6,'(6F15.6)') (G(i),i=1,n)
call RecPrt('H',' ',H,n,n)
call RecPrt('D',' ',D,n,n)
#endif
!----------------------------------------------------------------------*
! Compute start energy.                                                *
!----------------------------------------------------------------------*
Eref = sum(C(:)*G(:))
do i=1,n
  Eref = Eref+Half*C(i)*sum(C(:)*H(i,:))
end do
#ifdef _DEBUGPRINT_
write(u6,'(a,F15.6)') 'Eref = ',Eref
#endif
!----------------------------------------------------------------------*
! Do a scan in all tripletwise directions.                             *
!----------------------------------------------------------------------*
Iter = 0
DidChange = .true.
Step = 0.1_wp
call Abend()

! Below we explore changes to all triplets that are consistent with
! that the constrains are not broken.

do Iter=1,500

  DidChange = .false.

  do i=1,n-2

    do j=i+1,n-1

      do k=j+1,n

        ! Make sure that step trivially is within the range.
        ! The first constraint. Now we vary d(j) and d(k)
        ! such that the constraints are not broken.

        Step_pi = min(Step,One-C(i))
        Steps(1,1) = Step_pi
        Steps(1,2) = Step_pi

        ! Given the change d(i) compute d(k) and d(l),
        ! subjects to the simultaneous equations.
        !
        ! Sum_i d(i) = 0
        !
        ! Sum_ij (2(C(i)+d(i))*d(j)*D(i,j) = 0

        AI = Two*(C(i)+Step_pi)*Step_pi*D(i,i)+C(k)*Step_pi*D(k,i)+C(j)*Step_pi*D(j,i)
        AJ1 = (C(i)+Two*Step_pi)*D(i,j)+C(j)*D(j,j)+C(k)*D(k,j)
        AK1 = (C(i)+Two*Step_pi)*D(i,k)+C(k)*D(k,k)+C(j)*D(j,k)
        AJK = Two*D(j,k)
        AJ2 = D(j,j)
        AK2 = D(k,k)

        BI = AI-AJ1*Step_pi+AJ2*Step_pi**2
        BJ1 = -AJ1+AK1-(AJK-Two*AJ2)*Step_pi
        BJ2 = -AJK+AJ2+AK2

        A = BJ1/(Two*BJ2)
        B = BI/BJ2

        Step_pk = -A+sqrt(A**2-B)
        Step_mk = -A-sqrt(A**2-B)

        Step_pj = -(Step_pi+Step_pk)
        Step_mj = -(Step_pi+Step_mk)

        C(i) = C(i)+Step_pi

        C(j) = C(j)+Step_mj
        C(k) = C(k)+Step_mk
        Steps(2,1) = Step_mj
        Steps(3,1) = Step_mk

        Es(1) = optim_E(C,G,H,n)

        C(j) = C(j)-Step_mj
        C(k) = C(k)-Step_mk

        C(j) = C(j)+Step_pj
        C(k) = C(k)+Step_pk
        Steps(2,2) = Step_pj
        Steps(3,2) = Step_pk

        Es(2) = optim_E(C,G,H,n)

        C(j) = C(j)-Step_pj
        C(k) = C(k)-Step_pk

        C(i) = C(i)-Step_pi

        Step_mi = min(C(i),Step)
        Steps(1,3) = Step_mi
        Steps(1,4) = Step_mi

        AI = Two*(C(i)+Step_mi)*Step_mi*D(i,i)+C(k)*Step_mi*D(k,i)+C(j)*Step_mi*D(j,i)
        AJ1 = (C(i)+Two*Step_mi)*D(i,j)+C(j)*D(j,j)+C(k)*D(k,j)
        AK1 = (C(i)+Two*Step_mi)*D(i,k)+C(k)*D(k,k)+C(j)*D(j,k)
        AJK = Two*D(j,k)
        AJ2 = D(j,j)
        AK2 = D(k,k)

        BI = AI-AJ1*Step_mi+AJ2*Step_mi**2
        BJ1 = -AJ1+AK1-(AJK-Two*AJ2)*Step_mi
        BJ2 = -AJK+AJ2+AK2

        A = BJ1/(Two*BJ2)
        B = BI/BJ2

        Step_pk = -A+sqrt(A**2-B)
        Step_mk = -A-sqrt(A**2-B)

        Step_pj = -(Step_mi+Step_pk)
        Step_mj = -(Step_mi+Step_mk)

        C(i) = C(i)+Step_mi

        C(j) = C(j)+Step_mj
        C(k) = C(k)+Step_mk
        Steps(2,3) = Step_mj
        Steps(3,3) = Step_mk

        Es(3) = optim_E(C,G,H,n)

        C(j) = C(j)-Step_mj
        C(k) = C(k)-Step_mk

        C(j) = C(j)+Step_pj
        C(k) = C(k)+Step_pk
        Steps(2,4) = Step_pj
        Steps(3,4) = Step_pk

        Es(4) = optim_E(C,G,H,n)

        C(j) = C(j)-Step_pj
        C(k) = C(k)-Step_pk

        C(i) = C(i)-Step_mi

#       ifdef _DEBUGPRINT_
        write(u6,'(a,3I3,4F20.10)') 'i,j,k,Es(i),i=1,4)=',i,j,k,(Es(l),l=1,4)
#       endif
        nLow = 0
        do l=1,4
          if (Es(l) < Eref) then
            nLow = l
            Eref = Es(l)
            DidChange = .true.
          end if
        end do
        if (nLow /= 0) then
          C(i) = C(i)+Steps(1,nLow)
          C(j) = C(j)+Steps(2,nLow)
          C(k) = C(k)+Steps(3,nLow)
        end if

      end do
    end do
  end do

  if (.not. DidChange) then
    if (Step > 0.9e-4_wp) then
      Step = 0.1_wp*Step
      DidChange = .true.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Step is',Step
#     endif
    end if
  end if
# ifdef _DEBUGPRINT_

  ! Check that the constraint is fulfilled.

  rSum = sum(C(:))
  write(u6,*) 'optim: rSum-1',rSum-One

  rSum = Zero
  do i=1,n
    rSum = rSum+C(i)*sum(C(:)*D(i,:))
  end do
  write(u6,*) 'optim: rSum-r2',rSum-r2
# endif
  if (.not. DidChange) exit
end do
#ifdef _DEBUGPRINT_
write(u6,*) 'ERef=',ERef
#endif
E_Pred = ERef

end subroutine Optim2
