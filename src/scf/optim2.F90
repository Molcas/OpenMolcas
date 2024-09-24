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
subroutine Optim2(E_Pred,G,H,D,C,n,nDim,n0,n1,r2)
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
!       nDim    : leading dimension for H                              *
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

use Constants, only: Zero, Half, One, Two

implicit none
!----------------------------------------------------------------------*
! Dummy arguments.                                                     *
!----------------------------------------------------------------------*
integer n, nDim, n0, n1
real*8 G(nDim), H(nDim,nDim), C(nDim), D(nDim,nDim)
!----------------------------------------------------------------------*
! Local variables.                                                     *
!----------------------------------------------------------------------*
real*8 Step, Eref, E_Pred
real*8 Step_pi, Step_mi, Step_pj, Step_mj, Step_pk, Step_mk
real*8 Optim_E
#ifdef _DEBUGPRINT_
real*8 sum
#endif
logical DidChange
integer Iter
integer i, j, k, l, nLow
real*8 A, B, r2
real*8 AI, AJ1, AK1, AJK, AJ2, AK2
real*8 BI, BJ1, BJ2
real*8 Steps(3,4), Es(4)

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

do i=1,n
  C(i) = Zero
end do
A = Two*(D(n0,n1)-D(n1,n1))/(D(n0,n0)-Two*D(n0,n1)+D(n1,n1))
B = (D(n1,n1)-r2)/(D(n0,n0)-Two*D(n0,n1)+D(n1,n1))
C(n0) = -(A/Two)+sqrt((A/Two)**2-B)
if (C(n0) > One) C(n0) = -(A/Two)-sqrt((A/Two)**2-B)
if (C(n0) < Zero) then
  write(6,*) 'C(n0) < Zero'
  call Abend()
end if
C(n1) = One-C(n0)
if (C(n1) > One) then
  write(6,*) 'C(n1) > One'
  call Abend()
end if
if (C(n1) < Zero) then
  write(6,*) 'C(n1) < Zero'
  call Abend()
end if

! At this point we have a set of Cs which fulfil the two contraints.
#ifdef _DEBUGPRINT_
write(6,'(a)') 'Start C:'
write(6,'(6F15.6)') (C(i),i=1,n)
write(6,'(a)') 'Start G:'
write(6,'(6F15.6)') (G(i),i=1,n)
call RecPrt('H',' ',H,nDim,nDim)
call RecPrt('D',' ',D,nDim,nDim)
#endif
!----------------------------------------------------------------------*
! Compute start energy.                                                *
!----------------------------------------------------------------------*
Eref = Zero
do i=1,n
  Eref = Eref+C(i)*G(i)
  do j=1,n
    Eref = Eref+Half*C(i)*C(j)*H(i,j)
  end do
end do
#ifdef _DEBUGPRINT_
write(6,'(a,F15.6)') 'Eref = ',Eref
#endif
!----------------------------------------------------------------------*
! Do a scan in all tripletwise directions.                             *
!----------------------------------------------------------------------*
Iter = 0
DidChange = .true.
Step = 0.10d0
call Abend()

! Below we explore changes to all triplets that are consistent with
! that the constrains are not broken.

100 continue
if ((Iter < 500) .and. DidChange) then

  Iter = Iter+1
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

        Es(1) = optim_E(C,G,H,n,nDim)

        C(j) = C(j)-Step_mj
        C(k) = C(k)-Step_mk

        C(j) = C(j)+Step_pj
        C(k) = C(k)+Step_pk
        Steps(2,2) = Step_pj
        Steps(3,2) = Step_pk

        Es(2) = optim_E(C,G,H,n,nDim)

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

        Es(3) = optim_E(C,G,H,n,nDim)

        C(j) = C(j)-Step_mj
        C(k) = C(k)-Step_mk

        C(j) = C(j)+Step_pj
        C(k) = C(k)+Step_pk
        Steps(2,4) = Step_pj
        Steps(3,4) = Step_pk

        Es(4) = optim_E(C,G,H,n,nDim)

        C(j) = C(j)-Step_pj
        C(k) = C(k)-Step_pk

        C(i) = C(i)-Step_mi

#       ifdef _DEBUGPRINT_
        write(6,'(a,3I3,4F20.10)') 'i,j,k,Es(i),i=1,4)=',i,j,k,(Es(l),l=1,4)
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
    if (Step > 0.9d-4) then
      Step = 0.1d0*Step
      DidChange = .true.
#     ifdef _DEBUGPRINT_
      write(6,*) 'Step is',Step
#     endif
    end if
  end if
# ifdef _DEBUGPRINT_

  ! Check that the constraint is fullfilled.

  sum = Zero
  do i=1,n
    sum = sum+C(i)
  end do
  write(6,*) 'optim: sum-1',sum-One

  sum = Zero
  do i=1,n
    do j=1,n
      sum = sum+C(i)*C(j)*D(i,j)
    end do
  end do
  write(6,*) 'optim: sum-r2',sum-r2
# endif
  Go To 100
end if
#ifdef _DEBUGPRINT_
write(6,*) 'ERef=',ERef
#endif
E_Pred = ERef

end subroutine Optim2
