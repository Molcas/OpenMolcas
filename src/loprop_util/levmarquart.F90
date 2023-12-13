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

subroutine LevMarquart(Potte,nPick,Pick,EPCo,Coo,dMullig,lMax,A,iAtom,jAtom,chP,Thrs1,Thrs2,nThrs,Chi2B,iPrint,AboveMul)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Ten, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nPick, Pick(nPick), lMax, iAtom, jAtom, nThrs, iPrint
real(kind=wp), intent(in) :: Potte(nPick), EPCo(3,*), Coo(3), dMullig((lMax+1)*(lMax+2)*(lMax+3)/6), chP, Thrs1, Thrs2
real(kind=wp), intent(out) :: A(2), Chi2B
logical(kind=iwp), intent(inout) :: AboveMul(2)
integer(kind=iwp) :: i, ind, iP, Iter, j, nStep
real(kind=wp) :: AlfMat(4), AlfMatI(4), ARaw(2,2), B(2), BRaw(2), Chi2, dA(2), ddLower, ddUpper, Der1, Der2, DerMax, Diffo, &
                 dLambda, dLower, dUpper, r, rinv, x, y, z
logical(kind=iwp) :: lScreen1, lScreen2, lScreen3, lScreen4, lStop1, lStop2, lStop3, lStop4
character(len=60) :: UtChar
real(kind=wp), allocatable :: Pout(:), rinvStore(:), rStore(:), xStore(:), yStore(:), zStore(:)
real(kind=wp), external :: ElPot
#include "warnings.h"

! Set iteration count to zero here at the top of the
! Levenberg-Marquart thing. Initiate lambda and the exponents.
! The exponents are not God-given, but these rather low values
! have been found to work nicely in most cases. Set highest
! and lowest allowed values of exponents. Set highest and lowest
! allowed values of exponent steps.

Iter = 0
nStep = 0
dLambda = 1.0e-3_wp
A(:) = Half
dUpper = Three
dLower = 0.1_wp
ddUpper = 0.7_wp
ddLower = -0.7_wp

! Start loop, set zeros.

do
  Iter = Iter+1
  Chi2 = Zero
  Chi2B = Zero
  do i=1,2
    BRaw(i) = Zero
    dA(i) = Zero
    do j=1,2
      ARaw(i,j) = Zero
    end do
  end do

  ! Compute distances and store them.

  call mma_allocate(xStore,nPick,label='xStore')
  call mma_allocate(yStore,nPick,label='yStore')
  call mma_allocate(zStore,nPick,label='zStore')
  call mma_allocate(rStore,nPick,label='rStore')
  call mma_allocate(rinvStore,nPick,label='rinvStore')
  do iP=1,nPick
    ind = Pick(iP)
    x = EPCo(1,ind)-Coo(1)
    y = EPCo(2,ind)-Coo(2)
    z = EPCo(3,ind)-Coo(3)
    xStore(iP) = x
    yStore(iP) = y
    zStore(iP) = z
    rStore(iP) = sqrt(x**2+y**2+z**2)
    rinvStore(iP) = One/rStore(iP)
  end do

  ! Compute the electric potential and the difference between true and
  ! modelled potential. Also compute first derivatives. Formaly, there
  ! are second derivatives as well, but they are neglected for
  ! stability reasons. For this relatively simple non-linear problem,
  ! it is probably not meaningful to code a more intelligent
  ! routine that includes the second derivatives with some
  ! stabalizing modifications.

  DerMax = Zero
  do iP=1,nPick
    r = rStore(iP)
    rinv = rinvStore(iP)
    x = xStore(iP)
    y = yStore(iP)
    z = zStore(iP)
    Diffo = Potte(iP)-ElPot(r,rinv,x,y,z,dMullig,lMax,A,chP,.true.,.true.)
    Der1 = dMullig(1)*(One+Two*A(1)*r)*exp(-Two*A(1)*r)
    Der2 = (dMullig(2)*x+dMullig(3)*y+dMullig(4)*z)*(A(2)**2+Two*A(2)**3*r)*exp(-Two*A(2)*r)
    if (abs(Der1) > DerMax) Dermax = abs(Der1)
    if (abs(Der2) > DerMax) Dermax = abs(Der2)

    ! Assemble the alpha-matrix, the beta-vector and the chi2 scalar.

    Chi2 = Chi2+Diffo**2
    BRaw(1) = BRaw(1)+Diffo*Der1
    BRaw(2) = BRaw(2)+Diffo*Der2
    ARaw(1,1) = ARaw(1,1)+Der1*Der1
    ARaw(2,1) = ARaw(2,1)+Der2*Der1
    ARaw(2,2) = ARaw(2,2)+Der2*Der2
  end do
  Chi2 = Chi2/real(nPick,kind=wp)

  ! Check if numerical and statistical well-behaved calculations are
  ! expected or not.

  if (nPick < 5) then
    write(u6,*)
    write(u6,'(A)') ' Error in numerical fit for exponent!'
    write(u6,'(A)') '    Too few electric potential points have'
    write(u6,'(A,2I4)') '    been sampled for (At1,At2):',iAtom,jAtom
    write(u6,*)
    call Quit(_RC_GENERAL_ERROR_)
  end if
  if (DerMax < 1.0e-4_wp) then
    write(u6,*)
    write(u6,'(A)') ' Error in numerical fit for exponent!'
    write(u6,'(A)') '    To small derivative for parameter'
    write(u6,'(A,2I4)') '    in centre (At1,At2):',iAtom,jAtom
    write(u6,'(A,ES12.4)') '    Maximal derivative:',DerMax
    write(u6,*)
    write(u6,'(A)') '    Either include closer points, give a'
    write(u6,'(A)') '    better initial estimate, or increase'
    write(u6,'(A)') '    multipole magnitude threshold.'
    write(u6,*)
    call Quit(_RC_GENERAL_ERROR_)
  end if

  ! Solve the linear system to get step, dA.

  call SolveA(AlfMat,AlfMatI,dLambda,ARaw,BRaw,dA,iPrint,AboveMul,ddUpper,ddLower)

  ! Make a screening of dA: if the maximal or minimal exponent
  ! has been reached, but the optimization still pushes on, screen
  ! them by putting them as zero, and modifying AboveMul.

  do i=1,2
    if (AboveMul(i)) then
      lScreen1 = (A(i)+1.0e-8_wp) > dUpper
      lScreen2 = (dA(i)+1.0e-8_wp) > ddUpper
      lScreen3 = (A(i)-1.0e-8_wp) < dLower
      lScreen4 = (dA(i)-1.0e-8_wp) < ddLower
      if ((lScreen1 .and. lScreen2) .or. (lScreen3 .and. lScreen4)) then
        dA(i) = Zero
        AboveMul(i) = .false.
      end if
    end if
  end do
  if (iPrint >= 7) then
    call RecPrt('deltaA',' ',dA,2,1)
  end if

  ! Construct trial parameters and compute the chi2-scalar for
  ! the trial parameters. To keep things within reasonable bounds
  ! do not exceed certain values.

  B(1) = A(1)+dA(1)
  B(2) = A(2)+dA(2)
  if (B(1) < dLower) B(1) = dLower
  if (B(2) < dLower) B(2) = dLower
  if (B(1) > dUpper) B(1) = dUpper
  if (B(2) > dUpper) B(2) = dUpper
  do iP=1,nPick
    r = rStore(iP)
    rinv = rinvStore(iP)
    x = xStore(iP)
    y = yStore(iP)
    z = zStore(iP)
    Diffo = Potte(iP)-ElPot(r,rinv,x,y,z,dMullig,lMax,B,chP,.true.,.true.)
    Chi2B = Chi2B+Diffo**2
  end do
  Chi2B = Chi2B/real(nPick,kind=wp)

  ! Take appropriate action given the different Stop-criterion. They
  ! measure the following: (1) The error should not change too much
  ! between steps; observe that it is not meaningful to have too
  ! tight a threshold in this regard due to the statistical nature of
  ! the problem. (2) The modifier parameter has to stabalize and not
  ! be far out in the linear region, rather in the second order regime.
  ! (3) The last step should be a decrease. (4) A certain number of
  ! steps should preceed that decreases the error; this threshold
  ! 'overlaps' some with second threshold. Halt the optimization if
  ! no convergence is reached after 40 iterations. This is by far
  ! a generous limit.

  lStop1 = abs(Chi2-Chi2B) < Thrs1
  lStop2 = dLambda < Thrs2
  lStop3 = Chi2 > Chi2B
  lStop4 = nStep >= nThrs
  if (iPrint >= 6) then
    write(u6,*)
    write(u6,790) Iter
    write(u6,794) 'Chi2','Chi2, trial','Lambda'
    write(u6,791) Chi2,Chi2B,dLambda
    write(u6,795) 'Exponents:','Charge','Dipole'
    write(u6,792) A(1),A(2)
  end if
  if (lStop1 .and. lStop2 .and. lStop3 .and. lStop4) then
    A(1) = B(1)
    A(2) = B(2)
    exit
  end if
  if (Iter > 40) then
    write(u6,*)
    write(u6,*) 'No convergence in the Levenberg-Marquart section.'
    call Quit(_RC_NOT_CONVERGED_)
  end if
  if (Chi2B >= Chi2) then
    dLambda = dLambda*Ten
  else
    dLambda = dLambda/Ten
    nStep = nStep+1
    A(1) = B(1)
    A(2) = B(2)
  end if
end do

! Optional printing when convergence is reached.

if (iPrint >= 5) then
  write(u6,*)
  write(u6,*)
  write(u6,793) Iter
  write(u6,794) 'Chi2','Chi2, trial','Lambda'
  write(u6,791) Chi2,Chi2B,dLambda
  write(u6,795) 'Exponents:','Charge','Dipole'
  write(u6,792) A(1),A(2)
  call mma_allocate(Pout,nPick,label='Pout')
  do iP=1,nPick
    r = rStore(iP)
    rinv = rinvStore(iP)
    x = xStore(iP)
    y = yStore(iP)
    z = zStore(iP)
    Pout(iP) = ElPot(r,rinv,x,y,z,dMullig,lMax,A,chP,.true.,.true.)
  end do
  write(UtChar,'(A,2I3)') 'Approximate partial density potential, centre',iAtom,jAtom
  call RecPrt(UtChar,' ',Pout,nPick,1)
  call RecPrt('Distance to points',' ',rStore,nPick,1)
  call mma_deallocate(Pout)
end if

call mma_deallocate(xStore)
call mma_deallocate(yStore)
call mma_deallocate(zStore)
call mma_deallocate(rStore)
call mma_deallocate(rinvStore)

! And yes, there is an end here as well.

return

790 format('Levenberg-Marquart optimization. Iteration:',I3)
794 format('   ',A,'          ',A,'   ',A)
791 format(3(ES13.6,' '))
795 format(' ',A,' ',A,'      ',A)
792 format('       ',2F12.6)
793 format('Convergence reached in iteration ',I2)

end subroutine LevMarquart
