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

subroutine LevMarquart(iPotte,nPick,ipPick,nEPP,ipEPCo,Coo,dMullig,lMax,A,iAtom,jAtom,chP,Thrs1,Thrs2,nThrs,Chi2B,iPrint,AboveMul)

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
#include "warnings.h"
dimension Coo(3), A(2), B(2), dMullig((lMax*(lMax**2+6*lMax+11)+6)/6)
dimension AlfMat(4), AlfMatI(4), rStore(nPick), rinvStore(nPick)
dimension ARaw(2,2), BRaw(2), dA(2), Pout(nPick)
dimension xStore(nPick), yStore(nPick), zStore(nPick)
logical AboveMul(2)
logical lStop1, lStop2, lStop3, lStop4
logical lScreen1, lScreen2, lScreen3, lScreen4
character*60 UtChar

! Set iteration count to zero here at the top of the
! Levenberg-Marquart thing. Initiate lambda and the exponents.
! The exponents are not God-given, but these rather low values
! have been found to work nicely in most cases. Set highest
! and lowest allowed values of exponents. Set highest and lowest
! allowed values of exponent steps.

Iter = 0
nStep = 0
dLambda = 1d-3
A(1) = 0.5d0
A(2) = 0.5d0
dUpper = 3.0d0
dLower = 1d-1
ddUpper = 0.7d0
ddLower = -0.7d0

! Start loop, set zeros.

9901 continue
Iter = Iter+1
Chi2 = 0.0d0
Chi2B = 0.0d0
do i=1,2
  BRaw(i) = 0.0d0
  dA(i) = 0.0d0
  do j=1,2
    ARaw(i,j) = 0.0d0
  end do
end do

! Compute distances and store them.

do iP=1,nPick
  ind = iWork(ipPick+iP-1)
  x = Work(ipEPCo+(ind-1)*3+0)-Coo(1)
  y = Work(ipEPCo+(ind-1)*3+1)-Coo(2)
  z = Work(ipEPCo+(ind-1)*3+2)-Coo(3)
  xStore(iP) = x
  yStore(iP) = y
  zStore(iP) = z
  rStore(iP) = sqrt(x**2+y**2+z**2)
  rinvStore(iP) = 1.0d0/rStore(iP)
end do

! Compute the electric potential and the difference between true and
! modelled potential. Also compute first derivatives. Formaly, there
! are second derivatives as well, but they are neglected for
! stability reasons. For this relatively simple non-linear problem,
! it is probably not meaningful to code a more intelligent
! routine that includes the second derivatives with some
! stabalizing modifications.

DerMax = 0.0d0
do iP=1,nPick
  r = rStore(iP)
  rinv = rinvStore(iP)
  x = xStore(iP)
  y = yStore(iP)
  z = zStore(iP)
  Potte = ElPot(r,rinv,x,y,z,dMullig,lMax,A,chP,.true.,.true.)
  Diffo = Work(iPotte+iP-1)-Potte
  Der1 = dMullig(1)*(1.0d0+2.0d0*A(1)*r)*exp(-2.0d0*A(1)*r)
  Der2 = (dMullig(2)*x+dMullig(3)*y+dMullig(4)*z)*(A(2)**2+2.0d0*A(2)**3*r)*exp(-2.0d0*A(2)*r)
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
Chi2 = Chi2/dble(nPick)

! Check if numerical and statistical well-behaved calculations are
! expected or not.

if (nPick < 5) then
  write(6,*)
  write(6,'(A)') ' Error in numerical fit for exponent!'
  write(6,'(A)') '    Too few electric potential points have'
  write(6,'(A,2I4)') '    been sampled for (At1,At2):',iAtom,jAtom
  write(6,*)
  call Quit(_RC_GENERAL_ERROR_)
end if
if (DerMax < 1d-4) then
  write(6,*)
  write(6,'(A)') ' Error in numerical fit for exponent!'
  write(6,'(A)') '    To small derivative for parameter'
  write(6,'(A,2I4)') '    in centre (At1,At2):',iAtom,jAtom
  write(6,'(A,E12.4)') '    Maximal derivative:',DerMax
  write(6,*)
  write(6,'(A)') '    Either include closer points, give a'
  write(6,'(A)') '    better initial estimate, or increase'
  write(6,'(A)') '    multipole magnitude threshold.'
  write(6,*)
  call Quit(_RC_GENERAL_ERROR_)
end if

! Solve the linear system to get step, dA.

call SolveA(AlfMat,AlfMatI,dLambda,dMullig,lMax,ARaw,BRaw,dA,iPrint,AboveMul,ddUpper,ddLower)

! Make a screening of dA: if the maximal or minimal exponent
! has been reached, but the optimization still pushes on, screen
! them by putting them as zero, and modifying abovemul.

do i=1,2
  if (AboveMul(i)) then
    lScreen1 = (A(i)+1d-8) > dUpper
    lScreen2 = (dA(i)+1d-8) > ddUpper
    lScreen3 = (A(i)-1d-8) < dLower
    lScreen4 = (dA(i)-1d-8) < ddLower
    if ((lScreen1 .and. lScreen2) .or. (lScreen3 .and. lScreen4)) then
      dA(i) = 0.0d0
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
  Potte = ElPot(r,rinv,x,y,z,dMullig,lMax,B,chP,.true.,.true.)
  Diffo = Work(iPotte+iP-1)-Potte
  Chi2B = Chi2B+Diffo**2
end do
Chi2B = Chi2B/dble(nPick)

! Take appropriate action given the different Stop-criterion. They
! measure the following: (1) The error should not change too much
! between steps; observe that it is not meaningful to have a too
! tight threshold in this regard due to the statistical nature of
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
  write(6,*)
  write(6,790) Iter
  write(6,794) 'Chi2','Chi2, trial','Lambda'
  write(6,791) Chi2,Chi2B,dLambda
  write(6,795) 'Exponents:','Charge','Dipole'
  write(6,792) A(1),A(2)
end if
if (lStop1 .and. lStop2 .and. lStop3 .and. lStop4) then
  A(1) = B(1)
  A(2) = B(2)
  Go To 9902
end if
if (Iter > 40) then
  write(6,*)
  write(6,*) 'No convergence in the Levenberg-Marquart section.'
  call Quit(_RC_NOT_CONVERGED_)
end if
if (Chi2B >= Chi2) then
  dLambda = dLambda*10.0d0
  Go To 9901
else
  dLambda = dLambda*0.1d0
  nStep = nStep+1
  A(1) = B(1)
  A(2) = B(2)
  Go To 9901
end if

9902 continue

! Optional printing when convergence is reached.

if (iPrint >= 5) then
  write(6,*)
  write(6,*)
  write(6,793) Iter
  write(6,794) 'Chi2','Chi2, trial','Lambda'
  write(6,791) Chi2,Chi2B,dLambda
  write(6,795) 'Exponents:','Charge','Dipole'
  write(6,792) A(1),A(2)
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
end if

790 format('Levenberg-Marquart optimization. Iteration:',I3)
794 format('   ',A,'          ',A,'   ',A)
791 format(3(E13.6,' '))
795 format(' ',A,' ',A,'      ',A)
792 format('       ',2F12.6)
793 format('Convergence reached in iteration ',I2)

!nd yes, there is an end here as well.

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nEPP)

end subroutine LevMarquart
