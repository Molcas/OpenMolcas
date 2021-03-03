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

subroutine dqelg(n,epstab,reslt,abserr,res3la,nres)
!***begin prologue  dqelg
!***refer to  dqagie,dqagoe,dqagpe,dqagse
!***routines called  d1mach
!***revision date  830518   (yymmdd)
!***keywords  epsilon algorithm, convergence acceleration,
!             extrapolation
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine determines the limit of a given sequence of
!            approximations, by means of the epsilon algorithm of
!            p.wynn. an estimate of the absolute error is also given.
!            the condensed epsilon table is computed. only those
!            elements needed for the computation of the next diagonal
!            are preserved.
!***description
!
!           epsilon algorithm
!           standard fortran subroutine
!           real*8 version
!
!           parameters
!              n      - integer
!                       epstab(n) contains the new element in the
!                       first column of the epsilon table.
!
!              epstab - real*8
!                       vector of dimension 52 containing the elements
!                       of the two lower diagonals of the triangular
!                       epsilon table. the elements are numbered
!                       starting at the right-hand corner of the
!                       triangle.
!
!              reslt  - real*8
!                       resulting approximation to the integral
!
!              abserr - real*8
!                       estimate of the absolute error computed from
!                       result and the 3 previous results
!
!              res3la - real*8
!                       vector of dimension 3 containing the last 3
!                       results
!
!              nres   - integer
!                       number of calls to the routine
!                       (should be zero at first call)
!
!***end prologue  dqelg

use Constants, only: One, Five
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: n, nres
real(kind=wp), intent(inout) :: epstab(52), res3la(3)
real(kind=wp), intent(out) :: reslt, abserr
integer(kind=iwp) :: i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num
real(kind=wp) :: delta1, delta2, delta3, epmach, epsinf, error, err1, err2, err3, e0, e1, e1abs, e2, e3, oflow, res, ss, tol1, &
                 tol2, tol3
real(kind=wp), external :: d1mach

! list of major variables
! -----------------------
!
! e0     - the 4 elements on which the computation of a new
! e1       element in the epsilon table is based
! e2
! e3                 e0
!              e3    e1    new
!                    e2
! newelm - number of elements to be computed in the new
!          diagonal
! error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
! reslt  - the element in the new diagonal with least value
!          of error
!
! machine dependent constants
! ---------------------------
!
! epmach is the largest relative spacing.
! oflow is the largest positive magnitude.
! limexp is the maximum number of elements the epsilon
! table can contain. if this number is reached, the upper
! diagonal of the epsilon table is deleted.
!
!***first executable statement  dqelg

epmach = d1mach(4)
oflow = d1mach(2)
nres = nres+1
abserr = oflow
reslt = epstab(n)
if (n < 3) then
  abserr = max(abserr,Five*epmach*abs(reslt))
  return
end if
limexp = 50
epstab(n+2) = epstab(n)
newelm = (n-1)/2
epstab(n) = oflow
num = n
k1 = n
do i=1,newelm
  k2 = k1-1
  k3 = k1-2
  res = epstab(k1+2)
  e0 = epstab(k3)
  e1 = epstab(k2)
  e2 = res
  e1abs = abs(e1)
  delta2 = e2-e1
  err2 = abs(delta2)
  tol2 = max(abs(e2),e1abs)*epmach
  delta3 = e1-e0
  err3 = abs(delta3)
  tol3 = max(e1abs,abs(e0))*epmach
  if ((err2 < tol2) .and. (err3 < tol3)) then

    ! if e0, e1 and e2 are equal to within machine
    ! accuracy, convergence is assumed.
    ! reslt = e2
    ! abserr = abs(e1-e0)+abs(e2-e1)

    reslt = res
    abserr = err2+err3
    ! ***jump out of do-loop
    abserr = max(abserr,Five*epmach*abs(reslt))
    return
  end if
  e3 = epstab(k1)
  epstab(k1) = e1
  delta1 = e1-e3
  err1 = abs(delta1)
  tol1 = max(e1abs,abs(e3))*epmach

  ! if two elements are very close to each other, omit
  ! a part of the table by adjusting the value of n

  if ((err1 <= tol1) .or. (err2 <= tol2) .or. (err3 <= tol3)) then
    n = i+i-1
    ! ***jump out of do-loop
    exit
  end if
  ss = One/delta1+One/delta2-One/delta3
  epsinf = abs(ss*e1)

  ! test to detect irregular behaviour in the table, and
  ! eventually omit a part of the table adjusting the value
  ! of n.

  if (epsinf <= 1.0e-4_wp) then
    n = i+i-1
    ! ***jump out of do-loop
    exit
  end if

  ! compute a new element and eventually adjust
  ! the value of reslt.

  res = e1+One/ss
  epstab(k1) = res
  k1 = k1-2
  error = err2+abs(res-e2)+err3
  if (error <= abserr) cycle
  abserr = error
  reslt = res
end do

! shift the table.

if (n == limexp) n = 2*(limexp/2)-1
ib = 1
if ((num/2)*2 == num) ib = 2
ie = newelm+1
do i=1,ie
  ib2 = ib+2
  epstab(ib) = epstab(ib2)
  ib = ib2
end do
if (num /= n) then
  indx = num-n+1
  do i=1,n
    epstab(i) = epstab(indx)
    indx = indx+1
  end do
end if
if (nres < 4) then
  res3la(nres) = reslt
  abserr = oflow
  abserr = max(abserr,Five*epmach*abs(reslt))
  return
end if

! compute error estimate

abserr = abs(reslt-res3la(3))+abs(reslt-res3la(2))+abs(reslt-res3la(1))
res3la(1) = res3la(2)
res3la(2) = res3la(3)
res3la(3) = reslt
abserr = max(abserr,Five*epmach*abs(reslt))

return

end subroutine dqelg
