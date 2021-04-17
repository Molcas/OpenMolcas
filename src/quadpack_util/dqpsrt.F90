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

subroutine dqpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
!***begin prologue  dqpsrt
!***refer to  dqage,dqagie,dqagpe,dqawse
!***routines called  (none)
!***revision date  810101   (yymmdd)
!***keywords  sequential sorting
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  this routine maintains the descending ordering in the
!            list of the local error estimated resulting from the
!            interval subdivision process. at each call two error
!            estimates are inserted using the sequential search
!            method, top-down for the largest error estimate and
!            bottom-up for the smallest error estimate.
!***description
!
!           ordering routine
!           standard fortran subroutine
!           real*8 version
!
!           parameters (meaning at output)
!              limit  - integer
!                       maximum number of error estimates the list
!                       can contain
!
!              last   - integer
!                       number of error estimates currently in the list
!
!              maxerr - integer
!                       maxerr points to the nrmax-th largest error
!                       estimate currently in the list
!
!              ermax  - real*8
!                       nrmax-th largest error estimate
!                       ermax = elist(maxerr)
!
!              elist  - real*8
!                       vector of dimension last containing
!                       the error estimates
!
!              iord   - integer
!                       vector of dimension last, the first k elements
!                       of which contain pointers to the error
!                       estimates, such that
!                       elist(iord(1)),...,  elist(iord(k))
!                       form a decreasing sequence, with
!                       k = last if last <= (limit/2+2), and
!                       k = limit+1-last otherwise
!
!              nrmax  - integer
!                       maxerr = iord(nrmax)
!
!***end prologue  dqpsrt

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: limit, last
integer(kind=iwp), intent(inout) :: maxerr, iord(last), nrmax
real(kind=wp), intent(in) :: elist(last)
real(kind=wp), intent(out) :: ermax
real(kind=wp) :: errmax, errmin
integer(kind=iwp) :: i, ibeg, ido, isucc, j, jbnd, jupbn, k
logical(kind=iwp) :: jump

!***first executable statement  dqpsrt

! check whether the list contains more than
! two error estimates.

if (last <= 2) then
  iord(1) = 1
  iord(2) = 2
  maxerr = iord(nrmax)
  ermax = elist(maxerr)
  return
end if

! this part of the routine is only executed if, due to a
! difficult integrand, subdivision increased the error
! estimate. in the normal case the insert procedure should
! start after the nrmax-th largest error estimate.

errmax = elist(maxerr)
if (nrmax /= 1) then
  ido = nrmax-1
  do i=1,ido
    isucc = iord(nrmax-1)
    ! ***jump out of do-loop
    if (errmax <= elist(isucc)) exit
    iord(nrmax) = isucc
    nrmax = nrmax-1
  end do
end if

! compute the number of elements in the list to be maintained
! in descending order. this number depends on the number of
! subdivisions still allowed.

jupbn = last
if (last > (limit/2+2)) jupbn = limit+3-last
errmin = elist(last)

! insert errmax by traversing the list top-down,
! starting comparison from the element elist(iord(nrmax+1)).

jbnd = jupbn-1
ibeg = nrmax+1
jump = .false.
if (ibeg <= jbnd) then
  do i=ibeg,jbnd
    isucc = iord(i)
    ! ***jump out of do-loop
    if (errmax >= elist(isucc)) then
      jump = .true.
      exit
    end if
    iord(i-1) = isucc
  end do
end if
if (.not. jump) then
  iord(jbnd) = maxerr
  iord(jupbn) = last
else

  ! insert errmin by traversing the list bottom-up.

  iord(i-1) = maxerr
  k = jbnd
  do j=i,jbnd
    isucc = iord(k)
    ! ***jump out of do-loop
    if (errmin < elist(isucc)) then
      iord(k+1) = last
      exit
    end if
    iord(k+1) = isucc
    k = k-1
  end do
  iord(i) = last
end if

! set maxerr and ermax.

maxerr = iord(nrmax)
ermax = elist(maxerr)

return

end subroutine dqpsrt
