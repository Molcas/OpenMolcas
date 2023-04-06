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

subroutine multdot(wrk,wrksize,nind,a,ssa,b,ssb,scalar,rc)
! this routine does dot product
! scalar = A(ind)*B(ind)
!
! nind   - number of indices in both A B (I)
! a      - map type of A (I)
! ssa    - symmetry state of A (I)
! b      - map type of A (I)
! ssb    - symmetry state of B (I)
! scalar - final dot product (O)
! rc     - return (error) code (O)
!
! N.B. A and B must have same ordering of indices
!
! indA     indB     Implemented
! 1        1         Not yet
! 2        2           Yes
! 3        3           Yes
! 4        4           Yes
! more                 No

use ccsd_global, only: Map_Type
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, nind, ssa, ssb
real(kind=wp), intent(in) :: wrk(wrksize)
type(Map_Type), intent(in) :: a, b
real(kind=wp), intent(out) :: scalar
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: iia, iib, length, nhelp, posa, posb, symp, symq, symr
real(kind=wp) :: scal

rc = 0

!T some tests

do nhelp=1,nind
  if (a%d(0,nhelp) /= b%d(0,nhelp)) then
    ! RC =1 : nonidentical types of indices (NCI/Stup)
    rc = 1
    return
  end if
end do

if (a%d(0,5) /= b%d(0,5)) then
  ! RC =2 : nonidentical number of symmetry blocks in A and B (Stup)
  rc = 2
  return
end if

if (a%d(0,6) /= b%d(0,6)) then
  ! RC =3 : nonidentical type of A and B (Stup)
  rc = 3
  return
end if

if (ssa /= ssb) then
  ! RC =4 : nonidentical symmetry state of A and B (Stup)
  rc = 4
  return
end if

if (nind == 4) then

  !I 4 index matrices

  scalar = Zero
  do iia=1,a%d(0,5)

    !I.1 def parameters of A
    symp = a%d(iia,3)
    symq = a%d(iia,4)
    symr = a%d(iia,5)
    ! syms is redundant
    posa = a%d(iia,1)

    !I.2 def parameters of B
    iib = b%i(symp,symq,symr)
    posb = b%d(iib,1)

    !I.3 length must be common for both A and B
    length = a%d(iia,2)

    if (length > 0) then
      call mr0u3wt(length,length,length,1,1,wrk(posa),wrk(posb),scal)
      scalar = scalar+scal
    end if

  end do

else if (nind == 3) then

  !II 3 index matrices

  scalar = Zero
  do iia=1,a%d(0,5)

    !II.1 def parameters of A
    symp = a%d(iia,3)
    symq = a%d(iia,4)
    ! symr is redundant
    posa = a%d(iia,1)

    !II.2 def parameters of B
    iib = b%i(symp,symq,1)
    posb = b%d(iib,1)

    !II.3 length must be common for both A and B
    length = a%d(iia,2)

    if (length > 0) then
      call mr0u3wt(length,length,length,1,1,wrk(posa),wrk(posb),scal)
      scalar = scalar+scal
    end if

  end do

else if (nind == 2) then

  !III 2 index matrices

  scalar = Zero
  do iia=1,a%d(0,5)

    !III.1 def parameters of A
    symp = a%d(iia,3)
    ! symq is redundant
    posa = a%d(iia,1)

    !III.2 def parameters of B
    iib = b%i(symp,1,1)
    posb = b%d(iib,1)

    !III.3 length must be common for both A and B
    length = a%d(iia,2)

    if (length > 0) then
      call mr0u3wt(length,length,length,1,1,wrk(posa),wrk(posb),scal)
      scalar = scalar+scal
    end if

  end do

else if (nind == 1) then

  !IV 1 index matrices
  ! RC=5 : 1 index matrices (NCI)
  rc = 5
  return

else
  !V more than 4 index matrices
  ! RC=6 : more than 4 index matrices (NCI)
  rc = 6
  return
end if

return

end subroutine multdot
