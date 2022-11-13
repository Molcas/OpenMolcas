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

subroutine SolveA(AlfMat,AlfMatI,dLambda,ARaw,BRaw,dA,iPrint,AboveMul,ddUpper,ddLower)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: AlfMat(4), AlfMatI(4), dA(2)
integer(kind=iwp), intent(in) :: iPrint
real(kind=wp), intent(in) :: dLambda, ARaw(2,2), BRaw(2), ddUpper, ddLower
logical(kind=iwp), intent(in) :: AboveMul(2)
integer(kind=iwp) :: i, iMull, j, kaunt, nDim
real(kind=wp) :: Beta(2), Det, dtA(2)
logical(kind=iwp) :: Yeps(2)

! Which elements can be non-zero? Too low magnitude of multipoles
! and they are screened.

nDim = 0
do iMull=1,2
  if (AboveMul(iMull)) then
    Yeps(iMull) = .true.
    nDim = nDim+1
    Beta(nDim) = BRaw(iMull)
  else
    Yeps(iMull) = .false.
  end if
end do
if (iPrint >= 10) then
  call RecPrt('Beta',' ',Beta,nDim,1)
end if

! Shuffle around, and create a matrix for exponents with non-zero factors.

kaunt = 0
do i=1,2
  do j=1,2
    if (Yeps(i) .and. Yeps(j)) then
      kaunt = kaunt+1
      if (i == j) then
        AlfMat(kaunt) = ARaw(max(i,j),min(i,j))*(One+dLambda)
      else
        AlfMat(kaunt) = ARaw(max(i,j),min(i,j))
      end if
    end if
  end do
end do

! Invert and solve.

call MInv(AlfMat,AlfMatI,Det,nDim)
call dGeMV_('N',nDim,nDim,One,AlfMatI,nDim,Beta,1,Zero,dtA,1)

! Optional printing.

if (iPrint >= 10) then
  call RecPrt('Alfa',' ',AlfMat,nDim,nDim)
  call RecPrt('InverseA',' ',AlfMatI,nDim,nDim)
  call RecPrt('deltatA',' ',dtA,nDim,1)
end if

! Damp large steps since such steps can take the optimization
! too far away and put it in a region with very small derivatives,
! and there things turns into baloney.

if (dtA(1) < ddLower) dtA(1) = ddLower
if (dtA(2) < ddLower) dtA(2) = ddLower
if (dtA(1) > ddUpper) dtA(1) = ddUpper
if (dtA(2) > ddUpper) dtA(2) = ddUpper

! Extend to full dimension.

kaunt = 0
do i=1,2
  if (Yeps(i)) then
    kaunt = kaunt+1
    dA(i) = dtA(kaunt)
  else
    dA(i) = Zero
  end if
end do

return

end subroutine SolveA
