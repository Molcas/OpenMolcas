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

function D_Bend(Ind,iOp_,nSym)

use Slapaf_Info, only: jStab, nStab
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: D_Bend
integer(kind=iwp), intent(in) :: Ind(3), iOp_(3), nSym
integer(kind=iwp) :: iAtom, ideg, iOp_E, iOp_ER, iOp_ET, iOp_R, iOp_T, iU_A, iU_AB, iU_ABC, iU_AC, iU_B, iU_C, jAtom, mAtom, nU_A, &
                     nU_ABC, nU_B, nU_C
integer(kind=iwp), external :: iU, iUR, nU

!                                                                      *
!***********************************************************************
!                                                                      *
! B E N D  (iAtom,mAtom,jAtom)

D_Bend = Zero

iAtom = Ind(1)
mAtom = Ind(2)
jAtom = Ind(3)
iOp_E = iOp_(1)
iOp_R = iOp_(2)
iOp_T = iOp_(3)

! Pick up the dimension of the stabilizer and store the operator
! indices in iU_*

nU_A = nStab(iAtom)
iU_A = iU(jStab(0,iAtom),nU_A)

nU_B = nStab(mAtom)
iU_B = iU(jStab(0,mAtom),nU_B)

nU_C = nStab(jAtom)
iU_C = iU(jStab(0,jAtom),nU_C)

!write(u6,*) ' U_A',iU_A,nU_A
!write(u6,*) ' U_B',iU_B,nU_B
!write(u6,*) ' U_C',iU_C,nU_C

! Form the stabilizer for ((iAtom,mAtom),jAtom)

if ((iAtom == mAtom) .and. (iAtom == jAtom)) then

  ! A-R(A)-T(A)

  iOp_ER = ieor(iOp_E,iOp_R)
  iU_AB = ior(iU_A,iUR(iOp_ER,iU_A))
  iOp_ET = ieor(iOp_E,iOp_T)
  iU_ABC = ior(iU_AB,iUR(iOp_ET,iU_C))
else if (iAtom == jAtom) then

  ! A-R(B)-T(A)

  iOp_ET = ieor(iOp_E,iOp_T)
  iU_AC = ior(iU_A,iUR(iOp_ET,iU_C))
  iU_ABC = iand(iU_AC,iU_B)
else if (mAtom == jAtom) then

  ! A-R(B)-T(B)

  iU_ABC = iand(iU_A,iU_C)
else if (iAtom == mAtom) then

  ! A-R(A)-T(C)

  iU_ABC = iand(iU_A,iU_C)
else

  ! A-R(B)-T(C)

  iU_AB = iand(iU_A,iU_B)
  iU_ABC = iand(iU_AB,iU_C)
end if
nU_ABC = nU(iU_ABC)

! Compute the degeneracy of the angle

ideg = nSym/nU_ABC
D_Bend = real(iDeg,kind=wp)

!write(u6,*) ' D_Bend=',D_Bend

return

end function D_Bend
