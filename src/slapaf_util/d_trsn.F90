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

function D_Trsn(Ind,iOp_,nSym)

use Slapaf_Info, only: jStab, nStab
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: D_Trsn
integer(kind=iwp), intent(in) :: Ind(4), iOp_(4), nSym
integer(kind=iwp) :: iAtom, iDeg, iOp_E, iOp_ER, iOp_ES, iOp_ET, iOp_R, iOp_S, iOp_T, iOp_TS, iU_A, iU_AB, iU_ABCD, iU_B, iU_C, &
                     iU_CD, iU_D, jAtom, kAtom, lAtom, nU_A, nU_ABCD, nU_B, nU_C, nU_D
integer(kind=iwp), external :: iU, iUR, nU

!                                                                      *
!***********************************************************************
!                                                                      *
! T O R S I O N S ((iAtom,R(jAtom)),T(kAtom,S(lAtom)))

D_Trsn = Zero

iAtom = Ind(1)
jAtom = Ind(2)
kAtom = Ind(3)
lAtom = Ind(4)
iOp_E = iOp_(1)
iOp_R = iOp_(2)
iOp_T = iOp_(3)
iOp_TS = iOp_(4)
iOp_S = ieor(iOp_T,iOp_TS)
!write(u6,*) iAtom,jAtom,kAtom,lAtom
!write(u6,*) iOp_E,iOp_R,iOp_T,iOp_TS,iOp_S

nU_A = nStab(iAtom)
iU_A = iU(jStab(0,iAtom),nU_A)
nU_B = nStab(jAtom)
iU_B = iU(jStab(0,jAtom),nU_B)
nU_C = nStab(kAtom)
iU_C = iU(jStab(0,kAtom),nU_C)
nU_D = nStab(lAtom)
iU_D = iU(jStab(0,lAtom),nU_D)

!write(u6,*) ' U_A:',iU_A,nU_A
!write(u6,*) ' U_B:',iU_B,nU_B
!write(u6,*) ' U_C:',iU_C,nU_C
!write(u6,*) ' U_D:',iU_D,nU_D

! Form stabilizer for (iAtom,jAtom)

iOp_ER = ieor(iOp_E,iOp_R)
if (iAtom == jAtom) then
  iU_AB = ior(iU_A,iUR(iOp_ER,iU_B))
else
  iU_AB = iand(iU_A,iU_B)
end if
!write(u6,*) iAtom == jAtom
!write(u6,*) ' U_AB:',iU_AB,nU(iU_AB)

! Form stabilizer for (kAtom,lAtom)

iOp_ES = ieor(iOp_E,iOp_S)
if (kAtom == lAtom) then
  iU_CD = ior(iU_C,iUR(iOp_ES,iU_D))
else
  iU_CD = iand(iU_C,iU_D)
end if
!write(u6,*) kAtom == lAtom
!write(u6,*) ' U_CD:',iU_CD,nU(iU_CD)

! Form the stabilizer for the torsion

if ((iAtom /= lAtom) .or. (jAtom /= kAtom) .or. ((iAtom == lAtom) .and. (jAtom == kAtom) .and. (iOp_ER /= iOp_ES))) then
  iU_ABCD = iand(iU_AB,iU_CD)
  !write(u6,*) ' Ops!'
else
  iOp_ET = ieor(iOp_E,iOp_T)
  iU_ABCD = ior(iU_AB,iUR(iOp_ET,iU_CD))
end if
!write(u6,*) iAtom /= lAtom
!write(u6,*) jAtom /= kAtom
!write(u6,*) ((iAtom == lAtom) .and. (jAtom == kAtom) .and. (iOp_ER /= iOp_ES))
nU_ABCD = nU(iU_ABCD)
!write(u6,*) ' U_ABCD:',iU_ABCD,nU(iU_ABCD)

! Compute the degeneracy of the torsion

iDeg = nSym/nU_ABCD
D_Trsn = real(iDeg,kind=wp)

!write(u6,*) ' D_Trsn=',D_Trsn

return

end function D_Trsn
