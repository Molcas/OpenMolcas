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

subroutine SupSym(FrcCrt,nsAtom,Coor,nSupSy,Idntcl,iAtom)
!***********************************************************************
!                                                                      *
!     Object: To constrain forces to higher symmetry.                  *
!                                                                      *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nsAtom, nSupSy, Idntcl(nSupSy), iAtom(nsAtom)
real(kind=wp), intent(inout) :: FrcCrt(3,nsAtom)
real(kind=wp), intent(in) :: Coor(3,nsAtom)
integer(kind=iwp) :: I, ICEN, IDIV, J, K
real(kind=wp) :: cMass(3), E(3), FORCE, PROJ
integer(kind=iwp), external :: iDeg
real(kind=wp), external :: DDot_

call CofMss(Coor,nsAtom,cMass)

! Loop over groups of centers which are identical.

K = 0
do I=1,nSupSy
  FORCE = Zero

  ! loop over centers which are identical.

  IDIV = 0
  do J=1,Idntcl(I)
    K = K+1
    ! Get the number of this center
    ICEN = iAtom(K)

    ! Get an unit vector in the direction from center of mass to
    ! the ICEN'th center.

    E(:) = Coor(:,ICEN)-cMass(:)
    ! Normalize the vector.
    E(:) = E(:)/sqrt(sum(E(:)**2))

    ! Project the force of the ICEN'th center on this unit vector

    PROJ = DDot_(3,E,1,FrcCrt(1,ICEN),1)

    ! Sum forces which should be equal
    FORCE = FORCE+PROJ*real(iDeg(Coor(1,iCen)),kind=wp)
    IDIV = IDIV+iDeg(Coor(1,iCen))

    ! Save the unit vector of the new force of this center

    FrcCrt(:,ICEN) = E(:)
  end do

  ! Get the new force
  FORCE = FORCE/real(IDIV,kind=wp)
  K = K-Idntcl(I)

  ! Symmetrize the forces

  do J=1,Idntcl(I)
    K = K+1
    ICEN = iAtom(K)

    FrcCrt(:,ICEN) = FrcCrt(:,ICEN)*FORCE
  end do
end do

end subroutine SupSym
