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
integer(kind=iwp) :: nsAtom, nSupSy, Idntcl(nSupSy), iAtom(nsAtom)
real(kind=wp) :: FrcCrt(3,nsAtom), Coor(3,nsAtom)
integer(kind=iwp) :: I, ICEN, IDIV, J, K, L
real(kind=wp) :: ABSE, cMass(3), E(3), E2, FORCE, PROJ
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

    E2 = Zero
    do L=1,3
      E(L) = Coor(L,ICEN)-cMass(L)
      E2 = E2+E(L)*E(L)
    end do
    ABSE = sqrt(E2)
    ! Normalize the vector.
    do L=1,3
      E(L) = E(L)/ABSE
    end do

    ! Project the force of the ICEN'th center on this unit vector

    PROJ = DDot_(3,E,1,FrcCrt(1,ICEN),1)

    ! Sum forces which should be equal
    FORCE = FORCE+PROJ*real(iDeg(Coor(1,iCen)),kind=wp)
    IDIV = IDIV+iDeg(Coor(1,iCen))

    ! Save the unit vector of the new force of this center

    do L=1,3
      FrcCrt(L,ICEN) = E(L)
    end do
  end do

  ! Get the new force
  FORCE = FORCE/real(IDIV,kind=wp)
  K = K-Idntcl(I)

  ! Symmetrize the forces

  do J=1,Idntcl(I)
    K = K+1
    ICEN = iAtom(K)

    do L=1,3
      FrcCrt(L,ICEN) = FrcCrt(L,ICEN)*FORCE
    end do
  end do
end do

end subroutine SupSym
