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

use Slapaf_Info, only: dMass, Smmtrc

implicit real*8(a-h,o-z)
#include "real.fh"
! global arrays
real*8 FrcCrt(3,nsAtom), Coor(3,nsAtom)
integer Idntcl(nSupSy), iAtom(nsAtom)
! local arrays
real*8 cMass(3)
real*8 E(3)

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
    FORCE = FORCE+PROJ*dble(iDeg(Coor(1,iCen)))
    IDIV = IDIV+iDeg(Coor(1,iCen))

    ! Save the unit vector of the new force of this center

    do L=1,3
      FrcCrt(L,ICEN) = E(L)
    end do
  end do

  ! Get the new force
  FORCE = FORCE/dble(IDIV)
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

contains

subroutine CofMss(Coor,nsAtom,cMass)
  !*********************************************************************
  !   Object: To calculate the molecular mass, the center of mass and  *
  !           move the coordinates so origo is the center of mass.     *
  !*********************************************************************

  real*8 COOR(3,nsAtom), cMass(3)
  integer i, j

  ! Calculate the molecular mass.

  TMass = Zero
  do I=1,nsAtom
    TMass = TMass+dMass(I)*dble(iDeg(Coor(1,i)))
  end do
  iCOM = -1
  if (TMass >= 1.d99) then
    do i=1,nsAtom
      if (dMass(i) == 1.d99) then
        iCOM = i
        Go To 99
      end if
    end do
  end if
99 continue

  ! calculate the center of mass

  cMass(:) = Zero
  ! Loop over the unique centers
  do i=1,nsAtom
    do j=1,3
      ! Add contribution
      if (Smmtrc(j,i)) cMass(j) = cMass(j)+dMass(i)*Coor(j,i)*dble(iDeg(Coor(1,i)))
    end do
  end do

  do i=1,3
    cMass(i) = cMass(i)/TMass
  end do
  if ((iCOM >= 1) .and. (iCOM <= nsAtom)) cMass(:) = Coor(:,iCom)

# ifdef _DEBUGPRINT_
  if (LWrite) write(6,100) (cMass(i),i=1,3),TMass
100 format(//,' Center of Mass (Bohr) ',3F10.5,/,' Molecular Mass   (au) ',1F15.5)
# endif
# ifdef _DO_NOT_

  ! translate the center of mass to origo

  do i=1,nsAtom
    do j=1,3
      Coor(j,i) = Coor(j,i)-cMass(j)
    end do
  end do
# endif

  return

end subroutine CofMss

end subroutine SupSym
