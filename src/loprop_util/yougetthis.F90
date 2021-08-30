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

subroutine YouGetThis(EC,Pot_Expo,Pot_Point,Pot_Fac,Diffed,MP,lMax,lMaxF,nij,LuYou)

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lMax, lMaxF, nij, LuYou
real(kind=wp), intent(in) :: EC(3,nij), Pot_Expo(nij*4), Pot_Point(nij), Pot_Fac(nij*4), MP(nij,*)
logical(kind=iwp), intent(in) :: Diffed(nij*4)
integer(kind=iwp) :: i, k, kauntA, kk, l, nS, nT

! Number of centres and maximal angular momentum.

write(LuYou,101) nij
write(LuYou,102) lMax,lMaxF
kauntA = 0

! Loop over centres.

do i=1,nij
  kauntA = kauntA+1

  ! Coordinates of centre.

  write(LuYou,103) (EC(k,i),k=1,3)
  do l=0,lMax
    nS = l*(l+1)*(l+2)/6
    nT = (l+1)*(l+2)*(l+3)/6
    if (l <= lMaxF) then
      if (Diffed(2*(kauntA-1)+l+1)) then

        ! Factor and exponent.

        write(LuYou,104) Two*Pot_Expo(2*(kauntA-1)+l+1)
        write(LuYou,105) (Pot_Fac(4*(kauntA-1)+kk),kk=nS+1,nT)
      else

        ! Factor and dummy-exponent, if this multipole should not be
        ! made diffuse.

        write(LuYou,104) -7.91204_wp
        write(LuYou,105) (Pot_Fac(4*(kauntA-1)+kk),kk=nS+1,nT)
      end if
    else

      ! The pure multipole, which under no circumstance can be diffuse.

      write(LuYou,104) -7.91204_wp
      write(LuYou,105) (MP(kauntA,kk),kk=nS+1,nT)
    end if
  end do

  ! Throw in the point-part.

  write(LuYou,104) Pot_Point(kauntA)
end do

101 format(I5)
102 format(2I5)
103 format(3(F20.14))
104 format(F20.14)
105 format(3(F20.14))

return

end subroutine YouGetThis
