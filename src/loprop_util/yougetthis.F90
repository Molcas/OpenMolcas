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

subroutine YouGetThis(nAt,EC,Pot_Expo,Pot_Point,Pot_Fac,Diffed,ipMP,lMax,lMaxF,nij,LuYou)

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
dimension EC(3,nij), Pot_Expo(nij*4), Pot_Point(nij), Pot_Fac(nij*4)
logical Diffed(nij*4)

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

        write(LuYou,104) 2.0d0*Pot_Expo(2*(kauntA-1)+l+1)
        write(LuYou,105) (Pot_Fac(4*(kauntA-1)+kk),kk=nS+1,nT)
      else

        ! Factor and dummy-exponent, if this multipole should not be
        ! made diffuse.

        write(LuYou,104)-7.91204d0
        write(LuYou,105) (Pot_Fac(4*(kauntA-1)+kk),kk=nS+1,nT)
      end if
    else

      ! The pure multipole, which under no circumstance can be diffuse.

      write(LuYou,104)-7.91204d0
      write(LuYou,105) (Work(ipMP+nij*kk+kauntA-1),kk=nS,nT-1)
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
! Avoid unused argument warnings
if (.false.) call Unused_integer(nAt)

end subroutine YouGetThis
