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

function spin(l,i_dim,m1,m2)
! Returns the value of the < m1 | S_l | m2 >

use Constants, only: Zero, One, Half, cZero, cOne, Onei
use Definitions, only: wp, iwp, u6

implicit none
complex(kind=wp) :: spin
integer(kind=iwp), intent(in) :: l, i_dim, m1, m2
integer(kind=iwp) :: ipar
real(kind=wp) :: D, F, MM1, MM2, R, S

spin = cZero
ipar = mod(i_dim,2)
S = real(i_dim-1,kind=wp)*Half
R = Zero

#ifdef _DEBUGPRINT_
write(u6,'(A,4I3)') 'l,i_dim,m1,m2=',l,i_dim,m1,m2
write(u6,'(A,I3,F8.3)') ' ipar,  S  =',ipar,S
#endif

! set up the true values of the MM1 and MM2
! i.e. the projections of the spin momentum S
if (ipar == 0) then
  if (m1 < 0) then
    MM1 = real(m1,kind=wp)+Half
  else
    MM1 = real(m1,kind=wp)-Half
  end if

  if (m2 < 0) then
    MM2 = real(m2,kind=wp)+Half
  else
    MM2 = real(m2,kind=wp)-Half
  end if
else
  MM1 = real(m1,kind=wp)
  MM2 = real(m2,kind=wp)
end if
#ifdef _DEBUGPRINT_
write(u6,'(A,2F8.3)') 'mm1,mm2=',mm1,mm2
call xFlush(u6)
#endif
! compute the value of the matrix element, for each possible cartesian component.
! l=1  => X
! l=2  => Y
! l=3  => Z

if (l == 1) then

  if (MM1-One == MM2) then
    D = S+MM1
    F = S-MM1+One
    R = Half*sqrt(D*F)

    spin = R*cOne

  else if (MM1+One == MM2) then

    D = S-MM1
    F = S+MM1+One
    R = Half*sqrt(D*F)
    spin = R*cOne

  else
    spin = cZero
  end if
# ifdef _DEBUGPRINT_
  write(u6,'(A,2F8.3)') 'X, spin =',spin
  call xFlush(u6)
# endif

else if (l == 2) then

  if (MM1-One == MM2) then

    D = S+MM1
    F = S-MM1+One
    R = -Half*sqrt(D*F)
    spin = R*Onei

  else if (MM1+One == MM2) then

    D = S-MM1
    F = S+MM1+One
    R = Half*sqrt(D*F)
    spin = R*Onei

  else
    spin = cZero
  end if
# ifdef _DEBUGPRINT_
  write(u6,'(A,2F8.3)') 'Y, spin =',spin
  call xFlush(u6)
# endif

else if (l == 3) then

  if (MM1 /= MM2) then
    spin = cZero
  else
    spin = MM1*cOne
  end if
# ifdef _DEBUGPRINT_
  write(u6,'(A,2F8.3)') 'Z, spin =',spin
  call xFlush(u6)
# endif

else
  write(u6,'(A)') 'The spin function gives a wrong number'
  return
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'Upon Return:  spin =',spin
call xFlush(u6)
#endif

return

end function spin
