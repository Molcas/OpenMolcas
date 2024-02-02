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

complex*16 function spin(l,dim,m1,m2)
! Returns the value of the < m1 | S_l | m2 >

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: m1, m2, dim, l
integer :: ipar
real(kind=8) :: S, MM1, MM2, R, D, F
logical :: dbg

dbg = .false.

spin = 0.0_wp
ipar = mod(dim,2)
S = dble(dim-1)/2.0_wp
R = 0.0_wp

if (dbg) then
  write(6,'(A,4I3)') 'l,dim,m1,m2=',l,dim,m1,m2
  write(6,'(A,I3,F8.3)') ' ipar,  S  =',ipar,S
end if

! set up the true values of the MM1 and MM2
! i.e. the projections of the spin momentum S
if (ipar == 0) then
  if (m1 < 0) then
    MM1 = dble(m1)+0.5_wp
  else
    MM1 = dble(m1)-0.5_wp
  end if

  if (m2 < 0) then
    MM2 = dble(m2)+0.5_wp
  else
    MM2 = dble(m2)-0.5_wp
  end if
else
  MM1 = dble(m1)
  MM2 = dble(m2)
end if
if (dbg) write(6,'(A,2F8.3)') 'mm1,mm2=',mm1,mm2
call xFlush(6)
! compute the value of the matrix element, for each possible cartesian component.
! l=1  => X
! l=2  => Y
! l=3  => Z

if (l == 1) then

  if ((MM1-1.0_wp) == MM2) then
    D = S+MM1
    F = S-MM1+1.0_wp
    R = 0.5_wp*sqrt(D*F)

    spin = cmplx(R,0.0_wp,wp)

  else if ((MM1+1.0_wp) == MM2) then

    D = S-MM1
    F = S+MM1+1.0_wp
    R = 0.5_wp*sqrt(D*F)
    spin = cmplx(R,0.0_wp,wp)

  else
    spin = cmplx(0.0_wp,0.0_wp,wp)
  end if
  if (dbg) write(6,'(A,2F8.3)') 'X, spin =',spin
  call xFlush(6)

else if (l == 2) then

  if ((MM1-1.0_wp) == MM2) then

    D = S+MM1
    F = S-MM1+1.0_wp
    R = -0.5_wp*sqrt(D*F)
    spin = cmplx(0.0_wp,R,wp)

  else if ((MM1+1.0_wp) == MM2) then

    D = S-MM1
    F = S+MM1+1.0_wp
    R = 0.5_wp*sqrt(D*F)
    spin = cmplx(0.0_wp,R,wp)

  else
    spin = cmplx(0.0_wp,0.0_wp,wp)
  end if
  if (dbg) write(6,'(A,2F8.3)') 'Y, spin =',spin
  call xFlush(6)

else if (l == 3) then

  if (MM1 /= MM2) then
    spin = cmplx(0.0_wp,0.0_wp,wp)
  else
    spin = cmplx(MM1,0.0_wp,wp)
  end if
  if (dbg) write(6,'(A,2F8.3)') 'Z, spin =',spin
  call xFlush(6)

else
  write(6,'(A)') 'The spin function gives a wrong number'
  return
end if

if (dbg) write(6,*) 'Upon Return:  spin =',spin
call xFlush(6)

return

end function spin
