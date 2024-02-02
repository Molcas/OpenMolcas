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

subroutine prbar(ist,s1,s2,M)

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: ist
character(len=5), intent(in) :: s1, s2
complex(kind=8), intent(in) :: M(3)
! local
character(len=30) :: fx, fy, fz
character(len=40) :: f1, f2
real(kind=8) :: R

write(fx,'(i2,5a)') ist,'. | <',s1,' | mu_X |',s2,' > |'
write(fy,'(i2,5a)') ist,'. | <',s1,' | mu_Y |',s2,' > |'
write(fz,'(i2,5a)') ist,'. | <',s1,' | mu_Z |',s2,' > |'
R = 0.0_wp
R = (abs(M(1))+abs(M(2))+abs(M(3)))/3.0_wp

f1 = '(2x,a,2ES19.11,1x,A,      23x,A)'
f2 = '(2x,a,2ES19.11,1x,A,ES22.12,1x,A)'
write(6,f1) fx,M(1),'|','|'
write(6,f2) fy,M(2),'|',R,'|'
write(6,f1) fz,M(3),'|','|'

return

end subroutine prbar
