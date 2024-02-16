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

subroutine Dipol_Exchange(N1,N2,vec,dist,M1,M2,HDIP)
! this Subroutine computes the dipolar coupling between the two moments

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: N1, N2
real(kind=8), intent(in) :: vec(3), dist
complex(kind=8), intent(in) :: M1(3,N1,N1)
complex(kind=8), intent(in) :: M2(3,N2,N2)
complex(kind=8), intent(out) :: HDIP(N1,N1,N2,N2)
! local variables
integer :: m, i1, j1, i2, j2
real(kind=8) :: MB2
complex(kind=8) :: p2a, p2b, p1, HL, d3, mb2c, threeC, vec1(3)

MB2 = 0.4329701512063995_wp ! in cm-1*T-1   -- the value of (mu_Bohr*mu_Bohr)/(angstrom^3)  in cm-1

if ((N1 <= 0) .or. (N2 <= 0)) return
call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,HDIP,1)
if (dist == 0.0_wp) then
  write(6,'(A)') 'DIPOL_EXCHANGE::  dist = 0'
  write(6,'(A)') 'this is not normal. Stop.'
  return
end if
! switch to Complex arithmetic:
! kind=8, complex double precision
d3 = (0.0_wp,0.0_wp)
mb2c = (0.0_wp,0.0_wp)
threeC = (3.0_wp,0.0_wp)
d3 = cmplx(dist*dist*dist,0.0_wp,wp)
mb2c = cmplx(MB2,0.0_wp,wp)
do m=1,3
  vec1(m) = cmplx(vec(m),0.0_wp,wp)
end do
! calculate the dipolar coupling:
do i1=1,N1
  do j1=1,N1
    do i2=1,N2
      do j2=1,N2
        p2a = (0.0_wp,0.0_wp)
        p2b = (0.0_wp,0.0_wp)
        p1 = (0.0_wp,0.0_wp)
        HL = (0.0_wp,0.0_wp)
        do m=1,3
          p2a = p2a+M1(m,i1,j1)*vec1(m)
          p2b = p2b+M2(m,i2,j2)*vec1(m)
          p1 = p1+M1(m,i1,j1)*M2(m,i2,j2)
        end do
        HL = p1-threeC*p2a*p2b
        HDIP(i1,j1,i2,j2) = mb2c*HL/d3
      end do
    end do
  end do
end do

return

end subroutine Dipol_Exchange
