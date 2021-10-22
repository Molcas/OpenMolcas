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

subroutine Find_Dipole_Center(q_A,q_B,Dipole_A,Dipole_B,qn_A,qn_B,R_A,R_B,t,iPlot)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: q_A, q_B, Dipole_A, Dipole_B, qn_A, qn_B, R_A, R_B
real(kind=wp), intent(out) :: t
integer(kind=iwp), intent(in) :: iPlot
integer(kind=iwp) :: i
real(kind=wp) :: ax, bx, cx, Delta, E, fa, fb, fc, golden, R, r_best, r_best2, R_Test
integer(kind=iwp) :: nSteps = 100
real(kind=wp), parameter :: E_Threshold = 1.0e-12_wp, R_Threshold = 1.0e-12_wp
logical(kind=iwp), parameter :: Absolute = .false.
real(kind=wp), external :: Golden2, Multipole_Expansion

R = R_B-R_A

Delta = R/real(nSteps+1,kind=wp)
R_Test = R_A+Delta
if (iPlot == 1) then
  write(u6,*) 'Electronic contributions: q_A, q_B = ',q_a,q_b
  do i=1,nSteps
    R_Test = R_A+Delta*real(i,kind=wp)
    E = Multipole_Expansion(q_A,q_B,Dipole_A,Dipole_B,R_A,R_B,R_Test,Absolute)
    write(u6,'(1X,A,F6.3,1X,F20.12)') 'R, E = ',R_test,E
    call xFlush(u6)
  end do
  write(u6,*) 'Nuclear contributions: q_A, q_B = ',qn_A,qn_B
  do i=1,nSteps
    R_Test = R_A+Delta*real(i,kind=wp)
    E = Multipole_Expansion(qn_A,qn_B,Zero,Zero,R_A,R_B,R_Test,Absolute)
    write(u6,'(1X,A,F6.3,1X,F20.12)') 'R, E = ',R_test,E
    call xFlush(u6)
  end do
  write(u6,*) 'Total contributions: q_A, q_B = ',q_a+qn_a,q_b+qn_b
  do i=1,nSteps
    R_Test = R_A+Delta*real(i,kind=wp)
    E = Multipole_Expansion(q_A+qn_A,q_B+qn_B,Dipole_A,Dipole_B,R_A,R_B,R_Test,Absolute)
    write(u6,'(1X,A,F6.3,1X,F20.12)') 'R, E = ',R_test,E
    call xFlush(u6)
  end do
end if

ax = (R_A+R_B)*Half+Delta
bx = (R_A+R_B)*Half-Delta
! First make an initial bracketing of the minima
call mnBrak2(ax,bx,cx,fa,fb,fc,multipole_expansion,q_a,q_b,dipole_a,dipole_b,r_a,r_b)
! Use Golden to find the minima
golden = golden2(ax,bx,cx,multipole_expansion,r_threshold,e_threshold,r_best,q_a,q_b,dipole_a,dipole_b,r_a,r_b)
t = (r_best-half*R)/R
write(u6,'(A,3F18.10)') 't_el , r_best, golden = ',t,r_best,golden
call xflush(u6)

ax = (R_A+R_B)*Half+Delta
bx = (R_A+R_B)*Half-Delta
!ax = -0.245_wp
!bx = -0.263_wp
! First make an initial bracketing of the minima
call mnBrak2(ax,bx,cx,fa,fb,fc,multipole_expansion,qn_a,qn_b,zero,zero,r_a,r_b)
! Use Golden to find the minima
golden = golden2(ax,bx,cx,multipole_expansion,r_threshold,e_threshold,r_best2,qn_a,qn_b,zero,zero,r_a,r_b)
t = (r_best2-half*R)/R
write(u6,'(A,3F18.10)') 't_nuc, r_best, golden = ',t,r_best2,golden
call xflush(u6)
r_best = (abs(q_a+q_b)*R_best+abs(qn_a+qn_b)*R_best2)/(abs(q_a+q_b)+abs(qn_a+qn_b))
t = (r_best-half*R)/R
write(u6,'(A,3F18.10)') 't_fit, r_best, golden = ',t,r_best,golden
call xflush(u6)

!?? ax = (R_A+R_B)*Half + Delta
!?? bx = (R_A+R_B)*Half - Delta
!?? !ax = -0.245_wp
!?? !bx = -0.263_wp
!?? ! First make an initial bracketing of the minima
!?? call mnBrak2(ax,bx,cx,fa,fb,fc,multipole_expansion,q_a+qn_a,q_b+qn_b,dipole_a,dipole_b,r_a,r_b)
!?? ! Use Golden to find the minima
!?? golden = golden2(ax,bx,cx,multipole_expansion,r_threshold,e_threshold,r_best,q_a+qn_a,q_b+qn_b,dipole_a,dipole_b,r_a,r_b)
!?? t = (r_best-half*R)/R
!?? write(u6,'(A,3F18.10)') 't_sum, r_best, golden = ',t, r_best, golden
!?? call xflush(u6)

return

end subroutine Find_Dipole_Center
