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

subroutine ModGauss(A,Xi,w)

use Constants, only: Zero, One, Two, Three, Half, rBohr
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: A
real(kind=wp), intent(out) :: Xi, w
integer(kind=iwp) :: i, iNeg, Iter, MaxIter
real(kind=wp) :: A3, Delta_R, Delta_W, Det, Errors(0:12), g(2), H(2,2), HInv(2,2), R0, r_90, RMS, Step(2), T, Thr, W0
real(kind=wp), parameter :: Facts(2,0:12) = reshape([Zero,Zero,One,Zero,-One,Zero,Two,Zero,-Two,Zero,Zero,One,Zero,-One,Zero,Two, &
                                                     Zero,-Two,One,One,-One,One,One,-One,-One,-One],[2,13])

!                                                                      *
!***********************************************************************
!                                                                      *
!write(u6,*) 'A=',A
A3 = real(A,kind=wp)**(One/Three)
! DA Eq. 51 (See NucExp)
RMS = 0.836_wp*A3+0.570_wp  ! fm
RMS = RMS*1.0e-15_wp/rBohr  ! bohr
!write(u6,*) 'RMS:',RMS
w = Zero
Xi = One/R(w)**2
!if (A <= 9) write(u6,*) 'Use the Gaussian model!'
if (A <= 9) return
!write(u6,*) 'Xi:',Xi
!Xi = OneHalf/RMS**2
!write(u6,*) 'Xi:',Xi

T = 2.3_wp              ! fm
T = T*1.0e-15_wp/rBohr  ! bohr

! Start seeds

w = Half
r_90 = Half*RMS
!                                                                      *
!***********************************************************************
!                                                                      *
! Start of iterations

Thr = 1.0e-7_wp
MaxIter = 100

do Iter=1,MaxIter
  !write(u6,*)
  !write(u6,*) 'Iteration:', Iter
  !write(u6,*) 'w.r_90:',w,r_90
  Delta_w = 1.0e-4_wp*w
  Delta_r = 1.0e-4_wp*r_90
  !                                                                      *
  !***********************************************************************
  !                                                                      *
  ! The error is a function of r_90 and w

  !RMS_s = RMS
  !RMS = RMS*1.0e-15_wp/rBohr  ! bohr
  !write(u6,*) 'RMS:',RMS
  !Xi = One/R(w)**2
  !write(u6,*) 'Xi=',Xi
  !RMS = RMS_s
  do i=0,12
    w0 = w+Facts(1,i)*Delta_w
    r0 = r_90+Facts(2,i)*Delta_r

    Errors(i) = (f(r0/R(w0),w0)-0.9_wp)**2+(f((r0+T)/R(w0),w0)-0.1_wp)**2

    !if (i == 0) then
    !  write(u6,*) 'r_90,f(x_90)=',r0,f(r0/R(w0),w0)
    !  write(u6,*) 'r_10,f(x_10)=',r0+T,f((r0+T)/R(w0),w0)
    !  write(u6,*) 'w,R(w)      =',w0,R(w0)
    !end if

  end do
  !write(u6,*) 'Error=',Errors(0)
  !call RecPrt('Errors',' ',Errors,13,1)
  g(1) = (Errors(1)-Errors(2))/(Two*Delta_w)
  H(1,1) = (Errors(3)+Errors(4)-Two*Errors(0))/(Two*Delta_w)**2
  g(2) = (Errors(5)-Errors(6))/(Two*Delta_r)
  H(2,2) = (Errors(7)+Errors(8)-Two*Errors(0))/(Two*Delta_r)**2
  H(1,2) = (Errors(9)+Errors(12)-Errors(10)-Errors(11))/((Two*Delta_w)*(Two*Delta_r))
  H(2,1) = H(1,2)
  !call RecPrt('gradient',' ',g,1,2)
  !call RecPrt('Hessian ',' ',H,2,2)

  call DiagMtrx_x(H,2,iNeg)
  !write(u6,*) 'iNeg=',iNeg
  !call RecPrt('Hessian ',' ',H,2,2)
  call MInv(H,HInv,Det,2)
  !call RecPrt('HInv ',' ',HInv,2,2)
  Step(1) = HInv(1,1)*g(1)+HInv(1,2)*g(2)
  Step(2) = HInv(2,1)*g(1)+HInv(2,2)*g(2)
  !call RecPrt('Step',' ',Step,1,2)
  Step(1) = sign(min(abs(Step(1)),0.1_wp*w),Step(1))
  Step(2) = sign(min(abs(Step(2)),0.1_wp*r_90),Step(2))
  !call RecPrt('Step',' ',Step,1,2)
  w = w-Step(1)
  r_90 = r_90-Step(2)

  if (Errors(0) <= Thr) exit
end do
w0 = w
r0 = r_90
!write(u6,*)
!write(u6,*) 'Iterations:',Iter
!write(u6,*) 'Error=',Errors(0)
!write(u6,*) 'r_90,f(x_90)=',r0,f(r0/R(w0),w0)
!write(u6,*) 'r_10,f(x_10)=',r0+T,f((r0+T)/R(w0),w0)
!write(u6,*) 'w,R(w)      =',w0,R(w0)
!write(u6,*)
!                                                                      *
!***********************************************************************
!                                                                      *
Xi = One/R(w)**2
w = w*Xi
!write(u6,*) 'Xi:',Xi
!write(u6,*) ' w:',w
!                                                                      *
!***********************************************************************
!                                                                      *
return

contains

pure function f(x,w)

  real(kind=wp) :: f
  real(kind=wp), intent(in) :: x, w

  f = (One+w*x**2)*exp(-x**2)

end function

pure function R(w)

  use Constants, only: Five

  real(kind=wp) :: R
  real(kind=wp), intent(in) :: w

  R = sqrt(Two*RMS**2*(Three*w+Two)/(Three*(Two+Five*w)))

end function

end subroutine ModGauss
