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

subroutine ModGauss(Z,A,Xi,w)

use Constants, only: rBohr

implicit none
real*8 Z, Xi, w
integer A
real*8 Facts(2,0:12), Errors(0:12), g(2), H(2,2), HInv(2,2), Step(2)
data Facts/0.0d0,0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0,2.0d0,0.0d0,-2.0d0,0.0d0,0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0,2.0d0,0.0d0,-2.0d0,1.0d0, &
           1.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,-1.0d0,-1.0d0/
real*8 A3, RMS, T, r_90, Thr, X, Delta_W, Delta_R, W0, R0, Det, e, f, r, x1, x2
integer MaxIter, Iter, i, iNeg
! Statement functions
f(x,w) = (1.0d0+w*x**2)*exp(-x**2)
R(w) = sqrt(2.0d0*RMS**2*(3.0d0*w+2.0d0)/(3.0d0*(2.0d0+5.0d0*w)))
E(x1,x2,w) = (f(x1,w)-0.9d0)**2+(f(x2,w)-0.1d0)**2

!                                                                      *
!***********************************************************************
!                                                                      *
!write(6,*) 'A=',A
A3 = 1.0d0*dble(A)**(1.0d0/3.0d0)
! DA Eq. 51
RMS = 0.836d0*A3+0.570d0 ! fm
RMS = RMS*1.0D-15/rBohr  ! bohr
!write(6,*) 'RMS:',RMS
w = 0.0d0
Xi = 1.0d0/R(w)**2
!if (A <= 9) write(6,*) 'Use the Gaussian model!'
if (A <= 9) return
!write(6,*) 'Xi:',Xi
!Xi = 1.5D0/RMS**2
!write(6,*) 'Xi:',Xi

T = 2.30d0           ! fm
T = T*1.0D-15/rBohr  ! bohr

! Start seeds

w = 0.5d0
r_90 = RMS/2.0d0
!                                                                      *
!***********************************************************************
!                                                                      *
! Start of iterations

Thr = 1.0D-7
MaxIter = 100
Iter = 0

777 continue

Iter = Iter+1
!write(6,*)
!write(6,*) 'Iteration:', Iter
!write(6,*) 'w.r_90:',w,r_90
Delta_w = 0.0001d0*w
Delta_r = 0.0001d0*r_90
!                                                                      *
!***********************************************************************
!                                                                      *
! The error is a function of r_90 and w

!RMS_s = RMS
!RMS = RMS*1.0D-15/rBohr  ! bohr
!write(6,*) 'RMS:',RMS
!Xi = 1.0D0/R(w)**2
!write(6,*) 'Xi=',Xi
!RMS = RMS_s
do i=0,12
  w0 = w+Facts(1,i)*Delta_w
  r0 = r_90+Facts(2,i)*Delta_r

  Errors(i) = E(r0/R(w0),(r0+T)/R(w0),w0)

  !if (i == 0) then
  !  write(6,*) 'r_90,f(x_90)=',r0,f(r0/R(w0),w0)
  !  write(6,*) 'r_10,f(x_10)=',r0+T,f((r0+T)/R(w0),w0)
  !  write(6,*) 'w,R(w)      =',w0,R(w0)
  !end if

end do
!write(6,*) 'Error=',Errors(0)
!call RecPrt('Errors',' ',Errors,13,1)
g(1) = (Errors(1)-Errors(2))/(2.0d0*Delta_w)
H(1,1) = (Errors(3)+Errors(4)-2.0d0*Errors(0))/(2.0d0*Delta_w)**2
g(2) = (Errors(5)-Errors(6))/(2.0d0*Delta_r)
H(2,2) = (Errors(7)+Errors(8)-2.0d0*Errors(0))/(2.0d0*Delta_r)**2
H(1,2) = (Errors(9)+Errors(12)-Errors(10)-Errors(11))/((2.0d0*Delta_w)*(2.0d0*Delta_r))
H(2,1) = H(1,2)
!call RecPrt('gradient',' ',g,1,2)
!call RecPrt('Hessian ',' ',H,2,2)

call DiagMtrx_x(H,2,iNeg)
!write(6,*) 'iNeg=',iNeg
!call RecPrt('Hessian ',' ',H,2,2)
call MInv(H,HInv,Det,2)
!call RecPrt('HInv ',' ',HInv,2,2)
Step(1) = HInv(1,1)*g(1)+HInv(1,2)*g(2)
Step(2) = HInv(2,1)*g(1)+HInv(2,2)*g(2)
!call RecPrt('Step',' ',Step,1,2)
Step(1) = sign(min(abs(Step(1)),0.1d0*w),Step(1))
Step(2) = sign(min(abs(Step(2)),0.1d0*r_90),Step(2))
!call RecPrt('Step',' ',Step,1,2)
w = w-Step(1)
r_90 = r_90-Step(2)

if ((Iter < MaxIter) .and. (Errors(0) > Thr)) Go To 777
w0 = w
r0 = r_90
!write(6,*)
!write(6,*) 'Iterations:',Iter
!write(6,*) 'Error=',Errors(0)
!write(6,*) 'r_90,f(x_90)=',r0,f(r0/R(w0),w0)
!write(6,*) 'r_10,f(x_10)=',r0+T,f((r0+T)/R(w0),w0)
!write(6,*) 'w,R(w)      =',w0,R(w0)
!write(6,*)
!                                                                      *
!***********************************************************************
!                                                                      *
Xi = 1.0d0/R(w)**2
w = w*Xi
!write(6,*) 'Xi:',Xi
!write(6,*) ' w:',w
!                                                                      *
!***********************************************************************
!                                                                      *
return
! Avoid unused argument warnings
if (.false.) call Unused_real(Z)

end subroutine ModGauss
