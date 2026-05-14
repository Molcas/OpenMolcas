!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021-2023, Vladislav Kochetov                          *
!***********************************************************************
! set of routines implementing different integrators, mainly of
! Runge-Kutta family.

module integrators

use rhodyn_data, only: dt, d, ak1, ak2, ak3, ak4, ak5, ak6, timestep, len_sph, midk1, midk2, midk3, midk4
use Definitions, only: wp

implicit none
private

public :: classic_rk4, rk4, rk5, rk45, rkck, rk4_sph

contains

subroutine classic_rk4(t0,y)
  !*********************************************************************
  ! convenient Runge-Kutta method of 4th order
  !*********************************************************************

  use Constants, only: Six, Half

  real(kind=wp), intent(in) :: t0
  complex(kind=wp), intent(inout) :: y(d,d)

  call equation(t0,y,ak1)
  call equation(t0+Half*timestep,y+Half*timestep*ak1,ak2)
  call equation(t0+Half*timestep,y+Half*timestep*ak2,ak3)
  call equation(t0+timestep,y+timestep*ak3,ak4)
  y(:,:) = y+timestep/Six*(ak1+2*ak2+2*ak3+ak4)

end subroutine classic_rk4

subroutine rk4(t0,y)
  !*********************************************************************
  ! Runge-Kutta method of 4th order with proper adjusted midpoints
  !*********************************************************************

  real(kind=wp), intent(in) :: t0
  complex(kind=wp), intent(inout) :: y(d,d)
  real(kind=wp) :: x
  real(kind=wp), parameter :: a2 = 0.25_wp, &
                              a3 = 0.375_wp, &
                              a4 = 12.0_wp/13.0_wp, &
                              c21 = 0.25_wp, &
                              c31 = 0.09375_wp, &
                              c32 = 0.28125_wp, &
                              c41 = 1932.0_wp/2197.0_wp, &
                              c42 = -7200.0_wp/2197.0_wp, &
                              c43 = 7296.0_wp/2197.0_wp, &
                              c51 = 439.0_wp/216.0_wp, &
                              c52 = -8.0_wp, &
                              c53 = 3680.0_wp/513.0_wp, &
                              c54 = -845.0_wp/4104.0_wp, &
                              c1 = 25.0_wp/216.0_wp, &
                              c3 = 1408.0_wp/2565.0_wp, &
                              c4 = 2197.0_wp/4104.0_wp, &
                              c5 = -0.2_wp

  x = t0
  call equation(x,y,ak1)
  x = t0+a2*timestep
  call equation(x,y+ak1*timestep*c21,ak2)
  x = t0+a3*timestep
  call equation(x,y+timestep*(c31*ak1+c32*ak2),ak3)
  x = t0+a4*timestep
  call equation(x,y+timestep*(c41*ak1+c42*ak2+c43*ak3),ak4)
  x = t0+timestep
  call equation(x,y+timestep*(c51*ak1+c52*ak2+c53*ak3+c54*ak4),ak5)
  y(:,:) = y+(c1*ak1+c3*ak3+c4*ak4+c5*ak5)*timestep

end subroutine rk4

subroutine rk5(t0,y)
  !*********************************************************************
  ! Runge-Kutta method of 5th order
  !*********************************************************************

  use Constants, only: Half

  real(kind=wp), intent(in) :: t0
  complex(kind=wp), intent(inout) :: y(d,d)
  real(kind=wp) :: x
  real(kind=wp), parameter :: a2 = 0.25_wp, &
                              a3 = 0.375_wp, &
                              a4 = 12.0_wp/13.0_wp, &
                              c21 = 0.25_wp, &
                              c31 = 0.09375_wp, &
                              c32 = 0.28125_wp, &
                              c41 = 1923.0_wp/2197.0_wp, &
                              c42 = -7200.0_wp/2197.0_wp, &
                              c43 = 7296.0_wp/2197.0_wp, &
                              c51 = 439.0_wp/216.0_wp, &
                              c52 = -8.0_wp, &
                              c53 = 3680.0_wp/513.0_wp, &
                              c54 = -845.0_wp/4104.0_wp, &
                              c61 = -8.0_wp/27.0_wp, &
                              c62 = 2.0_wp, &
                              c63 = -3544.0_wp/2565.0_wp, &
                              c64 = 1859.0_wp/4104.0_wp, &
                              c65 = -0.275_wp, &
                              c1 = 16.0_wp/135.0_wp, &
                              c3 = 6656.0_wp/12825.0_wp, &
                              c4 = 28561.0_wp/56430.0_wp, &
                              c5 = -0.18_wp, &
                              c6 = 2.0_wp/55.0_wp

  x = t0
  call equation(x,y,ak1)
  x = t0+a2*timestep
  call equation(x,y+timestep*c21*ak1,ak2)
  x = t0+a3*timestep
  call equation(x,y+timestep*(c31*ak1+c32*ak2),ak3)
  x = t0+a4*timestep
  call equation(x,y+timestep*(c41*ak1+c42*ak2+c43*ak3),ak4)
  x = t0+timestep
  call equation(x,y+timestep*(c51*ak1+c52*ak2+c53*ak3+c54*ak4),ak5)
  x = t0+Half*timestep
  call equation(x,y+timestep*(c61*ak1+c62*ak2+c63*ak3+c64*ak4+c65*ak5),ak6)
  y(:,:) = y+(c1*ak1+c3*ak3+c4*ak4+c5*ak5+c6*ak6)*timestep

end subroutine rk5

subroutine rk45(t0,y,err)
  !*********************************************************************
  ! Runge-Kutta-Fehlberg 4(5) integration algorithm
  !*********************************************************************
  real(kind=wp), intent(in) :: t0
  complex(kind=wp), intent(inout) :: y(d,d)
  real(kind=wp), intent(out) :: err
  real(kind=wp) :: x
  real(kind=wp), parameter :: a2 = 0.25_wp, &
                              a3 = 3.0_wp/8.0_wp, &
                              a4 = 12.0_wp/13.0_wp, &
                              a6 = 0.5_wp, &
                              c21 = 0.25_wp, &
                              c31 = 3.0_wp/32.0_wp, &
                              c32 = 9.0_wp/32.0_wp, &
                              c41 = 1932.0_wp/2197.0_wp, &
                              c42 = -7200.0_wp/2197.0_wp, &
                              c43 = 7296.0_wp/2197.0_wp, &
                              c51 = 439.0_wp/216.0_wp, &
                              c52 = -8.0_wp, &
                              c53 = 3680.0_wp/513.0_wp, &
                              c54 = -845.0_wp/4104.0_wp, &
                              c61 = -8.0_wp/27.0_wp, &
                              c62 = 2.0_wp, &
                              c63 = -3544.0/2565.0_wp, &
                              c64 = 1859.0/4104.0_wp, &
                              c65 = -0.275_wp, &
                              c1 = 16.0_wp/135.0_wp, &
                              c3 = 6656.0_wp/12825.0_wp, &
                              c4 = 28561.0_wp/56430.0_wp, &
                              c5 = -0.18_wp, &
                              c6 = 2.0_wp/55.0_wp, &
                              dc1 = c1-25.0_wp/216.0_wp, &
                              dc3 = c3-1408.0_wp/2565.0_wp, &
                              dc4 = c4-2197.0_wp/4104.0_wp, &
                              dc5 = c5+0.2_wp

  ! 1st step
  x = t0
  call equation(x,y,ak1)
  ! 2nd step
  x = t0+a2*dt
  call equation(x,y+dt*c21*ak1,ak2)
  ! 3rd step
  x = t0+a3*dt
  call equation(x,y+dt*(c31*ak1+c32*ak2),ak3)
  ! 4th step
  x = t0+a4*dt
  call equation(x,y+dt*(c41*ak1+c42*ak2+c43*ak3),ak4)
  ! 5th step
  x = t0+dt
  call equation(x,y+dt*(c51*ak1+c52*ak2+c53*ak3+c54*ak4),ak5)
  ! 6th step
  x = t0+a6*dt
  call equation(x,y+dt*(c61*ak1+c62*ak2+c63*ak3+c64*ak4+c65*ak5),ak6)
  ! Accumulate increments with proper weights
  y(:,:) = y+dt*(c1*ak1+c3*ak3+c4*ak4+c5*ak5+c6*ak6)
  err = maxval(abs(dt*(dc1*ak1+dc3*ak3+dc4*ak4+dc5*ak5+c6*ak6)))

end subroutine rk45

subroutine rkck(t0,y,err)
  !*********************************************************************
  ! Runge-Kutta-Cash-Karp integration algorithm
  ! implementated following the Numerical Recipes in Fortran 90
  ! by W.H. Press, S.A. Teukolsky et al. (1997).
  ! Also, a clear ansatz can be found in the book
  ! Numerical computations with GPUs by V. Kindratenko
  !*********************************************************************

  real(kind=wp), intent(in) :: t0
  real(kind=wp), intent(out) :: err
  complex(kind=wp), intent(inout) :: y(d,d)
  real(kind=wp) :: x
  real(kind=wp), parameter :: a2 = 0.2_wp, &
                              a3 = 0.3_wp, &
                              a4 = 0.6_wp, &
                              a6 = 0.875_wp, &
                              b21 = 0.2_wp, &
                              b31 = 3.0_wp/40.0_wp, &
                              b32 = 9.0_wp/40.0_wp, &
                              b41 = 0.3_wp, &
                              b42 = -0.9_wp, &
                              b43 = 1.2_wp, &
                              b51 = -11.0_wp/54.0_wp, &
                              b52 = 2.5_wp, &
                              b53 = -70.0_wp/27.0_wp, &
                              b54 = 35.0_wp/27.0_wp, &
                              b61 = 1631.0_wp/55296.0_wp, &
                              b62 = 175.0_wp/512.0_wp, &
                              b63 = 575.0_wp/13824.0_wp, &
                              b64 = 44275.0_wp/110592.0_wp, &
                              b65 = 253.0_wp/4096.0_wp, &
                              c1 = 37.0_wp/378.0_wp, &
                              c3 = 250.0_wp/621.0_wp, &
                              c4 = 125.0_wp/594.0_wp, &
                              c6 = 512.0_wp/1771.0_wp, &
                              dc1 = c1-2825.0_wp/27648.0_wp, &
                              dc3 = c3-18575.0_wp/48384.0_wp, &
                              dc4 = c4-13525.0_wp/55296.0_wp, &
                              dc5 = -277.0_wp/14336.0_wp, &
                              dc6 = c6-0.25_wp
  ! 1st step
  x = t0
  call equation(x,y,ak1)
  ! 2nd step
  x = t0+a2*dt
  call equation(x,y+b21*dt*ak1,ak2)
  ! 3rd step
  x = t0+a3*dt
  call equation(x,y+dt*(b31*ak1+b32*ak2),ak3)
  ! 4th step
  x = t0+a4*dt
  call equation(x,y+dt*(b41*ak1+b42*ak2+b43*ak3),ak4)
  ! 5th step
  x = t0+dt
  call equation(x,y+dt*(b51*ak1+b52*ak2+b53*ak3+b54*ak4),ak5)
  ! 6th step
  x = t0+a6*dt
  call equation(x,y+dt*(b61*ak1+b62*ak2+b63*ak3+b64*ak4+b65*ak5),ak6)
  ! Accumulate increments with proper weights
  y(:,:) = y+dt*(c1*ak1+c3*ak3+c4*ak4+c6*ak6)
  err = maxval(abs(dt*(dc1*ak1+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)))

end subroutine rkck

subroutine rk4_sph(t0,y)
  !*********************************************************************
  ! convenient Runge-Kutta method of 4th order
  ! for integration of set of matrices in ITOs basis
  !*********************************************************************

  use Constants, only: Six, Half

  real(kind=wp), intent(in) :: t0
  complex(kind=wp), intent(inout) :: y(len_sph,d,d)

  call equation_sph(t0,y,midk1)
  call equation_sph(t0+Half*timestep,y+Half*timestep*midk1,midk2)
  call equation_sph(t0+Half*timestep,y+Half*timestep*midk2,midk3)
  call equation_sph(t0+timestep,y+timestep*midk3,midk4)
  y(:,:,:) = y+timestep/Six*(midk1+2*midk2+2*midk3+midk4)

end subroutine rk4_sph

end module integrators
