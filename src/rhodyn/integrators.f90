subroutine classic_rk4(t0,y)
  use rhodyn_data
  implicit none
  real(8) :: t0
  complex(8),dimension(:,:) :: y
  procedure(equation_func)  :: equation
  call equation(t0,y,ak1)
  call equation(t0+0.5*timestep,y+0.5*timestep*ak1,ak2)
  call equation(t0+0.5*timestep,y+0.5*timestep*ak2,ak3)
  call equation(t0+timestep,y+timestep*ak3,ak4)
  y = y + timestep/6.0*(ak1+2*ak2+2*ak3+ak4)
end

subroutine rk4(t0,y)
  use rhodyn_data
  implicit none
  real(8) :: t0, x
  complex(8),dimension(:,:) :: y
  procedure(equation_func)  :: equation
  real(8),parameter :: a2  = 0.25,&
                       a3  = 0.375d0,&
                       a4  = 0.92307692307692313d0,&
                       c21 = 0.25d0,&
                       c31 = 0.09375d0,&
                       c32 = 0.28125d0,&
                       c41 = 0.87938097405553028d0,&
                       c42 =-3.2771961766044608d0,&
                       c43 = 3.3208921256258535d0,&
                       c51 = 2.0324074074074074d0,&
                       c52 =-8,&
                       c53 = 7.1734892787524362d0,&
                       c54 =-0.20589668615984405d0,&
                       c1  = 0.11574074074074074d0,&
                       c3  = 0.54892787524366471d0,&
                       c4  = 0.53533138401559455d0,&
                       c5  =-0.2d0
  x=t0
  call equation(x,y,ak1)
  x=t0+a2*timestep
  call equation(x,y+ak1*timestep*c21,ak2)
  x=t0+a3*timestep
  call equation(x,y+timestep*(c31*ak1+c32*ak2),ak3)
  x=t0+a4*timestep
  call equation(x,y+timestep*(c41*ak1+c42*ak2+c43*ak3),ak4)
  x=t0+timestep
  call equation(x,y+timestep*(c51*ak1+c52*ak2+c53*ak3+c54*ak4),ak5)
  y = y + (c1*ak1+c3*ak3+c4*ak4+c5*ak5)*timestep
end

subroutine rk5(t0,y)
  use rhodyn_data
  implicit none
  real(8) :: t0, x
  complex(8),dimension(:,:) :: y
  procedure(equation_func)  :: equation
  real(8),parameter :: a2  = 0.25d0,&
                      a3  = 0.375d0,&
                      a4  = 0.92307692307692313d0,&
                      c21 = 0.25d0,&
                      c31 = 0.09375d0,&
                      c32 = 0.28125d0,&
                      c41 = 0.87938097405553028,&
                      c42 =-3.2771961766044608,&
                      c43 = 3.3208921256258535,&
                      c51 = 2.0324074074074074,&
                      c52 =-8.d0,&
                      c53 = 7.1734892787524362,&
                      c54 =-0.20589668615984405,&
                      c61 =-0.29629629629629628,&
                      c62 = 2d0,&
                      c63 =-1.3816764132553607,&
                      c64 = 0.45297270955165692,&
                      c65 =-0.275d0,&
                      c1 = 0.11851851851851852,&
                      c3 = 0.51898635477582844,&
                      c4 = 0.50613149034201665,&
                      c5 =-0.18d0,&
                      c6 = 0.036363636363636362
  x=t0
  call equation(x,y,ak1)
  x=t0+a2*timestep
  call equation(x,y+timestep*c21*ak1,ak2)
  x=t0+a3*timestep
  call equation(x,y+timestep*(c31*ak1+c32*ak2),ak3)
  x=t0+a4*timestep
  call equation(x,y+timestep*(c41*ak1+c42*ak2+c43*ak3),ak4)
  x=t0+timestep
  call equation(x,y+timestep*(c51*ak1+c52*ak2+c53*ak3+c54*ak4),ak5)
  x=t0+0.5d0*timestep
  call equation(x,y+timestep*(c61*ak1+c62*ak2+c63*ak3+&
                              c64*ak4+c65*ak5),ak6)
  y = y + (c1*ak1+c3*ak3+c4*ak4+c5*ak5+c6*ak6)*timestep
end

subroutine rk45(t0,y,err)
  use rhodyn_data
  implicit none
!
! Runge-Kutta-Fehlberg 4(5) integration algorithm
!
  real(8) :: t0,err,x
  complex(8),dimension(:,:) :: y
  procedure(equation_func)  :: equation
  real(8), parameter ::a2  = 0.25,&
                       a3  = 3.0/8.0,&
                       a4  = 12.0/13.0,&
                       a6  = 0.5,&
                       c21 = 0.25,&
                       c31 = 3.0/32.0,&
                       c32 = 9.0/32.0,&
                       c41 = 1932.0/2197.0,&
                       c42 =-7200.0/2197.0,&
                       c43 = 7296.0/2197.0,&
                       c51 = 439.0/216.0,&
                       c52 =-8.0,&
                       c53 = 3680.0/513.0,&
                       c54 =-845.0/4104.0,&
                       c61 =-8.0/27.0,&
                       c62 = 2.0,&
                       c63 =-3544.0/2565.0,&
                       c64 = 1859.0/4104.0,&
                       c65 =-0.275,&
                       c1  = 16.0/135.0,&
                       c3  = 6656.0/12825.0,&
                       c4  = 28561.0/56430.0,&
                       c5  =-0.18,&
                       c6  = 2.0/55.0,&
                       dc1 = c1-25.0/216.0,&
                       dc3 = c3-1408.0/2565.0,&
                       dc4 = c4-2197.0/4104.0,&
                       dc5 = c5+0.2
!     1st step
  x=t0
  call equation(x,y,ak1)
!     2nd step
  x=t0+a2*dt
  call equation(x,y+dt*c21*ak1,ak2)
!     3rd step
  x=t0+a3*dt
  call equation(x,y+dt*(c31*ak1+c32*ak2),ak3)
!     4th step
  x=t0+a4*dt
  call equation(x,y+dt*(c41*ak1+c42*ak2+c43*ak3),ak4)
!     5th step
  x=t0+dt
  call equation(x,y+dt*(c51*ak1+c52*ak2+c53*ak3+c54*ak4),ak5)
!     6th step
  x=t0+a6*dt
  call equation(x,y+dt*(c61*ak1+c62*ak2+c63*ak3+&
                                c64*ak4+c65*ak5),ak6)
!     Accumulate increments with proper weights
  y = y + dt*(c1*ak1+c3*ak3+c4*ak4+c5*ak5+c6*ak6)
  err = maxval(abs(dt*(dc1*ak1+dc3*ak3+dc4*ak4+dc5*ak5+c6*ak6)))
end

subroutine rkck(t0,y,err)
  use rhodyn_data
  implicit none
!
!     Runge-Kutta-Cash-Karp integration algorithm
!     implementated following the Numerical Recipes in Fortran 90
!     by W.H. Press, S.A. Teukolsky et al. (1997).
!     Also, a clear ansatz can be found in the book
!     Numerical computations with GPUs by V. Kindratenko
!
  real(8) :: err, t0, x
  complex(8),dimension(:,:) :: y
  procedure(equation_func)  :: equation
  real(8), parameter :: a2  = 0.2,&
                        a3  = 0.3,&
                        a4  = 0.6,&
                        a6  = 0.875,&
                        b21 = 0.2,&
                        b31 = 3.0/40.0,&
                        b32 = 9.0/40.0,&
                        b41 = 0.3,&
                        b42 =-0.9,&
                        b43 = 1.2,&
                        b51 =-11.0/54.0,&
                        b52 = 2.5,&
                        b53 =-70.0/27.0,&
                        b54 = 35.0/27.0,&
                        b61 = 1631.0/55296.0,&
                        b62 = 175.0/512.0,&
                        b63 = 575.0/13824.0,&
                        b64 = 44275.0/110592.0,&
                        b65 = 253.0/4096.0,&
                        c1  = 37.0/378.0,&
                        c3  = 250.0/621.0,&
                        c4  = 125.0/594.0,&
                        c6  = 512.0/1771.0,&
                        dc1 = c1-2825.0/27648.0,&
                        dc3 = c3-18575.0/48384.0,&
                        dc4 = c4-13525.0/55296.0,&
                        dc5 =-277.0/14336.0,&
                        dc6 = c6-0.25
! 1st step
  x=t0
  call equation(x,y,ak1)
! 2nd step
  x=t0+a2*dt
  call equation(x,y+b21*dt*ak1,ak2)
! 3rd step
  x=t0+a3*dt
  call equation(x,y+dt*(b31*ak1+b32*ak2),ak3)
! 4th step
  x=t0+a4*dt
  call equation(x,y+dt*(b41*ak1+b42*ak2+b43*ak3),ak4)
! 5th step
  x=t0+dt
  call equation(x,y+dt*(b51*ak1+b52*ak2+b53*ak3+b54*ak4),ak5)
! 6th step
  x=t0+a6*dt
  call equation(x,y+dt*(b61*ak1+b62*ak2+b63*ak3+&
                                b64*ak4+b65*ak5),ak6)
! Accumulate increments with proper weights
  y = y + dt*(c1*ak1+c3*ak3+c4*ak4+c6*ak6)
  err = maxval(abs(dt*(dc1*ak1+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)))
end
