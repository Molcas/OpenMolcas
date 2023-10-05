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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChkAcc(K,InitR,Error,R,Change)
!-----------------------------------------------------------------------
! Function : Check the convergence
!-----------------------------------------------------------------------

use ReMez_mod, only: IW
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: K
integer(kind=iwp), intent(inout) :: InitR
real(kind=wp), intent(in) :: Error, R
logical(kind=iwp), intent(inout) :: Change
real(kind=wp) :: ErrMax, ErrMin
logical(kind=iwp) :: Skip
real(kind=wp), parameter :: E2(7) = [2.128e-2_wp,2.080e-4_wp,1.833e-6_wp,1.542e-8_wp,1.261e-10_wp,1.012e-12_wp,8.020e-15_wp], &
                            E5(7) = [7.075e-2_wp,3.437e-3_wp,1.500e-4_wp,6.258e-6_wp,2.543e-7_wp,1.016e-8_wp,4.007e-10_wp], &
                            E10(16) = [8.556e-2_wp,8.752e-3_wp,7.145e-4_wp,5.577e-5_wp,4.243e-6_wp,3.173e-7_wp,2.344e-8_wp, &
                                       1.716e-9_wp,1.248e-10_wp,9.021e-12_wp,6.492e-13_wp,4.654e-14_wp,3.326e-15_wp,2.371e-16_wp, &
                                       1.708e-17_wp,5.136e-18_wp], &
                            E20(17) = [Zero,1.448e-2_wp,1.819e-3_wp,2.169e-4_wp,2.521e-5_wp,2.880e-6_wp,3.252e-7_wp,3.640e-8_wp, &
                                       4.045e-9_wp,4.470e-10_wp,4.918e-11_wp,5.389e-12_wp,5.888e-13_wp,6.414e-14_wp,6.972e-15_wp, &
                                       7.562e-16_wp,8.202e-17_wp], &
                            E30(17) = [Zero,1.699e-2_wp,2.627e-3_wp,3.795e-4_wp,5.336e-5_wp,7.379e-6_wp,1.008e-6_wp,1.365e-7_wp, &
                                       1.836e-8_wp,2.456e-9_wp,3.270e-10_wp,4.337e-11_wp,5.735e-12_wp,7.562e-13_wp,9.947e-14_wp, &
                                       1.306e-14_wp,1.711e-15_wp], &
                            E40(18) = [Zero,1.784e-2_wp,3.215e-3_wp,5.230e-4_wp,8.266e-5_wp,1.285e-5_wp,1.973e-6_wp,3.003e-7_wp, &
                                       4.539e-8_wp,6.822e-9_wp,1.021e-9_wp,1.522e-10_wp,2.262e-11_wp,3.352e-12_wp,4.956e-13_wp, &
                                       7.312e-14_wp,1.077e-14_wp,1.584e-15_wp], &
                            E50(19) = [Zero,1.785e-2_wp,3.659e-3_wp,6.469e-4_wp,1.110e-4_wp,1.872e-5_wp,3.121e-6_wp,5.155e-7_wp, &
                                       8.456e-8_wp,1.380e-8_wp,2.241e-9_wp,3.625e-10_wp,5.847e-11_wp,9.405e-12_wp,1.509e-12_wp, &
                                       2.417e-13_wp,3.863e-14_wp,6.165e-15_wp,9.826e-16_wp], &
                            E60(20) = [Zero,Zero,4.001e-3_wp,7.541e-4_wp,1.377e-4_wp,2.471e-5_wp,4.382e-6_wp,7.701e-7_wp, &
                                       1.344e-7_wp,2.333e-8_wp,4.031e-9_wp,6.939e-10_wp,1.191e-10_wp,2.038e-11_wp,3.479e-12_wp, &
                                       5.927e-13_wp,1.008e-13_wp,1.711e-14_wp,2.902e-15_wp,4.914e-16_wp], &
                            E70(20) = [Zero,Zero,4.271e-3_wp,8.474e-4_wp,1.626e-4_wp,3.066e-5_wp,5.711e-6_wp,1.054e-6_wp, &
                                       1.933e-7_wp,3.525e-8_wp,6.398e-9_wp,1.157e-9_wp,2.085e-10_wp,3.749e-11_wp,6.723e-12_wp, &
                                       1.203e-12_wp,2.150e-13_wp,3.834e-14_wp,6.829e-15_wp,1.215e-15_wp], &
                            E80(20) = [Zero,Zero,4.485e-3_wp,9.293e-4_wp,1.858e-4_wp,3.648e-5_wp,7.077e-6_wp,1.361e-6_wp, &
                                       2.598e-7_wp,4.932e-8_wp,9.323e-9_wp,1.756e-9_wp,3.295e-10_wp,6.169e-11_wp,1.152e-11_wp, &
                                       2.147e-12_wp,4.004e-13_wp,7.420e-14_wp,1.376e-14_wp,2.549e-15_wp], &
                            E90(20) = [Zero,Zero,4.655e-3_wp,1.002e-3_wp,2.074e-4_wp,4.213e-5_wp,8.458e-6_wp,1.683e-6_wp, &
                                       3.325e-7_wp,6.532e-8_wp,1.278e-8_wp,2.490e-9_wp,4.836e-10_wp,9.369e-11_wp,1.811e-11_wp, &
                                       3.492e-12_wp,6.723e-13_wp,1.292e-13_wp,2.480e-14_wp,4.754e-15_wp], &
                            E100(20) = [Zero,Zero,4.789e-3_wp,1.066e-3_wp,2.274e-4_wp,4.760e-5_wp,9.841e-6_wp,2.016e-6_wp, &
                                        4.103e-7_wp,8.303e-8_wp,1.673e-8_wp,3.357e-9_wp,6.716e-10_wp,1.340e-10_wp,2.667e-11_wp, &
                                        5.298e-12_wp,1.050e-12_wp,2.079e-13_wp,4.110e-14_wp,8.114e-15_wp], &
                            E200(20) = [Zero,Zero,5.052e-3_wp,1.456e-3_wp,3.707e-4_wp,9.217e-5_wp,2.261e-5_wp,5.498e-6_wp, &
                                        1.327e-6_wp,3.186e-7_wp,7.613e-8_wp,1.812e-8_wp,4.301e-9_wp,1.018e-9_wp,2.403e-10_wp, &
                                        5.663e-11_wp,1.332e-11_wp,3.128e-12_wp,7.333e-13_wp,1.717e-13_wp], &
                            E300(20) = [Zero,Zero,Zero,1.628e-3_wp,4.554e-4_wp,1.235e-4_wp,3.297e-5_wp,8.724e-6_wp,2.292e-6_wp, &
                                        5.986e-7_wp,1.557e-7_wp,4.032e-8_wp,1.041e-8_wp,2.681e-9_wp,6.888e-10_wp,1.766e-10_wp, &
                                        4.519e-11_wp,1.155e-11_wp,2.946e-12_wp,7.507e-13_wp], &
                            E400(20) = [Zero,Zero,Zero,1.695e-3_wp,5.117e-4_wp,1.467e-4_wp,4.141e-5_wp,1.157e-5_wp,3.209e-6_wp, &
                                        8.850e-7_wp,2.430e-7_wp,6.645e-8_wp,1.812e-8_wp,4.926e-9_wp,1.336e-9_wp,3.616e-10_wp, &
                                        9.771e-11_wp,2.636e-11_wp,7.099e-12_wp,1.910e-12_wp], &
                            E500(20) = [Zero,Zero,Zero,1.700e-3_wp,5.517e-4_wp,1.649e-4_wp,4.842e-5_wp,1.407e-5_wp,4.061e-6_wp, &
                                        1.165e-6_wp,3.326e-7_wp,9.463e-8_wp,2.683e-8_wp,7.587e-9_wp,2.140e-9_wp,6.025e-10_wp, &
                                        1.693e-10_wp,4.750e-11_wp,1.331e-11_wp,3.723e-12_wp], &
                            E600(20) = [Zero,Zero,Zero,Zero,5.811e-4_wp,1.795e-4_wp,5.438e-5_wp,1.630e-5_wp,4.848e-6_wp, &
                                        1.434e-6_wp,4.221e-7_wp,1.238e-7_wp,3.618e-8_wp,1.055e-8_wp,3.067e-9_wp,8.902e-10_wp, &
                                        2.579e-10_wp,7.459e-11_wp,2.154e-11_wp,6.213e-12_wp], &
                            E700(20) = [Zero,Zero,Zero,Zero,6.031e-4_wp,1.915e-4_wp,5.952e-5_wp,1.829e-5_wp,5.577e-6_wp, &
                                        1.691e-6_wp,5.101e-7_wp,1.533e-7_wp,4.594e-8_wp,1.372e-8_wp,4.091e-9_wp,1.217e-9_wp, &
                                        3.613e-10_wp,1.071e-10_wp,3.170e-11_wp,9.371e-12_wp], &
                            E800(20) = [Zero,Zero,Zero,Zero,6.193e-4_wp,2.016e-4_wp,6.401e-5_wp,2.008e-5_wp,6.252e-6_wp, &
                                        1.935e-6_wp,5.960e-7_wp,1.829e-7_wp,5.593e-8_wp,1.706e-8_wp,5.190e-9_wp,1.576e-9_wp, &
                                        4.776e-10_wp,1.445e-10_wp,4.367e-11_wp,1.318e-11_wp], &
                            E900(20) = [Zero,Zero,Zero,Zero,6.309e-4_wp,2.102e-4_wp,6.799e-5_wp,2.172e-5_wp,6.881e-6_wp, &
                                        2.167e-6_wp,6.795e-7_wp,2.122e-7_wp,6.605e-8_wp,2.050e-8_wp,6.349e-9_wp,1.962e-9_wp, &
                                        6.052e-10_wp,1.864e-10_wp,5.732e-11_wp,1.760e-11_wp], &
                            E1000(20) = [Zero,Zero,Zero,Zero,6.385e-4_wp,2.177e-4_wp,7.153e-5_wp,2.321e-5_wp,7.468e-6_wp, &
                                         2.389e-6_wp,7.605e-7_wp,2.412e-7_wp,7.623e-8_wp,2.403e-8_wp,7.555e-9_wp,2.371e-9_wp, &
                                         7.426e-10_wp,2.322e-10_wp,7.251e-11_wp,2.261e-11_wp], &
                            E2000(20) = [Zero,Zero,Zero,Zero,6.428e-4_wp,2.570e-4_wp,9.365e-5_wp,3.341e-5_wp,1.180e-5_wp, &
                                         4.140e-6_wp,1.445e-6_wp,5.025e-7_wp,1.741e-7_wp,6.017e-8_wp,2.074e-8_wp,7.134e-9_wp, &
                                         2.450e-9_wp,8.396e-10_wp,2.874e-10_wp,9.824e-11_wp], &
                            E3000(20) = [Zero,Zero,Zero,Zero,Zero,2.646e-4_wp,1.047e-4_wp,3.926e-5_wp,1.456e-5_wp,5.360e-6_wp, &
                                         1.963e-6_wp,7.157e-7_wp,2.601e-7_wp,9.426e-8_wp,3.407e-8_wp,1.229e-8_wp,4.425e-9_wp, &
                                         1.590e-9_wp,5.708e-10_wp,2.046e-10_wp], &
                            E4000(20) = [Zero,Zero,Zero,Zero,Zero,Zero,1.110e-4_wp,4.316e-5_wp,1.654e-5_wp,6.284e-6_wp, &
                                         2.375e-6_wp,8.936e-7_wp,3.351e-7_wp,1.253e-7_wp,4.673e-8_wp,1.739e-8_wp,6.461e-9_wp, &
                                         2.396e-9_wp,8.873e-10_wp,3.281e-10_wp], &
                            E5000(20) = [Zero,Zero,Zero,Zero,Zero,Zero,1.146e-4_wp,4.596e-5_wp,1.804e-5_wp,7.020e-6_wp, &
                                         2.715e-6_wp,1.046e-6_wp,4.013e-7_wp,1.536e-7_wp,5.861e-8_wp,2.232e-8_wp,8.484e-9_wp, &
                                         3.220e-9_wp,1.220e-9_wp,4.617e-10_wp], &
                            E6000(20) = [Zero,Zero,Zero,Zero,Zero,Zero,1.162e-4_wp,4.807e-5_wp,1.923e-5_wp,7.626e-6_wp, &
                                         3.004e-6_wp,1.178e-6_wp,4.604e-7_wp,1.794e-7_wp,6.972e-8_wp,2.704e-8_wp,1.046e-8_wp, &
                                         4.043e-9_wp,1.560e-9_wp,6.011e-10_wp], &
                            E7000(20) = [Zero,Zero,Zero,Zero,Zero,Zero,1.163e-4_wp,4.969e-5_wp,2.021e-5_wp,8.137e-6_wp, &
                                         3.254e-6_wp,1.295e-6_wp,5.138e-7_wp,2.032e-7_wp,8.014e-8_wp,3.154e-8_wp,1.239e-8_wp, &
                                         4.858e-9_wp,1.902e-9_wp,7.438e-10_wp], &
                            E8000(20) = [Zero,Zero,Zero,Zero,Zero,Zero,Zero,5.096e-5_wp,2.103e-5_wp,8.576e-6_wp,3.474e-6_wp, &
                                         1.400e-6_wp,5.624e-7_wp,2.252e-7_wp,8.993e-8_wp,3.584e-8_wp,1.425e-8_wp,5.659e-9_wp, &
                                         2.244e-9_wp,8.883e-10_wp], &
                            E9000(20) = [Zero,Zero,Zero,Zero,Zero,Zero,Zero,5.195e-5_wp,2.172e-5_wp,8.958e-6_wp,3.669e-6_wp, &
                                         1.495e-6_wp,6.070e-7_wp,2.457e-7_wp,9.917e-8_wp,3.995e-8_wp,1.606e-8_wp,6.445e-9_wp, &
                                         2.582e-9_wp,1.033e-9_wp], &
                            E10000(20) = [Zero,Zero,Zero,Zero,Zero,Zero,Zero,5.271e-5_wp,2.232e-5_wp,9.296e-6_wp,3.844e-6_wp, &
                                          1.582e-6_wp,6.481e-7_wp,2.648e-7_wp,1.079e-7_wp,4.388e-8_wp,1.780e-8_wp,7.213e-9_wp, &
                                          2.918e-9_wp,1.179e-9_wp], &
                            ERMax(20) = [Zero,Zero,Zero,Zero,Zero,Zero,Zero,5.392e-5_wp,2.611e-5_wp,1.312e-5_wp,6.807e-6_wp, &
                                         3.630e-6_wp,1.984e-6_wp,1.108e-6_wp,6.311e-7_wp,3.659e-7_wp,2.155e-7_wp,1.289e-7_wp, &
                                         7.811e-8_wp,4.794e-8_wp], &
                            RList(30) = [2.0e0_wp,5.0e0_wp,1.0e1_wp,2.0e1_wp,3.0e1_wp,4.0e1_wp,5.0e1_wp,6.0e1_wp,7.0e1_wp, &
                                         8.0e1_wp,9.0e1_wp,1.0e2_wp,2.0e2_wp,3.0e2_wp,4.0e2_wp,5.0e2_wp,6.0e2_wp,7.0e2_wp, &
                                         8.0e2_wp,9.0e2_wp,1.0e3_wp,2.0e3_wp,3.0e3_wp,4.0e3_wp,5.0e3_wp,6.0e3_wp,7.0e3_wp, &
                                         8.0e3_wp,9.0e3_wp,1.0e4_wp]

Skip = .false.
select case (InitR)
  case default ! (1)
    ! probably default should be 31...
    ! this is what you get when you use computed gotos
    ErrMin = E2(K)
    ErrMax = E5(K)
  case (2)
    ErrMin = E5(K)
    ErrMax = E10(K)
  case (3)
    ErrMin = E10(K)
    ErrMax = E20(K)
  case (4)
    ErrMin = E20(K)
    ErrMax = E30(K)
  case (5)
    ErrMin = E30(K)
    ErrMax = E40(K)
  case (6)
    ErrMin = E40(K)
    ErrMax = E50(K)
  case (7)
    ErrMin = E50(K)
    ErrMax = E60(K)
  case (8)
    ErrMin = E60(K)
    ErrMax = E70(K)
  case (9)
    ErrMin = E70(K)
    ErrMax = E80(K)
  case (10)
    ErrMin = E80(K)
    ErrMax = E90(K)
  case (11)
    ErrMin = E90(K)
    ErrMax = E100(K)
  case (12)
    ErrMin = E100(K)
    ErrMax = E200(K)
  case (13)
    ErrMin = E200(K)
    ErrMax = E300(K)
  case (14)
    ErrMin = E300(K)
    ErrMax = E400(K)
  case (15)
    ErrMin = E400(K)
    ErrMax = E500(K)
  case (16)
    ErrMin = E500(K)
    ErrMax = E600(K)
  case (17)
    ErrMin = E600(K)
    ErrMax = E700(K)
  case (18)
    ErrMin = E700(K)
    ErrMax = E800(K)
  case (19)
    ErrMin = E800(K)
    ErrMax = E900(K)
  case (20)
    ErrMin = E900(K)
    ErrMax = E1000(K)
  case (21)
    ErrMin = E1000(K)
    ErrMax = E2000(K)
  case (22)
    ErrMin = E2000(K)
    ErrMax = E3000(K)
  case (23)
    ErrMin = E3000(K)
    ErrMax = E4000(K)
  case (24)
    ErrMin = E4000(K)
    ErrMax = E5000(K)
  case (25)
    ErrMin = E5000(K)
    ErrMax = E6000(K)
  case (26)
    ErrMin = E6000(K)
    ErrMax = E7000(K)
  case (27)
    ErrMin = E7000(K)
    ErrMax = E8000(K)
  case (28)
    ErrMin = E8000(K)
    ErrMax = E9000(K)
  case (29)
    ErrMin = E9000(K)
    ErrMax = E10000(K)
  case (30)
    ErrMin = E10000(K)
    ErrMax = ERMax(K)
  case (31)
    Skip = .true.
end select

if (.not. Skip) then
  write(IW,'(/A/)') ' Check the accuracy of the convergence'
  write(IW,'(A,F10.3,2X,A,E18.9E2)') ' R =',RList(InitR),'Maximum error = ',ErrMin
  write(IW,'(A,F10.3,2X,A,E18.9E2)') ' R =',R,'Maximum error = ',Error
  write(IW,'(A,F10.3,2X,A,E18.9E2/)') ' R =',RList(InitR+1),'Maximum error = ',ErrMax
  if ((Error > ErrMin) .and. (Error < ErrMax)) then
    write(IW,'(A)') ' Convergence is GOOD.'
    Change = .false.
  else
    write(IW,'(A)') ' Convergence is not good.'
    InitR = InitR+1
    Change = .true.
  end if
end if

return

end subroutine ChkAcc
