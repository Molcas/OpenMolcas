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
! Copyright (C) 2020, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Random_Vector
!
!> @brief
!>   Generate a random vector in an \f$ N \f$ dimensional hypersphere
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Generate a random \f$ N \f$ dimensional vector of unit length or less,
!> by using random normal-distributed numbers. \cite Mul1959-CACM-2-19
!> The normal-distributed random numbers are generated with the Box--Muller
!> transform. \cite Box1958-AMS-29-610
!>
!> @param[in]  N    Dimension of the generated vector
!> @param[out] Vec  Generated random vector
!> @param[in]  UVec Specifies whether the vector should be of unit length
!***********************************************************************

subroutine Random_Vector(N,Vec,UVec)

use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(out) :: Vec(N)
logical(kind=iwp), intent(in) :: UVec
integer(kind=iwp) :: i, iSeed = 0
real(kind=wp) :: m, sm, tot_m, U, V, X, Y
real(kind=wp), parameter :: Thr = 1.0e-8_wp
real(kind=wp), external :: Random_Molcas

! Initialize random seed
if (iSeed == 0) call GetSeed(iSeed)

! To reduce numerical errors, repeat until the size is reasonable
tot_m = Zero
do while ((tot_m < Thr) .or. (tot_m > One/Thr))
  tot_m = Zero
  do i=1,N,2
    ! Get two independent normal distributed-variales, X and Y
    ! See doi:10.1214/aoms/1177706645
    U = Random_Molcas(iSeed)
    V = Two*Pi*Random_Molcas(iSeed)
    m = -Two*log(U)
    sm = sqrt(m)
    X = sm*cos(V)
    Y = sm*sin(V)
    ! Add them to the vector,
    ! being careful with the last one if N is odd
    ! See doi:10.1145/377939.377946
    if (i == N) then
      Vec(i) = X
      tot_m = tot_m+X*X
    else
      Vec(i:i+1) = [X,Y]
      tot_m = tot_m+m
    end if
  end do
end do
! Normalize the vector
! and scale by a random (0,1) number if no unit vector desired
if (UVec) then
  sm = One
else
  sm = Random_Molcas(iSeed)
end if
Vec(:) = sm/sqrt(tot_m)*Vec(:)

return

end subroutine Random_Vector
