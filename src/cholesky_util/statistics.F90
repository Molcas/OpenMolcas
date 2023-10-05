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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  Statistics
!
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Computes the following statistics for a data set \p X(i), \p i = ``1``, ``2``, ..., ``n``:
!>
!> - \p Stat(ip_Mean):      mean value
!> - \p Stat(ip_MeanAbs):   mean absolute value
!> - \p Stat(ip_Min):       minimum value
!> - \p Stat(ip_Max):       maximum value
!> - \p Stat(ip_MaxAbs):    maximum absolute value
!> - \p Stat(ip_Variance):  biased variance (i.e. RMS w.r.t. mean value)
!> - \p Stat(ip_VarianceU): unbiased variance
!>
!> @param[in]  X            Data array
!> @param[in]  n            Dimension of \p X
!> @param[out] Stat         Statistics
!> @param[in]  ip_Mean      Pointer to mean value in \p Stat
!> @param[in]  ip_MeanAbs   Pointer to mean abs. value in \p Stat
!> @param[in]  ip_Min       Pointer to min. value in \p Stat
!> @param[in]  ip_Max       Pointer to max. value in \p Stat
!> @param[in]  ip_MaxAbs    Pointer to max. abs. value in \p Stat
!> @param[in]  ip_Variance  Pointer to biased variance in \p Stat
!> @param[in]  ip_VarianceU Pointer to unbiased variance in \p Stat
!***********************************************************************

subroutine Statistics(X,n,Stat,ip_Mean,ip_MeanAbs,ip_Min,ip_Max,ip_MaxAbs,ip_Variance,ip_VarianceU)

use Constants, only: One
use Definitions, only: iwp, wp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: n, ip_Mean, ip_MeanAbs, ip_Min, ip_Max, ip_MaxAbs, ip_Variance, ip_VarianceU
real(kind=wp), intent(in) :: X(n)
real(kind=wp), intent(_OUT_) :: Stat(*)
integer(kind=iwp) :: i
real(kind=wp) :: xMax, xMean, xMeanAbs, xMin, xn, xn1, xVariance

if (n < 1) return

xn = One/real(n,kind=wp)
if (n == 1) then
  xn1 = 9.99e15_wp
else
  xn1 = One/real(n-1,kind=wp)
end if

xMean = sum(X(1:n))
xMeanAbs = sum(abs(X(1:n)))
xMax = X(1)
xMin = X(1)
do i=2,n
  xMax = max(xMax,X(i))
  xMin = min(xMin,X(i))
end do
xMean = xMean*xn

if (ip_Mean > 0) Stat(ip_Mean) = xMean
if (ip_MeanAbs > 0) Stat(ip_MeanAbs) = xMeanAbs*xn
if (ip_Min > 0) Stat(ip_Min) = xMin
if (ip_Max > 0) Stat(ip_Max) = xMax
if (ip_MaxAbs > 0) Stat(ip_MaxAbs) = max(abs(xMax),abs(xMin))

if ((ip_Variance > 0) .or. (ip_VarianceU > 0)) then
  xVariance = sum((X(1:n)-xMean)**2)
  if (ip_VarianceU > 0) Stat(ip_VarianceU) = sqrt(xVariance*xn1)
  if (ip_Variance > 0) Stat(ip_Variance) = sqrt(xVariance*xn)
end if

end subroutine Statistics
