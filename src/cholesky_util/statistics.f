************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  Statistics
*
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Computes the following statistics for a data set \p X(i), \p i = ``1``, ``2``, ..., ``n``:
*>
*> - \p Stat(ip_Mean):      mean value
*> - \p Stat(ip_MeanAbs):   mean absolute value
*> - \p Stat(ip_Min):       minimum value
*> - \p Stat(ip_Max):       maximum value
*> - \p Stat(ip_MaxAbs):    maximum absolute value
*> - \p Stat(ip_Variance):  biased variance (i.e. RMS w.r.t. mean value)
*> - \p Stat(ip_VarianceU): unbiased variance
*>
*> @param[in]  X            Data array
*> @param[in]  n            Dimension of \p X
*> @param[out] Stat         Statistics
*> @param[in]  ip_Mean      Pointer to mean value in \p Stat
*> @param[in]  ip_MeanAbs   Pointer to mean abs. value in \p Stat
*> @param[in]  ip_Min       Pointer to min. value in \p Stat
*> @param[in]  ip_Max       Pointer to max. value in \p Stat
*> @param[in]  ip_MaxAbs    Pointer to max. abs. value in \p Stat
*> @param[in]  ip_Variance  Pointer to biased variance in \p Stat
*> @param[in]  ip_VarianceU Pointer to unbiased variance in \p Stat
************************************************************************
      SubRoutine Statistics(X,n,Stat,ip_Mean,ip_MeanAbs,ip_Min,ip_Max,
     &                      ip_MaxAbs,ip_Variance,ip_VarianceU)
      Implicit None
      Integer n
      Real*8  X(n), Stat(*)
      Integer ip_Mean, ip_MeanAbs, ip_Min, ip_Max, ip_MaxAbs
      Integer ip_Variance, ip_VarianceU

      Integer i
      Real*8  xn, xn1
      Real*8  xMean, xMeanAbs, xMin, xMax
      Real*8  xVariance

      If (n .lt. 1) Return

      xn = 1.0d0/dble(n)
      If (n .eq. 1) Then
         xn1 = 9.99D15
      Else
         xn1 = 1.0d0/dble(n-1)
      End If

      xMean    = X(1)
      xMeanAbs = abs(X(1))
      xMax     = X(1)
      xMin     = X(1)
      Do i = 2,n
         xMean    = xMean + X(i)
         xMeanAbs = xMeanAbs + abs(X(i))
         xMax     = max(xMax,X(i))
         xMin     = min(xMin,X(i))
      End Do
      xMean = xMean*xn

      If (ip_Mean    .gt. 0) Stat(ip_Mean)    = xMean
      If (ip_MeanAbs .gt. 0) Stat(ip_MeanAbs) = xMeanAbs*xn
      If (ip_Min     .gt. 0) Stat(ip_Min)     = xMin
      If (ip_Max     .gt. 0) Stat(ip_Max)     = xMax
      If (ip_MaxAbs  .gt. 0) Stat(ip_MaxAbs)  = max(abs(xMax),abs(xMin))

      If (ip_Variance.gt.0 .or. ip_VarianceU.gt.0) Then
         xVariance = (X(1)-xMean)**2
         Do i = 2,n
            xVariance = xVariance + (X(i)-xMean)**2
         End Do
         If (ip_VarianceU .gt. 0) Then
            Stat(ip_VarianceU) = sqrt(xVariance*xn1)
         End If
         If (ip_Variance .gt. 0) Then
            Stat(ip_Variance) = sqrt(xVariance*xn)
         End If
      End If

      End
