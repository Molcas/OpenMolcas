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
      SubRoutine Statistics(X,n,Stat,ip_Mean,ip_MeanAbs,ip_Min,ip_Max,
     &                      ip_MaxAbs,ip_Variance,ip_VarianceU)
************************************************************
*
*   <DOC>
*     <Name>Statistics</Name>
*     <Syntax>Call Statistics(X,n,Stat,ip\_Mean,ip\_MeanAbs,
*                             ip\_Min,ip\_Max,
*                             ip\_MaxAbs,ip\_Variance,ip\_VarianceU)
*     </Syntax>
*     <Arguments>
*     \Argument{X}{Data array}{Real*8}{in}
*     \Argument{n}{Dimension of X}{Integer}{in}
*     \Argument{Stat}{Statistics}{Real*8}{out}
*     \Argument{ip\_Mean}{Pointer to mean value in Stat}{Integer}{in}
*     \Argument{ip\_MeanAbs}{Pointer to mean abs. value in Stat}
*              {Integer}{in}
*     \Argument{ip\_Min}{Pointer to min. value in Stat}{Integer}{in}
*     \Argument{ip\_Max}{Pointer to max. value in Stat}{Integer}{in}
*     \Argument{ip\_MaxAbs}{Pointer to max. abs. value in Stat}
*              {Integer}{in}
*     \Argument{ip\_Variance}{Pointer to biased variance in Stat}
*              {Integer}{in}
*     \Argument{ip\_VarianceU}{Pointer to unbiased variance in Stat}
*              {Integer}{in}
*     </Arguments>
*     <Purpose></Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     Computes the following statistics for a data set X(i),i=1,2,..,n:
*     Stat(ip\_Mean)     : mean value
*     Stat(ip\_MeanAbs)  : mean absolute value
*     Stat(ip\_Min)      : minimum value
*     Stat(ip\_Max)      : maximum value
*     Stat(ip\_MaxAbs)   : maximum absolute value
*     Stat(ip\_Variance) : biased variance (i.e. RMS w.r.t. mean value)
*     Stat(ip\_VarianceU): unbiased variance
*     </Description>
*    </DOC>
*
************************************************************

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
