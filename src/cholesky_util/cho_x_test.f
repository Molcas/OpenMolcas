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
      SubRoutine Cho_X_Test(X,n,Square,Vec,nVec,xf,Y,lY,Thr,irc)
************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_Test</Name>
*     <Syntax>
*       Call Cho\_X\_Test(X,n,Square,Vec,nVec,xf,Y,lY,Thr,irc)
*     </Syntax>
*     <Arguments>
*     \Argument{X}{The n-by-n matrix to be tested}{Real*8}{in}
*       \Argument{n}{Dimension of matrix X}{Integer}{in}
*       \Argument{Square}{Flag for packing of X}{Logical}{in}
*       \Argument{Vec}{Cholesky vectors representing X}{Real*8}{in}
*       \Argument{nVec}{Number of Cholesky vectors}{Integer}{in}
*       \Argument{xf}{Factor for the scaling of the vectors}{Real*8}{in}
*       \Argument{Y}{Y = X - xf*Vec*VecT}{Real*8}{out}
*       \Argument{lY}{Dimension of array Y}{Integer}{in}
*       \Argument{Thr}{Threshold allowed for RMS error}{Real*8}{in}
*       \Argument{irc}{Return code}{}{out}
*     </Arguments>
*     <Purpose>Check Cholesky decomposition of matrix X</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        The difference between the matrix X and its Cholesky
*        representation is calculated and returned in array Y. If the
*        RMS error is less than Thr, irc=0 is returned, else a positive
*        value is returned (a negative value signals an input error).
*        The input matrix X may be stored as a lower triangle
*        (Square=.false.) or as a full square matrix (Square=.true.).
*        The dimension of Y will be the same as that of X
*        (i.e. n(n+1)/2 for lower triangular storage or n*n for full
*        square storage).
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
      Dimension X(*), Vec(n,nVec), Y(lY)
      Logical   Square

      external ddot_

C     Check input.
C     ------------

      If (n .lt. 1) Then ! nothing to do
         irc = 0
         Return
      Else If (nVec.lt.0 .or. Thr.lt.0.0d0) Then ! input error
         irc = -1
         Return
      End If
      If (Square) Then
         lX = n*n
      Else
         lX = n*(n+1)/2
      End If
      If (lY .lt. lX) Then ! insufficient memory
         irc = -2
         Return
      End If

C     Compute Y(ij) = X(ij) - xf * sum_J L(iJ) * L(jJ).
C     -------------------------------------------------

      Call dCopy_(lX,X,1,Y,1)
      If (Square) Then
         Call DGEMM_('N','T',n,n,nVec,
     &              -xf,Vec,n,Vec,n,1.0d0,Y,n)
      Else
         Call dGeMM_Tri('N','T',n,n,nVec,
     &                  -xf,Vec,n,Vec,n,1.0d0,Y,n)
      End If

C     Check RMS error.
C     ----------------

      RMS = sqrt(dDot_(lX,Y,1,Y,1)/dble(lX))
      If (RMS .gt. Thr) Then
         irc = 1
      Else
         irc = 0
      End If

      End
