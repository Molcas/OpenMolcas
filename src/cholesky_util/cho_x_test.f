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
*  Cho_X_Test
*
*> @brief
*>   Check Cholesky decomposition of matrix \p X
*> @author Thomas Bondo Pedersen
*>
*> @details
*> The difference between the matrix \p X and its Cholesky
*> representation is calculated and returned in array \p Y. If the
*> RMS error is less than \p Thr, \p irc = ``0`` is returned, else a positive
*> value is returned (a negative value signals an input error).
*> The input matrix \p X may be stored as a lower triangle
*> (\p Square = ``.false.``) or as a full square matrix (\p Square = ``.true.``).
*> The dimension of \p Y will be the same as that of \p X
*> (i.e. ``n(n+1)/2`` for lower triangular storage or ``n*n`` for full
*> square storage).
*>
*> @param[in]  X      The \p n &times; \p n matrix to be tested
*> @param[in]  n      Dimension of matrix \p X
*> @param[in]  Square Flag for packing of \p X
*> @param[in]  Vec    Cholesky vectors representing \p X
*> @param[in]  nVec   Number of Cholesky vectors
*> @param[in]  xf     Factor for the scaling of the vectors
*> @param[out] Y      ``Y = X - xf*Vec*VecT``
*> @param[in]  lY     Dimension of array \p Y
*> @param[in]  Thr    Threshold allowed for RMS error
*> @param[out] irc    Return code
************************************************************************
      SubRoutine Cho_X_Test(X,n,Square,Vec,nVec,xf,Y,lY,Thr,irc)
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
