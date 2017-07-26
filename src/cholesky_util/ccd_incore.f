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
*  CCD_InCore
*
*> @brief
*>   Complete Cholesky decomposition of the symmetric positive definite matrix \p X
*> @author Thomas Bondo Pedersen
*>
*> @details
*> The \p n &times; \p n matrix \p X is Cholesky decomposed and the resulting
*> Cholesky vectors are returned in the lower triangle of \p X.
*> A non-zero return code signals
*> that an error has occured (\p X might e.g. be non-positive
*> definite) and the output is ill-defined.
*>
*> @param[in,out] X   Matrix to be Cholesky decomposed;
*>                    on exit, lower triangle contains the vectors
*> @param[in]     n   Linear dimension of \p X
*> @param[out]    irc Return code
************************************************************************
      Subroutine CCD_InCore(X,n,irc)
      Implicit None
      Integer n
      Real*8  X(n,n)
      Integer irc

      Integer i, j, k
      Real*8  Fac

      irc = 0
      If (n .lt. 1) Return ! return (nothing to do)

      Do j=1,n
         ! Check for negative diagonal
         If (X(j,j).gt.0.0d0) Then
            Fac=1.0d0/sqrt(X(j,j))
         Else
            irc=1
            Return
         End If
         ! Compute vector j
         Do i=j,n
            X(i,j)=Fac*X(i,j)
         End Do
         ! Subtract from remaining columns
         Do k=j+1,n
            Fac=X(k,j)
            Do i=k,n
               X(i,k)=X(i,k)-X(i,j)*Fac
            End Do
         End DO
      End Do

      End
