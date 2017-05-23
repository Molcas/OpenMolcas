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
      Subroutine CCD_InCore(X,n,irc)
************************************************************
*
*   <DOC>
*     <Name>CCD\_InCore</Name>
*     <Syntax>Call CCD\_InCore(X,n,irc)</Syntax>
*     <Arguments>
*       \Argument{X}{Matrix to be Cholesky decomposed, dimension
*                    X(n,n); on exit, lower triangle contains
*                    the vectors}{Real*8}{inout}
*       \Argument{n}{Linear dimension of X}{Integer}{in}
*       \Argument{irc}{Return code}{Integer}{out}
*     </Arguments>
*     <Purpose>Complete Cholesky decomposition of the symmetric positive
*              definite matrix X</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        The n-by-n matrix X is Cholesky decomposed and the resulting
*        Cholesky vectors are returned in the lower triangle of X.
*        A non-zero return code signals
*        that an error has occured (X might f.ex. be non-positive
*        definite) and the output is ill-defined.
*     </Description>
*    </DOC>
*
************************************************************
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
