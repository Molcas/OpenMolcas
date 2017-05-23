************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine ISWAP(N,X,incX,Y,incY)
************************************************************************
*                                                                      *
*     Interchange vectors X and Y                                      *
*                                                                      *
*     calling arguments:                                               *
*     N       : Integer, input.                                        *
*               Number of input elements.                              *
*     X       : Array of Integer                                       *
*               Vector X                                               *
*     incX    : Integer, input.                                        *
*               Stride of vector X.                                    *
*     Y       : Array of Integer                                       *
*               Vector Y                                               *
*     incY    : Integer, input.                                        *
*               Stride of vector Y.                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Dimension X(*),Y(*)

      If ( N.lt.0 ) then
         Write (6,*)
         Write (6,*) '  *** Error in subroutine ISWAP ***'
         Write (6,*) '  Invalid number of elements in vectors X and Y :'
         Write (6,*) '  N must be larger than zero'
         Write (6,*)
         Call Abend
      End If

      iX=1
      If ( incX.lt.0 ) iX=1+(1-N)*incX
      iY=1
      If ( incY.lt.0 ) iY=1+(1-N)*incY

      Do i=0,N-1
         Temp=Y(iY+i*incY)
         Y(iY+i*incY)=X(iX+i*incX)
         X(iX+i*incX)=Temp
      End Do

      Return
      End
