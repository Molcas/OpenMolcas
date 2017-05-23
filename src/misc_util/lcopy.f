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
* Copyright (C) 1992, Markus P. Fuelscher                              *
************************************************************************
      Subroutine lCopy(N,X,incX,Y,incY)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Copy vectors of integers                                         *
*                                                                      *
*     calling arguments:                                               *
*     N      : Number of elements                                      *
*     X      : input vector                                            *
*     incX   : stride for vector X                                     *
*     X      : output vector                                           *
*     incY   : stride for vector Y                                     *
*     LuCom  : FORTRAN unit number                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher                                                  *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Integer (A-Z)
      Logical X,Y
      Dimension X(1+incX*(N-1)),Y(1+incY*(N-1))
*
*----------------------------------------------------------------------*
*     Start procedure                                                  *
*----------------------------------------------------------------------*
      If ( N.eq.0 ) Return
      If ( N.lt.0 ) Then
         Write (6,*) 'lcopy: N.lt.0'
         Write (6,*) 'N=',N
         Call QTrace
         Call Abend()
      End If
      If ( incX.eq.1 .and. incY.eq.1 ) then
         M=Mod(N,4)
         If ( M.gt.0 ) then
            Do 10 i=1,M
               Y(i)=X(i)
10          Continue
         End If
         Do 20 i=M+1,N,4
            Y(i  )=X(i  )
            Y(i+1)=X(i+1)
            Y(i+2)=X(i+2)
            Y(i+3)=X(i+3)
20       Continue
      Else
        iX = 1
        iY = 1
        If (incX.lt.0 ) iX=(-N+1)*incX+1
        If (incY.lt.0 ) iY=(-N+1)*incY+1
        Do 30 i=1,N
          Y(iY)=X(iX)
          iX=iX+incX
          iY=iY+incY
30      Continue
      End If
*----------------------------------------------------------------------*
*     End procedure                                                    *
*----------------------------------------------------------------------*
      Return
      End
