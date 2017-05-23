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
* Copyright (C) 1996, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Fold(nSym,nBas,A,B)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Fold up symmetry blocked matrix A stored as square blocks        *
*     into triangular matrices and scale the off-diagonal elements     *
*     by a factor of 2 (two). The input and output matrices may be     *
*     identical.                                                       *
*                                                                      *
*     calling arguments:                                               *
*     nSym    : input, integer                                         *
*               number of symmetry blocks                              *
*     nBas    : input, array of integers                               *
*               matrix dimension per symmetry block                    *
*     A       : input, array of real*8                                 *
*               Unfolded input matrix                                  *
*     B       : output, array of real*8                                *
*               Folded output matrix                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Dimension nBas(*) , A(*) , B(*)

      Parameter ( Two=2.0d0 )

      iOff1 = 0
      iOff2 = 0
      Do iSym = 1, nSym
        mBas = nBas(iSym)
        Do iBas= 1, mBas
          Do jBas = 1 , iBas-1
            B(iOff2+jBas) =  Two * A(iOff1+jBas)
          End Do
          B(iOff2+iBas) =  A(iOff1+iBas)
          iOff1 = iOff1 + mBas
          iOff2 = iOff2 + iBas
        End Do
      End Do

      Return
      End
