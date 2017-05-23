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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      Subroutine Fold_Mat(nSym,nBas,A,B)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Fold up symmetry blocked matrix A stored as SQUARE blocks        *
*     into TRIANGULAR matrices. The matrix A can be non-symmetric.     *
*     The input and output matrices CANNOT BE IDENTICAL !!!            *
*                                                                      *
*     calling arguments:                                               *
*     nSym    : input, integer                                         *
*               number of symmetry blocks                              *
*     nBas    : input, array of integers                               *
*               matrix dimension per symmetry block                    *
*     A       : input, array of real*8  (square storage)               *
*               Unfolded input matrix                                  *
*     B       : output, array of real*8 (triangular storage)           *
*               Folded output matrix                                   *
*                                                                      *
*----------------------------------------------------------------------*
*     Author:   F. Aquilante                                           *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Dimension nBas(*) , A(*) , B(*)
*************************************
      iit(j) = j*(j+1)/2
*************************************
      ijt(i,j) = i*(i-1)/2 + j
*************************************
      iis(n,j) = n*(j-1) + j
*************************************
      ijs(n,i,j) = n*(j-1) + i
*************************************


      iOff1 = 0
      iOff2 = 0

      Do iSym=1,nSym

        Do j=1,nBas(iSym)

          B( iOff1 + iit(j) ) =  A( iOff2 + iis(nBas(iSym),j) )

          Do i=j+1,nBas(iSym)

            B( iOff1 + ijt(i,j) ) = A( iOff2 + ijs(nBas(iSym),i,j) )
     &                            + A( iOff2 + ijs(nBas(iSym),j,i) )

          End Do

        End Do

        iOff1 = iOff1 + iit(nBas(iSym))
        iOff2 = iOff2 + nBas(iSym)**2

      End Do

      Return
      End
