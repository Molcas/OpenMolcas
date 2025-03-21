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
      Subroutine Dmat_MCLR(CMO,OCC,D)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Compute the active one-body density                              *
*                                                                      *
*     calling arguments:                                               *
*     CMO     : input, array of real*8                                 *
*               MO-coefficients                                        *
*     OCC     : input, array of real*8                                 *
*               occupation numbers                                     *
*     D       : output, array of real*8                                *
*               total one-body density                                 *
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
      use Constants, only: Zero
      use input_mclr, only: nSym,nBas
      Implicit None

      Real*8 CMO(*) , OCC(*) , D(*)

      Integer iOff1,iOff2,iOff3,iSym,iBas,i,ii,j,k
      Real*8 Sum

      iOff1 = 0
      iOff2 = 0
      iOff3 = 0
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        If ( iBas.ne.0 ) then
          Do i = 1,iBas
            ii = (i*i-i)/2
            Do j = 1,i
              Sum = Zero
              Do k = 1,ibas
                Sum = Sum + OCC(iOff3+k)
     &                    * CMO(iOff1+(k-1)*iBas+i)
     &                    * CMO(iOff1+(k-1)*iBas+j)
              End Do
              D(iOff2+ii+j) = 2.0d0*Sum
              If (j.eq.i) D(iOff2+ii+j) = Sum
            End Do
          End Do
        End If
        iOff1 = iOff1 + iBas*iBas
        iOff2 = iOff2 + (iBas*iBas+iBas)/2
        iOff3 = iOff3 + iBas
      End Do

      End Subroutine Dmat_MCLR
