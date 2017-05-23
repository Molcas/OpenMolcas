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
* Copyright (C) 1999, Roland Lindh                                     *
*               2012, Andy May                                         *
************************************************************************
************************************************************************
*                                                                      *
*----- Statement functions to access the k2 data in the k2 data        *
*      memory block                                                    *
*                                                                      *
*      These functions should eventually be used globally for any      *
*      seward utility, be it integrals, gradients, direct Fock,        *
*      direct integrals, and second order derivatives.                 *
*                                                                      *
*                          W A R N I N G !                             *
*                                                                      *
*                                                                      *
*      Observe that the index array (IndZ) should always be placed     *
*      after all arrays with real!                                     *
*                                                                      *
*      Real Arrays:                                                    *
*      alpha  : alpha Gaussian exponent                                *
*      beta   : beta  Gaussian exponent                                *
*      Z      : alpha+beta Gaussian exponent                           *
*      Kappa  : Common factor from Gaussian product theorem            *
*      Pcoor  : Coordinates of Gaussian product                        *
*      ab     : Max(|(ab|ab)|^{1/2})    (primitive basis)              *
*               over all angular components                            *
*      abCon  : Max(|(ab|ab)|^{1/2})* C ("contracted basis")           *
*               over all angular components                            *
*                                                                      *
*      Integer Arrays:                                                 *
*      IndZ   : Pair index (rectangular indexation), the last          *
*               entry in the array contains the effective number       *
*               of elements                                            *
*                                                                      *
*      Scalars:                                                        *
*      EstI   : Estimated largest contracted integral |(ab|ab)|^{1/2}  *
*      ZtMax  : Z of the largest abCon                                 *
*      abMax  : largest abCon                                          *
*      ZetaM  : largest Zeta value                                     *
*      ZtMaxD : Z of the largest ab * D                                *
*      abMaxD : largest ab * D                                         *
*                                                                      *
*      Auxiliary arrays:                                               *
*      HrrMtrx: matrices to use for transformation from intermediate   *
*               integrals to real spherical harmonics                  *
*                                                                      *
*      Author: R. Lindh                                                *
*              Dept. of Chemical Physics                               *
*              University of Lund, Sweden                              *
*              April 1999                                              *
*                                                                      *
*      Converted from statement functions by A. May June 2012          *
************************************************************************
*                                                                      *
      Integer Function ip_Z(iZeta,nZeta)
      Implicit None
      Integer iZeta,nZeta
      ip_Z = iZeta
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nZeta)
      End

      Integer Function ip_Kappa(iZeta,nZeta)
      Implicit None
      Integer iZeta,nZeta
      ip_Kappa = iZeta + nZeta
      Return
      End

      Integer Function ip_Pcoor(iZeta,nZeta)
      Implicit None
      Integer iZeta,nZeta
      ip_Pcoor = iZeta + 2*nZeta
      Return
      End

      Integer Function ip_ZInv(iZeta,nZeta)
      Implicit None
      Integer iZeta,nZeta
      ip_ZInv = iZeta + 5*nZeta
      Return
      End

      Integer Function ip_ab(iZeta,nZeta)
      Implicit None
      Integer iZeta,nZeta
      ip_ab = iZeta + 6*nZeta
      Return
      End

      Integer Function ip_abCon(iZeta,nZeta)
      Implicit None
      Integer iZeta,nZeta
      ip_abCon = iZeta + 7*nZeta
      Return
      End

      Integer Function ip_Alpha(iZeta,nZeta,iAB)
      Implicit None
      Integer iZeta,nZeta,iAB
      ip_Alpha = iZeta + (8+iAB-1)*nZeta
      Return
      End

      Integer Function ip_Beta(iZeta,nZeta,iAB)
      Implicit None
      Integer iZeta,nZeta,iAB
      ip_Beta = iZeta + (8+iAB-1)*nZeta
      Return
      End

      Integer Function ip_IndZ(iZeta,nZeta)
      Implicit None
      Integer iZeta,nZeta
      ip_IndZ = iZeta + 10*nZeta
      Return
      End

      Integer Function ip_EstI(nZeta)
      Implicit None
      Integer nZeta
      ip_EstI = 1 + 11*nZeta + 1
      Return
      End

      Integer Function ip_ZtMax(nZeta)
      Implicit None
      Integer nZeta
      ip_ZtMax = 2 + 11*nZeta + 1
      Return
      End

      Integer Function ip_abMax(nZeta)
      Implicit None
      Integer nZeta
      ip_abMax = 3 + 11*nZeta + 1
      Return
      End

      Integer Function ip_ZetaM(nZeta)
      Implicit None
      Integer nZeta
      ip_ZetaM = 4 + 11*nZeta + 1
      Return
      End

      Integer Function ip_ZtMaxD(nZeta)
      Implicit None
      Integer nZeta
      ip_ZtMaxD = 5 + 11*nZeta + 1
      Return
      End

      Integer Function ip_abMaxD(nZeta)
      Implicit None
      Integer nZeta
      ip_abMaxD = 6 + 11*nZeta + 1
      Return
      End

      Integer Function ip_nHm(nZeta)
      Implicit None
      Integer nZeta
      ip_nHm = 7 + 11*nZeta + 1
      Return
      End

      Integer Function ip_HrrMtrx(nZeta)
      Implicit None
      Integer nZeta
      ip_HrrMtrx = 8 + 11*nZeta + 1
      Return
      End

      Integer Function ip_abG(nZeta,nHm)
      Implicit None
      Integer nZeta,nHm
      ip_abG = 8+nHm + 11*nZeta + 1
      Return
      End
*                                                                      *
************************************************************************
