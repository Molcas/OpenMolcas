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
* Copyright (C) 2023, Ignacio Fdez. Galvan                             *
************************************************************************

      Subroutine Prepare_Kriging(Model_E,Model_G,nData,nDim,iFirst)
*     Prepare energy and gradients in the right format for Setup_Kriging

      Use Kriging_mod, Only: nSet
      Use Slapaf_Info, Only: dqInt, dqInt_Aux, Energy, Energy0
      Use Constants, Only: Zero
      Use Definitions, Only: wp, iwp

      Implicit None
      Integer(kind=iwp), Intent(in) :: nData, nDim, iFirst
      Real(kind=wp), Intent(out) :: Model_E(nData,nSet),
     &                              Model_G(nDim,nData,nSet)
      Integer(kind=iwp) :: i, iLast

      iLast = iFirst+nData-1
*
*     Trivial case, single surface: just copy energy and gradient
*
      Model_E(:,1) = Energy(iFirst:iLast)
      Model_G(:,:,1) = -dqInt(:,iFirst:iLast)
*
*     Get data for additional surfaces
*
      Do i=2,nSet
        If (i == 2) Then
          Model_E(:,i) = Energy0(iFirst:iLast)
        Else
          Model_E(:,i) = Zero
        End If
        Model_G(:,:,i) = -dqInt_Aux(:,iFirst:iLast,i-1)
      End Do

      End Subroutine Prepare_Kriging
