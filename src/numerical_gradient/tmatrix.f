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
* Copyright (C) 2015, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine TMatrix(TMx)
      Implicit None
************************************************************************
*     subroutine to get the T matrix that defines the constrained and  *
*     unconstrained subspaces.                                         *
************************************************************************
#include "info_slapaf.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8, Intent(InOut) :: TMx(mInt,mInt)
*
      Integer Lambda1,Lambda2
      Logical Invert
      Real*8, Allocatable, Dimension(:,:) :: C1,C2,CT
*
*     Get the global constraint vectors
*
      Call mma_Allocate(C1,mInt,nLambda)
      Call get_drdq(C1,Lambda1)
*
*     Get the NG constraint vectors
*
      Call Merge_Constraints('UDC.NG','','UDC',nLambda,iRow_c)
      Call Fix_UDC(iRow_c,nLambda,AtomLbl,nsAtom,nStab,.True.)
      Call mma_Allocate(C2,mInt,nLambda)
      Call get_drdq(C2,Lambda2)
*
*     Combine both sets of constraints and get the T matrix
*
      nLambda=Lambda1+Lambda2
      Call mma_Allocate(CT,mInt,mInt)
      Call dCopy_(mInt*Lambda1,C1,1,CT(1,1),1)
      Call dCopy_(mInt*Lambda2,C2,1,CT(1,Lambda1+1),1)
      Call FZero(CT(1,nLambda+1),mInt*(mInt-nLambda))
      If (nLambda.gt.0) Then
        Call GS(CT,nLambda,TMx,mInt,.False.,.True.)
      Else
        Call FZero(TMx,mInt**2)
        Call dCopy_(mInt,[One],0,TMx,mInt+1)
      End If
*
*     If NG constraints are to be inverted, combine the complement
*     with the global constraints instead
*
      Call Qpg_iScalar('Invert constraints',Invert)
      If (Invert) Call Get_lScalar('Invert constraints',Invert)
      If (Invert) Then
        Lambda2=mInt-nLambda
        Call dCopy_(mInt*Lambda1,C1,1,CT(1,1),1)
        Call dCopy_(mInt*Lambda2,TMx(1,nLambda+1),1,CT(1,Lambda1+1),1)
        nLambda=Lambda1+Lambda2
        Call FZero(CT(1,nLambda+1),mInt*(mInt-nLambda))
        If (nLambda.gt.0) Then
          Call GS(CT,nLambda,TMx,mInt,.False.,.True.)
        Else
          Call FZero(TMx,mInt**2)
          Call dCopy_(mInt,[One],0,TMx,mInt+1)
        End If
      End If
*
      Call mma_Deallocate(C1)
      Call mma_Deallocate(C2)
      Call mma_Deallocate(CT)
*
      End Subroutine TMatrix
