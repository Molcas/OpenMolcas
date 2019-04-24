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
*  Store_Not_Grad
*
*> @brief Mark a vector as non-computable in the gradients file
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Marks a gradient or coupling vector in the gradients file (GRADS)
*> as non-computable. This is useful, for example, to allow SLAPAF
*> to proceed with an approximate coupling vector if it cannot be
*> computed with the current method.
*> For a gradient use \p iNAC, \p jNAC = 0. For a coupling vector use
*> \p iRoot = 0.
*>
*> @param[in] iRoot Root number to which the gradient belongs
*> @param[in] iNAC  First root of a coupling vector
*> @param[in] jNAC  Second root of a coupling vector
************************************************************************
      Subroutine Store_Not_Grad(iRoot,iNAC,jNAC)
      Implicit None
#include "stdalloc.fh"
      Integer :: nGrad,iRoot,iNAC,jNAC
      Integer, Dimension(5) :: TOC
      Integer, Dimension(1) :: Length
      Integer, Dimension(:), Allocatable :: i_grad,i_nac
      Integer :: nRoots,nCoup,LuGrad,iAd,iSt,jSt,idx
      Logical :: Found
      Character(Len=5) :: Filename
*
* Create GRADS file if it does not exist
*
      Call Get_iScalar('Number of roots',nRoots)
      Call Get_iScalar('Unique atoms',nGrad)
      nGrad=3*nGrad
      Filename='GRADS'
      LuGrad=20
      Call f_Inquire(Filename,Found)
      If (.Not.Found) Call Create_Grads(Filename,nRoots,nGrad)
*
* Read the header
*
      Call DaName(LuGrad,Filename)
      iAd=0
      Call iDaFile(LuGrad,2,TOC,Size(TOC),iAd)
      Call iDaFile(LuGrad,2,Length,1,iAd)
      If (Length(1).ne.nRoots) Then
        Call WarningMessage(2,'Bad number of roots in GRADS file')
        Call Abend()
      End If
      Call iDaFile(LuGrad,2,Length,1,iAd)
      If (Length(1).ne.nGrad) Then
        Call WarningMessage(2,'Bad length in GRADS file')
        Call Abend()
      End If
      nCoup=Max(1,nRoots*(nRoots-1)/2)
      Call mma_Allocate(i_grad,nRoots)
      Call mma_Allocate(i_nac,nCoup)
      Call iDaFile(LuGrad,2,i_grad,nRoots,iAd)
      Call iDaFile(LuGrad,2,i_nac,nCoup,iAd)
*
* Write the negative index that marks it as non-computable
*
      If (iRoot.eq.0) Then
        If ((iNAC.ne.0).and.(jNAC.ne.0)) Then
          iSt=Max(iNAC,jNAC)-1
          jSt=Min(iNAC,jNAC)
          idx=iSt*(iSt-1)/2+jSt
          i_nac(idx)=-1
          iAd=TOC(4)
          Call iDaFile(LuGrad,1,i_nac,nCoup,iAd)
        End If
      Else
        idx=iRoot
        i_grad(idx)=-1
        iAd=TOC(3)
        Call iDaFile(LuGrad,1,i_grad,nRoots,iAd)
      End If
*
      Call DaClos(LuGrad)
      Call mma_Deallocate(i_grad)
      Call mma_Deallocate(i_nac)
*
      End Subroutine Store_Not_Grad
