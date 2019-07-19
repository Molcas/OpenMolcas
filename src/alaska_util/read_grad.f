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
*  Read_Grad
*
*> @brief Read a vector from the gradients file
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Reads a gradient or coupling vector from the gradients file (GRADS).
*> For a gradient use \p iNAC, \p jNAC = 0. For a coupling vector use
*> \p iRoot = 0.
*>
*> @param[out] Grad  Gradient vector to read
*> @param[in]  nGrad Length of \p Grad
*> @param[in]  iRoot Root number to which the gradient belongs
*> @param[in]  iNAC  First root of a coupling vector
*> @param[in]  jNAC  Second root of a coupling vector
*>
*> @return \p 0 if the vector was not found, \p 1 if the vector was
*>         found, \p -1 if the vector was marked as non-computable
************************************************************************
      Function Read_Grad(Grad,nGrad,iRoot,iNAC,jNAC)
      Implicit None
#include "real.fh"
#include "stdalloc.fh"
      Integer :: Read_Grad,nGrad,iRoot,iNAC,jNAC
      Real*8 :: Grad(nGrad)
      Integer, Dimension(5) :: TOC
      Integer, Dimension(1) :: iDum
      Integer, Dimension(:), Allocatable :: i_grad,i_nac
      Integer :: nRoots,nCoup,LuGrad,iAd,iSt,jSt,idx
      Logical :: Found
      Character(Len=5) :: Filename
*
* If the GRADS file does not exist, there is no gradient
*
      Filename='GRADS'
      Call f_Inquire(Filename,Found)
      If (.Not.Found) Then
        Read_Grad=0
      Else
*
* Read the header
*
        LuGrad=20
        Call DaName(LuGrad,Filename)
        iAd=0
        Call iDaFile(LuGrad,2,TOC,Size(TOC),iAd)
        Call iDaFile(LuGrad,2,iDum,1,iAd)
        nRoots=iDum(1)
        If (Max(iRoot,iNAC,jNAC).gt.nRoots) Then
          Call WarningMessage(2,'Bad number of roots in GRADS file')
          Call Abend()
        End If
        Call iDaFile(LuGrad,2,iDum,1,iAd)
        If (iDum(1).ne.nGrad) Then
          Call WarningMessage(2,'Bad length in GRADS file')
          Call Abend()
        End If
        nCoup=Max(1,nRoots*(nRoots-1)/2)
        Call mma_Allocate(i_grad,nRoots)
        Call mma_Allocate(i_nac,nCoup)
        Call iDaFile(LuGrad,2,i_grad,nRoots,iAd)
        Call iDaFile(LuGrad,2,i_nac,nCoup,iAd)
*
* Read the gradient or NAC vector
*
        If (iRoot.eq.0) Then
          If ((iNAC.ne.0).and.(jNAC.ne.0)) Then
            iSt=Max(iNAC,jNAC)-1
            jSt=Min(iNAC,jNAC)
            idx=iSt*(iSt-1)/2+jSt
            iAd=i_nac(idx)
          Else
            iAd=-1
          End If
        Else
          idx=iRoot
          iAd=i_grad(idx)
        End If
*
        If (iAd.eq.0) Then
          Read_Grad=0
        Else If (iAd.lt.0) Then
          Read_Grad=-1
        Else
          Call dDaFile(LuGrad,2,Grad,nGrad,iAd)
          Read_Grad=1
        End If
*
        Call DaClos(LuGrad)
        Call mma_Deallocate(i_grad)
        Call mma_Deallocate(i_nac)
*
      End If
      If (Read_Grad.le.0) Then
        Call FZero(Grad,nGrad)
      End If
*
      End Function Read_Grad
