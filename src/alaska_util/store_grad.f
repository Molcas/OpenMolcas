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
*  Store_Grad
*
*> @brief Store a vector in the gradients file
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Stores a gradient or coupling vector in the gradients file (GRADS).
*> For a gradient use \p iNAC, \p jNAC = 0. For a coupling vector use
*> \p iRoot = 0.
*>
*> @param[in] Grad  Gradient vector to store
*> @param[in] nGrad Length of \p Grad
*> @param[in] iRoot Root number to which the gradient belongs
*> @param[in] iNAC  First root of a coupling vector
*> @param[in] jNAC  Second root of a coupling vector
************************************************************************
      Subroutine Store_Grad(Grad,nGrad,iRoot,iNAC,jNAC)
      Implicit None
#include "real.fh"
#include "stdalloc.fh"
      Integer :: nGrad,iRoot,iNAC,jNAC
      Real*8 :: Grad(nGrad)
      Real*8 , Dimension(:), Allocatable :: TmpGrad
      Integer, Dimension(5) :: TOC
      Integer, Dimension(:), Allocatable :: i_grad,i_nac
      Integer :: nRoots,nCoup,LuGrad,iAd,Length,iSt,jSt,idx
      Integer, External :: AixRm
      Logical :: Found, BadFile
      Character(Len=5) :: Filename
*
*define _DEBUG_
*
* Create GRADS file if it does not exist
*
      Call Get_iScalar('Number of roots',nRoots)
      Filename='GRADS'
      LuGrad=20
      Call f_Inquire(Filename,Found)
      If (.Not.Found) Call Create_Grads(Filename,nRoots,nGrad)
*
* Read the header
*
      BadFile=.False.
      Call DaName(LuGrad,Filename)
      iAd=0
      Call iDaFile(LuGrad,2,TOC,Size(TOC),iAd)
      iAd=TOC(1)
      Call iDaFile(LuGrad,2,Length,1,iAd)
      If (Length.ne.nRoots) BadFile=.True.
      iAd=TOC(2)
      Call iDaFile(LuGrad,2,Length,1,iAd)
      If (Length.ne.nGrad) BadFile=.True.
      If (BadFile) Then
        Call DaClos(LuGrad)
        If (AixRm('GRADS').ne.0) Call Abend()
        Call WarningMessage(1,'Number of roots and/or length of '//
     &       'gradients do not match, re-creating GRADS file')
        Call Create_Grads(Filename,nRoots,nGrad)
        Call DaName(LuGrad,Filename)
        iAd=0
        Call iDaFile(LuGrad,1,TOC,Size(TOC),iAd)
      End If
      nCoup=Max(1,nRoots*(nRoots-1)/2)
      Call mma_Allocate(i_grad,nRoots)
      Call mma_Allocate(i_nac,nCoup)
      iAd=TOC(3)
      Call iDaFile(LuGrad,2,i_grad,nRoots,iAd)
      iAd=TOC(4)
      Call iDaFile(LuGrad,2,i_nac,nCoup,iAd)
*
* Write the gradient or NAC vector
*
      If (iRoot.eq.0) Then
        If ((iNAC.ne.0).and.(jNAC.ne.0)) Then
          Call mma_Allocate(TmpGrad,nGrad)
          Call dCopy_(nGrad,Grad,1,TmpGrad,1)
          If (iNAC.lt.jNAC) Call dScal_(nGrad,-One,TmpGrad,1)
          iSt=Max(iNAC,jNAC)-1
          jSt=Min(iNAC,jNAC)
          idx=iSt*(iSt-1)/2+jSt
          If (i_nac(idx).eq.0) Then
            i_nac(idx)=TOC(5)
            Call dDaFile(LuGrad,1,TmpGrad,nGrad,TOC(5))
            iAd=0
            Call iDaFile(LuGrad,1,TOC,Size(TOC),iAd)
            iAd=TOC(4)
            Call iDaFile(LuGrad,1,i_nac,nCoup,iAd)
          Else
            iAd=i_nac(idx)
            Call dDaFile(LuGrad,1,TmpGrad,nGrad,iAd)
          End If
          Call mma_Deallocate(TmpGrad)
        End If
      Else
        idx=iRoot
        If (i_grad(idx).eq.0) Then
          i_grad(idx)=TOC(5)
          Call dDaFile(LuGrad,1,Grad,nGrad,TOC(5))
          iAd=0
          Call iDaFile(LuGrad,1,TOC,Size(TOC),iAd)
          iAd=TOC(3)
          Call iDaFile(LuGrad,1,i_grad,nRoots,iAd)
        Else
          iAd=i_grad(idx)
          Call dDaFile(LuGrad,1,Grad,nGrad,iAd)
        End If
      End If
*
#ifdef _DEBUG_
      write(6,*) 'iRoot, iNAC, jNAC:',iRoot,iNAC,jNAC
      write(6,*) 'grads:',i_grad
      write(6,*) 'nacs:',i_nac
#endif
      Call DaClos(LuGrad)
      Call mma_Deallocate(i_grad)
      Call mma_Deallocate(i_nac)
*
      End Subroutine Store_Grad
