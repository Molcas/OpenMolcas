************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Hess_Vec(nAtoms,Hess,EVec,mAtoms,nDim)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Hess(3*nAtoms*(3*nAtoms+1)/2),EVec((3*mAtoms)**2)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the eigenvalues and the eigenvectors of the Hessian
*
*define _DEBUG_
*
*---- Set up a unit matrix
*
      nQ = nDim
      call dcopy_(nQ*nQ,[Zero],0,EVec,1)
      call dcopy_(nQ,[One],0,EVec,nQ+1)
*
*---- Compute eigenvalues and eigenvectors
*
      Call NIDiag_new(Hess,EVec,nQ,nQ,0)
      Call JacOrd(Hess,EVec,nQ,nQ)
*
      ThrD=1.0D-10
      Do iQ = 1, nQ
*        Fix standard direction.
         rZ = 0.0D0
         Do iElem = 1, nQ
            ij=(iQ-1)*nQ + iElem
            If (Abs(EVec(ij)).gt.Abs(rZ)+ThrD)
     &         rZ = EVec(ij)
         End Do
         ij=(iQ-1)*nQ + 1
         If (rZ.lt.0.0D0) Call DScal_(nQ,-One,EVec(ij),1)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Call RecPrt(' Eigenvectors','(12f6.2)',EVec,nDim,nDim)
      Call TriPrt(' Eigenvalues','(12E8.2)',Hess,nDim)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
