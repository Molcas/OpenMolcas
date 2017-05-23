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
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine Remove_TR(nQ,nX,nQQ,KMat,nK,TRVec,nTR,BM,iBM,nqB,nB)
      Implicit None
      Integer :: nQ,nX,nQQ,nK,nTR,nB
      Real*8 :: KMat(nQ,nK),TRVec(nX,nTR),BM(nB)
      Integer :: iBM(nB),nqB(nQ)
#include "real.fh"
#include "stdalloc.fh"
      Integer :: i,iK,iX,iQ,iB,iTR,iV
      Real*8, Dimension(:), Allocatable :: VecInt,TR
      Real*8, External :: DDot_

      Call mma_allocate(TR,nK)
      Call mma_allocate(VecInt,nX)

      Call FZero(TR,nK)

      Do iK=1,nK
*
*       Get the normalized Cartesian vector for this internal coordinate
        Call FZero(VecInt,nX)
        iB=0
        Do iQ=1,nQ
          Do i=1,nqB(iQ)
            iB=iB+1
            iX=iBM(iB)
            VecInt(iX)=VecInt(iX)+KMat(iQ,iK)*BM(iB)
          End Do
        End Do
        Call dScal_(nX,One/Sqrt(dDot_(nX,VecInt,1,VecInt,1)),VecInt,1)
*
*       Compute the overlap with the external translations and rotations
        Do iTR=1,nTR
          TR(iK)=TR(iK)+dDot_(nX,VecInt,1,TRvec(1,iTR),1)**2
        End Do
        write(6,*) 'vibration',iK,TR(iK)
        Do i=1,nX,3
          write(6,*) VecInt(i:i+2)
        End Do
      End Do
*
*     Put the nK-nQQ vectors with largest overlap at the end
      Do iK=nK,nQQ+1,-1
        iV=MaxLoc(TR(1:iK),1)
        write(6,*) 'IFG ',iK,iV
        If (iV.ne.iK) Call dSwap_(nQ,KMat(1,iK),1,KMat(1,iV),1)
      End Do

      Return
      End

