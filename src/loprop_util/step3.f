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
      Subroutine Step3(iCenter,SMatrix,nDim,TMatrix,iType)
*                                                                      *
************************************************************************
*                                                                      *
*     Step 3. GS S2 ->S3
*
      Implicit ReaL*8 (A-H,O-Z)
      Real*8 SMatrix(nDim*nDim),TMatrix(nDim*nDim)
      Integer iCenter(nDim),iType(nDim)
#include "real.fh"
*
clg   write (*,*) 'Step 3', nDim
clg   Call RecPrt('T before GS 3',' ',TMatrix,nDim,nDim)
clg   write (*,*)
      call dcopy_(nDim**2,Zero,0,TMatrix,1)
      call dcopy_(nDim,One,0,TMatrix,nDim+1)
      Call GramSchmidt(SMatrix,TMatrix,nDim,iType,iCenter,1)
*
      Return
      End
