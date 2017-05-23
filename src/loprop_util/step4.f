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
      Subroutine Step4(SMatrix,nDim,TMatrix,iType)
*                                                                      *
************************************************************************
*                                                                      *
*     Step 4. LW S3 ->S4
*
      Implicit ReaL*8 (A-H,O-Z)
      Real*8 SMatrix(nDim*nDim),TMatrix(nDim*nDim)
      Integer iType(nDim)
#include "real.fh"
*
clg   write (*,*) 'Step 4', nDim
clg   Call RecPrt('T before LW 4',' ',TMatrix,nDim,nDim)
clg   Call RecPrt('S in step4 ',' ',SMatrix,nDim,nDim)
clg   write (*,*)
      k=0
      Do i=1,nDim
         Do j=1,nDim
            k=k+1
            If (i.ne.j .and.(iType(i).ne.iType(j))) then
               SMatrix(k)=Zero
            EndIf
         End Do
      End Do
clg   Call RecPrt('S before LW 4',' ',SMatrix,nDim,nDim)

      call dcopy_(nDim**2,Zero,0,TMatrix,1)
      call dcopy_(nDim,One,0,TMatrix,nDim+1)
      Call Lowdin(SMatrix,TMatrix,nDim)
*
      Return
      End
