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
      Subroutine Step2(iMatrix,SMatrix,nDim,TMatrix,iType,
     &                 SMatrix_Save,Temp)
*                                                                      *
************************************************************************
*                                                                      *
*     Step 2. LO S0 ->S1
*     occupied virtual on the same center
*     and orthonormalization of the original AO basis
*     The call to lowdin gives me a transforation A1, and a transformed
*     S matrix
*
      Implicit ReaL*8 (A-H,O-Z)
      Real*8 SMatrix(nDim*nDim),TMatrix(nDim*nDim),
     &       Temp(nDim*nDim), SMatrix_Save(nDim*nDim)
      Integer iMatrix(nDim),iType(nDim)
#include "real.fh"
*
clg     write (*,*) 'Step 2', nDim
clg     Call RecPrt('Save before LW 2',' ',SMatrix_Save,nDim,nDim)
clg     write (*,*)
         k=0
      Do i=1,nDim
         Do j=1,nDim
             k=k+1
c               write (*,*) iMatrix(i),iMatrix(j), SMatrix(k)
             if ((iMatrix(i).ne.iMatrix(j)) .and.
     &        (iType(i).ne.iType(j)))  then
                SMatrix(k)=Zero
            endif
         End Do
      End Do
clg   Call RecPrt('SMatrix before LW 2',' ',SMatrix,nDim,nDim)
clg   Call RecPrt('SMatrix_Save before LW 2',' ',SMatrix_Save,nDim,nDim)


      call dcopy_(nDim**2,[Zero],0,TMatrix,1)
      call dcopy_(nDim,[One],0,TMatrix,nDim+1)
      Call Lowdin(SMatrix,TMatrix,nDim)
*     Pick up S2
clg   Call RecPrt('SMatrix after LW 2',' ',SMatrix,nDim,nDim)
clg   Call RecPrt('TMatrix after LW 2',' ',TMatrix,nDim,nDim)
      call dcopy_(nDim**2,SMatrix_Save,1,SMatrix,1)
*
*     Now apply T2 to S2:  S3=T2(T)*S2*T2
*
      Call FZero(Temp,nDim**2)
      Call DGEMM_('N','N',
     &            nDim,nDim,nDim,
     &            1.0d0,SMatrix,nDim,
     &            TMatrix,nDim,
     &            0.0d0,Temp,nDim)
      Call DGEMM_('T','N',
     &            nDim,nDim,nDim,
     &            1.0d0,TMatrix,nDim,
     &            Temp,nDim,
     &            0.0d0,SMatrix,nDim)
C     Call RecPrt('S3',' ',Work(ip_s),nBas(1),nBas(1))
      call dcopy_(nDim**2,SMatrix,1,SMatrix_Save,1)


      Return
      End
