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
      Subroutine Step1(iCenter,Matrix,nDim,TMatrix,iType,Matrix0,
     &                 Temp)
*                                                                      *
************************************************************************
*                                                                      *
*     Step 1. LO S0 ->S1
*     occupied virtual on the same center
*     and orthonormalization of the original AO basis
*     The call to lowdin gives me a transforation A1, and a transformed
*     S matrix
*
      Implicit ReaL*8 (A-H,O-Z)
      Real*8 Matrix(nDim*nDim),TMatrix(nDim*nDim),Matrix0(nDim*nDim),
     &       Temp(nDim*nDim)
      Integer iCenter(nDim),iType(nDim)
#include "real.fh"
*
      k=0
      Do i=1,nDim
         Do j=1,nDim
             k=k+1
c               write (*,*) iCenter(i),iCenter(j), Matrix(k)
             if (iCenter(i).ne.iCenter(j)) then
                Matrix(k)=Zero
            endif
         End Do
      End Do
      Call GramSchmidt(Matrix,TMatrix,nDim,iType,iCenter,0)
      call dcopy_(nDim**2,Matrix0,1,Matrix,1)
*
*     Now apply T1 to original S: S2=T1(T)*S*T1
*
      Call DGEMM_('N','N',
     &            nDim,nDim,nDim,
     &            1.0d0,Matrix,nDim,
     &            TMatrix,nDim,
     &            0.0d0,Temp,nDim)
      Call DGEMM_('T','N',
     &            nDim,nDim,nDim,
     &            1.0d0,TMatrix,nDim,
     &            Temp,nDim,
     &            0.0d0,Matrix,nDim)
*
      Return
      End
