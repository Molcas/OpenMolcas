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
      SubRoutine WrHDsk(Hess,ngrad)
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "disp2.fh"
#include "real.fh"
      real*8 Hess(nGrad*(nGrad+1)/2)
      Character*8 Label
      Real*8, Allocatable:: Temp(:), HStat(:), EVec(:), EVal(:)
*
      Call mma_allocate(Temp,nGrad**2,Label='Temp')
      nH=0
      Do iIrrep=0,nIrrep-1
         nH=nH+lDisp(iIrrep)
      End Do
      Call mma_allocate(HStat,nH,Label='HStat')
*
*---- Reorder Hessian to lower triangular form
*
      iGrad1=1
      iGrad2=0
      iG=0
      ip_Acc=1
      Do iIrrep=0,nIrrep-1
         iGrad2=iGrad2+lDisp(iIrrep)
*
         Do iG1=iGrad1,iGrad2
           Do iG2=iGrad1,iG1
              iG=iG+1
              Temp(iG)=Hess(iG1*(iG1-1)/2+IG2)
           End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*------- Diagonalize and keep eigen values for check facility
*
         mH=lDisp(iIrrep)
         Call mma_allocate(EVal,mH*(mH+1)/2,Label='EVal')
         Call mma_allocate(EVec,mH*mH,Label='EVec')
*
         call dcopy_(mH*(mH+1)/2,Temp,1,EVal,1)
         call dcopy_(mH*mH,[Zero],0,EVec,1)
         call dcopy_(mH,[One],0,EVec,mH+1)
*
*------- Compute eigenvalues and eigenvectors
*
         Call Jacob(EVal,EVec,mH,mH)
         Call Jacord(EVal,EVec,mH,mH)
*
         Do i = 1, mH
            HStat(ip_Acc)=EVal(i*(i+1)/2)
            ip_Acc=ip_Acc+1
         End Do
*
         Call mma_deallocate(EVec)
         Call mma_deallocate(EVal)
*                                                                      *
************************************************************************
*                                                                      *
         iGrad1=iGrad1+lDisp(iIrrep)
      End Do
*
*---- Write eigen values to the check file.
*
      Call Add_Info('HStat',HStat,nH,5)
*
      iRc=-1
      iOpt=0
      Label='StatHess'
      Call dWrMck(iRC,iOpt,Label,idum,Temp,idum)
      If (iRc.ne.0) Then
         Write (6,*) 'WrHDsk: Error writing to MCKINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
      Call mma_deallocate(HStat)
      Call mma_deallocate(Temp)
      Return
      End
