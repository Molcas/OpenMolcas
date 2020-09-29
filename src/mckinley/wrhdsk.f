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
#include "WrkSpc.fh"
#include "disp.fh"
#include "disp2.fh"
#include "real.fh"
      real*8 Hess(nGrad*(nGrad+1)/2)
      Character*8 Label
*
      Call Getmem('Temp','ALLO','REAL',ipTemp,nGrad**2)
      nH=0
      Do iIrrep=0,nIrrep-1
         nH=nH+lDisp(iIrrep)
      End Do
      Call Getmem('HStat','ALLO','REAL',ipHStat,nH)
*
*---- Reorder Hessian to lower triangular form
*
      iGrad1=1
      iGrad2=0
      iG=0
      ip_Acc=ipHStat
      Do iIrrep=0,nIrrep-1
         iGrad2=iGrad2+lDisp(iIrrep)
*
         ipH=ipTemp
         Do iG1=iGrad1,iGrad2
           Do iG2=iGrad1,iG1
              iG=iG+1
              Work(ipTemp+iG-1)=Hess(iG1*(iG1-1)/2+IG2)
           End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*------- Diagonalize and keep eigen values for check facility
*
         mH=lDisp(iIrrep)
         Call GetMem('EVal','Allo','Real',ipEVal,mH*(mH+1)/2)
         Call GetMem('EVec','Allo','Real',ipEVec,mH*mH)
*
         call dcopy_(mH*(mH+1)/2,Work(ipH),1,Work(ipEVal),1)
         call dcopy_(mH*mH,[Zero],0,Work(ipEVec),1)
         call dcopy_(mH,[One],0,Work(ipEVec),mH+1)
*
*------- Compute eigenvalues and eigenvectors
*
         Call Jacob(Work(ipEVal),Work(ipEVec),mH,mH)
         Call Jacord(Work(ipEVal),Work(ipEVec),mH,mH)
*
         Do i = 1, mH
            Work(ip_Acc)=Work(i*(i+1)/2+ipEVal-1)
            ip_Acc=ip_Acc+1
         End Do
*
         Call GetMem('EVec','Free','Real',ipEVec,mH*mH)
         Call GetMem('EVal','Free','Real',ipEVal,mH*(mH+1)/2)
*                                                                      *
************************************************************************
*                                                                      *
         iGrad1=iGrad1+lDisp(iIrrep)
      End Do
*
*---- Write eigen values to the check file.
*
      Call Add_Info('HStat',Work(ipHStat),nH,5)
*
      iRc=-1
      iOpt=0
      Label='StatHess'
      Call dWrMck(iRC,iOpt,Label,idum,Work(ipTemp),idum)
      If (iRc.ne.0) Then
         Write (6,*) 'WrHDsk: Error writing to MCKINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
      Call Getmem('HStat','Free','REAL',ipHStat,nH)
      Call Getmem('Temp','Free','REAL',ipTemp,nGrad**2)
      Return
      End
