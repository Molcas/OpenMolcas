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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      Subroutine Update_sl(iter,NmIter,nQQ,Step_Trunc,nWndw,kIter)
************************************************************************
*                                                                      *
*     Object: to update coordinates                                    *
*                                                                      *
*    Input:                                                            *
*      iter           : iteration counter                              *
*      NmIter         : number of iteration in numerical approach      *
*      nQQ         : total number of internal coordinates           *
*      Beta           : damping factor                                 *
*                                                                      *
*    OutPut:                                                           *
*      Step_Trunc     : character label to denote truncation type      *
*                                                                      *
*                                                                      *
*     Author: Roland Lindh                                             *
*             2000                                                     *
************************************************************************
      use Slapaf_Info, only: Shift, qInt
      use Slapaf_Parameters, only: Beta, Beta_disp
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "Molcas.fh"
      Character Step_Trunc
      Real*8 Dummy(1)
      Real*8, Allocatable:: t_Shift(:,:), t_qInt(:,:), tmp(:)
*
      Logical Kriging_Hessian
*                                                                      *
************************************************************************
*                                                                      *
      iRout=153
      iPrint=nPrint(iRout)
*
#ifdef _DEBUGPRINT_
      Call RecPrt('Update_sl: qInt',' ',qInt,nQQ,Iter)
      Call RecPrt('Update_sl: Shift',' ',Shift,nQQ,Iter-1)
#endif
*
      iOpt_RS=0
      qBeta=Beta
      qBeta_Disp=Beta_Disp
*                                                                      *
************************************************************************
*                                                                      *
      Call Mk_Hss_Q()
      Kriging_Hessian =.FALSE.
*                                                                      *
************************************************************************
*                                                                      *
*     Select between numerical evaluation of the Hessian or a molcular
*     geometry optimization.
*
      Step_Trunc=' '
      If (iter.eq.NmIter.and.NmIter.ne.1) Then
*                                                                      *
************************************************************************
*                                                                      *
*------- On the first iteration after a numerical evaluation of the
*        Hessian we like the step to be relative to the initial
*        structure.
*
#ifdef _DEBUGPRINT_
         Write(6,*)'UpDate_SL: first iteration'
#endif
         iter_=1
         Call mma_Allocate(t_Shift,nQQ,SIZE(Shift,2),Label='t_Shift')
         t_Shift(:,:)=Shift(:,:)
         Shift(:,:)=Zero

         Call mma_Allocate(t_qInt,nQQ,SIZE(qInt,2),Label='t_qInt')
         t_qInt(:,:)=qInt(:,:)
         qInt(:,:)=Zero
         qInt(:,1)=t_qInt(:,1)
*
         Call Update_inner(iter_,nQQ,qInt,Shift,Beta,Beta_Disp,
     &                    Step_Trunc,nWndw,kIter,Kriging_Hessian,qBeta,
     &                    iOpt_RS,.True.,iter_,qBeta_Disp)
*                                                                      *
************************************************************************
*                                                                      *
*------- Move new coordinates to the correct position and compute the
*        corresponding shift.
*
         Call mma_allocate(Tmp,SIZE(qInt,1),Label='tmp')
         tmp(:)=qInt(:,2)
         qInt(:,:)=t_qInt(:,:)
         qInt(:,iter+1)=tmp(:)
         Shift(:,:)=t_Shift(:,:)
         Shift(:,iter)=tmp(:)-qInt(:,iter)

*
         Call mma_deallocate(tmp)
         Call mma_deallocate(t_qInt)
         Call mma_deallocate(t_Shift)
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
*        Conventional optimization.
*
         Call Update_inner(iter,nQQ,qInt,Shift,Beta,Beta_Disp,
     &                     Step_Trunc,nWndw,kIter,Kriging_Hessian,qBeta,
     &                     iOpt_RS,.True.,iter,qBeta_Disp)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Write out the shift in internal coordinates basis.
*
#ifdef _DEBUGPRINT_
      Call RecPrt('Shifts in internal coordinate basis / au or rad',
     &      ' ',Shift,nQQ,Iter)
      Call RecPrt('qInt in internal coordinate basis / au or rad',
     &      ' ',qInt,nQQ,Iter+1)
#endif

*
*---- Remove unneeded fields from the runfile
      Dummy(1)=-Zero
      Call Put_dArray('BMxOld',Dummy(1),0)
      Call Put_dArray('TROld',Dummy(1),0)
*
      Return
      End
