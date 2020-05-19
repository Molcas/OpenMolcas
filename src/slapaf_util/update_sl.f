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
      Subroutine Update_sl(iter,MaxItr,NmIter,iInt,nFix,nInter,qInt,
     &                     Shift,
     &                     Grad,iOptC,Beta,Beta_Disp,Lbl,GNrm,
     &                     Energy,UpMeth,ed,Line_Search,Step_Trunc,
     &                     nLambda,iRow_c,nsAtom,AtomLbl,nSym,iOper,
     &                     mxdc,jStab,nStab,BMx,Smmtrc,nDimBC,
     &                     rLambda,ipCx,GrdMax,StpMax,GrdLbl,StpLbl,
     &                     iNeg,nLbl,Labels,nLabels,FindTS,TSC,nRowH,
     &                     nWndw,Mode,MF,
     &                     iOptH,HUpMet,kIter,GNrm_Threshold,IRC,
     &                     dMass,HrmFrq_Show,CnstWght,Curvilinear,
     &                     Degen)
************************************************************************
*                                                                      *
*     Object: to update coordinates                                    *
*                                                                      *
*    Input:                                                            *
*      iter           : iteration counter                              *
*      MaxItr         : max number of iterations                       *
*      NmIter         : number of iteration in numerical approach      *
*      iInt           : number of internal coordinates to vary         *
*      nFix           : number of frozen internal coordinates          *
*      nInter         : total number of internal coordinates           *
*      qInt(*,iter )  : the internal coordinates                       *
*      Grad(*,iter )  : the gradient in the internal coordinates       *
*      iOptC          : option flag for update methods                 *
*      Beta           : damping factor                                 *
*      Lbl            : character labels for internal coordinates      *
*      nLbl           : length of Lbl                                  *
*      GNrm           : the norm of the gradient in each iteration     *
*      Energy         : the energy of each iteration                   *
*      Line_Search    : logical flag for line search                   *
*      nLambda        : number of contraints                           *
*      iRow_c         : number of lines on the UDC file                *
*      nsAtom         : number of symmetry unique atoms                *
*      AtomLbl        : character string with atom labels              *
*      nSym           : number of irreps                               *
*      iOper          : integer representations of symmetry operators  *
*      mxdc           : max number of nsAtom                           *
*      jStab          : integer list of stabilizers                    *
*      nStab          : number of stabilizers                          *
*      BMx            : the so-called Wilson B matrix                  *
*      Smmtrc         : logical flag for symmetry properties           *
*      nDimBC         : dimension of redundant coordinates(?)          *
*      rLambda        : vector for Lagrange multipliers                *
*      ipCx           : pointer to cartesian coordinates               *
*      iNeg           : Hessian index                                  *
*      Labels         : character string of primitive int. coord.      *
*      nLabels        : length of Labels                               *
*      CnstWght       : constraints weight                             *
*                                                                      *
*    OutPut:                                                           *
*      Shift(*,iter ) : the shift of the internal coordinates          *
*      qInt(*,iter +1): the internal coordinates to be used in the     *
*                       next iteration                                 *
*      UpMeth         : character label with update method abrivation  *
*      ed             : estimated energy change to the next point      *
*      Step_Trunc     : character label to denote truncation type      *
*      GrdMax         : largest gradient component                     *
*      StpMax         : largest step component                         *
*      GrdLbl         : character label of component with GrdMax       *
*      StpLbl         : character label of component with StpLbl       *
*                                                                      *
*                                                                      *
*     Author: Roland Lindh                                             *
*             2000                                                     *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "Molcas.fh"
      Real*8 qInt(nInter,MaxItr), Shift(nInter,MaxItr),
     &       Grad(nInter,MaxItr), GNrm(MaxItr), Energy(MaxItr),
     &       BMx(3*nsAtom,3*nsAtom), rLambda(nLambda,MaxItr),
     &       dMass(nsAtom), Degen(3*nsAtom), MF(3*nsAtom)
      Integer iOper(0:nSym-1), jStab(0:7,nsAtom), nStab(nsAtom),
     &        iNeg(2)
      Logical Line_Search, Smmtrc(3*nsAtom),
     &        FindTS, TSC, HrmFrq_Show, Curvilinear
      Character Lbl(nLbl)*8, GrdLbl*8, StpLbl*8, Step_Trunc,
     &          Labels(nLabels)*8, AtomLbl(nsAtom)*(LENIN), UpMeth*6,
     &          HUpMet*6
      Real*8, Allocatable:: t_Shift(:,:), t_qInt(:,:)
*
      Logical Kriging_Hessian
*
      iRout=153
      iPrint=nPrint(iRout)
      Call QEnter('Update')
*
      If (iPrint.ge.99) Then
         Call RecPrt('Update: qInt',' ',qInt,nInter,Iter)
         Call RecPrt('Update: Energy',' ',Energy,1,Iter)
         Call RecPrt('Update: Grad',' ',Grad,nInter,Iter)
*        Call RecPrt('Update: Shift',' ',Shift,nInter,Iter-1)
*        Call RecPrt('Update: GNrm',' ',GNrm,Iter,1)
      End If
*
      iOpt_RS=0
      Kriging_Hessian =.FALSE.
      qBeta=Beta
      qBeta_Disp=Beta_Disp
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
         If (iPrint.ge.99) Write(6,*)'UpDate_SL: first iteration'
         iter_=1
         Call mma_Allocate(t_Shift,nInter,iter_,Label='t_Shift')
         Call mma_Allocate(t_qInt,nInter,iter_+1,Label='t_qInt')
*
         t_qInt(:,1)=qInt(:,1)
*
         Call Update_sl_(iter_,iInt,nFix,nInter,t_qInt,
     &                   t_Shift,Grad,iOptC,Beta,Beta_Disp,
     &                   Lbl,GNrm,Energy,UpMeth,ed,Line_Search,
     &                   Step_Trunc,nLambda,iRow_c,nsAtom,AtomLbl,nSym,
     &                   iOper,mxdc,jStab,nStab,BMx,Smmtrc,nDimBC,
     &                   rLambda,ipCx,GrdMax,StpMax,GrdLbl,StpLbl,
     &                   iNeg,nLbl,Labels,nLabels,FindTS,TSC,nRowH,
     &                   nWndw,Mode,MF,
     &                   iOptH,HUpMet,kIter,GNrm_Threshold,IRC,dMass,
     &                   HrmFrq_Show,CnstWght,Curvilinear,Degen,
     &                   Kriging_Hessian,qBeta,iOpt_RS,.True.,iter_,
     &                   qBeta_Disp)
*                                                                      *
************************************************************************
*                                                                      *
*------- Move new coordinates to the correct position and compute the
*        corresponding shift.
*
         qInt(:,iter+1)=t_qInt(:,2)
         Shift(:,iter)=t_qInt(:,2)-qInt(:,iter)
*
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
         Call Update_sl_(iter,iInt,nFix,nInter,qInt,Shift,
     &                Grad,iOptC,Beta,Beta_Disp,Lbl,GNrm,Energy,
     &                UpMeth,ed,Line_Search,Step_Trunc,nLambda,
     &                iRow_c,nsAtom,AtomLbl,nSym,iOper,mxdc,jStab,
     &                nStab,BMx,Smmtrc,nDimBC,rLambda,ipCx,
     &                GrdMax,StpMax,GrdLbl,StpLbl,iNeg,nLbl,
     &                Labels,nLabels,FindTS,TSC,nRowH,
     &                nWndw,Mode,MF,
     &                iOptH,HUpMet,kIter,GNrm_Threshold,IRC,dMass,
     &                HrmFrq_Show,CnstWght,Curvilinear,Degen,
     &                Kriging_Hessian,qBeta,iOpt_RS,.True.,iter,
     &                qBeta_Disp)
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
      If (iPrint.ge.99) Then
         Call RecPrt(
     &      'Shifts in internal coordinate basis / au or rad',
     &      ' ',Shift,nInter,Iter)
         Call RecPrt(
     &      'qInt in internal coordinate basis / au or rad',
     &      ' ',qInt,nInter,Iter+1)
      End If
*
*---- Remove unneeded fields from the runfile
      Dummy=Zero
      Call Put_dArray('BMxOld',Dummy,0)
      Call Put_dArray('TROld',Dummy,0)
*
      Call QExit('Update')
      Return
      End
