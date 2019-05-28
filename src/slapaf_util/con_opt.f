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
* Copyright (C) 2003, Roland Lindh                                     *
************************************************************************
      Subroutine Con_Opt(r,drdq,T,dEdq,rLambda,q,dq,dy,dx,dEdq_,du,x,
     &                   dEdx,W,GNrm,nWndw,
     &                   Hess,nInter,nIter,iOptC,Mode,ipMF,
     &                   iOptH,HUpMet,jPrint_,Energy,nLambda,
     &                   mIter,nRowH,Err,EMx,RHS,iPvt,dg,A,nA,ed,
     &                   Beta,nFix,iP,UpMeth,
     &                   Line_Search,Step_Trunc,Lbl,GrdLbl,StpLbl,
     &                   GrdMax,StpMax,d2rdq2,nsAtom,IRC,CnstWght,
     &                   Restriction)
************************************************************************
*                                                                      *
*     Object: to perform an constrained optimization. The constraints  *
*             are optimized directly as a linear problem in the        *
*             subspace defined by the gradients of the constraints.    *
*             The minimization is performed in the complemental        *
*             subspace with the ordinary routines.                     *
*                                                                      *
*             L(q,l) = E(q) - l^T r(q)                                 *
*                                                                      *
*                                                                      *
*                                                                      *
*             Ref:                                                     *
*             "A Reduced-Restricted-Quasi-Newton-Raphson Method for    *
*              Locating and Optimizing Energy Crossing Points Between  *
*              Two Potential Energy Surfaces",                         *
*             J. M. Anglada and J. M. Bofill,                          *
*             J. Comput. Chem., 18, 992-1003 (1997)                    *
*                                                                      *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemical Physics,                 *
*             University of Lund, SWEDEN                               *
*             July '03                                                 *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      External Restriction
      Real*8 Restriction
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 r(nLambda,nIter), drdq(nInter,nLambda,nIter),
     &       T(nInter,nInter), dEdq(nInter,nIter),
     &       rLambda(nLambda,nIter+1), q(nInter,nIter+1),
     &       dq(nInter,nIter), dy(nLambda), dx(nInter-nLambda,nIter),
     &       dEdq_(nInter,nIter), du(nInter),
     &       x(nInter-nLambda,nIter+1),
     &       dEdx(nInter-nLambda,nIter),
     &       W(nInter-nLambda,nInter-nLambda),
     &       Hess(nInter,nInter),
     &       Energy(nIter),
     &       Err(nInter,nIter+1), EMx((nIter+1)**2), RHS(nIter+1),
     &       dg(mIter), A(nA), d2rdq2(nInter,nInter,nLambda)
      Integer iPvt(nInter+1), iP(nInter), iNeg(2)
      Logical Line_Search, Found, IRC_setup
      Character HUpMet*6, UpMeth*6, Step_Trunc*1, Lbl(nInter+nLambda)*8,
     &          GrdLbl*8, StpLbl*8, StpLbl_Save*8
*                                                                      *
************************************************************************
*                                                                      *
      jPrint=jPrint_
      If (jPrint.ge.99) Then
         Call RecPrt('r',' ',r,nLambda,nIter)
         Call RecPrt('drdq(orig)',' ',drdq,nInter,nLambda*nIter)
         Call RecPrt('dEdq',' ',dEdq,nInter,nIter)
         Call RecPrt('Energy',' ',Energy,nIter,1)
         Call RecPrt('q',' ',q,nInter,nIter+1)
         Call RecPrt('T',' ',T,nInter,nInter)
         Do iLambda = 1, nLambda
            Call RecPrt('d2rdq2',' ',d2rdq2(1,1,iLambda),nInter,nInter)
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
      ipTb=1
      ipTti=ipTb+nLambda
*
      yBeta=One
      gBeta=One
      xBeta=One
      dydy_last=Beta
      gg_last=Beta
      dxdx_last=Beta
      Sf=Sqrt(Two)
      dxdx=Zero
      Thr=1.0D-6
      Call GetMem('dqtmp','Allo','Real',ipdq,nInter)
      Call Get_iScalar('iOff_Iter',iOff_Iter)
      Do iIter = iOff_Iter+1, nIter
*                                                                      *
************************************************************************
*                                                                      *
*------- If we have that the derivative of the constraint is a null
*        vector replace it with an arbitrary vector which is orthogonal
*        to the gradient.
*
*        In case of the first macro iteration of an IRC search replace
*        it with the reaction vector.
*
         Do iLambda = 1, nLambda
            RR=Sqrt(DDot_(nInter,drdq(1,iLambda,iIter),1,
     &                          drdq(1,iLambda,iIter),1))
            If (RR.lt.1.0D-12) Then
*
               xBeta = xBeta*Half
*
               If (Abs(IRC).ne.0) Then
                  iMEP=0
                  Call Qpg_iScalar('nMEP',Found)
                  If (Found) Call Get_iScalar('nMEP',iMEP)
                  IRC_Setup=iIter.eq.1.and.iMEP.eq.0
               Else
                  IRC_Setup=.False.
               End If
*
               If (iIter.eq.nIter) Then
*
                  Write (6,*) 'Warning: constraint ',iLambda,
     &                        ' has a null vector, I''ll fix it!'
*
                  If (IRC_Setup.and.IRC.eq.1) Then
                     Write (6,*) ' IRC forward direction.'
                  Else If (IRC_Setup.and.IRC.eq.-1) Then
                     Write (6,*) ' IRC backward direction.'
                  End If
               End If
               r(iLambda,iIter)=Zero
*
               If (IRC_SetUp) Then
                  Call ReacQ(Work(ipMF),3*nsAtom,
     &                       dEdq(1,        iIter),nInter)
               End If
*
*              Try to use the transverse direction if the gradient is too small
*
               RR=Sqrt(DDot_(nInter,dEdq(1,iIter),1,dEdq(1,iIter),1))
               If (RR.lt.1.0D-12) Then
                  Call qpg_dArray('Transverse',Found,nTrans)
                  If (Found.and.(nTrans.eq.3*nsAtom)) Then
                     Call Allocate_Work(ipTrans,3*nsAtom)
                     Call Get_dArray('Transverse',Work(ipTrans),nTrans)
                     Call ReacQ(Work(ipTrans),3*nsAtom,
     &                          dEdq(1,iIter),nInter)
                     RR=Sqrt(DDot_(nInter,dEdq(1,iIter),1,
     &                                   dEdq(1,iIter),1))
                     Call DScal_(nInter,One/RR,dEdq(1,iIter),1)
                     Call Free_Work(ipTrans)
                  End If
               Else
                  Call DScal_(nInter,One/RR,dEdq(1,iIter),1)
               End If
*
               Do iInter = 1, nInter
c                 drdq(iInter,iLambda,iIter) =
c    &              dEdq(iInter,iIter)
                  drdq(iInter,iLambda,iIter) =
     &              Sign(One,dEdq(iInter,iIter))
                  If (Abs(dEdq(iInter,iIter)).le.1.0D-6)
     &               drdq(iInter,iLambda,iIter) = Zero
               End Do
               If (jPrint.ge.99)
     &         Call RecPrt('drdq(1)',' ',drdq,nInter,nLambda*nIter)
*
*------------- Orthogonalize against the gradient and all other
*              constraints
*
               RG=DDot_(nInter,drdq(1,iLambda,iIter),1,
     &                        dEdq(1,        iIter),1)
               Call DaXpY_(nInter,-RG,dEdq(1,        iIter),1,
     &                                drdq(1,iLambda,iIter),1)
               RR=Sqrt(DDot_(nInter,drdq(1,iLambda,iIter),1,
     &                             drdq(1,iLambda,iIter),1))
C              Call DScal_(nInter,One/RR,drdq(1,iLambda,iIter),1)
               Do jLambda = 1, nLambda
                  If (jLambda.ne.iLambda) Then
                     RG=DDot_(nInter,drdq(1,iLambda,iIter),1,
     &                              drdq(1,jLambda,iIter),1)
                     Call DaXpY_(nInter,-RG,drdq(1,jLambda,iIter),1,
     &                                     drdq(1,iLambda,iIter),1)
                     RR=Sqrt(DDot_(nInter,drdq(1,iLambda,iIter),1,
     &                                   drdq(1,iLambda,iIter),1))
                     Call DScal_(nInter,One/RR,drdq(1,iLambda,iIter),1)
                  End If
               End Do
*
               Continue
               If (jPrint.ge.99)
     &         Call RecPrt('drdq(2)',' ',drdq,nInter,nLambda*nIter)
            End If
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*------- Set up the T-matrix which separates the space into the
*        subspace in which the constraint is accomplished and
*        the complemental subspace.
*
         Call GS(drdq(1,1,iIter),nLambda,T,nInter,.False.,.False.)
C        Call RecPrt('T-Matrix',' ',T,nInter,nInter)
*                                                                      *
************************************************************************
*                                                                      *
*------- Compute delta y  (Eqn. 11)
*
         Call Allocate_Work(ipRTInv,nLambda**2)
         Call Allocate_Work(ipRT,nLambda**2)
         Call FZero(Work(ipRT),nLambda**2)
         Call DGEMM_('T','N',
     &               nLambda,nLambda,nInter,
     &               1.0d0,drdq(1,1,iIter),nInter,
     &               T(1,ipTb),nInter,
     &               0.0d0,Work(ipRT),nLambda)
         Call MInv(Work(ipRT),Work(ipRTInv),iSing,Det,nLambda)
         Call Free_Work(ipRT)
         Call FZero(dy,nLambda)
         Call DGEMM_('N','N',
     &               nLambda,1,nLambda,
     &               1.0d0,Work(ipRTInv),nLambda,
     &               r(1,iIter),nLambda,
     &               0.0d0,dy,nLambda)
         Call Free_Work(ipRTInv)
*
         Call DScal_(nLambda,-One,dy,1)
         If (jPrint.ge.99) Call RecPrt('dy',' ',dy,nLambda,1)
*
*                                                                      *
************************************************************************
*                                                                      *
*------- Compute new lambdas (Eqn. 17)
*
         Call Allocate_Work(ipTRT,nLambda*nInter)
         Call FZero(Work(ipTRT),nLambda*nInter)
         Call Allocate_Work(ipTRInv,nLambda**2)
         Call Allocate_Work(ipTR,nLambda**2)
         Call FZero(Work(ipTR),nLambda**2)
*
         Call DGEMM_('T','N',
     &               nLambda,nLambda,nInter,
     &               1.0d0,drdq(1,1,iIter),nInter,
     &               drdq(1,1,iIter),nInter,
     &               0.0d0,Work(ipTR),nLambda)
*
         Call MInv(Work(ipTR),Work(ipTRInv),iSing,Det,nLambda)
         Call Free_Work(ipTR)
         Call DGEMM_('N','T',
     &               nLambda,nInter,nLambda,
     &               1.0d0,Work(ipTRInv),nLambda,
     &               drdq(1,1,iIter),nInter,
     &               0.0d0,Work(ipTRT),nLambda)
         Call Free_Work(ipTRInv)
         Call FZero(rLambda(1,iIter),nLambda)
         Call DGEMM_('N','N',
     &               nLambda,1,nInter,
     &               1.0d0,Work(ipTRT),nLambda,
     &               dEdq(1,iIter),nInter,
     &               0.0d0,rLambda(1,iIter),nLambda)
         Call Free_Work(ipTRT)
*
         Call DScal_(nLambda,-One,rLambda(1,iIter),1) ! Sign
*                                                                      *
************************************************************************
*                                                                      *
*        Add contributions from constraints according to the last
*        iterations. See Eqn. 7
*
*        W^0 = H^0 - Sum(i) l_i d2rdq2(q_0)
*
         If (iIter.eq.nIter) Then
            Do iLambda = 1, nLambda
               Call DaXpY_(nInter**2,-rLambda(iLambda,nIter),
     &                    d2rdq2(1,1,iLambda),1,Hess,1)
            End Do
         End If
*                                                                      *
************************************************************************
*                                                                      *
*------- Transform various vectors and matrices to the nInter-nLambda
*        subspace. See Eqn. 10 and text following.
*
*------- Compute x
*
         If (nInter-nLambda.gt.0) Then
*
            Call FZero(x(1,iIter),nInter-nLambda)
            Call DGEMM_('T','N',
     &                  nInter-nLambda,1,nInter,
     &                  1.0d0,T(1,ipTti),nInter,
     &                  q(1,iIter),nInter,
     &                  0.0d0,x(1,iIter),nInter-nLambda)
*
*---------- Compute dx
*
            Call FZero(dx(1,iIter),nInter-nLambda)
            Call DGEMM_('T','N',
     &                  nInter-nLambda,1,nInter,
     &                  1.0d0,T(1,ipTti),nInter,
     &                  dq(1,iIter),nInter,
     &                  0.0d0,dx(1,iIter),nInter-nLambda)
*
*           Compute step restriction based on information from the
*           minimization in the x subspace. This restriction is based
*           on the step length in the x subspace.
*
            If (iIter.ne.nIter) Then
*
*              Compute dq step in the x subspace
*
               Call FZero(du,nInter)
               call dcopy_(nInter-nLambda,dx(1,iIter),1,du(1+nLambda),1)
               Call DGEMM_('N','N',
     &                     nInter,1,nInter,
     &                     One,T,nInter,
     &                     du,nInter,
     &                     Zero,Work(ipdq),nInter)
               dxdx=Sqrt(DDot_(nInter,Work(ipdq),1,Work(ipdq),1))
*              dxdx=Sqrt(DDot_(nInter-nLambda,dx(1,iIter),1,
*    &                                       dx(1,iIter),1))
*
               If (dxdx.lt.0.75D0*dxdx_last.and.
     &             dxdx.lt.(Beta-Thr)) Then
*                 Increase trust radius
                  xBeta=Min(One,xBeta*Sf)
C                 xBeta=xBeta*Sf
C
C
               Else If (dxdx.gt.1.25D0*dxdx_last.or.
     &                  dxdx.ge.(Beta+Thr)) Then
*                 Reduce trust radius
                  xBeta=Max(One/Five,xBeta/Sf)
               End If
               dxdx_last=dxdx
C              Write (6,*) 'dxdx=',dxdx
C              Write (6,*) 'xBeta=',xBeta
            End If
*
*---------- Compute the reduced gradient, observe again that we store
*           the forces rather than the gradients.
*
*           See text after Eqn. 13.
*
            Call Allocate_Work(ipTmp1,nInter)
            Call Allocate_Work(ipTmp3,nInter)
            Call FZero(Work(ipTmp3),nInter)
            Call Allocate_Work(ipTmp2,nInter*nLambda)
            Call FZero(Work(ipTmp2),nInter*nLambda)
            call dcopy_(nInter,dEdq(1,iIter),1,Work(ipTmp1),1)
*
            Call DGEMM_('N','N',
     &                  nInter,nLambda,nInter,
     &                  1.0d0,Hess,nInter,
     &                  T(1,ipTb),nInter,
     &                  0.0d0,Work(ipTmp2),nInter)
*
            Call DGEMM_('N','N',
     &                  nInter,1,nLambda,
     &                  1.0d0,Work(ipTmp2),nInter,
     &                  dy,nLambda,
     &                  0.0d0,Work(ipTmp3),nInter)
            Call Free_Work(ipTmp2)
!           Sign
            Call DaXpY_(nInter,-One,Work(ipTmp3),1,Work(ipTmp1),1)
*
            Call Free_Work(ipTmp3)
            Call FZero(dEdx(1,iIter),nInter-nLambda)
            Call DGEMM_('T','N',
     &                  nInter-nLambda,1,nInter,
     &                  1.0d0,T(1,ipTti),nInter,
     &                  Work(ipTmp1),nInter,
     &                  0.0d0,dEdx(1,iIter),nInter-nLambda)
            Call Free_Work(ipTmp1)
*
*           Compute step restriction based on information from the
*           minimization in the x subspace. This restriction is based
*           on the gradient in the x subspace.
*
            gg=Sqrt(DDot_(nInter-nLambda,dEdx(1,iIter),1,
     &                                  dEdx(1,iIter),1))
            If (gg.lt.0.75D0*gg_last.and.gg.lt.(Beta-Thr)) Then
*              Increase trust radius
               gBeta=Min(One,gBeta*Sf)
C              gBeta=gBeta*Sf
            Else If (gg.gt.1.25D0*gg_last.or.gg.ge.(Beta+Thr)) Then
*              Reduce trust radius
               gBeta=Max(One/Five,gBeta/Sf)
            End If
            gg_last=gg
C           Write (6,*) 'gg=',gg
C           Write (6,*) 'gBeta=',gBeta
*
         End If
*                                                                      *
************************************************************************
*                                                                      *
*------- Add contributions to the energy change
*
         If (iIter.eq.nIter) Then
*
*---------- Term due to constraint
*
            Call Allocate_Work(ipTmp,nLambda)
            Call FZero(Work(ipTmp),nLambda)
            Call Allocate_Work(ip_Tdy,nInter)
            Call FZero(Work(ip_Tdy),nInter)
            Call DGEMM_('N','N',
     &                  nInter,1,nLambda,
     &                  1.0d0,T(1,ipTb),nInter,
     &                  dy,nLambda,
     &                  0.0d0,Work(ip_Tdy),nInter)
            Call DGEMM_('T','N',
     &                  nLambda,1,nInter,
     &                  1.0d0,drdq(1,1,iIter),nInter,
     &                  Work(ip_Tdy),nInter,
     &                  0.0d0,Work(ipTmp),nLambda)
            Call Free_Work(ip_Tdy)
            Call DaXpY_(nLambda,One,r(1,iIter),1,Work(ipTmp),1)
            ed = ed - DDot_(nLambda,rLambda(1,iIter),1,Work(ipTmp),1)
            Call Free_Work(ipTmp)
*
*---------- Term due to coupling
*
            Call Allocate_Work(ip_Tr,nInter)
            Call FZero(Work(ip_Tr),nInter)
            Call DGEMM_('N','N',
     &                  nInter,1,nLambda,
     &                  1.0d0,T(1,ipTb),nInter,
     &                  dy,nLambda,
     &                  0.0d0,Work(ip_Tr),nInter)
            ed = ed + DDot_(nInter,Work(ip_Tr),1,dEdq,1) ! Sign
            Call Allocate_Work(ip_WTr,nInter)
            Call FZero(Work(ip_WTr),nInter)
            Call DGEMM_('N','N',
     &                  nInter,1,nInter,
     &                  1.0d0,Hess,nInter,
     &                  Work(ip_Tr),nInter,
     &                  0.0d0,Work(ip_WTr),nInter)
            ed = ed + Half * DDot_(nInter,Work(ip_Tr),1,Work(ip_WTr),1)
            Call Free_Work(ip_WTr)
            Call Free_Work(ip_Tr)
         End If
*                                                                      *
************************************************************************
*                                                                      *
*------- Restrict actual step, dy
*
*        Step in direction which fulfills the restriction is given
*        priority. This step is only reduced if it is larger than
*        half the overall step restriction.
*
*        Compute dq step in the y subspace
*
         Call FZero(du,nInter)
         call dcopy_(nLambda,dy,1,du,1)
         Call DGEMM_('N','N',
     &               nInter,1,nInter,
     &               One,T,nInter,
     &               du,nInter,
     &               Zero,Work(ipdq),nInter)
         dydy=Sqrt(DDot_(nInter,Work(ipdq),1,Work(ipdq),1))
*
*        Reduce y step size if larger than some maximum size
*        (the x step or half the total max step times a factor)
*
         dydymax=CnstWght*max(dxdx,Half*Beta)
         If (dydy.gt.dydymax) Then
*           Write (6,*) 'Reduce dydy!',dydy,' -> ',dydymax
            Call DScal_(nLambda,dydymax/dydy,dy,1)
            dydy=dydymax
         End If
*
*        The step reduction in the space which we minimize is such that
*        while the fulfillment of the constraint is not improving from
*        step to step we reduce the step length in the subspace in which
*        we do the minimization.

C        If (dydy.lt.0.75D0*dydy_last.or.dydy.lt.1.0D-5) Then
         If (dydy.lt.0.75D0*dydy_last.or.dydy.lt.1.0D-2) Then
*---------- Recent step in the space for the restriction is smaller than
*           the previous, or the recent step is smaller than some
*           threshold. Then increase step length for x, however not
*           more than the overall step restriction.
C           yBeta=Min(One,yBeta*Sf)
C           yBeta=yBeta*Sf
            yBeta=Min(Two,yBeta*Sf)
         Else If (dydy.gt.1.25D0*dydy_last.and.dydy.ge.1.0D-5) Then
*           Otherwise decrease step direction.
            yBeta=Max(One/Ten,yBeta/Sf)
         End If
*
*        Twist for MEP optimizations.
*
         If (iIter.eq.iOff_iter+1 .and. dydy.lt.1.0D-4
     &       .and. iIter.ne.1) Then
*           yBeta=Beta/Ten
            xBeta=xBeta*Half
         End If
C        Write (6,*) 'dydy_last=',dydy_last
C        Write (6,*) 'dydy=',dydy
         dydy_last=dydy
C        Write (6,*) 'yBeta=',yBeta
*
         If (jPrint.ge.99) Call RecPrt('dy(actual)',' ',dy,nLambda,1)
*                                                                      *
************************************************************************
*                                                                      *
      End Do
      Call GetMem('dqtmp','Free','Real',ipdq,nInter)
*                                                                      *
************************************************************************
*                                                                      *
*
      If (jPrint.ge.99) Then
         Call RecPrt('dEdx',' ',dEdx,nInter-nLambda,nIter)
         Call RecPrt('Lambda',' ',rLambda,nLambda,nIter)
         Call RecPrt('x',' ',x,nInter-nLambda,nIter)
         Call RecPrt('dx',' ',dx,nInter-nLambda,nIter)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the h vectors.   (Eqn. 16)
*     Observe the sign again since we store the force rather than the
*     gradient.
*
      call dcopy_(nInter*nIter,dEdq,1,dEdq_,1)
      If (jPrint.ge.99) Call RecPrt('dEdq',' ',dEdq,nInter,nIter)
      Do iIter = iOff_iter+1, nIter
         Do iLambda = 1, nLambda
            Call DaXpY_(nInter,rLambda(iLambda,nIter), ! Sign
     &                 drdq(1,iLambda,iIter),1,
     &                 dEdq_(1,iIter),1)
         End Do
      End Do
      If (jPrint.ge.99) Call RecPrt('dEdq_',' ',dEdq_,nInter,nIter)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the value of the Lagrangian
*
      Do iIter = iOff_iter+1, nIter
         Temp = Energy(iIter)
         Do iLambda = 1, nLambda
            Temp = Temp + rLambda(iLambda,iIter)*r(iLambda,iIter)
         End Do
         Energy(iIter)=Temp
      End Do
      If (jPrint.ge.99) Call RecPrt('Energy',' ',Energy,nIter,1)
*                                                                      *
************************************************************************
*                                                                      *
*---- Update the Hessian in the n subspace
*     There should be no negative eigenvalue, if so change the sign.
*
      If (jPrint.ge.6) Then
         Write (6,*)
         Write (6,*)
         Write (6,*) ' *** Updating the reduced Hessian ***'
         Write (6,*)
      End If
*
      If (iAnd(iOptC,4096).eq.4096) Then
*
*        If FINDTS option force minimization option during
*        the Hessian update.
*
         iOptC_Temp=128
      Else
         iOptC_Temp=iOptC
      End If
*
      Dummy = 0.0D0
      Call Update_H(nWndw,Hess,nInter,
     &              nIter,iOptC_Temp,Mode,ipMF,
     &              dq,dEdq_,iNeg,iOptH,HUpMet,nRowH,
     &              jPrint,Dummy,Dummy,nsAtom,IRC,.False.)

      If (jPrint.ge.99) Call RecPrt('Hess',' ',Hess,nInter,nInter)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the reduced Hessian
*
*     See text after Eqn. 13.
*
      If (nInter-nLambda.gt.0) Then
*
         Call Allocate_Work(ipTmp,(nInter-nLambda)*nInter)
         Call FZero(Work(ipTmp),(nInter-nLambda)*nInter)
         Call DGEMM_('T','N',
     &               nInter-nLambda,nInter,nInter,
     &               1.0d0,T(1,ipTti),nInter,
     &               Hess,nInter,
     &               0.0d0,Work(ipTmp),nInter-nLambda)
         Call FZero(W,(nInter-nLambda)**2)
         Call DGEMM_('N','N',
     &               nInter-nLambda,nInter-nLambda,nInter,
     &               1.0d0,Work(ipTmp),nInter-nLambda,
     &               T(1,ipTti),nInter,
     &               0.0d0,W,nInter-nLambda)
         Call Free_Work(ipTmp)
         If (jPrint.ge.99)
     &      Call RecPrt('W',' ',W,nInter-nLambda,nInter-nLambda)
*                                                                      *
************************************************************************
*                                                                      *
*----    Update dx
*
C        Write (6,*) 'yBeta=',yBeta
C        Write (6,*) 'gBeta=',gBeta
C        Write (6,*) 'xBeta=',xBeta
         tBeta= Max(Beta*yBeta*Min(xBeta,gBeta),Beta/Ten)
C        Write (6,*) 'tBeta=',tBeta
         Call Newq(x,nInter-nLambda,nIter,dx,W,dEdx,Err,EMx,
     &             RHS,iPvt,dg,A,nA,ed,iOptC,tBeta,
     &             nFix,ip,UpMeth,Energy,Line_Search,Step_Trunc,
     &             Restriction)
         GNrm=
     &    Sqrt(DDot_(nInter-nLambda,dEdx(1,nIter),1,dEdx(1,nIter),1))
*
      Else
*        Negative norm should serve as sign that there is no gradient
         GNrm=-One
      End If
*
*     A small trick here to avoid a false convergence.
*
c     If (GNrm.lt.2.0D-4.or.nIter.eq.iOff_Iter+1) Then
c        Do iLambda = 1, nLambda
c           GNrm=GNrm+Sqrt((rLambda(iLambda,nIter)*
c    &                            r(iLambda,nIter))**2)
c        End Do
c        GNrm=Min(1.0D-2,GNrm)
c        If (nIter.eq.iOff_Iter+1) GNrm=1.0D-2
c     End If
      If (jPrint.ge.99) Call RecPrt('dx',' ',dx,nInter-nLambda,nIter)
*                                                                      *
************************************************************************
*                                                                      *
*---- Back transform dy and dx to dq.
*
*     See Eqn. 10.
*
*temporary code
*
*     dy only, constraint
*
C     Call FZero(du,nInter)
C     call dcopy_(nLambda,dy,1,du,1)
C     Call DGEMM_('N','N',
C    &            nInter,1,nInter,
C    &            One,T,nInter,
C    &            du,nInter,
C    &            Zero,dq(1,nIter),nInter)
C     dydy=Sqrt(DDot_(nInter,dq,1,dq,1))
C     Call RecPrt('dq(dy)',' ',dq,nInter,nIter)
C     Write (6,*) '<R(q_0)|dy>=',DDot_(nInter,
C    &              dRdq(1,1,nIter),1,dq(1,nIter),1)
*
*     dx only, constraint minimization
*
C     Call FZero(du,nInter)
C     call dcopy_(nInter-nLambda,dx(1,nIter),1,du(1+nLambda),1)
C     Call DGEMM_('N','N',
C    &            nInter,1,nInter,
C    &            One,T,nInter,
C    &            du,nInter,
C    &            Zero,dq(1,nIter),nInter)
C     dxdx=Sqrt(DDot_(nInter,dq,1,dq,1))
C     Call RecPrt('dq(dx)',' ',dq,nInter,nIter)
C     Write (6,*) '<R(q_0)|dx>=',DDot_(nInter,
C    &              dRdq(1,1,nIter),1,dq(1,nIter),1)
*temporary code
*
*     dy+dx, full step
*
      Call FZero(du,nInter)
      call dcopy_(nLambda,dy,1,du,1)
      call dcopy_(nInter-nLambda,dx(1,nIter),1,du(1+nLambda),1)
      Call DGEMM_('N','N',
     &            nInter,1,nInter,
     &            One,T,nInter,
     &            du,nInter,
     &            Zero,dq(1,nIter),nInter)
*
*     Compute q for the next iteration
*
      call dcopy_(nInter,q(1,nIter),1,q(1,nIter+1),1)
      Call DaXpY_(nInter,One,dq(1,nIter),1,q(1,nIter+1),1)
      If (jPrint.ge.99) Call RecPrt('q',' ',q,nInter,nIter+1)
*                                                                      *
************************************************************************
*                                                                      *
*     StpMax from q
*
      Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,nInter,
     &            dEdq_(1,nIter),dq(1,nIter),Lbl)
*
*     GrdMax for dEdx
*
      If (nInter-nLambda.gt.0) Then
*
         StpMax_Save=StpMax
         StpLbl_Save=StpLbl
         Call Allocate_Work(ipLbl,nInter-nLambda)
         Call Char2Real(Lbl,Work(ipLbl),(nInter-nLambda)*8)
         GrdMax=Zero
         Do i = 1, nInter-nLambda
            Write (Lbl(i),'(A,I3.3)') 'dEdx',i
         End Do
         Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,nInter-nLambda,
     &               dEdx(1,nIter),dx(1,nIter),Lbl)
         StpMax=StpMax_Save
         StpLbl=StpLbl_Save
         Call Real2Char(Work(ipLbl),Lbl,(nInter-nLambda)*8)
         Call Free_Work(ipLbl)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
