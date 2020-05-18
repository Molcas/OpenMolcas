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
* Copyright (C) 2003, 2020 Roland Lindh                                *
************************************************************************
      Subroutine Con_Opt(r,drdq,T,dEdq,rLambda,q,dq,dy,dx,dEdq_,du,x,
     &                   dEdx,W,GNrm,nWndw,
     &                   Hess,nInter,nIter,iOptC,Mode,MF,
     &                   iOptH,HUpMet,jPrint_,Energy,nLambda,
     &                   mIter,nRowH,Err,EMx,RHS,iPvt,dg,A,nA,ed,
     &                   Beta,Beta_Disp,nFix,iP,UpMeth,
     &                   Line_Search,Step_Trunc,Lbl,GrdLbl,StpLbl,
     &                   GrdMax,StpMax,d2rdq2,nsAtom,IRC,CnstWght,
     &                   iOpt_RS,Thr_RS,iter_)
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
*             Ref:                                                     *
*             "A Reduced-Restricted-Quasi-Newton-Raphson Method for    *
*              Locating and Optimizing Energy Crossing Points Between  *
*              Two Potential Energy Surfaces",                         *
*             J. M. Anglada and J. M. Bofill,                          *
*             J. Comput. Chem., 18, 992-1003 (1997)                    *
*                                                                      *
*             For more details consult:                                *
*             "New General Tools for Constrained Geometry Optimization"*
*             L. De Vico, M. Olivucci, and R. Lindh,                   *
*             JCTC, 1:1029-1037, 2005                                  *
*             DOI: 10.1021/ct500949                                    *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemical Physics,                 *
*             University of Lund, SWEDEN                               *
*             July '03                                                 *
************************************************************************
      use RDC
      Implicit Real*8 (a-h,o-z)
      External Restriction_Step, Restriction_Disp_Con
      Real*8 Restriction_Step, Restriction_Disp_Con
#include "real.fh"
#include "stdalloc.fh"
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
     &       dg(mIter), A(nA), d2rdq2(nInter,nInter,nLambda),
     &       MF(3*nsAtom)
      Integer iPvt(nInter+1), iP(nInter), iNeg(2)
      Logical Line_Search, Found, IRC_setup, Corrected
      Character HUpMet*6, UpMeth*6, Step_Trunc*1, Lbl(nInter+nLambda)*8,
     &          GrdLbl*8, StpLbl*8, StpLbl_Save*8
      Real*8, Allocatable:: dq_xy(:), Trans(:), Tmp1(:), Tmp2(:,:)
      Real*8, Allocatable:: RT(:,:), RTInv(:,:), RRR(:,:), RRInv(:,:),
     &                      RR(:,:), Tdy(:), Tr(:), WTr(:),
     &                      Hessian(:,:)
*                                                                      *
************************************************************************
*                                                                      *
      jPrint=jPrint_
*#define _DEBUG_
#ifdef _DEBUG_
      Write (6,*)
      Write (6,*) '****************************************************'
      Write (6,*) '**** Con_Opt: Input data ***************************'
      Write (6,*) '****************************************************'
      Write (6,*)
      Write (6,*) 'iOpt_RS=',iOpt_RS
      Write (6,*)
      Call RecPrt('Con_Opt: Energy',' ',Energy,nIter,1)
      Call RecPrt('Con_Opt: q',' ',q,nInter,nIter)
      Call RecPrt('Con_Opt: dEdq',' ',dEdq,nInter,nIter)
      Call RecPrt('Con_Opt: Hess(in)',' ',Hess,nInter,nInter)
      Call RecPrt('Con_Opt: r',' ',r,nLambda,nIter)
      Do iIter = 1, nIter
         Write (6,*)' iIter=',iIter
         Call RecPrt('Con_Opt: drdq(orig)',' ',drdq(1,1,iIter),
     &               nInter,nLambda)
      End Do
      Do iLambda = 1, nLambda
         Call RecPrt('Con_Opt: d2rdq2(iLambda)',' ',
     &                            d2rdq2(1,1,iLambda),nInter,nInter)
      End Do
#endif
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
      Beta_Disp_Min=1.0D-10
      Call Get_iScalar('iOff_Iter',iOff_Iter)
      Call mma_allocate(dq_xy,nInter,Label='dq_xy')
      Call mma_allocate(Hessian,nInter,nInter,Label='Hessian')
      If (iOpt_RS.ne.0) Then
         nInter_=nInter
         nLambda_=nLambda
         Call mma_allocate(q_,nInter,Label='q_')
         Call mma_allocate(T_,nInter,nInter,Label='T_')
         Call mma_allocate(dy_,nLambda,Label='dy_')
      End If
*
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      iSt=Max(iOff_Iter+1,nIter-nWndw+1)
      Do iIter = iSt, nIter
*     Do iIter = iOff_Iter+1, nIter
#ifdef _DEBUG_
         Write (6,*)
         Write (6,*) '>>>>>> iIter=',iIter
         Write (6,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*------- Compute the lambdas values (Eqn. 17)
*
*        dEdq=drdq l ; from  h(q_0,l)=dEdq - drdq l = 0
*
*        drdq^T dEdq= (drdq^T drdq) l
*
*        l =  (drdq^T drdq)^{-1} drdq^T dEdq
*
         Call mma_allocate(RRR,nLambda,nInter,Label='RRR')
         Call mma_allocate(RRInv,nLambda,nLambda,Label='RRInv')
         Call mma_allocate(RR,nLambda,nLambda,Label='RR')
*
*        drdq^T drdq
*
         RR(:,:)=Zero
         Call DGEMM_('T','N',nLambda,nLambda,nInter,
     &               1.0d0,drdq(1,1,iIter),nInter,
     &                     drdq(1,1,iIter),nInter,
     &               0.0d0,RR,nLambda)
*
*        (drdq^T drdq)^{-1}
*
         Call MInv(RR,RRInv,iSing,Det,nLambda)
         Call mma_deallocate(RR)
*
*        (drdq^T drdq)^{-1} drdq^T
*
         RRR(:,:)=Zero
         Call DGEMM_('N','T',nLambda,nInter,nLambda,
     &               1.0d0,RRInv,nLambda,
     &                     drdq(1,1,iIter),nInter,
     &               0.0d0,RRR,nLambda)
         Call mma_deallocate(RRInv)
*
*        l = (T_b^T drdq)^{-1} drdq^T dEdq
*
*        Note the sign conflict due to dEdq stored as a force.
*
         rLambda(:,iIter)=Zero
         Call DGEMM_('N','N',nLambda,1,nInter,
     &              -1.0d0,RRR,nLambda,    ! Sign conflict
     &                     dEdq(1,iIter),nInter,
     &               0.0d0,rLambda(1,iIter),nLambda)
         Call mma_deallocate(RRR)
#ifdef _DEBUG_
         Call RecPrt('rLambda(iIter)',' ',rLambda(1,iIter),nLambda,1)
#endif
      End Do
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iIter = iSt, nIter
*     Do iIter = iOff_Iter+1, nIter
#ifdef _DEBUG_
         Write (6,*)
         Write (6,*) '>>>>>> iIter=',iIter
         Write (6,*)
#endif
*                                                                      *
************************************************************************
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
            RR_=Sqrt(DDot_(nInter,drdq(1,iLambda,iIter),1,
     &                          drdq(1,iLambda,iIter),1))
*
*           Make sure that we dont mess up gradients which are zero vectors.
*
            If ( RR_.lt.1.0D-12 ) Then
*           If ( RR_.lt.1.0D-12 .and.
*    &           Abs(r(iLambda,iIter)).gt.1.0D-12 ) Then
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
               If (IRC_SetUp) Call ReacQ(MF,3*nsAtom,dEdq(1,iIter),
     &                                   nInter)
*
*              Try to use the transverse direction if the gradient is
*              too small
*
               RR_=Sqrt(DDot_(nInter,dEdq(1,iIter),1,dEdq(1,iIter),1))
               If (RR_.lt.1.0D-12) Then
                  Call qpg_dArray('Transverse',Found,nTrans)
                  If (Found.and.(nTrans.eq.3*nsAtom)) Then
                     Call mma_Allocate(Trans,3*nsAtom,Label='Trans')
                     Call Get_dArray('Transverse',Trans,nTrans)
                     Call ReacQ(Trans,3*nsAtom,dEdq(1,iIter),nInter)
                     RR_=Sqrt(DDot_(nInter,dEdq(1,iIter),1,
     &                                   dEdq(1,iIter),1))
                     Call DScal_(nInter,One/RR_,dEdq(1,iIter),1)
                     Call mma_deallocate(Trans)
                  End If
               Else
                  Call DScal_(nInter,One/RR_,dEdq(1,iIter),1)
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
#ifdef _DEBUG_
               Call RecPrt('Con_Opt: drdq(1)',' ',drdq,nInter,
     &                                            nLambda*nIter)
#endif
*
*------------- Orthogonalize against the gradient and all other
*              constraints
*
               RG=DDot_(nInter,drdq(1,iLambda,iIter),1,
     &                        dEdq(1,        iIter),1)
               Call DaXpY_(nInter,-RG,dEdq(1,        iIter),1,
     &                                drdq(1,iLambda,iIter),1)
               RR_=Sqrt(DDot_(nInter,drdq(1,iLambda,iIter),1,
     &                             drdq(1,iLambda,iIter),1))
C              Call DScal_(nInter,One/RR_,drdq(1,iLambda,iIter),1)
               Do jLambda = 1, nLambda
                  If (jLambda.ne.iLambda) Then
                     RG=DDot_(nInter,drdq(1,iLambda,iIter),1,
     &                              drdq(1,jLambda,iIter),1)
                     Call DaXpY_(nInter,-RG,drdq(1,jLambda,iIter),1,
     &                                     drdq(1,iLambda,iIter),1)
                     RR_=Sqrt(DDot_(nInter,drdq(1,iLambda,iIter),1,
     &                                   drdq(1,iLambda,iIter),1))
                     Call DScal_(nInter,One/RR_,drdq(1,iLambda,iIter),1)
                  End If
               End Do
*
               Continue
#ifdef _DEBUG_
               Call RecPrt('Con_Opt; drdq(2)',' ',drdq,nInter,
     &                                            nLambda*nIter)
#endif
            End If
         End Do
************************************************************************
*                                                                      *
*        NOTE: for historical reasons the code stores the force rather *
*              than the gradient. Hence, be careful of the sign when   *
*              ever the term dEdq, dEdq_h show up.                     *
*                                                                      *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*------- Set up the T-matrix which separates the space into the
*        subspace in which the constraint is accomplished and
*        the complemental subspace.
*
*        [T_b,T_{ti}]
*
*        With the properties
*
*        drdq^T T_b =/= 0 (but not = 1, see note below)
*
*        drdq^T T_{ti} = T_b^T T_{ti} = 0
*
*        T_{ti}^T T_{ti} = 1
*
         Call GS(drdq(1,1,iIter),nLambda,T,nInter,.False.,.False.)
         If (iOpt_RS.ne.0) T_(:,:)=T(:,:)
#ifdef _DEBUG_
         Call RecPrt('Con_Opt: T-Matrix',' ',T,nInter,nInter)
         Call RecPrt('Con_Opt: T_b',' ',T(1,ipTB),nInter,nLambda)
         Call RecPrt('Con_Opt: T_ti',' ',T(1,ipTti),nInter,
     &                                              nInter-nLambda)
#endif
*                                                                      *
*       Note that the paper by Anglada and Bofill has some errors on   *
*       the properties of T. Especially with respect to the properties *
*       of T_b.                                                        *
*                                                                      *
************************************************************************
*                                                                      *
*------- Compute delta y  (Eqn. 11) -- in case of full relaxation --
*
*        r(q)+drdq^T dq=0
*
*        dq = [T_b, T_{ti}] (dy,dx)^T
*
*        r(q) = - drdq^T T_b dy - drdq^T T_{ti} dx
*
*        r(q) = - drdq^T T_b dy
*
*        dy = - (drdq^T T_b)^{-1} r(q)
*
         Call mma_Allocate(RTInv,nLambda,nLambda,Label='RTInv')
         Call mma_Allocate(RT   ,nLambda,nLambda,Label='RT   ')
         RT(:,:) = Zero
*
*        drdq^T T_b
*
         Call DGEMM_('T','N',nLambda,nLambda,nInter,
     &               1.0d0,drdq(1,1,iIter),nInter,
     &                     T(1,ipTb),nInter,
     &               0.0d0,RT,nLambda)
#ifdef _DEBUG_
         Call RecPrt('Con_Opt: RT',' ',RT,nLambda,nLambda)
#endif
*
         Call MInv(RT,RTInv,iSing,Det,nLambda)
         Call mma_deallocate(RT)
#ifdef _DEBUG_
         Call RecPrt('Con_Opt: RTInv',' ',RTInv,nLambda,nLambda)
#endif
*
*        dy = - (drdq^T T_b)^{-1} r(q)
*
         dy(:)=Zero
         Call DGEMM_('N','N',nLambda,1,nLambda,
     &              -1.0d0,RTInv,nLambda,
     &                     r(1,iIter),nLambda,
     &               0.0d0,dy,nLambda)
         Call mma_deallocate(RTInv)
#ifdef _DEBUG_
         Call RecPrt('Con_Opt: dy(full)',' ',dy,nLambda,1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        Add contributions from constraints according to the last
*        iterations. See Eqn. 7,
*
*        The Hessian of the Lagranian expressed in the full space.
*
*        W^0 = H^0 - Sum(i) l_{0,i} d2rdq2(q_0)
         If (iOpt_RS.eq.0) Then
            Hessian(:,:) = Hess(:,:)
            If (iIter.eq.nIter) Then
               Do iLambda = 1, nLambda
                  Call DaXpY_(nInter**2,-rLambda(iLambda,nIter),
     &                                  d2rdq2(1,1,iLambda),1,Hessian,1)
               End Do
            End If
         Else
            Call Hessian_Kriging_Layer(q(1,iIter),Hessian,nInter)
            Do iLambda = 1, nLambda
               Call DaXpY_(nInter**2,-rLambda(iLambda,iIter),
     &                               d2rdq2(1,1,iLambda),1,Hessian,1)
            End Do
         End If
#ifdef _DEBUG_
         Call RecPrt('W',' ',Hessian,nInter,nInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*------- Transform various vectors and matrices to the nInter-nLambda
*        subspace. See Eqn. 10 and text following.
*
*------- Compute x, the coordinates of the reduced space.
*
         If (nInter-nLambda.gt.0) Then
*
*           q = T_b y + T_{ti} x = T [y,x]^T
*
*           x = T_{ti}^T q, since T_{ti}^T T_b = 0
*
            x(:,iIter)=Zero
            Call DGEMM_('T','N',nInter-nLambda,1,nInter,
     &                  1.0d0,T(1,ipTti),nInter,
     &                        q(1,iIter),nInter,
     &                  0.0d0,x(1,iIter),nInter-nLambda)
*
*---------- Compute dx
*
*           dx = T_{ti}^T dq
*
            dx(:,iIter)=Zero
            Call DGEMM_('T','N',nInter-nLambda,1,nInter,
     &                  1.0d0,T(1,ipTti),nInter,
     &                        dq(1,iIter),nInter,
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
*              du = [0,dx]^T
*
*              dq = T [0,dx]^T
*
               du(:)=Zero
               call dcopy_(nInter-nLambda,dx(1,iIter),1,du(1+nLambda),1)
               dq_xy(:)=Zero
               Call DGEMM_('N','N',nInter,1,nInter,
     &                     One,T,nInter,
     &                         du,nInter,
     &                     Zero,dq_xy,nInter)
               dxdx=Sqrt(DDot_(nInter,dq_xy,1,dq_xy,1))
*              dxdx=Sqrt(DDot_(nInter-nLambda,dx(1,iIter),1,
*    &                                        dx(1,iIter),1))
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
*                                                                      *
************************************************************************
*                                                                      *
*---------- Compute the reduced gradient, observe again that we store
*           the forces rather than the gradients.
*
*           See text after Eqn. 13.
*
*           dEdx = T_{ti}^T (dEdq + W_{ex}T_b dy)
*
            Call mma_allocate(Tmp1,nInter,Label='Tmp1')
            Call mma_allocate(Tmp2,nInter,nLambda,Label='Tmp2')
*
#ifdef _DEBUG_
            Call RecPrt('Con_Opt: dEdq',' ',dEdq(1,iIter),1,nInter)
            Call RecPrt('Con_Opt: W',' ',Hessian,nInter,nInter)
            Call RecPrt('Con_Opt: T',' ',T,nInter,nInter)
            Write (6,*) 'ipTb,ipTti=',ipTb,ipTti
            Call RecPrt('Con_Opt: dy',' ',dy,1,nLambda)
#endif
*
*           W_{ex} T_b
*
            Tmp2(:,:)=Zero
            Call DGEMM_('N','N',nInter,nLambda,nInter,
     &                  1.0d0,Hessian,nInter,
     &                        T(1,ipTb),nInter,
     &                  0.0d0,Tmp2,nInter)
#ifdef _DEBUG_
            Call RecPrt('W_{ex} T_b',' ',Tmp2,nInter,nLambda)
#endif
*
*           dEdq + W_{ex} T_b dy
*
*           Since we are computing the force we have the conflicting
*           sign below.
*
            Tmp1(:)=dEdq(:,iIter)
            Call DGEMM_('N','N',nInter,1,nLambda,
     &                 -1.0d0,Tmp2,nInter,   ! Sign conflict
     &                        dy,nLambda,
     &                  1.0d0,Tmp1,nInter)
            Call mma_deallocate(Tmp2)
#ifdef _DEBUG_
            Call RecPrt('dEdq + W_{ex} T_b dy',' ',Tmp1,1,nInter)
#endif
*
*           dEdx = T^t_{ti} (dEdq + W_{ex} T_b dy)
*
            dEdx(:,iIter)=Zero
            Call DGEMM_('T','N',nInter-nLambda,1,nInter,
     &                  1.0d0,T(1,ipTti),nInter,
     &                        Tmp1,nInter,
     &                  0.0d0,dEdx(1,iIter),nInter-nLambda)
#ifdef _DEBUG_
            Write (6,*) 'iIter=',iIter
            Call RecPrt('dEdx(1,iIter)',' ',dEdx(1,iIter),1,
     &                  nInter-nLambda)
#endif
            Call mma_deallocate(Tmp1)
*
*           Compute step restriction based on information from the
*           minimization in the x subspace. This restriction is based
*           on the gradient in the x subspace.
*
            gg=Sqrt(DDot_(nInter-nLambda,dEdx(1,iIter),1,
     &                                   dEdx(1,iIter),1))
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
************************************************************************
*                                                                      *
*------- Add contributions to the Lagranian energy change
*
         If (iIter.eq.nIter) Then
*
*---------- Term due to constraint
*
            Call mma_allocate(Tmp1,nLambda,Label='Tmp1')
            Call mma_allocate(Tdy,nInter,Label='Tdy')
*
*           T_b dy
*
            Tdy(:)=Zero
            Call DGEMM_('N','N',nInter,1,nLambda,
     &                  1.0d0,T(1,ipTb),nInter,
     &                        dy,nLambda,
     &                  0.0d0,Tdy,nInter)
*
*           drdq^T dq = drdq^T T_b dy
*
            Tmp1(:)=Zero
            Call DGEMM_('T','N',nLambda,1,nInter,
     &                  1.0d0,drdq(1,1,iIter),nInter,
     &                        Tdy,nInter,
     &                  0.0d0,Tmp1,nLambda)
            Call mma_deallocate(Tdy)
*
*           r + drdq^T dq
*
            Call DaXpY_(nLambda,One,r(1,iIter),1,Tmp1,1)
*
*           l^T (r + drdq^T dq)
*
            ed = ed - DDot_(nLambda,rLambda(1,iIter),1,Tmp1,1)
            Call mma_deallocate(Tmp1)
*
*---------- Term due to coupling
*
            Call mma_allocate(Tr,nInter,Label='Tr')
*
*           T_b dy
*
            Tr(:)=Zero
            Call DGEMM_('N','N',nInter,1,nLambda,
     &                  1.0d0,T(1,ipTb),nInter,
     &                        dy,nLambda,
     &                  0.0d0,Tr,nInter)
*
*           dEdq^T T_b dy
*
*           Note the sign conflict.
*
            ed = ed - DDot_(nInter,Tr,1,dEdq,1)
*
            Call mma_allocate(WTr,nInter,Label='WTr')
*
*           W T_b dy
*
            WTr(:)=Zero
            Call DGEMM_('N','N',nInter,1,nInter,
     &                  1.0d0,Hessian,nInter,
     &                        Tr,nInter,
     &                  0.0d0,WTr,nInter)
*
*            dy^T T_b^T W T_b dy
            ed = ed + Half * DDot_(nInter,Tr,1,WTr,1)
            Call mma_deallocate(WTr)
            Call mma_deallocate(Tr)
         End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*------- Restrict actual step, dy
*
*        Step in direction which fulfills the restriction is given
*        priority. This step is only reduced if it is larger than
*        half the overall step restriction.
*
         If (iOpt_RS.eq.0) Then
*
*           Compute dq step in the y subspace
*
*           dq_y = T [dy, 0]^T

            du(:)=Zero
            du(1:nLambda)=dy(:)
            dq_xy(:)=Zero
            Call DGEMM_('N','N',nInter,1,nInter,
     &                  One,T,nInter,
     &                      du,nInter,
     &                  Zero,dq_xy,nInter)
*
            dydy=Sqrt(DDot_(nInter,dq_xy,1,dq_xy,1))
*
*
*           Reduce y step size if larger than some maximum size
*           (the x step or half the total max step times a factor)
*
            dydymax=CnstWght*max(dxdx,Half*Beta)
            If (dydy.gt.dydymax) Then
*              Write (6,*) 'Reduce dydy!',dydy,' -> ',dydymax
               Call DScal_(nLambda,dydymax/dydy,dy,1)
               dydy=dydymax
               Step_Trunc='*'
            Else
*              Write (6,*) 'No reduce dydy!',dydy,' < ',dydymax
               Step_Trunc=' '
            End If
*
*           The step reduction in the space which we minimize is such
*           that while the fulfillment of the constraint is not
*           improving from step to step we reduce the step length in
*           the subspace in which we do the minimization.

            If (dydy.lt.0.75D0*dydy_last.or.dydy.lt.1.0D-2) Then
*------------- Recent step in the space for the restriction is smaller
*              than the previous, or the recent step is smaller than
*              some threshold. Then increase step length for x, however
*              not more than the overall step restriction.
               yBeta=Min(Two,yBeta*Sf)
            Else If (dydy.gt.1.25D0*dydy_last.and.dydy.ge.1.0D-5) Then
*              Otherwise decrease step direction.
               yBeta=Max(One/Ten,yBeta/Sf)
            End If
         Else
*
*           Here in the case of kriging and restricted-variance
*           optimization.
*
            iCount=1
            iCount_Max=100
            If (Step_Trunc.eq.'N') Step_Trunc=' '
            Fact=One
*
*           We only need this for the last point.
*
            If (iIter.ne.nIter) Go to 667
*
            tmp=0.0D0
            Do i = 1, nLambda
               Do j = 1, nInter
                  tmp = Max(tmp,Abs(drdq(j,i,iIter)))
               End Do
            End Do
            tmp=Min(tmp,0.30D0) ! Some constraints can have huge
                                ! gradients. So be a bit careful.
            Beta_Disp_=Max(Beta_Disp_Min,tmp*Half*Beta_Disp)
*
            q_(:)=q(:,iIter)
            du(1:nInter-nLambda)=Zero ! Fake dx(:)=Zero
*
#ifdef _DEBUG_
            Write (6,*) 'Step_trunc=',Step_trunc
            Write (6,*) 'Beta_Disp_=',Beta_Disp_
            Write (6,*) 'Start: dy(:)=',dy(:)
#endif
*
            If (DDot_(nLambda,dy,1,dy,1).lt.1.0D-12) Go To 667
*
 666        Continue
            du(1:nLambda)=(One/Fact)*dy(:)
            du(1+nLambda:nInter)=Zero
            Call RecPrt('du(1)',' ',du,1,nInter)
            Call Backtrans_K(T,du,dq_xy,nInter,1)
            q(:,iIter+1)=q(:,iIter)+dq_xy(:)
*
            Call Dispersion_Kriging_Layer(q(1,iIter+1),dydy,nInter)
            Write (6,*) 'dydy(1)=',dydy
*
            du(1:nInter-nLambda)=Zero ! Fake dx(:)=Zero
            dy_(:)=(One/Fact)*dy(:)
*
            dydy=Restriction_Disp_Con(x(1,iIter),du,nInter-nLambda)
            Write (6,*) 'dydy(2)=',dydy
            If (iCount.eq.1) Then
               Fact_long=Fact
               dydy_long=dydy
               Fact_short=Zero
               dydy_short=dydy_long+One
            End If
#ifdef _DEBUG_
            Write (6,*) 'dydy,Fact,iCount=', dydy,Fact,iCount
#endif
            If (dydy.gt.Beta_Disp_ .or. iCount.gt.1) Then
               If (Abs(Beta_Disp_-dydy).lt.Thr_RS) Go To 667
               iCount=iCount+1
               If (iCount.gt.iCount_Max) Then
                  Write (6,*) 'iCount.gt.iCount_Max'
                  Call Abend()
               End If
               Call Find_RFO_Root(Fact_long,dydy_long,
     &                            Fact_short,dydy_short,
     &                            Fact,dydy,Beta_Disp_)
               Step_Trunc='*'
               Go To 666
            End If
 667        Continue
            dy(:)=(One/Fact)*dy(:)
#ifdef _DEBUG_
            Write (6,*) 'Step_trunc=',Step_trunc
            Write (6,*) 'Final: dy(:)=',dy(:)
#endif
            Write (6,*) 'Done dy!'
*
         End If
*
*        Twist for MEP optimizations.
*
         If (iIter.eq.iOff_iter+1 .and. dydy.lt.1.0D-4
     &       .and. iIter.ne.1) xBeta=xBeta*Half
         dydy_last=dydy
*
#ifdef _DEBUG_
         Call RecPrt('Con_Opt: dy(actual)',' ',dy,nLambda,1)
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End Do ! Do iIter = iOff_Iter+1, nIter
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
#ifdef _DEBUG_
      Call RecPrt('Con_Opt: dEdx',' ',dEdx,nInter-nLambda,nIter)
      Call RecPrt('Con_Opt: Lambda',' ',rLambda,nLambda,nIter)
      Call RecPrt('Con_Opt: x',' ',x,nInter-nLambda,nIter)
      Call RecPrt('Con_Opt: dx',' ',dx,nInter-nLambda,nIter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the h vectors.   (Eqn. 16)
*
*     h(q,l) = dEdq - drdq l_0^T
*
*     Note the sign conflict due to storage of the force rather than the
*     gradient.
#ifdef _DEBUG_
      Call RecPrt('Con_Opt: dEdq',' ',dEdq,nInter,nIter)
#endif
      dEdq_(:,:) = dEdq(:,:)
      Do iIter = iOff_iter+1, nIter
         Do iLambda = 1, nLambda
            Call DaXpY_(nInter,rLambda(iLambda,nIter),   ! Sign conflict
     &                          drdq(1,iLambda,iIter),1,
     &                          dEdq_(1,iIter),1)
         End Do
      End Do
#ifdef _DEBUG_
      Call RecPrt('Con_Opt: dEdq_',' ',dEdq_,nInter,nIter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the value of the Lagrangian
*
*     L = E + l r(q)
*
      Do iIter = iOff_iter+1, nIter
         Temp = Energy(iIter)
         Do iLambda = 1, nLambda
            Temp = Temp + rLambda(iLambda,iIter)*r(iLambda,iIter)
         End Do
         Energy(iIter)=Temp
      End Do
#ifdef _DEBUG_
      Call RecPrt('Con_Opt: Lagrangian',' ',Energy,nIter,1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Update the Hessian in the m subspace
*     There should be no negative eigenvalue, if so change the sign.
*
#ifdef _DEBUG_
      Write (6,*)
      Write (6,*)
      Write (6,*) ' *** Updating the reduced Hessian ***'
      Write (6,*)
#endif
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
#ifdef _DEBUG_
      Call RecPrt('Con_Opt: Hessian(raw)',' ',Hessian,nInter,nInter)
      Write (6,*) 'iOptH=',iOptH
#endif
      If (Step_Trunc.eq.'N') Step_Trunc=' '
      Call Update_H(nWndw,Hessian,nInter,
     &              nIter,iOptC_Temp,Mode,MF,
     &              dq,dEdq_,iNeg,iOptH,HUpMet,nRowH,
     &              jPrint,Dummy,Dummy,nsAtom,IRC,.False.,Corrected)
      If (Corrected) Step_Trunc='#'

#ifdef _DEBUG_
      Call RecPrt('Con_Opt: Hessian(updated)',' ',Hessian,nInter,nInter)
      Write (6,*) 'Step_Trunc=',Step_trunc
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the reduced Hessian
*
*     See text after Eqn. 13.
*
*     T_{ti}^T W T_{ti}
*
      If (nInter-nLambda.gt.0) Then
*
         Call mma_allocate(Tmp2,nInter-nLambda,nInter,Label='Tmp2')
         Tmp2(:,:)=Zero
         Call DGEMM_('T','N',nInter-nLambda,nInter,nInter,
     &               1.0d0,T(1,ipTti),nInter,
     &                     Hessian,nInter,
     &               0.0d0,Tmp2,nInter-nLambda)
         W(:,:)=Zero
         Call DGEMM_('N','N',nInter-nLambda,nInter-nLambda,nInter,
     &               1.0d0,Tmp2,nInter-nLambda,
     &                     T(1,ipTti),nInter,
     &               0.0d0,W,nInter-nLambda)
         Call mma_deallocate(Tmp2)
#ifdef _DEBUG_
         Call RecPrt('Con_Opt: W',' ',W,nInter-nLambda,nInter-nLambda)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*----    Update dx
*
*        Set threshold depending on if restriction is w.r.t. step size
*        or variance.
*
         fact=One
         If (iOpt_RS.eq.0) Then
            Beta_Disp_= 1.0D0 ! Dummy assign
         Else
*
*           Note that we use the dEdx data for the last point on the
*           real PES.
*
#ifdef _DEBUG_
            Write (6,*) 'Beta_Disp_=',Beta_Disp_
#endif
            Beta_Disp_Min=1.0D-10
            tmp=0.0D0
            Do i = 1, nInter-nLambda
               tmp = Max(tmp,Abs(dEdx(i,iter_)))
            End Do
            Beta_Disp_=Max(Beta_Disp_,Beta_Disp_Min,tmp*Beta_Disp)
#ifdef _DEBUG_
            Write (6,*) 'tmp,Beta_Disp_=',tmp,Beta_Disp_
#endif
*
*           Copy stuff so that the restricted_disp_Cons routine
*           can figure out how to compute the dispersion.
*
            q_(:)=q(:,nIter)
            dy_(:)=dy(:)
         End If
         GNrm=
     &    Sqrt(DDot_(nInter-nLambda,dEdx(1,iter_),1,dEdx(1,iter_),1))
         tBeta= Min(1.0D3*GNrm,Beta)
         Thr_RS=1.0D-7
#ifdef _DEBUG_
            Write (6,*) 'Step_Trunc(0)=',Step_Trunc
#endif
         Do
            Call Newq(x,nInter-nLambda,nIter,dx,W,dEdx,Err,EMx,
     &                RHS,iPvt,dg,A,nA,ed,iOptC,tBeta,
     &                nFix,ip,UpMeth,Energy,Line_Search,Step_Trunc,
     &                Restriction_Step,Thr_RS)
            If (Step_Trunc.eq.'N') Step_Trunc=' '
            If (iOpt_RS.eq.0) Exit
*
            du(1:nLambda)=dy(:)
            du(1+nLambda:nInter)=dx(:,nIter)
            Call RecPrt('du(1)',' ',du,1,nInter)
            Call Backtrans_K(T,du,dq_xy,nInter,1)
            q(:,nIter+1)=q(:,nIter)+dq_xy(:)
*
            Call Dispersion_Kriging_Layer(q(1,nIter+1),disp,nInter)
            Write (6,*) 'Disp(1)=',Disp

            disp=Restriction_Disp_Con(x(1,nIter),dx(1,nIter),
     &                                nInter-nLambda)
            Write (6,*) 'Disp(2)=',Disp
#ifdef _DEBUG_
            Write (6,*) 'disp=',disp
#endif
            fact=Half*fact
            tBeta=Half*tBeta
            If ((One-disp/Beta_Disp_.gt.1.0D-3)) Exit
            If ((fact.lt.1.0D-5) .or. (disp.lt.Beta_Disp_)) Exit
            Step_Trunc='*'
         End Do
#ifdef _DEBUG_
            Write (6,*) 'Step_Trunc(n)=',Step_Trunc
#endif
         GNrm=
     &    Sqrt(DDot_(nInter-nLambda,dEdx(1,nIter),1,dEdx(1,nIter),1))
*
      Else
*        Negative norm should serve as sign that there is no gradient
         GNrm=-One
      End If
*
#ifdef _DEBUG_
      Write (6,*) 'Step_Trunc=',Step_trunc
      Call RecPrt('Con_Opt: dx',' ',dx,nInter-nLambda,nIter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Back transform dy and dx to dq.
*
*     See Eqn. 10.
*
#ifdef _DEBUG_
*
*     dy only, constraint
*
      du(:)=Zero
      du(1:nLambda)=dy(:)
      Call DGEMM_('N','N',
     &           nInter,1,nInter,
     &            One,T,nInter,
     &            du,nInter,
     &            Zero,dq(1,nIter),nInter)
      dydy=Sqrt(DDot_(nInter,dq,1,dq,1))
      Call RecPrt('dq(dy)',' ',dq,nInter,nIter)
      Write (6,*) '<R(q_0)|dy>=',DDot_(nInter,
     &              dRdq(1,1,nIter),1,dq(1,nIter),1)
*
*     dx only, constraint minimization
*
      du(:)=Zero
      du(1+nLambda:nInter)=dx(:,nIter)
      Call DGEMM_('N','N',
     &            nInter,1,nInter,
     &            One,T,nInter,
     &            du,nInter,
     &            Zero,dq(1,nIter),nInter)
      dxdx=Sqrt(DDot_(nInter,dq,1,dq,1))
      Call RecPrt('dq(dx)',' ',dq,nInter,nIter)
      Write (6,*) '<R(q_0)|dx>=',DDot_(nInter,
     &              dRdq(1,1,nIter),1,dq(1,nIter),1)
#endif
*
*     dy+dx, full step
*
*     dq = T [dy, dx]
*
      du(:)=Zero
      du(1:nLambda)=dy(1:nLambda)
      du(nLambda+1:nInter)=dx(:,nIter)
      Call DGEMM_('N','N',nInter,1,nInter,
     &            One,T,nInter,
     &                du,nInter,
     &            Zero,dq(1,nIter),nInter)
*
*     Compute q for the next iteration
*
      q(:,nIter+1) = q(:,nIter) + dq(:,nIter)
*
#ifdef _DEBUG_
      Call RecPrt('Con_Opt: q',' ',q,nInter,nIter+1)
#endif
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
         Call mma_allocate(Tmp1,nInter-nLambda,Label='Tmp1')
         Call Char2Real(Lbl,Tmp1,(nInter-nLambda)*8)
         GrdMax=Zero
         Do i = 1, nInter-nLambda
            Write (Lbl(i),'(A,I3.3)') 'dEdx',i
         End Do
         Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,nInter-nLambda,
     &               dEdx(1,nIter),dx(1,nIter),Lbl)
         StpMax=StpMax_Save
         StpLbl=StpLbl_Save
         Call Real2Char(Tmp1,Lbl,(nInter-nLambda)*8)
         Call mma_deallocate(Tmp1)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iOpt_RS.ne.0) Then
         Call mma_deallocate(dy_)
         Call mma_deallocate(T_)
         Call mma_deallocate(q_)
      End If
      Call mma_deallocate(dq_xy)
      Call mma_deallocate(Hessian)
      Return
      End
