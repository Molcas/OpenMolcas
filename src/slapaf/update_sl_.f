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
      Subroutine Update_sl_(kIter,iInt,nFix,nInter,qInt,Shift,
     &                     Grad,iOptC,Beta,Beta_Disp,Lbl,GNrm,
     &                     Energy,UpMeth,ed,Line_Search,Step_Trunc,
     &                     nLambda,iRow_c,nsAtom,AtomLbl,nSym,iOper,
     &                     mxdc,jStab,nStab,BMx,Smmtrc,nDimBC,
     &                     rLambda,Cx,GrdMax,StpMax,GrdLbl,StpLbl,
     &                     iNeg,nLbl,Labels,nLabels,FindTS,TSC,nRowH,
     &                     nWndw,Mode,MF,
     &                     iOptH,HUpMet,mIter,GNrm_Threshold,IRC,
     &                     dMass,HrmFrq_Show,CnstWght,Curvilinear,
     &                     Degen,Kriging_Hessian,qBeta,iOpt_RS,
     &                     First_MicroIteration,Iter,qBeta_Disp)
************************************************************************
*     Object: to update coordinates                                    *
*                                                                      *
*    Input:                                                            *
*      kIter          : iteration counter                              *
*      iInt           : number of internal coordinates to vary         *
*      nFix           : number of frozen internal coordinates          *
*      nInter         : total number of internal coordinates           *
*      qInt(*,kIter)  : the internal coordinates                       *
*      Shift(*,kIter) : the shift of the internal coordinates          *
*      Grad(*,kIter)  : the gradient in the internal coordinates       *
*      iOptC          : option flag for update methods                 *
*      Beta           : damping factor step length                     *
*      Beta_Disp      : damping factor variance                        *
*      Lbl            : character labels for internal coordinates      *
*      nLbl           : length of Lbl                                  *
*      GNrm           : the norm of the gradient in each iteration     *
*      Energy         : the energy of each iteration                   *
*      Line_Search    : logical flag for line search                   *
*      nLambda        : number of constraints                          *
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
*      Cx             : Cartesian coordinates                          *
*      iNeg           : Hessian index                                  *
*      Labels         : character string of primitive int. coord.      *
*      nLabels        : length of Labels                               *
*      CnstWght       : constraints weight                             *
*                                                                      *
*    OutPut:                                                           *
*      qInt(*,kIter+1): the internal coordinates to be used in the     *
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
#include "WrkSpc.fh"
#include "print.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
      Real*8 qInt(nInter,kIter+1), Shift(nInter,kIter),
     &       Grad(nInter,kIter), GNrm(kIter), Energy(kIter),
     &       dMass(nsAtom), BMx(3*nsAtom,3*nsAtom),
     &       rLambda(nLambda,kIter+1), Degen(3*nsAtom), MF(3*nsAtom),
     &       Cx(3*nsAtom,kIter+1)
      Integer iOper(0:nSym-1), jStab(0:7,nsAtom), nStab(nsAtom),
     &        iNeg(2)
*    &        iNeg(2), jNeg(2)
      Logical Line_Search, Smmtrc(3*nsAtom),FindTS, TSC, HrmFrq_Show,
     &        Found, Curvilinear, Kriging_Hessian, First_MicroIteration,
     &        Corrected
      Character Lbl(nLbl)*8, GrdLbl*8, StpLbl*8, Step_Trunc,
     &          Labels(nLabels)*8, AtomLbl(nsAtom)*(LENIN), UpMeth*6,
     &          HUpMet*6, File1*8, File2*8
      Real*8, Allocatable:: Hessian(:,:), Wess(:,:), AMat(:), dg(:),
     &                      RHS(:), ErrVec(:,:), EMtrx(:,:)
      Integer, Allocatable:: Pvt(:), Index(:)
      Real*8, Allocatable:: R(:,:), dRdq(:,:,:), CInt(:), CInt0(:)
      iRout=153
      iPrint=nPrint(iRout)
      Lu=6
*#define _DEBUG_
#ifdef _DEBUG_
      Write (Lu,*)'Update_sl_:iOpt_RS,Beta,Beta_Disp=',
     &                     iOpt_RS,Beta,Beta_Disp
      Call RecPrt('Update_sl_: qInt',' ',qInt,nInter,kIter)
      Call RecPrt('Update_sl_: Shift',' ',Shift,nInter,kIter-1)
      Call RecPrt('Update_sl_: GNrm',' ',GNrm,kIter,1)
#endif
*
      GrdMax=Zero
      StpMax=Zero
*
      mInter=nInter
      Step_Trunc='N'
      nA = (Max(mInter,kIter)+1)**2
      Call mma_Allocate(AMat,nA,Label='AMat')
      Call mma_Allocate(dg,mInter,Label='dg')
      Call mma_Allocate(Pvt,kIter+1,Label='Pvt')
      Call mma_Allocate(Index,kIter,Label='Index')
      Call mma_Allocate(ErrVec,mInter,(kIter+1),Label='ErrVec')
      Call mma_Allocate(EMtrx,kIter+1,kIter+1,Label='EMtrx')
      Call mma_Allocate(RHS,kIter+1)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Generate the Hessian in internal coordinates and then update it
*     according to some Hessian update method (BFGS, MSP, etc.)
*
*
      Call mma_Allocate(Hessian,nInter,nInter,Label='Hessian')
      If (Kriging_Hessian) Then
         iOptH_ = iOr(8,iAnd(iOptH,32))
      Else
         Call Mk_Hss_Q()
         iOptH_ = iOptH
      End If
      Call Get_dArray('Hss_Q',Hessian,nInter**2)
*
*     Perform the Hessian update, in case of GEK it essentially will
*     modify the Hessian if it is needed to guide 2nd order
*     optimization towards a minimum or a TS.
*
#ifdef _DEBUG_
      Write (Lu,*)
      Write (Lu,*)
      Write (Lu,*) ' *** Updating the molecular Hessian ***'
      Write (Lu,*)
#endif
      iRout=154
      jPrint=nPrint(iRout)
      Call Update_H(nWndw,Hessian,nInter,
     &              mIter,iOptC,Mode,MF,
     &              Shift(1,kIter-mIter+1),Grad(1,kIter-mIter+1),
     &              iNeg,iOptH_,HUpMet,nRowH,jPrint,GNrm(kIter),
     &              GNrm_Threshold,nsAtom,IRC,.True.,Corrected)
      If (Corrected) Step_Trunc='#'
*
*     Call RecPrt('Update_sl_: Hessian',' ',Hessian,nInter,nInter)
*     Write (6,*) 'After corrections'
*     Call DiagMtrx(Hessian,nInter,jNeg)
*
*     Save the number of internal coordinates on the runfile.
*
      Call Put_iScalar('No of Internal coordinates',mInter)
*
*     Save the force constant matrix on the runfile.
*
      Call Put_dArray('Hess',Hessian,mInter**2)
*
*     Optional frequency analysis
*
      iDo_dDipM=0
      If (HrmFrq_Show) Call GF_on_the_fly(iDo_DipM)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---- If this is a constrained optimization with the goal of finding a
*     transition state and we have one negative eigenvalue of the
*     Hessian then shift over to Mode Following RF to find the TS.
*
*     If TSConstraints were given, merge them with global constraints if
*     they exist, unless we are in TS regime.
*     If TS constraints were not given, remove global constraints when in
*     TS regime.
*
      Call qpg_darray('TanVec',Found,nRP)
      If (FindTS.and.First_MicroIteration) Then
         File1='UDC'
         File2='TSC'
         If (.not.TSC) File2=''
         If (iNeg(1).ge.1 ) Then
            If ((GNrm(kIter).le.GNrm_Threshold).or.Found) Then
*              Change to MFRF optimization.
               Mask=1+2+4+8+16+32+64+256+512+1024+2048+8192
               iOptC=iAnd(iOptC,Mask)
               Call Put_lScalar('TS Search',.True.)
               If (TSC) Then
                  File2=''
               Else
                  File1=''
               End If
            End If
         End If
         Call Merge_Constraints(File1,File2,'UDC',nLambda,iRow_c)
         Call Fix_UDC(iRow_c,nLambda,AtomLbl,nsAtom,nStab,.False.)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If ( nLambda.eq.0) Then
         M=3*nsAtom
         N=nInter
         NRHS=1
         Call Allocate_Work(ip_Tmp,M)
*
*        Iterate to avoid too large displacement in cartesians.
*        If the maximum displacement is more than 2*Beta, reduce the step.
*
*        Initial setup to ensure fCart=1.0 at first iteration.
         rInter=Beta
         fCart=Ten
         rCart=fCart*rInter
         nLoop=0
         Do While (rCart.ge.Two*Beta)
            nLoop=nLoop+1
            If (nLoop.gt.10) Exit
            If (rCart.gt.rInter) Then
              fCart=fCart*rInter/rCart
            Else
              fCart=fCart*0.9D0
            End If
            mInter=nInter+nLambda
            nA = (Max(mInter,kIter)+1)**2
*                                                                      *
************************************************************************
*                                                                      *
            gBeta=One
            xBeta=One
            gg_last=Beta
            dxdx_last=Beta
            Sf=Sqrt(Two)
            If (iOpt_RS.eq.0) Then
               kStart=Max(1,kIter-4)
               iEnd=kIter
            Else
               kStart=Max(1,Iter-4)
               iEnd=Iter
            End If
            Thr=1.0D-6
            Do iIter = kStart, iEnd
*
               If (iIter.ne.kIter) Then
                  dxdx=
     &             Sqrt(DDot_(mInter,Shift(1,iIter),1,Shift(1,iIter),1))
*
                  If (dxdx.lt.0.75D0*dxdx_last.and.
     &                dxdx.lt.(Beta-Thr)) Then
*                    Increase trust radius
C                    xBeta=xBeta*Sf
                     xBeta=Min(Two,xBeta*Sf)
                  Else If (dxdx.gt.1.25D0*dxdx_last.or.
     &                     dxdx.ge.(Beta+Thr)) Then
*                    Reduce trust radius
                     xBeta=Max(One/Five,xBeta/Sf)
                  End If
                  dxdx_last=dxdx
               End If
*
               gg=Sqrt(DDot_(nInter-nLambda,Grad(1,iIter),1,
     &                                     Grad(1,iIter),1))
               If (gg.lt.0.75D0*gg_last.and.gg.lt.(Beta-Thr)) Then
*                 Increase trust radius
C                 gBeta=gBeta*Sf
                  gBeta=Min(Two,gBeta*Sf)
               Else If (gg.gt.1.25D0*gg_last.or.gg.ge.(Beta+Thr)) Then
*                 Reduce trust radius
                  gBeta=Max(One/Five,gBeta/Sf)
               End If
               gg_last=gg
*
            End Do
            tBeta= Max(Beta*Min(xBeta,gBeta),Beta/Ten)
C           Write (6,*) 'tBeta=',tBeta
*                                                                      *
************************************************************************
*                                                                      *
*---------- Compute updated geometry in Internal coordinates
*
            fact=One
            qBeta=fCart*tBeta
            Thr_RS=1.0D-7
            Do
               Call Newq(qInt,mInter,kIter,Shift,Hessian,Grad,
     &                   ErrVec,EMtrx,RHS,
     &                   Pvt,dg,AMat,nA,
     &                   ed,iOptC,qBeta,nFix,Index,UpMeth,
     &                   Energy,Line_Search,Step_Trunc,Thr_RS)
               If (Step_Trunc.eq.'N') Step_Trunc=' '
               If (iOpt_RS.eq.0) Exit
*
               qInt(:,kIter+1)=qInt(:,kIter)+Shift(:,kIter)
               Call Dispersion_Kriging_Layer(qInt(1,kIter+1),Disp,
     &                                       nInter)
#ifdef _DEBUG_
               Write (6,*) 'Disp,Beta_Disp=',Disp,Beta_Disp
#endif
               fact=Half*fact
               qBeta=Half*qBeta
               If ((One-disp/Beta_Disp.gt.1.0D-3)) Exit
               If ((fact.lt.1.0D-5) .or. (disp.lt.Beta_Disp)) Exit
               Step_Trunc='*'
            End Do
#ifdef _DEBUG_
               Write (6,*) 'Step_Trunc=',Step_Trunc
#endif
*
            Call MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,mInter,
     &                  Grad(1,kIter),Shift(1,kIter),Lbl)
*
*---------- Set shift vector to zero for frozen internal coordinates.
*
            If (nFix.gt.0)
     &         Call dcopy_(nFix,[Zero],0,Shift(iInt+1,kIter),1)
*
*           Rough conversion to Cartesians
*
            Call Eq_Solver('T',M,N,NRHS,BMx,Curvilinear,Degen,
     &                     Shift(1,kIter),Work(ip_Tmp))
            rInter=Sqrt(dDot_(N,Shift(1,kIter),1,Shift(1,kIter),1))
            rCart=Zero
            Do i=ip_Tmp,ip_Tmp+nsAtom-1
               rCart=Max(rCart,
     &                   Sqrt(dDot_(3,Work(i),nsAtom,Work(i),nsAtom)))
            End Do
         End Do
         Call Free_Work(ip_Tmp)
         Call Put_dScalar('Max error',Zero)
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
*        Optimization with constraints using Lagrangian technique.
*
*        Allocate memory for the new arrays
*
         Call mma_allocate(R,nLambda,kIter,Label='R')
         Call mma_allocate(dRdq,nInter,nLambda,kIter,Label='dRdq')
         Call Getmem('ipT',   'Allo','Real',ipT,nInter*nInter)
         Call FZero(Work(ipT),nInter**2)
         Call Getmem('ipdy',  'Allo','Real',ipdy,nLambda)
         Call GetMem('dx',    'Allo','Real',ipdx,(nInter-nLambda)*kIter)
         Call FZero(Work(ipdx),(nInter-nLambda)*kIter)
         Call GetMem('dEdq_', 'Allo','Real',ipdEdq_,nInter*kIter)
         Call GetMem('du',    'Allo','Real',ipdu,nInter)
         Call GetMem('x','Allo','Real',ipx,(nInter-nLambda)*(kIter+1))
         Call FZero(Work(ipx),(nInter-nLambda)*(kIter+1))
         Call GetMem('dEdx','Allo','Real',ipdEdx,(nInter-nLambda)*kIter)
         Call mma_allocate(Wess,nInter-nLambda,nInter-nLambda,
     &                     Label='Wess')
         Wess(:,:)=0.0D0
         Call GetMem('Energy','Allo','Real',ipEnergy,kIter)
*
         call dcopy_(kIter,Energy,1,Work(ipEnergy),1)
*
         If (mInter.gt.nLbl) Then
            Call WarningMessage(2,'Update_Sl: mInter.gt.nLbl')
            Write (6,*) 'mInter=',mInter
            Write (6,*) 'nLbl=',nLbl
            Call Abend()
         End If
*
         nBVec=iRow_c-nLambda-1
         If (nBVec.gt.nLabels) Then
            Call WarningMessage(2,'Update_Sl: nBVec.gt.nLabels')
            Write (6,*) 'nBVec=',nBVec
            Write (6,*) 'nLabels=',nLabels
            Call Abend()
         End If
         Call GetMem('dBMx', 'Allo','Real',ipdBMx,nLambda*(3*nsAtom)**2)
         Call GetMem('BMx',  'Allo','Real',ipBMx,3*nsAtom*nLambda)
*
         Call GetMem('BVec', 'Allo','Real',ipBVec,3*nsAtom*nBVec)
         Call GetMem('Value','Allo','Real',ipValue,nBVec)
         Call GetMem('Value0','Allo','Real',ipValue0,nBVec)
         Call mma_allocate(CInt,nLambda,Label='CInt')
         Call mma_allocate(CInt0,nLambda,Label='CInt0')
         Call GetMem('Mult', 'Allo','Real',ipMult,nBVec**2)
         Call GetMem('iFlip','Allo','Inte',ip_iFlip,nBVec)
         Call GetMem('dBVec','Allo','Real',ipdBVec,nBVec*(3*nsAtom)**2)
*
*        Compute the constraints
*
         n1=3*nsAtom
         n2=n1**2
         Do lIter = 1, kIter
            Call DefInt2(Work(ipBVec),Work(ipdBVec),nBVec,Labels,
     &                   Work(ipBMx),nLambda,nsAtom,iRow_c,
     &                   Work(ipValue),cInt,cInt0,
     &                   Lbl(nInter+1),AtomLbl,Cx(1,lIter),
     &                   (lIter.eq.kIter).and.First_MicroIteration,
     &                   nSym,iOper,jStab,nStab,mxdc,
     &                   Work(ipMult),Smmtrc,nDimBC,Work(ipdBMx),
     &                   Work(ipValue0),lIter,iWork(ip_iFlip),dMass)
*
*           Assemble r
*
            R(:,lIter)=cInt(:)-cInt0(:)
*
*           Assemble dr/dq: Solve  B dr/dq = dr/dx
*
            dRdq(:,:,lIter)=Zero
*           Call RecPrt('BMx',' ',BMx,n1,nInter)
*           Call RecPrt('Work(ipBMx)',' ',Work(ipBMx),n1,nLambda)
*
            M=n1
            N=nInter
            NRHS=nLambda
*
*           Temporary fix of the dC/dx vector which always
*           is propted up with the full degeneracy factor.
*
            If (.NOT.Curvilinear) Then
               Do iLambda=1,nLambda
                  Do i = 1, n1
                     ij = (iLambda-1)*n1 + i -1
                     Work(ij+ipBMx)=Work(ij+ipBMx)/Degen(i)
                  End Do
               End Do
            End If
            LudRdX=30
            Call DaName(LudRdX,'dRdX')
            iAd=0
            Call iDaFile(LudRdX,1,[nLambda],1,iAd)
            Call iDaFile(LudRdX,1,[n1],1,iAd)
            Call dDaFile(LudRdX,1,Work(ipBMx),nLambda*n1,iAd)
            Call DaClos(LudRdX)
            Call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,
     &                     Work(ipBMx),dRdq(1,1,lIter))
*           Call RecPrt('drdq(1,1,lIter)',' ',
*    &                   drdq(1,1,lIter),nInter,nLambda)
*
         End Do     ! lIter
         Call GetMem('dBVec','Free','Real',ipdBVec,nBVec*(3*nsAtom)**2)
         Call GetMem('iFlip','Free','Inte',ip_iFlip,nBVec)
         Call GetMem('Mult', 'Free','Real',ipMult,nBVec**2)
         Call mma_deallocate(cInt)
         Call mma_deallocate(cInt0)
         Call GetMem('Value0','Free','Real',ipValue0,nBVec)
         Call GetMem('Value','Free','Real',ipValue,nBVec)
         Call GetMem('BVec', 'Free','Real',ipBVec,3*nsAtom*nBVec)
*
         If (iPrint.ge.99) Then
            Call RecPrt('R',' ',R,nLambda,kIter)
            Call RecPrt('drdq',' ',dRdq,nInter*nLambda,kIter)
            Call RecPrt('dC/dx',' ',Work(ipBMx),n1,nLambda)
            Call RecPrt('dQ/dx',' ',BMx,n1,nInter)
            Do i = 1, nLambda
               jpx=ipdBMx+(i-1)*n2
               Write (6,*) ' iLambda=',i
               Call RecPrt('d2C/dx2',' ',Work(jpx),n1,n1)
            End Do
            If (Curvilinear) Call dBPrint(nInter,nDimBC)
         End If
         Call GetMem('BMx',  'Free','Real',ipBMx,3*nsAtom*nLambda)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*------- Assemble d^2C/dQ^2
*
*        The expression is written as
*
*        dQ/dx * d^2C/dQ^2 * dQ/dx = d^2C/dx^2 - Sum (d^2Q/dx^2 * dC/dQ)
*
*        Solve first for y in
*                dQ/dx y = d^2C/dx^2 - Sum (d^2Q/dx^2 x dC/dQ)
*
#ifdef _NEW_CODE_
*        Compute Sum (d^2Q/dx^2 * dC/dQ)
*
         Call GetMem('d2Qdx2_dCdQ','Allo','Real',ip_QC,
     &                nDimBC**2*nLambda)
         If (Curvilinear) Call dBMult(dRdq(1,1,kIter),
     &                                Work(ip_QC),nInter,nDimBC,nLambda)
         Call RecPrt('drdq(1,1,kIter)',' ',drdq(1,1,kIter),n1,1)
         Call RecPrt('Work(ip_QC)',' ',Work(ip_QC),nDimBC**2,nLambda)
*
*        Subtract the term d^2C/dx^2
*
         Do k = 1, nLambda
            j = 0
            Do jx = 1, n1
               If (Smmtrc(jx)) Then
                  j = j + 1
*
                  i=0
                  Do ix = 1, n1
                     If (Smmtrc(ix)) Then
                        i = i + 1
                        ijk = (k-1)*n2 + (jx-1)*n1 + ix-1 + ipdBMx
                        ijk1= (k-1)*nDimBC**2 + (j-1)*nDimBC
     &                      + i-1 + ip_QC
                        Work(ijk) = Work(ijk) - Work(ijk1)
                     End If
                  End Do
               End If
            End Do
         End Do
         Call Free_Work(ip_QC)
#endif
*
         Call GetMem('Scr2','Allo','Real',ipScr2,nInter*n1*nLambda)
         Call GetMem('Scr1','Allo','Real',ipScr1,nInter*n1*nLambda)
*        Call RecPrt('d^2C/dx^2',' ',Work(ipdBMx),n2,nLambda)
*
*        Temporary fix of the d^2C/dx^2 vector which always
*        is propted up with the full degeneracy factor.
*
         If (.NOT.Curvilinear) Then
            Do k = 1, nLambda
               Do j = 1, n1
                  Do i = 1, n1
                     ijk = (k-1)*n2 + (j-1)*n1 + i-1 + ipdBMx
                     Work(ijk)=Work(ijk)/(Degen(i)*Degen(j))
                  End Do
               End Do
            End Do
         End If
*
         M=3*nsAtom
         N=nInter
         NRHS=n1*nLambda
*        Call RecPrt('dQ/dx',' ',BMx,n1,nInter)
         Call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,
     &                  Work(ipdBMx),Work(ipScr1))
         Call GetMem('dBMx', 'Free','Real',ipdBMx,nLambda*(3*nsAtom)**2)
         Call TRNSPS(nInter,n1*nLambda,Work(ipScr1),Work(ipScr2))
*        Call RecPrt('Scr1',' ',Work(ipScr2),n1*nLambda,nInter)
*
*        Followed by solve for x in B x = y^T for which x = d^2r/dq^2
*
         M=3*nsAtom
         N=nInter
         NRHS=nLambda*nInter
         Call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,
     &                  Work(ipScr2),Work(ipScr1))
         Call Free_Work(ipScr2)
         Call GetMem('d2L','Allo','Real',ipd2L,nInter*nInter*nLambda)
         Call TRNSPS(nInter*nLambda,nInter,Work(ipScr1),Work(ipd2L))
         Call Free_Work(ipScr1)
*        Call RecPrt('ipd2L',' ',Work(ipd2L),nInter**2,nLambda)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*------- Compute updated geometry in Internal coordinates
*
         Mode_=0
         M=3*nsAtom
         N=nInter
         NRHS=1
         Call Allocate_Work(ip_Tmp,M)
*
*        Iterate to avoid too large displacement in cartesians.
*        If the maximum displacement is more than 2*Beta, reduce the step.
*
*        Initial setup to ensure fCart=1.0 at first iteration.
         Thr_RS=1.0D-7
         rInter=Beta
         fCart=Ten
         rCart=fCart*rInter
         nLoop=0
         Do While (rCart.ge.Two*Beta)
            nLoop=nLoop+1
            If (nLoop.gt.100) Exit
            If (rCart.gt.rInter) Then
              fCart=fCart*rInter/rCart
            Else
              fCart=fCart*0.9D0
            End If
            qBeta=fCart*Beta
            Call Con_Opt(R,dRdq,Work(ipT),Grad,
     &                rLambda,qInt,Shift,Work(ipdy),Work(ipdx),
     &                Work(ipdEdq_),Work(ipdu),Work(ipx),Work(ipdEdx),
     &                Wess,GNrm(kIter),
     &                nWndw,Hessian,nInter,kIter,
     &                iOptC,Mode_,MF,iOptH_,HUpMet,jPrint,
     &                Work(ipEnergy),nLambda,mIter,nRowH,
     &                ErrVec,EMtrx,RHS,Pvt,
     &                dg,AMat,nA,ed,qBeta,qBeta_Disp,nFix,
     &                Index,UpMeth,Line_Search,Step_Trunc,Lbl,
     &                GrdLbl,StpLbl,GrdMax,StpMax,Work(ipd2L),nsAtom,
     &                IRC,CnstWght,iOpt_RS,Thr_RS,iter,
     &                First_Microiteration)
            If (iOpt_RS.eq.1) Exit
*
*           Rough conversion to Cartesians
*
            Call Eq_Solver('T',M,N,NRHS,BMx,Curvilinear,Degen,
     &                     Shift(1,kIter),Work(ip_Tmp))
            rInter=Sqrt(dDot_(N,Shift(1,kIter),1,Shift(1,kIter),1))
            rCart=Zero
            Do i=ip_Tmp,ip_Tmp+nsAtom-1
               rCart=Max(rCart,
     &                   Sqrt(dDot_(3,Work(i),nsAtom,Work(i),nsAtom)))
            End Do
         End Do
         Call Free_Work(ip_Tmp)
*
         If (iPrint.ge.99) Then
            Write (Lu,*)
            Write (Lu,*) '********************************************'
            Write (Lu,*) '* Lagrange multipliers for the constraints *'
            Write (Lu,*) '********************************************'
            Write (Lu,'(1X,A,2X,ES13.6)')
     &          (Lbl(nInter+iInt),-One*rLambda(iInt,mIter),
     &           iInt=1,nLambda)
            Write (Lu,*)
         End If
*
         Call Free_Work(ipd2L   )
         Call Free_Work(ipEnergy)
         Call mma_Deallocate(Wess)
         Call Free_Work(ipdEdx  )
         Call Free_Work(ipx     )
         Call Free_Work(ipdu    )
         Call Free_Work(ipdEdq_ )
         Call Free_Work(ipdx    )
         Call Free_Work(ipdy    )
         Call Free_Work(ipT     )
         Call mma_deallocate(dRdq)
         Call mma_deallocate(R)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_Deallocate(Hessian)
      Call mma_Deallocate(RHS)
      Call mma_Deallocate(EMtrx)
      Call mma_Deallocate(ErrVec)
      Call mma_Deallocate(Index)
      Call mma_Deallocate(Pvt)
      Call mma_Deallocate(dg)
      Call mma_Deallocate(AMat)
*
      Return
      End
