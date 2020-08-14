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
* Copyright (C) 1990, Roland Lindh                                     *
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Cnt1El2(Kernel,KrnlMm,Label,
     &                 iDCnt,iDCar,loper,rHrmt,DiffOp,
     &                  Lab_Dsk,iadd,isym,kcar,nordop)
************************************************************************
*                                                                      *
* Object: to compute the one-electron integrals. The method employed at*
*         this point is not necessarily the fastest. However, the total*
*         time for the computation of integrals will depend on the time*
*         spent in computing the two-electron integrals.               *
*         The memory at this point is assumed to be large enough to do *
*         the computation in core.                                     *
*         The data is structured with respect to four indices, two (my *
*         ny or i j) refer to primitives or basis functions and two (a *
*         b) refer to the components of the cartesian or spherical     *
*         harmonic gaussians.                                          *
*                                                                      *
* Called from: Drv1El                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              ICopy                                                   *
*              GetMem                                                  *
*              DCopy    (ESSL)                                         *
*              KrnlMm                                                  *
*              ZXia                                                    *
*              MemSO1                                                  *
*              DCR                                                     *
*              Inter                                                   *
*              SetUp1                                                  *
*              Kernel                                                  *
*              DGEMM_   (ESSL)                                         *
*              DGeTMO   (ESSL)                                         *
*              CarSph                                                  *
*              SymAd1                                                  *
*              DScal    (ESSL)                                         *
*              SOSctt                                                  *
*              PrMtrx                                                  *
*              XProp                                                   *
*              WrOne                                                   *
*              ErrOne                                                  *
*              Prop                                                    *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
*             Rewritten for gradients needed in hessian calculations   *
*             and general operator treatment                           *
*             May '95 By:                                              *
*             Anders Bernhardsson , Dept. of Theoretical Chemistry,    *
*             University  of Lund, SWEDEN.                             *
************************************************************************
      use Real_Spherical
      use iSD_data
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
c#include "print.fh"
#include "disp.fh"
#include "disp2.fh"
#include "nsd.fh"
#include "setup.fh"
* log trans   integer dcent
      Real*8 A(3), B(3), RB(3),CCoor(3)
      Character Label*8
      Integer nOp(2), ip(8),ipc(0:7),
     &          iDCRR(0:7), iDCRT(0:7), iStabM(0:7), iStabO(0:7),
     &          IndGrd(0:7)
      Logical AeqB,TstFnc,TF,IfGrd(3,2),EQ,DiffOP,DiffCnt,Trans(2)
      Integer iTwoj(0:7)
      Character*8 Lab_dsk
      Data iTwoj/1,2,4,8,16,32,64,128/
*
*     Statement functions
*
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                       nIrrep/nStab(mdc),iChTbl,iIrrep,iComp,
     &                       nStab(mdc))
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
*-----Compute the number of blocks from each component of the operator
*     and the irreps it will span.
*
*
      Call FZero(CCoor,3)
      Call iCopy(nIrrep,[0],0,IndGrd,1)
      loper=0
      ii=1
      Do i=0,nirrep-1
         ipC(i)=ii
         ii=ii+nBas(i)**2
      End Do
      nnIrrep=nIrrep
      If (sIrrep) nnIrrep=1
      Do iIrrep=0,nnIrrep-1
         jIrrep=nropr(ieor(ioper(iIrrep), ioper(isym)),
     &          ioper,nirrep)
         nDisp = IndDsp(iDcnt,iIrrep)
         Do iCar=1,3
            iComp = 2**(iCar-1)
            If ( TF(iDCnt,iIrrep,iComp)) Then
               ndisp=ndisp+1
               If (iDCar.eq.icar) Then
                  loper=loper+2**jIrrep
                  IndGrd(jIrrep) = nDisp
               End If
            End If
         End Do
      End Do
      nIC=0
      If (loper.eq.0) Return
      Call ICopy(nIrrep,[0],0,ip,1)
      Do iIrrep =0,nIrrep-1
         If (iAnd(2**iIrrep,loper).ne.0) Then
            LenInt=nFck(iIrrep)
            nIc=nIC+1
            Call GetMem(Label,'ALLO','REAL',ip(NIC),LenInt)
            call dcopy_(LenInt,[Zero],0,Work(ip(nIC)),1)
         End If
      End Do
      Call SOS(iStabO,nStabO,1)
*
*-----Auxiliary memory allocation.
*
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Nr_Shells(nSkal)
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Double loop over shells.
*
      Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
         A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
         Do jS = 1, iS
            jShll  = iSD( 0,jS)
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jPrim  = iSD( 5,jS)
            jAO    = iSD( 7,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
            jCnttp = iSD(13,jS)
            jCnt   = iSD(14,jS)
            B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
*-------Call kernel routine to get memory requirement. Observe, however
*       that kernels which will use the HRR will allocate that
*       memory internally.
*
        maxi=maxPrm(iAng)*maxprm(jang)
        Call GetMem('Zeta','ALLO','REAL',iZeta,maxi)
        Call GetMem('Zeta','ALLO','REAL',ipZI ,Maxi)
        Call GetMem('Kappa','ALLO','REAL',iKappa,Maxi)
        Call GetMem('PCoor','ALLO','REAL',iPCoor,Maxi*3)
        Call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
*
*       Memory requirements for contraction and Symmetry
*       adoption of derivatives.
*
        MaxP= Max(MaxPrm(iAng),MaxPrm(jAng))
        MaxZeta=MaxPrm(iAng)*MaxPrm(jAng)
        MaxB= Max(MaxBas(iAng),MaxBas(jAng))
        lFinal = MaxPrm(iAng) * MaxPrm(jAng) *
     &           nElem(iAng)*nElem(jAng)*nIrrep
*
        MemKrn=Max(MemKer*Maxi,lFinal)
        Call GetMem('Kernel','ALLO','REAL',iKern,MemKrn)
*
*            Save some memory and use Scrt area for
*            transformation
*
*       Allocate memory for the final integrals all in the
*       primitive basis.
*
        Call GetMem('Final','ALLO','REAL',ipFnl,lFinal)
*
*       Scratch area for the transformation to spherical gaussians
*
        nScr1=MaxBas(iAng)*MaxBas(jAng)*nElem(iAng)*nElem(jAng)*nIC
        Call GetMem('ScrSph','ALLO','REAL',iScrt1,nScr1)
*
*         At this point we can compute Zeta.
*         This is now computed in the ij or ji order.
*
          Call ZXia(Work(iZeta),Work(ipZI),
     &              iPrim,jPrim,Shells(iShll)%Exp,
     &                          Shells(jShll)%Exp)
*
            DiffCnt=(mdci.eq.iDCnt).or.(mdcj.eq.iDCnt)
            If ((.not.DiffCnt).and.(.not.DiffOp)) Goto 131
            AeqB = iS.eq.jS
            Call lCopy(6,[.false.],0,IfGrd,1)
            Call lCopy(2,[.false.],0,trans,1)
            If (mdci.eq.iDCnt) Then
                IfGrd(idCar,1)=.true.
            End If
            If (mdcj.eq.iDCnt) Then
                IfGrd(idCar,2)=.true.
            End If
*
            If (IfGrd(iDCar,1).and.IfGrd(iDCar,2).and.
     &          (.not.DiffOp)) Then
              IfGrd(iDCar,2)=.false.
              Trans(2)=.true.
            End If
            If (Label.eq.'CONNECTI') Trans(2)=.false.
*
*           Allocate memory for SO integrals that will be generated by
*           this batch of AO integrals.
*
            nSO=0
            Do iIrrep=0,nIrrep-1
                If (iAnd(loper,2**iIrrep).ne.0) Then
                 iSmLbl=2**iIrrep
                 nSO=nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
               End If
            End Do
c           If (iPrint.ge.29) Write (*,*) ' nSO=',nSO
            If (nSO.eq.0) Go To 131
            Call GetMem(' SO ','ALLO','REAL',ipSO,nSO*iBas*jBas)
            call dcopy_(nSO*iBas*jBas,[Zero],0,Work(ipSO),1)
*
*           Find the DCR for A and B
*
            Call DCR(LmbdR,iOper,nIrrep,
     &               jStab(0,mdci),nStab(mdci),
     &               jStab(0,mdcj),nStab(mdcj),iDCRR,nDCRR)
*
*           Find the stabilizer for A and B
*
            Call Inter(jStab(0,mdci),nStab(mdci),
     &                 jStab(0,mdcj),nStab(mdcj),
     &                 iStabM,nStabM)
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &               iDCRT,nDCRT)
*
*           Compute normalization factor
*
            iuv = nStab(mdci)*nStab(mdcj)
            Fact = DBLE(iuv*nStabO) / DBLE(nIrrep**2 * LmbdT)
            If (MolWgh.eq.1) Then
               Fact = Fact * DBLE(nIrrep)**2 / DBLE(iuv)
            Else If (MolWgh.eq.2) Then
               Fact = sqrt(DBLE(iuv))*DBLE(nStabO)/DBLE(nIrrep*LmbdT)
            End If
*
*           Loops over symmetry operations acting on the basis.
*
            nOp(1) = NrOpr(0,iOper,nIrrep)
            if(jBas.lt.-999999) write(6,*) 'gcc overoptimization',nDCRR
            Do 140 lDCRR = 0, nDCRR-1
             RB(1) = DBLE(iPhase(1,iDCRR(lDCRR)))*B(1)
             RB(2) = DBLE(iPhase(2,iDCRR(lDCRR)))*B(2)
             RB(3) = DBLE(iPhase(3,iDCRR(lDCRR)))*B(3)
             nOp(2) = NrOpr(iDCRR(lDCRR),iOper,nIrrep)
             If (Label.ne.'CONNECTI'
     &           .and.EQ(A,RB).and. (.Not.DiffOp)) Go To 140
*
*            Compute kappa and P.
*
             Call Setup1(Shells(iShll)%Exp,iPrim,
     &                   Shells(jShll)%Exp,jPrim,
     &                   A,RB,Work(iKappa),Work(iPCoor),Work(ipZI))
*
*            Compute AO integrals.
*            for easy implementation of NA integrals.
*
             Call Kernel(Shells(iShll)%Exp,iPrim,
     &                   Shells(jShll)%Exp,jPrim,
     &                   Work(iZeta),Work(ipZI),
     &                   Work(iKappa),Work(iPCoor),
     &                   Work(ipFnl),iPrim*jPrim,
     &                   iAng,jAng,A,RB,nOrder,Work(iKern),
     &                   MemKrn,Ccoor,nOrdOp,IfGrd,IndGrd,nop,
     &                   nStab(mdci),
     &                   nStab(mdcj),nic,idcar,idcnt,
     &                   iStabM,nStabM,trans,kcar,isym)
*
*
*        Transform from primitive to contracted basis functions.
*        Order of transformation is fixed. It has been shown through
*        testing that the index order ij,ab will give a performance
*        that is up to 20% faster than the ab,ij index order.
*
*
*            Transform i,jabx to jabx,I
             kk=nElem(iAng)*nElem(jAng)*nIC
             Call DGEMM_('T','N',
     &                   jPrim*kk,iBas,iPrim,
     &                   1.0d0,Work(ipFnl),iPrim,
     &                         Shells(iShll)%pCff,iPrim,
     &                   0.0d0,Work(iKern),jPrim*kk)
*
*            Transform j,abxI to abxI,J
*
             Call DGEMM_('T','N',
     &                   kk*iBas,jBas,jPrim,
     &                   1.0d0,Work(iKern),jPrim,
     &                         Shells(jShll)%pCff,jPrim,
     &                   0.0d0,Work(ipFnl),kk*iBas)
*
*            Transform to spherical gaussians if needed.
*
                 kk=nElem(iAng)*nElem(jAng)
*
                 If (Shells(iShll)%Transf.or.Shells(jShll)%Transf) Then
*
*             Result comes back as IJAB or IJAb
*
                   Call CarSph(Work(ipFnl),kk,iBas*jBas*nIC,
     &                    Work(iKern),nScr1,
     &                    RSph(ipSph(iAng)),iAng,
     &                    Shells(iShll)%Transf,Prjct(iShll),
     &                    RSph(ipSph(jAng)),
     &                    jAng,Shells(jShll)%Transf,
     &                    Prjct(jShll),Work(iScrt1),iCmp*jCmp)
*
                  Call DGeTmO(Work(iScrt1),nIC,nIC,
     &                    iBas*jBas*iCmp*jCmp,
     &                    Work(iKern),iBas*jBas*iCmp*jCmp)

*
                Else
*
*             Transpose abx,IJ back to IJ,abx
*
                    Call DGeTmO(Work(ipFnl),kk*nIC,kk*nIC,
     &                   iBas*jBas,Work(iKern),iBas*jBas)
                End If
*
*            At this point accumulate the batch of integrals onto the
*            final symmetry adapted integrals.
*
*               If (iPrint.ge.99) Then
*                 Call RecPrt (' Accumulated SO integrals, so far...',
*    &                               ' ',Work(ipSO),iBas*jBas,nSO)
*               End If
*
*------------Symmetry adapt component by component
*
             iSOBlk = ipSO
             iIC=1
             Do iIrrep = 0, nIrrep-1
                iSmLbl=iAnd(lOper,iTwoj(iIrrep))
                mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
                If (mSO.eq.0) Then
                   Do jIrrep = 0, nIrrep-1
                      If (iAnd(iSmLbl,iTwoj(jIrrep)).ne.0) iIC = iIC + 1
                   End Do
                Else
                   Call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,
     &                         iShell,jShell,iShll,jShll,
     &                         Work(iKern),iBas,jBas,nIC,iIC,
     &                         Work(iSOBlk),mSO,nOp)
                   iSOBlk = iSOBlk + mSO*iBas*jBas
                End If
             End Do
*
 140        Continue
*
*           Multiply with factors due to projection operators
*
           If (Fact.ne.One) Call DScal_(nSO*iBas*jBas,Fact,Work(ipSO),1)
*
*           Scatter the SO's on to the non-zero blocks of the
*           lower triangle.
*
             iSOBlk=ipSO
             iIC=0
             Do  iIrrep = 0, nIrrep-1
               If (iAnd(lOper,2**iIrrep).ne.0) Then
                 iSmlbl=2**iIrrep
                 iiC=iiC+1
                 mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
                 If (nfck(iirrep).ne.0.and.mSO.ne.0)
     &            Call SOSctt(Work(iSOBlk),iBas,jBas,mSO,
     &                    Work(ip(iIC)),
     &                    nFck(iIrrep),iSmLbl,
     &                    iCmp,jCmp,iShell,jShell,iAO,jAO,
     &                    nIC,Label,2**iIrrep,rHrmt)
                 iSOBlk = iSOBlk + mSO*iBas*jBas
               End If
             End Do
*
            Call GetMem('  SO ','FREE','REAL',ipSO,nSO*iBas*jBas)
 131        Continue
        Call GetMem('Kappa','FREE','REAL',iKappa,Maxi)
        Call GetMem('PCoor','FREE','REAL',iPCoor,Maxi*3)
        Call GetMem('Zeta','FREE','REAL',ipZI ,Maxi)
        Call GetMem('Zeta','FREE','REAL',iZeta,Maxi)
        Call GetMem('ScrSph','Free','REAL',iScrt1,nScr1)
        Call GetMem('Final','FREE','REAL',ipFnl,lFinal)
        Call GetMem('Kernel','FREE','REAL',iKern,MemKrn)
         End Do
      End Do
      Call Free_iSD()
*
*     Compute properties or write integrals to disc and
*     deallocate core.
*
      ipOut = 0
      mDim = 0
      nDens=0
      ipNuc = 0
      ndenssq=0
      Do iI=0,nIrrep-1
         ndenssq=ndenssq+nbas(ii)**2
         nDens=nDens+nBas(iI)*(nBas(iI)+1)/2
      End Do
      nrOp=0

      Call Getmem('Temp','ALLO','REAL',ipscr,ndenssq)
      Do 16 iIrrep = 0, nIrrep-1
         iSmLbl = 2**iIrrep
         If (iAnd(ismLbl,loper).ne.0) Then
            nrOp=nrOp+1
            jdisp=indgrd(iirrep)
            kOper=2**iIrrep
*           Write(*,*) koper,isym,jdisp,iirrep
            If (iadd.ne.0) Then
               irc=-1
               iopt=0
               iipscr=ip_of_iWork_d(work(ipscr))
               call rdmck(irc,iOpt,Lab_dsk,jdisp,iwork(iipscr),koper)
               If (irc.ne.0)
     & Call SysAbendMsg('cnt1el2','error during read in rdmck',' ')
               call daxpy_(nfck(iirrep),one,
     &                    work(ipscr),1,
     &                   work(ip(nrop)),1)
            End If
            irc=-1
            iopt=0
*           Write(*,*) Lab_dsk,jdisp,koper
            call dwrmck(irc,iOpt,Lab_dsk,jdisp,work(ip(nrop)),koper)
            If (irc.ne.0)
     & Call SysAbendMsg('cnt1el2','error during write in dwrmck',' ')
         Call GetMem(Label,'FREE','REAL',ip(nrOp),nFck(iIrrep))
         End If
 16   Continue
*

      Call Getmem('Temp','FREE','REAL',ipscr,2*ii)
*
      Return
      End
