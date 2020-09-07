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
* Copyright (C) 1990,1991,1993,1999, Roland Lindh                      *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drv_Fck(Label,ip,lOper,nComp,CCoor,
     &                   nOrdOp,rNuc,rHrmt,iChO,
     &                   opmol,ipad,opnuc,iopadr,idirect,isyop,
     &                   PtChrg,nGrid,iAddPot)
      Use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "real.fh"
#include "warnings.fh"
      Character Label*8
      Real*8 CCoor(3,nComp), rNuc(nComp), PtChrg(nGrid)
      Integer ip(nComp), lOper(nComp), iChO(nComp), iStabO(0:7)
      Dimension opmol(*),opnuc(*),iopadr(*)
      Real*8, Dimension(:), Allocatable :: Int1El
      Integer iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 112
      iPrint = nPrint(iRout)
      Call qEnter('Drv_Fck')
      If (iPrint.ge.19) Then
         Write (6,*) ' In OneEl: Label', Label
         Write (6,*) ' In OneEl: nComp'
         Write (6,'(1X,8I5)') nComp
         Write (6,*) ' In OneEl: lOper'
         Write (6,'(1X,8I5)') lOper
         Write (6,*) ' In OneEl: n2Tri'
         Do iComp = 1, nComp
            ip(iComp) = n2Tri(lOper(iComp))
         End Do
         Write (6,'(1X,8I5)') (ip(iComp),iComp=1,nComp)
         Call RecPrt(' CCoor',' ',CCoor,3,nComp)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the number of blocks from each component of the operator
*     and the irreps it will span.
*
      nIC = 0
      llOper = 0
      Do iComp = 1, nComp
         llOper = iOr(llOper,lOper(iComp))
         Do iIrrep = 0, nIrrep-1
            If (iAnd(lOper(iComp),iTwoj(iIrrep)).ne.0) nIC = nIC + 1
         End Do
      End Do
      If (iPrint.ge.20) Write (6,*) ' nIC =',nIC
      If (nIC.eq.0) Go To 999
      Call SOS(iStabO,nStabO,llOper)
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate memory for symmetry adapted one electron integrals.
*     Will just store the unique elements, i.e. low triangular blocks
*     and lower triangular elements in the diagonal blocks.
*
      Call ICopy(nComp,[-1],0,ip,1)
      LenTot=0
      Do iComp = 1, nComp
         LenInt=n2Tri(lOper(iComp))
         LenTot=LenTot+LenInt+4
      End Do
      Call mma_allocate(Int1El,LenTot)
      ip(1)=1
      Call DCopy_(LenTot,[Zero],0,Int1El(ip(1)),1)
      iadr=ip(1)
      do iComp = 1, nComp
         LenInt=n2Tri(lOper(iComp))
         ip(icomp)=iadr
         iadr=iadr+LenInt+4
*        Copy center of operator to work area.
         Call DCopy_(3,Ccoor(1,iComp),1,Int1El(ip(iComp)+LenInt),1)
*        Copy nuclear contribution to work area.
         Int1El(ip(iComp)+LenInt+3) = rNuc(iComp)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute all SO integrals for all components of the operator.
*
      Call Drv_Fck_Internal(Label,ip,Int1El,LenTot,lOper,nComp,CCoor,
     &                      nOrdOp,rHrmt,iChO,
     &                      opmol,opnuc,ipad,iopadr,idirect,isyop,
     &                      iStabO,nStabO,nIC,
     &                      PtChrg,nGrid,iAddPot)
*                                                                      *
************************************************************************
*                                                                      *
*                    P O S T P R O C E S S I N G                       *
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.10)    Call PrMtrx(Label,lOper,nComp,ip,Int1El)
*                                                                      *
************************************************************************
*                                                                      *
*---- Write integrals to disc.
*
      mpp_state=1
      Do iComp = 1, nComp
         iSmLbl = lOper(iComp)
*                                                                      *
************************************************************************
*                                                                      *
*------- Write integrals to disc
*
         iOpt = 0
         iRC = -1
         If (Label(1:3).eq.'PAM')
     &      Write(Label,'(A5,I3.3)') 'PAM  ',iPAMcount
c        Write(6,*) ' oneel *',Label,'*'


         Call WrOne(iRC,iOpt,Label,iComp,Int1El(ip(iComp)),iSmLbl)

         If (Label(1:3).eq.'PAM')
     &      Call WrOne(iRC,iOpt,Label,1,Int1El(ip(iComp)),iSmLbl)
         iPAMcount=iPAMcount+1

         If (iRC.ne.0) then
            Call qTrace
            Write(6,*) ' *** Error in subroutine ONEEL ***'
            Write(6,*) '     Abend in subroutine WrOne'
            Call Quit(_RC_IO_ERROR_WRITE_)
         End If
      End Do  ! iComp
*                                                                      *
************************************************************************
*                                                                      *
*---- Deallocate memory for integral
*
      Call mma_deallocate(Int1El)
*                                                                      *
************************************************************************
*                                                                      *
 999  Continue
      Call qExit('Drv_Fck')
      Return
      End
      Subroutine Drv_Fck_Internal(Label,ip,Int1El,LenTot,lOper,nComp,
     &                            CCoor,nOrdOp,rHrmt,iChO,opmol,opnuc,
     &                            ipad,iopadr,idirect,isyop,iStabO,
     &                            nStabO,nIC,PtChrg,nGrid,iAddPot)
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
*              DCopy    (ESSL)                                         *
*              KrnlMm                                                  *
*              ZXia                                                    *
*              MemSO1                                                  *
*              DCR                                                     *
*              Inter                                                   *
*              SetUp1                                                  *
*              Kernel                                                  *
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
*             Modified for Hermite-Gauss quadrature November '90       *
*             Modified for Rys quadrature November '90                 *
*             Modified for multipole moments November '90              *
*                                                                      *
*             Modified for general kernel routines January  91         *
*             Modified for nonsymmetrical operators February  91       *
*             Modified for better symmetry treatement October  93      *
*             Modified loop structure April 99                         *
************************************************************************
      use Real_Spherical
      use iSD_data
      use Basis_Info
      use Phase_Info
      Implicit Real*8 (A-H,O-Z)
#include "angtp.fh"
#include "info.fh"
#include "real.fh"
#include "rmat_option.fh"
#include "stdalloc.fh"
#include "lundio.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 A(3), B(3), RB(3), CCoor(3,nComp), PtChrg(nGrid)
      Character ChOper(0:7)*3, Label*8
      Integer nOp(2), ip(nComp), lOper(nComp), iChO(nComp),
     &        iDCRR(0:7), iDCRT(0:7), iStabM(0:7), iStabO(0:7)
      Integer iTwoj(0:7)
      Dimension opmol(*),opnuc(*),iopadr(*)
      Real*8, Dimension(:), Allocatable :: Zeta, ZI, SO, Fnl
      Real*8 Int1El(LenTot)
      Data iTwoj/1,2,4,8,16,32,64,128/
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement functions
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 112
      iPrint = nPrint(iRout)
*     iPrint = 99
      Call qEnter('Drv_Fck_')
*
*-----Auxiliary memory allocation.
*
      Call mma_allocate(Zeta,m2Max)
      Call mma_allocate(ZI,m2Max)
*                                                                      *
************************************************************************
*                                                                      *
      Call Nr_Shells(nSkal)
*                                                                      *
************************************************************************
*                                                                      *
*-----Double loop over shells. These loops decide the integral type
*
      Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         IndShl = iSD( 8,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
         A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
         Do jS = iS, iS
            jShll  = iSD( 0,jS)
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jPrim  = iSD( 5,jS)
            jAO    = iSD( 7,jS)
            JndShl = iSD( 8,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
            jCnttp = iSD(13,jS)
            jCnt   = iSD(14,jS)
            B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
************************************************************************
*                                                                      *
*           Allocate memory for SO integrals that will be generated by
*           this batch of AO integrals.
*
            nSO=0
            Do iComp = 1, nComp
               iSmLbl=lOper(iComp)
               nSO=nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,
     &                        IndShl,JndShl)
            End Do
            If (iPrint.ge.29) Write (6,*) ' nSO=',nSO
            If (nSO.eq.0) Go To 131
            Call mma_allocate(SO,nSO*iBas*jBas)
            Call DCopy_(nSO*iBas*jBas,[Zero],0,SO,1)
*                                                                      *
************************************************************************
*                                                                      *
            If (iPrint.ge.19) Write (6,'(A,A,A,A,A)')
     &        ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
*                                                                      *
************************************************************************
*                                                                      *
*           Allocate memory for the final integrals all in the
*           primitive basis.
            lFinal = nIC * MaxPrm(iAng) * MaxPrm(jAng) *
     &               nElem(iAng)*nElem(jAng)
            Call mma_allocate(Fnl,lFinal)
            Call dCopy_(lFinal,[Zero],0,Fnl,1)
*                                                                      *
************************************************************************
*                                                                      *
*           At this point we can compute Zeta.
*           This is now computed in the ij or ji order.
*
            Call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,
     &                                    Shells(jShll)%Exp)
*                                                                      *
************************************************************************
*                                                                      *
*           Find the DCR for A and B
*
            Call DCR(LmbdR,iOper,nIrrep,jStab(0,mdci),
     &               nStab(mdci),jStab(0,mdcj),
     &               nStab(mdcj),iDCRR,nDCRR)
*
*           Find the stabilizer for A and B
*
            Call Inter(jStab(0,mdci),nStab(mdci),
     &                 jStab(0,mdcj),nStab(mdcj),
     &                 iStabM,nStabM)
*
            Call DCR(LambdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &               iDCRT,nDCRT)
*
            If (iPrint.ge.19) Then
               Write (6,*)
               Write (6,*) ' g      =',nIrrep
               Write (6,*) ' u      =',nStab(mdci)
               Write (6,'(9A)') '(U)=',(ChOper(jStab(ii,mdci)),
     &               ii = 0, nStab(mdci)-1)
               Write (6,*) ' v      =',nStab(mdcj)
               Write (6,'(9A)') '(V)=',(ChOper(jStab(ii,mdcj)),
     &               ii = 0, nStab(mdcj)-1)
               Write (6,*) ' LambdaR=',LmbdR
               Write (6,*) ' r      =',nDCRR
               Write (6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),
     &               ii = 0, nDCRR-1)
               Write (6,*) ' m      =',nStabM
               Write (6,'(9A)') '(M)=',(ChOper(iStabM(ii)),
     &               ii = 0, nStabM-1)
            End If
*                                                                      *
************************************************************************
*                                                                      *
*           Compute normalization factor
*
            iuv = nStab(mdci)*nStab(mdcj)
            If (MolWgh.eq.1) Then
               Fact = DBLE(nStabO) / DBLE(LambdT)
            Else If (MolWgh.eq.0) Then
               Fact = DBLE(iuv*nStabO) / DBLE(nIrrep**2 * LambdT)
            Else
               Fact = Sqrt(DBLE(iuv))*DBLE(nStabO)/
     &                DBLE(nirrep*LambdT)
            End If
            Fact = One / Fact
*                                                                      *
************************************************************************
*                                                                      *
*           Loops over symmetry operations acting on the basis.
*
            nOp(1) = NrOpr(0,iOper,nIrrep)
*           Do lDCRR = 0, nDCRR-1
            Do lDCRR = 0, 0
             RB(1) = DBLE(iPhase(1,iDCRR(lDCRR)))*B(1)
             RB(2) = DBLE(iPhase(2,iDCRR(lDCRR)))*B(2)
             RB(3) = DBLE(iPhase(3,iDCRR(lDCRR)))*B(3)
             nOp(2) = NrOpr(iDCRR(lDCRR),iOper,nIrrep)
             If (iPrint.ge.49) Write (6,'(A,3F6.2,2X,3F6.2)') '*',
     &             (A(i),i=1,3),(RB(i),i=1,3)
*                                                                      *
************************************************************************
*                                                                      *
*            Pick up epsilon from memory
*
            Call FZero(Fnl,iBas*jBas*iCmp*jCmp*nIC)
            Do iB = 1, iBas
               Do jB = 1, iBas
                  ijB=(jB-1)*iBas+iB
                  Do iC = 1, iCmp
                     ijC=(iC-1)*iCmp+iC
                     iTo= + (ijC-1)*iBas**2+ijB
#ifdef _DEBUG_
                     Write (6,*) 'ijB,ijC=',ijB,ijC
                     Write (6,*) 'Fnl(iTo),Shell(iShll)%FockOp(iB,jB)=',
     &                            Fnl(iTo),Shell(iShll)%FockOp(iB,jB)
#endif
                     Fnl(iTo)=Shells(iShll)%FockOp(iB,jB)
                  End Do
               End Do
            End Do
#ifdef _DEBUG_
            Call RecPrt('EOrb',' ',Shells(iShll)%FockOp,iBas,1)
            Call RecPrt('EOrb',' ',Shells(iShll)%FockOp,iBas,iBas)
            Call RecPrt('FckInt',' ',Fnl,iBas*jBas,iCmp*jCmp*nIC)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*            At this point accumulate the batch of integrals onto the
*            final symmetry adapted integrals.
*
             If (iPrint.ge.99) Then
                Call RecPrt (' Accumulated SO integrals, so far...',
     &                               ' ',SO,iBas*jBas,nSO)
             End If
*
*------------Symmetry adapt component by component
*
             iSOBlk = 1
             iIC = 1
             Do iComp = 1, nComp
              iSmLbl=lOper(iComp)
              mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,IndShl,JndShl)
              If (mSO.eq.0) Then
                 Do iIrrep = 0, nIrrep-1
                    If (iAnd(lOper(iComp),iTwoj(iIrrep)).ne.0)
     &                  iIC = iIC + 1
                 End Do
              Else
                 Call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,
     &                       iShell,jShell,iShll,jShll,
     &                       IndShl,JndShl,Fnl,
     &                       iBas,jBas,nIC,iIC,SO(iSOBlk),mSO,nOp)
                 iSOBlk = iSOBlk + mSO*iBas*jBas
              End If
             End Do
*
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Multiply with factors due to projection operators
*
            If (Fact.ne.One) Call DScal_(nSO*iBas*jBas,Fact,SO,1)
            If (iPrint.ge.99) Then
               Write (6,*) ' Scaling SO''s', Fact
               Call RecPrt(' Accumulated SO integrals',' ',
     &                     SO,iBas*jBas,nSO)
            End If
*                                                                      *
************************************************************************
*                                                                      *
*           Scatter the SO's on to the non-zero blocks of the
*           lower triangle.
*
            iSOBlk = 1
            Do iComp = 1, nComp
              iSmLbl=lOper(iComp)
              If (n2Tri(iSmLbl).ne.0) Then
                 mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,
     &                      IndShl,JndShl)
              Else
                 mSO=0
              End If
              If (mSO.ne.0) Then
                 Call SOSctt(SO(iSOBlk),iBas,jBas,mSO,Int1El(ip(iComp)),
     &                       n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,
     &                       jShell,IndShl,JndShl,
     &                       iAO,jAO,nComp,Label,lOper,rHrmt)
                 iSOBlk = iSOBlk + mSO*iBas*jBas
              End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
            Call mma_deallocate(Fnl)
            Call mma_deallocate(SO)
*                                                                      *
************************************************************************
*                                                                      *
 131        Continue
         End Do
      End Do
*
      Call mma_deallocate(ZI)
      Call mma_deallocate(Zeta)
*
      Call qExit('Drv_Fck_')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(CCoor)
         Call Unused_integer(nOrdOp)
         Call Unused_integer_array(iChO)
         Call Unused_real_array(opmol)
         Call Unused_real_array(opnuc)
         Call Unused_integer(ipad)
         Call Unused_integer_array(iopadr)
         Call Unused_integer(idirect)
         Call Unused_integer(isyop)
         Call Unused_real_array(PtChrg)
         Call Unused_integer(iAddPot)
      End If
      End
