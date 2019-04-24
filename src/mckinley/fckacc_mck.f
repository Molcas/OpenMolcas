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
* Copyright (C) 1993,1998, Roland Lindh                                *
************************************************************************
#define _USE_OLD_CODE_
#ifdef _USE_OLD_CODE_
      Subroutine FckAcc_mck(iAng, iCmp,jCmp,kCmp,lCmp,Shijij,
     &                  iShll, iShell, kOp, nijkl,
     &                  AOInt,TwoHam,nDens,Scrt,nScrt,
     &                  iAO,iAOst,iBas,jBas,kBas,lBas,
     &                  Dij,ij1,ij2,ij3,ij4,
     &                  Dkl,kl1,kl2,kl3,kl4,
     &                  Dik,ik1,ik2,ik3,ik4,
     &                  Dil,il1,il2,il3,il4,
     &                  Djk,jk1,jk2,jk3,jk4,
     &                  Djl,jl1,jl2,jl3,jl4,
     &                  FT,nFT,
     &                  tfact,iCar,iCent,pert,indgrd,ipdisp)
************************************************************************
*                                                                      *
*  Object: to accumulate contibutions from the AO integrals directly   *
*          to the symmatry adapted Fock matrix.                        *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*          In addition to this complication we have that the order of  *
*          indicies in the integrals are not ordered canonically but   *
*          rather in an order such that the contraction step will be   *
*          optimal. Hence, special care has to be taken when tracing   *
*          the density with the integrals so that both entities have   *
*          the same order.                                             *
*                                                                      *
*          The Fock matrix is computed in lower triangular form.       *
*                                                                      *
*          The density matrix is not folded if the shell indices and   *
*          the angular indices are identical.                          *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              DNrm2_  (ESSL)                                          *
*              DGeTMO  (ESSL)                                          *
*              DGeMV   (ESSL)                                          *
*              FckDst                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden. February '93                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "disp2.fh"
      Real*8 AOInt(nijkl,iCmp,jCmp,kCmp,lCmp), TwoHam(nDens),
     &       Scrt(nScrt), FT(nFT),
     &       Dij(ij1*ij2+1,ij3,ij4),
     &       Dkl(kl1*kl2+1,kl3,kl4),
     &       Dik(ik1*ik2+1,ik3,ik4),
     &       Dil(il1*il2+1,il3,il4),
     &       Djk(jk1*jk2+1,jk3,jk4),
     &       Djl(jl1*jl2+1,jl3,jl4)
c     Logical Qij, Qkl
      Logical Shij, Shkl, Shijij, Qijij,
     &        iShij, iShkl, iQij, iQkl,
     &        iQik, iShik, iQil, iShil, iQjk, iShjk, iQjl, iShjl,
     &        lFij, lFkl, lFik, lFjl, lFil, lFjk
      Integer iAng(4), iShell(4), iShll(4), kOp(4), kOp2(4),
     &        iAO(4), iAOst(4),
     &        iCmpa(4)
      Logical Pert(0:nIrrep-1)
      integer indgrd(3,4,0:nirrep-1),ipdisp(*)

*     Local Arrays
      Integer iSym(4,0:7), iTwoj(0:7)
      Real*8 Prmt(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement Function
*
      iOff(ixyz)  = ixyz*(ixyz+1)*(ixyz+2)/6
      xPrmt(i,j) = Prmt(iAnd(i,j))
*
      iprint=0
*
*     Write (*,*) DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1),
*    &            DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One  ,0)
*     If (DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1).gt.Zero
*    &    .or.
*    &    DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One,0).gt.Zero) Then
*     If (iPrint.ge.49) Then
*     Call RecPrt('Dik','(5G20.10)',Dik,ik1*ik2+1,ik3*ik4)
*        Write (*,'(A,2G20.10)')
*    &                  ' FckAcc:AOIn',DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
*    &               AOInt,1,AOInt,1), DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
*    &               AOInt,1,One,0)
*     End If
*     End If
*     If (iPrint.ge.99) Then
*        Call RecPrt('FckAcc:AOInt','(5G20.10)',
*    &                       AOInt,nijkl,iCmp*jCmp*kCmp*lCmp)
*        Write (*,'(A,G20.10)')'Dij=',XDot(Dij,ij1,ij2,ij3,ij4)
*        Write (*,'(A,G20.10)')'Dkl=',XDot(Dkl,kl1,kl2,kl3,kl4)
*        Write (*,'(A,G20.10)')'Dik=',XDot(Dik,ik1,ik2,ik3,ik4)
*        Write (*,'(A,G20.10)')'Dil=',XDot(Dil,il1,il2,il3,il4)
*        Write (*,'(A,G20.10)')'Djk=',XDot(Djk,jk1,jk2,jk3,jk4)
*        Write (*,'(A,G20.10)')'Djl=',XDot(Djl,jl1,jl2,jl3,jl4)
*     End If
*
*     Write (*,'(A,8L1)') 'Pert=',Pert
      If (iBas*jBas*kBas*lBas.gt.nScrt) Then
         Write (6,*) 'FckAcc_McK: iBas*jBas*kBas*lBas.gt.nScrt'
         Write (6,*) 'iBas,jBas,kBas,lBas,nScrt=',
     &         iBas,jBas,kBas,lBas,nScrt
         Call QTrace
         Call Abend()
      End If
      ii = iOff(iAng(1))
      jj = iOff(iAng(2))
      kk = iOff(iAng(3))
      ll = iOff(iAng(4))
      kOp2(1) = iOper(kOp(1))
      kOp2(2) = iOper(kOp(2))
      kOp2(3) = iOper(kOp(3))
      kOp2(4) = iOper(kOp(4))
      iCmpa(1) = iCmp
      iCmpa(2) = jCmp
      iCmpa(3) = kCmp
      iCmpa(4) = lCmp
      lFij = .False.
      lFkl = .False.
      lFik = .False.
      lFjl = .False.
      lFil = .False.
      lFjk = .False.
*
      ipFij = 1
      nFij  = iBas*jBas*iCmpa(1)*iCmpa(2)
*
      ipFkl = ipFij + nFij
      nFkl  = kBas*lBas*iCmpa(3)*iCmpa(4)
*
      ipFik = ipFkl + nFkl
      nFik  = iBas*kBas*iCmpa(1)*iCmpa(3)
*
      ipFjl = ipFik + nFik
      nFjl  = jBas*lBas*iCmpa(2)*iCmpa(4)
*
      ipFil = ipFjl + nFjl
      nFil  = iBas*lBas*iCmpa(1)*iCmpa(4)
*
      ipFjk = ipFil + nFil
      nFjk  = jBas*kBas*iCmpa(2)*iCmpa(3)
*
      call dcopy_(nFij+nFkl+nFik+nFjl+nFil+nFjk,[Zero],0,FT(ipFij),1)
*
*     Quadruple loop over elements of the basis functions angular
*     description. Loops are reduced to just produce unique SO integrals
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
      iShij = iShell(1).eq.iShell(2)
      iShkl = iShell(3).eq.iShell(4)
      iShik = iShell(1).eq.iShell(3)
      iShil = iShell(1).eq.iShell(4)
      iShjk = iShell(2).eq.iShell(3)
      iShjl = iShell(2).eq.iShell(4)
      Do 100 i1 = 1, iCmp
         Do 101 j = 0, nIrrep-1
            iSym(1,j) =
     &        iAnd(IrrCmp(IndS(iShell(1))+i1),iTwoj(j))
101      Continue
         jCmpMx = jCmp
         If (Shij) jCmpMx = i1
         iChBs = iChBas(ii+i1)
         If (Transf(iShll(1))) iChBs = iChBas(iSphCr(ii+i1))
         pEa = xPrmt(iOper(kOp(1)),iChBs)
         Do 200 i2 = 1, jCmpMx
            Do 201 j = 0, nIrrep-1
               iSym(2,j) = iAnd(IrrCmp(IndS(iShell(2))+i2),iTwoj(j))
201         Continue
            jChBs = iChBas(jj+i2)
            If (Transf(iShll(2))) jChBs = iChBas(iSphCr(jj+i2))
            pRb = xPrmt(iOper(kOp(2)),jChBs)
*           Qij = i1.eq.i2
            If (iShell(2).gt.iShell(1)) Then
               i12 = jCmp*(i1-1) + i2
            Else
               i12 = iCmp*(i2-1) + i1
            End If
            Do 300 i3 = 1, kCmp
               Do 301 j = 0, nIrrep-1
                  iSym(3,j) = iAnd(IrrCmp(IndS(iShell(3))+i3),iTwoj(j))
301            Continue
               lCmpMx = lCmp
               If (Shkl) lCmpMx = i3
               kChBs = iChBas(kk+i3)
               If (Transf(iShll(3))) kChBs = iChBas(iSphCr(kk+i3))
               pTc = xPrmt(iOper(kOp(3)),kChBs)
               Do 400 i4 = 1, lCmpMx
                  Do 401 j = 0, nIrrep-1
                     iSym(4,j) =
     &                 iAnd(IrrCmp(IndS(iShell(4))+i4),iTwoj(j))
401               Continue
*                 Qkl = i3.eq.i4
                  lChBs = iChBas(ll+i4)
                  If (Transf(iShll(4))) lChBs = iChBas(iSphCr(ll+i4))
                  pTSd= xPrmt(iOper(kOp(4)),lChBs)
                  If (iShell(4).gt.iShell(3)) Then
                     i34 = lCmp*(i3-1) + i4
                  Else
                     i34 = kCmp*(i4-1) + i3
                  End If
                  If (Shijij .and. i34.gt.i12) Go To 400
                  Qijij = Shijij .and. i12.eq.i34
                  iQij = iShij.and.i1.eq.i2
                  iQkl = iShkl.and.i3.eq.i4
                  iQik = iShik.and.i1.eq.i3
                  iQil = iShil.and.i1.eq.i4
                  iQjk = iShjk.and.i2.eq.i3
                  iQjl = iShjl.and.i2.eq.i4
                  pFctr=pEa*pRb*pTc*pTSd
*
                  mFij = 0
                  mFkl = 0
                  mFik = 0
                  mFjl = 0
                  mFil = 0
                  mFjk = 0
                  Do 1000 iIrrep = 0, nIrrep-1
                     If (iSym(1,iIrrep).ne.0 .and.
     &                   iSym(2,iIrrep).ne.0) mFkl=mFkl+1
                     If (iSym(3,iIrrep).ne.0 .and.
     &                   iSym(4,iIrrep).ne.0) mFij=mFij+1
                     If (iSym(1,iIrrep).ne.0 .and.
     &                   iSym(3,iIrrep).ne.0) mFjl=mFjl+1
                     If (iSym(2,iIrrep).ne.0 .and.
     &                   iSym(4,iIrrep).ne.0) mFik=mFik+1
                     If (iSym(1,iIrrep).ne.0 .and.
     &                   iSym(4,iIrrep).ne.0) mFjk=mFjk+1
                     If (iSym(2,iIrrep).ne.0 .and.
     &                   iSym(3,iIrrep).ne.0) mFil=mFil+1
 1000             Continue
                  If (mFij+mFkl+mFik+mFjl+mFil+mFjk.eq.0) Go To 400
*
                  vijkl = DNrm2_(iBas*jBas*kBas*lBas,
     &                          AOInt(1,i1,i2,i3,i4),1)
                  If (vijkl.lt.CutInt) Go To 400
************************************************************************
*                                                                      *
*-----------------Fij = hij + Sum(kl) Dkl {(ij|kl)-1/2(ik|jl)}         *
*                                                                      *
*-----------------or                                                   *
*                                                                      *
*-----------------Fij = hij + Sum(k=>l) Dkl {(2-d(kl)} P(ij|kl)        *
*                                                                      *
*-----------------where P(ij|kl)=(ij|kl)-1/4(ik|jl)-1/4(il|jk)         *
*                                                                      *
*-----------------or in the case of no sum restriction                 *
*                                                                      *
*-----------------P(ij|kl) = (ij|kl) - 1/2(ik|jl)                      *
*                                                                      *
*-----------------C o l o u m b  C o n t r i b u t i o n s             *
*                                                                      *
*---------------- Fij = Dkl * (ij|kl)                                  *
*                                                                      *
*-----------------Order density matrix in accordance with the integrals*
*                                                                      *
************************************************************************
                  If (mFij.eq.0) Go To 3203
                  If (iShell(3).lt.iShell(4)) Then
                     vkl = Dkl(kBas*lBas+1,i4,i3)
                  Else
                     vkl = Dkl(kBas*lBas+1,i3,i4)
                  End If
*
*-----------------Pickup the right column of the density matrix and
*                 change order if not canonical.
*
                  qFctr = One
                  ipFij1 = ((i2-1)*iCmpa(1)+i1-1)*iBas*jBas+ipFij
                  Fac = One
                  If (iQij) Fac = Half
                  D_kl=Two
                  If (iQkl) D_kl=One
                  Fac=Fac*D_kl
                  If (vkl*vijkl*Abs(Fac*qFctr*pFctr).lt.CutInt)
     &               Go To 3203
                  If (iShell(3).lt.iShell(4)) Then
                     Call DGetMO(Dkl(1,i4,i3),kl1,kl1,
     &                           kl2,Scrt,kl2)
                     Call dGeMV_('N',iBas*jBas,kBas*lBas,
     &                 Fac * qFctr*pFctr,AOInt(1,i1,i2,i3,i4),iBas*jBas,
     &                          Scrt,1,
     &                          One,FT(ipFij1),1)
                  Else
                     Call dGeMV_('N',iBas*jBas,kBas*lBas,
     &                 Fac * qFctr*pFctr,AOInt(1,i1,i2,i3,i4),iBas*jBas,
     &                          Dkl(1,i3,i4),1,
     &                          One,FT(ipFij1),1)
                  End If
C                 Call RecPrt('Fij',' ',FT(ipFij1),iBas,jBas)
                  lFij = .True.
 3203             Continue
                  If (Qijij) Go To 1200
*---------------- Fkl = Dij * (ij|kl)
                  If (mFkl.eq.0) Go To 1200
                  If (iShell(1).lt.iShell(2)) Then
                     vij = Dij(iBas*jBas+1,i2,i1)
                  Else
                     vij = Dij(iBas*jBas+1,i1,i2)
                  End  If
                  qFctr = One
                  ipFkl1 = ((i4-1)*iCmpa(3)+i3-1)*kBas*lBas+ipFkl
                  Fac = One
                  If (iQkl) Fac = Half
                  D_ij=Two
                  If (iQij) D_ij=One
                  Fac=Fac*D_ij
                  If (vij*vijkl*Abs(Fac*qFctr*pFctr).lt.CutInt)
     &               Go To 1200
                  If (iShell(1).lt.iShell(2)) Then
                     Call DGeTMO(Dij(1,i2,i1),ij1,ij1,
     &                           ij2,Scrt,ij2)
                     Call dGeMV_('T',iBas*jBas,kBas*lBas,
     &                  Fac* qFctr*pFctr,AOInt(1,i1,i2,i3,i4),iBas*jBas,
     &                          Scrt,1,
     &                          One,FT(ipFkl1),1)
                  Else
                     Call dGeMV_('T',iBas*jBas,kBas*lBas,
     &                  Fac* qFctr*pFctr,AOInt(1,i1,i2,i3,i4),iBas*jBas,
     &                          Dij(1,i1,i2),1,
     &                          One,FT(ipFkl1),1)
                  End If
C                 Call RecPrt('Fkl',' ',FT(ipFkl1),kBas,lBas)
                  lFkl = .True.
 1200             Continue
*
*-----------------E x c h a n g e   c o n t r i b u t i o n s
*
*-----------------Change the order ijkl to ikjl. Make sure also that
*                 the index pairs are canonical.
*
                  ipD = 1 + iBas*jBas*kBas*lBas
                  np = ipD-1 + Max(iBas*lBas,jBas*lBas,
     &                             iBas*kBas,jBas*kBas)
                  If (np.gt.nScrt) Then
                     Write (6,*) 'FckAcc_McK: np.gt.nScrt'
                     Write (6,*) 'np,nScrt=',np,nScrt
                     Call QTrace
                     Call Abend()
                  End If
                  If (mFik+mFjl.eq.0) Go To 1210
                  If (mFik.ne.0.and.iShell(2).lt.iShell(4)) Then
                     vjl=Djl(jBas*lBas+1,i4,i2)
                  Else If (mFik.ne.0) Then
                     vjl=Djl(jBas*lBas+1,i2,i4)
                  Else
                     vjl=Zero
                  End If
                  If (mFjl.ne.0.and.iShell(1).lt.iShell(3)) Then
                     vik=Dik(iBas*kBas+1,i3,i1)
                  Else If (mFjl.ne.0) Then
                     vik=Dik(iBas*kBas+1,i1,i3)
                  Else
                    vik=Zero
                  End If
                  If (vik*vijkl/Four .lt. CutInt .and.
     &                vjl*vijkl/Four .lt. CutInt) Go To 1210
                  Do 600 j = 1, jBas
                     nij = (j-1)*iBas + 1
                     Do 601 k = 1, kBas
                        nik = (k-1)*iBas + 1
                        Do 602 l = 1, lBas
                           nkl = (l-1)*kBas + k
                           njl = (l-1)*jBas + j
                           iOut = (nkl-1)*iBas*jBas + nij
                           iIn  = (njl-1)*iBas*kBas + nik
                           call dcopy_(iBas,AOInt(iOut,i1,i2,i3,i4),1,
     &                                Scrt(iIn),1)
 602                    Continue
 601                 Continue
 600              Continue
************************************************************************
*                                                                      *
*-----------------Fik = - 1/4 * Djl * P(jl|ik)                         *
*                                                                      *
*                 P(jl|ik) = (jl|ik) - 1/4(ji|lk) - 1/4(jk|il)         *
*                                                                      *
*                 P(jl|ik) = (jl|ik) - 1/2(ji|lk)                      *
*                                                                      *
*-----------------Change factor if                                     *
*                 a) asymmetrical P matrix is implied                  *
*                 b) if the two exchange integrals in the symmetrical  *
*                    P matrix are identical.                           *
*                                                                      *
************************************************************************
                  If (mFik.eq.0) Go To 2213
                  qFctr = One
                  ipFik1 = ((i3-1)*iCmpa(1)+i1-1)*iBas*kBas+ipFik
                  Fac = -Quart
                  If (iQjl.and. .Not.iQik) Fac = -Half
                  D_jl=Two
                  If (iQjl) D_jl=One
                  Fac=Fac*D_jl
                  If (vjl*vijkl*Abs(Fac*qFctr*pFctr).lt.CutInt)
     &               Go To 2213
                  If (iShell(2).lt.iShell(4)) Then
                     Call DGeTMO(Djl(1,i4,i2),jl1,jl1,
     &                           jl2,Scrt(ipD),jl2)
                     Call dGeMV_('N',iBas*kBas,jBas*lBas,
     &                          Fac*qFctr*pFctr,Scrt,iBas*kBas,
     &                          Scrt(ipD),1,
     &                          One,FT(ipFik1),1)
                  Else
                     Call dGeMV_('N',iBas*kBas,jBas*lBas,
     &                          Fac*qFctr*pFctr,Scrt,iBas*kBas,
     &                          Djl(1,i2,i4),1,
     &                          One,FT(ipFik1),1)
                  End If
C                 Call RecPrt('Fik',' ',FT(ipFik1),iBas,kBas)
                  lFik = .True.
 2213             Continue
                  If (iQij.and.iQkl) Go To 2220
************************************************************************
*                                                                      *
*-----------------Fjl = - 1/4 * Dik * P(jl|ik)                         *
*                                                                      *
*                 P(jl|ik) = (jl|ik) - 1/4(ji|lk) - 1/4(jk|il)         *
*                                                                      *
*                 P(jl|ik) = (jl|ik) - 1/2(ji|lk)                      *
*                                                                      *
************************************************************************
                  If (mFjl.eq.0) Go To 1210
                  qFctr = One
                  ipFjl1 = ((i4-1)*iCmpa(2)+i2-1)*jBas*lBas+ipFjl
                  Fac = -Quart
                  If (iQik.and. .Not.iQjl) Fac = -Half
                  D_ik=Two
                  If (iQik) D_ik=One
                  Fac=Fac*D_ik
                  If (vik*vijkl*Abs(Fac*qFctr*pFctr).lt.CutInt)
     &               Go To 1210
                  If (iShell(1).lt.iShell(3)) Then
                     Call DGeTMO(Dik(1,i3,i1),ik1,ik1,
     &                           ik2,Scrt(ipD),ik2)
                     Call dGeMV_('T',iBas*kBas,jBas*lBas,
     &                          Fac*qFctr*pFctr,Scrt,iBas*kBas,
     &                          Scrt(ipD),1,
     &                          One,FT(ipFjl1),1)
                  Else
                     Call dGeMV_('T',iBas*kBas,jBas*lBas,
     &                          Fac*qFctr*pFctr,Scrt,iBas*kBas,
     &                          Dik(1,i1,i3),1,
     &                          One,FT(ipFjl1),1)
                  End If
C                 Call RecPrt('Fjl',' ',FT(ipFjl1),jBas,lBas)
                  lFjl = .True.
 1210             Continue
*
*-----------------Change order ijkl to iljk
*
                  If (iQij.or.iQkl) Go To 2220
                  If (mFil+mFjk.eq.0) Go To 2220
                  If (mFil.ne.0.and.iShell(2).lt.iShell(3)) Then
                     vjk = Djk(jBas*kBas+1,i3,i2)
                  Else If (mFil.ne.0) Then
                     vjk = Djk(jBas*kBas+1,i2,i3)
                  Else
                     vjk=Zero
                  End If
                  If (mFjk.ne.0.and.iShell(1).lt.iShell(4)) Then
                     vil = Dil(iBas*lBas+1,i4,i1)
                  Else If (mFjk.ne.0) Then
                     vil = Dil(iBas*lBas+1,i1,i4)
                  Else
                     vil=Zero
                  End If
                  If (vil*vijkl/Four .lt. CutInt .and.
     &                vjk*vijkl/Four .lt. CutInt) Go To 2220
                  i = 1
                  Do 610 j = 1, jBas
                     nij = (j-1)*iBas + i
                     Do 611 k = 1, kBas
                        njk = (k-1)*jBas + j
                         Do 612 l = 1, lBas
                            nkl = (l-1)*kBas + k
                            nil = (l-1)*iBas + i
                            ijkl = (nkl-1)*iBas*jBas + nij
                            iljk = (njk-1)*iBas*lBas + nil
                            call dcopy_(iBas,AOInt(ijkl,i1,i2,i3,i4),1,
     &                                 Scrt(iljk),1)
 612                     Continue
 611                  Continue
 610               Continue
*-----------------Fil = - 1/4 * Djk * (ij|kl)
                  If (mFil.eq.0) Go To 1220
                  qFctr = One
                  ipFil1 = ((i4-1)*iCmpa(1)+i1-1)*iBas*lBas+ipFil
                  Fac = -Quart
                  If (iQjk.and. .Not.iQil) Fac = -Half
                  D_jk=Two
                  If (iQjk) D_jk=One
                  Fac=Fac*D_jk
                  If (vjk*vijkl*Abs(Fac*qFctr*pFctr).lt.CutInt)
     &               Go To 1220
                  If (iShell(2).lt.iShell(3)) Then
                     Call DGeTMO(Djk(1,i3,i2),jk1,jk1,
     &                           jk2,Scrt(ipD),jk2)
                     Call dGeMV_('N',iBas*lBas,jBas*kBas,
     &                          Fac*qFctr*pFctr,Scrt,iBas*lBas,
     &                          Scrt(ipD),1,
     &                          One,FT(ipFil1),1)
                  Else
                     Call dGeMV_('N',iBas*lBas,jBas*kBas,
     &                          Fac*qFctr*pFctr,Scrt,iBas*lBas,
     &                          Djk(1,i2,i3),1,
     &                          One,FT(ipFil1),1)
                  End If
C                 Call RecPrt('Fil',' ',FT(ipFil1),iBas,lBas)
                  lFil = .True.
 1220             Continue
                  If (Qijij) Go To 2220
*-----------------Fjk = - 1/4 * Dil * (ij|kl)
                  If (mFjk.eq.0) Go To 2220
                  qFctr = One
                  ipFjk1 = ((i3-1)*iCmpa(2)+i2-1)*jBas*kBas+ipFjk
                  Fac = -Quart
                  If (iQil.and. .Not.iQjk) Fac = -Half
                  D_il=Two
                  If (iQil) D_il=One
                  Fac=Fac*D_il
                  If (vil*vijkl*Abs(Fac*qFctr*pFctr).lt.CutInt)
     &               Go To 2220
                  If (iShell(1).lt.iShell(4)) Then
                     Call DGeTMO(Dil(1,i4,i1),il1,il1,
     &                           il2,Scrt(ipD),il2)
                     Call dGeMV_('T',iBas*lBas,jBas*kBas,
     &                          Fac*qFctr*pFctr,Scrt,iBas*lBas,
     &                          Scrt(ipD),1,
     &                          One,FT(ipFjk1),1)
                  Else
                     Call dGeMV_('T',iBas*lBas,jBas*kBas,
     &                          Fac*qFctr*pFctr,Scrt,iBas*lBas,
     &                          Dil(1,i1,i4),1,
     &                          One,FT(ipFjk1),1)
                  End If
C                 Call RecPrt('Fjk',' ',FT(ipFjk1),jBas,kBas)
                  lFjk = .True.
 2220             Continue
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
*
*     Write (*,*) ' Fij'
      nnIrrep=nIrrep
      If (sIrrep) nnIrrep=1
      Do iIrrep=0,nnIrrep-1
*        Write (*,'(I2,L1)') iIrrep, pert(iIrrep)
         If (pert(iIrrep)) Then
            ip=ipDisp(abs(indgrd(iCar,iCent,iIrrep)))
            rCh=xPrmt(iOper(kOp(iCent)),iChBas(1+iCar))*
     &                rChTbl(iIrrep,kOp(iCent))
            Fact=tfact*rCh
*           Write (*,*) 'Level ij'
            If (lFij) Call FckDst(TwoHam(ip),ndens,FT(ipFij),
     &                             iBas,jBas,iCmpa(1),iCmpa(2),
     &                             kOp2(1),kOp2(2),
     &                             iIrrep,iShell(1),iShell(2),iShij,
     &                             iAO(1),iAO(2),iAOst(1),iAOst(2),
     &                             fact)
*           Write (*,*) 'Level kl'
            If (lFkl) Call FckDst(TwoHam(ip),ndens,FT(ipFkl),
     &                             kBas,lBas,iCmpa(3),iCmpa(4),
     &                             kOp2(3),kOp2(4),
     &                             iIrrep,iShell(3),iShell(4),iShkl,
     &                             iAO(3),iAO(4),iAOst(3),iAOst(4),
     &                             fact)
*           Write (*,*) 'Level ik'
            If (lFik) Call FckDst(TwoHam(ip),ndens,FT(ipFik),
     &                             iBas,kBas,iCmpa(1),iCmpa(3),
     &                             kOp2(1),kOp2(3),
     &                             iIrrep,iShell(1),iShell(3),iShik,
     &                             iAO(1),iAO(3),iAOst(1),iAOst(3),
     &                             fact)
*           Write (*,*) 'Level jl'
            If (lFjl) Call FckDst(TwoHam(ip),ndens,FT(ipFjl),
     &                             jBas,lBas,iCmpa(2),iCmpa(4),
     &                             kOp2(2),kOp2(4),
     &                             iIrrep,iShell(2),iShell(4),iShjl,
     &                             iAO(2),iAO(4),iAOst(2),iAOst(4),
     &                             fact)
*           Write (*,*) 'Level il'
            If (lFil) Call FckDst(TwoHam(ip),ndens,FT(ipFil),
     &                             iBas,lBas,iCmpa(1),iCmpa(4),
     &                             kOp2(1),kOp2(4),
     &                             iIrrep,iShell(1),iShell(4),iShil,
     &                             iAO(1),iAO(4),iAOst(1),iAOst(4),
     &                             fact)
*           Write (*,*) 'Level jk'
            If (lFjk) Call FckDst(TwoHam(ip),ndens,FT(ipFjk),
     &                             jBas,kBas,iCmpa(2),iCmpa(3),
     &                             kOp2(2),kOp2(3),
     &                             iIrrep,iShell(2),iShell(3),iShjk,
     &                             iAO(2),iAO(3),iAOst(2),iAOst(3),
     &                             fact)
*           If (DDot_(3468,TwoHam,1,TwoHam,1).gt.Zero) Then
*           If (Abs(DDot_(3468,TwoHam,1,One,0)).gt.1.0D-16) Then
*           Write (*,'(A,G20.6,G20.7)') 'TwoHam=',
*    &             DDot_(3468,TwoHam,1,TwoHam,1),
*    &             DDot_(3468,TwoHam,1,One,0)
*           Else
*           Write (*,'(A,G20.6)') 'TwoHam=',
*    &             DDot_(3468,TwoHam,1,TwoHam,1)
*           End If
*           Call FZero(TwoHam,3468)
*           End If
         End If
      End Do
c     Call GetMem(' Exit FckAcc','CHECK','REAL',iDum,iDum)
c     Call QExit('FckAcc')
      Return
      End
#else
      Subroutine FckAcc_Mck(iAng, iCmp, jCmp, kCmp, lCmp, Shijij,
     &                  iShll, iShell, kOp, nijkl,
     &                  AOInt,TwoHam,nDens,Scrt,nScrt,
     &                  iAO,iAOst,iBas,jBas,kBas,lBas,
     &                  Dij,ij1,ij2,ij3,ij4,
     &                  Dkl,kl1,kl2,kl3,kl4,
     &                  Dik,ik1,ik2,ik3,ik4,
     &                  Dil,il1,il2,il3,il4,
     &                  Djk,jk1,jk2,jk3,jk4,
     &                  Djl,jl1,jl2,jl3,jl4,
     &                  FT,nFT,
     &                  tfact,iCar,iCent,pert,indgrd,ipdisp)
************************************************************************
*                                                                      *
*  Object: to accumulate contibutions from the AO integrals directly   *
*          to the symmatry adapted Fock matrix.                        *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*          In addition to this complication we have that the order of  *
*          indicies in the integrals are not ordered canonically but   *
*          rather in an order such that the contraction step will be   *
*          optimal. Hence, special care has to be taken when tracing   *
*          the density with the integrals so that both entities have   *
*          the same order.                                             *
*                                                                      *
*          The Fock matrix is computed in lower triangular form.       *
*                                                                      *
*          The density matrix is not folded if the shell indices and   *
*          the angular indices are identical.                          *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              DNrm2_  (ESSL)                                          *
*              DGeTMO  (ESSL)                                          *
*              DGeMV   (ESSL)                                          *
*              FckDst                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden. February '93                            *
*                                                                      *
*     Modified July '98 in Tokyo by R. Lindh                           *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "disp2.fh"
#include "print.fh"
      Real*8 AOInt(nijkl,iCmp,jCmp,kCmp,lCmp), TwoHam(nDens),
     &       Scrt(nScrt), FT(nFT),
     &       Dij(ij1*ij2+1,ij3,ij4),
     &       Dkl(kl1*kl2+1,kl3,kl4),
     &       Dik(ik1*ik2+1,ik3,ik4),
     &       Dil(il1*il2+1,il3,il4),
     &       Djk(jk1*jk2+1,jk3,jk4),
     &       Djl(jl1*jl2+1,jl3,jl4)
      Logical Shijij, Qijij,
     &        iShij, iShkl, iQij, iQkl,
     &        iQik, iShik, iQil, iShil, iQjk, iShjk, iQjl, iShjl,
     &        lFij, lFkl, lFik, lFjl, lFil, lFjk
      Integer iAng(4), iShell(4), iShll(4), kOp(4), kOp2(4),
     &        iAO(4), iAOst(4),
     &        iCmpa(4)
      Logical Pert(0:nIrrep-1)
      integer indgrd(3,4,0:nirrep-1),ipdisp(*)
*     Local Arrays
      Integer iSym(4)
      Real*8 Prmt(0:7)
c     Character*72 Label
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement Function
*
      iOff(ixyz)  = ixyz*(ixyz+1)*(ixyz+2)/6
      xPrmt(i,j) = Prmt(iAnd(i,j))
c     iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*
      iRout = 38
      iPrint = nPrint(iRout)
*
*     If (iPrint.ge.49) Then
*        Write (*,*) ' FckAcc:AOIn',DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
*    &               AOInt,1,AOInt,1), DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
*    &               AOInt,1,One,0)
*     End If
*     If (iPrint.ge.99) Then
*        Call RecPrt('FckAcc:AOInt',' ',AOInt,nijkl,iCmp*jCmp*kCmp*lCmp)
*        Write (*,*)'Dij=',XDot(Dij,ij1,ij2,ij3,ij4)
*        Write (*,*)'Dkl=',XDot(Dkl,kl1,kl2,kl3,kl4)
*        Write (*,*)'Dik=',XDot(Dik,ik1,ik2,ik3,ik4)
*        Write (*,*)'Dil=',XDot(Dil,il1,il2,il3,il4)
*        Write (*,*)'Djk=',XDot(Djk,jk1,jk2,jk3,jk4)
*        Write (*,*)'Djl=',XDot(Djl,jl1,jl2,jl3,jl4)
*     End If
*
      ExFac=One
      ThrInt=0.0D0
      If (iBas*jBas*kBas*lBas.gt.nScrt) Then
         Write (6,*) 'FckAcc: nScrt too small!'
         Call Abend
      End If
      ii = iOff(iAng(1))
      jj = iOff(iAng(2))
      kk = iOff(iAng(3))
      ll = iOff(iAng(4))

      kOp2(1) = iOper(kOp(1))
      kOp2(2) = iOper(kOp(2))
      kOp2(3) = iOper(kOp(3))
      kOp2(4) = iOper(kOp(4))
      iCmpa(1) = iCmp
      iCmpa(2) = jCmp
      iCmpa(3) = kCmp
      iCmpa(4) = lCmp
      lFij = .False.
      lFkl = .False.
      lFik = .False.
      lFjl = .False.
      lFil = .False.
      lFjk = .False.
*
      ipFij = 1
      nFij  = iBas*jBas*iCmpa(1)*iCmpa(2)
*
      ipFkl = ipFij + nFij
      nFkl  = kBas*lBas*iCmpa(3)*iCmpa(4)
*
      ipFik = ipFkl + nFkl
      nFik  = iBas*kBas*iCmpa(1)*iCmpa(3)
*
      ipFjl = ipFik + nFik
      nFjl  = jBas*lBas*iCmpa(2)*iCmpa(4)
*
      ipFil = ipFjl + nFjl
      nFil  = iBas*lBas*iCmpa(1)*iCmpa(4)
*
      ipFjk = ipFil + nFil
      nFjk  = jBas*kBas*iCmpa(2)*iCmpa(3)
*
      call dcopy_(nFij+nFkl+nFik+nFjl+nFil+nFjk,Zero,0,FT(ipFij),1)
*
      ipDij = 1
      ipDkl = 1
      ipDik = 1
      ipDil = 1
      ipDjk = 1
      ipDjl = 1
      ipFij1= 1
      ipFkl1= 1
      ipFik1= 1
      ipFil1= 1
      ipFjk1= 1
      ipFjl1= 1
*
*     Quadruple loop over elements of the basis functions angular
*     description. Loops are reduced to just produce unique SO integrals
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      iShij = iShell(1).eq.iShell(2)
      iShkl = iShell(3).eq.iShell(4)
      iShik = iShell(1).eq.iShell(3)
      iShil = iShell(1).eq.iShell(4)
      iShjk = iShell(2).eq.iShell(3)
      iShjl = iShell(2).eq.iShell(4)
      mijkl=iBas*jBas*kBas*lBas
      Do 100 i1 = 1, iCmp
         iSym(1)=IrrCmp(IndS(iShell(1))+i1)
         jCmpMx = jCmp
         If (iShij) jCmpMx = i1
         iChBs = iChBas(ii+i1)
         If (Transf(iShll(1))) iChBs = iChBas(iSphCr(ii+i1))
         pEa = xPrmt(iOper(kOp(1)),iChBs)
         Do 200 i2 = 1, jCmpMx
            iSym(2) =IrrCmp(IndS(iShell(2))+i2)
            jChBs = iChBas(jj+i2)
            If (Transf(iShll(2))) jChBs = iChBas(iSphCr(jj+i2))
            pRb = xPrmt(iOper(kOp(2)),jChBs)
            If (iShell(2).gt.iShell(1)) Then
               i12 = jCmp*(i1-1) + i2
            Else
               i12 = iCmp*(i2-1) + i1
            End If
            Do 300 i3 = 1, kCmp
               iSym(3) =IrrCmp(IndS(iShell(3))+i3)
               lCmpMx = lCmp
               If (iShkl) lCmpMx = i3
               kChBs = iChBas(kk+i3)
               If (Transf(iShll(3))) kChBs = iChBas(iSphCr(kk+i3))
               pTc = xPrmt(iOper(kOp(3)),kChBs)
               Do 400 i4 = 1, lCmpMx
                  iSym(4) =IrrCmp(IndS(iShell(4))+i4)
                  lChBs = iChBas(ll+i4)
                  If (Transf(iShll(4))) lChBs = iChBas(iSphCr(ll+i4))
                  pTSd= xPrmt(iOper(kOp(4)),lChBs)
                  If (iShell(4).gt.iShell(3)) Then
                     i34 = lCmp*(i3-1) + i4
                  Else
                     i34 = kCmp*(i4-1) + i3
                  End If
                  If (Shijij .and. i34.gt.i12) Go To 400
                  vijkl=0
                  do ijkl=1,mijkl
                    vijkl=max(vijkl,abs(AOInt(ijkl,i1,i2,i3,i4)))
                  end do
c                 vijkl = DNrm2_(iBas*jBas*kBas*lBas,
c    &                          AOInt(1,i1,i2,i3,i4),1)
C                 If (vijkl.lt.CutInt) Go To 400
*
                  Qijij = Shijij .and. i12.eq.i34
                  iQij = iShij.and.i1.eq.i2
                  iQkl = iShkl.and.i3.eq.i4
                  iQik = iShik.and.i1.eq.i3
                  iQil = iShil.and.i1.eq.i4
                  iQjk = iShjk.and.i2.eq.i3
                  iQjl = iShjl.and.i2.eq.i4
                  pFctr=pEa*pRb*pTc*pTSd
************************************************************************
*
                  Fac_ij = pFctr
                  Fac_kl = pFctr
                  Fac_ik =-Quart * pFctr
                  Fac_jl =-Quart * pFctr
                  Fac_il =-Quart * pFctr
                  Fac_jk =-Quart * pFctr
                  If (iQij) Fac_ij = Fac_ij * Half
                  If (iQkl) Fac_kl = Fac_kl * Half
                  If (iQjl.and. .Not.iQik) Fac_ik = Fac_ik * Two
                  If (iQik.and. .Not.iQjl) Fac_jl = Fac_jl * Two
                  If (iQjk.and. .Not.iQil) Fac_il = Fac_il * Two
                  If (iQil.and. .Not.iQjk) Fac_jk = Fac_jk * Two
*
                  D_ij=Two
                  If (iQij) D_ij=One
                  D_kl=Two
                  If (iQkl) D_kl=One
                  D_ik=Two
                  If (iQik) D_ik=One
                  D_jl=Two
                  If (iQjl) D_jl=One
                  D_il=Two
                  If (iQil) D_il=One
                  D_jk=Two
                  If (iQjk) D_jk=One
*
                  Fac_ij=Fac_ij*D_kl
                  Fac_kl=Fac_kl*D_ij
                  Fac_jl=Fac_jl*D_ik
                  Fac_ik=Fac_ik*D_jl
                  Fac_jk=Fac_jk*D_il
                  Fac_il=Fac_il*D_jk
*
C                 Write (*,*)
C                 Write (*,*) 'iShell(1),iShell(2),i1,i2=',
C    &                         iShell(1),iShell(2),i1,i2
C                 Write (*,*) 'Dij=',Dij(1,i1,i2)
C                 Write (*,*) 'Fac_ij,iQij=',Fac_ij,iQij
C                 Write (*,*)
C                 Write (*,*) 'iShell(3),iShell(4),i3,i4=',
C    &                         iShell(3),iShell(4),i3,i4
C                 Write (*,*) 'Dkl=',Dkl(1,i3,i4)
C                 Write (*,*) 'Fac_kl,iQkl=',Fac_kl,iQkl
C                 Write (*,*)
                  If (Qijij) Then
                     Fac_kl = Zero
                     Fac_jk = Zero
                  End If
                  If (iQij.and.iQkl) Then
                     Fac_jl = Zero
                     Fac_il = Zero
                     Fac_jk = Zero
                  End If
                  If (iQij.or.iQkl) Then
                     Fac_il = Zero
                     Fac_jk = Zero
                  End If
************************************************************************
*
                  iOpt=0
                  ip = 1
                  If (iAnd(iSym(1),iSym(2)).ne.0 .or.
     &                iAnd(iSym(3),iSym(4)).ne.0) Then
                     iOpt = iOpt + 1
*
                     If (iShell(3).lt.iShell(4)) Then
                        vkl = Dkl(kBas*lBas+1,i4,i3)
                        ipDkl=ip
                        ip = ip + kBas*lBas
                        Call DGetMO(Dkl(1,i4,i3),kl1,kl1,
     &                              kl2,Scrt(ipDkl),kl2)
                     Else
                        vkl = Dkl(kBas*lBas+1,i3,i4)
                        loc1=(idLoc(Dkl(1,i3,i4))-idLoc(Scrt))
                        loc2=(idLoc(Scrt(2))-idLoc(Scrt(1)))
                        ipDkl = 1 + loc1/loc2
                     End If
                     ipFij1 = ((i2-1)*iCmpa(1)+i1-1)*iBas*jBas
     &                      + ipFij
                     If (iShell(1).lt.iShell(2)) Then
                        vij = Dij(iBas*jBas+1,i2,i1)
                        ipDij=ip
                        ip = ip + iBas*jBas
                        Call DGeTMO(Dij(1,i2,i1),ij1,ij1,
     &                              ij2,Scrt(ipDij),ij2)
                     Else
                        vij = Dij(iBas*jBas+1,i1,i2)
                        loc1=(idLoc(Dij(1,i1,i2))-idLoc(Scrt))
                        loc2=(idLoc(Scrt(2))-idLoc(Scrt(1)))
                        ipDij = 1 + loc1/loc2
                     End  If
                     ipFkl1 = ((i4-1)*iCmpa(3)+i3-1)*kBas*lBas
     &                      + ipFkl
                     If (vkl*vijkl*Abs(Fac_ij).lt.ThrInt .and.
     &                   vij*vijkl*Abs(Fac_kl).lt.ThrInt) Then
                        iOpt = iOpt -1
                     Else
                        lFij=.True.
                        lFkl=.True.
                     End If
                  End If
*
                  If (iAnd(iSym(1),iSym(3)).ne.0 .or.
     &                iAnd(iSym(2),iSym(4)).ne.0) Then
                     iOpt = iOpt + 2
*
                     If (iShell(2).lt.iShell(4)) Then
                        vjl=Djl(jBas*lBas+1,i4,i2)
                        ipDjl=ip
                        ip = ip + jBas*lBas
                        Call DGeTMO(Djl(1,i4,i2),jl1,jl1,
     &                              jl2,Scrt(ipDjl),jl2)
                     Else
                        vjl=Djl(jBas*lBas+1,i2,i4)
                        loc1= (idLoc(Djl(1,i2,i4))-idLoc(Scrt))
                        loc2=(idLoc(Scrt(2))-idLoc(Scrt(1)))
                        ipDjl = 1 + loc1/loc2
                     End If
                     ipFik1 = ((i3-1)*iCmpa(1)+i1-1)*iBas*kBas
     &                      + ipFik
                     If (iShell(1).lt.iShell(3)) Then
                        vik=Dik(iBas*kBas+1,i3,i1)
                        ipDik = ip
                        ip = ip + iBas*kBas
                        Call DGeTMO(Dik(1,i3,i1),ik1,ik1,
     &                              ik2,Scrt(ipDik),ik2)
                     Else
                        vik=Dik(iBas*kBas+1,i1,i3)
                        loc1=(idLoc(Dik(1,i1,i3))-idLoc(Scrt))
                        loc2=(idLoc(Scrt(2))-idLoc(Scrt(1)))
                        ipDik = 1 + loc1/loc2
                     End If
                     ipFjl1 = ((i4-1)*iCmpa(2)+i2-1)*jBas*lBas
     &                      + ipFjl
                     If (vik*vijkl*Abs(Fac_jl) .lt. ThrInt .and.
     &                   vjl*vijkl*Abs(Fac_ik) .lt. ThrInt) Then
                        iOpt = iOpt - 2
                     Else
                        lFik = .True.
                        lFjl = .True.
                     End If
                  End If
*
                  If (iAnd(iSym(1),iSym(4)).ne.0 .or.
     &                iAnd(iSym(2),iSym(3)).ne.0) Then
                     iOpt = iOpt + 4
*
                     If (iShell(2).lt.iShell(3)) Then
                        vjk = Djk(jBas*kBas+1,i3,i2)
                        ipDjk = ip
                        ip = ip + jBas*kBas
                        Call DGeTMO(Djk(1,i3,i2),jk1,jk1,
     &                              jk2,Scrt(ipDjk),jk2)
                     Else
                        vjk = Djk(jBas*kBas+1,i2,i3)
                        loc1=(idLoc(Djk(1,i2,i3))-idLoc(Scrt))
                        loc2= (idLoc(Scrt(2))-idLoc(Scrt(1)))
                        ipDjk = 1 + loc1/loc2
                     End If
                     ipFil1 = ((i4-1)*iCmpa(1)+i1-1)*iBas*lBas
     &                      + ipFil
                     If (iShell(1).lt.iShell(4)) Then
                        vil = Dil(iBas*lBas+1,i4,i1)
                        ipDil = ip
                        ip = ip + iBas*lBas
                        Call DGeTMO(Dil(1,i4,i1),il1,il1,
     &                              il2,Scrt(ipDil),il2)
                     Else
                        vil = Dil(iBas*lBas+1,i1,i4)
                        loc1= (idLoc(Dil(1,i1,i4))-idLoc(Scrt))
                        loc2= (idLoc(Scrt(2))-idLoc(Scrt(1)))
                        ipDil = 1 + loc1/loc2
                     End If
                     ipFjk1 = ((i3-1)*iCmpa(2)+i2-1)*jBas*kBas
     &                     + ipFjk
                     If (vil*vijkl*Abs(Fac_jk) .lt. ThrInt .and.
     &                   vjk*vijkl*Abs(Fac_il) .lt. ThrInt) Then
                        iOpt = iOpt -4
                     Else
                        lFil = .True.
                        lFjk = .True.
                     End If
                  End If
                  If (ip-1.gt.nScrt) Then
                     Write (6,*) 'FckAcc: nScrt too small!'
                     Call Abend
                  End If
                  Go To ( 1, 2, 3, 4, 5, 6, 7) iOpt
                  Go To 400
*                                                                      *
************************************************************************
*                                                                      *
 1                Continue
                  Call Fck1(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,
     &                      Scrt(ipDij),FT(ipFij1),Fac_ij,
     &                      Scrt(ipDkl),FT(ipFkl1),Fac_kl,ExFac)
                  Go To 400
 2                Continue
                  Call Fck2(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,
     &                      Scrt(ipDik),FT(ipFik1),Fac_ik,
     &                      Scrt(ipDjl),FT(ipFjl1),Fac_jl,ExFac)
                  Go To 400
 3                Continue
                  Call Fck3(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,
     &                      Scrt(ipDij),FT(ipFij1),Fac_ij,
     &                      Scrt(ipDkl),FT(ipFkl1),Fac_kl,
     &                      Scrt(ipDik),FT(ipFik1),Fac_ik,
     &                      Scrt(ipDjl),FT(ipFjl1),Fac_jl,ExFac)
                  Go To 400
 4                Continue
                  Call Fck4(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,
     &                      Scrt(ipDil),FT(ipFil1),Fac_il,
     &                      Scrt(ipDjk),FT(ipFjk1),Fac_jk,ExFac)
                  Go To 400
 5                Continue
                  Call Fck5(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,
     &                      Scrt(ipDij),FT(ipFij1),Fac_ij,
     &                      Scrt(ipDkl),FT(ipFkl1),Fac_kl,
     &                      Scrt(ipDil),FT(ipFil1),Fac_il,
     &                      Scrt(ipDjk),FT(ipFjk1),Fac_jk,ExFac)
                  Go To 400
 6                Continue
                  Call Fck6(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,
     &                      Scrt(ipDik),FT(ipFik1),Fac_ik,
     &                      Scrt(ipDjl),FT(ipFjl1),Fac_jl,
     &                      Scrt(ipDil),FT(ipFil1),Fac_il,
     &                      Scrt(ipDjk),FT(ipFjk1),Fac_jk,ExFac)
                  Go To 400
 7                Continue
                  Call Fck7(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,
     &                      Scrt(ipDij),FT(ipFij1),Fac_ij,
     &                      Scrt(ipDkl),FT(ipFkl1),Fac_kl,
     &                      Scrt(ipDik),FT(ipFik1),Fac_ik,
     &                      Scrt(ipDjl),FT(ipFjl1),Fac_jl,
     &                      Scrt(ipDil),FT(ipFil1),Fac_il,
     &                      Scrt(ipDjk),FT(ipFjk1),Fac_jk,ExFac)
                  Go To 400
************************************************************************
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
*
      nnIrrep=nIrrep
      If (sIrrep) nnIrrep=1
      Do iIrrep=0,nnIrrep-1
*
         If (pert(iIrrep)) Then
            ip=ipDisp(abs(indgrd(iCar,iCent,iIrrep)))
            rCh=xPrmt(iOper(kOp(iCent)),iChBas(1+iCar))*
     &                rChTbl(iIrrep,kOp(iCent))
            Fact=tfact*rCh
*           Write (*,*) 'Level ij'
            If (lFij) Call FckDst(TwoHam(ip),ndens,FT(ipFij),
     &                             iBas,jBas,iCmpa(1),iCmpa(2),
     &                             kOp2(1),kOp2(2),
     &                             iIrrep,iShell(1),iShell(2),iShij,
     &                             iAO(1),iAO(2),iAOst(1),iAOst(2),
     &                             fact)
*           Write (*,*) 'Level kl'
            If (lFkl) Call FckDst(TwoHam(ip),ndens,FT(ipFkl),
     &                             kBas,lBas,iCmpa(3),iCmpa(4),
     &                             kOp2(3),kOp2(4),
     &                             iIrrep,iShell(3),iShell(4),iShkl,
     &                             iAO(3),iAO(4),iAOst(3),iAOst(4),
     &                             fact)
*           Write (*,*) 'Level ik'
            If (lFik) Call FckDst(TwoHam(ip),ndens,FT(ipFik),
     &                             iBas,kBas,iCmpa(1),iCmpa(3),
     &                             kOp2(1),kOp2(3),
     &                             iIrrep,iShell(1),iShell(3),iShik,
     &                             iAO(1),iAO(3),iAOst(1),iAOst(3),
     &                             fact)
*           Write (*,*) 'Level jl'
            If (lFjl) Call FckDst(TwoHam(ip),ndens,FT(ipFjl),
     &                             jBas,lBas,iCmpa(2),iCmpa(4),
     &                             kOp2(2),kOp2(4),
     &                             iIrrep,iShell(2),iShell(4),iShjl,
     &                             iAO(2),iAO(4),iAOst(2),iAOst(4),
     &                             fact)
*           Write (*,*) 'Level il'
            If (lFil) Call FckDst(TwoHam(ip),ndens,FT(ipFil),
     &                             iBas,lBas,iCmpa(1),iCmpa(4),
     &                             kOp2(1),kOp2(4),
     &                             iIrrep,iShell(1),iShell(4),iShil,
     &                             iAO(1),iAO(4),iAOst(1),iAOst(4),
     &                             fact)
*           Write (*,*) 'Level jk'
            If (lFjk) Call FckDst(TwoHam(ip),ndens,FT(ipFjk),
     &                             jBas,kBas,iCmpa(2),iCmpa(3),
     &                             kOp2(2),kOp2(3),
     &                             iIrrep,iShell(2),iShell(3),iShjk,
     &                             iAO(2),iAO(3),iAOst(2),iAOst(3),
     &                             fact)
         End If
      End Do
*
      Return
      End
#endif
