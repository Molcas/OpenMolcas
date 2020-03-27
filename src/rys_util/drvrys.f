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
* Copyright (C) 2015, Roland Lindh                                     *
************************************************************************
      Subroutine DrvRys(iZeta,iEta,nZeta,nEta,mZeta,mEta,
     &                  nZeta_Tot,nEta_Tot,Data1,mData1,Data2,mData2,
     &                  nAlpha,nBeta,nGamma,nDelta,
     &                  IndZ,Zeta,ZInv,P,KappAB,IndZet,
     &                  IndE, Eta,EInv,Q,KappCD,IndEta,
     &                  ix1,iy1,iz1,ix2,iy2,iz2,ThrInt,CutInt,
     &                  vij,vkl,vik,vil,vjk,vjl,
     &                  Prescreen_On_Int_Only,NoInts,iAnga,Coor,CoorAC,
     &                  mabMin,mabMax,mcdMin,mcdMax,nijkl,nabcd,mabcd,
     &                  Wrk,iW2,iW4,nWork2,mWork2,
     &                  HMtrxAB,HMtrxCD,la,lb,lc,ld,iCmp,iShll,
     &                  NoPInts,Dij,mDij,Dkl,mDkl,Do_TnsCtl,kabcd,
     &                  Coeff1,iBasi,Coeff2,jBasj,
     &                  Coeff3,kBask,Coeff4,lBasl)
************************************************************************
* Routine for the computation of primitive integrals and accumulation  *
* to the (ab|cd) or the (e0|f0) set of integrals. If the primitive     *
* set of integrals are smaller than the set of contracted integrals    *
* the code selects to apply the HRR recursion {e0|f0} -> {ab|cd} here  *
* before the contraction generating the (ab|cd) set of integrals, if   *
* not the {e0|f0} set is contracted to the (e0|f0) set directly and    *
* HRR recursion is applied outside this routine.                       *
*                                                                      *
* For the contraction we have that either all primitive integrals      *
* can be computed in a single step, otherwise subsets of primitive     *
* integrals are computed and acculumulated to the contracted set.      *
*                                                                      *
* The Wrk array is subdiveded into 2 or 3 blocks depending on if the   *
* calling code iterate over subsets of primitive integrals.            *
*                                                                      *
* Memory blocking                                                      *
* ===============                                                      *
* For an interative use:                                               *
*      iW4 points to the start of Wkr, length nWork2-mWork2            *
*      iW2 points at nWork2-mWork+1, length mWork2                     *
*      iW3 points at nWork2, length nWork3                             *
*                                                                      *
* For single iteration use:                                            *
*      iW4 and iW2 points at the start of Wkr, length nWork2           *
*      iW3 points at nWork2, length nWork3                             *
*                                                                      *
* Usage of memory                                                      *
*      Screen: does not use Wrk                                        *
*      Rys:    use iW2 secrtion                                        *
*      HRR:    use the aggregated iW2 and iW3 section                  *
*      Cntrct: use the iW2, iW3, and iW4 sections seperately           *
*                                                                      *
* Author: Roland Lindh                                                 *
*         Dept Chemistry - Angstrom, the Theoretical Chem. Prog.       *
*         Uppsala University, Uppsala, Sweden                          *
*         2015                                                         *
************************************************************************
      Implicit None
      External TERI,ModU2,vCff2D,vRys2D
      External ip_ZtMax, ip_abMax, ip_ZtMaxD, ip_abMaxD
      Integer iZeta,iEta,nZeta,nEta,mZeta,mEta,nZeta_Tot,nEta_Tot,
     &        mData1,mData2,nAlpha,nBeta,nGamma,nDelta,
     &        IndZ(nZeta), IndE(nEta), IndZet(nZeta), IndEta(nEta),
     &        ix1, iy1, iz1, ix2, iy2, iz2, mDij, mDkl, iOffZ, iOffE,
     &        iAnga(4),iCmp(4), iShll(4), la,lb,lc,ld,
     &        mabMin,mabMax,mcdMin,mcdMax,nijkl,nabcd,mabcd,
     &        nWork2, mWork2, iW2, iW4,
     &        iBasi,jBasj,kBask,lBasl, kabcd,
     &        ip_ZtMax, ip_abMax, ip_ZtMaxD, ip_abMaxD
      Real*8  Data1(mData1), Data2(mData2),
     &        Zeta(nZeta), ZInv(nZeta), P(nZeta,3), KappAB(nZeta),
     &         Eta(nEta),  EInv(nEta),  Q(nEta, 3), KappCD(nEta),
     &        ThrInt,CutInt,vij,vkl,vik,vil,vjk,vjl,
     &        Coor(3,4), CoorAC(3,2), Wrk(nWork2),
     &        HMtrxAB(*), HMtrxCD(*), Dij(mDij), Dkl(mDkl),
     &        Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj),
     &        Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl)
      Logical Prescreen_On_Int_Only, NoInts,NoPInts, Do_TnsCtl
*     Local arrays
      Integer lZeta, lEta, i_Int, n1, n2, n3, n4, iW3, nWork3, nW2
      Logical Nospecial
*define _DEBUG_
#ifdef _DEBUG_
      Write (6,*) 'Enter DrvRys'
      Write (6,*) 'iZeta, nZeta, mZeta, nZeta_Tot=',
     &             iZeta, nZeta, mZeta, nZeta_Tot
      Write (6,*) 'iEta , nEta , mEta , nEta_Tot=',
     &             iEta , nEta , mEta , nEta_Tot
      Call RecPrt('Coeff1',' ',Coeff1,nAlpha,iBasi)
      Call RecPrt('Coeff2',' ',Coeff2,nBeta ,jBasj)
      Call RecPrt('Coeff3',' ',Coeff3,nGamma,kBask)
      Call RecPrt('Coeff4',' ',Coeff4,nDelta,lBasl)
      Call RecPrt('KappAB',' ',KappAB,1,nZeta)
      Call RecPrt('KappCD',' ',KappCD,1,nEta)
#endif
*
      NoSpecial=.False. ! Use special code if possible
*
*-----Transfer k2 data and prescreen
*
      iOffZ = mDij-nZeta
      iOffE = mDkl-nEta
      Call Screen(nZeta,nEta,mZeta,mEta,lZeta,lEta,
     &            Zeta,ZInv,P,KappAB,IndZet,
     &            Data1(iZeta),nAlpha,nBeta,IndZ(iZeta),
     &            Data1(ip_ZtMax (nZeta)),
     &            Data1(ip_abMax (nZeta)),
     &            Data1(ip_ZtMaxD(nZeta)),
     &            Data1(ip_abMaxD(nZeta)),
     &            Eta,EInv,Q,KappCD,IndEta,
     &            Data2(iEta),nGamma,nDelta,IndE(iEta),
     &            Data2(ip_ZtMax ( nEta)),
     &            Data2(ip_abMax ( nEta)),
     &            Data2(ip_ZtMaxD( nEta)),
     &            Data2(ip_abMaxD( nEta)),
     &            Dij(iOffZ),Dkl(iOffE),
     &            ix1,iy1,iz1,ix2,iy2,iz2,ThrInt,CutInt,
     &            vij,vkl,vik,vil,vjk,vjl,
     &            Prescreen_On_Int_Only)
*     Write (6,*) 'lZeta,lEta:', lZeta, lEta
      If (lZeta*lEta.eq.0) Then
         Call FZero(Wrk(iW2),mWork2)
         Go To 99
      End If
      NoInts=.False.
*
*-----Compute [a0|c0], ijkl,a,c
*
      Call Rys(iAnga,lZeta*lEta,
     &         Zeta,ZInv,lZeta,
     &         Eta,EInv,lEta,
     &         P,nZeta,Q,nEta,
     &         KappAB,KappCD,
     &         Coor,Coor,CoorAC,
     &         mabMin,mabMax,mcdMin,mcdMax,
     &         Wrk(iW2),mWork2,TERI,ModU2,
     &         vCff2D,vRys2D,NoSpecial)
*
*-----Select between HRR before contraction or to contract
*     and perform the HRR later once the complete set of
*     contracted integrals have been generated.
*
      If (lZeta*lEta.lt.nijkl .and.
     &    mZeta.eq.nZeta_tot .and. mEta.eq.nEta_tot ) Then
*
*--------Apply the HRR recusions first. Note that this is only
*        executed if used in single iteration mode. Hence,
*        iW2 and iW4 are identical.
*
         iW3=iW2+lZeta*lEta*mabcd
         Call DGeTMO(Wrk(iW2),lZeta*lEta,lZeta*lEta,
     &               mabcd,Wrk(iW3),mabcd)
         call dcopy_(mabcd*lZeta*lEta,Wrk(iW3),1,Wrk(iW2),1)
         Call TnsCtl(Wrk(iW2),nWork2,Coor,
     &               mabcd,lZeta*lEta,
     &               mabMax,mabMin,mcdMax,mcdMin,
     &               HMtrxAB,HMtrxCD,la,lb,lc,ld,
     &               iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &               iShll(1),iShll(2),iShll(3),iShll(4),i_Int)
         If (i_Int.ne.iW2) call dcopy_(lZeta*lEta*nabcd,Wrk(i_int),1,
     &                                                 Wrk(iW2),1)
         Do_TnsCtl=.False.
         n1=1
         n2=iCmp(1)*iCmp(2)
         n3=1
         n4=iCmp(3)*iCmp(4)
         kabcd=nabcd
      Else
*
*--------Postpone application of the HRR recusions until later.
*
         Do_TnsCtl=.True.
         n1=mabMin
         n2=mabMax
         n3=mcdMin
         n4=mcdMax
         kabcd=mabcd
      End If
*
*-----Accumulate to the contracted integrals
*
      If (iW4.ne.iW2) Then
*        Account for size of the integrals in
         nW2=lZeta*lEta * kabcd
      Else ! iW4.eq.iW2
*        Account for size of the integrals in and out
         nW2=Max(iBasi*jBasj*kBask*lBasl,lZeta*lEta) * kabcd
      End If
      iW3=iW2+nW2
      nWork3=mWork2-nW2
*     Write (6,*) 'iW4,iW2,iW3:',iW4,iW2,iW3
*     Write (6,*) 'nWork3:',nWork3
      Call Cntrct(NoPInts,
     &            Coeff1,nAlpha,iBasi,Coeff2,nBeta ,jBasj,
     &            Coeff3,nGamma,kBask,Coeff4,nDelta,lBasl,
     &            Wrk(iW2),n1,n2,n3,n4,Wrk(iW3),nWork3,Wrk(iW4),
     &            IndZet,lZeta,IndEta,lEta)
*
 99   Continue
#ifdef _DEBUG_
      Write (6,*) 'iW4,iW2,iW3:',iW4,iW2,iW3
      Call RecPrt('DrvRys:(e0|0f)',' ',Wrk(iW4),kabcd,
     &            iBasi*jBasj*kBask*lBasl)
#endif
*
      Return
      End
