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
* Copyright (C) 1990,1991,1994, Roland Lindh                           *
*               1990, IBM                                              *
*               2017, Ignacio Fdez. Galvan                             *
************************************************************************
      SubRoutine Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &               Eta,EInv,nEta,
     &               P,lP,Q,lQ,rKapab,rKapcd,Coori,Coora,CoorAC,
     &               mabMin,mabMax,mcdMin,mcdMax,Array,nArray,
     &               Tvalue,ModU2,Cff2D,Rys2D,NoSpecial)
************************************************************************
*                                                                      *
* Object: to compute the source integrals for the transfer equation    *
*         with the Rys quadrature, i.e. the integrals [e0|f0] will be  *
*         computed (e=Max(la,lb),la+lb,f=Max(lc,ld),lc+ld).            *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              Tvalue                                                  *
*              RtsWgh                                                  *
*              vRysRW                                                  *
*              ModU2                                                   *
*              Cff2D                                                   *
*              Rys2D                                                   *
*              RysEF                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Modified for k2 loop. August '91                         *
*             Modified for decreased memory access January '94         *
*             Modified for special routines Jan-Mar '94                *
************************************************************************
      use vRys_RW
      use Real_Info, only: ChiI2
      use Temporary_parameters, only: IsChi
      Implicit Real*8 (A-H,O-Z)
      External Tvalue, ModU2, Cff2D, Rys2D
#include "notab.fh"
#include "print.fh"
#include "real.fh"
cgh - stuff for short range integrals
#include "FMM.fh"
      logical secondpass
#include "srint.fh"
      Real*8 Zeta(nZeta), ZInv(nZeta), P(lP,3), rKapab(nZeta),
     &       Eta(nEta),   EInv(nEta),  Q(lQ,3), rKapcd(nEta),
     &       CoorAC(3,2), Coora(3,4), Coori(3,4), Array(nArray)
      Integer iAnga(4)
      Logical AeqB, CeqD, AeqC, EQ, NoSpecial
*
*     Statement function for canonical index, etc.
*
      iTri(i,j) = (Max(i,j)*(Max(i,j)-1))/2 + Min(i,j)
*

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
         Write (6,*) 'NoSpecial=',NoSpecial
         Call RecPrt(' In Rys:P','(10G15.5)',P,lP,3)
         Call RecPrt(' In Rys:Q','(10G15.5)',Q,lQ,3)
         Call RecPrt(' In Rys:Zeta','(10G15.5)',Zeta,nZeta,1)
         Call RecPrt(' In Rys:Eta','(10G15.5)',Eta,nEta,1)
         Write (6,*) ' In Rys: iAnga=',iAnga
         Call RecPrt('CoorAC',' ',CoorAC,3,2)
         Call RecPrt('Coora',' ',Coora,3,4)
         Call RecPrt('Coori',' ',Coori,3,4)
         Call RecPrt('rKapab',' ',rKapab,1,nZeta)
         Call RecPrt('rKapcd',' ',rKapcd,1,nEta)
#endif
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      AeqB = EQ(Coori(1,1),Coori(1,2))
      CeqD = EQ(Coori(1,3),Coori(1,4))
      AeqC = EQ(Coori(1,1),Coori(1,3))
      nRys = (la+lb+lc+ld+2)/2
      nabMax = la+lb
      nabMin = Max(la,lb)
      ncdMax = lc+ld
      ncdMin = Max(lc,ld)
      nabcd = (nabMax+1)*(ncdMax+1)
*
*     In some cases a pointer to Array will not be used. However, the
*     subroutine call still have the same number of arguments. In this
*     cases a "dummy pointer" is used. The pointer could point to any
*     odd element in Array. I will use the last one for some unknown reason.
*
      ip_Array_Dummy=nArray
*
      ijkl = 0
      If (NoSpecial) ijkl = -1
*
*     For FMM, compute short-range integrals disabling special cases
cgh - disable special cases anyway for the short range integrals
*
      If (shortrange.or.FMM_shortrange) ijkl = -1
*
      ij = iTri(la+1,lb+1)
      kl = iTri(lc+1,ld+1)
      If (ijkl.eq.0) ijkl = iTri(ij,kl)
      Select Case (ijkl)
*
       Case (1)
*
*--------(ss|ss)
*
         Call ssss(Array,
     &             Zeta,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),
     &             Eta, nEta,Q,lQ,rKapCD,Coori(1,3),Coori(1,4),
     &             Tmax(1),Map(iMap(1)),nMap(1),
     &             x0(ix0(1)),nx0(1),Cff(iCffW(6,1)),
     &             Cff(iCffW(5,1)),Cff(iCffW(4,1)),Cff(iCffW(3,1)),
     &             Cff(iCffW(2,1)),Cff(iCffW(1,1)),Cff(iCffW(0,1)),
     &             ddx(1),HerW2(iHerW2(1)),IsChi,ChiI2)
*
       Case (2)
*
*--------(ps|ss) & (ss|ps)
*
         If (kl.eq.2)
     &   Call sssp(Array,Zeta,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),
     &                    Eta, nEta,Q,lQ,rKapCD,Coori(1,3),Coori(1,4),
     &             CoorAC,
     &             Tmax(1),Map(iMap(1)),nMap(1),
     &             x0(ix0(1)),nx0(1),Cff(iCffW(6,1)),
     &             Cff(iCffW(5,1)),Cff(iCffW(4,1)),Cff(iCffW(3,1)),
     &             Cff(iCffW(2,1)),Cff(iCffW(1,1)),Cff(iCffW(0,1)),
     &             Cff(iCffR(6,1)),
     &             Cff(iCffR(5,1)),Cff(iCffR(4,1)),Cff(iCffR(3,1)),
     &             Cff(iCffR(2,1)),Cff(iCffR(1,1)),Cff(iCffR(0,1)),
     &             ddx(1),HerW2(iHerW2(1)),HerR2(iHerR2(1)),IsChi,ChiI2)
         If (ij.eq.2)
     &   Call psss(Array,
     &             Zeta,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),
     &              Eta, nEta,Q,lQ,rKapCD,Coori(1,3),Coori(1,4),
     &             CoorAC,
     &             Tmax(1),Map(iMap(1)),nMap(1),
     &             x0(ix0(1)),nx0(1),Cff(iCffW(6,1)),
     &             Cff(iCffW(5,1)),Cff(iCffW(4,1)),Cff(iCffW(3,1)),
     &             Cff(iCffW(2,1)),Cff(iCffW(1,1)),Cff(iCffW(0,1)),
     &             Cff(iCffR(6,1)),
     &             Cff(iCffR(5,1)),Cff(iCffR(4,1)),Cff(iCffR(3,1)),
     &             Cff(iCffR(2,1)),Cff(iCffR(1,1)),Cff(iCffR(0,1)),
     &             ddx(1),HerW2(iHerW2(1)),HerR2(iHerR2(1)),IsChi,ChiI2)
*
       Case (3)
*
*--------(ps|ps)
*
         Call psps(Array,Zeta,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),
     &                    Eta, nEta,Q,lQ,rKapCD,Coori(1,3),Coori(1,4),
     &             CoorAC,
     &             Tmax(2),Map(iMap(2)),nMap(2),
     &             x0(ix0(2)),nx0(2),Cff(iCffW(6,2)),
     &             Cff(iCffW(5,2)),Cff(iCffW(4,2)),Cff(iCffW(3,2)),
     &             Cff(iCffW(2,2)),Cff(iCffW(1,2)),Cff(iCffW(0,2)),
     &             Cff(iCffR(6,2)),
     &             Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)),
     &             Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),
     &             ddx(2),HerW2(iHerW2(2)),HerR2(iHerR2(2)),IsChi,ChiI2)
*
       Case (4)
*
*--------(pp|ss) & (ss|pp)
*
      If (ij.eq.3)
     &   Call ppss(Array,Zeta,ZInv,nZeta,P,lP,rKapAB,
     &             Coori(1,1),Coori(1,2),
     &                    Eta,EInv, nEta,Q,lQ,rKapCD,
     &             Coori(1,3),Coori(1,4),
     &             CoorAC,
     &             Tmax(2),Map(iMap(2)),nMap(2),
     &             x0(ix0(2)),nx0(2),Cff(iCffW(6,2)),
     &             Cff(iCffW(5,2)),Cff(iCffW(4,2)),Cff(iCffW(3,2)),
     &             Cff(iCffW(2,2)),Cff(iCffW(1,2)),Cff(iCffW(0,2)),
     &             Cff(iCffR(6,2)),
     &             Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)),
     &             Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),
     &             ddx(2),HerW2(iHerW2(2)),HerR2(iHerR2(2)),IsChi,ChiI2)
      If (kl.eq.3)
     &   Call sspp(Array,Zeta,ZInv,nZeta,P,lP,rKapAB,
     &             Coori(1,1),Coori(1,2),
     &                    Eta,EInv, nEta,Q,lQ,rKapCD,
     &             Coori(1,3),Coori(1,4),
     &             CoorAC,
     &             Tmax(2),Map(iMap(2)),nMap(2),
     &             x0(ix0(2)),nx0(2),Cff(iCffW(6,2)),
     &             Cff(iCffW(5,2)),Cff(iCffW(4,2)),Cff(iCffW(3,2)),
     &             Cff(iCffW(2,2)),Cff(iCffW(1,2)),Cff(iCffW(0,2)),
     &             Cff(iCffR(6,2)),
     &             Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)),
     &             Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),
     &             ddx(2),HerW2(iHerW2(2)),HerR2(iHerR2(2)),IsChi,ChiI2)
*
       Case (5)
*
*--------(pp|ps) & (sp|pp)
*
         If (ij.eq.3)
     &   Call ppps(Array,Zeta,ZInv,nZeta,P,lP,rKapAB,
     &             Coori(1,1),Coori(1,2),
     &                    Eta,EInv, nEta,Q,lQ,rKapCD,
     &             Coori(1,3),Coori(1,4),
     &             CoorAC,
     &             Tmax(2),Map(iMap(2)),nMap(2),
     &             x0(ix0(2)),nx0(2),Cff(iCffW(6,2)),
     &             Cff(iCffW(5,2)),Cff(iCffW(4,2)),Cff(iCffW(3,2)),
     &             Cff(iCffW(2,2)),Cff(iCffW(1,2)),Cff(iCffW(0,2)),
     &             Cff(iCffR(6,2)),
     &             Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)),
     &             Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),
     &             ddx(2),HerW2(iHerW2(2)),HerR2(iHerR2(2)),IsChi,ChiI2)
         If (kl.eq.3)
     &   Call sppp(Array,Zeta,ZInv,nZeta,P,lP,rKapAB,
     &             Coori(1,1),Coori(1,2),
     &                    Eta,EInv, nEta,Q,lQ,rKapCD,
     &             Coori(1,3),Coori(1,4),
     &             CoorAC,
     &             Tmax(2),Map(iMap(2)),nMap(2),
     &             x0(ix0(2)),nx0(2),Cff(iCffW(6,2)),
     &             Cff(iCffW(5,2)),Cff(iCffW(4,2)),Cff(iCffW(3,2)),
     &             Cff(iCffW(2,2)),Cff(iCffW(1,2)),Cff(iCffW(0,2)),
     &             Cff(iCffR(6,2)),
     &             Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)),
     &             Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),
     &             ddx(2),HerW2(iHerW2(2)),HerR2(iHerR2(2)),IsChi,ChiI2)
*
       Case (6)
*
*--------(pp|pp)
*
         Call pppp(Array,Zeta,ZInv,nZeta,P,lP,rKapAB,
     &             Coori(1,1),Coori(1,2),
     &                    Eta,EInv, nEta,Q,lQ,rKapCD,
     &             Coori(1,3),Coori(1,4),
     &             CoorAC,
     &             Tmax(3),Map(iMap(3)),nMap(3),
     &             x0(ix0(3)),nx0(3),Cff(iCffW(6,3)),
     &             Cff(iCffW(5,3)),Cff(iCffW(4,3)),Cff(iCffW(3,3)),
     &             Cff(iCffW(2,3)),Cff(iCffW(1,3)),Cff(iCffW(0,3)),
     &             Cff(iCffR(6,3)),
     &             Cff(iCffR(5,3)),Cff(iCffR(4,3)),Cff(iCffR(3,3)),
     &             Cff(iCffR(2,3)),Cff(iCffR(1,3)),Cff(iCffR(0,3)),
     &             ddx(3),HerW2(iHerW2(3)),HerR2(iHerR2(3)),IsChi,ChiI2)
*
       Case Default
*
*-----General code
*
*     Allocate memory for integrals of [a0|c0] type
      ip = 1
      ipAC = ip
      ip = ip + nT*(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
cgh - in order to produce the short range integrals, two arrays of
cgh - this type are needed - one for the ordinary full range, one for
cgh - the long range integrals
cgh - (additional memory has been declared in MemRys)
      ipAC_long = ipAC
      If (shortrange.or.FMM_shortrange) Then
        ipAC_long = ip
        ip = ip + nT*(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
      EndIf
*     Allocate memory for the normalization factors
      ipFact = ip
      ip = ip + nT
*     Allocate memory for the 2D-integrals.
      ipxyz = ip
      ip = ip + nabcd*3*nT*nRys
      secondpass=.false.
*   - jump mark for second pass:
 2304 continue
*     Allocate memory for the coefficients in the recurrence relation
*     of the 2D-integrals.
      nTR=nT*nRys
      If (NoSpecial) Then
         ipPAQP = ip
      Else
         If (nabMax.ge.1) Then
            iab = 2
            icd = 1
            iabcd = (nabMax+1)*(icd-1) + iab - 1
            ipPAQP = ipxyz + 3*nT*nRys*iabcd
         Else
            ipPAQP = ip_Array_Dummy
         End If
      End If
      ip = ip + nTR*3
      If (NoSpecial) Then
         ipQCPQ = ip
      Else
         If (ncdMax.ge.1) Then
            iab = 1
            icd = 2
            iabcd = (nabMax+1)*(icd-1) + iab - 1
            ipQCPQ = ipxyz + 3*nT*nRys*iabcd
         Else
            ipQCPQ = ip_Array_Dummy
         End If
      End If
      ip = ip + nTR*3
      lB10=Max(Min(nabMax-1,1),0)
      If (lB10.ge.1) Then
         ipB10 = ip
      Else
         ipB10 = ip_Array_Dummy
      End If
      ip = ip + nTR*3*lB10
      labMax = Min(nabMax,ncdMax)
      lB00=Max(Min(labMax,1),0)
      If (lB00.ge.1) Then
         ipB00 = ip
      Else
         ipB00 = ip_Array_Dummy
      End If
      ip = ip + nTR*3*lB00
      lB01=Max(Min(ncdMax-1,1),0)
      If (lB01.ge.1) Then
         ipB01 = ip
      Else
         ipB01 = ip_Array_Dummy
      End If
      ip = ip + nTR*3*lB01
*     Allocate memory for the roots.
      ipU2 = ip
      ip = ip + nT*nRys
*     Allocate memory for Zeta, ZInv, Eta, EInv
      ipZeta = ip
      ip = ip + nT
      ipEta  = ip
      ip = ip + nT
      ipZInv = ip
      ip = ip + nT
      ipEInv = ip
      ip = ip + nT
*     Allocate memory for P and Q
      ipP = ip
      ip = ip + 3*nT
      ipQ = ip
      ip = ip + 3*nT
*     Allocate memory for the inverse.
      ipDiv = ip
      ip = ip + nT
*     Allocate memory for the arguments.
      ipTv = ip
      ip = ip + nT
*     Allocate memory for rKapab and rKapcd
      iprKapab = ip
      ip = ip + nT
      iprKapcd = ip
      ip = ip + nT
*define _CHECK_
#ifdef _CHECK_
      If (ip-1.gt.nArray) Then
         Call WarningMessage(2,'Rys: ip-1 =/= nArray (pos.1)')
         Write (6,*) ' nArray=',nArray
         Write (6,*) ' ip-1  =',ip-1
         Write (6,*) ' nRys  =',nRys
         Write (6,*) ' nZeta =',nZeta
         Write (6,*) ' nEta  =',nEta
         Call Abend()
      End If
#endif
*
*     Expand Zeta, ZInv, Eta ,EInv, rKapab, rKapcd, P, and Q
*
      If (nEta*nZeta.ne.nT) Then
         If (nEta.ne.nT .and.nZeta.ne.nT) Then
            Write (6,*) 'Corrupted parameters!'
            Call Abend()
         End If
         iOff = 0
         call dcopy_(nZeta,Zeta,  1,Array(iOff+ipZeta  ),1)
         call dcopy_(nZeta,ZInv,  1,Array(iOff+ipZInv  ),1)
         call dcopy_(nZeta,rKapab,1,Array(iOff+iprKapab),1)
         call dcopy_(nZeta,P(1,1),1,Array(iOff+ipP     ),1)
         iOff = iOff + nT
         call dcopy_(nZeta,P(1,2),1,Array(iOff+ipP     ),1)
         iOff = iOff + nT
         call dcopy_(nZeta,P(1,3),1,Array(iOff+ipP     ),1)
         iOff = 0
         call dcopy_(nEta, Eta,  1,Array(iOff+ipEta) ,1)
         call dcopy_(nEta,EInv,  1,Array(iOff+ipEInv),1)
         call dcopy_(nEta,rKapcd,1,Array(iOff+iprKapcd),1)
         call dcopy_(nEta,Q(1,1),1,Array(iOff+ipQ     ),1)
         iOff = iOff + nT
         call dcopy_(nEta,Q(1,2),1,Array(iOff+ipQ     ),1)
         iOff = iOff + nT
         call dcopy_(nEta,Q(1,3),1,Array(iOff+ipQ     ),1)
      Else
      Do iEta = 1, nEta
         iOff = (iEta-1)*nZeta
         call dcopy_(nZeta,Zeta,  1,Array(iOff+ipZeta  ),1)
         call dcopy_(nZeta,ZInv,  1,Array(iOff+ipZInv  ),1)
         call dcopy_(nZeta,rKapab,1,Array(iOff+iprKapab),1)
         call dcopy_(nZeta,P(1,1),1,Array(iOff+ipP     ),1)
         iOff = iOff + nZeta*nEta
         call dcopy_(nZeta,P(1,2),1,Array(iOff+ipP     ),1)
         iOff = iOff + nZeta*nEta
         call dcopy_(nZeta,P(1,3),1,Array(iOff+ipP     ),1)
      End Do
      Do iZeta = 1, nZeta
         iOff = iZeta-1
         call dcopy_(nEta, Eta,  1,Array(iOff+ipEta) ,nZeta)
         call dcopy_(nEta,EInv,  1,Array(iOff+ipEInv),nZeta)
         call dcopy_(nEta,rKapcd,1,Array(iOff+iprKapcd),nZeta)
         call dcopy_(nEta,Q(1,1),1,Array(iOff+ipQ     ),nZeta)
         iOff = iOff + nZeta*nEta
         call dcopy_(nEta,Q(1,2),1,Array(iOff+ipQ     ),nZeta)
         iOff = iOff + nZeta*nEta
         call dcopy_(nEta,Q(1,3),1,Array(iOff+ipQ     ),nZeta)
      End Do
      End If
*
*     Compute tha arguments for which we will compute the roots and
*     the weights.
*
      Call Tvalue(Array(ipZeta),Array(ipEta),Array(ipP),Array(ipQ),
     &            Array(iprKapab),Array(iprKapcd),Array(ipTv),
     &            Array(ipFact), Array(ipDiv),nT,IsChi,ChiI2)
*     Let go of rKapab and rKapcd
      ip = ip - 2*nT
*
*     Compute roots and weights. Make sure that the weights ends up in
*     the array where the z component of the 2D integrals will be.
*     Call vRysRW if roots and weights are tabulated in various Taylor
*     expansions. If not tabulated call RtsWgh. If from scratch
*     (no table at all), call RysRtsWgh
*
*     Pointer to z-component of 2D-integrals where the weights will be
*     put directly. This corresponds to xyz2D(1,1,3,0,0).
      ipWgh = ipxyz + 2*nT*nRys
#ifdef _RYS_SCRATCH_
#ifdef _CHECK_
      If (ip-1.gt.nArray) Then
         Call WarningMessage(2,'Rys: ip-1 =/= nArray (pos.2)')
         Write (6,*) ' nArray=',nArray
         Write (6,*) ' ip-1  =',ip-1
         Call Abend()
      End If
#endif
      Call RysRtsWgh(Array(ipTv),nT,Array(ipU2),Array(ipWgh),nRys)
#else
      If (nRys.gt.nMxRys .or. NoTab) Then
#ifdef _CHECK_
         If (ip-1.gt.nArray) Then
            Call WarningMessage(2,'Rys: ip-1 =/= nArray (pos.2)')
            Write (6,*) ' nArray=',nArray
            Write (6,*) ' ip-1  =',ip-1
            Call Abend()
         End If
#endif
         Call RtsWgh(Array(ipTv),nT,Array(ipU2),Array(ipWgh),nRys)
      Else
#ifdef _CHECK_
         If (ip-1.gt.nArray) Then
            Call WarningMessage(2,'Rys: ip-1 =/= nArray (pos.3)')
            Write (6,*) ' nArray=',nArray
            Write (6,*) ' ip-1  =',ip-1
            Call Abend()
         End If
#endif
         Call vRysRW(la,lb,lc,ld,Array(ipTv),Array(ipU2),Array(ipWgh),
     &               nT,nRys)
      End If
#endif
*     Let go of arguments
      ip = ip - nT
*
*     Compute coefficients for the recurrence relations of the
*     2D-integrals
*
      If (la+lb+lc+ld.gt.0) Call ModU2(Array(ipU2),nT,nRys,Array(ipDiv))
*     Let go of inverse
      ip = ip - nT

      Call Cff2D(Max(nabMax-1,0),Max(ncdMax-1,0),nRys,
     &           Array(ipZeta),Array(ipZInv),Array(ipEta),Array(ipEInv),
     &           nT,Coori,CoorAC,Array(ipP),Array(ipQ),la,lb,lc,ld,
     &           Array(ipU2),Array(ipPAQP),Array(ipQCPQ),
     &           Array(ipB10),Array(ipB00),labMax,Array(ipB01))
*     Let go of roots
      ip = ip - nT*nRys
*     Let go of Zeta, ZInv, Eta, and EInv
      ip = ip - nT*4
*     Let go of P and Q
      ip = ip - 6*nT
*
*     Compute the 2D-integrals from the roots and weights
*
      Call Rys2D(Array(ipxyz),nT,nRys,nabMax,
     &           ncdMax,Array(ipPAQP),Array(ipQCPQ),
     &           Array(ipB10),Max(nabMax-1,0),
     &           Array(ipB00),labMax,
     &           Array(ipB01),Max(ncdMax-1,0))
      ip = ip - nTR*3*lB01
      ip = ip - nTR*3*lB00
      ip = ip - nTR*3*lB10
      ip = ip - nTR*3
      ip = ip - nTR*3
*
*     Compute [a0|c0] integrals
*
      ipScr = ip
      ip = ip + nT*nRys
      AeqB = EQ(Coora(1,1),Coora(1,2))
      CeqD = EQ(Coora(1,3),Coora(1,4))
*                                                                      *
************************************************************************
*                                                                      *
*     Use Molpro Coulomb attenuation driver for the
*     FMM short-range integrals
*
      If ((shortrange.and..not.(isr_simulate.gt.1)) .or.
     &     FMM_shortrange) Then
*                                                                      *
************************************************************************
*                                                                      *
         If (.not.secondpass) Then
*
*           [in the first pass, the ordinary full integrals are created
*            in Array(ipScr)]
            Call RysEF(Array(ipxyz),nT,nT,nRys,
     &                 nabMin,nabMax,ncdMin,ncdMax,
     &                 Array(ipAC),mabMin,mabMax,mcdMin,mcdMax,
     &                 Array(ipScr),Array(ipFact),AeqB,CeqD)
*           [release memory at ipScr]
            ip = ip - nT*nRys
*           [in the second pass we will make the long range integrals:]
            If (FMM_shortrange) Then
               asymptotic_Rys = .True.
            Else
               IsChi = 1
            End If
*           [set flag for 2nd pass, then go ahead and do 2nd pass]
            secondpass=.true.
            goto 2304
*
         Else
*
*           [in the second run, the long range integrals are created
*            in Array(ipScr_long)]
            Call RysEF(Array(ipxyz),nT,nT,nRys,
     &                 nabMin,nabMax,ncdMin,ncdMax,
     &                 Array(ipAC_long),mabMin,mabMax,mcdMin,mcdMax,
     &                 Array(ipScr),Array(ipFact),AeqB,CeqD)
*           [make difference to produce the desired short range
*            integrals]
            If((.not.(isr_simulate.gt.0)) .or.
     &         FMM_shortrange) then
              call daxpy_(nT*(mabMax-mabMin+1)*(mcdMax-mcdMin+1),
     &                    -One,Array(ipAC_long),1,Array(ipAC),1)
            End If
*
*           [reset IsChi for ordinary full integrals]
            If (FMM_shortrange) Then
               asymptotic_Rys = .False.
            Else
               IsChi=0
            End If
*
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
         Call RysEF(Array(ipxyz),nT,nT,nRys,
     &              nabMin,nabMax,ncdMin,ncdMax,
     &              Array(ipAC),mabMin,mabMax,mcdMin,mcdMax,
     &              Array(ipScr),Array(ipFact),AeqB,CeqD)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      ip = ip - nT*nRys
      ip = ip - nabcd*3*nT*nRys
      ip = ip - nT
* -   release additional memory allocated for long range integrals
      If (shortrange.or.FMM_shortrange)
     &   ip = ip - nT*(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
      End Select
#ifdef _DEBUGPRINT_
      mabcd=(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
      Call RecPrt('{e0|f0}',' ',Array,nT,mabcd)
#endif
*
      Return
      End
