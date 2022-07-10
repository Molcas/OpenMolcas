!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990,1991,1994, Roland Lindh                           *
!               1990, IBM                                              *
!               2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Rys(iAnga,nT,Zeta,ZInv,nZeta,Eta,EInv,nEta,P,lP,Q,lQ,rKapab,rKapcd,Coori,Coora,CoorAC,mabMin,mabMax,mcdMin,mcdMax, &
               Array,nArray,Tvalue,ModU2,Cff2D,Rys2D,NoSpecial)
!***********************************************************************
!                                                                      *
! Object: to compute the source integrals for the transfer equation    *
!         with the Rys quadrature, i.e. the integrals [e0|f0] will be  *
!         computed (e=Max(la,lb),la+lb,f=Max(lc,ld),lc+ld).            *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Modified for k2 loop. August '91                         *
!             Modified for decreased memory access January '94         *
!             Modified for special routines Jan-Mar '94                *
!***********************************************************************

use vRys_RW, only: Cff, ddx, HerR2, HerW2, iCffR, iCffW, iHerR2, iHerW2, iMap, ix0, Map, nMap, nx0, TMax, x0
use Gateway_Info, only: ChiI2
use Gateway_global, only: asymptotic_Rys, FMM_shortrange, IsChi, NoTab
#ifdef _RYS_SCRATCH_
use RysScratch, only: RysRtsWgh
#else
use vRys_RW, only: nMxRys
#endif
use Index_Functions, only: iTri
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iAnga(4), nT, nZeta, nEta, lP, lQ, mabMin, mabMax, mcdMin, mcdMax, nArray
real(kind=wp), intent(in) :: Zeta(nZeta), ZInv(nZeta), Eta(nEta), EInv(nEta), P(lP,3), Q(lQ,3), rKapab(nZeta), rKapcd(nEta), &
                             Coori(3,4), Coora(3,4), CoorAC(3,2)
real(kind=wp), intent(inout) :: Array(nArray)
external :: Tvalue, ModU2, Cff2D, Rys2D
logical(kind=iwp), intent(in) :: NoSpecial
integer(kind=iwp) :: iab, iabcd, icd, iEta, ij, ijkl, iOff, ip, ip_Array_Dummy, ipAC, ipAC_long, ipB00, ipB01, ipB10, ipDiv, &
                     ipEInv, ipEta, ipFact, ipP, ipPAQP, ipQ, ipQCPQ, iprKapab, iprKapcd, ipScr, ipTv, ipU2, ipWgh, ipxyz, ipZeta, &
                     ipZInv, iZeta, kl, la, labMax, lb, lB00, lB01, lB10, lc, ld, nabcd, nabMax, nabMin, ncdMax, ncdMin, nRys, &
                     ntmp, nTR
logical(kind=iwp) :: AeqB, CeqD, secondpass
logical(kind=iwp), external :: EQ

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(u6,*) 'NoSpecial=',NoSpecial
call RecPrt(' In Rys:P','(10G15.5)',P,lP,3)
call RecPrt(' In Rys:Q','(10G15.5)',Q,lQ,3)
call RecPrt(' In Rys:Zeta','(10G15.5)',Zeta,nZeta,1)
call RecPrt(' In Rys:Eta','(10G15.5)',Eta,nEta,1)
write(u6,*) ' In Rys: iAnga=',iAnga
call RecPrt('CoorAC',' ',CoorAC,3,2)
call RecPrt('Coora',' ',Coora,3,4)
call RecPrt('Coori',' ',Coori,3,4)
call RecPrt('rKapab',' ',rKapab,1,nZeta)
call RecPrt('rKapcd',' ',rKapcd,1,nEta)
#endif
la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)
AeqB = EQ(Coori(1,1),Coori(1,2))
CeqD = EQ(Coori(1,3),Coori(1,4))
nRys = (la+lb+lc+ld+2)/2
nabMax = la+lb
nabMin = max(la,lb)
ncdMax = lc+ld
ncdMin = max(lc,ld)
nabcd = (nabMax+1)*(ncdMax+1)

! In some cases a pointer to Array will not be used. However, the
! subroutine call still have the same number of arguments. In this
! cases a "dummy pointer" is used. The pointer could point to any
! odd element in Array. I will use the last one for some unknown reason.

ip_Array_Dummy = nArray

ijkl = 0
if (NoSpecial) ijkl = -1

! For FMM, compute short-range integrals disabling special cases
!gh - disable special cases anyway for the short range integrals

if (FMM_shortrange) ijkl = -1

ij = iTri(la+1,lb+1)
kl = iTri(lc+1,ld+1)
if (ijkl == 0) ijkl = iTri(ij,kl)
select case (ijkl)

  case (1)

    ! (ss|ss)

    call ssss(Array,Zeta,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),Eta,nEta,Q,lQ,rKapCD,Coori(1,3),Coori(1,4),Tmax(1),Map(iMap(1)), &
              nMap(1),x0(ix0(1)),nx0(1),Cff(iCffW(6,1)),Cff(iCffW(5,1)),Cff(iCffW(4,1)),Cff(iCffW(3,1)),Cff(iCffW(2,1)), &
              Cff(iCffW(1,1)),Cff(iCffW(0,1)),ddx(1),HerW2(iHerW2(1)),IsChi,ChiI2)

  case (2)

    ! (ps|ss) & (ss|ps)

    if (kl == 2) call sssp(Array,Zeta,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),Eta,nEta,Q,lQ,rKapCD,Coori(1,3),Coori(1,4),CoorAC, &
                           Tmax(1),Map(iMap(1)),nMap(1),x0(ix0(1)),nx0(1),Cff(iCffW(6,1)),Cff(iCffW(5,1)),Cff(iCffW(4,1)), &
                           Cff(iCffW(3,1)),Cff(iCffW(2,1)),Cff(iCffW(1,1)),Cff(iCffW(0,1)),Cff(iCffR(6,1)),Cff(iCffR(5,1)), &
                           Cff(iCffR(4,1)),Cff(iCffR(3,1)),Cff(iCffR(2,1)),Cff(iCffR(1,1)),Cff(iCffR(0,1)),ddx(1), &
                           HerW2(iHerW2(1)),HerR2(iHerR2(1)),IsChi,ChiI2)
    if (ij == 2) call psss(Array,Zeta,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),Eta,nEta,Q,lQ,rKapCD,Coori(1,3),Coori(1,4),CoorAC, &
                           Tmax(1),Map(iMap(1)),nMap(1),x0(ix0(1)),nx0(1),Cff(iCffW(6,1)),Cff(iCffW(5,1)),Cff(iCffW(4,1)), &
                           Cff(iCffW(3,1)),Cff(iCffW(2,1)),Cff(iCffW(1,1)),Cff(iCffW(0,1)),Cff(iCffR(6,1)),Cff(iCffR(5,1)), &
                           Cff(iCffR(4,1)),Cff(iCffR(3,1)),Cff(iCffR(2,1)),Cff(iCffR(1,1)),Cff(iCffR(0,1)),ddx(1), &
                           HerW2(iHerW2(1)),HerR2(iHerR2(1)),IsChi,ChiI2)

  case (3)

    ! (ps|ps)

    call psps(Array,Zeta,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),Eta,nEta,Q,lQ,rKapCD,Coori(1,3),Coori(1,4),CoorAC,Tmax(2), &
              Map(iMap(2)),nMap(2),x0(ix0(2)),nx0(2),Cff(iCffW(6,2)),Cff(iCffW(5,2)),Cff(iCffW(4,2)),Cff(iCffW(3,2)), &
              Cff(iCffW(2,2)),Cff(iCffW(1,2)),Cff(iCffW(0,2)),Cff(iCffR(6,2)),Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)), &
              Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),ddx(2),HerW2(iHerW2(2)),HerR2(iHerR2(2)),IsChi,ChiI2)

  case (4)

    ! (pp|ss) & (ss|pp)

    if (ij == 3) call ppss(Array,Zeta,ZInv,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),Eta,EInv,nEta,Q,lQ,rKapCD,Coori(1,3), &
                           Coori(1,4),CoorAC,Tmax(2),Map(iMap(2)),nMap(2),x0(ix0(2)),nx0(2),Cff(iCffW(6,2)),Cff(iCffW(5,2)), &
                           Cff(iCffW(4,2)),Cff(iCffW(3,2)),Cff(iCffW(2,2)),Cff(iCffW(1,2)),Cff(iCffW(0,2)),Cff(iCffR(6,2)), &
                           Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)),Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),ddx(2), &
                           HerW2(iHerW2(2)),HerR2(iHerR2(2)),IsChi,ChiI2)
    if (kl == 3) call sspp(Array,Zeta,ZInv,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),Eta,EInv,nEta,Q,lQ,rKapCD,Coori(1,3), &
                           Coori(1,4),CoorAC,Tmax(2),Map(iMap(2)),nMap(2),x0(ix0(2)),nx0(2),Cff(iCffW(6,2)),Cff(iCffW(5,2)), &
                           Cff(iCffW(4,2)),Cff(iCffW(3,2)),Cff(iCffW(2,2)),Cff(iCffW(1,2)),Cff(iCffW(0,2)),Cff(iCffR(6,2)), &
                           Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)),Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),ddx(2), &
                           HerW2(iHerW2(2)),HerR2(iHerR2(2)),IsChi,ChiI2)

  case (5)

    ! (pp|ps) & (sp|pp)

    if (ij == 3) call ppps(Array,Zeta,ZInv,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),Eta,EInv,nEta,Q,lQ,rKapCD,Coori(1,3), &
                           Coori(1,4),CoorAC,Tmax(2),Map(iMap(2)),nMap(2),x0(ix0(2)),nx0(2),Cff(iCffW(6,2)),Cff(iCffW(5,2)), &
                           Cff(iCffW(4,2)),Cff(iCffW(3,2)),Cff(iCffW(2,2)),Cff(iCffW(1,2)),Cff(iCffW(0,2)),Cff(iCffR(6,2)), &
                           Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)),Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),ddx(2), &
                           HerW2(iHerW2(2)),HerR2(iHerR2(2)),IsChi,ChiI2)
    if (kl == 3) call sppp(Array,Zeta,ZInv,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),Eta,EInv,nEta,Q,lQ,rKapCD,Coori(1,3), &
                           Coori(1,4),CoorAC,Tmax(2),Map(iMap(2)),nMap(2),x0(ix0(2)),nx0(2),Cff(iCffW(6,2)),Cff(iCffW(5,2)), &
                           Cff(iCffW(4,2)),Cff(iCffW(3,2)),Cff(iCffW(2,2)),Cff(iCffW(1,2)),Cff(iCffW(0,2)),Cff(iCffR(6,2)), &
                           Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)),Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),ddx(2), &
                           HerW2(iHerW2(2)),HerR2(iHerR2(2)),IsChi,ChiI2)

  case (6)

    ! (pp|pp)

    call pppp(Array,Zeta,ZInv,nZeta,P,lP,rKapAB,Coori(1,1),Coori(1,2),Eta,EInv,nEta,Q,lQ,rKapCD,Coori(1,3),Coori(1,4),CoorAC, &
              Tmax(3),Map(iMap(3)),nMap(3),x0(ix0(3)),nx0(3),Cff(iCffW(6,3)),Cff(iCffW(5,3)),Cff(iCffW(4,3)),Cff(iCffW(3,3)), &
              Cff(iCffW(2,3)),Cff(iCffW(1,3)),Cff(iCffW(0,3)),Cff(iCffR(6,3)),Cff(iCffR(5,3)),Cff(iCffR(4,3)),Cff(iCffR(3,3)), &
              Cff(iCffR(2,3)),Cff(iCffR(1,3)),Cff(iCffR(0,3)),ddx(3),HerW2(iHerW2(3)),HerR2(iHerR2(3)),IsChi,ChiI2)

  case default

    ! General code

    ! Allocate memory for integrals of [a0|c0] type
    ip = 1
    ipAC = ip
    ip = ip+nT*(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
    !gh - in order to produce the short range integrals, two arrays of
    !gh - this type are needed - one for the ordinary full range, one for
    !gh - the long range integrals
    !gh - (additional memory has been declared in MemRys)
    ipAC_long = ipAC
    if (FMM_shortrange) then
      ipAC_long = ip
      ip = ip+nT*(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
    end if
    ! Allocate memory for the normalization factors
    ipFact = ip
    ip = ip+nT
    ! Allocate memory for the 2D-integrals.
    ipxyz = ip
    ip = ip+nabcd*3*nT*nRys
    secondpass = .false.
    ! jump mark for second pass:
    do
      ! Allocate memory for the coefficients in the recurrence relation
      ! of the 2D-integrals.
      nTR = nT*nRys
      if (NoSpecial) then
        ipPAQP = ip
      else
        if (nabMax >= 1) then
          iab = 2
          icd = 1
          iabcd = (nabMax+1)*(icd-1)+iab-1
          ipPAQP = ipxyz+3*nT*nRys*iabcd
        else
          ipPAQP = ip_Array_Dummy
        end if
      end if
      ip = ip+nTR*3
      if (NoSpecial) then
        ipQCPQ = ip
      else
        if (ncdMax >= 1) then
          iab = 1
          icd = 2
          iabcd = (nabMax+1)*(icd-1)+iab-1
          ipQCPQ = ipxyz+3*nT*nRys*iabcd
        else
          ipQCPQ = ip_Array_Dummy
        end if
      end if
      ip = ip+nTR*3
      lB10 = max(min(nabMax-1,1),0)
      if (lB10 >= 1) then
        ipB10 = ip
      else
        ipB10 = ip_Array_Dummy
      end if
      ip = ip+nTR*3*lB10
      labMax = min(nabMax,ncdMax)
      lB00 = max(min(labMax,1),0)
      if (lB00 >= 1) then
        ipB00 = ip
      else
        ipB00 = ip_Array_Dummy
      end if
      ip = ip+nTR*3*lB00
      lB01 = max(min(ncdMax-1,1),0)
      if (lB01 >= 1) then
        ipB01 = ip
      else
        ipB01 = ip_Array_Dummy
      end if
      ip = ip+nTR*3*lB01
      ! Allocate memory for the roots.
      ipU2 = ip
      ip = ip+nT*nRys
      ! Allocate memory for Zeta, ZInv, Eta, EInv
      ipZeta = ip
      ip = ip+nT
      ipEta = ip
      ip = ip+nT
      ipZInv = ip
      ip = ip+nT
      ipEInv = ip
      ip = ip+nT
      ! Allocate memory for P and Q
      ipP = ip
      ip = ip+3*nT
      ipQ = ip
      ip = ip+3*nT
      ! Allocate memory for the inverse.
      ipDiv = ip
      ip = ip+nT
      ! Allocate memory for the arguments.
      ipTv = ip
      ip = ip+nT
      ! Allocate memory for rKapab and rKapcd
      iprKapab = ip
      ip = ip+nT
      iprKapcd = ip
      ip = ip+nT
!#     define _CHECK_
#     ifdef _CHECK_
      if (ip-1 > nArray) then
        call WarningMessage(2,'Rys: ip-1 =/= nArray (pos.1)')
        write(u6,*) ' nArray=',nArray
        write(u6,*) ' ip-1  =',ip-1
        write(u6,*) ' nRys  =',nRys
        write(u6,*) ' nZeta =',nZeta
        write(u6,*) ' nEta  =',nEta
        call Abend()
      end if
#     endif

      ! Expand Zeta, ZInv, Eta ,EInv, rKapab, rKapcd, P, and Q

      if (nEta*nZeta /= nT) then
        if ((nEta /= nT) .and. (nZeta /= nT)) then
          write(u6,*) 'Corrupted parameters!'
          call Abend()
        end if
        iOff = 0
        call dcopy_(nZeta,Zeta,1,Array(iOff+ipZeta),1)
        call dcopy_(nZeta,ZInv,1,Array(iOff+ipZInv),1)
        call dcopy_(nZeta,rKapab,1,Array(iOff+iprKapab),1)
        call dcopy_(nZeta,P(1,1),1,Array(iOff+ipP),1)
        iOff = iOff+nT
        call dcopy_(nZeta,P(1,2),1,Array(iOff+ipP),1)
        iOff = iOff+nT
        call dcopy_(nZeta,P(1,3),1,Array(iOff+ipP),1)
        iOff = 0
        call dcopy_(nEta,Eta,1,Array(iOff+ipEta),1)
        call dcopy_(nEta,EInv,1,Array(iOff+ipEInv),1)
        call dcopy_(nEta,rKapcd,1,Array(iOff+iprKapcd),1)
        call dcopy_(nEta,Q(1,1),1,Array(iOff+ipQ),1)
        iOff = iOff+nT
        call dcopy_(nEta,Q(1,2),1,Array(iOff+ipQ),1)
        iOff = iOff+nT
        call dcopy_(nEta,Q(1,3),1,Array(iOff+ipQ),1)
      else
        do iEta=1,nEta
          iOff = (iEta-1)*nZeta
          call dcopy_(nZeta,Zeta,1,Array(iOff+ipZeta),1)
          call dcopy_(nZeta,ZInv,1,Array(iOff+ipZInv),1)
          call dcopy_(nZeta,rKapab,1,Array(iOff+iprKapab),1)
          call dcopy_(nZeta,P(1,1),1,Array(iOff+ipP),1)
          iOff = iOff+nZeta*nEta
          call dcopy_(nZeta,P(1,2),1,Array(iOff+ipP),1)
          iOff = iOff+nZeta*nEta
          call dcopy_(nZeta,P(1,3),1,Array(iOff+ipP),1)
        end do
        do iZeta=1,nZeta
          iOff = iZeta-1
          call dcopy_(nEta,Eta,1,Array(iOff+ipEta),nZeta)
          call dcopy_(nEta,EInv,1,Array(iOff+ipEInv),nZeta)
          call dcopy_(nEta,rKapcd,1,Array(iOff+iprKapcd),nZeta)
          call dcopy_(nEta,Q(1,1),1,Array(iOff+ipQ),nZeta)
          iOff = iOff+nZeta*nEta
          call dcopy_(nEta,Q(1,2),1,Array(iOff+ipQ),nZeta)
          iOff = iOff+nZeta*nEta
          call dcopy_(nEta,Q(1,3),1,Array(iOff+ipQ),nZeta)
        end do
      end if

      ! Compute the arguments for which we will compute the roots and the weights.

      call Tvalue(Array(ipZeta),Array(ipEta),Array(ipP),Array(ipQ),Array(iprKapab),Array(iprKapcd),Array(ipTv),Array(ipFact), &
                  Array(ipDiv),nT,IsChi,ChiI2)
      ! Let go of rKapab and rKapcd
      ip = ip-2*nT

      ! Compute roots and weights. Make sure that the weights ends up in
      ! the array where the z component of the 2D integrals will be.
      ! Call vRysRW if roots and weights are tabulated in various Taylor
      ! expansions. If not tabulated call RtsWgh. If from scratch
      ! (no table at all), call RysRtsWgh

      ! Pointer to z-component of 2D-integrals where the weights will be
      ! put directly. This corresponds to xyz2D(1,1,3,0,0).
      ipWgh = ipxyz+2*nT*nRys
#     ifdef _RYS_SCRATCH_
#     ifdef _CHECK_
      if (ip-1 > nArray) then
        call WarningMessage(2,'Rys: ip-1 =/= nArray (pos.2)')
        write(u6,*) ' nArray=',nArray
        write(u6,*) ' ip-1  =',ip-1
        call Abend()
      end if
#     endif
      call RysRtsWgh(Array(ipTv),nT,Array(ipU2),Array(ipWgh),nRys)
#     else
      if ((nRys > nMxRys) .or. NoTab) then
#       ifdef _CHECK_
        if (ip-1 > nArray) then
          call WarningMessage(2,'Rys: ip-1 =/= nArray (pos.2)')
          write(u6,*) ' nArray=',nArray
          write(u6,*) ' ip-1  =',ip-1
          call Abend()
        end if
#       endif
        call RtsWgh(Array(ipTv),nT,Array(ipU2),Array(ipWgh),nRys)
      else
#       ifdef _CHECK_
        if (ip-1 > nArray) then
          call WarningMessage(2,'Rys: ip-1 =/= nArray (pos.3)')
          write(u6,*) ' nArray=',nArray
          write(u6,*) ' ip-1  =',ip-1
          call Abend()
        end if
#       endif
        call vRysRW(la,lb,lc,ld,Array(ipTv),Array(ipU2),Array(ipWgh),nT,nRys)
      end if
#     endif
      ! Let go of arguments
      ip = ip-nT

      ! Compute coefficients for the recurrence relations of the 2D-integrals

      if (la+lb+lc+ld > 0) call ModU2(Array(ipU2),nT,nRys,Array(ipDiv))
      ! Let go of inverse
      ip = ip-nT

      call Cff2D(max(nabMax-1,0),max(ncdMax-1,0),nRys,Array(ipZeta),Array(ipZInv),Array(ipEta),Array(ipEInv),nT,Coori,CoorAC, &
                 Array(ipP),Array(ipQ),la,lb,lc,ld,Array(ipU2),Array(ipPAQP),Array(ipQCPQ),Array(ipB10),Array(ipB00),labMax, &
                 Array(ipB01))
      ! Let go of roots
      ip = ip-nT*nRys
      ! Let go of Zeta, ZInv, Eta, and EInv
      ip = ip-nT*4
      ! Let go of P and Q
      ip = ip-6*nT

      ! Compute the 2D-integrals from the roots and weights

      call Rys2D(Array(ipxyz),nT,nRys,nabMax,ncdMax,Array(ipPAQP),Array(ipQCPQ),Array(ipB10),Array(ipB00),Array(ipB01))
      ip = ip-nTR*3*lB01
      ip = ip-nTR*3*lB00
      ip = ip-nTR*3*lB10
      ip = ip-nTR*3
      ip = ip-nTR*3

      ! Compute [a0|c0] integrals

      ipScr = ip
      ip = ip+nT*nRys
      AeqB = EQ(Coora(1,1),Coora(1,2))
      CeqD = EQ(Coora(1,3),Coora(1,4))
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Use Molpro Coulomb attenuation driver for the
      ! FMM short-range integrals

      if (FMM_shortrange) then
        !                                                              *
        !***************************************************************
        !                                                              *
        if (.not. secondpass) then

          ! [in the first pass, the ordinary full integrals are created in Array(ipScr)]
          call RysEF(Array(ipxyz),nT,nT,nRys,nabMin,nabMax,ncdMin,ncdMax,Array(ipAC),mabMin,mabMax,mcdMin,mcdMax,Array(ipScr), &
                     Array(ipFact),AeqB,CeqD)
          ! [release memory at ipScr]
          ip = ip-nT*nRys
          ! [in the second pass we will make the long range integrals:]
          if (FMM_shortrange) then
            asymptotic_Rys = .true.
          else
            IsChi = 1
          end if
          ! [set flag for 2nd pass, then go ahead and do 2nd pass]
          secondpass = .true.

        else

          ! [in the second run, the long range integrals are created in Array(ipScr_long)]
          call RysEF(Array(ipxyz),nT,nT,nRys,nabMin,nabMax,ncdMin,ncdMax,Array(ipAC_long),mabMin,mabMax,mcdMin,mcdMax, &
                     Array(ipScr),Array(ipFact),AeqB,CeqD)
          ! [make difference to produce the desired short range integrals]
          if (FMM_shortrange) then
            ntmp = nT*(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
            Array(ipAC:ipAC+ntmp-1) = Array(ipAC:ipAC+ntmp-1)-Array(ipAC_long:ipAC_long+ntmp-1)
          end if

          ! [reset IsChi for ordinary full integrals]
          if (FMM_shortrange) then
            asymptotic_Rys = .false.
          else
            IsChi = 0
          end if
          exit

        end if
        !                                                              *
        !***************************************************************
        !                                                              *
      else
        !                                                              *
        !***************************************************************
        !                                                              *
        call RysEF(Array(ipxyz),nT,nT,nRys,nabMin,nabMax,ncdMin,ncdMax,Array(ipAC),mabMin,mabMax,mcdMin,mcdMax,Array(ipScr), &
                   Array(ipFact),AeqB,CeqD)
        exit
        !                                                              *
        !***************************************************************
        !                                                              *
      end if
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ip = ip-nT*nRys
    ip = ip-nabcd*3*nT*nRys
    ip = ip-nT
    ! - release additional memory allocated for long range integrals
    if (FMM_shortrange) ip = ip-nT*(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
end select
#ifdef _DEBUGPRINT_
mabcd = (mabMax-mabMin+1)*(mcdMax-mcdMin+1)
call RecPrt('{e0|f0}',' ',Array,nT,mabcd)
#endif

return

end subroutine Rys
