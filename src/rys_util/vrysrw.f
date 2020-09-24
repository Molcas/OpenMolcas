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
*               1990, IBM                                              *
************************************************************************
      SubRoutine vRysRW(la,lb,lc,ld,Arg,Root,Weight,nArg,nRys)
************************************************************************
*                                                                      *
*  Object: to compute the roots and weights of the Rys polynomials.    *
*          This is done with two approximations. For low arguments     *
*          we will use a 6'th order polynomial and for high arguments  *
*          we will use the asymptotic formulas which are based on the  *
*          roots and weight of Hermite polynomials.                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             September '90                                            *
************************************************************************
      use vRys_RW
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "FMM.fh"
      Real*8 Arg(nArg), Root(nRys,nArg), Weight(nRys,nArg), Tmax_
*
      iRout = 78
      iPrint = nPrint(iRout)
*     Call qEnter('vRysRW')
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Call RecPrt('In vRysRW:Arg',' ',Arg,nArg,1)
#endif
      labcd=1
*
      If (nRys.gt.nMxRys) Then
          Call WarningMessage(2,
     &            'vRysrw: nRys in vRysRW is larger than nMxRys!')
          Write (6,*) ' nRys  =',nRys
          Write (6,*) ' nMxRys=',nMxRys
          Call Abend()
      End If
*
*  For the FMM we use the asymptotic limit to compute the
*  multipole-component of the integrals
*
      TMax_=TMax(nRys)
      If (asymptotic_Rys) TMax_=1.0D99
*
      If (nRys.eq.1) Then
      labcd=la+lb+lc+ld
      If (labcd.eq.0) Then
         Call Rys01(Arg,nArg,Weight,Map(iMap(1)),nMap(1),
     &              x0(ix0(1)),nx0(1),Cff(iCffW(6,1)),
     &              Cff(iCffW(5,1)),Cff(iCffW(4,1)),Cff(iCffW(3,1)),
     &              Cff(iCffW(2,1)),
     &              Cff(iCffW(1,1)),Cff(iCffW(0,1)),ddx(nRys),
     &              HerW2(iHerW2(1)),TMax_)
      Else
         Call Rys11(Arg,nArg,Root,Weight,
     &              Map(iMap(1)),nMap(1),
     &              x0(ix0(1)),nx0(1),Cff(iCffR(6,1)),
     &              Cff(iCffR(5,1)),Cff(iCffR(4,1)),Cff(iCffR(3,1)),
     &              Cff(iCffR(2,1)),Cff(iCffR(1,1)),Cff(iCffR(0,1)),
     &              Cff(iCffW(6,1)),Cff(iCffW(5,1)),Cff(iCffW(4,1)),
     &              Cff(iCffW(3,1)),Cff(iCffW(2,1)),Cff(iCffW(1,1)),
     &              Cff(iCffW(0,1)),ddx(nRys),
     &              HerW2(iHerW2(1)),HerR2(iHerR2(1)),TMax_)
      End If
*
      Else If (nRys.eq.2) Then
      Call Rys22(Arg,nArg,Root,Weight,
     &           Map(iMap(2)),nMap(2),
     &           x0(ix0(2)),nx0(2),Cff(iCffR(6,2)),
     &           Cff(iCffR(5,2)),Cff(iCffR(4,2)),Cff(iCffR(3,2)),
     &           Cff(iCffR(2,2)),Cff(iCffR(1,2)),Cff(iCffR(0,2)),
     &           Cff(iCffW(6,2)),Cff(iCffW(5,2)),Cff(iCffW(4,2)),
     &           Cff(iCffW(3,2)),Cff(iCffW(2,2)),Cff(iCffW(1,2)),
     &           Cff(iCffW(0,2)),ddx(nRys),
     &           HerW2(iHerW2(2)),HerR2(iHerR2(2)),TMax_)
*
      Else If (nRys.eq.3) Then
      Call Rys33(Arg,nArg,Root,Weight,
     &           Map(iMap(3)),nMap(3),
     &           x0(ix0(3)),nx0(3),Cff(iCffR(6,3)),
     &           Cff(iCffR(5,3)),Cff(iCffR(4,3)),Cff(iCffR(3,3)),
     &           Cff(iCffR(2,3)),Cff(iCffR(1,3)),Cff(iCffR(0,3)),
     &           Cff(iCffW(6,3)),Cff(iCffW(5,3)),Cff(iCffW(4,3)),
     &           Cff(iCffW(3,3)),Cff(iCffW(2,3)),Cff(iCffW(1,3)),
     &           Cff(iCffW(0,3)),ddx(nRys),
     &           HerW2(iHerW2(3)),HerR2(iHerR2(3)),TMax_)
*
      Else If (nRys.eq.4) Then
      Call Rys44(Arg,nArg,Root,Weight,
     &           Map(iMap(4)),nMap(4),
     &           x0(ix0(4)),nx0(4),Cff(iCffR(6,4)),
     &           Cff(iCffR(5,4)),Cff(iCffR(4,4)),Cff(iCffR(3,4)),
     &           Cff(iCffR(2,4)),Cff(iCffR(1,4)),Cff(iCffR(0,4)),
     &           Cff(iCffW(6,4)),Cff(iCffW(5,4)),Cff(iCffW(4,4)),
     &           Cff(iCffW(3,4)),Cff(iCffW(2,4)),Cff(iCffW(1,4)),
     &           Cff(iCffW(0,4)),ddx(nRys),
     &           HerW2(iHerW2(4)),HerR2(iHerR2(4)),TMax_)
*
      Else If (nRys.eq.5) Then
      Call Rys55(Arg,nArg,Root,Weight,
     &           Map(iMap(5)),nMap(5),
     &           x0(ix0(5)),nx0(5),Cff(iCffR(6,5)),
     &           Cff(iCffR(5,5)),Cff(iCffR(4,5)),Cff(iCffR(3,5)),
     &           Cff(iCffR(2,5)),Cff(iCffR(1,5)),Cff(iCffR(0,5)),
     &           Cff(iCffW(6,5)),Cff(iCffW(5,5)),Cff(iCffW(4,5)),
     &           Cff(iCffW(3,5)),Cff(iCffW(2,5)),Cff(iCffW(1,5)),
     &           Cff(iCffW(0,5)),ddx(nRys),
     &           HerW2(iHerW2(5)),HerR2(iHerR2(5)),TMax_)
*
      Else If (nRys.eq.6) Then
      Call Rys66(Arg,nArg,Root,Weight,
     &           Map(iMap(6)),nMap(6),
     &           x0(ix0(6)),nx0(6),Cff(iCffR(6,6)),
     &           Cff(iCffR(5,6)),Cff(iCffR(4,6)),Cff(iCffR(3,6)),
     &           Cff(iCffR(2,6)),Cff(iCffR(1,6)),Cff(iCffR(0,6)),
     &           Cff(iCffW(6,6)),Cff(iCffW(5,6)),Cff(iCffW(4,6)),
     &           Cff(iCffW(3,6)),Cff(iCffW(2,6)),Cff(iCffW(1,6)),
     &           Cff(iCffW(0,6)),ddx(nRys),
     &           HerW2(iHerW2(6)),HerR2(iHerR2(6)),TMax_)
*
      Else If (nRys.eq.7) Then
      Call Rys77(Arg,nArg,Root,Weight,
     &           Map(iMap(7)),nMap(7),
     &           x0(ix0(7)),nx0(7),Cff(iCffR(6,7)),
     &           Cff(iCffR(5,7)),Cff(iCffR(4,7)),Cff(iCffR(3,7)),
     &           Cff(iCffR(2,7)),Cff(iCffR(1,7)),Cff(iCffR(0,7)),
     &           Cff(iCffW(6,7)),Cff(iCffW(5,7)),Cff(iCffW(4,7)),
     &           Cff(iCffW(3,7)),Cff(iCffW(2,7)),Cff(iCffW(1,7)),
     &           Cff(iCffW(0,7)),ddx(nRys),
     &           HerW2(iHerW2(7)),HerR2(iHerR2(7)),TMax_)
*
      Else If (nRys.eq.8) Then
      Call Rys88(Arg,nArg,Root,Weight,
     &           Map(iMap(8)),nMap(8),
     &           x0(ix0(8)),nx0(8),Cff(iCffR(6,8)),
     &           Cff(iCffR(5,8)),Cff(iCffR(4,8)),Cff(iCffR(3,8)),
     &           Cff(iCffR(2,8)),Cff(iCffR(1,8)),Cff(iCffR(0,8)),
     &           Cff(iCffW(6,8)),Cff(iCffW(5,8)),Cff(iCffW(4,8)),
     &           Cff(iCffW(3,8)),Cff(iCffW(2,8)),Cff(iCffW(1,8)),
     &           Cff(iCffW(0,8)),ddx(nRys),
     &           HerW2(iHerW2(8)),HerR2(iHerR2(8)),TMax_)
*
      Else If (nRys.eq.9) Then
      Call Rys99(Arg,nArg,Root,Weight,
     &           Map(iMap(9)),nMap(9),
     &           x0(ix0(9)),nx0(9),Cff(iCffR(6,9)),
     &           Cff(iCffR(5,9)),Cff(iCffR(4,9)),Cff(iCffR(3,9)),
     &           Cff(iCffR(2,9)),Cff(iCffR(1,9)),Cff(iCffR(0,9)),
     &           Cff(iCffW(6,9)),Cff(iCffW(5,9)),Cff(iCffW(4,9)),
     &           Cff(iCffW(3,9)),Cff(iCffW(2,9)),Cff(iCffW(1,9)),
     &           Cff(iCffW(0,9)),ddx(nRys),
     &           HerW2(iHerW2(9)),HerR2(iHerR2(9)),TMax_)
*
      Else
      Call WarningMessage(2,
     &            ' vRysRW: nRys in vRysRW is larger than MaxRys!')
      Call Abend()
*
      End If
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         If (labcd.ne.0)
     &      Call Recprt(' In vRysRW: Roots ',' ',Root  ,nRys,nArg)
         Call Recprt(' In vRysRW: Weight',' ',Weight,nRys,nArg)
      End If
#endif
*     Call qExit('vRysRW')
      Return
      End
