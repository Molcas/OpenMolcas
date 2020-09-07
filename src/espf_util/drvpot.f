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
* Copyright (C) 1991, Roland Lindh                                     *
*               2001, Hans-Joachim Werner                              *
************************************************************************
      SubRoutine Drvpot(Ccoor,opnuc,ncmp,ptchrg,ngrid,iaddpot)
************************************************************************
*                                                                      *
* Object: driver for computation of one-electron property matrices     *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              OneEl                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '91                                              *
*                                                                      *
*     Modified for Properties only by HJW Aug 2001                     *
*     Restricted to POT: Ignacio Fdez. Galvan, March 2019              *
************************************************************************
      use Basis_Info
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
      External PotInt, NAMem
#include "stdalloc.fh"
#include "espf.fh"
      Character*8 Label
      Real*8 Ccoor(3),opnuc(*),ptchrg(*)
      Real*8, Allocatable :: Centr(:,:)
      Logical Do_ESPF
      Dimension dummy(1),iopadr(1)
*                                                                      *
************************************************************************
*                                                                      *
      Call qEnter('DrvPot')
*
      Call IniSewM('mltpl',0)
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      ntdg = 0
      Do iIrrep = 0, nIrrep - 1
         ntdg = ntdg + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
      Call DecideOnESPF(Do_ESPF)
c
      Call mma_allocate(Centr,3,mCentr)
      ndc = 0
      nc = 1
      Do jCnttp = 1, nCnttp
         mCnt = dbsc(jCnttp)%nCntr
         If (dbsc(jCnttp)%Aux) mCnt = 0
         Do jCnt = 1, mCnt
            ndc = ndc + 1
            Do i = 0, nIrrep/dc(ndc)%nStab - 1
               Call OA(iCoset(i,0,ndc),dbsc(jCnttp)%Coor(1:3,jCnt),
     &                 Centr(1:3,nc))
               nc = nc + 1
            End Do
            jxyz = jxyz + 3
         End Do
      End Do
      nc = nc-1
c
      nComp=1
      nOrdOp = 0
      Call GetMem('ip    ','ALLO','INTE',ip1,nComp)
      Call GetMem('lOper ','ALLO','INTE',ip2,nComp)
      Call GetMem('kOper ','ALLO','INTE',ip3,nComp)
      Label='Pot '
      If (iaddpot.le.0.and..not.Do_ESPF) Then
         Call GetMem('Nuc ','ALLO','REAL',ipNuc,ngrid)
         Call Pot_nuc(CCoor,work(ipnuc),ngrid)
      Else
        ipnuc=ip_Dummy
      End if
      If (iaddpot.lt.0) Then
         If (iaddpot.eq.-1) Then
            Call Get_D1ao_Var(ipdens,Length)
         Else
            Call Get_D1ao(ipdens,Length)
         End If
         call Drv1_Pot(work(ipdens),CCoor,ptchrg,ngrid,1,0)
         Call GetMem('DENS','FREE','REAL',ipdens,ntdg)
         If (.not.Do_ESPF) Then
            Call AddVec(ptchrg,ptchrg,work(ipnuc),ngrid)
            Call dCopy_(ngrid,work(ipnuc),1,opnuc,1)
         End If
      Else
        iWork(ip2) = 2**nirrep-1
        Call OneEl(PotInt,NAMem,Label,iWork(ip1),iWork(ip2),ncmp,
     &             Ccoor,nOrdOp,work(ipnuc),rHrmt,iWork(ip3),
     &             dummy,1,opnuc,iopadr,1,1,
     &             ptchrg,ngrid,iaddpot)
         If (iaddpot.eq.0.and..not.Do_ESPF)
     &      opnuc(1)=work(ipnuc)
      End If
      If (iaddpot.le.0.and..not.Do_ESPF)
     &   Call GetMem('Nuc ','FREE','REAL',ipNuc,ngrid)
      Call GetMem('kOper ','FREE','INTE',ip3,nComp)
      Call GetMem('lOper ','FREE','INTE',ip2,nComp)
      Call GetMem('ip    ','FREE','INTE',ip1,nComp)
*
      Call mma_deallocate(Centr)
      Call Free_iSD()
*
      Call QExit('DrvPot')
*
      Return
      End
