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
* Copyright (C) 2002, Roland Lindh                                     *
************************************************************************
      SubRoutine DrvDFTg(Grad,Temp,nGrad)
************************************************************************
*                                                                      *
* Object: driver for computation of gradient with respect to the DFT   *
*         energy.                                                      *
*                                                                      *
* Called from: Alaska or Drvg1                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              OneEl                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chem. Phys.                       *
*             University of Lund, SWEDEN                               *
*             August 2002                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "rctfld.fh"
#include "disp.fh"
#include "nq_info.fh"
      Character Label*80, KSDFT*16
      Real*8 Grad(nGrad), Temp(nGrad)
      Logical First, Dff, Do_Grad, King, l_casdft
      Character*4 DFTFOCK
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu1,TWall1)
*                                                                      *
************************************************************************
*                                                                      *
*...  Prologue
      DFTFOCK='SCF '
      iRout = 131
      iPrint = nPrint(iRout)
      Call qEnter('DrvDFTg')
      LuWr=6
*
      nDens = 0
      Do iIrrep = 0, nIrrep - 1
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
*
*     D F T - g r a d i e n t                                          *
************************************************************************
*8)                                                                    *
*     D F T - g r a d i e n t                                          *
*                                                                      *
*     Call Get_iOption(iDFT)

      Call Get_cArray('DFT functional',KSDFT,16)
      l_casdft = KSDFT(1:5).eq.'TLSDA'   .or.
     &           KSDFT(1:6).eq.'TLSDA5'  .or.
     &           KSDFT(1:5).eq.'TBLYP'   .or.
     &           KSDFT(1:6).eq.'TSSBSW'  .or.
     &           KSDFT(1:5).eq.'TSSBD'  .or.
     &           KSDFT(1:5).eq.'TS12G'  .or.
     &           KSDFT(1:4).eq.'TPBE'    .or.
     &           KSDFT(1:5).eq.'FTPBE'   .or.
     &           KSDFT(1:7).eq.'TREVPBE' .or.
     &           KSDFT(1:8).eq.'FTREVPBE'.or.
     &           KSDFT(1:6).eq.'FTLSDA'  .or.
     &           KSDFT(1:6).eq.'FTBLYP'

      Call Get_iScalar('System BitSwitch',iDFT)
      If( l_casdft ) then
        DFTFOCK='ROKS'
        iDFT=iOr(iDFT,2**6)
        Call Put_iScalar('System BitSwitch',iDFT)
      End IF

      Call Get_iScalar('System BitSwitch',iDFT)
      If (iAnd(iDFT,2**6).ne.0) Then
*
         Call StatusLine(' Alaska:',' Computing DFT gradients')
*
         First=.True.
         Dff  =.False.
         Call Get_cArray('DFT functional',KSDFT,16)
         ExFac=Zero ! Set to proper value at retrun!
         Do_Grad=.True.
         Call Get_iScalar('Multiplicity',iSpin)
*        Write (LuWr,*) 'DrvDFTg: KSDFT=',KSDFT
*        Write (LuWr,*) 'DrvDFTg: ExFac=',ExFac
         iDumm=1
         Call DrvDFT(Dummy1,Dummy2,Dummy3,Dummy4,nDens,First,Dff,
     &               lRF,KSDFT,ExFac,Do_Grad,Temp,nGrad,iSpin,
     &               Dumm0,Dumm1,iDumm,DFTFOCK)
*
         iEnd=1
 99      Continue
         If (KSDFT(iEnd:iEnd).eq.' ') Then
            iEnd = iEnd - 1
         Else
            iEnd = iEnd + 1
            Go To 99
         End If
         Label='The DFT('//KSDFT(1:iEnd)//') contribution'
         jPrint=nPrint(112)
!AMS
!        jprint=15
         If (jPrint.ge.15) Call PrGrad(Label,Temp,nGrad,lIrrep,ChDisp,5)
         If (king()) Call DaXpY_(nGrad,One,Temp,1,Grad,1)
         If (iPrint.lt.6) Go To 777
         Write (LuWr,*)
         If (Grid_Type.eq.Moving_Grid) Then
            Write(LuWr,*) 'DFT contribution computed for a moving grid.'
         Else
            Write(LuWr,*) 'DFT contribution computed for a fixed grid.'
         End If
         Write (LuWr,*)
 777     Continue

*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu2,TWall2)
      Call SavTim(5,TCpu2-TCpu1,TWall2-TWall1)
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('DrvDFTg')
      Return
      End
