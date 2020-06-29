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
      SubRoutine Cmbn_EF_DPnt(EF,nTs,DPnt,MxAto,DCntr,nS,iSph,Q,
     &                        Grad,nGrad)
********************************************************************
*                                                                  *
*      Combine EF with DPnt array.                                 *
*                                                                  *
*      Roland Lindh                                                *
*      020117                                                      *
*                                                                  *
********************************************************************
      use Basis_Info
      Implicit Real*8(a-h,o-z)
      parameter (tol=1d-8)
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 EF(3,2,nTs),DPnt(nTs,MxAto,3,3), Grad(nGrad),
     &       DCntr(nS,MxAto,3,3), Q(2,nTs)
      Integer iSph(nTs)
      Logical TF,TstFnc
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                       nIrrep/nStab(mdc),iChTbl,iIrrep,iComp,
     &                       nStab(mdc))
*
      iIrrep=0
*
      mdc=0
      iCen=1
      Do iCnttp=1,nCnttp
         If (AuxCnttp(iCnttp)) Cycle
         Do iCnt=1,dbsc(iCnttp)%nCntr
            mdc=mdc+1
            nDispS = IndDsp(mdc,iIrrep)

            Do iCar=0,2
               iComp = 2**iCar
               If ( TF(mdc,iIrrep,iComp)) Then
                  nDispS = nDispS + 1
                  Do iTs = 1, nTs
                     jSph=iSph(iTs)
                     QTot=Q(1,iTs)+Q(2,iTs)
                     Grad(nDispS) = Grad(nDispS)
     &                  + QTot * (
     &                    (EF(1,1,iTs)-EF(1,2,iTs))
     &                  *(DPnt(iTs,iCen,iCar+1,1)
     &                   +DCntr(jSph,iCen,iCar+1,1))
     &                  + (EF(2,1,iTs)-EF(2,2,iTs))
     &                  *(DPnt(iTs,iCen,iCar+1,2)
     &                   +DCntr(jSph,iCen,iCar+1,2))
     &                  + (EF(3,1,iTs)-EF(3,2,iTs))
     &                  *(DPnt(iTs,iCen,iCar+1,3)
     &                   +DCntr(jSph,iCen,iCar+1,3))
     &                           )
                  End Do
               End If
            End Do
            iCen = iCen + nIrrep/nStab(mdc)
         End Do
      End Do
*
      Return
      End
