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
* Copyright (C) 1996, Anders Bernhardsson                              *
*               2002, Roland Lindh                                     *
************************************************************************
      SubRoutine GrdTr_Alaska(GradIn,MxAto,GradOut,nGrad)
********************************************************************
*                                                                  *
*      The inverse of                                              *
*      Transforms a symmetry adapted gradient to unsymmetric  form *
*                                                                  *
*       Written by Anders Bernhardsson                             *
*       960427                                                     *
*       Modified by Roland Lindh                                   *
*       020115                                                     *
*                                                                  *
********************************************************************
      Implicit Real*8(a-h,o-z)
      parameter (tol=1d-8)
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 GradIn(3,MxAto), GradOut(nGrad)
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
         If (AuxCnttp(iCnttp)) Return
         Do iCnt=1,nCntr(iCnttp)
            mdc=mdc+1
            nDispS = IndDsp(mdc,iIrrep)

            Do iCar=0,2
               iComp = 2**iCar
               If ( TF(mdc,iIrrep,iComp)) Then
                  nDispS = nDispS + 1
                  GradOut(nDispS)=GradIn(iCar+1,iCen)
               End If
            End Do
            iCen = iCen + nIrrep/nStab(mdc)
         End Do
      End Do
*
      Return
      End
