************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Drvel1(Grad)
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
      Logical TF,tstfnc
      Real*8 Grad(*)
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                       nIrrep/nStab(mdc),iChTbl,iIrrep,iComp,
     &                       nStab(mdc))
      idisp=0
      do jIrrep=0,nirrep-1
       Do Jcar=1,3
        iirrep=irrfnc(2**(jcar-1))
        If (jirrep.eq.iirrep) Then
         mdc=0
         Do  iCnttp = 1, nCnttp
          ZA = Charge(iCnttp)
          Do  iCnt = 1, dbsc(iCnttp)%nCntr
           mdc=mdc+1
           Do iCar=1,3
            iComp = 2**(iCar-1)
            If ( TF(mdc,jIrrep,iComp)) Then
              idisp=idisp+1
              If (icar.eq.jcar) Grad(idisp)=ZA
            End If
           End Do
          End Do
         End Do
        End If
       End Do
      End Do
      Return
      End
