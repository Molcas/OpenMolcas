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
      Subroutine MSP(B,Bd,Gamma,Delta,nDim)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 B(nDim,nDim), Bd(nDim),Gamma(nDim),Delta(nDim)
*
*                              T       T            ( T)
*                    |(1-phi)/d g phi/d d|        | (g )
*                    |     T       T    T         | ( T)
*     B = B + (g  d )|phi/d d    -Phi*d g/(d d)**2| (d )
*
*
      iRout=212
      iPrint=nPrint(iRout)
*
*
      gd= DDot_(nDim,Gamma,1,Delta,1)
      dd= DDot_(nDim,Delta,1,Delta,1)
      gg= DDot_(nDim,Gamma,1,Gamma,1)
      phi=(One-((gd**2)/(dd*gg)))
      e_msp=(gd/dd)**2*((Two/(One-Phi*Sqrt(Phi)))-One)
      If (iPrint.ge.99) Then
         Call RecPrt(' MSP: Hessian',' ',B,nDim,nDim)
         Call RecPrt(' MSP: Delta',' ',Delta,nDim,1)
         Call RecPrt(' MSP: Gamma',' ',Gamma,nDim,1)
         Write (6,*) 'MSP: Phi=',Phi
         Write (6,*) 'gd,dd,gg=', gd,dd,gg
         Write (6,*) 'MSP: a=',Sqrt(Phi)
         Write (6,*) 'MSP: E_msp=',E_msp
      End If
      Do i = 1, nDim
         Do j = 1, nDim
            B(i,j) = B(i,j)
     &             + ((One-phi)/gd)*Gamma(i)*Gamma(j)
     &             + phi*( (Gamma(i)*Delta(j)+Delta(i)*Gamma(j))/dd
     &                   - gd*Delta(i)*Delta(j)/dd**2 )
         End Do
      End Do
*
      If (iPrint.ge.99)
     &   Call RecPrt(' MSP: Updated Hessian',' ',B,nDim,nDim)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Bd)
      End
