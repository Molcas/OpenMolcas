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
      Subroutine CiSelector(nEqState,nState,iSTC,nCIRef,iCIInd,dCIRef)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "WrkSpc.fh"

      Dimension iCIInd(MxState)
      Dimension dCIRef(MxState)

*
*-- Initial stuff
*
      dScalMAX=0.0d0
      indMAX=1
*
*-- Compute relevant scalar products
*
      Do 501, iState=1,nState
        dScal=0.0d0
        Do 502, iRef=1,nCIRef
          indBase=nState*(iState-1)
          indx=indBase+iCIInd(iRef)-1
          dScal=dScal+Work(iSTC+indx)*dCIRef(iRef)
502     Continue
        dScal=abs(dScal)
*
*---- Test if largest
*
        If(dScal.gt.dScalMAX) then
          dScalMAX=dScal
          indMAX=iState
        Endif
501   Continue

*
*-- If maximum overlap is small, scream!
*
      If(dScalMAX.lt.0.7071067811d0) then
        Write(6,*)
        Write(6,*)'   WARNING! Less than 50% of CISElect reference'
     &//'found. Consider to redefine reference!'
      Endif

*
*-- Now set nEqState
*
      nEqState=indMAX

*
*-- Auf Wiedersehen
*
      Return
      End
