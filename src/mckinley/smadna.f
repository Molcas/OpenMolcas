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
      SubRoutine SmAdNa(ArrIn,nb,ArrOut,nop,
     &                  lOper,IndGrd,
     &                  iuv,IfGrd,Index,iDCar,rf,IFG,tr)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
c#include "print.fh"
      Real*8 ArrIn (nb,*),
     &       ArrOut(nb,*)
      Integer  lOper,
     &          IndGrd(3,4,0:nIrrep-1),iuv(3),Index(3,4),nOp(3)
      Logical   IfGrd(3,4),IFG(4),tr(4)
*
*     Statement function for Cartesian index
*
*
c     iRout = 200
c     iPrint = nPrint(iRout)
c     Call qEnter('SymAdO')
*
*--------Accumulate contributions
*
      iComp=0
      Do 102 iIrrep=0,nIrrep-1
        If (iAnd(lOper,2**iIrrep).ne.0) Then
          iComp=iComp+1
          Do 103 iCn=1,3
*           If (Index(idCar,iCn).ne.0) Then
            If ( (Indgrd(idCar,iCn,iIrrep).ne.0) .and.
     &         ( (index(idcar,icn).gt.0).or.tr(icn)))
     &      Then
*              Accumulate contribution to the gradient
               i1=0
               i2=0
               If (iCn.eq.1) Then
                   ps = DBLE( iPrmt( nOp(1), iChBas(1+idCar) ) )
                   Fact = rf*DBLE(iuv(1))/DBLE(nIrrep)
                   If (.not.tr(iCn)) Then
                    i1=Index(idCar,iCn)
                   Else
                    If (index(idcar,2).gt.0) i1=Index(idCar,2)
                    If (index(idCar,3).gt.0) i2=Index(idCar,3)
                    Fact=-Fact
                   End If
               Else If (iCn.eq.2) Then
                   ps=DBLE(iChTbl(iIrrep,nOp(2)))
                   ps = ps*DBLE( iPrmt( nOp(2), iChBas(1+idCar) ) )
                   Fact = rf*ps *
     &                    DBLE(iuv(2))/DBLE(nIrrep)
                   If (.not.tr(iCn)) Then
                    i1=Index(idCar,iCn)
                   Else
                    If (index(idcar,1).gt.0) i1=Index(idCar,1)
                    If (index(idCar,3).gt.0) i2=Index(idCar,3)
                    Fact=-Fact
                   End If
               Else
                   ps=DBLE(iChTbl(iIrrep,nOp(3)))
                   ps = ps*DBLE( iPrmt( nOp(3), iChBas(1+idCar) ) )
                   Fact = rf*ps *
     &                    DBLE(iuv(3))/DBLE(nIrrep)
                   If (.not.tr(iCn)) Then
                    i1=Index(idCar,iCn)
                   Else
                    If (index(idcar,1).gt.0) i1=Index(idCar,1)
                    If (index(idCar,2).gt.0) i2=Index(idCar,2)
                    Fact=-Fact
                   End If
               End if
            If (i1.ne.0)
     &          Call DaXpY_(nb,Fact,
     &                 ArrIn(1,i1),1,ArrOut(1,iComp),1)
            If (i2.ne.0)
     &          Call DaXpY_(nb,Fact,
     &                     ArrIn(1,i2),1,ArrOut(1,iComp),1)
           End If
 103  Continue
      End If
 102  Continue
*
*     Call GetMem(' Exit SymAdO','LIST','REAL',iDum,iDum)
c     Call qExit('SymAdO')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical_array(IfGrd)
         Call Unused_logical_array(IFG)
      End If
      End
