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
      SubRoutine SymAdO_mck2(ArrIn,nB,ArrOut,nrOp,nop,
     &                  IndGrd,ksym,iu,iv,ifgrd,idCar,trans)
      use Symmetry_Info, only: iChTbl, iOper, iChBas
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
      Real*8 ArrIn (nB,2), ArrOut(nB,nrOp)
      Integer IndGrd(0:nIrrep-1),nop(2)
      Logical IfGrd(3,2),trans(2)
*
*--------Accumulate contributions
*
      n=0
      Do  jIrrep=0,nIrrep-1
       iirrep=nropr(ieor(ioper(jirrep),ioper(ksym)))
       If (Indgrd(jIrrep).ne.0) Then
        n=n+1
        Do iCn=1,2
          If ((Trans(iCn).or.IfGrd(idCar,iCn)).and.
     &          (IndGrd(jIrrep).ne.0)) Then
*
*              Accumulate contribution to the gradient
*
            If (iCn.eq.1) Then
                   ps = DBLE( iPrmt( nOp(1),iChBas(1+idCar) ) )
                   Fact = DBLE(iu)/DBLE(nIrrep)
                   If (trans(1)) Then
                     Fact=-Fact
                   End If
            Else
                   ps=DBLE(iChTbl(iIrrep,nOp(2)))
                   ps = ps*DBLE( iPrmt( nOp(2),iChBas(1+idCar) ) )
                   Fact = ps * DBLE(iv)/DBLE(nIrrep)
                   If (trans(2)) Then
                    Fact=-Fact
                   End If
            End if
            Call DaXpY_(nB,Fact,ArrIn(1,icn),1,ArrOut(1,n),1)
          End If
        End Do
       End If
      End Do
*
      Return
      End
