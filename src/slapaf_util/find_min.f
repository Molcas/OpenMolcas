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
      Subroutine Find_Min(nOrder,XStart,A,XMin,RC,XLow,XHi,ENew)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 A(0:nOrder)
      Logical RC
*
      iRout=117
      iPrint=nPrint(iRout)
      If (iPrint.ge.99) Then
          Call RecPrt('Find_Min: A',' ',A,1,nOrder+1)
      End If
      Thr=1.0D-12
      XValue=XStart
      RC=.True.
      MaxIter=100
      Do i = 1, MaxIter
         X=XValue
         fnc  = Zero
         XX   = One
         Do j=0,nOrder
            fnc=fnc+A(j)*XX
            XX = XX*X
         End Do
         dfnc = Zero
         XX   = One
         Do j=1,nOrder
            tmp = DBLE(j)
            dfnc=dfnc+A(j)*tmp*XX
            XX = XX*X
         End Do
         ddfnc = Zero
         XX    = One
         Do j=2,nOrder
            tmp = DBLE(j*j-j)
            ddfnc=ddfnc+A(j)*tmp*XX
            XX = XX*X
         End Do
         XInc=dfnc/ddfnc
         XValue=XValue-XInc
         If (iPrint.eq.99) Then
            Write (6,*) 'Fnc,dFnc,ddFnc=',Fnc,dFnc,ddFnc
         End If
         If (Abs(XInc).lt.Thr) Then
            ENew=fnc
            XMin=XValue
            Return
         End If
         If (XValue.gt.XHi) XValue=XHi
         If (XValue.lt.XLow) XValue=XLow
      End Do
      If (iPrint.ge.6)
     &   Write (6,*) '-- Too many iterations in Find_Min'
      RC=.False.
*
      Return
      End
