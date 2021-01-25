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
************************************************************************
      SubRoutine HRR1(ab1,nab1,a1b,na1b,cffAB,ab,nab,
     &           na,nb,na1,nb1,nPrim,la,lb)
************************************************************************
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             June '91                                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
*     Character*72 Label
      Real*8 ab1(nPrim,nab1), a1b(nPrim,na1b), cffAB(3),
     &       ab(nPrim,nab)
*
      Ind(iy,iz) = (iy+iz)*(iy+iz+1)/2 + iz + 1
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
*     iRout = 99
*     iPrint = nPrint(iRout)
*     If (iPrint.ge.99) Then
*        Write(Label,'(A,i1,A,i1,A)') ' Source: (',na1,',',nb,'|'
*        Call RecPrt(Label,' ',a1b,nPrim,na1b)
*        Call RecPrt(' Coordinates (A-B)',' ',cffAB,1,3)
*        Write(Label,'(A,i1,A,i1,A)') ' Source: (',na,',',nb,'|'
*        Call RecPrt(Label,' ',ab,nPrim,nab)
*     End If
*
*     Loop over indices of the target batch
*
      Do 20 ixb = nb1, 0, -1
         Do 25 iyb = nb1-ixb, 0, -1
            izb = nb1-ixb-iyb
            ixyzb1=Ind(iyb,izb)
*
            Do 30 ixa = na, 0, -1
               Do 35 iya = na-ixa, 0, -1
                  iza = na-ixa-iya
*
                  ixyza =Ind(iya,iza)
*
*                 Find a angular index which can be decremented
*
                  If (ixb.ne.0) Then
                     ipxyz = 1
                     ixyza1=Ind(iya,iza)
                     ixyzb =Ind(iyb,izb)
                  Else If (iyb.ne.0) Then
                     ipxyz = 2
                     ixyza1=Ind(iya+1,iza)
                     ixyzb =Ind(iyb-1,izb)
                  Else
                     ipxyz = 3
                     ixyza1=Ind(iya,iza+1)
                     ixyzb =Ind(iyb,izb-1)
                  End If
                  If (la.ge.lb) Then
                     ipab1 = ixyza  + nElem(na )*(ixyzb1-1)
                     ipa1b = ixyza1 + nElem(na1)*(ixyzb -1)
                     ipab  = ixyza  + nElem(na )*(ixyzb -1)
                  Else
                     ipab1 = ixyzb1 + nElem(nb1)*(ixyza -1)
                     ipa1b = ixyzb  + nElem(nb )*(ixyza1-1)
                     ipab  = ixyzb  + nElem(nb )*(ixyza -1)
                  End If
                  If (cffAB(ipxyz).ne.Zero) Then
                     Call DZaXpY(nPrim,cffAB(ipxyz),ab(1,ipab),1,
     &                                             a1b(1,ipa1b),1,
     &                                             ab1(1,ipab1),1)
                  Else
*                    call dcopy_(nPrim,a1b(1,ipa1b),1,ab1(1,ipab1),1)
                     Do 40 i = 1, nPrim
                        ab1(i,ipab1) = a1b(i,ipa1b)
 40                  Continue
                  End If
 35            Continue
 30         Continue
 25      Continue
 20   Continue
*
*     If (iPrint.ge.99) Then
*        Write(Label,'(A,i1,A,i1,A)') ' Target: (',na,',',nb1,'|'
*        Call RecPrt(Label,' ',ab1,nPrim,nab1)
*     End If
*
      Return
      End
