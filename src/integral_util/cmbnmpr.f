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
* Copyright (C) Kurt Pfingst                                           *
************************************************************************
      SubRoutine CmbnMPr(Rnr,nZeta,la,lb,lr,Zeta,Final,nComp)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: MltInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: K.Pfingst                                                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      External gammat,gammaf
#include "print.fh"
#include "real.fh"
#include "rmat.fh"
#include "gam.fh"
#include "nrmf.fh"
      Real*8 Final(nZeta,nComp,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2),
     *       Zeta(nZeta), Rnr(nZeta,0:(la+lb+lr))
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
      iRout = 134
      iPrint = nPrint(iRout)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     iPrint = 99
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*     Call QEnter('CmbnMP')
*     Call GetMem(' Enter CmbnMP','LIST','REAL',iDum,iDum)
*
      Do 10 ixa = 0, la
         iyaMax=la-ixa
      Do 11 ixb = 0, lb
         iybMax=lb-ixb
         Do 20 iya = 0, iyaMax
            iza = la-ixa-iya
            ipa= Ind(la,ixa,iza)
         Do 21 iyb = 0, iybMax
            izb = lb-ixb-iyb
            ipb= Ind(lb,ixb,izb)
            If (iPrint.ge.99) Then
               Write (6,*) ixa,iya,iza,ixb,iyb,izb
               Write (6,*) ipa,ipb
            End If
*
*           Combine multipole moment integrals
*
            iComp = 0
            Do 41 ix = lr, 0, -1
               Do 42 iy = lr-ix, 0, -1
                  iz = lr-ix-iy
                  iComp=iComp+1
*                 Write (*,*) ix, iy, iz, iComp
                  lrs=ixa+ixb+ix+iya+iyb+iy+iza+izb+iz
                  lcost=iza+izb+iz
                  lsint=ixa+ixb+ix+iya+iyb+iy
                  lsinf=iya+iyb+iy
                  lcosf=ixa+ixb+ix
                  Fact=gammath(lsint,lcost)*gammaph(lsinf,lcosf)
*                 Fact1=gammat(x)*gammaf(x)
*                 write(*,*) '  fact',fact
*                 write(*,*) ' fact1',fact1
                  Do 30 iZeta = 1, nZeta
                     Final(iZeta,iComp,ipa,ipb) = Fact *
     *                       Rnr(iZeta,lrs)
30                Continue
42             Continue
41          Continue
*
21       Continue
20       Continue
11    Continue
10    Continue
*
*     Call GetMem(' Exit CmbnMP','LIST','REAL',iDum,iDum)
*     Call QExit('CmbnMP')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Zeta)
      End
