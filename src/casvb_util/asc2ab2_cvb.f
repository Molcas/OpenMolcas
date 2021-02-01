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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine asc2ab2_cvb(detvec,nvec,nel,nalf,
     >  nbet,ndet,
     >  mindet,maxdet,nkdet,xdet,locc)
      implicit real*8 (a-h,o-w,y-z),integer(x)
      dimension detvec(ndet,nvec)
      dimension mindet(0:nel),maxdet(0:nel),nkdet(0:nel),
     >  xdet(0:nel,0:nalf)
      dimension locc(nel)

      do 100 iorb=0,nel
      mindet(iorb)=max(iorb-nbet,0)
      maxdet(iorb)=min(iorb,nalf)
100   continue
      call weight_cvb(xdet,mindet,maxdet,nalf,nel)
      call imove_cvb(maxdet,nkdet,nel+1)
      call occupy_cvb(nkdet,nel,locc,locc(nalf+1))
      inddet=1
200   continue
      call dscal_(nvec,party_cvb(locc,nel),detvec(inddet,1),ndet)
      call loind_cvb(nel,nalf,nkdet,mindet,maxdet,
     >               locc,locc(nalf+1),inddet,xdet,*200)
      return
      end
c  ********************************
c  ** VB determinant information **
c  ********************************
