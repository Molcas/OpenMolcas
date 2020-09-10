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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
c  *********************************************************************
c  *                                                                   *
c  *  CNFINI   := Set NVBR, NDETVB, NDETVB2, MNION, MXION,             *
c  *              NCONFION, and IFSC.                                  *
c  *                                                                   *
c  *********************************************************************
      subroutine cnfini_cvb(iconfs,nconf1,nel1,
     >  nS,i2s,nMs,nalf1,nbet1,
     >  nvbr1,ndetvb1,ndetvb21,mnion1,mxion1,nconfion,ifsc1)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


      dimension iconfs(noe,nconf1)
      dimension i2s(nS),nalf1(nMs),nbet1(nMs)
      dimension nconfion(0:*)

c  Main loop over configurations :
      mnion1=nel1/2
      mxion1=0
      call izero(nconfion,1+nel1/2)
      ndetvb1=0
      ndetvb21=0
      nvbr1=0
      do 100 iconf=1,nconf1
      ion=0
      do 200 iorb=1,norb
      if(iconfs(iorb,iconf).eq.2)ion=ion+1
200   continue
      if(ion.lt.mnion1)mnion1=ion
      if(ion.gt.mxion1)mxion1=ion
      nconfion(ion)=nconfion(ion)+1
      do 300 iS=1,nS
      call icomb_cvb(nel1-2*ion,(nel1-i2s(iS))/2-ion,iretval1)
      call icomb_cvb(nel1-2*ion,(nel1-i2s(iS))/2-ion-1,iretval2)
      nvbr1=nvbr1+iretval1-iretval2
300   continue
      do 400 iMs=1,nMs
      call icomb_cvb(nel1-2*ion,nalf1(iMs)-ion,iretval)
      ndetvb1=ndetvb1+iretval
      ndetvb21=ndetvb21+(iretval+1)/2
400   continue
100   continue
      if(norb.eq.nel1.and.nconf1.eq.1)then
        ifsc1=1
        do 500 i=1,nel1
        if(iconfs(i,nconf1).ne.1)ifsc1=0
500     continue
      else
        ifsc1=0
      endif
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(nbet1)
      end
