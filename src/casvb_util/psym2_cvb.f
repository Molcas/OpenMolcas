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
      subroutine psym2_cvb(civec1,civec2,
     >  isymalf,isymbet,iasyind,ibsyind,osym,ips)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension civec1(nda,ndb),civec2(nda,ndb)
      dimension isymalf(nda),isymbet(ndb)
      dimension iasyind(0:mxirrep),ibsyind(0:mxirrep)
      dimension osym(mxirrep)

      if(ips.eq.1)then
        do 1000 irp=1,nirrep
        if(isympr(irp).eq.1)goto 1000
        do 100 jrpa=1,nirrep
        jrpb=md2h(irp,jrpa)
        do 101 ida=iasyind(jrpa-1)+1,iasyind(jrpa)
        inda=isymalf(ida)
        do 102 idb=ibsyind(jrpb-1)+1,ibsyind(jrpb)
        civec1(inda,isymbet(idb))=zero
102     continue
101     continue
100     continue
1000    continue
      elseif(ips.eq.2)then
        do 2000 irp=1,nirrep
        osym(irp)=zero
        do 200 jrpa=1,nirrep
        jrpb=md2h(irp,jrpa)
        do 201 ida=iasyind(jrpa-1)+1,iasyind(jrpa)
        inda=isymalf(ida)
        do 202 idb=ibsyind(jrpb-1)+1,ibsyind(jrpb)
        osym(irp)=osym(irp)+
     >    civec1(inda,isymbet(idb))*civec2(inda,isymbet(idb))
202     continue
201     continue
200     continue
2000    continue
      endif
      return
      end
