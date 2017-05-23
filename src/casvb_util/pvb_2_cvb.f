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
      subroutine pvb_2_cvb(cfrom,cto,csk,
     >  iapr,ixapr,mult)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension cfrom(nda,ndb),cto(nda,ndb),csk(ndetvb)
      dimension iapr(ndetvb),ixapr(nda+1)

      if(mult.eq.-1)then
        idetvb=0
        do 50 ia=1,nda
        do 50 ixa=ixapr(ia),ixapr(ia+1)-1
        idetvb=idetvb+1
        ib=iapr(ixa)
50      csk(idetvb)=cfrom(ia,ib)
      elseif(mult.eq.0)then
        call fzero(cto,nda*ndb)
        idetvb=0
        do 100 ia=1,nda
        do 100 ixa=ixapr(ia),ixapr(ia+1)-1
        idetvb=idetvb+1
        ib=iapr(ixa)
        cto(ia,ib)=cfrom(ia,ib)
100     csk(idetvb)=cfrom(ia,ib)
      elseif(mult.eq.1)then
        csk(1)=zero
        do 200 ia=1,nda
        do 200 ixa=ixapr(ia),ixapr(ia+1)-1
200     csk(1)=csk(1)+cto(ia,iapr(ixa))*cfrom(ia,iapr(ixa))
      elseif(mult.eq.2)then
        call fzero(cto,nda*ndb)
        idetvb=0
        do 300 ia=1,nda
        do 300 ixa=ixapr(ia),ixapr(ia+1)-1
        idetvb=idetvb+1
        ib=iapr(ixa)
300     cto(ia,ib)=csk(idetvb)
      elseif(mult.eq.3)then
        csk(1)=zero
        idetvb=0
        do 400 ia=1,nda
        do 400 ixa=ixapr(ia),ixapr(ia+1)-1
        idetvb=idetvb+1
c CFROM is really CDETVB
400     csk(1)=csk(1)+cto(ia,iapr(ixa))*cfrom(idetvb,1)
      endif
      return
      end
