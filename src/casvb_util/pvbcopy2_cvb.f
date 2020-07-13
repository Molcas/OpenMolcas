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
      subroutine pvbcopy2_cvb(cfrom,cto,
     >  iapr,ixapr,ret,ic)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


      dimension cfrom(nda,ndb),cto(nda,ndb)
      dimension iapr(ndetvb),ixapr(nda+1)

      if(ic.eq.0)then
        call fzero(cto,nda*ndb)
        idetvb=0
        do 100 ia=1,nda
        do 101 ixa=ixapr(ia),ixapr(ia+1)-1
        idetvb=idetvb+1
        ib=iapr(ixa)
        cto(ia,ib)=cfrom(ia,ib)
101     continue
100     continue
      elseif(ic.eq.1)then
        ret=zero
        do 200 ia=1,nda
        do 201 ixa=ixapr(ia),ixapr(ia+1)-1
        ret=ret+cto(ia,iapr(ixa))*cfrom(ia,iapr(ixa))
201     continue
200     continue
      endif
      return
      end
