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
c  *  PVB    := Zero parts of CI vector not in VB wfn.                 *
c  *                                                                   *
c  *********************************************************************
      subroutine ci2vb2_cvb(civec,cvbdet,
     >  iapr,ixapr,ret,ic)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension civec(nda,ndb),cvbdet(ndetvb)
      dimension iapr(ndetvb),ixapr(nda+1)

      if(ic.eq.0)then
        idetvb=0
        do 100 ia=1,nda
        do 100 ixa=ixapr(ia),ixapr(ia+1)-1
        idetvb=idetvb+1
        ib=iapr(ixa)
100     cvbdet(idetvb)=civec(ia,ib)
      elseif(ic.eq.1)then
        call fzero(civec,nda*ndb)
        idetvb=0
        do 200 ia=1,nda
        do 200 ixa=ixapr(ia),ixapr(ia+1)-1
        idetvb=idetvb+1
        ib=iapr(ixa)
200     civec(ia,ib)=cvbdet(idetvb)
      elseif(ic.eq.2)then
        idetvb=0
        do 300 ia=1,nda
        do 300 ixa=ixapr(ia),ixapr(ia+1)-1
        idetvb=idetvb+1
        ib=iapr(ixa)
300     civec(ia,ib)=civec(ia,ib)+cvbdet(idetvb)
      elseif(ic.eq.3)then
        ret=zero
        idetvb=0
        do 400 ia=1,nda
        do 400 ixa=ixapr(ia),ixapr(ia+1)-1
        idetvb=idetvb+1
        ib=iapr(ixa)
400     ret=ret+civec(ia,ib)*cvbdet(idetvb)
      endif
      return
      end
