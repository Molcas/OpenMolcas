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
      subroutine mol2vb2_cvb(vecvb,vecmol,isyml,fac,iwr,
     >  indxa,indxb,nstra,nstrb,nsa,nsb)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension vecvb(ndet),vecmol(*)
      dimension indxa(nsa),indxb(nsb)
      dimension nstra(mxirrep),nstrb(mxirrep)

      call indxab_cvb(indxa,indxb,nstra,nstrb,nsa,nsb)

c  Now loop casvb -> molcas
      idet=0
      do 500 isyma=1,mxirrep
      isymb=md2h(isyma,isyml)
      nnsa=nstra(isyma)
      nnsb=nstrb(isymb)
      if(nnsa.le.0 .or. nnsb.le.0)goto 500

      ioffsa=0
      do 510 is=1,isyma-1
      ioffsa=ioffsa+nstra(is)
510   continue
      ioffsb=0
      do 520 is=1,isymb-1
      ioffsb=ioffsb+nstrb(is)
520   continue

      do 600 isb=1,nnsb
      indbet=indxb(isb+ioffsb)
      do 601 isa=1,nnsa
      index=indxa(isa+ioffsa)+(indbet-1)*nda
      idet=idet+1
      if(iwr.eq.0)then
        vecmol(idet)=vecvb(index)
      elseif(iwr.eq.1)then
        vecvb(index)=vecmol(idet)
      elseif(iwr.eq.2)then
        vecvb(index)=vecvb(index)+fac*vecmol(idet)
      endif
601   continue
600   continue
500   continue
      return
      end
