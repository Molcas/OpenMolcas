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
      subroutine sminus2_cvb(bikfrom,bikto,
     >  nel,nalffrom,ndetfrom,nalfto,ndetto,nvec,
     >  xdetto,ioccfrom,ioccto)
      implicit real*8 (a-h,o-w,y-z),integer(x)
      dimension bikfrom(ndetfrom,nvec),bikto(ndetto,nvec)
      dimension xdetto(0:nel,0:nalfto)
      dimension ioccfrom(nalffrom),ioccto(nalfto)
#include "malloc_cvb.fh"

      call fzero(bikto,ndetto*nvec)

c Determinant (to) weight array:
      call weightfl_cvb(xdetto,nalfto,nel)
      if(ndetto.ne.xdetto(nel,nalfto))then
        write(6,*) ' Discrepancy in NDET:',ndetto,xdetto(nel,nalfto)
        call abend_cvb()
      endif

      call loopstr0_cvb(ioccfrom,indfrom,nalffrom,nel)
100   continue
      call imove_cvb(ioccfrom(2),ioccto,nalfto)
      do 200 iexc=1,nalffrom
      indto=minind_cvb(ioccto,nalfto,nel,xdetto)
      call daxpy_(nvec,1d0,bikfrom(indfrom,1),ndetfrom,
     >  bikto(indto,1),ndetto)
      if(iexc.lt.nalffrom)ioccto(iexc)=ioccfrom(iexc)
200   continue
      call loopstr_cvb(ioccfrom,indfrom,nalffrom,nel)
      if(indfrom.ne.1)goto 100
      return
      end
