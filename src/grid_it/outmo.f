************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine outmo(imo,ipower,cmo,clincomb,cout,nbas,nmo)
************************************************************************
* Adapted from SAGIT to work with OpenMolcas (October 2020)            *
************************************************************************
      implicit real*8 (a-h,o-z)
      dimension cmo(nbas,nmo),clincomb(nmo),cout(nbas)
#include "WrkSpc.fh"

      if(imo.ne.0)then
         call fmove(cmo(1,imo),cout,nbas)
         call power(cout,nbas,ipower)
         return
      else
        call fzero(cout,nbas)
        Call GetMem('TmpMo','ALLO','REAL',ipTmpMo,nbas)
        do 100 i=1,nmo
        if(clincomb(i).ne.0d0)then
          call fmove(cmo(1,i),Work(ipTmpMo),nbas)
          call power(Work(ipTmpMO),nbas,ipower)
          call daxpy_(nbas,clincomb(i),Work(ipTmpMo),1,cout,1)
        endif
100     continue
        Call GetMem('TmpMo','FREE','REAL',ipTmpMo,nbas)
      endif
      return
      end

      subroutine power(c,n,ipower)
      implicit real*8 (a-h,o-z)
      dimension c(n)

      if(ipower.eq.1)then
        return
      elseif(ipower.eq.2)then
        do 100 i=1,n
         c(i)=c(i)*c(i)
100     continue
      elseif(ipower.eq.-2)then
        do 200 i=1,n
          c(i)=sign(c(i)*c(i),c(i))
200     continue
      endif
      return
      end
      subroutine MakePab(cmo,occ,cout,nCMO,nMOs,nIrrep,nBas)
      implicit real*8 (a-h,o-z)
      dimension cmo(nCMO),occ(nMOs),cout(nMOs), nBas(0:7)
#include "WrkSpc.fh"
        Call GetMem('List','LIST','REAL',id,id)
        id=0
        id2=0
        call fzero(cout,nMOs)
        do iIrr=0,nIrrep-1
        nd=nBas(iIrr)
        nd2=nd*nd
        do 100 i=1,nd
c        if(occ(i+id).ne.0d0)then
        do 200 j=1,nd
        cout(i+id)=cout(i+id)+occ(i+id)*
     &                ( cmo(j+i*(nd-1)+id2) ** 2 )
200     continue
100     continue
        id=id+nd
        id2=id2+nd2
        enddo

      return
      end
