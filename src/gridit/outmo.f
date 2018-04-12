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


      subroutine outmo_nosupport(imo,ipower,cmo,clincomb,cout,nbas,nmo)
      implicit real*8 (a-h,o-z)
      dimension cmo(nbas,nmo),clincomb(nmo),cout(nbas)
      Include 'WrkSpc.fh'

      if(imo.ne.0)then
         call fmove(cmo(1,imo),cout,nbas)
         call power_nosupport(cout,nbas,ipower)
         return
      else
        call fzero(cout,nbas)
        Call GetMem('TmpMo','Allo','Real',ipTmpMo,nbas)
        do 100 i=1,nmo
        if(clincomb(i).ne.0d0)then
          call fmove(cmo(1,i),Work(ipTmpMo),nbas)
          call power_nosupport(Work(ipTmpMO),nbas,ipower)
          call daxpy(nbas,clincomb(i),work(ipTmpMo),1,cout,1)
        endif
100     continue
        Call GetMem('TmpMo','Free','Real',ipTmpMo,nbas)
      endif
      return
      end
      subroutine power_nosupport(c,n,ipower)
      implicit real*8 (a-h,o-z)
      dimension c(n)

      if(ipower.eq.1)then
        return
      elseif(ipower.eq.2)then
        do 100 i=1,n
100     c(i)=c(i)*c(i)
      elseif(ipower.eq.-2)then
        do 200 i=1,n
200     c(i)=sign(c(i)*c(i),c(i))
      endif
      return
      end
      subroutine MakePab_nosupport(cmo,occ,cout,nCMO,nMOs,nIrrep,nBas)
      implicit real*8 (a-h,o-z)
      dimension cmo(nCMO),occ(nMOs),cout(nMOs), nBas(*)
      Include 'WrkSpc.fh'
        Call GetMem('List','List','Real',id,id)
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
