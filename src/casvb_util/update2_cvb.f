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
      subroutine update2_cvb(orbs,cvb,orbsp,cvbp,
     >  sorbs,dxorg,
     >  ic,
     >  norb,nvb,nprorb,npr,orbopt,strucopt,sym,
     >  dx,iorts,nort,sorbsinv)
      implicit real*8 (a-h,o-z)
      logical orbopt,strucopt,sym
#include "print_cvb.fh"
      dimension orbs(norb,norb),cvb(nvb)
      dimension orbsp(norb,norb),cvbp(nvb)
      dimension sorbs(norb,norb)
      dimension dxorg(npr),dx(npr)
      dimension iorts(2,nort)
      dimension sorbsinv(norb,norb)
      save zero
      data zero/0d0/

      call free2all_cvb(dxorg,dx,1)
      if(ip(3).ge.3.and.ic.eq.1)then
        write(6,'(/,a)')' Update vector :'
        call vecprint_cvb(dx,npr)
      endif

      call fmove(orbsp,orbs,norb*norb)
      call fmove(cvbp,cvb,nvb)

      if(orbopt)then
        call mxattb_cvb(orbsp,orbsp,norb,norb,norb,sorbs)

        ij=0
        do 200 iorb=1,norb
        do 200 jorb=1,norb
        if(jorb.ne.iorb)then
          ij=ij+1
          do 400 i=1,norb
400       orbs(i,iorb)=orbs(i,iorb)+dx(ij)*orbsp(i,jorb)
        endif
200     continue

c  2nd-order correction for orthogonality constraints:
        call fmove(sorbs,sorbsinv,norb*norb)
        call mxinv_cvb(sorbsinv,norb)
        do 500 iort=1,nort
        iorb=iorts(1,iort)
        jorb=iorts(2,iort)
        sdidj=zero
        do 600 k=1,norb-1
        korb=k
        if(korb.ge.iorb)korb=korb+1
        do 600 l=1,norb-1
        lorb=l
        if(lorb.ge.jorb)lorb=lorb+1
600     sdidj=sdidj+sorbs(korb,lorb)*dx(k+(iorb-1)*(norb-1))
     >    *dx(l+(jorb-1)*(norb-1))
        fac=-.5d0*sdidj
        do 700 i=1,norb
        do 700 j=1,norb
        orbs(i,iorb)=orbs(i,iorb)+fac*orbsp(i,j)*sorbsinv(j,jorb)
700     orbs(i,jorb)=orbs(i,jorb)+fac*orbsp(i,j)*sorbsinv(j,iorb)
500     continue
      endif

      if(strucopt)then
        call addvec(cvb,cvb,dx(nprorb+1),nvb)
        call scalstruc_cvb(orbs,cvb)
        call cvbnrm_cvb(cvb)
      endif

      call nize_cvb(orbs,norb,dum,norb,0,0)
      if(sym)call symtriz_cvb(orbs,cvb)
      return
      end
