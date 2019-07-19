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
      subroutine getmo2_cvb(cmo,cmo2,cmoblk,
     >  ic,ic2)
      implicit real*8 (a-h,o-z)
#include "mo_cvb.fh"
      dimension cmo(nbas_mo,nbas_mo),cmo2(nbas_mo,nbas_mo)
      dimension cmoblk(nbasisq_mo)

c  Construct full matrix of MOs in symmetry-adapted AO basis :
      call getmoblk_cvb(cmoblk,ic2)
      call fzero(cmo,nbas_mo*nbas_mo)
      do 100 isk=1,nsym_mo
      do 100 jbas=1,nbasi_mo(isk)
100   call fmove_cvb(cmoblk(1+nbassqf_mo(isk)+(jbas-1)*nbasi_mo(isk)),
     >  cmo(nbasf_mo(isk)+1,jbas+nbasf_mo(isk)),nbasi_mo(isk))

      if(mod(ic,2).eq.1)then
        call mxinv_cvb(cmo,nbas_mo)
        call transp_cvb(cmo,cmo,nbas_mo,nbas_mo)
      endif

      if(ic.ge.2)then
        do 200 iorb=1,nact_mo
200     call fmove_cvb(cmo(1,iact_mo(iorb)),cmo2(1,iorb),nbas_mo)
      endif
      return
      end
