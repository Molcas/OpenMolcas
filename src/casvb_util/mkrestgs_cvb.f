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
      subroutine mkrestgs_cvb(orbsao,irdorbs,cvb,
     >  cvbdet,iapr,ixapr,iabind,cvbdet1)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "mo_cvb.fh"
      dimension orbsao(nbas_mo,norb),irdorbs(norb),cvb(nvb)
      dimension cvbdet(ndetvb),iapr(ndetvb),ixapr(nda+1)
      dimension iabind(*),cvbdet1(*)
      dimension idum(1)

      ioffs=0
      call rdis_cvb(idum,1,recn_tmp04,ioffs)
      ndetvb1=idum(1)
      call rdis_cvb(idum,1,recn_tmp04,ioffs)
      norb1=idum(1)
      call rdis_cvb(idum,1,recn_tmp04,ioffs)
      nalf1=idum(1)
      call rdis_cvb(idum,1,recn_tmp04,ioffs)
      nbet1=idum(1)
      if(norb1.ne.norb.or.nalf1.ne.nalf.or.nbet1.ne.nbet)then
        write(6,'(2a)')' Inconsistency between previous and current',
     >    ' VB wavefunction definitions.'
        write(6,*)' NORB now ',norb,' before ',norb1
        write(6,*)' NALF now ',nalf,' before ',nalf1
        write(6,*)' NBET now ',nbet,' before ',nbet1
        call abend_cvb()
      endif
      do 100 iorb=1,norb
      irdorbs(iorb)=1
      call rdrs_cvb(orbsao(1,iorb),norb,recn_tmp04,ioffs)
100   continue
      call rdis_cvb(iabind,ndetvb1,recn_tmp04,ioffs)
      call rdrs_cvb(cvbdet1,ndetvb1,recn_tmp04,ioffs)

      call fzero(cvbdet,ndetvb)
      do 200 idetvb1=1,ndetvb1
c  NDA & string definitions assumed the same :
      ib=(iabind(idetvb1)-1)/nda+1
      ia=iabind(idetvb1)-(ib-1)*nda
      do 300 ixa=ixapr(ia),ixapr(ia+1)-1
      if(ib.eq.iapr(ixa))cvbdet(ixa)=cvbdet1(idetvb1)
300   continue
200   continue
      kbasiscvb=kbasis
      call vb2strc_cvb(cvbdet,cvb)
      return
      end
