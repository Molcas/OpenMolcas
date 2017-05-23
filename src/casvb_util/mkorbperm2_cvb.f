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
      subroutine mkorbperm2_cvb(orbs,cvb,owrk,cvbdet)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "spinb_cvb.fh"
      dimension orbs(norb,norb),cvb(nvb)
      dimension owrk(norb,norb),cvbdet(ndetvb)

      if(ip(1).ge.1)then
        write(6,'(/,a)')' Permuting orbitals :'
        write(6,'(1x,30i4)')(iorbprm(iorb),iorb=1,norb)
      endif
      do 100 iorb=1,norb
      jorb=abs(iorbprm(iorb))
      sgn=dble(sign(1,iorbprm(iorb)))
      call fmove(orbs(1,jorb),owrk(1,iorb),norb)
100   call dscal_(norb,sgn,owrk(1,iorb),1)
      call fmove(owrk,orbs,norb*norb)
      call str2vbc_cvb(cvb,cvbdet)
      call permvb_cvb(cvbdet,iorbprm)
      call vb2strc_cvb(cvbdet,cvb)
      return
      end
