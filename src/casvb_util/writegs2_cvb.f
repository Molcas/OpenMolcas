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
      subroutine writegs2_cvb(orbs,cvb,
     >  cvbdet,iapr,ixapr,iabind)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension orbs(norb,norb),cvb(nvb)
      dimension cvbdet(ndetvb),iapr(ndetvb),ixapr(nda+1),iabind(ndetvb)

      call str2vbc_cvb(cvb,cvbdet)
      ioffs=0
      call wris_cvb([ndetvb],1,recn_tmp04,ioffs)
      call wris_cvb([norb],1,recn_tmp04,ioffs)
      call wris_cvb([nalf],1,recn_tmp04,ioffs)
      call wris_cvb([nbet],1,recn_tmp04,ioffs)
      call wrrs_cvb(orbs,norb*norb,recn_tmp04,ioffs)
      idetvb=0
      do 100 ia=1,nda
      do 101 ixa=ixapr(ia),ixapr(ia+1)-1
      idetvb=idetvb+1
      ib=iapr(ixa)
      iabind(idetvb)=ia+(ib-1)*nda
101   continue
100   continue
      call wris_cvb(iabind,ndetvb,recn_tmp04,ioffs)
      call wrrs_cvb(cvbdet,ndetvb,recn_tmp04,ioffs)
      call make_cvb('WRITEGS')
      return
      end
