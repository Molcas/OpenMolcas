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
      subroutine ppgs2_cvb(cvb,cvbdet,ifnss)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"
#include "malloc_cvb.fh"
      dimension cvb(nvb),cvbdet(ndetvb),ifnss(0:nel,0:nel)

c  First applicable configuration with first possible spin in
c  each fragment is set to perfect-pairing:
      call dfill(nvb,1d-2,cvb,1)
      ioffs_cvb=0
      icoffs_nconf=0
      do 100 ifrag=1,nfrag
      nelsing=nel_fr(ifrag)-2*mnion_fr(ifrag)
      do 200 iS=1,nS_fr(ifrag)
      if(i2s_fr(iS,ifrag).le.nelsing)then
        cvb(ifnss(nelsing,i2s_fr(iS,ifrag))+ioffs_cvb)=1d0
        goto 300
      endif
200   continue
300   ioffs_cvb=ioffs_cvb+nvb_fr(ifrag)
      icoffs_nconf=icoffs_nconf+nconf_fr(ifrag)
100   continue
      kbasiscvb_kp=kbasiscvb
      kbasiscvb=1
      call str2vbc_cvb(cvb,cvbdet)
      kbasiscvb=kbasiscvb_kp
      call vb2strc_cvb(cvbdet,cvb)
      return
      end
c  Changes phases between alpha-beta separated determinants, and
c  determinants with increasing orbital numbers:
