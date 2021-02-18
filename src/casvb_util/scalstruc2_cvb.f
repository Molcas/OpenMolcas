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
      subroutine scalstruc2_cvb(orbs,cvb,iconfs,ifnss)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "frag_cvb.fh"
      dimension orbs(norb,norb),cvb(nvb)
      dimension iconfs(noe,nconf),ifnss(0:nel,0:nel)

      if(sc)then
        fac=one
        do 100 iorb=1,norb
        fac2=ddot_(norb,orbs(1,iorb),1,orbs(1,iorb),1)
        fac=fac*sqrt(fac2)
100     continue
        call dscal_(nvb,fac,cvb,1)
      else
        do 200 iorb=1,norb
        fac2=ddot_(norb,orbs(1,iorb),1,orbs(1,iorb),1)
        fac1=sqrt(fac2)
        istr=0
        iconf_off=0
        do 300 ifrag=1,nfrag
        do 301 iS=1,nS_fr(ifrag)
        do 302 ion=0,nel/2
        nelsing=nel-2*ion
        do 400 i=iconf_off+1,iconf_off+nconfion_fr(ion,ifrag)
        if(iconfs(iorb,i).eq.1)then
          call dscal_(ifnss(nelsing,i2s_fr(iS,ifrag)),
     >      fac1,cvb(istr+1),1)
        elseif(iconfs(iorb,i).eq.2)then
          call dscal_(ifnss(nelsing,i2s_fr(iS,ifrag)),
     >      fac2,cvb(istr+1),1)
        endif
        istr=istr+ifnss(nelsing,i2s_fr(iS,ifrag))
400     continue
        iconf_off=iconf_off+nconfion_fr(ion,ifrag)
302     continue
301     continue
300     continue
        if(istr.ne.nvb)then
          write(6,*)' ISTR not equal to NVB in SCALSTRUC! ',istr,nvb
          call abend_cvb()
        endif
200     continue
      endif
      return
      end
