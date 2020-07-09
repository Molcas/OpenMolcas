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
      subroutine input2_cvb(
     >  iorbrel,mxdimrel,ifxorb,
     >  iorts,irots,izeta,orbs,irdorbs)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "inpmod_cvb.fh"
#include "spinb_cvb.fh"
#include "frag_cvb.fh"
#include "malloc_cvb.fh"
      dimension iorbrel(mxdimrel),ifxorb(mxorb)
      dimension iorts(2,*),irots(2,*),izeta(*)
      dimension orbs(mxaobf,mxorb),irdorbs(mxorb)

      ibase = mstacki_cvb(0)
      ip_iconfs = mheapiz_cvb(0)
      ip_cvb = mheaprz_cvb(0)
      ip_symelm = mheaprz_cvb(0)
      ifxstr = mheapiz_cvb(0)
      idelstr = mheapiz_cvb(0)
      noe=2*mxorb

      call hini_cvb()

      call defs_cvb()
      call casinfodef_cvb()

c  Counters
      nconf=0
      nvbinp=0
      nsyme=0
      norbrel=0
      ndimrel=0
      nijrel=0
      nfxorb=0
      nfxvb=0
      nzrvb=0
      nort=0
      ndrot=0
      lfxvb=0
      lzrvb=0
      call maxdims0_cvb()
      call izero(ifxorb,mxorb)
      call izero(izeta,mxsyme)
      call fraginit_cvb()

      call input3_cvb(
     >  iorbrel,mxdimrel,ifxorb,ifxstr,
     >  idelstr,iorts,irots,izeta,
     >  ip_iconfs,orbs,irdorbs,ip_cvb,ip_symelm,kbasiscvb_inp)

      if(inputmode.eq.2)then
c  Input parsing complete for this step ...
c  ... Work out NORB, NEL, S ...
        noe1=noe
        call casinfoset_cvb()
c  ... Do ICONFS before others to get NVB and related info ...
        do 100 iconf=1,nconf
        call imove_cvb(iw((iconf-1)*noe1+ip_iconfs),
     >    iw((iconf-1)*noe+ip_iconfs),noe)
100     continue
        call mrealloci_cvb(ip_iconfs,noe*nconf)

        if(nfrag.le.1)then
          nMs_fr(1)=1
          nalf_fr(1,1)=nalf
          nbet_fr(1,1)=nbet
        else
          do ifrag=1,nfrag
          nMs_fr(ifrag)=1
          nalf_fr(1,ifrag)=(nel_fr(ifrag)+i2s_fr(1,ifrag))/2
          nbet_fr(1,ifrag)=nel_fr(ifrag)-nalf_fr(1,ifrag)
          enddo
        endif
        if(nfrag.eq.0)then
          nfrag=1
          nel_fr(1)=nel
          nconf_fr(1)=nconf
          nS_fr(1)=1
          i2s_fr(1,1)=nalf-nbet
        endif
        do ifrag=1,nfrag
        if(nS_fr(ifrag).eq.0)then
          nS_fr(ifrag)=1
          i2s_fr(1,ifrag)=nalf-nbet
        endif
        enddo
        iconf_add=0
        do 1001 ifrag=1,nfrag
        if(nel_fr(ifrag).eq.0)then
          nel_fr(ifrag)=nel
          nalf_fr(1,ifrag)=nalf
          nbet_fr(1,ifrag)=nbet
        endif
        if(nS_fr(ifrag).eq.0)then
          nS_fr(ifrag)=1
          i2s_fr(1,ifrag)=nalf-nbet
        endif
        if(nconf_fr(ifrag).eq.0)then
          nconf=nconf+1
          nconf_fr(ifrag)=1
          call mrealloci_cvb(ip_iconfs,noe*nconf)
          do jconf=nconf,iconf_add+2,-1
          call imove_cvb(iw((jconf-2)*noe+ip_iconfs),
     >      iw((jconf-1)*noe+ip_iconfs),noe)
          enddo
          call izero(iw(iconf_add*noe+ip_iconfs),noe)
          do 1201 i=1,min(nel_fr(ifrag),norb)
          iw(i+iconf_add*noe+ip_iconfs-1)=1
1201      continue
          do 1301 i=1,nel_fr(ifrag)-norb
          iw(i+iconf_add*noe+ip_iconfs-1)=2
1301      continue
        endif
        call cnfcheck_cvb(iw(iconf_add*noe+ip_iconfs),nconf_fr(ifrag),
     >    nel_fr(ifrag))
        call cnfini_cvb(iw(iconf_add*noe+ip_iconfs),nconf_fr(ifrag),
     >    nel_fr(ifrag),
     >    nS_fr(ifrag),i2s_fr(1,ifrag),
     >    nMs_fr(ifrag),nalf_fr(1,ifrag),nbet_fr(1,ifrag),
     >    nvbr_fr(ifrag),ndetvb_fr(ifrag),ndetvb2_fr(ifrag),
     >    mnion_fr(ifrag),mxion_fr(ifrag),nconfion_fr(0,ifrag),ifsc)
        iconf_add=iconf_add+nconf_fr(ifrag)
1001    continue
        ndetvb=0
        ndetvb2=0
        nvbr=0
        nelcheck=0
        do i=1,nfrag
        ndetvb=ndetvb+ndetvb_fr(i)
        ndetvb2=ndetvb2+ndetvb2_fr(i)
        nvbr=nvbr+nvbr_fr(i)
        nelcheck=nelcheck+nel_fr(i)
        enddo
        if(nelcheck.ne.nel)then
          write(6,*)' Error: total number of electrons in fragment ',
     >      'wavefunctions :',nelcheck,
     >      ' not equal to number of electrons ',nel
          call abend_cvb()
        endif
        sc=(nfrag.eq.1.and.ifsc.eq.1)
c  Set absym and use just lowest spin value if spinbas=determinants :
        absym(1)=(nalf.eq.nbet)
        do ifrag=1,nfrag
        i2s_min=nel_fr(ifrag)
        do iS=1,nS_fr(ifrag)
        if(i2s_fr(iS,ifrag).ne.0)absym(1)=.false.
        is2_min=min(i2s_min,i2s_fr(iS,ifrag))
        enddo
        if(kbasis.eq.6)then
          nS_fr(ifrag)=1
          i2s_fr(1,ifrag)=i2s_min
        endif
        enddo
        do i=2,5
        absym(i)=absym(1)
        enddo
        nvb=nvb_cvb(kbasis)
        mnion=mnion_fr(1)
        mxion=mxion_fr(1)
        do i=2,nfrag
        mnion=min(mnion,mnion_fr(i))
        mxion=max(mxion,mxion_fr(i))
        enddo
c  ... Now remaining quantities that depend on NORB or NVB ...
c  SYMELM
        ip_from=ip_symelm
        ip_to=ip_symelm
        do 400 isyme=1,nsyme
        do 500 iorb=1,norb
        if(ip_from.ne.ip_to)call fmove_cvb(w(ip_from),w(ip_to),norb)
        ip_from=ip_from+mxorb
        ip_to=ip_to+norb
500     continue
        ip_from=ip_from+(mxorb-norb)*mxorb
400     continue
c  IORBREL
        ifrom=1
        ito=1
600     continue
        if(ifrom.le.ndimrel)then
          iorb=iorbrel(ifrom)
          jorb=iorbrel(ifrom+1)
          nmov=3+iorbrel(ifrom+2)
          if(iorb.le.norb.and.jorb.le.norb)then
            if(ifrom.ne.ito)
     *            call imove_cvb(iorbrel(ifrom),iorbrel(ito),nmov)
            ito=ito+nmov
          endif
          ifrom=ifrom+nmov
          goto 600
        endif
        ndimrel=ito-1
c  IFXSTR
        ito=0
        do 700 ifrom=1,nfxvb
        if(iw(ifrom+ifxstr-1).le.nvb)then
          ito=ito+1
          iw(ito+ifxstr-1)=iw(ifrom+ifxstr-1)
        endif
700     continue
        nfxvb=ito
c  IDELSTR
        ito=0
        do 800 ifrom=1,nzrvb
        if(iw(ifrom+idelstr-1).le.nvb)then
          ito=ito+1
          iw(ito+idelstr-1)=iw(ifrom+idelstr-1)
        endif
800     continue
        nzrvb=ito
c  IORTS
        ito=0
        do 900 ifrom=1,nort
        if(iorts(1,ifrom).le.norb.and.iorts(2,ifrom).le.norb)then
          ito=ito+1
          iorts(1,ito)=iorts(1,ifrom)
          iorts(2,ito)=iorts(2,ifrom)
        endif
900     continue
        nort=ito
c  IROTS
        ito=0
        do 1000 ifrom=1,ndrot
        if(irots(1,ifrom).le.norb.and.irots(2,ifrom).le.norb)then
          ito=ito+1
          irots(1,ito)=irots(1,ifrom)
          irots(2,ito)=irots(2,ifrom)
        endif
1000    continue
        ndrot=ito
c  Calling DEFS2 before INITOPT is required in order to set things
c  such as ICRIT :
        call defs2_cvb(ifxorb)
        call initopt_cvb(icrit,lfxvb,nfxvb,iorts,nort,norb)

        call defs2_cvb(ifxorb)

c Try for new record
        call rdioff1_cvb(need)
        need=need+3*ihlf_cvb(1)+ihlf_cvb(noe*nconf)+
     >    mxaobf*norb+ihlf_cvb(norb)+nvbinp+nsyme*norb*norb+
     >    ihlf_cvb(ndimrel)+
     >    ihlf_cvb(norb)+ihlf_cvb(nfxvb)+ihlf_cvb(nzrvb)+
     >    ihlf_cvb(2*nort)+ihlf_cvb(2*ndrot)+ihlf_cvb(2*ndrot)+
     >    ihlf_cvb(nsyme)
        if(recinp.eq.0d0)then
          recinp=recn_tmp01
        elseif(recinp_old.eq.0d0)then
          recinp_old=recn_tmp01
          recinp=recn_tmp02
        else
          swap=recinp_old
          recinp_old=recinp
          recinp=swap
        endif
        call reserv_cvb(need,recinp)
        call rdioff1_cvb(ioffs)
        call wrioff_cvb(1,recinp,ioffs)
        call wris_cvb([noe],1,recinp,ioffs)
        call wrioff_cvb(2,recinp,ioffs)
        call wris_cvb([nconf],1,recinp,ioffs)
        call wrioff_cvb(3,recinp,ioffs)
        call wris_cvb([kbasiscvb_inp],1,recinp,ioffs)
        call wrioff_cvb(4,recinp,ioffs)
        call wris_cvb(iw(ip_iconfs),noe*nconf,recinp,ioffs)
        call wrioff_cvb(5,recinp,ioffs)
        call wrrs_cvb(orbs,mxaobf*norb,recinp,ioffs)
        call wrioff_cvb(6,recinp,ioffs)
        call wris_cvb(irdorbs,norb,recinp,ioffs)
        call wrioff_cvb(7,recinp,ioffs)
        call wrrs_cvb(w(ip_cvb),nvbinp,recinp,ioffs)
        call wrioff_cvb(8,recinp,ioffs)
        call wrrs_cvb(w(ip_symelm),nsyme*norb*norb,recinp,ioffs)
        call wrioff_cvb(9,recinp,ioffs)

        call dset_cvb(iorbrel,ifxorb,ifxstr,
     >    idelstr,iorts,irots,izeta)
      else
        call maxdims_cvb()
      endif
      call mhpfreei_cvb(ip_iconfs)
      call mhpfreer_cvb(ip_cvb)
      call mhpfreer_cvb(ip_symelm)
      call mhpfreei_cvb(ifxstr)
      call mhpfreei_cvb(idelstr)
      call hend_cvb()
      call mfreei_cvb(ibase)

      return
      end
