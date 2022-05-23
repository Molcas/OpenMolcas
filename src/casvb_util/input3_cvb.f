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
      subroutine input3_cvb(
     >  iorbrel,mxdimrel,ifxorb,ifxstr,
     >  izrstr,iorts,irots,izeta,
     >  ip_iconfs,orbs,irdorbs,ip_cvb,ip_symelm,kbasiscvb_inp)
      implicit real*8 (a-h,o-z)
c ... Files/Hamiltonian available ...
      logical, external :: valid_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "inpmod_cvb.fh"
#include "spinb_cvb.fh"
#include "malloc_cvb.fh"
      parameter (nglob=5,nstrin=51,nendvb=3,nspec=3,
     >  ncrit=2,nmeth=12,nwkw=5,ncmp=4)
      character*8 global,string,endvb,specl
      character*8 crit,weightkw,methkw
      character*50 inpstr
      logical firsttime_cvb
      logical DoCholesky
      external firsttime_cvb
      dimension global(nglob),string(nstrin),endvb(nendvb),specl(nspec)
      dimension crit(ncrit)
      dimension weightkw(nwkw),methkw(nmeth)
      dimension iorbrel(mxdimrel),ifxorb(mxorb)
      dimension iorts(*),irots(*),izeta(*)
      dimension orbs(mxaobf,mxorb),irdorbs(mxorb)
      dimension idum(1)
      save global,string,endvb,specl
      save crit,weightkw,methkw
      data global/'XXXXxxxx','START   ','GUESS   ','PRINT   ',
     >            'PREC    '/
      data string/'XXXXxxxx','XXXXxxxx','SAVE    ','XXXXxxxx',
     >            'ORBPERM ','COUPLE  ','MAXITER ','CRIT    ',
     >            'CASPROJ ','PROJCAS ','NOCASPRO','NOPROJCA',
     >            'XXXXxxxx','XXXXxxxx','SYMELM  ','ORBREL  ',
     >            'XXXXxxxx','SYMPROJ ','NOSYMPRO','FIXORB  ',
     >            'FIXSTRUC','DELSTRUC','FREORB  ','FRESTRUC',
     >            'ORTHCON ','SADDLE  ','SHSTRUC ','VBWEIGHT',
     >            'CIWEIGHT','SCORR   ','NOSCORR ','METHOD  ',
     >            'OPTIM   ','OPT     ','ENDOPTIM','REPORT  ',
     >            'ENDREPOR','XXXXxxxx','TUNE    ','XXXXxxxx',
     >            'OPPOSITE','XXXXxxxx','XXXXxxxx','STAT    ',
     >            'INIT    ','NOINIT  ','TIDY    ','PLOC    ',
     >            'NOPLOC  ','ALTERNAT','ENDALTER'/
      data endvb/ 'ENDVB   ','ENDCASVB','END     '/
      data specl/ 'SERVICE ','MOSCOW  ','PERFLOC '/
      data crit/  'OVERLAP ','ENERGY  '/
      data weightkw/'CHIRGWIN','LOWDIN  ','INVERSE ','NONE    ',
     >              'ALL     '/
      data methkw/'FLETCHER','TRIM    ','TRUSTOPT','DAVIDSON',
     >            'STEEP   ','VB2CAS  ','AUGHESS ','AUG2    ',
     >            'CHECK   ','DFLETCH ','NONE    ','SUPER   '/


      Call DecideOnCholesky(DoCholesky)
      If (DoCholesky) Then
        write(6,*)'** Cholesky or RI/DF not yet implemented in CASVB **'
        call abend_cvb()
      EndIf

      call fstring_cvb(specl,nspec,istr,ncmp,2)
      if(istr.eq.1)then
c 'SERVICE'
        service=.true.
        write(6,'(1x,a,/)') '**** Service mode **** '
        call service_cvb()
        return
      elseif(istr.eq.2)then
c 'MOSCOW'
        service=.true.
        write(6,'(1x,a,/)') '**** MOSCOW mode **** '
        call moscow_cvb()
        return
      elseif(istr.eq.3)then
        service=.true.
        write(6,'(1x,a,/)') '**** PERFLOC mode **** '
        !call perfloc_plc(3)
        call perfloc_plc()
        return
      endif

1000  continue

c  CASSCF wavefunction information :
      call casinfoinp_cvb()
c  VB wavefunction information :
      call fraginp_cvb(ip_iconfs)

      igroup=0
      call fstring_cvb(global,nglob,istr,ncmp,2)
      if(istr.ne.0)then
        igroup=1
        goto 1110
      endif
      call fstring_cvb(string,nstrin,istr,ncmp,2)
      if(istr.ne.0)then
        igroup=2
        goto 1110
      endif
      call fstring_cvb(endvb,nendvb,istr,ncmp,2)
      if(istr.ne.0)then
        igroup=3
        goto 1110
      endif
1110  continue
      if(igroup.eq.3)then
c 'ENDVB', 'ENDCASVB' or 'END'
        istr=0
      endif

      if(igroup.eq.2)goto 1111
      if(istr.eq.2)then
c 'START'
        strtvb=zero
520     call string_cvb(inpstr,1,nread,1)
        if(nread.eq.1)then
          if(inpstr(1:3).eq.'CI=')then
            call setfn_cvb(strtci,inpstr(4:50))
            goto 520
          elseif(inpstr(1:3).eq.'VB=')then
            call setfn_cvb(strtvb,inpstr(4:50))
            goto 520
          elseif(inpstr(1:3).eq.'MO=')then
            call setfn_cvb(strtmo,inpstr(4:50))
            goto 520
          elseif(inpstr(1:4).eq.'INT=')then
            call setfn_cvb(strtint,inpstr(5:50))
            goto 520
          endif
        endif
        if(valid_cvb(strtvb).and.firsttime_cvb())
     >    call touch_cvb('STRTGS')
      elseif(istr.eq.3)then
c 'GUESS'
        call gsinp_cvb(
     >    orbs,irdorbs,ip_cvb,nvbinp,kbasiscvb_inp,
     >    mxaobf,mxorb,kbasis,strtvb)
      elseif(istr.eq.4)then
c 'PRINT'
        call int_cvb(ip,10,nread,1)
      elseif(istr.eq.5)then
c 'PREC'
        call int_cvb(idum,1,nread,1)
        iprec=idum(1)
        if(iprec.lt.0)then
          write(6,*)' Illegal precision :',iprec
          call abend_cvb()
        endif
        call int_cvb(idum,1,nread,1)
        iwidth=idum(1)
        call formats_cvb()
      endif

c 'ENDVB', 'ENDCASVB' , 'END' or unrecognized keyword -- end of input :
      if(istr.ne.0)goto 1000
1111  continue

      if(istr.eq.1)then
      elseif(istr.eq.2)then
      elseif(istr.eq.3)then
c 'SAVE'
1520    call string_cvb(inpstr,1,nread,1)
        if(nread.eq.1)then
          if(inpstr(1:5).eq.'VBCI=')then
            call setfn_cvb(savvbci,inpstr(6:50))
            goto 1520
          elseif(inpstr(1:3).eq.'VB=')then
            call setfn_cvb(savvb,inpstr(4:50))
            goto 1520
          endif
        endif
      elseif(istr.eq.4)then
      elseif(istr.eq.5)then
c 'ORBPERM'
        if(firsttime_cvb())call touch_cvb('ORBPERM')
        call int_cvb(iorbprm,mxorb,nread,0)
        if(nread.gt.mxorb)then
          write(6,*)' Too many orbitals in ORBPERM keyword!'
          call abend_cvb()
        endif
        do 13350 iorb=1,nread
        if(abs(iorbprm(iorb)).lt.1.or.abs(iorbprm(iorb)).gt.mxorb)then
          write(6,'(a,40i3)')' Illegal orbital label(s) in ORBPERM:',
     >      (iorbprm(ior),ior=1,nread)
          call abend_cvb()
        endif
13350   continue
      elseif(istr.eq.6)then
c 'COUPLE'
        kbasis_old=kbasis
        call fstring_cvb(spinbkw,nspinb,kbasis,ncmp,1)
        if(kbasis.eq.0)kbasis=kbasis_old
        if(kbasis.eq.7)kbasis=6
      elseif(istr.eq.7)then
c  'MAXITER'
        call int_cvb(idum,1,nread,0)
        mxiter=idum(1)
      elseif(istr.eq.8)then
c 'CRIT'
        call fstring_cvb(crit,ncrit,icrit,ncmp,1)
        if(icrit.ne.1.and.icrit.ne.2)then
          write(6,*)' Unrecognized CRIT keyword!'
          call abend_cvb()
        endif
      elseif(istr.eq.9.or.istr.eq.10)then
c 'CASPROJ' or 'PROJCAS'
        projcas=.true.
      elseif(istr.eq.11.or.istr.eq.12)then
c 'NOCASPROJ' or 'NOPROJCAS'
        projcas=.false.
      elseif(istr.eq.15)then
c 'SYMELM'
        call symelminp_cvb(ip_symelm,nsyme,tags,izeta,
     >    mxirrep,mxorb,mxsyme,ityp)
      elseif(istr.eq.16)then
c 'ORBREL'
        iorb=0
        jorb=0
        call int_cvb(idum,1,nread,1)
        iorb=idum(1)
        call int_cvb(idum,1,nread,1)
        jorb=idum(1)
        if(iorb.lt.1.or.iorb.gt.mxorb.or.jorb.lt.1.or.jorb.gt.mxorb)then
          write(6,*)' Illegal orbital number(s) in ORBREL:',iorb,jorb
          call abend_cvb()
        endif
        iorbrel(1+ndimrel)=iorb
        iorbrel(2+ndimrel)=jorb
        nops=0
15300   call fstring_cvb(tags,nsyme,itag,3,1)
        if(itag.ne.0)then
          nops=nops+1
          if(ndimrel+3+nops.gt.mxdimrel)then
            write(6,*)' Too many symmetry elements in ORBREL keyword!'
            call abend_cvb()
          endif
          iorbrel(nops+3+ndimrel)=itag
          goto 15300
        endif
        iorbrel(3+ndimrel)=nops
        norbrel=norbrel+1
        ndimrel=ndimrel+3+nops
      elseif(istr.eq.18)then
c 'SYMPROJ'
        projsym=.true.
        call izero(isympr,mxirrep)
        call int_cvb(idum,1,nread,1)
        isymput=idum(1)
        if(nread.eq.1)then
          isympr(isymput)=1
15320     call int_cvb(idum,1,nread,1)
          isymput=idum(1)
          if(nread.eq.1)then
            isympr(isymput)=1
            goto 15320
          endif
        else
          call imove_cvb(isymv,isympr,mxirrep)
        endif
      elseif(istr.eq.19)then
c 'NOSYMPROJ'
        projsym=.false.
      elseif(istr.eq.20)then
c 'FIXORB'
        itmp = mstacki_cvb(mxorb)
        call intchk_cvb(iw(itmp),mxorb,nfxorb,0,'FIXORB',-1)
        call izero(ifxorb,mxorb)
        do 15340 i=1,nfxorb
        ifxorb(iw(i+itmp-1))=1
15340   continue
        call mfreei_cvb(itmp)
      elseif(istr.eq.21)then
c 'FIXSTRUC'
        lfxvb=0
        call mhpfreei_cvb(ifxstr)
        mxread=mavaili_cvb()/2
        ifxstr=mheapi_cvb(mxread)
        call intchk_cvb(iw(ifxstr),mxread,nfxvb,0,'FIXSTRUC',lfxvb)
        call mrealloci_cvb(ifxstr,nfxvb)
      elseif(istr.eq.22)then
c 'DELSTRUC'
        lzrvb=0
        call mhpfreei_cvb(izrstr)
        mxread=mavaili_cvb()/2
        izrstr=mheapi_cvb(mxread)
        call intchk_cvb(iw(izrstr),mxread,nzrvb,0,'DELSTRUC',lzrvb)
        call mrealloci_cvb(izrstr,nzrvb)
      elseif(istr.eq.23)then
c 'FREORB' - not implemented
        itmp = mstacki_cvb(mxorb)
        call intchk_cvb(iw(itmp),mxorb,nfrorb1,0,'FREORB',-1)
        itmp2 = mstackiz_cvb(mxorb)
        do i=1,nfrorb1
        iw(iw(i+itmp-1)+itmp2-1)=1
        enddo
        nfxorb1=0
        do i=1,mxorb
        if(iw(i+itmp2-1).eq.1)then
          nfxorb1=nfxorb1+1
          iw(nfxorb1+itmp-1)=i
        endif
        enddo
        call mfreei_cvb(itmp)
        nfxorb=max(nfxorb,nfxorb1)
      elseif(istr.eq.24)then
c 'FRESTRUC' - not implemented (and code incomplete)
        lfxvb=1
        call mhpfreei_cvb(ifxstr)
        mxread=mavaili_cvb()/2
        ifxstr=mheapi_cvb(mxread)
        call intchk_cvb(iw(ifxstr),mxread,nfxvb,0,'FRESTRUC',lfxvb)
        call mrealloci_cvb(ifxstr,nfxvb)
      elseif(istr.eq.25)then
c 'ORTHCON'
        mxgroup=40
        mxortl=40
        mxpair=mxorb*(mxorb-1)/2
        itmpa = mstacki_cvb(mxorb*mxorb)
        itmpb = mstacki_cvb(mxorb*mxgroup)
        itmpc = mstacki_cvb(mxgroup)
        itmpd = mstacki_cvb(mxortl)
        call orthcon_cvb(iorts,iw(itmpa),iw(itmpb),iw(itmpc),
     >  iw(itmpd),mxortl,mxpair)
        call mfreei_cvb(itmpa)
      elseif(istr.eq.26)then
c 'SADDLE'
        call int_cvb(idum,1,nread,1)
        isaddle=idum(1)
      elseif(istr.eq.27)then
c 'SHSTRUC'
        ishstruc=1
      elseif(istr.eq.28)then
c 'VBWEIGHT'
        ivbweights=0
15600   call fstring_cvb(weightkw,nwkw,istr2,ncmp,1)
        if(istr2.eq.1)then
          if(mod(ivbweights,2).eq.0)ivbweights=ivbweights+1
        elseif(istr2.eq.2)then
          if(mod(ivbweights,4).le.1)ivbweights=ivbweights+2
        elseif(istr2.eq.3)then
          if(mod(ivbweights,8).le.3)ivbweights=ivbweights+4
        elseif(istr2.eq.4)then
          ivbweights=0
        elseif(istr2.eq.5)then
          ivbweights=7
        endif
        if(istr2.gt.0)goto 15600
      elseif(istr.eq.29)then
c 'CIWEIGHT'
        npcf=10
        iciweights=0
15700   call fstring_cvb(weightkw,nwkw,istr2,ncmp,1)
        if(istr2.eq.1)then
          if(mod(iciweights,2).eq.0)iciweights=iciweights+1
        elseif(istr2.eq.2)then
          if(mod(iciweights,4).le.1)iciweights=iciweights+2
        elseif(istr2.eq.3)then
          if(mod(iciweights,8).le.3)iciweights=iciweights+4
        elseif(istr2.eq.4)then
          iciweights=0
        elseif(istr2.eq.5)then
          iciweights=7
        endif
        if(istr2.gt.0)goto 15700
        call int_cvb(idum,1,nread,1)
        npcf=idum(1)
      elseif(istr.eq.30)then
c 'SCORR'
        sij=.true.
      elseif(istr.eq.31)then
c 'NOSCORR'
        sij=.false.
      elseif(istr.eq.32)then
c 'METHOD'
        call fstring_cvb(methkw,nmeth,istr2,ncmp,1)
        if(istr2.ne.0)then
          imethod=istr2
        endif
      elseif(istr.eq.33.or.istr.eq.34)then
c 'OPTIM' or 'OPT'
        call maxdims_cvb()
        call loopcntr_cvb(1)
      elseif(istr.eq.35)then
c 'ENDOPTIM'
        call maxdims_cvb()
        call loopcntr_cvb(2)
      elseif(istr.eq.36)then
c 'REPORT '
        call maxdims_cvb()
        call loopcntr_cvb(3)
      elseif(istr.eq.37)then
c 'ENDREPOR'
        call maxdims_cvb()
        call loopcntr_cvb(4)
      elseif(istr.eq.39)then
c 'TUNE'
        call tuneinp_cvb()
      elseif(istr.eq.41)then
c 'OPPOSITE'
        opposite=.true.
      elseif(istr.eq.44)then
c 'STAT'
        if(firsttime_cvb())call touch_cvb('STAT')
      elseif(istr.eq.45)then
c 'INIT'
        initial=1
      elseif(istr.eq.46)then
c 'NOINIT'
        initial=0
      elseif(istr.eq.47)then
c 'TIDY'
      elseif(istr.eq.48)then
c 'PLOC'
        ploc=.true.
      elseif(istr.eq.49)then
c 'NOPLOC'
        ploc=.false.
      elseif(istr.eq.50)then
c 'ALTERN'
        mxalter=50
        call int_cvb(idum,1,nread,1)
        mxalter=idum(1)
        call loopcntr2_cvb(5,mxalter)
      elseif(istr.eq.51)then
c 'ENDALTER'
        call loopcntr_cvb(6)
      endif

c 'ENDVB', 'ENDCASVB' , 'END' or unrecognized keyword -- end of input :
      if(istr.ne.0)goto 1000
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(irots)
      end
