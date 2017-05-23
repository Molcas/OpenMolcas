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
      subroutine ciweight2_cvb(civec,civbs,civb,citmp,civec5,
     >  orbs,sorbs,orbinv,owrk,gjorb,gjorb2,gjorb3,
     >  vec1,vec2,vec3,
     >  vec4,vec5,
     >  wghtion1,wghtion2,wghtion3,wghtion4,wghtion5,wghtion6,
     >  mingrph,maxgrph,xalf,xbet,iaocc,ibocc,
     >  mingion,maxgion,nkion,xion,locion,lunion,
     >  mingsng,maxgsng,nksng,xsng,locsng,lunsng,
     >  mingasg,maxgasg,nkasg,xasg,locasg,lunasg,
     >  gal1,gal2,indavec,indbvec,
     >  ionmin,ionmax,mxrem,mxsng,mxasg,ncnfcas,mxdetcas)
      implicit real*8 (a-h,o-w,y-z),integer(x)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "formats_cvb.fh"
      character*240 line
      dimension civec(ndet),civbs(ndet),civb(ndet),citmp(ndet),
     >  civec5(ndet)
      dimension orbs(norb,norb),sorbs(norb,norb)
      dimension orbinv(norb,norb),owrk(norb,norb)
      integer gjorb,gjorb2,gjorb3
      dimension gjorb(*),gjorb2(*),gjorb3(*)

      dimension vec1(ndet),vec2(ndet),vec3(ndet),vec4(ndet),vec5(ndet)
      dimension wghtion1(ionmin:ionmax),wghtion2(ionmin:ionmax)
      dimension wghtion3(ionmin:ionmax),wghtion4(ionmin:ionmax)
      dimension wghtion5(ionmin:ionmax),wghtion6(ionmin:ionmax)
      dimension mingrph(0:norb),maxgrph(0:norb)
      dimension xalf(0:norb,0:nalf),xbet(0:norb,0:nbet)
      dimension iaocc(norb),ibocc(norb)
      dimension mingion(0:norb),maxgion(0:norb),nkion(0:norb)
      dimension xion((norb+1)*(ionmax+1)),locion(norb),lunion(norb)
      dimension mingsng(0:norb),maxgsng(0:norb),nksng(0:norb)
      dimension xsng((mxrem+1)*(mxsng+1)),locsng(norb),lunsng(norb)
      dimension mingasg(0:norb),maxgasg(0:norb),nkasg(0:norb)
      dimension xasg((mxsng+1)*(mxasg+1)),locasg(norb),lunasg(norb)
      dimension gal1(ncnfcas),gal2(ncnfcas)
      dimension indavec(mxdetcas),indbvec(mxdetcas)

      dimension cprint(6)

      call cidot_cvb(civb,civbs,cnrm)
      fac=svb/sqrt(cnrm)

      call cicopy_cvb(civec,citmp)
      call fmove(orbs,orbinv,norb*norb)
      call mxinv_cvb(orbinv,norb)
      call gaussj_cvb(orbinv,gjorb)
      call applyt_cvb(civec,gjorb)
c  Chirgwin-Coulson weights
      if(mod(iciweights,2).eq.1)then
        call transp_cvb(orbs,owrk,norb,norb)
        call gaussj_cvb(owrk,gjorb2)
        call applyt_cvb(citmp,gjorb2)
        do 200 idet=1,ndet
200     vec2(idet)=(vec1(idet)-fac*vec2(idet))*
     >              (vec3(idet)-fac*vec4(idet))
        do 300 idet=1,ndet
300     vec1(idet)=vec1(idet)*vec3(idet)
      endif
      do 400 idet=1,ndet
400   vec4(idet)=vec3(idet)-fac*vec4(idet)

c  Inverse-overlap weights
      if(.not.mod(iciweights,8).gt.3)goto 4010
      call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)
      call fmove(sorbs,orbinv,norb*norb)
      call mxinv_cvb(orbinv,norb)
      call gaussj_cvb(orbinv,gjorb)
c Alpha weight array:
      do 1100 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nalf,0)
1100  maxgrph(iorb)=min(iorb,nalf)
      call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
c Beta weight array:
      do 1200 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nbet,0)
1200  maxgrph(iorb)=min(iorb,nbet)
      call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)

      nc=0
      do 2000 ion=ionmin,ionmax
      nsing=nel-2*ion
      mrem=norb-ion
      nalfsng=nalf-ion
      nbetsng=nbet-ion
c  Initialise loop for ionic orbitals
      do 2100 iorb=0,norb
      mingion(iorb)=max(iorb-norb+ion,0)
2100  maxgion(iorb)=min(iorb,ion)
      call weight_cvb(xion,mingion,maxgion,ion,norb)
      call imove_cvb(maxgion,nkion,norb+1)
      call occupy_cvb(nkion,norb,locion,lunion)
c  Initialise loop for singly occupied orbitals
      do 2200 iorb=0,mrem
      mingsng(iorb)=max(iorb-mrem+nsing,0)
2200  maxgsng(iorb)=min(iorb,nsing)
      call weight_cvb(xsng,mingsng,maxgsng,nsing,mrem)
      call imove_cvb(maxgsng,nksng,mrem+1)
      call occupy_cvb(nksng,mrem,locsng,lunsng)
c  Initialise loop for singly occupied alpha orbitals
      do 2300 iorb=0,nsing
      mingasg(iorb)=max(iorb-nsing+nalfsng,0)
2300  maxgasg(iorb)=min(iorb,nalfsng)
      call weight_cvb(xasg,mingasg,maxgasg,nalfsng,nsing)
      call imove_cvb(maxgasg,nkasg,nsing+1)
      call occupy_cvb(nkasg,nsing,locasg,lunasg)

c  Loop ionic
      indion=1
3000  continue
c  Loop singly occupied
      indsng=1
3100  continue
      call fzero(vec5,ndet)
      s11=zero
      s22=zero
      s12=zero
c  Loop singly occupied alpha
      indasg=1
3200  continue

      call izero(iaocc,norb)
      call izero(ibocc,norb)
      do 3300 i=1,ion
      iaocc(locion(i))=1
3300  ibocc(locion(i))=1
      do 3400 ia=1,nalfsng
      iaorb=lunion(locsng(locasg(ia)))
3400  iaocc(iaorb)=1
      do 3500 ib=1,nbetsng
      iborb=lunion(locsng(lunasg(ib)))
3500  ibocc(iborb)=1
      inda=indget_cvb(iaocc,nalf,norb,xalf)
      indb=indget_cvb(ibocc,nbet,norb,xbet)
      indavec(indasg)=inda
      indbvec(indasg)=indb

      indab=(indb-1)*nda+inda
      vec5(indab)=vec3(indab)
      s11=s11+vec3(indab)*vec3(indab)
      s22=s22+vec4(indab)*vec4(indab)
      s12=s12+vec3(indab)*vec4(indab)

      call loind_cvb(nsing,nalfsng,nkasg,mingasg,maxgasg,
     >                       locasg,lunasg,indasg,xasg,*3200)

      call applyt_cvb(civec5,gjorb)

      call icomb_cvb(nsing,nalfsng,nindasg)
      sm1=zero
      do 3600 indasg=1,nindasg
      inda=indavec(indasg)
      indb=indbvec(indasg)
      indab=(indb-1)*nda+inda
3600  sm1=sm1+vec3(indab)*vec5(indab)

      if(abs(sm1).gt.1.d-20)then
        sm1=s11*s11/sm1
      elseif(abs(sm1).le.1.d-20)then
        sm1=zero
      endif

      call fzero(vec5,ndet)
      do 3700 indasg=1,nindasg
      inda=indavec(indasg)
      indb=indbvec(indasg)
      indab=(indb-1)*nda+inda
3700  vec5(indab)=vec4(indab)

      call applyt_cvb(civec5,gjorb)

      sm2=zero
      do 3800 indasg=1,nindasg
      inda=indavec(indasg)
      indb=indbvec(indasg)
      indab=(indb-1)*nda+inda
3800  sm2=sm2+vec4(indab)*vec5(indab)

      if(abs(sm2).gt.1.d-20)then
        sm2=s22*s22/sm2
      elseif(abs(sm2).le.1.d-20)then
        sm2=zero
      endif

      nc=nc+1
      gal1(nc)=sm1
      gal2(nc)=sm2

      call loind_cvb(mrem,nsing,nksng,mingsng,maxgsng,
     >                       locsng,lunsng,indsng,xsng,*3100)
      call loind_cvb(norb,ion,nkion,mingion,maxgion,
     >                       locion,lunion,indion,xion,*3000)
2000  continue
      sum1=zero
      sum2=zero
      do 3900 ic=1,ncnfcas
      sum1=sum1+gal1(ic)
3900  sum2=sum2+gal2(ic)
      fac1=one/sum1
      if(abs(one-svb*svb).lt.1.d-20.and.abs(sum2).lt.1.d-20)then
        fac2=one
      else
        fac2=(one-svb*svb)/sum2
      endif
      do 4000 ic=1,ncnfcas
      gal1(ic)=fac1*gal1(ic)
4000  gal2(ic)=fac2*gal2(ic)
4010  continue

c  Weights of Lowdin orthonormalized structures
      if(mod(iciweights,4).gt.1)then
        call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)
        call mxsqrt_cvb(sorbs,norb,1)
        call gaussj_cvb(sorbs,gjorb3)
        call applyt_cvb(civec,gjorb3)
        call applyt_cvb(civb,gjorb3)
        call cidot_cvb(civb,civb,cnrm)
        fac=svb/sqrt(cnrm)
        do 4100 idet=1,ndet
        vec4(idet)=vec4(idet)*vec4(idet)
4100    vec3(idet)=vec3(idet)*vec3(idet)
      endif

      write(6,'(/,2a)')' Weights of CASSCF configurations ',
     >                 'in VB basis (c_res=c_cas-Svb*c_vb) :'
      write(6,'(2a)')  ' ---------------------------------',
     >                 '------------------------------------'
      if(mod(iciweights,8).gt.3)then
        write(6,'(a)')' Sum of inverse-overlap weights :'
        write(6,form2AD)' c_cas :',sum1,' expected :',one
        write(6,form2AD)' c_res :',sum2,' expected :',one-svb*svb
        write(6,'(a)')' '
      endif
      lenfld=8+iprec
      ix1=max(0,min(3,2*lenfld-16))
      call cblank_cvb(line,240)
      if(npcf.gt.0.or.npcf.eq.-1)then
        line(1:21)='  Conf. =>  Orbitals '
        ibeg=max(14+3*nel,22)
        ibegt=ibeg
        if(mod(iciweights,2).eq.1)then
          line(ibegt+ix1:ibegt+2*lenfld-1)='Chirgwin-Coulson'
          ibegt=ibegt+2*lenfld
        endif
        if(mod(iciweights,4).gt.1)then
          line(ibegt+ix1:ibegt+2*lenfld-1)='Lowdin'
          ibegt=ibegt+2*lenfld
        endif
        if(mod(iciweights,8).gt.3)then
          line(ibegt+ix1:ibegt+2*lenfld-1)='Inverse'
          ibegt=ibegt+2*lenfld
        endif
        write(6,'(a)')line(1:len_trim_cvb(line))
        call cblank_cvb(line,240)
        if(mod(iciweights,2).eq.1)then
          line(ibeg+ix1:ibeg+lenfld-1)='c_cas'
          ibeg=ibeg+lenfld
          line(ibeg+ix1:ibeg+lenfld-1)='c_res'
          ibeg=ibeg+lenfld
        endif
        if(mod(iciweights,4).gt.1)then
          line(ibeg+ix1:ibeg+lenfld-1)='c_cas'
          ibeg=ibeg+lenfld
          line(ibeg+ix1:ibeg+lenfld-1)='c_res'
          ibeg=ibeg+lenfld
        endif
        if(mod(iciweights,8).gt.3)then
          line(ibeg+ix1:ibeg+lenfld-1)='c_cas'
          ibeg=ibeg+lenfld
          line(ibeg+ix1:ibeg+lenfld-1)='c_res'
          ibeg=ibeg+lenfld
        endif
        write(6,'(a)')line(1:len_trim_cvb(line))
      endif
      call cblank_cvb(line,240)

c Alpha weight array:
      do 5100 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nalf,0)
5100  maxgrph(iorb)=min(iorb,nalf)
      call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
c Beta weight array:
      do 5200 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nbet,0)
5200  maxgrph(iorb)=min(iorb,nbet)
      call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)

      nc=0
      call fzero(wghtion1,ionmax-ionmin+1)
      call fzero(wghtion2,ionmax-ionmin+1)
      call fzero(wghtion3,ionmax-ionmin+1)
      call fzero(wghtion4,ionmax-ionmin+1)
      call fzero(wghtion5,ionmax-ionmin+1)
      call fzero(wghtion6,ionmax-ionmin+1)
      do 6000 ion=ionmin,ionmax
      nsing=nel-2*ion
      mrem=norb-ion
      nalfsng=nalf-ion
      nbetsng=nbet-ion
c  Initialise loop for ionic orbitals
      do 6100 iorb=0,norb
      mingion(iorb)=max(iorb-norb+ion,0)
6100  maxgion(iorb)=min(iorb,ion)
      call weight_cvb(xion,mingion,maxgion,ion,norb)
      call imove_cvb(maxgion,nkion,norb+1)
      call occupy_cvb(nkion,norb,locion,lunion)
c  Initialise loop for singly occupied orbitals
      do 6200 iorb=0,mrem
      mingsng(iorb)=max(iorb-mrem+nsing,0)
6200  maxgsng(iorb)=min(iorb,nsing)
      call weight_cvb(xsng,mingsng,maxgsng,nsing,mrem)
      call imove_cvb(maxgsng,nksng,mrem+1)
      call occupy_cvb(nksng,mrem,locsng,lunsng)
c  Initialise loop for singly occupied alpha orbitals
      do 6300 iorb=0,nsing
      mingasg(iorb)=max(iorb-nsing+nalfsng,0)
6300  maxgasg(iorb)=min(iorb,nalfsng)
      call weight_cvb(xasg,mingasg,maxgasg,nalfsng,nsing)
      call imove_cvb(maxgasg,nkasg,nsing+1)
      call occupy_cvb(nkasg,nsing,locasg,lunasg)

c  Loop ionic
      indion=1
7000  continue
c  Loop singly occupied
      indsng=1
7100  continue
      c1=zero
      c2=zero
      c3=zero
      c4=zero
c  Loop singly occupied alpha
      indasg=1
7200  continue

      call izero(iaocc,norb)
      call izero(ibocc,norb)
      do 7300 i=1,ion
      iaocc(locion(i))=1
7300  ibocc(locion(i))=1
      do 7400 ia=1,nalfsng
      iaorb=lunion(locsng(locasg(ia)))
7400  iaocc(iaorb)=1
      do 7500 ib=1,nbetsng
      iborb=lunion(locsng(lunasg(ib)))
7500  ibocc(iborb)=1
      inda=indget_cvb(iaocc,nalf,norb,xalf)
      indb=indget_cvb(ibocc,nbet,norb,xbet)

      indab=(indb-1)*nda+inda
      c1=c1+vec1(indab)
      c2=c2+vec2(indab)
      c3=c3+vec3(indab)
      c4=c4+vec4(indab)

      call loind_cvb(nsing,nalfsng,nkasg,mingasg,maxgasg,
     >                       locasg,lunasg,indasg,xasg,*7200)
      nc=nc+1
      wghtion1(ion)=wghtion1(ion)+c1
      wghtion2(ion)=wghtion2(ion)+c2
      wghtion3(ion)=wghtion3(ion)+c3
      wghtion4(ion)=wghtion4(ion)+c4
      wghtion5(ion)=wghtion5(ion)+gal1(nc)
      wghtion6(ion)=wghtion6(ion)+gal2(nc)
      if(nc.le.npcf.or.npcf.eq.-1)then
        call int2char_cvb(line,nc,7)
        line(9:10)='=>'
        ilin=11
        do 7600 i=1,ion
        call int2char_cvb(line(ilin:ilin+2),locion(i),3)
        ilin=ilin+3
        call int2char_cvb(line(ilin:ilin+2),locion(i),3)
7600    ilin=ilin+3
        do 7700 i=1,nsing
        call int2char_cvb(line(ilin:ilin+2),lunion(locsng(i)),3)
7700    ilin=ilin+3
        ilin=max(ilin,19)
        nprint=0
        if(mod(iciweights,2).eq.1)then
          cprint(nprint+1)=c1
          cprint(nprint+2)=c2
          nprint=nprint+2
        endif
        if(mod(iciweights,4).gt.1)then
          cprint(nprint+1)=c3
          cprint(nprint+2)=c4
          nprint=nprint+2
        endif
        if(mod(iciweights,8).gt.3)then
          cprint(nprint+1)=gal1(nc)
          cprint(nprint+2)=gal2(nc)
          nprint=nprint+2
        endif
        write(6,formAD)line(1:ilin),(cprint(mp),mp=1,nprint)
      endif
      call loind_cvb(mrem,nsing,nksng,mingsng,maxgsng,
     >                       locsng,lunsng,indsng,xsng,*7100)
      call loind_cvb(norb,ion,nkion,mingion,maxgion,
     >                       locion,lunion,indion,xion,*7000)
6000  continue

      call cblank_cvb(line,240)
      line(1:21)=' Accumulated weights:'
      ibeg=22
      if(mod(iciweights,2).eq.1)then
        line(ibeg+ix1:ibeg+2*lenfld-1)='Chirgwin-Coulson'
        ibeg=ibeg+2*lenfld
      endif
      if(mod(iciweights,4).gt.1)then
        line(ibeg+ix1:ibeg+2*lenfld-1)='Lowdin'
        ibeg=ibeg+2*lenfld
      endif
      if(mod(iciweights,8).gt.3)then
        line(ibeg+ix1:ibeg+2*lenfld-1)='Inverse'
        ibeg=ibeg+2*lenfld
      endif
      write(6,'(/,a)')line(1:len_trim_cvb(line))
      call cblank_cvb(line,240)
      ibeg=22
      if(mod(iciweights,2).eq.1)then
        line(ibeg+ix1:ibeg+lenfld-1)='c_cas'
        ibeg=ibeg+lenfld
        line(ibeg+ix1:ibeg+lenfld-1)='c_res'
        ibeg=ibeg+lenfld
      endif
      if(mod(iciweights,4).gt.1)then
        line(ibeg+ix1:ibeg+lenfld-1)='c_cas'
        ibeg=ibeg+lenfld
        line(ibeg+ix1:ibeg+lenfld-1)='c_res'
        ibeg=ibeg+lenfld
      endif
      if(mod(iciweights,8).gt.3)then
        line(ibeg+ix1:ibeg+lenfld-1)='c_cas'
        ibeg=ibeg+lenfld
        line(ibeg+ix1:ibeg+lenfld-1)='c_res'
        ibeg=ibeg+lenfld
      endif
      write(6,'(a)')line(1:len_trim_cvb(line))
      call cblank_cvb(line,240)
      total1=zero
      total2=zero
      total3=zero
      total4=zero
      total5=zero
      total6=zero
      do 8000 ion=ionmin,ionmax
      total1=total1+wghtion1(ion)
      total2=total2+wghtion2(ion)
      total3=total3+wghtion3(ion)
      total4=total4+wghtion4(ion)
      total5=total5+wghtion5(ion)
      total6=total6+wghtion6(ion)
      line(1:19)=' For ionicity    : '
      call int2char_cvb(line(14:16),ion,3)
      nprint=0
      if(mod(iciweights,2).eq.1)then
        cprint(nprint+1)=wghtion1(ion)
        cprint(nprint+2)=wghtion2(ion)
        nprint=nprint+2
      endif
      if(mod(iciweights,4).gt.1)then
        cprint(nprint+1)=wghtion3(ion)
        cprint(nprint+2)=wghtion4(ion)
        nprint=nprint+2
      endif
      if(mod(iciweights,8).gt.3)then
        cprint(nprint+1)=wghtion5(ion)
        cprint(nprint+2)=wghtion6(ion)
        nprint=nprint+2
      endif
      write(6,formAD)line(1:19),(cprint(mp),mp=1,nprint)
8000  continue
      nprint=0
      if(mod(iciweights,2).eq.1)then
        cprint(nprint+1)=total1
        cprint(nprint+2)=total2
        nprint=nprint+2
      endif
      if(mod(iciweights,4).gt.1)then
        cprint(nprint+1)=total3
        cprint(nprint+2)=total4
        nprint=nprint+2
      endif
      if(mod(iciweights,8).gt.3)then
        cprint(nprint+1)=total5
        cprint(nprint+2)=total6
        nprint=nprint+2
      endif
      write(6,formAD)' Total all       : ',(cprint(mp),mp=1,nprint)
      write(6,'(a)')' '
      return
      end
