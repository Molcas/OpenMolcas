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
      subroutine orthcon_cvb(ipairs,ipair,igroups,ngroup,iorthlst,
     >  mxortl,mxpair)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      parameter (nstrin=7,ncmp=4,mxgroup=40)
      character*8 string(nstrin)
      character*3 glabel(mxgroup)
      dimension ipairs(2,mxpair),ipair(mxorb,mxorb)
      dimension igroups(mxorb,mxgroup),ngroup(mxgroup)
      dimension iorthlst(mxortl)
      save string
      data string/'GROUP   ','ORTH    ','PAIRS   ','STRONG  ',
     >            'FULL    ','END     ','ENDORTHC'/

      call izero(ipair,mxorb*mxorb)
      ngrp=0
2000  call fstring_cvb(string,nstrin,istr,ncmp,2)
      if(istr.eq.1)then
c 'GROUP'
      ngrp=ngrp+1
      if(ngrp.gt.mxgroup)then
        write(6,*)' Too many GROUP keywords in input!',mxgroup
        call abend_cvb()
      endif
      glabel(ngrp)=' '
      call string_cvb(glabel(ngrp),1,nread,1)
      if(glabel(ngrp)(1:1).lt.'A'.or.glabel(ngrp)(1:1).gt.'Z')then
        write(6,*)' Group label must begin with a character A-Z: ',
     >    glabel(ngrp)
        call abend_cvb()
      endif
      call int_cvb(igroups(1,ngrp),mxorb,ngroup(ngrp),0)
      if(ngroup(ngrp).eq.-1)then
        write(6,*)' Too many elements for group ',glabel(ngrp)
        call abend_cvb()
      endif
      do 100 i=1,ngroup(ngrp)
      if(igroups(i,ngrp).lt.1.or.igroups(i,ngrp).gt.mxorb)then
        write(6,*)' Illegal orbital number in group ',glabel(ngrp),
     >    ' :',igroups(i,ngrp)
        call abend_cvb()
      endif
100   continue
      do 150 i=1,ngrp-1
      if(glabel(ngrp).eq.glabel(i))then
        write(6,*)' Repeated label in GROUP keywords : ',glabel(ngrp)
        call abend_cvb()
      endif
150   continue
      elseif(istr.eq.2)then
c 'ORTH'
      nsp=0
175   continue
      call int_cvb(iorthlst(1+nsp),mxortl-nsp,nread,0)
      nsp=nsp+nread
      if(mxortl-nsp.gt.0)then
        call fstring_cvb(glabel,ngrp,igrp,3,0)
        if(igrp.gt.0)then
          nsp=nsp+1
          iorthlst(nsp)=-igrp
          if(mxortl-nsp.gt.0)goto 175
        endif
      endif
      do 500 isp=1,nsp
      do 500 jsp=isp+1,nsp
      ior=iorthlst(isp)
      jor=iorthlst(jsp)
      if(ior.gt.0.and.jor.gt.0)then
        ipair(ior,jor)=1
      elseif(ior.gt.0.and.jor.lt.0)then
        jor2=-jor
        do 600 jo=1,ngroup(jor2)
600     ipair(ior,igroups(jo,jor2))=1
      elseif(ior.lt.0.and.jor.gt.0)then
        ior2=-ior
        do 700 io=1,ngroup(ior2)
700     ipair(jor,igroups(io,ior2))=1
      elseif(ior.lt.0.and.jor.lt.0)then
        ior2=-ior
        jor2=-jor
        do 800 io=1,ngroup(ior2)
        do 800 jo=1,ngroup(jor2)
800     ipair(igroups(io,ior2),igroups(jo,jor2))=1
      endif
500   continue
      elseif(istr.eq.3)then
c 'PAIRS'
      nsp=0
975   continue
      call int_cvb(ipairs(1+nsp,1),2*mxpair-nsp,nread,0)
      nsp=nsp+nread
      if(2*mxpair-nsp.gt.0)then
        call fstring_cvb(glabel,ngrp,igrp,3,0)
        if(igrp.gt.0)then
          nsp=nsp+1
          ipairs(nsp,1)=-igrp
          if(mxortl-nsp.gt.0)goto 975
        endif
      endif
      if(mod(nsp,2).eq.1)then
        write(6,*)' Odd number of orthogonalization numbers in PAIRS!'
        call abend_cvb()
      endif
      npairs=nsp/2
      do 1300 ipar=1,npairs
      ior=ipairs(1,ipar)
      jor=ipairs(2,ipar)
      if(ior.gt.0.and.jor.gt.0)then
        ipair(ior,jor)=1
      elseif(ior.gt.0.and.jor.lt.0)then
        jor2=-jor
        do 1400 jo=1,ngroup(jor2)
1400    ipair(ior,igroups(jo,jor2))=1
      elseif(ior.lt.0.and.jor.gt.0)then
        ior2=-ior
        do 1500 io=1,ngroup(ior2)
1500    ipair(jor,igroups(io,ior2))=1
      elseif(ior.lt.0.and.jor.lt.0)then
        ior2=-ior
        jor2=-jor
        do 1600 io=1,ngroup(ior2)
        do 1600 jo=1,ngroup(jor2)
1600    ipair(igroups(io,ior2),igroups(jo,jor2))=1
      endif
1300  continue
      elseif(istr.eq.4)then
c 'STRONG'
      do 1700 i=1,mxorb
      do 1700 j=i+1,mxorb
1700  if(.not.(mod(i,2).eq.1.and.j.eq.i+1))ipair(i,j)=1
      elseif(istr.eq.5)then
c 'FULL'
      do 1800 i=1,mxorb
      do 1800 j=i+1,mxorb
1800  ipair(i,j)=1
      endif
c 'END' , 'ENDORTHC' or unrecognized keyword -- end of ORTHCON input :
      if(.not.(istr.eq.6.or.istr.eq.7.or.istr.eq.0))goto 2000
      call izero(ipairs,2*mxpair)
      nort=0
      do 1900 i=1,mxorb
      do 1900 j=i+1,mxorb
      if(ipair(i,j).eq.1.or.ipair(j,i).eq.1)then
        nort=nort+1
        ipairs(1,nort)=i
        ipairs(2,nort)=j
      endif
1900  continue
      return
      end
