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
      subroutine symelminp_cvb(ip_symelm,nsyme,tags,izeta,
     >  mxirrep,mxorb,mxsyme,ityp)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "malloc_cvb.fh"
      parameter (nsymelm=5,nsign=2,ncmp=4)
      character*8 symelm(nsymelm),sign(nsign)
      character*3 tags(mxsyme)
      dimension izeta(*)
      dimension ityp(mxorb)
      dimension iaux(1),daux(1)
      save symelm,sign
      data symelm/'IRREPS  ','COEFFS  ','TRANS   ','END     ',
     >            'ENDSYMEL'/
      data sign/'+       ','-       '/
      save zero,one
      data zero/0d0/,one/1d0/

      nsyme=nsyme+1
      if(nsyme.gt.mxsyme)then
        write(6,*)' Too many symmetry elements found :',nsyme,mxsyme
        call abend_cvb()
      endif
      tags(nsyme)=' '
      call string_cvb(tags(nsyme),1,nread,1)
      call fstring_cvb(sign,nsign,isign,ncmp,1)
      if(isign.eq.1)then
        izeta(nsyme)=1
      elseif(isign.eq.2)then
        izeta(nsyme)=-1
      else
        izeta(nsyme)=0
      endif
      call mreallocr_cvb(ip_symelm,mxorb*mxorb*nsyme)
      ishft=mxorb*mxorb*(nsyme-1)
      call mxunit_cvb(w(ishft+ip_symelm),mxorb)

1000  call fstring_cvb(symelm,nsymelm,istr2,ncmp,2)
      if(istr2.eq.1)then
c    'IRREPS'
        do 1100 i=1,mxirrep
        iaux=0
        call int_cvb(iaux,1,nread,0)
        irrep=iaux(1)
        if(irrep.ne.0)then
          do 1200 iorb=1,mxorb
1200      if(irrep.eq.ityp(iorb))
     >      w(iorb+(iorb-1)*mxorb+ishft+ip_symelm-1)=-one
        endif
1100    continue
      elseif(istr2.eq.2)then
c    'COEFFS'
        do 1300 i=1,mxorb
        iaux=0
        call int_cvb(iaux,1,nread,0)
        iorb=iaux(1)
        if(iorb.ne.0)then
          w(iorb+(iorb-1)*mxorb+ishft+ip_symelm-1)=-one
        else
          goto 1301
        endif
1300    continue
1301    continue
      elseif(istr2.eq.3)then
c    'TRANS'
        iaux=0
        call int_cvb(iaux,1,nread,0)
        idim=iaux(1)
        if(idim.lt.1.or.idim.gt.mxorb)then
          write(6,*)' Illegal dimension in TRANS:',idim,mxorb
          call abend_cvb()
        endif
        itmp = mstacki_cvb(idim)
        do 1400 i=1,idim
        call int_cvb(iaux,1,nread,0)
        iorb=iaux(1)
        if(iorb.lt.1.or.iorb.gt.mxorb)then
          write(6,*)' Illegal orbital number in TRANS:',iorb
          call abend_cvb()
        endif
1400    iw(i+itmp-1)=iorb
        do 1500 ior=1,idim
        iorb=iw(ior+itmp-1)
        do 1500 jor=1,idim
        jorb=iw(jor+itmp-1)
        daux=zero
        call real_cvb(daux(1),1,nread,0)
1500    w(iorb+(jorb-1)*mxorb+ishft+ip_symelm-1)=daux(1)
        call mfreei_cvb(itmp)
      endif
c    'END' , 'ENDSYMEL' or unrecognized keyword -- end SYMELM input :
      if(.not.(istr2.eq.4.or.istr2.eq.5.or.istr2.eq.0))goto 1000
      if(.not.mxorth_cvb(w(ishft+ip_symelm),mxorb))then
        write(6,*)' Symmetry element ',tags(nsyme),' not orthogonal!'
        write(6,*)' Check usage of TRANS keyword.'
        call abend_cvb()
      endif
      return
      end
