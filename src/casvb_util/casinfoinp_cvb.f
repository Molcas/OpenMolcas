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
      subroutine casinfoinp_cvb()
      implicit real*8 (a-h,o-z)
#include "casinfo_cvb.fh"
      parameter (nstrin=6,ncmp=4)
      character*8 string
      dimension string(nstrin)
      save string
      data string/'RAS2    ','INACTIVE','FROZEN  ','NACTEL  ',
     >  'SPIN    ','SYMMETRY'/
      dimension nactel(3)
      save one
      data one/1d0/

1000  call fstring_cvb(string,nstrin,istr,ncmp,2)

      if(istr.eq.1)then
c 'RAS2'
        call izero(iorocc_d,mxirrep_ci)
        call int_cvb(iorocc_d,mxirrep_ci,nread,1)
      elseif(istr.eq.2)then
c 'INACTIVE'
        call izero(iorclos_d,mxirrep_ci)
        call int_cvb(iorclos_d,mxirrep_ci,nread,1)
      elseif(istr.eq.3)then
c 'FROZEN'
        call izero(iorcore_d,mxirrep_ci)
        call int_cvb(iorcore_d,mxirrep_ci,nread,1)
      endif
      if(istr.eq.4.or.istr.eq.5.or.istr.eq.6)then
c 'NACTEL' or 'SPIN' or 'SYMMETRY'
        if(nstsym_d.eq.0)then
          call izero(istnel_d,mxstsy_ci)
          call izero(istsy_d,mxstsy_ci)
          call izero(istms2_d,mxstsy_ci)
          call izero(nstats_d,mxstsy_ci)
          call fzero(weight_d,mxstt_ci*mxstsy_ci)
          nstsym_d=1
          istsy_d(nstsym_d)=1
          nstats_d(nstsym_d)=1
          weight_d(1,nstsym_d)=one
        endif
      endif
      if(istr.eq.4)then
c 'NACTEL'
        call izero(nactel,3)
        call int_cvb(nactel,3,nread,1)
        if(nactel(2).ne.0.or.nactel(3).ne.0)then
          write(6,*)' Illegal NACTEL read :',nactel
          write(6,*)' Use CASVB only with CASSCF wavefunctions!'
          call abend_cvb()
        endif
        istnel_d(nstsym_d)=nactel(1)
      elseif(istr.eq.5)then
c 'SPIN'
        call int_cvb(istms2_d(nstsym_d),1,nread,1)
        istms2_d(nstsym_d)=istms2_d(nstsym_d)-1
      elseif(istr.eq.6)then
c 'SYMMETRY'
        call int_cvb(istsy_d(nstsym_d),1,nread,1)
      endif
c Unrecognized keyword -- end of input :
      if(istr.ne.0)goto 1000
      return
      end
