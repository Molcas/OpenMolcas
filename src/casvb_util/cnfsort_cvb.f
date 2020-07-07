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
      subroutine cnfsort_cvb(iconfs,nconf1,nel1,ioncty,iconfs2)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


      dimension iconfs(noe,nconf1),ioncty(nconf1),iconfs2(noe,nconf1)

      mnion1=nel1/2
      mxion1=0
      do 100 iconf=1,nconf1
      ion=0
      do 200 iorb=1,norb
      if(iconfs(iorb,iconf).eq.2)ion=ion+1
200   continue
      ioncty(iconf)=ion
      if(ion.lt.mnion1)mnion1=ion
      if(ion.gt.mxion1)mxion1=ion
100   continue
      jconf=0
      do 300 ion=mnion1,mxion1
      do 301 iconf=1,nconf1
      if(ioncty(iconf).eq.ion)then
        jconf=jconf+1
        call imove_cvb(iconfs(1,iconf),iconfs2(1,jconf),noe)
      endif
301   continue
300   continue
      if(jconf.ne.nconf1)then
        write(6,*)' Error in cnfsort - jconf not same as nconf1 :',
     >    jconf,nconf1
        call abend_cvb()
      endif
        call imove_cvb(iconfs2,iconfs,noe*nconf1)
      return
      end
