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
      subroutine mxprint_cvb(a,nrow,ncol,itype)
c Prints matrix A, stored according to ITYPE
      implicit real*8 (a-h,o-z)
#include "print_cvb.fh"
#include "formats_cvb.fh"
      parameter (mxbuf=8)
      dimension buffer(mxbuf),ibuf(mxbuf),a(*)

      nbuf=min((iwidth-4)/(iprec+4),mxbuf)
      if(nbuf.eq.7)nbuf=6
      iform=0
      jin=1
100   jend=jin+nbuf-1
      if(ncol.le.nbuf) jend=ncol
      if(jend.gt.ncol+nbuf-1) return
      jend=min(ncol,jend)
      k=0
      do 200 j=jin,jend
      k=k+1
      ibuf(k)=j
200   continue
      if(iform.eq.0)then
        write(6,formMXP1)(ibuf(i),i=1,jend-jin+1)
      else
        write(6,formMXP2)(ibuf(i),i=1,jend-jin+1)
      endif
      do 300 i=1,nrow
      k=0
      do 400 j=jin,jend
      k=k+1
      if(itype.eq.0)then
        ind=(j-1)*nrow+i
      elseif(itype.eq.1)then
        if(i.ge.j)then
          ind=i*(i-1)/2+j
        else
          ind=j*(j-1)/2+i
        endif
      else
        ind=(i-1)*ncol+j
      endif
      buffer(k)=a(ind)
400   continue
      if(iform.eq.0)then
        write(6,formMXP3)i,(buffer(ii),ii=1,jend-jin+1)
      else
        write(6,formMXP4)i,(buffer(ii),ii=1,jend-jin+1)
      endif
300   continue
      jin=jend+1
      if(ncol.gt.nbuf)goto 100
      return
      end
      subroutine mxprintd_cvb(a,nrow,ncol,itype)
      implicit real*8 (a-h,o-z)
#include "print_cvb.fh"
#include "formats_cvb.fh"
      parameter (mxbuf=8)
      dimension buffer(mxbuf),ibuf(mxbuf),a(*)
      nbuf=min((iwidth-4)/(iprec+8),mxbuf)
      if(nbuf.eq.7)nbuf=6
      iform=1
      jin=1
100   jend=jin+nbuf-1
      if(ncol.le.nbuf) jend=ncol
      if(jend.gt.ncol+nbuf-1) return
      jend=min(ncol,jend)
      k=0
      do 200 j=jin,jend
      k=k+1
      ibuf(k)=j
200   continue
      if(iform.eq.0)then
        write(6,formMXP1)(ibuf(i),i=1,jend-jin+1)
      else
        write(6,formMXP2)(ibuf(i),i=1,jend-jin+1)
      endif
      do 300 i=1,nrow
      k=0
      do 400 j=jin,jend
      k=k+1
      if(itype.eq.0)then
        ind=(j-1)*nrow+i
      elseif(itype.eq.1)then
        if(i.ge.j)then
          ind=i*(i-1)/2+j
        else
          ind=j*(j-1)/2+i
        endif
      else
        ind=(i-1)*ncol+j
      endif
      buffer(k)=a(ind)
400   continue
      if(iform.eq.0)then
        write(6,formMXP3)i,(buffer(ii),ii=1,jend-jin+1)
      else
        write(6,formMXP4)i,(buffer(ii),ii=1,jend-jin+1)
      endif
300   continue
      jin=jend+1
      if(ncol.gt.nbuf)goto 100
      return
      end
