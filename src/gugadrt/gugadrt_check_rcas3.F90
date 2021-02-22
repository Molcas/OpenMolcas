!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      subroutine gugadrt_check_rcas3(jk,ind,inb,ndj,locu)
#include "gendrt.fh"
#include "casrst_drt.fh"
      dimension ind(8,max_node),lsym(8),iexcit(ndj),locu(8,ndj)
      inb=0
      nsumel=0
      do i=1,8
        lsym(i)=ind(i,jk)
        nsumel=nsumel+lsym(i)
      enddo
      do i=1,ndj
        iexcit(i)=0
        do m=1,8
          iex=lsym(m)-locu(m,i)
         if(iex.gt.0) then
           iexcit(i)=iexcit(i)+iex
          endif
        enddo
      enddo
      inb=iexcit(1)
      do i=2,ndj
         inb=min(inb,iexcit(i))
      enddo
      inb=inb+ja(jk)*2+jb(jk)
      return
      end
