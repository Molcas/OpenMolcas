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
      subroutine biks_cvb(aikcof,bikcof,ikcoff,
     >  nel,kbasis,share,iprint)
      implicit real*8 (a-h,o-z)
      logical share
      character*10 basis(7)
      dimension aikcof(*),bikcof(*),ikcoff(0:nel,0:nel,0:nel)
      save basis
      data basis/'Kotani','Serber','Rumer','Rumer (LT)',
     >  'projected','Determ','Determ'/

      aikcof(1)=DBLE(kbasis)
      bikcof(1)=DBLE(kbasis)
      if(kbasis.eq.6)return

      if(iprint.ge.1)write(6,6100)
     >  basis(kbasis)(1:len_trim_cvb(basis(kbasis)))

      do 100 nel1=0,nel
      do 100 nalf1=0,nel
      do 100 i2s1=0,nel
      if(ikcoff(nel1,nalf1,i2s1).ne.-1)then
        ifns=ifns_cvb(nel1,(nel1+i2s1)/2,kbasis)
        ndet=ndet_cvb(nel1,nalf1)
        call bikset_cvb(
     >    aikcof(2+ikcoff(nel1,nalf1,i2s1)),
     >    bikcof(2+ikcoff(nel1,nalf1,i2s1)),
     >    nel1,nalf1,i2s1,ndet,ifns,kbasis,share,iprint)
      endif
100   continue
      return
6100  format(/,' Generate ',a,' spin functions.')
      end
