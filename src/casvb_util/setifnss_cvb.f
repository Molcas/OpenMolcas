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
      subroutine setifnss_cvb(ifnss1,ifnss2,ndetvbs)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


      dimension ifnss1(0:nel,0:nel),ifnss2(0:nel,0:nel)
      dimension ndetvbs(0:nel,0:nel)

      call izero(ifnss1,(nel+1)*(nel+1))
      call izero(ifnss2,(nel+1)*(nel+1))
      call izero(ndetvbs,(nel+1)*(nel+1))

      do 100 n=0,nel
      do 101 nalfa=(n+1)/2,n
      nbeta=n-nalfa
      if(nbeta.gt.nalfa)goto 101
      call icomb_cvb(n,nbeta,iretval1)
      call icomb_cvb(n,nbeta-1,iretval2)
      ifnss1(n,nalfa-nbeta)=iretval1-iretval2
      call icomb_cvb(n,nalfa,ifnss2(n,nalfa-nbeta))
      if(nalfa.eq.nbeta)ifnss2(n,nalfa-nbeta)=
     >  (ifnss2(n,nalfa-nbeta)+1)/2
      call icomb_cvb(n,nalfa,ndetvbs(n,nalfa))
101   continue
100   continue
      return
      end
