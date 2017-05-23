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
      subroutine mmstringen_cvb(norb,nel,locc,lunocc,nstring,
     >  nkmin,nkmax)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension locc(*),lunocc(*)
      dimension nkmin(0:norb),nkmax(0:norb)

      i_nk    = mstacki_cvb(norb+1)
c Spin string loop initialization (use xdet as graph storage) :
      call imove_cvb(nkmax,iw(i_nk),norb+1)
c  Spin string loop starts here :
      index=0
100   index=index+1
      i_locc=(index-1)*nel+1
      i_lunocc=(index-1)*(norb-nel)+1
      call occupy_cvb(iw(i_nk),norb,locc(i_locc),lunocc(i_lunocc))
      call loop_cvb(norb,iw(i_nk),nkmin,nkmax,*100)
      call mfreei_cvb(i_nk)
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(nstring)
      end
      real*8 function party_cvb(iperm,n)
c  Returns parity of permutation
      implicit real*8(a-h,o-z)
      dimension iperm(n)
#include "malloc_cvb.fh"

      i1 = mstacki_cvb(n)
      call imove_cvb(iperm,iw(i1),n)
      call party2_cvb(iw(i1),n,party)
      call mfreei_cvb(i1)
      party_cvb=party
      return
      end
