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
      subroutine setipermzeta_cvb(ipermzeta,
     >  orbs,symelm,izeta,
     >  orbinv,owrk,owrk2)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension ipermzeta(norb,nzeta)
      dimension orbs(norb,norb)
      dimension symelm(norb*norb,nsyme),izeta(nsyme)
      dimension orbinv(norb,norb),owrk(norb,norb),owrk2(norb,norb)
      save thresh
      data thresh/1.d-8/

      if(nzeta.gt.0)then
        call fmove_cvb(orbs,orbinv,norb*norb)
        call mxinv_cvb(orbinv,norb)
      endif

      izeta1=0
      do 100 isyme=1,nsyme
      if(izeta(isyme).ne.0)then
      izeta1=izeta1+1
c  Determine orbital permutation for sym. operation ISYME :
        call mxatb_cvb(symelm(1,isyme),orbs,norb,norb,norb,owrk2)
        call mxatb_cvb(orbinv,owrk2,norb,norb,norb,owrk)
        do 200 iorb=1,norb
        do 201 jorb=1,norb
        if(abs(abs(owrk(jorb,iorb))-one).lt.thresh)then
          ipermzeta(iorb,izeta1)=nint(owrk(jorb,iorb))*jorb
        elseif(abs(owrk(jorb,iorb)).gt.thresh)then
          write(6,*)' Fatal error! Symmetry operation ',tags(isyme),
     >      ' does not permute the VB orbitals!'
          call mxprint_cvb(owrk,norb,norb,0)
          call abend_cvb()
        endif
201     continue
200     continue
      endif
100   continue
      return
      end
