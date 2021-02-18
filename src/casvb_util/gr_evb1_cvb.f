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
      subroutine gr_evb1_cvb(civbh,civbs,civb,dvbdet,
     >   grad,grad1,grad2,gradx,vec1)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "fx_cvb.fh"
      dimension civbh(ndet),civbs(ndet),civb(ndet)
      dimension dvbdet(ndetvb)
      dimension grad(npr),grad1(npr),grad2(npr),gradx(norb,norb)

c  VEC1 dimension is MAX(NPRORB,NDETVB)
      dimension vec1(*)

      f1=one/ovraa
      f2=-two*f1*f1*ww
      f1=two*f1
      f3=-f1*f1
      f4=-f1*f3*ww

      call fzero(gradx,norb*norb)
      call onedens_cvb(civb,civbs,gradx,.true.,1)

      call mkgrd_cvb(civb,civbs,grad1,dvbdet,npr,.true.)
      call mkgrd_cvb(civb,civbh,grad2,dvbdet,npr,.true.)

c  Use VEC1 as work:
      do 100 i=1,npr
      vec1(i)=f1*grad2(i)+f2*grad1(i)
100   continue
      call prgrad_cvb(vec1,npr)

      call make_cvb('ORBFREE')
      call make_cvb('CIFREE')
      call all2free_cvb(vec1,grad,1)
      return
      end
