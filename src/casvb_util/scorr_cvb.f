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
      subroutine scorr_cvb(cvbdet,dvbdet,evbdet)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "malloc_cvb.fh"
      dimension cvbdet(ndetvb),dvbdet(ndetvb),evbdet(ndetvb)

      k1 = mstackr_cvb(norb*norb)
      k2 = mstackr_cvb(ndetvb)
      k3 = mstacki_cvb(norb)
      call scorr2_cvb(cvbdet,dvbdet,evbdet,
     >      w(k1),w(k2),iw(k3))
      call mfreer_cvb(k1)
      return
      end
