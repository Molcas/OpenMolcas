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
      subroutine mkvbinfo_cvb()
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"

      if(nfrag.gt.1)then
        call dpgendet_cvb()
      else
        call vbgendet_cvb(
     >    iw(ll(11)),iw(ll(12)),iw(ll(13)),iw(ll(14)),
     >    iw(ll(15)),iw(ll(17)),
     >    nconf,nconfion_fr(0,1),
     >    nda,ndb,ndetvb,nel,
     >    noe,nalf,nbet,norb)
      endif
      return
      end
