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
      subroutine vbgendet_cvb(
     >  iapr,ixapr,ibpr,ixbpr,
     >  iconfs,idetvb,
     >  nconf,nconfion,
     >  nda,ndb,ndetvb,nel,
     >  noe,nalf,nbet,norb)
      implicit real*8 (a-h,o-z)
      dimension iapr(ndetvb),ixapr(nda+1),ibpr(ndetvb),ixbpr(ndb+1)
      dimension iconfs(noe,nconf),idetvb(ndetvb)
      dimension nconfion(0:nel)
#include "malloc_cvb.fh"

      k1= mstacki_cvb(ndetvb)
      k2= mstacki_cvb(ndetvb)
      k3= mstacki_cvb(ndetvb)
      k4= mstacki_cvb(ndetvb)
      call vbgendet2_cvb(
     >  iapr,ixapr,ibpr,ixbpr,
     >  iconfs,idetvb,
     >  nconf,nconfion,
     >  nda,ndb,ndetvb,nel,
     >  noe,nalf,nbet,norb,
     >  iw(k1),iw(k2),iw(k3),iw(k4))
      call mfreei_cvb(k1)
      return
      end
