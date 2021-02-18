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
      subroutine ao2mo_cvb(orbsao,orbs,norb1)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "mo_cvb.fh"
      dimension orbsao(nbas_mo,norb1),orbs(norb,norb1)

      if(norb1.eq.0)return
      i1=mstackr_cvb(nbas_mo*norb)
      call getmo_cvb(w(i1),3,0)
      call mxattb_cvb(w(i1),orbsao,norb,nbas_mo,norb1,orbs)
      call mfreer_cvb(i1)
      return
      end
      subroutine mo2ao_cvb(orbs,orbsao,norb1)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "mo_cvb.fh"
      dimension orbsao(nbas_mo,norb1),orbs(norb,norb1)
      if(norb1.eq.0)return
      i1=mstackr_cvb(nbas_mo*norb)
      call getmo_cvb(w(i1),2,0)
      call mxatb_cvb(w(i1),orbs,nbas_mo,norb,norb1,orbsao)
      call mfreer_cvb(i1)
      return
      end
      subroutine lmo2ao_cvb(orbs,orbsao,norb1)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "mo_cvb.fh"
      dimension orbsao(nbas_mo,norb1),orbs(norb,norb1)
      if(norb1.eq.0)return
      i1=mstackr_cvb(nbas_mo*norb)
      call getmo_cvb(w(i1),2,1)
      call mxatb_cvb(w(i1),orbs,nbas_mo,norb,norb1,orbsao)
      call mfreer_cvb(i1)
      return
      end
