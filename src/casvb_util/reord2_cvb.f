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
      subroutine reord2_cvb(cfrom,cto,imode)
      implicit real*8(a-h,o-z)
c Front-end routine for molcas reord2, transforms
c from SGA CSFs to split-graph-GUGA CSFs.
      dimension cfrom(*),cto(*)
#include "WrkSpc.fh"
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "csfbas.fh"

c NAC      rasscf.fh
c NACTEL   general.fh
c STSYM    general.fh
c IPR      rasscf.fh
c KICONF   csfbas.fh
c KCFTP    csfbas.fh
      call getmem('kcnf','allo','inte',ivkcnf,nactel)
      call reord2(nac,nactel,stsym,imode,
     >  iwork(kiconf(1)),iwork(kcftp),
     >  cfrom,cto,iWork(ivkcnf))
      call getmem('kcnf','free','inte',ivkcnf,nactel)
      return
      end
