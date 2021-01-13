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
      subroutine initopt_cvb(icrit,lfxvb,nfxvb,iorts,nort,norb)
      implicit real*8(a-h,o-z)
#include "inpmod_cvb.fh"
#include "seth_cvb.fh"
#include "initopt_cvb.fh"
      dimension iorts(2,*)

c  IOPTCODE    +1  = REPORT
c              +2  = OPTIM
c              +4  = Svb
c              +8  = freeze structure coefficients
c              +16 = strong-orthogonality constraints

      if(ioptim.eq.0) return

      if(mod(ioptcode(ioptim),4).ge.2)then
         call setifinish_cvb(1)
      elseif(mod(ioptcode(ioptim),2).eq.1)then
         call setifinish_cvb(3)
      endif
      if(mod(ioptcode(ioptim),8).ge.4)icrit=1
      if(mod(ioptcode(ioptim),16).ge.8)then
        lfxvb=1
        nfxvb=0
      endif
      if(mod(ioptcode(ioptim),32).ge.16)then
        nort=0
        do 100 iorb=1,norb
        do 101 jorb=iorb+1,norb
        if(.not.(jorb.eq.iorb+1.and.mod(iorb,2).eq.1))then
          nort=nort+1
          iorts(1,nort)=iorb
          iorts(2,nort)=jorb
        endif
101     continue
100     continue
      endif
      return
      end
