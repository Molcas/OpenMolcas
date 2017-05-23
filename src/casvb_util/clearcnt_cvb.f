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
      subroutine clearcnt_cvb(icode)
c  ICODE=1 : Orbitals changed
c  ICODE=2 : CI coefficients changed
c  ICODE=3 : Everything changed
      implicit real*8 (a-h,o-z)
      logical initialize
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
      dimension iunset(mxciobj,2)
      data initialize/.true./
      save iunset,initialize
      if (initialize) then
       iunset(1,1)=0
       iunset(1,2)=0
       do 50 i=2,mxciobj
       iunset(i,1)=1
50     iunset(i,2)=1
       initialize=.false.
      endif

      if(icode.eq.3)then
        do 100 i=1,mxciobj
100     icnt_ci(i)=0
      else
        ipow1=2
        ipow2=1
        do 200 ichg=1,2
        if(mod(icode,ipow1).ge.ipow2)then
          do 300 i=1,mxciobj
300       if(iunset(i,ichg).eq.1)icnt_ci(i)=0
        endif
        ipow1=2*ipow1
200     ipow2=2*ipow2
      endif
      return
      end
