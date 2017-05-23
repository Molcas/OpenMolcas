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
      subroutine casinfoprint_cvb()
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"

      if(ip(1).ge.0.and.(.not.up2date_cvb('CASPRINT')))then
        write(6,'(/,a,i4)')' Number of active electrons :',nel
        write(6,'(a,i4)')  ' Number of active orbitals  :',norb
        write(6,'(a,f4.1)')' Total spin                 :',
     >    DBLE(nalf-nbet)/two
        if(nsym.eq.1)then
          write(6,'(a,i4)')' State symmetry             :',isym
        else
          iqisym=mstacki_cvb(nsym)
          incr=0
          do i=1,mxirrep
          if(isymv(i).eq.1)then
            incr=incr+1
            iw(incr+iqisym-1)=i
          endif
          enddo
          write(6,'(a,i4,7i3)')' State symmetries           :',
     >      (iw(ii+iqisym-1),ii=1,nsym)
          call mfreei_cvb(iqisym)
        endif
        write(6,'(/,a,100i3)')' Symmetries of active MOs   : ',
     >    (ityp(ii),ii=1,norb)
        call make_cvb('CASPRINT')
      endif
      return
      end
