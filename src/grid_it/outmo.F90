!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      subroutine outmo(imo,ipower,cmo,clincomb,cout,nbas,nmo)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
      implicit real*8 (a-h,o-z)
      dimension cmo(nbas,nmo),clincomb(nmo),cout(nbas)
#include "WrkSpc.fh"

      if(imo.ne.0)then
         call fmove(cmo(1,imo),cout,nbas)
         call power(cout,nbas,ipower)
         return
      else
        call fzero(cout,nbas)
        Call GetMem('TmpMo','ALLO','REAL',ipTmpMo,nbas)
        do 100 i=1,nmo
        if(clincomb(i).ne.0d0)then
          call fmove(cmo(1,i),Work(ipTmpMo),nbas)
          call power(Work(ipTmpMO),nbas,ipower)
          call daxpy_(nbas,clincomb(i),Work(ipTmpMo),1,cout,1)
        endif
100     continue
        Call GetMem('TmpMo','FREE','REAL',ipTmpMo,nbas)
      endif
      return
      end
