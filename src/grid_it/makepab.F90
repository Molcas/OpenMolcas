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
      subroutine MakePab(cmo,occ,cout,nCMO,nMOs,nIrrep,nBas)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
      implicit real*8 (a-h,o-z)
      dimension cmo(nCMO),occ(nMOs),cout(nMOs), nBas(0:7)
#include "WrkSpc.fh"
        Call GetMem('List','LIST','REAL',id,id)
        id=0
        id2=0
        call fzero(cout,nMOs)
        do iIrr=0,nIrrep-1
        nd=nBas(iIrr)
        nd2=nd*nd
        do 100 i=1,nd
!        if(occ(i+id).ne.0d0)then
        do 200 j=1,nd
        cout(i+id)=cout(i+id)+occ(i+id)*                                &
     &                ( cmo(j+i*(nd-1)+id2) ** 2 )
200     continue
100     continue
        id=id+nd
        id2=id2+nd2
        enddo

      return
      end
