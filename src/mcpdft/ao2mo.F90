!**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine ao2mo(cmo,d_ao,d_mo)
  use definitions,only:iwp,wp
  use constants,only:one,zero
  use stdalloc,only:mma_allocate,mma_deallocate
  implicit none

  real(kind=wp),intent(in) :: cmo(*),d_ao(*)
  real(kind=wp),intent(out) :: d_mo(*)

#include "rasdim.fh"
#include "general.fh"

  integer(kind=iwp) :: ioff1,ioff2,ioff3,isym,ibas,iorb,ifro
  real(kind=wp),allocatable :: tmp1(:),tmp2(:)

!     transform FA from AO to MO basis
  iOff1 = 1
  iOff2 = 1
  iOff3 = 1
  Do iSym = 1,nSym
    iBas = nBas(iSym)
    iOrb = nOrb(iSym)
    If(iBas == 0 .or. iOrb == 0) then
      Cycle
    endif
    iFro = nFro(iSym)
    Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
    Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
    Call Square(d_ao(iOff1),Tmp1,1,iBas,iBas)
    Call DGEMM_('N','N',iBas,iOrb,iBas, &
                1.0d0,Tmp1,iBas, &
                CMO(iOff2+(iFro*iBas)),iBas, &
                0.0d0,Tmp2,iBas)
    Call DGEMM_Tri('T','N',iOrb,iOrb,iBas, &
                   1.0D0,Tmp2,iBas, &
                   CMO(iOff2+(iFro*iBas)),iBas, &
                   0.0D0,d_mo(iOff3),iOrb)
    Call mma_deallocate(Tmp2)
    Call mma_deallocate(Tmp1)
    iOff1 = iOff1+(iBas*iBas+iBas)/2
    iOff2 = iOff2+iBas*iBas
    iOff3 = iOff3+(iOrb*iOrb+iOrb)/2
  EndDo
endsubroutine
