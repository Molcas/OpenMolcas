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
!***********************************************************************
! SubRoutine GATskL(create,nTsk,igaTsk)                                *
!  -> initialize or kill global task list                              *
!     create:   Logical: .TRUE. -> create / .FALSE. -> kill            *
!     nTsk:     # of tasks                                             *
!     igaTsk:   global array handle to global task list (on return)    *
!***********************************************************************

subroutine GATskL(create,nTsk,igaTsk)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use Definitions, only: u6
#endif
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: create
integer(kind=iwp) :: nTsk, igaTsk
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: Chk, nProcs
logical(kind=iwp) :: ok
#include "global.fh"
#include "mafdecls.fh"

if (.not. Is_Real_Par()) return
if (create) then
  nProcs = ga_nnodes()
  Chk = (nTsk+nProcs-1)/nProcs
  ok = ga_create(MT_INT,nTsk,1,'GlTskL',Chk,1,igaTsk)
  if (.not. ok) then
    write(u6,*) 'GATskL: ga_create not OK!'
    call GAStp('GATskL',42)
    call Abend()
  end if
  call ga_zero(igaTsk)
else
  ok = ga_destroy(igaTsk)
end if
#else
#include "macros.fh"
unused_var(create)
unused_var(nTsk)
unused_var(igaTsk)
#endif

end subroutine GATskL
