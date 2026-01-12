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

subroutine FreeStr()

use Str_Info, only: ITYP_Dummy, Str, Str_Hidden
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ITYP

do ITYP=1,ITYP_Dummy
  call mma_deallocate(Str_Hidden(ITYP)%OCSTR,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%STREO,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%STSM,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%STCL,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%NSTSO,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%ISTSO,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%EL1,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%EL3,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%EL123,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%Z,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%STSTMI,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%STSTMN,safe='*')
  call mma_deallocate(Str_Hidden(ITYP)%STSTM,safe='*')
  nullify(Str(ITYP)%OCSTR,Str(ITYP)%STREO,Str(ITYP)%STSM,Str(ITYP)%STCL,Str(ITYP)%NSTSO,Str(ITYP)%ISTSO,Str(ITYP)%EL1, &
          Str(ITYP)%EL3,Str(ITYP)%EL123,Str(ITYP)%Z,Str(ITYP)%STSTMI,Str(ITYP)%STSTMN,Str(ITYP)%STSTM)
end do

end subroutine
