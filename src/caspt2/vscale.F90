!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

subroutine V_SCALE(EIG,SCA,V,nRows,NAS,LDV,NIN,COND)

use caspt2_module, only: ThrShS
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nRows, NAS, LDV, NIN
real(kind=wp), intent(in) :: EIG(NAS), SCA(NAS)
real(kind=wp), intent(inout) :: V(LDV,*)
real(kind=wp), intent(out) :: COND(NIN)
integer(kind=iwp) :: I, iVec, J, jVEC
real(kind=wp) :: SZ

jVEC = 0
do J=1,NAS
  if (EIG(J) >= THRSHS) then
    jVEC = jVEC+1
    V(:,jVEC) = V(:,J)/sqrt(EIG(J))
  end if
end do
if (jVEC /= NIN) then
  write(u6,*) 'V_SCALE: inconsitency in linear dependence removal, ABORT'
  call AbEnd()
end if
! Addition, for the scaled symmetric ON.
do I=1,NIN
  V(1:nRows,I) = SCA(1:nRows)*V(1:nRows,I)
end do
! The condition number, after scaling, disregarding linear dep.
if (NIN >= 2) then
  do jVEC=1,NIN
    SZ = Zero
    do iVEC=1,nRows
      SZ = SZ+V(iVEC,jVEC)**2
    end do
    COND(jVEC) = SZ
  end do
end if
end subroutine V_SCALE

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(V_SCALE)

#endif
