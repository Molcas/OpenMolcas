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
! Copyright (C) 2011, Steven Vancoillie                                *
!***********************************************************************
!***********************************************************************
! New VECTOR UTILITIES, written by Steven Vancoillie, May 2011
! A set of subroutines that can transform RHS arrays using the parallel
! aware subroutines
!***********************************************************************

subroutine PSCAVEC(FACT,IVEC,JVEC)
! Scale vector nr IVEC with scale factor FACT and put the result in
! vector nr JVEC: |JVEC> <- FACT * |IVEC>

use caspt2_global, only: iPrGlb
use PrintLevel, only: USUAL
use caspt2_module, only: CPUSCA, nCases, nInDep, niSup, nSym, TIOSCA
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: FACT
integer(kind=iwp), intent(in) :: IVEC, JVEC
integer(kind=iwp) :: ICASE, ISYM, lg_V, NIN, NIS
real(kind=wp) :: CPU, CPU0, CPU1, SIGMA2, TIO, TIO0, TIO1
real(kind=wp), external :: RHS_DDOT

call TIMING(CPU0,CPU,TIO0,TIO)

if ((FACT == One) .and. (IVEC == JVEC)) return
SIGMA2 = Zero
do ICASE=1,NCASES
  do ISYM=1,NSYM
    NIN = NINDEP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (NIN*NIS /= 0) then
      call RHS_ALLO(NIN,NIS,lg_V)
      call RHS_READ(NIN,NIS,lg_V,ICASE,ISYM,IVEC)
      call RHS_SCAL(NIN,NIS,lg_V,FACT)
      if (FACT == -One) SIGMA2 = SIGMA2+RHS_DDOT(NIN,NIS,lg_V,lg_V)
      call RHS_SAVE(NIN,NIS,lg_V,ICASE,ISYM,JVEC)
      call RHS_FREE(lg_V)
    end if
  end do
end do
if ((IPRGLB >= USUAL) .and. (FACT == -One)) then ! it is at ITER=0
  write(u6,*)
  write(u6,'(1x,a,f18.10)') 'Variance of |WF0>: ',SIGMA2
end if

call TIMING(CPU1,CPU,TIO1,TIO)
CPUSCA = CPUSCA+(CPU1-CPU0)
TIOSCA = TIOSCA+(TIO1-TIO0)

end subroutine PSCAVEC
