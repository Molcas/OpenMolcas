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

subroutine SBDIAG()

use caspt2_global, only: iPrGlb
use PrintLevel, only: USUAL, VERBOSE
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use caspt2_module, only: nSym, ThrShn, ThrShs, Cases, nASup, nISup, nInDep
use definitions, only: iwp, wp, u6

implicit none
real(kind=wp) CondNr, CPU
integer(kind=iwp) iCase, iPar0, iPar1, iSym

if (IPRGLB >= VERBOSE) then
  write(u6,*)
  write(u6,*) ' Find transformation matrices to eigenbasis of block-diagonal part of H0.'
  write(u6,*) ' Eliminate linear dependency. Thresholds for:'
  write(u6,'(A,G12.4)') '   Initial squared norm  :',THRSHN
  write(u6,'(A,G12.4)') '   Eigenvalue of scaled S:',THRSHS
end if

! SVC.20100904: there are now two SBDIAG versions: a replicate
! subroutine expecting replicate S and B matrices (in upper triangular
! column-wise storage) and performing transformations on each process,
! and a global array subroutine expecting S and B matrices in global
! arrays (full column-wise storage) spread over processes.  The latter
! is currently used only for cases 1 (A) and 4 (C), for which global
! array mksmat and mkbmat routines have been implemented, as these can
! grow very big with increasing size of the active space.  For now, we
! still use the replicate routines for the other cases as they have more
! modest array sizes.

if (IPRGLB >= VERBOSE) then
  write(u6,*)
  write(u6,*) ' Condition numbers are computed after diagonal scaling and after removal of'
  write(u6,*) ' linear dependency. Resulting sizes, condition numbers, and times:'
  write(u6,'(3X,A10,4A12,A9)') 'CASE(SYM)','NASUP','NISUP','NINDEP','COND NR','CPU (s)'
end if

do iCASE=1,11
  do ISYM=1,NSYM
#   ifdef _MOLCAS_MPP_
    if (IS_REAL_PAR() .and. ((ICASE == 1) .or. (ICASE == 4))) then
      call SBDIAG_MPP(ISYM,ICASE,CONDNR,CPU)
    else
#   endif
      call SBDIAG_SER(ISYM,ICASE,CONDNR,CPU)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
    if (IPRGLB >= VERBOSE) write(u6,'(3X,A6,A1,I1,A1,1X,3I12,G11.2,I9)') &
      CASES(ICASE),'(',ISYM,')',NASUP(ISYM,ICASE),NISUP(ISYM,ICASE),NINDEP(ISYM,ICASE),CONDNR,nint(CPU)

  end do
end do

! usually print info on the total number of parameters
IPAR0 = 0
IPAR1 = 0
do ICASE=1,13
  do ISYM=1,NSYM
    IPAR0 = IPAR0+NASUP(ISYM,ICASE)*NISUP(ISYM,ICASE)
    IPAR1 = IPAR1+NINDEP(ISYM,ICASE)*NISUP(ISYM,ICASE)
  end do
end do
if (IPRGLB >= USUAL) then
  write(u6,*)
  write(u6,*) ' Total nr of CASPT2 parameters:'
  write(u6,'(a,i12)') '   Before reduction:',IPAR0
  write(u6,'(a,i12)') '   After  reduction:',IPAR1
end if

end subroutine SBDIAG
