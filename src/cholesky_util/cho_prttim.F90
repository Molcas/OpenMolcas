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

subroutine CHO_PRTTIM(SECTION,TCPU2,TCPU1,TWALL2,TWALL1,IOPT)
!
! Purpose: print timing for a section.

use Cholesky, only: LuPri
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: SECTION
real(kind=wp), intent(in) :: TCPU1, TCPU2, TWALL1, TWALL2
integer(kind=iwp), intent(in) :: IOPT
integer(kind=iwp) :: IHRC, IHRW, IMNC, IMNW, LENSEC
real(kind=wp) :: SECC, SECW, TCPUT, TWALLT
character(len=80) :: STRING

TCPUT = TCPU2-TCPU1
TWALLT = TWALL2-TWALL1
call CHO_CNVTIM(TCPUT,IHRC,IMNC,SECC)
call CHO_CNVTIM(TWALLT,IHRW,IMNW,SECW)

if (IOPT == 0) then
  LENSEC = len(SECTION)
  write(LUPRI,'(/,A,A,A)') '***** ',SECTION(1:LENSEC),' completed *****'
  write(LUPRI,'(A,I8,A,I2,A,F6.2,A)') 'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,' seconds'
  write(LUPRI,'(A,I8,A,I2,A,F6.2,A,/)') 'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,' seconds'
else if (IOPT == 1) then
  LENSEC = len(SECTION)
  write(LUPRI,'(///,A,A,A)') '***** ',SECTION(1:LENSEC),' completed *****'
  write(LUPRI,'(A,I8,A,I2,A,F6.2,A)') 'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,' seconds'
  write(LUPRI,'(A,I8,A,I2,A,F6.2,A,//)') 'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,' seconds'
else if (IOPT == 2) then
  LENSEC = min(len(SECTION),70)
  write(STRING,'(A10,A)') 'Timing of ',SECTION(1:LENSEC)
  LENSEC = LENSEC+10
  call CHO_HEAD(STRING(1:LENSEC),'=',80,LUPRI)
  write(LUPRI,'(/,A,I8,A,I2,A,F6.2,A)') 'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,' seconds'
  write(LUPRI,'(A,I8,A,I2,A,F6.2,A)') 'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,' seconds'
else
  write(LUPRI,'(/,A,I8,A,I2,A,F6.2,A)') 'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,' seconds'
  write(LUPRI,'(A,I8,A,I2,A,F6.2,A)') 'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,' seconds'
end if

call XFLUSH(LUPRI)

end subroutine CHO_PRTTIM
