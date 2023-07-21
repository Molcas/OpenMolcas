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

subroutine CHO_ANASIZE(VEC,LVEC,BIN,LBIN,LUPRI)
!
! Purpose: analyse vector (histogram).

implicit real*8(a-h,o-z)
real*8 VEC(LVEC), BIN(LBIN)
parameter(ZERO=0.0d0)
parameter(MBIN=20)
integer ICOUNT(MBIN)
logical FOUND

! Return if nothing to do.
! ------------------------

if ((LVEC < 1) .or. (LBIN < 1)) return

! Ensure that BIN is in descending order.
! ---------------------------------------

IJOB = -1
call CHO_ORDER(BIN,LBIN,IJOB)

! Test that BIN is positive.
! --------------------------

if (BIN(1) <= ZERO) return

! Analysis.
! ---------

NBIN = min(LBIN,MBIN)
call IZERO(ICOUNT,NBIN)
NLOW = 0
NZER = 0
NNEG = 0
XNEG = ZERO

do I=1,LVEC

  TEST = VEC(I)

  if (TEST < ZERO) then
    NNEG = NNEG+1
    XNEG = min(XNEG,TEST)
  else if (TEST == ZERO) then
    NZER = NZER+1
  end if

  IBIN = 0
  FOUND = .false.
  do while ((.not. FOUND) .and. (IBIN < NBIN))
    IBIN = IBIN+1
    if (TEST >= BIN(IBIN)) then
      ICOUNT(IBIN) = ICOUNT(IBIN)+1
      FOUND = .true.
    end if
  end do
  if (.not. FOUND) NLOW = NLOW+1

end do

! Print.
! ------

TOPCT = 1.0d2/dble(LVEC)

JCOUNT = ICOUNT(1)
write(LUPRI,'(/,1X,A,11X,D11.4,A,I12,1X,F7.2,A,3X,A,F7.2,A)') 'Larger than ',BIN(1),':',ICOUNT(1),dble(ICOUNT(1))*TOPCT,'%', &
                                                              'Accumulated: ',dble(JCOUNT)*TOPCT,'%'
do IBIN=2,NBIN
  JCOUNT = JCOUNT+ICOUNT(IBIN)
  write(LUPRI,'(1X,A,D11.4,A,D11.4,A,I12,1X,F7.2,A,3X,A,F7.2,A)') 'Between ',BIN(IBIN-1),' and ',BIN(IBIN),':',ICOUNT(IBIN), &
                                                                  dble(ICOUNT(IBIN))*TOPCT,'%','Accumulated: ', &
                                                                  dble(JCOUNT)*TOPCT,'%'
end do
JCOUNT = JCOUNT+NLOW
write(LUPRI,'(1X,A,10X,D11.4,A,I12,1X,F7.2,A,3X,A,F7.2,A)') 'Smaller than ',BIN(NBIN),':',NLOW,dble(NLOW)*TOPCT,'%', &
                                                            'Accumulated: ',dble(JCOUNT)*TOPCT,'%'

write(LUPRI,'(/,1X,A,I12,1X,F7.2,A)') 'Number of elements exactly 0.0D0 :',NZER,dble(NZER)*TOPCT,'%'
write(LUPRI,'(1X,A,I12,1X,F7.2,A)') 'Number of negative elements      :',NNEG,dble(NNEG)*TOPCT,'%'
if (NNEG > 0) write(LUPRI,'(1X,A,D12.4)') ' - numerically largest           :',XNEG

end subroutine CHO_ANASIZE
