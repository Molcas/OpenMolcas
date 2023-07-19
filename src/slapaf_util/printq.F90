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

subroutine PrintQ(rK,qLbl,nq,nQQ,LuIC,rMult)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nq, nQQ, LuIC
real(kind=wp), intent(in) :: rK(nq,nQQ), rMult(nq)
character(len=14), intent(in) :: qLbl(nq)
#include "print.fh"
integer(kind=iwp) :: i, iE, i_F, iiQQ, IncQQ, iPrint, iq, iQQ, iRout, istatus, jq, Lu, LuTmp, mQQ
real(kind=wp) :: temp
logical(kind=iwp) :: Start
character(len=80) :: Line
character(len=16) :: filnam
real(kind=wp), parameter :: Thr = 0.001_wp ! Threshold for printout.
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
iPrint = 99
#else
iRout = 122
iPrint = nPrint(iRout)
#endif

Lu = u6

if (iPrint > 5) then
  write(Lu,*)
  call CollapseOutput(1,'Internal coordinate section')

  if (iPrint >= 6) then
    write(Lu,*)
    write(Lu,'(80A)') ('*',i=1,80)
    write(Lu,*) ' Auto-Defined Internal coordinates'
    write(Lu,'(80A)') ('-',i=1,80)
    write(Lu,'(A)') '  Primitive Internal Coordinates:'
  else
    write(Lu,*)
    write(Lu,'(A)') '  Redundant Internal Coordinates:'
    write(Lu,*)
  end if

  rewind(LuIC)
  do
    read(LuIC,'(A)',iostat=istatus) Line
    if (istatus < 0) exit
    write(Lu,'(A)') Line
  end do
  rewind(LuIC)
  if (iPrint >= 6) then

    write(Lu,'(A)') '  Internal Coordinates:'
    do iQQ=1,nQQ
      write(Line,'(A,I3.3,A)') 'q',iQQ,' ='
      i_F = 7
      jq = 0
      Start = .true.
      do iq=1,nq
        temp = abs(rK(iq,iQQ))
        if (temp > Thr) then
          jq = jq+1
          if (jq > 4) then
            Line(80:80) = '&'
            write(Lu,'(A)') Line
            Line = ' '
            i_F = 6
            jq = 1
            Start = .false.
          end if
          if ((jq == 1) .and. Start) then
            iE = i_F+16
            write(Line(i_F:iE),'(A,F10.8,4A)') ' ',rK(iq,iQQ),' ',qLbl(iq)(1:4),' '
          else
            iE = i_F+17
            write(Line(i_F:iE),'(A,F10.8,4A)') '+ ',rK(iq,iQQ),' ',qLbl(iq)(1:4),' '
          end if
          i_F = iE+1
        end if
      end do
      write(Lu,'(A)') Line
    end do
    write(Lu,'(80A)') ('*',i=1,80)
    call CollapseOutput(0,'Internal coordinate section')
  end if
end if

! Write linear combinations to disc

LuTmp = 11
filnam = 'SPCINX'
call molcas_binaryopen_vanilla(luTmp,filnam)
!open(luTmp,File=filnam,Form='unformatted',Status='unknown')
rewind(LuTmp)

! put in degeneracy factor so that ddot will work.

write(LuTmp) nq,nQQ
do iq=1,nq
  write(LuTmp) qLbl(iq),rMult(iq)*rK(iq,:)
end do

close(LuTmp)

if ((iPrint >= 10) .and. (nQQ <= 12)) then
  write(Lu,*)
  write(Lu,*) ' Nonredundant internal coordinates'
end if
if ((iPrint >= 6) .and. (nQQ <= 12)) then
  write(Lu,*)
  write(Lu,*) ' Number of redundant coordinates:',nq
  write(Lu,*)
end if
if ((iPrint >= 10) .and. (nQQ <= 12)) then
  write(Lu,'(A,E10.3)') ' Threshold for printout:',Thr
  IncQQ = 8
  do iiQQ=1,nQQ,IncQQ
    mQQ = min(nQQ,iiQQ+IncQQ-1)
    write(Lu,*)
    write(Lu,'(14X,8I10)') (iQQ,iQQ=iiQQ,mQQ)
    do iq=1,nq
      temp = sqrt(DDot_(nQQ,rK(iq,1),nq,rK(iq,1),nq))
      if (temp > Thr) write(Lu,'(A,8F10.6)') qLbl(iq),rK(iq,iiQQ:mQQ)
    end do
    write(Lu,*)
  end do
end if

return

end subroutine PrintQ
