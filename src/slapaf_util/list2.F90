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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

subroutine List2(Title,Lbl,gq,nAtom,nInter,Smmtrc)
!***********************************************************************
!                                                                      *
! Object: to print cartesian internal coordinates.                     *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             1993                                                     *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtom, nInter
character(len=*), intent(in) :: Title, Lbl(nAtom)
real(kind=wp), intent(in) :: gq(3*nAtom,nInter)
logical(kind=iwp), intent(in) :: Smmtrc(3*nAtom)
integer(kind=iwp) :: i, iE, i_F, igq, ii, inc, iq, iQQ, jq, Lu, LuTmp, mInt, MxWdth, nLbl, nRow
real(kind=wp) :: temp
logical(kind=iwp) :: Start
character(len=80) :: Line
character(len=72) :: Frmt
character(len=16) :: filnam
character(len=14) :: qLbl_tmp
character(len=4), allocatable :: qLbl(:)
real(kind=wp), parameter :: Thr = 0.001_wp ! Threshold for printout.

Lu = u6

mInt = 3*nAtom
call mma_allocate(qLbl,mInt,Label='qLbl')

write(Lu,*)
call CollapseOutput(1,'Internal coordinates')
write(Lu,*)
write(Lu,*) ' Specification of the internal coordinates according to the user-defined internal'
write(Lu,*) ' coordinate format.'
write(Lu,*)
write(Lu,'(A)') 'Internal Coordinates'
iq = 0
do igq=1,mInt,3
  if (Smmtrc(igq)) then
    iq = iq+1
    write(qLbl(igq),'(A,I3.3)') 'c',iq
    write(Lu,'(3A)') qLbl(igq),' = Cartesian x ',Lbl((igq+2)/3)
  end if
  if (Smmtrc(igq+1)) then
    iq = iq+1
    write(qLbl(igq+1),'(A,I3.3)') 'c',iq
    write(Lu,'(3A)') qLbl(igq+1),' = Cartesian y ',Lbl((igq+2)/3)
  end if
  if (Smmtrc(igq+2)) then
    iq = iq+1
    write(qLbl(igq+2),'(A,I3.3)') 'c',iq
    write(Lu,'(3A)') qLbl(igq+2),' = Cartesian z ',Lbl((igq+2)/3)
  end if
end do
write(Lu,'(A)') 'Vary'
do iQQ=1,nInter
  write(Line,'(A,I3.3,A)') 'q',iQQ,' ='
  i_F = 7
  jq = 0
  Start = .true.
  do iq=1,mInt
    temp = abs(gq(iq,iQQ))
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
        write(Line(i_F:iE),'(A,F10.8,4A)') ' ',gq(iq,iQQ),' ',qLbl(iq),' '
      else
        iE = i_F+17
        write(Line(i_F:iE),'(A,F10.8,4A)') '+ ',gq(iq,iQQ),' ',qLbl(iq),' '
      end if
      i_F = iE+1
    end if
  end do
  write(Lu,'(A)') Line
end do
write(Lu,'(A)') 'End Of Internal Coordinates'
call CollapseOutput(0,'Internal coordinates')

! Write linear combinations to disc

LuTmp = 11
filnam = 'SPCINX'
call molcas_binaryopen_vanilla(luTmp,filnam)
!open(luTmp,File=filnam,Form='unformatted',Status='unknown')
rewind(LuTmp)

write(LuTmp) mInt,nInter
do iq=1,mInt
  qLbl_tmp = qLbl(iq)
  write(LuTmp) qLbl_tmp,(gq(iq,iQQ),iQQ=1,nInter)
end do
call mma_deallocate(qLbl)

close(LuTmp)

write(Lu,*)
call CollapseOutput(1,Title)

MxWdth = 132
nLbl = 8+1
nRow = 9
inc = min((MxWdth-nLbl)/nRow,nInter)

do ii=1,nInter,inc
  write(Lu,*)
  write(Frmt,'(A,I2,A)') '(A,1X,',inc,'(I5,4X))'
  write(Lu,Frmt) 'Internal',(i,i=ii,min(ii+inc-1,nInter))
  write(Lu,*)
  write(Frmt,'(A,I2,A)') '(A4,A4,1X,',inc,'(F8.5,1X))'
  do igq=1,mInt,3
    write(Lu,Frmt) Lbl((igq+2)/3),' x  ',(gq(igq,i),i=ii,min(ii+inc-1,nInter))
    write(Lu,Frmt) Lbl((igq+2)/3),' y  ',(gq(igq+1,i),i=ii,min(ii+inc-1,nInter))
    write(Lu,Frmt) Lbl((igq+2)/3),' z  ',(gq(igq+2,i),i=ii,min(ii+inc-1,nInter))
  end do
  write(Lu,*)
end do
call CollapseOutput(0,Title)

return

end subroutine List2
