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
! Copyright (C) 1991, Per-Olof Widmark                                 *
!               1993,1996,1997, Markus P. Fuelscher                    *
!               1996, Luis Serrano-Andres                              *
!               2002, Roland Lindh                                     *
!               2012, Victor P. Vysotskiy                              *
!***********************************************************************

subroutine FastIO(String)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Manage all I/O operations within the MOLCAS programs.            *
!                                                                      *
!     calling arguments:                                               *
!     String  : Character string                                       *
!               Command liner. Known commands are                      *
!               TRACE  : enable/disable trace output                   *
!               QUERY  : enable/disable updating of the calling tree   *
!               STATUS : print I/O statistics                          *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, IBM Sweden, 1991                                   *
!     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
!     L. Serrano-Andres, University of Lund, Sweden, 1996              *
!                                                                      *
!     Modified to MB, by RL June, 2002, in Tokyo, Japan.               *
!     New I/O Stat, V.P. Vysotskiy, University of Lund, Sweden, 2012   *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)
#include "fio.fh"
#ifdef _OLD_IO_STAT_
#include "ofio.fh"
#else
#include "pfio.fh"
#endif
#include "blksize.fh"
#include "SysDef.fh"
character*(*) String
real*8 tot_Rec_in, tot_Rec_out, file_size
real*8 tot_MB_in, tot_MB_out
#ifdef _OLD_IO_STAT_
real*8 Rec_in, Rec_out, MB_in, MB_out
#else
real*8 tot_Time_out, tot_Time_in, rwr, rrd
#endif
character*100 Fmt

lString = len(String)
if (lString >= 8) then
  ! Enable/disable tracing output
  if (String(1:8) == 'TRACE=ON') Trace = .true.
  if (String(1:9) == 'TRACE=OFF') Trace = .false.

  ! Enable/disable updating of the calling tree
  ! (This option may become very expensive in production runs!)
  if (String(1:8) == 'QUERY=ON') Query = .true.
  if (String(1:9) == 'QUERY=OFF') Query = .false.
end if
! Print I/O statistics
if ((String(1:6) == 'STATUS') .and. (iPrintLevel(-1) >= 3)) then
# ifdef _OLD_IO_STAT_

  call CollapseOutput(1,'I/O statistics')
  write(6,'(3X,A)') '--------------'
  write(6,'(1X,A)')
  write(6,'(1X,A)') ' Unit  Name          Xtent      pages      MBytes     pages      MBytes'
  write(6,'(1X,A)') '                    (MBytes)     in          in        out        out  '
  write(6,'(1X,A)') ' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  tot_Rec_in = 0.0d0
  tot_Rec_out = 0.0d0
  tot_MB_in = 0.0d0
  tot_MB_out = 0.0d0
  do i=1,MxFile
    if (FSCB(i) /= 0) then
      Rec_in = count(1,i)/1024.0d0
      Rec_out = count(3,i)/1024.0d0
      MB_in = count(2,i)/1024.0d0
      MB_out = count(4,i)/1024.0d0
      tot_Rec_in = tot_Rec_in+Rec_in
      tot_Rec_out = tot_Rec_out+Rec_out
      tot_MB_in = tot_MB_in+MB_in
      tot_MB_out = tot_MB_out+MB_out
      ! min_Block_length was undefined long ago...
      file_size = (dble(MxAddr(i))*dble(min_Block_length))/(1024.0d0)**2
      write(6,'(2X,I2,3X,A8,1X,5F11.2)') i,LuName(i),file_size,Rec_in,MB_in,Rec_out,MB_out
    end if
  end do
  write(6,'(1X,A)') ' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  write(6,'(1X,A,20X,4F11.2)') ' Total',tot_Rec_in,tot_MB_in,tot_Rec_out,tot_MB_out
  write(6,'(1X,A)') ' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  call CollapseOutput(0,'I/O statistics')
# else
  tot_Rec_in = 0.0d0
  tot_Rec_out = 0.0d0
  tot_MB_in = 0.0d0
  tot_MB_out = 0.0d0
  tot_Time_in = 0.0d0
  tot_Time_out = 0.0d0
  file_size = 0.0d0
  call CollapseOutput(1,'I/O STATISTICS')

  ! Part I: general I/O information
  write(6,*) ''
  write(6,'(1X,A)') ' I. General I/O information'
  write(6,'(1X,A)') ' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  write(6,'(1X,A)') ' Unit  Name          Flsize      Write/Read            MBytes           Write/Read'
  write(6,'(1X,A)') '                     (MBytes)       Calls              In/Out           Time, sec.'
  write(6,'(1X,A)') ' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

  do i=1,NProfFiles

    tot_Rec_in = tot_Rec_in+PRofData(1,i)
    tot_Rec_out = tot_Rec_out+PRofData(4,i)
    tot_MB_in = tot_MB_in+PRofData(2,i)
    tot_MB_out = tot_MB_out+PRofData(5,i)
    tot_Time_in = tot_Time_in+PRofData(3,i)
    tot_Time_out = tot_Time_out+PRofData(6,i)
    file_size = file_size+FlsSize(i)

    Fmt = '(2X,I2,2X,A8,3X,F11.2,A2,I8,A1,I8,A2,F9.1,A1,F9.1,A2,I8,A1,I8)'
    write(6,Fmt) i,LuNameProf(i),FlsSize(i)/1024.d0/1024.d0,' .',int(PRofData(1,i)),'/',int(PRofData(4,i)),' .', &
                 PRofData(2,i)/1024**2,'/',PRofData(5,i)/1024**2,' .',int(PRofData(3,i)),'/',int(PRofData(6,i))
  end do
  write(6,'(1X,A)') ' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

  Fmt = '(2X,A10,5X,F11.2,A2,I8,A1,I8,A2,F9.1,A1,F9.1,A2,I8,A1,I8)'
  write(6,Fmt) '*  TOTAL ',file_size/1024.d0/1024.d0,' .',int(tot_Rec_in),'/',int(tot_Rec_out),' .',tot_MB_in/1024**2,'/', &
               tot_MB_out/1024**2,' .',int(tot_Time_in),'/',int(tot_Time_out)
  write(6,'(1X,A)') ' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

  !Part II: a little bit more about I/O Access Patterns
  ! Percentage of write/read activity that was random writes/reads.
  write(6,*) ''
  write(6,'(1X,A)') ' II. I/O Access Patterns'
  write(6,'(1X,A)') ' - - - - - - - - - - - - - - - - - - - -'
  write(6,'(1X,A)') ' Unit  Name               % of random'
  write(6,'(1X,A)') '                        Write/Read calls'
  write(6,'(1X,A)') ' - - - - - - - - - - - - - - - - - - - -'
  do i=1,NProfFiles

    if (PRofData(1,i) > 0) then
      rwr = 100*PRofData(7,i)/PRofData(1,i)
    else
      rwr = 0.d0
    end if

    if (PRofData(4,i) > 0) then
      rrd = 100*PRofData(8,i)/PRofData(4,i)
    else
      rrd = 0.d0
    end if

    write(6,'(2X,I2,2X,A8,7X,F9.1,A1,F6.1)') i,LuNameProf(i),rwr,'/',rrd

  end do
  write(6,'(1X,A)') ' - - - - - - - - - - - - - - - - - - - -'
  call CollapseOutput(0,'I/O STATISTICS')
# endif
end if

return

end subroutine FastIO
