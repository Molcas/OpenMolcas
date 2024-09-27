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

subroutine plot_XT_no_Exp(label,nT,T,XTcalc)

use Constants, only: Zero, Five, Six
use Definitions, only: wp, iwp, u6

implicit none
character(len=50), intent(in) :: label
integer(kind=iwp), intent(in) :: nT
real(kind=wp), intent(in) :: T(nT), XTcalc(nT)
integer(kind=iwp) :: file_number, file_size, iErr, iT, Length, LuData, LuPlt
real(kind=wp) :: gnuplot_version, tmax, tmin, XTmax, XTmax_calc, XTmin, XTmin_calc
logical(kind=iwp) :: execute_gnuplot_cmd, file_exist, is_file_open
character(len=1023) :: gnuplot_CMD, realname_dat, realname_eps, realname_plt, realname_png
character(len=300) :: cdummy, datafile, epsfile, imagefile, line1, line2, plotfile
integer(kind=iwp), external :: AixRm, IsFreeUnit

#include "macros.fh"

iErr = 0
tmin = Zero
tmax = Zero
XTmin_calc = Zero
XTmax_calc = Zero
XTmin = Zero
XTmax = Zero
tmin = minval(T(:))-0.02_wp*maxval(T(:))
tmax = maxval(T(:))+0.02_wp*maxval(T(:))
XTmin_calc = minval(XTcalc)
XTmax_calc = maxval(XTcalc)
XTmin = XTmin_calc-0.01_wp*XTmax_calc
XTmax = XTmax_calc+0.01_wp*XTmax_calc

#ifdef _DEBUGPRINT_
write(u6,*) 'nT        = ',nT
write(u6,*) 'tmin      = ',tmin
write(u6,*) 'tmax      = ',tmax
write(u6,*) 'XTmin_calc= ',XTmin_calc
write(u6,*) 'XTmax_calc= ',XTmax_calc
write(u6,*) 'XTmin     = ',XTmin
write(u6,*) 'XTmax     = ',XTmax
do iT=1,nT
  write(u6,*) T(iT),XTcalc(iT)
end do
#endif

! generate the GNUPLOT script in the $WorkDir
line1 = ' '
line2 = ' '
gnuplot_CMD = ' '
file_exist = .false.
is_file_open = .false.
file_size = 0
execute_gnuplot_cmd = .false.
gnuplot_version = Zero

! check if the file lineOUT exists
inquire(file='lineOUT',exist=file_exist,opened=is_file_open,number=file_number)

if (file_exist) then
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'file "lineOUT" exists in WorkDir'
# endif
  if (is_file_open) then
#   ifdef _DEBUGPRINT_
    write(u6,'(A)') 'file "lineOUT" is opened'
#   endif
    ! close the file:
    close(unit=file_number,status='DELETE')
  end if
  ! delete the file
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'deleting the file...'
# endif
  iErr = AixRm('lineOUT')
# ifdef _DEBUGPRINT_
  write(u6,*) 'iErr = ',iErr
else
  write(u6,'(A)') 'file "lineOUT" does not exist in WorkDir'
# endif
end if

! find the gnuplot
#ifdef _DEBUGPRINT_
write(u6,'(A)') 'inquire which GNUPLOT'
#endif

!#ifdef __INTEL_COMPILER
call systemf('which gnuplot >> lineOUT',iErr)
#ifdef _DEBUGPRINT_
write(u6,*) 'iErr = ',iErr
#endif
!#else
!call execute_command_line('which gnuplot >> lineOUT')
!#endif

inquire(file='lineOUT',exist=file_exist,opened=is_file_open,number=file_number,size=file_size)

#ifdef _DEBUGPRINT_
write(u6,*) 'File_number =',file_number
write(u6,*) 'Is_file_open=',is_file_open
write(u6,*) 'File_exist  =',file_exist
write(u6,*) 'File_size   =',file_size
#endif

if (file_exist) then
  if (file_size > 0) then
#   ifdef _DEBUGPRINT_
    write(u6,'(A)') 'new file "lineOUT" exists in WorkDir'
#   endif

    file_number = IsFreeUnit(103)
    call molcas_open(file_number,'lineOUT')

    read(file_number,'(A)') line1

#   ifdef _DEBUGPRINT_
    write(u6,*) 'line1=',line1
    write(u6,*) trim(line1)
#   endif
    line2 = trim(line1)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'line2=',line2
#   endif

    close(file_number)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Closing lineOUT file'
    call xFlush(u6)
#   endif
    execute_gnuplot_cmd = .true.
  else
    ! file_size =0
    write(u6,'(A)') 'file "lineOUT" has a size=0. gnuplot was not found on the system.'
    write(u6,'(A)') 'plots will not be created.'
  end if
else
  write(u6,'(A)') 'file "lineOUT" does not exist in WorkDir'
end if
! remove file "lineOUT"
iErr = AixRm('lineOUT')
#ifdef _DEBUGPRINT_
write(u6,*) 'iErr = ',iErr
#endif
!-----------------------------------------------------------------------

! check the version of the gnuplot:
if (execute_gnuplot_cmd) then
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'inquire which version of GNUPLOT is installed'
# endif
  ! attempt to execute the script
  write(gnuplot_CMD,'(2A)') trim(line2),' --version > lineOUT'
# ifdef _DEBUGPRINT_
  write(u6,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
# endif
!# ifdef __INTEL_COMPILER
  call systemf(gnuplot_CMD,iErr)
# ifdef _DEBUGPRINT_
  write(u6,*) 'iErr = ',iErr
# endif
!# else
!  call execute_command_line(gnuplot_CMD)
!# endif
  file_number = IsFreeUnit(102)
  call molcas_open(file_number,'lineOUT')
  read(file_number,*) cdummy,gnuplot_version
  unused_var(cdummy)
# ifdef _DEBUGPRINT_
  write(u6,'(A,F4.1)') 'gnuplot_version = ',gnuplot_version
# endif
  if (abs(gnuplot_version) < 0.1_wp) execute_gnuplot_cmd = .false.
  close(file_number)
  ! remove file "lineOUT"
  iErr = AixRm('lineOUT')
# ifdef _DEBUGPRINT_
  write(u6,*) 'iErr = ',iErr
# endif
end if

!-----------------------------------------------------------------------
! get the true real names of the files on disk:
write(datafile,'(3A)') 'XT_'//trim(label)//'.dat'
write(imagefile,'(3A)') 'XT_'//trim(label)//'.png'
write(epsfile,'(3A)') 'XT_'//trim(label)//'.eps'
write(plotfile,'(3A)') 'XT_'//trim(label)//'.plt'
call prgmtranslate(datafile,realname_dat,Length)
call prgmtranslate(imagefile,realname_png,Length)
call prgmtranslate(epsfile,realname_eps,Length)
call prgmtranslate(plotfile,realname_plt,Length)
#ifdef _DEBUGPRINT_
write(u6,'(3A)') 'realname_dat=',trim(realname_dat)
write(u6,'(3A)') 'realname_png=',trim(realname_png)
write(u6,'(3A)') 'realname_eps=',trim(realname_eps)
write(u6,'(3A)') 'realname_plt=',trim(realname_plt)
#endif

!-----------------------------------------------------------------------
! create the file "XT.dat"
inquire(file=datafile,exist=file_exist,opened=is_file_open,number=file_number)
if (file_exist) iErr = AixRm(trim(datafile))
#ifdef _DEBUGPRINT_
write(u6,*) 'iErr = ',iErr
#endif
LuData = IsFreeUnit(104)
call molcas_open(LuData,datafile)
#ifdef _DEBUGPRINT_
write(u6,*) 'Opening "'//trim(datafile)//'" file'
write(u6,*) 'Opening "'//trim(realname_dat)//'" file'
#endif
do iT=1,nT
  write(LuData,'(2ES24.14)') T(iT),XTcalc(iT)
end do
#ifdef _DEBUGPRINT_
write(u6,*) 'Writing into the "'//trim(datafile)//'" file'
write(u6,*) 'Writing into the "'//trim(realname_dat)//'" file'
#endif
close(LuData)
#ifdef _DEBUGPRINT_
write(u6,*) 'Closing the "'//trim(datafile)//'" file'
write(u6,*) 'Closing the "'//trim(realname_dat)//'" file'
call xFlush(u6)
#endif

!-----------------------------------------------------------------------
! generate the GNUPLOT script in the $WorkDir
inquire(file=plotfile,exist=file_exist,opened=is_file_open,number=file_number)
if (file_exist) iErr = AixRm(trim(plotfile))
#ifdef _DEBUGPRINT_
write(u6,*) 'iErr = ',iErr
#endif
LuPlt = IsFreeUnit(105)
call molcas_open(LuPlt,plotfile)
#ifdef _DEBUGPRINT_
write(u6,*) 'Opening "'//trim(plotfile)//'" file'
write(u6,*) 'Opening "'//trim(realname_plt)//'" file'
#endif

if (gnuplot_version < Five) then
  !===  GNUPLOT VERSION 4 and below ==>>  generate EPS
  write(LuPlt,'(A)') 'set terminal postscript eps enhanced color  size 3.0, 2.0 font "arial, 10"'
  write(LuPlt,'(A)') 'set output "'//trim(realname_eps)//'" '
  write(LuPlt,'(A)') 'set grid'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# Axes'
  write(LuPlt,'(A)') '  set xlabel "Temperature / kelvin" font "Arial,14"'
  write(LuPlt,'(A)') '  set ylabel "{/Symbol c}T / cm^3 K mol^{-1}" font "Arial,14"'
  write(LuPlt,'(A)')
  write(LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[',tmin,':',tmax,']'
  write(LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',XTmin,':',XTmax,']'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# Tics for axes'
  write(LuPlt,'(A)') '  set xtics nomirror'
  write(LuPlt,'(A)') '  set ytics nomirror'
  write(LuPlt,'(A)') '  set mxtics'
  write(LuPlt,'(A)') '  set mytics'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# Margins'
  write(LuPlt,'(A)') '  set bmargin at screen 0.20'
  write(LuPlt,'(A)') '  set lmargin at screen 0.15'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# legend'
  write(LuPlt,'(A)') '  set key right bottom'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# actual plotting'
  write(LuPlt,'(A)') 'plot "'//trim(realname_dat)//'" using 1:2 with lines lt 1  lw  8 lc rgb "red" title "Calculation"'

else if ((gnuplot_version >= Five) .and. (gnuplot_version < Six)) then
  !===  GNUPLOT VERSION 5 and above ==>>  generate PNG
  write(LuPlt,'(A)') 'set terminal pngcairo transparent enhanced font "arial,10" fontscale 4.0 size 1800, 1200'
  write(LuPlt,'(A)') 'set output "'//trim(realname_png)//'" '
  write(LuPlt,'(A)') 'set grid'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# Axes'
  write(LuPlt,'(A)') '  set xlabel "Temperature / kelvin" font "Arial,14"'
  write(LuPlt,'(A)') '  set ylabel "{/Symbol c}T / cm^3 K mol^{-1}" font "Arial,14"'
  write(LuPlt,'(A)')
  write(LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[',tmin,':',tmax,']'
  write(LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',XTmin,':',XTmax,']'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# Tics for axes'
  write(LuPlt,'(A)') '  set xtics nomirror'
  write(LuPlt,'(A)') '  set ytics nomirror'
  write(LuPlt,'(A)') '  set mxtics'
  write(LuPlt,'(A)') '  set mytics'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# Margins'
  write(LuPlt,'(A)') '  set bmargin at screen 0.20'
  write(LuPlt,'(A)') '  set lmargin at screen 0.15'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# legend'
  write(LuPlt,'(A)') '  set key right bottom'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# actual plotting'
  write(LuPlt,'(A)') 'plot "'//trim(realname_dat)//'" using 1:2 with lines lt 1  lw  8 lc rgb "red" title "Calculation"'

else
  write(u6,*) 'GNUPLOT has version: ',gnuplot_version
  write(u6,*) 'This version of GNUPLOT is not known and thus, unsupported.'
  execute_gnuplot_cmd = .false.
end if
write(LuPlt,'(A)')
close(LuPlt)

if (execute_gnuplot_cmd) then
  ! attempt to execute the script
# ifdef _DEBUGPRINT_
  file_exist = .false.
  write(u6,*) trim(realname_plt)
  inquire(file=trim(realname_plt),exist=file_exist,opened=is_file_open,number=file_number)
  if (file_exist) then
    write(u6,'(A,i0,A)') 'File "'//trim(realname_plt)//'" exists.'
  else
    write(u6,'(A,i0,A)') 'File "'//trim(realname_plt)//'" does not exist.'
  end if
  file_exist = .false.
  write(u6,*) trim(realname_dat)
  inquire(file=trim(realname_dat),exist=file_exist,opened=is_file_open,number=file_number)
  if (file_exist) then
    write(u6,'(A,i0,A)') 'File "'//trim(realname_dat)//'" exists.'
  else
    write(u6,'(A,i0,A)') 'File "'//trim(realname_dat)//'" does not exist.'
  end if
# endif

  write(gnuplot_CMD,'(5A)') trim(line2),' < ',trim(realname_plt)
# ifdef _DEBUGPRINT_
  write(u6,'(A,A)') 'gnuplot_CMD=',trim(gnuplot_CMD)
# endif

!# ifdef __INTEL_COMPILER
  call systemf(gnuplot_CMD,iErr)
# ifdef _DEBUGPRINT_
  write(u6,*) 'iErr = ',iErr
# endif
!# else
!  call execute_command_line(gnuplot_CMD)
!# endif

  if (gnuplot_version < Five) then
    file_exist = .false.
    inquire(file=trim(realname_eps),exist=file_exist,opened=is_file_open,number=file_number)
    if (file_exist) then
      write(u6,'(A,i0,A)') 'File "'//trim(realname_eps)//'" was created in Working directory.'
    else
      write(u6,'(A,i0,A)') 'File "'//trim(realname_eps)//'" was NOT created in Working directory.'
    end if
  else
    file_exist = .false.
    inquire(file=trim(realname_png),exist=file_exist,opened=is_file_open,number=file_number)
    if (file_exist) then
      write(u6,'(A,i0,A)') 'File "'//trim(realname_png)//'" was created in Working directory.'
    else
      write(u6,'(A,i0,A)') 'File "'//trim(realname_png)//'" was NOT created in Working directory.'
    end if
  end if
end if

return

end subroutine plot_XT_no_Exp
