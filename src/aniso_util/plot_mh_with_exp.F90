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

subroutine plot_MH_with_Exp(nH,H,nTempMagn,TempMagn,MHcalc,MHexp,zJ)

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: nH, nTempMagn
real(wp), intent(in) :: H(nH), TempMagn(nTempMagn)
real(wp), intent(in) :: MHexp(nH,nTempMagn)
real(wp), intent(in) :: MHcalc(nH,nTempMagn)
real(wp), intent(in) :: zJ
! local variables
real(wp) :: hmin, hmax, MHmin_exp, MHmax_exp, MHmin_calc, MHmax_calc, MHmin, MHmax
real(wp) :: gnuplot_version
integer :: file_number, iH, iTempMagn, LuPlt, LuData, file_size, StdOut
logical :: file_exist, is_file_open, execute_gnuplot_cmd, dbg
character(len=300) :: line1, line2, cdummy
character(len=300) :: datafile, plotfile, imagefile, epsfile
integer, external :: AixRm
integer :: Length
character(len=1023) :: realname_plt, realname_dat, realname_png, realname_eps, gnuplot_CMD
integer :: iErr, ifilenumber
integer, external :: IsFreeUnit

dbg = .false.
iErr = 0
StdOut = 6
hmin = 0.0_wp
hmax = 0.0_wp
MHmin_exp = 0.0_wp
MHmax_exp = 0.0_wp
MHmin_calc = 0.0_wp
MHmax_calc = 0.0_wp
MHmin = 0.0_wp
MHmax = 0.0_wp
hmin = minval(H)-0.05_wp*maxval(H)
hmax = maxval(H)+0.05_wp*maxval(H)
MHmin_exp = minval(MHexp)
MHmax_exp = maxval(MHexp)
MHmin_calc = minval(MHcalc)
MHmax_calc = maxval(MHcalc)
MHmin = min(MHmin_exp,MHmin_calc)-0.08_wp*max(MHmax_exp,MHmax_calc)
MHmax = max(MHmax_exp,MHmax_calc)+0.08_wp*max(MHmax_exp,MHmax_calc)

if (dbg) then
  write(StdOut,*) 'nH        = ',nH
  write(StdOut,*) 'hmin      = ',hmin
  write(StdOut,*) 'hmax      = ',hmax
  write(StdOut,*) 'MHmin_exp = ',MHmin_exp
  write(StdOut,*) 'MHmax_exp = ',MHmax_exp
  write(StdOut,*) 'MHmin_calc= ',MHmin_calc
  write(StdOut,*) 'MHmax_calc= ',MHmax_calc
  write(StdOut,*) 'MHmin     = ',MHmin
  write(StdOut,*) 'MHmax     = ',MHmax
  write(StdOut,*) 'zJ        = ',zJ
  do iH=1,nH
    write(StdOut,*) H(iH),(MHcalc(iH,iTempMagn),MHexp(iH,iTempMagn),iTempMagn=1,nTempMagn)
  end do
end if

! generate the GNUPLOT script in the $WorkDir
line1 = ' '
line2 = ' '
gnuplot_CMD = ' '
file_exist = .false.
is_file_open = .false.
file_size = 0
execute_gnuplot_cmd = .false.
gnuplot_version = 0.0_wp

! check if the file lineOUT exists
inquire(file='lineOUT',exist=file_exist,opened=is_file_open,number=file_number)

if (file_exist) then
  if (dbg) write(StdOut,'(A)') 'file "lineOUT" exists in WorkDir'
  if (is_file_open) then
    if (dbg) write(StdOut,'(A)') 'file "lineOUT" is opened'
    ! close the file:
    close(unit=file_number,status='DELETE')
  end if
  ! delete the file
  if (dbg) write(StdOut,'(A)') 'deleting the file...'
  iErr = AixRm('lineOUT')
  if (dbg) write(StdOut,*) 'iErr = ',iErr
else
  if (dbg) write(StdOut,'(A)') 'file "lineOUT" does not exist in WorkDir'
end if

! find the gnuplot
if (dbg) write(StdOut,'(A)') 'inquire which GNUPLOT'

!#ifdef __INTEL_COMPILER
call systemf('which gnuplot >> lineOUT',iErr)
if (dbg) write(StdOut,*) 'iErr = ',iErr
!#else
!call execute_command_line('which gnuplot >> lineOUT')
!#endif

inquire(file='lineOUT',exist=file_exist,opened=is_file_open,number=file_number,size=file_size)

if (dbg) then
  write(StdOut,*) 'File_number =',file_number
  write(StdOut,*) 'Is_file_open=',is_file_open
  write(StdOut,*) 'File_exist  =',file_exist
  write(StdOut,*) 'File_size   =',file_size
end if

if (file_exist) then
  if (file_size > 0) then
    if (dbg) write(StdOut,'(A)') 'new file "lineOUT" exists in WorkDir'

    file_number = IsFreeUnit(93)
    call molcas_open(file_number,'lineOUT')

    read(file_number,'(A)') line1

    if (dbg) then
      write(StdOut,*) 'line1=',line1
      write(StdOut,*) trim(line1)
    end if
    line2 = trim(line1)
    if (dbg) write(StdOut,*) 'line2=',line2

    close(file_number)
    if (dbg) write(StdOut,*) 'Closing lineOUT file'
    execute_gnuplot_cmd = .true.
  else
    ! file_size =0
    write(StdOut,'(A)') 'file "lineOUT" has a size=0. gnuplot was not found on the system.'
    write(StdOut,'(A)') 'plots will not be created.'
  end if
else
  write(StdOut,'(A)') 'file "lineOUT" does not exist in WorkDir'
end if
! remove file "lineOUT"
iErr = AixRm('lineOUT')
if (dbg) write(StdOut,*) 'iErr = ',iErr
!-----------------------------------------------------------------------

! check the version of the gnuplot:
if (execute_gnuplot_cmd) then
  if (dbg) write(StdOut,'(A)') 'inquire which version of GNUPLOT is installed'
  ! attempt to execute the script
  write(gnuplot_CMD,'(2A)') trim(line2),' --version > lineOUT'
  if (dbg) write(StdOut,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
!# ifdef __INTEL_COMPILER
  call systemf(gnuplot_CMD,iErr)
  if (dbg) write(StdOut,*) 'iErr = ',iErr
!# else
!  call execute_command_line(gnuplot_CMD)
!# endif
  file_number = IsFreeUnit(94)
  call molcas_open(file_number,'lineOUT')
  read(file_number,*) cdummy,gnuplot_version
  if (dbg) write(StdOut,'(A,F4.1)') 'gnuplot_version = ',gnuplot_version
  if (abs(gnuplot_version) < 0.1_wp) execute_gnuplot_cmd = .false.
  close(file_number)
  ! remove file "lineOUT"
  iErr = AixRm('lineOUT')
  if (dbg) write(StdOut,*) 'iErr = ',iErr
end if

!write(datafile,'(3A)') 'MH.dat'
!write(imagefile,'(3A)') 'MH.png'
!write(epsfile,'(3A)') 'MH.eps'
!write(plotfile ,'(3A)') 'MH.plt'
!call prgmtranslate(datafile,realname_dat,Length)
!call prgmtranslate(imagefile,realname_png,Length)
!call prgmtranslate(epsfile,realname_eps,Length)
!call prgmtranslate(plotfile,realname_plt,Length)
!if (dbg) then
!  write(StdOut,'(3A)') 'realname_dat=',trim(realname_dat)
!  write(StdOut,'(3A)') 'realname_png=',trim(realname_png)
!  write(StdOut,'(3A)') 'realname_eps=',trim(realname_eps)
!  write(StdOut,'(3A)') 'realname_plt=',trim(realname_plt)
!end if
!-----------------------------------------------------------------------
do iTempMagn=1,nTempMagn

  ! generate the file "MH.dat":
  write(datafile,'(A,I0,A)') 'MH_T_',iTempMagn,'.dat'
  write(imagefile,'(A,I0,A)') 'MH_T_',iTempMagn,'.png'
  write(epsfile,'(A,I0,A)') 'MH_T_',iTempMagn,'.eps'
  write(plotfile,'(A,I0,A)') 'MH_T_',iTempMagn,'.plt'
  call prgmtranslate(datafile,realname_dat,Length)
  call prgmtranslate(imagefile,realname_png,Length)
  call prgmtranslate(epsfile,realname_eps,Length)
  call prgmtranslate(plotfile,realname_plt,Length)

  inquire(file=datafile,exist=file_exist,opened=is_file_open,number=file_number)
  if (file_exist) iErr = AixRm(trim(datafile))
  if (dbg) write(StdOut,*) 'iErr = ',iErr
  ifilenumber = 0
  ifilenumber = 95+iTempMagn
  LuData = IsFreeUnit(ifilenumber)
  call molcas_open(LuData,datafile)
  if (dbg) then
    write(StdOut,*) 'Opening "'//trim(datafile)//'" file. iTempMagn=',iTempMagn
    write(StdOut,*) 'Opening "'//trim(realname_dat)//'" file. iTempMagn=',iTempMagn
  end if
  do iH=1,nH
    write(LuData,'(3ES24.14)') H(iH),MHexp(iH,iTempMagn),MHcalc(iH,iTempMagn)
  end do
  if (dbg) then
    write(StdOut,*) 'Writing into the "'//trim(datafile)//'" file. iTempMagn=',iTempMagn
    write(StdOut,*) 'Writing into the "'//trim(realname_dat)//'" file. iTempMagn=',iTempMagn
  end if
  close(LuData)
  if (dbg) then
    write(StdOut,*) 'Closing the "'//trim(datafile)//'" file. iTempMagn=',iTempMagn
    write(StdOut,*) 'Closing the "'//trim(realname_dat)//'" file. iTempMagn=',iTempMagn
  end if
  flush(StdOut)

  ! generate the GNUPLOT script in the $WorkDir

  inquire(file=plotfile,exist=file_exist,opened=is_file_open,number=file_number)
  if (file_exist) iErr = AixRm(trim(plotfile))
  if (dbg) write(StdOut,*) 'iErr = ',iErr
  ifilenumber = 0
  ifilenumber = 105+iTempMagn
  LuPlt = IsFreeUnit(ifilenumber)
  call molcas_open(LuPlt,plotfile)
  if (dbg) then
    write(StdOut,*) 'Opening "'//trim(plotfile)//'" file'
    write(StdOut,*) 'Opening "'//trim(realname_plt)//'" file'
  end if

  ! EPS or PNG images:

  if (gnuplot_version < 5.0_wp) then
    !===  GNUPLOT VERSION 4 and below ==>>  generate EPS
    write(LuPlt,'(A)') 'set terminal postscript eps enhanced color  size 3.0, 2.0 font "arial, 10"'
    write(LuPlt,'(A)') 'set output "'//trim(realname_eps)//'" '
    write(LuPlt,'(A)') 'set grid'
    write(LuPlt,'(A)')
    write(LuPlt,'(A)') '# Axes'
    write(LuPlt,'(A)') '  set xlabel "Magnetic Field / Tesla" font "Arial,14"'
    write(LuPlt,'(A)') '  set ylabel "M / {/Symbol m}_{B}" font "Arial,14"'
    write(LuPlt,'(A)')
    write(LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[',hmin,':',hmax,']'
    write(LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',MHmin,':',MHmax,']'
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

    write(LuPlt,'(A,F7.3,A)') 'plot "'//trim(realname_dat)//'" using 1:2 with points  lt 1  lw 3 lc rgb "black"  title "Exp. T=', &
                              TempMagn(iTempMagn),' K.", \'
    write(LuPlt,'(A,F7.3,A)') '     "'//trim(realname_dat)//'" using 1:3 with lines   lt 1  lw 8 lc rgb "red"    title "Calc.T=', &
                              TempMagn(iTempMagn),' K."'
    write(LuPlt,'(A)')

  else if ((gnuplot_version >= 5.0_wp) .and. (gnuplot_version < 6.0_wp)) then
    !===  GNUPLOT VERSION 5 and above ==>>  generate PNG
    write(LuPlt,'(A)') 'set terminal pngcairo transparent enhanced font "arial,10" fontscale 4.0 size 1800, 1200'
    write(LuPlt,'(A)') 'set output "'//trim(realname_png)//'" '
    write(LuPlt,'(A)') 'set grid'
    write(LuPlt,'(A)')
    write(LuPlt,'(A)') '# Axes'
    write(LuPlt,'(A)') '  set xlabel "Magnetic Field / Tesla" font "Arial,14"'
    write(LuPlt,'(A)') '  set ylabel "M / {/Symbol m}_{B}" font "Arial,14"'
    write(LuPlt,'(A)')
    write(LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[',hmin,':',hmax,']'
    write(LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',MHmin,':',MHmax,']'
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

    write(LuPlt,'(A,F7.3,A)') 'plot "'//trim(realname_dat)//'" using 1:2 with circles lt 1  lw 3 lc rgb "black"  title "Exp. T=', &
                              TempMagn(iTempMagn),' K.", \'
    write(LuPlt,'(A,F7.3,A)') '     "'//trim(realname_dat)//'" using 1:3 with lines   lt 1  lw 8 lc rgb "red"    title "Calc.T=', &
                              TempMagn(iTempMagn),' K."'
    write(LuPlt,'(A)')
  else
    write(StdOut,*) 'GNUPLOT has version: ',gnuplot_version
    write(StdOut,*) 'This version of GNUPLOT is not known and thus, unsupported.'
    execute_gnuplot_cmd = .false.
  end if
  close(LuPlt)

  if (execute_gnuplot_cmd) then
    ! attempt to execute the script
    if (dbg) then
      file_exist = .false.
      write(StdOut,*) trim(realname_plt)
      inquire(file=trim(realname_plt),exist=file_exist,opened=is_file_open,number=file_number)
      if (file_exist) then
        write(StdOut,'(A,i0,A)') 'File "'//trim(realname_plt)//'" exists.'
      else
        write(StdOut,'(A,i0,A)') 'File "'//trim(realname_plt)//'" does not exist.'
      end if
      file_exist = .false.
      write(StdOut,*) trim(realname_dat)
      inquire(file=trim(realname_dat),exist=file_exist,opened=is_file_open,number=file_number)
      if (file_exist) then
        write(StdOut,'(A,i0,A)') 'File "'//trim(realname_dat)//'" exists.'
      else
        write(StdOut,'(A,i0,A)') 'File "'//trim(realname_dat)//'" does not exist.'
      end if
    end if

    write(gnuplot_CMD,'(5A)') trim(line2),' < ',trim(realname_plt)
    if (dbg) write(StdOut,'(A,A)') 'gnuplot_CMD=',trim(gnuplot_CMD)

!#   ifdef __INTEL_COMPILER
    call systemf(gnuplot_CMD,iErr)
    if (dbg) write(StdOut,*) 'iErr = ',iErr
!#   else
!    call execute_command_line(gnuplot_CMD)
!#   endif

    if (gnuplot_version < 5.0_wp) then
      file_exist = .false.
      inquire(file=trim(realname_eps),exist=file_exist,opened=is_file_open,number=file_number)
      if (file_exist) then
        write(StdOut,'(A,i0,A)') 'File "'//trim(realname_eps)//'" was created in Working directory.'
      else
        write(StdOut,'(A,i0,A)') 'File "'//trim(realname_eps)//'" was NOT created in Working directory.'
      end if
    else
      file_exist = .false.
      inquire(file=trim(realname_png),exist=file_exist,opened=is_file_open,number=file_number)
      if (file_exist) then
        write(StdOut,'(A,i0,A)') 'File "'//trim(realname_png)//'" was created in Working directory.'
      else
        write(StdOut,'(A,i0,A)') 'File "'//trim(realname_png)//'" was NOT created in Working directory.'
      end if
    end if
  end if

end do ! iTempMagn

return
#ifdef _WARNING_WORKAROUND_
if (.false.) call UNUSED_CHARACTER(cdummy)
#endif

end subroutine plot_MH_with_Exp
