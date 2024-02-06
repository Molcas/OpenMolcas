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

subroutine plot_MH_no_Exp(nH,H,nTempMagn,TempMagn,MHcalc,zJ)

use Constants, only: Zero, Five, Six
use Definitions, only: wp, u6

implicit none
integer, intent(in) :: nH, nTempMagn
real(wp), intent(in) :: H(nH), TempMagn(nTempMagn)
real(wp), intent(in) :: MHcalc(nH,nTempMagn)
real(wp), intent(in) :: zJ
! local variables
real(wp) :: hmin, hmax, MHmin_calc, MHmax_calc, MHmin, MHmax, r
real(wp) :: gnuplot_version
integer :: file_number, iH, iTempMagn, LuPlt, LuData, ik, ic, file_size
logical :: file_exist, is_file_open, execute_gnuplot_cmd, dbg
character(len=300) :: line1, line2, fmtline, cdummy
character(len=300) :: datafile, plotfile, imagefile, epsfile
character(len=7) :: color(111)
integer, external :: AixRm
integer :: Length
character(len=1023) :: realname_plt, realname_dat, realname_png, realname_eps, gnuplot_CMD
integer :: iErr, iSeed = 0
integer, external :: IsFreeUnit
real*8, external :: Random_Molcas

color(1) = '#ffffff'
color(2) = '#000000'
color(3) = '#a0a0a0'
color(4) = '#ff0000'
color(5) = '#00c000'
color(6) = '#0080ff'
color(7) = '#c000ff'
color(8) = '#00eeee'
color(9) = '#c04000'
color(10) = '#c8c800'
color(11) = '#4169e1'
color(12) = '#ffc020'
color(13) = '#008040'
color(14) = '#c080ff'
color(15) = '#306080'
color(16) = '#8b0000'
color(17) = '#408000'
color(18) = '#ff80ff'
color(19) = '#7fffd4'
color(20) = '#a52a2a'
color(21) = '#ffff00'
color(22) = '#40e0d0'
color(23) = '#000000'
color(24) = '#1a1a1a'
color(25) = '#333333'
color(26) = '#4d4d4d'
color(27) = '#666666'
color(28) = '#7f7f7f'
color(29) = '#999999'
color(30) = '#b3b3b3'
color(31) = '#c0c0c0'
color(32) = '#cccccc'
color(33) = '#e5e5e5'
color(34) = '#ffffff'
color(35) = '#f03232'
color(36) = '#90ee90'
color(37) = '#add8e6'
color(38) = '#f055f0'
color(39) = '#e0ffff'
color(40) = '#eedd82'
color(41) = '#ffb6c1'
color(42) = '#afeeee'
color(43) = '#ffd700'
color(44) = '#00ff00'
color(45) = '#006400'
color(46) = '#00ff7f'
color(47) = '#228b22'
color(48) = '#2e8b57'
color(49) = '#0000ff'
color(50) = '#00008b'
color(51) = '#191970'
color(52) = '#000080'
color(53) = '#0000cd'
color(54) = '#87ceeb'
color(55) = '#00ffff'
color(56) = '#ff00ff'
color(57) = '#00ced1'
color(58) = '#ff1493'
color(59) = '#ff7f50'
color(60) = '#f08080'
color(61) = '#ff4500'
color(62) = '#fa8072'
color(63) = '#e9967a'
color(64) = '#f0e68c'
color(65) = '#bdb76b'
color(66) = '#b8860b'
color(67) = '#f5f5dc'
color(68) = '#a08020'
color(69) = '#ffa500'
color(70) = '#ee82ee'
color(71) = '#9400d3'
color(72) = '#dda0dd'
color(73) = '#905040'
color(74) = '#556b2f'
color(75) = '#801400'
color(76) = '#801414'
color(77) = '#804014'
color(78) = '#804080'
color(79) = '#8060c0'
color(80) = '#8060ff'
color(81) = '#808000'
color(82) = '#ff8040'
color(83) = '#ffa040'
color(84) = '#ffa060'
color(85) = '#ffa070'
color(86) = '#ffc0c0'
color(87) = '#ffff80'
color(88) = '#ffffc0'
color(89) = '#cdb79e'
color(90) = '#f0fff0'
color(91) = '#a0b6cd'
color(92) = '#c1ffc1'
color(93) = '#cdc0b0'
color(94) = '#7cff40'
color(95) = '#a0ff20'
color(96) = '#bebebe'
color(97) = '#d3d3d3'
color(98) = '#d3d3d3'
color(99) = '#a0a0a0'
color(100) = '#a0b6cd'
color(101) = '#000000'
color(102) = '#1a1a1a'
color(103) = '#333333'
color(104) = '#4d4d4d'
color(105) = '#666666'
color(106) = '#7f7f7f'
color(107) = '#999999'
color(108) = '#b3b3b3'
color(109) = '#cccccc'
color(110) = '#e5e5e5'
color(111) = '#ffffff'

if (iSeed == 0) call GetSeed(iSeed)

dbg = .false.
iErr = 0
hmin = Zero
hmax = Zero
MHmin_calc = Zero
MHmax_calc = Zero
MHmin = Zero
MHmax = Zero
hmin = minval(H(:))-0.02_wp*maxval(H(:))
hmax = maxval(H(:))+0.02_wp*maxval(H(:))
MHmin_calc = minval(MHcalc)
MHmax_calc = maxval(MHcalc)
MHmin = MHmin_calc-0.01_wp*MHmax_calc
MHmax = MHmax_calc+0.01_wp*MHmax_calc

if (dbg) then
  write(u6,*) 'nH        = ',nH
  write(u6,*) 'hmin      = ',hmin
  write(u6,*) 'hmax      = ',hmax
  write(u6,*) 'MHmin_calc= ',MHmin_calc
  write(u6,*) 'MHmax_calc= ',MHmax_calc
  write(u6,*) 'MHmin     = ',MHmin
  write(u6,*) 'MHmax     = ',MHmax
  write(u6,*) 'zJ        = ',zJ
  do iH=1,nH
    write(u6,*) H(iH),(MHcalc(iH,iTempMagn),iTempMagn=1,nTempMagn)
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
gnuplot_version = Zero

! check if the file lineOUT exists
inquire(file='lineOUT',exist=file_exist,opened=is_file_open,number=file_number)

if (file_exist) then
  if (dbg) write(u6,'(A)') 'file "lineOUT" exists in WorkDir'
  if (is_file_open) then
    if (dbg) write(u6,'(A)') 'file "lineOUT" is opened'
    ! close the file:
    close(unit=file_number,status='DELETE')
  end if
  ! delete the file
  if (dbg) write(u6,'(A)') 'deleting the file...'
  iErr = AixRm('lineOUT')
  if (dbg) write(u6,*) 'iErr = ',iErr
else
  if (dbg) write(u6,'(A)') 'file "lineOUT" does not exist in WorkDir'
end if

! find the gnuplot
if (dbg) write(u6,'(A)') 'inquire which GNUPLOT'

!#ifdef __INTEL_COMPILER
call systemf('which gnuplot >> lineOUT',iErr)
if (dbg) write(u6,*) 'iErr = ',iErr
!#else
!call execute_command_line('which gnuplot >> lineOUT')
!#endif

inquire(file='lineOUT',exist=file_exist,opened=is_file_open,number=file_number,size=file_size)

if (dbg) then
  write(u6,*) 'File_number =',file_number
  write(u6,*) 'Is_file_open=',is_file_open
  write(u6,*) 'File_exist  =',file_exist
  write(u6,*) 'File_size   =',file_size
end if

if (file_exist) then
  if (file_size > 0) then
    if (dbg) write(u6,'(A)') 'new file "lineOUT" exists in WorkDir'

    file_number = IsFreeUnit(103)
    call molcas_open(file_number,'lineOUT')

    read(file_number,'(A)') line1

    if (dbg) then
      write(u6,*) 'line1=',line1
      write(u6,*) trim(line1)
    end if
    line2 = trim(line1)
    if (dbg) write(u6,*) 'line2=',line2

    close(file_number)
    if (dbg) write(u6,*) 'Closing lineOUT file'
    flush(u6)
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
if (dbg) write(u6,*) 'iErr = ',iErr
!-----------------------------------------------------------------------

! check the version of the gnuplot:
if (execute_gnuplot_cmd) then
  if (dbg) write(u6,'(A)') 'inquire which version of GNUPLOT is installed'
  ! attempt to execute the script
  write(gnuplot_CMD,'(2A)') trim(line2),' --version > lineOUT'
  if (dbg) write(u6,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
!# ifdef __INTEL_COMPILER
  call systemf(gnuplot_CMD,iErr)
  if (dbg) write(u6,*) 'iErr = ',iErr
!# else
!  call execute_command_line(gnuplot_CMD)
!# endif
  file_number = IsFreeUnit(102)
  call molcas_open(file_number,'lineOUT')
  read(file_number,*) cdummy,gnuplot_version
  if (dbg) write(u6,'(A,F4.1)') 'gnuplot_version = ',gnuplot_version
  if (abs(gnuplot_version) < 0.1_wp) execute_gnuplot_cmd = .false.
  close(file_number)
  ! remove file "lineOUT"
  iErr = AixRm('lineOUT')
  if (dbg) write(u6,*) 'iErr = ',iErr
end if

!-----------------------------------------------------------------------
! get the true real names of the files on disk:
write(datafile,'(3A)') 'MH.dat'
write(imagefile,'(3A)') 'MH.png'
write(epsfile,'(3A)') 'MH.eps'
write(plotfile,'(3A)') 'MH.plt'
call prgmtranslate(datafile,realname_dat,Length)
call prgmtranslate(imagefile,realname_png,Length)
call prgmtranslate(epsfile,realname_eps,Length)
call prgmtranslate(plotfile,realname_plt,Length)
if (dbg) then
  write(u6,'(3A)') 'realname_dat=',trim(realname_dat)
  write(u6,'(3A)') 'realname_png=',trim(realname_png)
  write(u6,'(3A)') 'realname_eps=',trim(realname_eps)
  write(u6,'(3A)') 'realname_plt=',trim(realname_plt)
end if

!-----------------------------------------------------------------------
! create the file "MH.dat"
inquire(file=datafile,exist=file_exist,opened=is_file_open,number=file_number)
if (file_exist) iErr = AixRm(trim(datafile))
if (dbg) write(u6,*) 'iErr = ',iErr
LuData = IsFreeUnit(104)
call molcas_open(LuData,datafile)
if (dbg) then
  write(u6,*) 'Opening "'//trim(datafile)//'" file'
  write(u6,*) 'Opening "'//trim(realname_dat)//'" file'
end if
write(fmtline,'(A,i0,A)') '(',nTempMagn+1,'ES24.14)'
do iH=1,nH
  write(LuData,fmtline) H(iH),(MHcalc(iH,iTempMagn),iTempMagn=1,nTempMagn)
end do
if (dbg) then
  write(u6,*) 'Writing into the "'//trim(datafile)//'"file'
  write(u6,*) 'Writing into the "'//trim(realname_dat)//'"file'
end if
close(LuData)
if (dbg) then
  write(u6,*) 'Closing the "'//trim(datafile)//'" file'
  write(u6,*) 'Closing the "'//trim(realname_dat)//'" file'
end if
flush(u6)

!-----------------------------------------------------------------------
! generate the GNUPLOT script in the $WorkDir
inquire(file=plotfile,exist=file_exist,opened=is_file_open,number=file_number)
if (file_exist) iErr = AixRm(trim(plotfile))
if (dbg) write(u6,*) 'iErr = ',iErr
LuPlt = IsFreeUnit(105)
call molcas_open(LuPlt,plotfile)
if (dbg) then
  write(u6,*) 'Opening "'//trim(plotfile)//'" file'
  write(u6,*) 'Opening "'//trim(realname_plt)//'" file'
end if

if (gnuplot_version < Five) then
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
  do iTempMagn=1,nTempMagn
    ik = iTempMagn+1
    r = Random_Molcas(iSeed)
    ic = floor(size(color)*r)+lbound(color,1)
    if ((iTempMagn == 1) .and. (nTempMagn > 1)) then
      write(LuPlt,'(A,i0,A,F7.3,A)') 'plot "'//trim(realname_dat)//'" using 1:',ik, &
        ' with lines lt 1  lw  8 lc rgb "'//color(ic)//'"  title "Calc. M at T=',TempMagn(iTempMagn),'K.", \'
    else if ((iTempMagn == 1) .and. (nTempMagn == 1)) then
      write(LuPlt,'(A,i0,A,F7.3,A)') 'plot "'//trim(realname_dat)//'" using 1:',ik, &
        ' with lines lt 1  lw  8 lc rgb "'//color(ic)//'"  title "Calc. M at T=',TempMagn(iTempMagn),'K."'
    end if
    if ((iTempMagn < nTempMagn) .and. (iTempMagn /= 1)) &
      write(LuPlt,'(A,i0,A,F7.3,A)') '     "'//trim(realname_dat)//'" using 1:',ik, &
        ' with lines lt 1  lw  8 lc rgb "'//color(ic)//'"  title "Calc. M at T=',TempMagn(iTempMagn),'K.", \'
    if ((iTempMagn == nTempMagn) .and. (nTempMagn > 1)) &
      write(LuPlt,'(A,i0,A,F7.3,A)') '     "'//trim(realname_dat)//'" using 1:',ik, &
        ' with lines lt 1  lw  8 lc rgb "'//color(ic)//'"  title "Calc. M at T=',TempMagn(iTempMagn),'K."'
  end do

else if ((gnuplot_version >= Five) .and. (gnuplot_version < Six)) then
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
  do iTempMagn=1,nTempMagn
    ik = iTempMagn+1
    r = Random_Molcas(iSeed)
    ic = floor(size(color)*r)+lbound(color,1)
    if ((iTempMagn == 1) .and. (nTempMagn > 1)) then
      write(LuPlt,'(A,i0,A,F7.3,A)') 'plot "'//trim(realname_dat)//'" using 1:',ik, &
        ' with lines lt 1  lw  8 lc rgb "'//color(ic)//'"  title "Calc. M at T=',TempMagn(iTempMagn),'K.", \'
    else if ((iTempMagn == 1) .and. (nTempMagn == 1)) then
      write(LuPlt,'(A,i0,A,F7.3,A)') 'plot "'//trim(realname_dat)//'" using 1:',ik, &
        ' with lines lt 1  lw  8 lc rgb "'//color(ic)//'"  title "Calc. M at T=',TempMagn(iTempMagn),'K."'
    end if
    if ((iTempMagn < nTempMagn) .and. (iTempMagn /= 1)) &
      write(LuPlt,'(A,i0,A,F7.3,A)') '     "'//trim(realname_dat)//'" using 1:',ik, &
        ' with lines lt 1  lw  8 lc rgb "'//color(ic)//'"  title "Calc. M at T=',TempMagn(iTempMagn),'K.", \'
    if ((iTempMagn == nTempMagn) .and. (nTempMagn > 1)) &
      write(LuPlt,'(A,i0,A,F7.3,A)') '     "'//trim(realname_dat)//'" using 1:',ik, &
        ' with lines lt 1  lw  8 lc rgb "'//color(ic)//'"  title "Calc. M at T=',TempMagn(iTempMagn),'K."'
  end do

else
  write(u6,*) 'GNUPLOT has version: ',gnuplot_version
  write(u6,*) 'This version of GNUPLOT is not known and thus, unsupported.'
  execute_gnuplot_cmd = .false.
end if
write(LuPlt,'(A)')
close(LuPlt)

if (execute_gnuplot_cmd) then
  ! attempt to execute the script
  if (dbg) then
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
  end if

  write(gnuplot_CMD,'(5A)') trim(line2),' < ',trim(realname_plt)
  if (dbg) write(u6,'(A,A)') 'gnuplot_CMD=',trim(gnuplot_CMD)

!# ifdef __INTEL_COMPILER
  call systemf(gnuplot_CMD,iErr)
  if (dbg) write(u6,*) 'iErr = ',iErr
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
#ifdef _WARNING_WORKAROUND_
if (.false.) call UNUSED_CHARACTER(cdummy)
#endif

end subroutine plot_MH_no_Exp
