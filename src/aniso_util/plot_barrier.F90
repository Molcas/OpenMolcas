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

subroutine plot_barrier(nBlock,nMult,nDIM,E,M)

use Constants, only: Zero, One, Three, Five, Six, Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBlock, nMult, nDIM(nMult)
real(kind=wp), intent(in) :: E(nBlock)
complex(kind=wp), intent(in) :: M(3,nMult,10,nMult,10)
real(kind=wp) :: dlt, emax, emin, F, fact, RAVE, RAVEMAX, RAVEMIN, X, xend, xmax, xmin, xstart, Y, yend, ystart, Z, gnuplot_version
integer(kind=iwp) :: file_number, file_size, i, i1, iErr, j, j1, k, l, Length, LuData, LuPlt
logical(kind=iwp) :: execute_gnuplot_cmd, file_exist, is_file_open
character(len=1023) :: gnuplot_CMD, realname_e_dat, realname_eps, realname_m_dat, realname_plt, realname_png
character(len=100) :: cdummy, datafile_e, datafile_m, epsfile, fmtx, imagefile, line1, line2, plotfile
integer(kind=iwp), external :: AixRm, IsFreeUnit

#include "macros.fh"

iErr = 0

xmin = minval(real(M(3,:,:,:,:),kind=wp))-0.1_wp*maxval(real(M(3,:,:,:,:),kind=wp))
xmax = maxval(real(M(3,:,:,:,:),kind=wp))+0.1_wp*maxval(real(M(3,:,:,:,:),kind=wp))
dlt = 0.01_wp*(xmax-xmin)

F = maxval(E(:))
fact = One
if (F < Ten) then
  fact = One
!else if ((F >= Ten) .and. (F < 100.0_wp)) then
!  fact = Ten
!else if ((F >= 1.0e2_wp) .and. (F < 1.0e3_wp)) then
!  fact = 1.0e2_wp
!else if ((F >= 1.0e3_wp) .and. (F < 1.0e4_wp)) then
!  fact = 1.0e3_wp
else if ((F >= 1.0e4_wp) .and. (F < 1.0e5_wp)) then
  fact = 1.0e4_wp
else if ((F >= 1.0e5_wp) .and. (F < 1.0e6_wp)) then
  fact = 1.0e5_wp
end if
Emin = minval(E(:))/fact-0.05_wp*maxval(E(:))/fact
Emax = maxval(E(:))/fact+0.05_wp*maxval(E(:))/fact

RAVEMIN = 9999999999999999.0_wp
RAVEMAX = -9999999999999999.0_wp

l = 0
do i=1,nMult
  do i1=1,ndim(i)
    l = l+1
    k = 0
    do j=1,nMult
      do j1=1,ndim(j)
        k = k+1
        if ((i == j) .and. (i1 == j1)) cycle
        X = abs(M(1,i,i1,j,j1))
        Y = abs(M(2,i,i1,j,j1))
        Z = abs(M(3,i,i1,j,j1))
        RAVE = (X+Y+Z)/Three
        if (RAVE < RAVEMIN) RAVEMIN = RAVE
        if (RAVE > RAVEMAX) RAVEMAX = RAVE
      end do
    end do
  end do
end do

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

call systemf('which gnuplot >> lineOUT',iErr)
#ifdef _DEBUGPRINT_
write(u6,*) 'iErr = ',iErr
#endif

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

    file_number = IsFreeUnit(45)
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

!!!!! check the version of the gnuplot:
if (execute_gnuplot_cmd) then
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'inquire which version of GNUPLOT is installed'
# endif
  ! attempt to execute the script
  write(gnuplot_CMD,'(2A)') trim(line2),' --version > lineOUT'
# ifdef _DEBUGPRINT_
  write(u6,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
# endif
  call systemf(gnuplot_CMD,iErr)
# ifdef _DEBUGPRINT_
  write(u6,*) 'iErr = ',iErr
# endif
  file_number = IsFreeUnit(42)
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
write(datafile_e,'(3A)') 'BARRIER_ENE.dat'
write(datafile_m,'(3A)') 'BARRIER_TME.dat'
write(imagefile,'(3A)') 'BARRIER.png'
write(epsfile,'(3A)') 'BARRIER.eps'
write(plotfile,'(3A)') 'BARRIER.plt'
call prgmtranslate(datafile_e,realname_e_dat,Length)
call prgmtranslate(datafile_m,realname_m_dat,Length)
call prgmtranslate(imagefile,realname_png,Length)
call prgmtranslate(epsfile,realname_eps,Length)
call prgmtranslate(plotfile,realname_plt,Length)
#ifdef _DEBUGPRINT_
write(u6,'(3A)') 'realname_e_dat=',trim(realname_e_dat)
write(u6,'(3A)') 'realname_m_dat=',trim(realname_m_dat)
write(u6,'(3A)') 'realname_png=',trim(realname_png)
write(u6,'(3A)') 'realname_eps=',trim(realname_eps)
write(u6,'(3A)') 'realname_plt=',trim(realname_plt)
#endif

!!!!! prepare the energy data file:
inquire(file=datafile_e,exist=file_exist,opened=is_file_open,number=file_number)
if (file_exist) iErr = AixRm(trim(datafile_e))
#ifdef _DEBUGPRINT_
write(u6,*) 'iErr = ',iErr
#endif
LuData = IsFreeUnit(75)
call molcas_open(LuData,datafile_e)
#ifdef _DEBUGPRINT_
write(u6,*) 'Opening "'//trim(datafile_e)//'" file'
write(u6,*) 'Opening "'//trim(realname_e_dat)//'" file'
call xFlush(u6)
#endif
! write energies
l = 0
do i=1,nMult
  do i1=1,ndim(i)
    l = l+1
    xstart = real(M(3,i,i1,i,i1),kind=wp)-dlt
    xend = real(M(3,i,i1,i,i1),kind=wp)+dlt

    write(LuData,'(3ES24.14)') xstart,E(l)/fact
    write(LuData,'(3ES24.14)') xend,E(l)/fact
    write(LuData,'(A)') ' '
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*) 'Writing into the "'//trim(datafile_e)//'" file'
write(u6,*) 'Writing into the "'//trim(realname_e_dat)//'" file'
#endif
close(LuData)
#ifdef _DEBUGPRINT_
write(u6,*) 'Closing the "'//trim(datafile_e)//'" file'
write(u6,*) 'Closing the "'//trim(realname_e_dat)//'" file'
call xFlush(u6)
#endif

!!!!!===================================================================
!!!!! prepare the transition matrix elements file:
inquire(file=datafile_m,exist=file_exist,opened=is_file_open,number=file_number)
if (file_exist) iErr = AixRm(trim(datafile_m))
#ifdef _DEBUGPRINT_
write(u6,*) 'iErr = ',iErr
#endif
LuData = IsFreeUnit(76)
call molcas_open(LuData,datafile_m)
#ifdef _DEBUGPRINT_
write(u6,*) 'Opening "'//trim(datafile_m)//'" file'
write(u6,*) 'Opening "'//trim(realname_m_dat)//'" file'
call xFlush(u6)
#endif
! write energies,
l = 0
do i=1,nMult
  do i1=1,ndim(i)
    l = l+1
    k = 0
    do j=1,nMult
      do j1=1,ndim(j)
        k = k+1

        if ((i == j) .and. (i1 == j1)) cycle

        X = abs(M(1,i,i1,j,j1))
        Y = abs(M(2,i,i1,j,j1))
        Z = abs(M(3,i,i1,j,j1))
        RAVE = (X+Y+Z)/Three

        xstart = real(M(3,i,i1,i,i1),kind=wp)
        xend = real(M(3,j,j1,j,j1),kind=wp)

        ystart = E(l)/fact
        yend = E(k)/fact

        fmtx = '(ES24.14,ES24.14,ES24.14)'
        if (RAVE > RAVEMAX*5.0e-3_wp) then
          write(LuData,fmtx) xstart,ystart,RAVE
          write(LuData,fmtx) xend,yend,RAVE
          write(LuData,'(A)') ' '
        end if
      end do
    end do
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*) 'Writing into the "'//trim(datafile_m)//'" file'
#endif
close(LuData)
#ifdef _DEBUGPRINT_
write(u6,*) 'Closing the "'//trim(datafile_m)//'" file'
call xFlush(u6)
#endif

!!!!!===================================================================
!!!!! generate the GNUPLOT script in the $WorkDir

inquire(file=plotfile,exist=file_exist,opened=is_file_open,number=file_number)
if (file_exist) iErr = AixRm(trim(plotfile))
#ifdef _DEBUGPRINT_
write(u6,*) 'iErr = ',iErr
#endif
LuPlt = IsFreeUnit(85)
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
  write(LuPlt,'(A)') '  set xlabel "Momentum / {/Symbol m}_{B}" font "Arial,12"'
  if (nint(fact) == 1) then
    write(LuPlt,'(A)') '  set ylabel "Energy / cm^{-1}" font "Arial,12"'
  else if (nint(fact) == 10) then
    write(LuPlt,'(A)') '  set ylabel "Energy (x10) / cm^{-1}" font "Arial,12"'
  else if (nint(fact) == 100) then
    write(LuPlt,'(A)') '  set ylabel "Energy (x100) / cm^{-1}" font "Arial,12"'
  else if (nint(fact) == 1000) then
    write(LuPlt,'(A)') '  set ylabel "Energy (x1000)/ cm^{-1}" font "Arial,12"'
  else if (nint(fact) == 10000) then
    write(LuPlt,'(A)') '  set ylabel "Energy (x10000)/ cm^{-1}" font "Arial,12"'
  else if (nint(fact) == 100000) then
    write(LuPlt,'(A)') '  set ylabel "Energy (x100000) / cm^{-1}" font "Arial,12"'
  end if
  write(LuPlt,'(A)')
  write(LuPlt,'(A,F24.14,A,F24.14,A)') '  set xrange[',xmin,':',xmax,']'
  write(LuPlt,'(A,F24.14,A,F24.14,A)') '  set yrange[',Emin,':',Emax,']'
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
  write(LuPlt,'(A)') '  set key off' !right bottom'
  write(LuPlt,'(A)')

  write(LuPlt,'(A)')
  write(LuPlt,'(A,ES24.14,A,ES24.14,A)') 'set palette defined ( ',RAVEMIN," '#FFF0F0', ",RAVEMAX," 'red' )"

  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# actual plotting'
  write(LuPlt,'(A)') 'plot  "'//trim(realname_m_dat)//'" using 1:2:3 with lines lt 1  lw 1 lc palette , \'
  write(LuPlt,'(A)') '      "'//trim(realname_e_dat)//'" using 1:2   with lines lt 1  lw 8 lc rgb "black" '
  write(LuPlt,'(A)')

else if ((gnuplot_version >= Five) .and. (gnuplot_version < Six)) then
  !===  GNUPLOT VERSION 5 and above ==>>  generate PNG
  write(LuPlt,'(A)') 'set terminal pngcairo transparent enhanced font "arial,10" fontscale 4.0 size 1800, 1200'
  write(LuPlt,'(A)') 'set output "'//trim(realname_png)//'" '
  write(LuPlt,'(A)') 'set grid'
  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# Axes'
  write(LuPlt,'(A)') '  set xlabel "Momentum / {/Symbol m}_{B}" font "Arial,12"'
  if (nint(fact) == 1) then
    write(LuPlt,'(A)') '  set ylabel "Energy / cm^{-1}" font "Arial,12"'
  else if (nint(fact) == 10) then
    write(LuPlt,'(A)') '  set ylabel "Energy (x10) / cm^{-1}" font "Arial,12"'
  else if (nint(fact) == 100) then
    write(LuPlt,'(A)') '  set ylabel "Energy (x100) / cm^{-1}" font "Arial,12"'
  else if (nint(fact) == 1000) then
    write(LuPlt,'(A)') '  set ylabel "Energy (x1000)/ cm^{-1}" font "Arial,12"'
  else if (nint(fact) == 10000) then
    write(LuPlt,'(A)') '  set ylabel "Energy (x10000)/ cm^{-1}" font "Arial,12"'
  else if (nint(fact) == 100000) then
    write(LuPlt,'(A)') '  set ylabel "Energy (x100000) / cm^{-1}" font "Arial,12"'
  end if
  write(LuPlt,'(A)')
  write(LuPlt,'(A,F24.14,A,F24.14,A)') '  set xrange[',xmin,':',xmax,']'
  write(LuPlt,'(A,F24.14,A,F24.14,A)') '  set yrange[',Emin,':',Emax,']'
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
  write(LuPlt,'(A)') '  set key off' !right bottom'
  write(LuPlt,'(A)')

  write(LuPlt,'(A)')
  write(LuPlt,'(A,ES24.14,A,ES24.14,A)') 'set palette defined ( ',RAVEMIN," '#FFF0F0', ",RAVEMAX," 'red' )"

  write(LuPlt,'(A)')
  write(LuPlt,'(A)') '# actual plotting'
  write(LuPlt,'(A)') 'plot  "'//trim(realname_m_dat)//'" using 1:2:3 with lines lt 1  lw 1 lc palette , \'
  write(LuPlt,'(A)') '      "'//trim(realname_e_dat)//'" using 1:2   with lines lt 1  lw 8 lc rgb "black" '
  write(LuPlt,'(A)')
else
  write(u6,*) 'GNUPLOT has version: ',gnuplot_version
  write(u6,*) 'This version of GNUPLOT is not known and thus, unsupported.'
  execute_gnuplot_cmd = .false.
end if
close(LuPlt)

if (execute_gnuplot_cmd) then
  ! attempt to execute the script
  write(gnuplot_CMD,'(5A)') trim(line2),' < ',trim(realname_plt)
# ifdef _DEBUGPRINT_
  write(u6,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
# endif
  call systemf(gnuplot_CMD,iErr)
# ifdef _DEBUGPRINT_
  write(u6,*) 'iErr = ',iErr
# endif

  if (gnuplot_version < Five) then
    inquire(file=trim(realname_eps),exist=file_exist,opened=is_file_open,number=file_number)
    if (file_exist) then
      write(u6,'(A,i0,A)') 'File "'//trim(realname_eps)//'" was created in Working directory.'
    else
      write(u6,'(A,i0,A)') 'File "'//trim(realname_eps)//'" was NOT created in Working directory.'
    end if
  else
    inquire(file=trim(realname_png),exist=file_exist,opened=is_file_open,number=file_number)
    if (file_exist) then
      write(u6,'(A,i0,A)') 'File "'//trim(realname_png)//'" was created in Working directory.'
    else
      write(u6,'(A,i0,A)') 'File "'//trim(realname_png)//'" was NOT created in Working directory.'
    end if
  end if
end if

return

end subroutine plot_barrier
