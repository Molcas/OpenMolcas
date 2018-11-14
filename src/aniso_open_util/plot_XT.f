************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine plot_XT( nT, T, XTcalc, XTexp, zJ )
#ifdef NAGFOR
      use f90_unix_env
      use f90_unix_proc
#endif

      Implicit None

      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)        :: nT
      Real (kind=wp), intent(in) :: T(nT)
      Real (kind=wp), intent(in) :: XTexp(nT)
      Real (kind=wp), intent(in) :: XTcalc(nT)
      Real (kind=wp), intent(in) :: zJ
      ! local variables
      Real (kind=wp) :: tmin, tmax, XTmin_exp, XTmax_exp, XTmin_calc,
     &                  XTmax_calc, XTmin, XTmax

      Integer           :: file_number, istat, iT, LuPlt, LuData
      Integer, external :: IsFreeUnit
      Logical           :: file_exist, is_file_open
      Character(100)    :: line1, line2, lineOut
      Character(100)    :: gnuplot_CMD

      Logical :: dbg

      Call qEnter('plot_XT')
      dbg=.false.

      tmin=0.0_wp
      tmax=0.0_wp
      XTmin_exp=0.0_wp
      XTmax_exp=0.0_wp
      XTmin_calc=0.0_wp
      XTmax_calc=0.0_wp
      XTmin=0.0_wp
      XTmax=0.0_wp
      tmin=MINVAL(T)-0.03_wp*MAXVAL(T)
      tmax=MAXVAL(T)+0.03_wp*MAXVAL(T)
      XTmin_exp=MINVAL(XTexp)
      XTmax_exp=MAXVAL(XTexp)
      XTmin_calc=MINVAL(XTcalc)
      XTmax_calc=MAXVAL(XTcalc)
      XTmin=MIN(XTmin_exp,XTmin_calc)-0.05_wp*MAX(XTmax_exp,XTmax_calc)
      XTmax=MAX(XTmax_exp,XTmax_calc)+0.05_wp*MAX(XTmax_exp,XTmax_calc)


      If(dbg) Write(6,*) 'nT        = ',nT
      If(dbg) Write(6,*) 'tmin      = ',tmin
      If(dbg) Write(6,*) 'tmax      = ',tmax
      If(dbg) Write(6,*) 'XTmin_exp = ',XTmin_exp
      If(dbg) Write(6,*) 'XTmax_exp = ',XTmax_exp
      If(dbg) Write(6,*) 'XTmin_calc= ',XTmin_calc
      If(dbg) Write(6,*) 'XTmax_calc= ',XTmax_calc
      If(dbg) Write(6,*) 'XTmin     = ',XTmin
      If(dbg) Write(6,*) 'XTmax     = ',XTmax
      If(dbg) Write(6,*) 'zJ        = ',zJ

      ! generate the GNUPLOT script in the $WorkDir
      lineOut=' '
      line1=' '
      line2=' '
      gnuplot_CMD=' '
      istat=0
      file_exist=.false.
      is_file_open=.false.

      ! check if the file lineOUT exists
      Inquire(FILE="lineOUT",EXIST=file_exist,OPENED=is_file_open,
     &         NUMBER=file_number)

      If (file_exist) then
         If(dbg) Write(6,'(A)') 'file  "lineOUT" exists in WorkDir'
         If(is_file_open) then
           If(dbg) Write(6,'(A)') 'file  "lineOUT" is opened'
           ! close the file:
           CLOSE (UNIT=file_number,STATUS='DELETE')
         End If
         ! delete the file
         If(dbg) Write(6,'(A)') 'deleting the file...'

         Call system ( "rm -rf lineOUT" )

      Else
         If(dbg) Write(6,'(A)') 'file  "lineOUT" does not exist in '//
     &                           'WorkDir'
      End If


      ! find the gnuplot
      If(dbg) Write(6,'(A)') 'inquire which GNUPLOT'

      Call system ( "which gnuplot >> lineOUT" )

      Inquire(FILE="lineOUT",EXIST=file_exist,OPENED=is_file_open,
     &         NUMBER=file_number)

      If(dbg) Write(6,*) 'File_number =',file_number
      If(dbg) Write(6,*) 'Is_file_open=',is_file_open
      If(dbg) Write(6,*) 'File_exist  =',file_exist

      If (file_exist) then
         If(dbg) Write(6,'(A)') 'new file  "lineOUT" exists in WorkDir'

         file_number=IsFreeUnit(53)
         Call molcas_open(file_number,'lineOUT')

         READ (file_number,'(A)') line1

         If(dbg) Write(6,*) 'line1=',line1
         If(dbg) Write(6,*) trim(line1)
         line2=trim(line1)
         If(dbg) Write(6,*) 'line2=',line2

         CLOSE(file_number)
         If(dbg) Write(6,*) 'Closing lineOUT file'
         Call xFlush(6)
      Else
         Write(6,'(A)') 'file  "lineOUT" does not exist in WorkDir'
      End If
      ! remove file "lineOUT"
      Call system ( "rm -rf lineOUT" )

!  generate the file "XT.dat":
      LuData=IsFreeUnit(54)

      Call molcas_open(LuData,"XT.dat")
      If(dbg) Write(6,*) 'Opening "XT.dat" file'
      Do iT=1,nT
         Write(LuData,'(3ES24.14)') T(iT), XTexp(iT), XTcalc(iT)
      End Do
      If(dbg) Write(6,*) 'Writing into the "XT.dat" file'
      Close(LuData)
      If(dbg) Write(6,*) 'Closing the "XT.dat" file'
      Call xFlush(6)

      ! generate the GNUPLOT script in the $WorkDir
      LuPlt=IsFreeUnit(55)
      Call molcas_open(LuPlt,"XT.plt")
      If(dbg) Write(6,*) 'Opening "XT.plt" file'

      Write(LuPlt,'(A)') 'set terminal pngcairo transparent '//
     &                   'enhanced font "arial,10" fontscale '//
     &                   '4.0 size 1800, 1200'
      Write(LuPlt,'(A)') 'set output "XT.png"'
      Write(LuPlt,'(A)') 'set grid'
      Write(LuPlt,'(A)')
      Write(LuPlt,'(A)') '# Axes'
      Write(LuPlt,'(A)') '  set xlabel "Temperature / Kelvin" font '//
     &               '"Arial,16"'
      Write(LuPlt,'(A)') '  set ylabel "{/Symbol c}T / '//
     &               'cm^3 K mol^{-1}" font "Arial,16"'
      Write(LuPlt,'(A)')
      Write(LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[',
     &                     tmin,':',tmax,']'
      Write(LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',
     &                     XTmin,':',XTmax,']'
      Write(LuPlt,'(A)')
      Write(LuPlt,'(A)') '# Tics for axes'
      Write(LuPlt,'(A)') '  set xtics nomirror'
      Write(LuPlt,'(A)') '  set ytics nomirror'
      Write(LuPlt,'(A)') '  set mxtics'
      Write(LuPlt,'(A)') '  set mytics'
      Write(LuPlt,'(A)')
      Write(LuPlt,'(A)') '# Margins'
      Write(LuPlt,'(A)') '  set bmargin at screen 0.15'
      Write(LuPlt,'(A)') '  set lmargin at screen 0.12'
      Write(LuPlt,'(A)')
      Write(LuPlt,'(A)') '# legend'
      Write(LuPlt,'(A)') '  set key'
      Write(LuPlt,'(A)')
      Write(LuPlt,'(A)') '# actual plotting'
      Write(LuPlt,'(A)') 'plot "XT.dat" using 1:2:(2) '//
     &                   'with circles lt 1  lw 3 lc rgb "black" '//
     &                   ' title "Experiment", \'
      Write(LuPlt,'(A)') '     "XT.dat" using 1:3 with'//
     &                   ' lines       lt 1  lw 10 lc rgb "red" '//
     &                   ' title "Calculation"'
      Write(LuPlt,'(A)')

      ! attempt to execute the script
      Write(gnuplot_CMD,'(A,A)') trim(line2),'  XT.plt'
      If(dbg) Write(6,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
      Close(LuPlt)

      Call system ( gnuplot_CMD )

      Call qExit('plot_XT')
      Return
      End Subroutine plot_XT
