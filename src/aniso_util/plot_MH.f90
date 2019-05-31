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
  Subroutine plot_MH_with_Exp( nH, H, nTempMagn, TempMagn, MHcalc, MHexp, zJ )

#ifdef NAGFOR
      use f90_unix_env
      use f90_unix_proc
#endif

  IMPLICIT NONE

  Integer, parameter    :: wp=SELECTED_REAL_KIND(p=15,r=307)
  INTEGER, INTENT(in)   :: nH, nTempMagn
  REAL (wp), INTENT(in) :: H(nH), TempMagn(nTempMagn)
  REAL (wp), INTENT(in) :: MHexp(nH,nTempMagn)
  REAL (wp), INTENT(in) :: MHcalc(nH,nTempMagn)
  REAL (wp), INTENT(in) :: zJ
  ! local variables
  REAL (wp) :: hmin, hmax, MHmin_exp, MHmax_exp, MHmin_calc, MHmax_calc, MHmin, MHmax
  REAL (wp) :: gnuplot_version
  INTEGER           :: file_number, istat, iH, iTempMagn, LuPlt, LuData, file_size, StdOut, iErr
  LOGICAL           :: file_exist, is_file_open, execute_gnuplot_cmd, dbg
  CHARACTER(LEN=100):: line1, line2, lineOut, cdummy
  CHARACTER(LEN=100):: gnuplot_CMD, filedat, fileplt
  CHARACTER(LEN=7) :: color(111)
  INTEGER, EXTERNAL:: AixRm

  color(  1)="#ffffff"; color(  2)="#000000"; color(  3)="#a0a0a0"; color(  4)="#ff0000"; color(  5)="#00c000"
  color(  6)="#0080ff"; color(  7)="#c000ff"; color(  8)="#00eeee"; color(  9)="#c04000"; color( 10)="#c8c800"
  color( 11)="#4169e1"; color( 12)="#ffc020"; color( 13)="#008040"; color( 14)="#c080ff"; color( 15)="#306080"
  color( 16)="#8b0000"; color( 17)="#408000"; color( 18)="#ff80ff"; color( 19)="#7fffd4"; color( 20)="#a52a2a"
  color( 21)="#ffff00"; color( 22)="#40e0d0"; color( 23)="#000000"; color( 24)="#1a1a1a"; color( 25)="#333333"
  color( 26)="#4d4d4d"; color( 27)="#666666"; color( 28)="#7f7f7f"; color( 29)="#999999"; color( 30)="#b3b3b3"
  color( 31)="#c0c0c0"; color( 32)="#cccccc"; color( 33)="#e5e5e5"; color( 34)="#ffffff"; color( 35)="#f03232"
  color( 36)="#90ee90"; color( 37)="#add8e6"; color( 38)="#f055f0"; color( 39)="#e0ffff"; color( 40)="#eedd82"
  color( 41)="#ffb6c1"; color( 42)="#afeeee"; color( 43)="#ffd700"; color( 44)="#00ff00"; color( 45)="#006400"
  color( 46)="#00ff7f"; color( 47)="#228b22"; color( 48)="#2e8b57"; color( 49)="#0000ff"; color( 50)="#00008b"
  color( 51)="#191970"; color( 52)="#000080"; color( 53)="#0000cd"; color( 54)="#87ceeb"; color( 55)="#00ffff"
  color( 56)="#ff00ff"; color( 57)="#00ced1"; color( 58)="#ff1493"; color( 59)="#ff7f50"; color( 60)="#f08080"
  color( 61)="#ff4500"; color( 62)="#fa8072"; color( 63)="#e9967a"; color( 64)="#f0e68c"; color( 65)="#bdb76b"
  color( 66)="#b8860b"; color( 67)="#f5f5dc"; color( 68)="#a08020"; color( 69)="#ffa500"; color( 70)="#ee82ee"
  color( 71)="#9400d3"; color( 72)="#dda0dd"; color( 73)="#905040"; color( 74)="#556b2f"; color( 75)="#801400"
  color( 76)="#801414"; color( 77)="#804014"; color( 78)="#804080"; color( 79)="#8060c0"; color( 80)="#8060ff"
  color( 81)="#808000"; color( 82)="#ff8040"; color( 83)="#ffa040"; color( 84)="#ffa060"; color( 85)="#ffa070"
  color( 86)="#ffc0c0"; color( 87)="#ffff80"; color( 88)="#ffffc0"; color( 89)="#cdb79e"; color( 90)="#f0fff0"
  color( 91)="#a0b6cd"; color( 92)="#c1ffc1"; color( 93)="#cdc0b0"; color( 94)="#7cff40"; color( 95)="#a0ff20"
  color( 96)="#bebebe"; color( 97)="#d3d3d3"; color( 98)="#d3d3d3"; color( 99)="#a0a0a0"; color(100)="#a0b6cd"
  color(101)="#000000"; color(102)="#1a1a1a"; color(103)="#333333"; color(104)="#4d4d4d"; color(105)="#666666"
  color(106)="#7f7f7f"; color(107)="#999999"; color(108)="#b3b3b3"; color(109)="#cccccc"; color(110)="#e5e5e5"
  color(111)="#ffffff"

  dbg=.false.
  iErr=0
  StdOut=6
  hmin=0.0_wp
  hmax=0.0_wp
  MHmin_exp=0.0_wp
  MHmax_exp=0.0_wp
  MHmin_calc=0.0_wp
  MHmax_calc=0.0_wp
  MHmin=0.0_wp
  MHmax=0.0_wp
  hmin=MINVAL(H)-0.05_wp*MAXVAL(H)
  hmax=MAXVAL(H)+0.05_wp*MAXVAL(H)
  MHmin_exp=MINVAL(MHexp)
  MHmax_exp=MAXVAL(MHexp)
  MHmin_calc=MINVAL(MHcalc)
  MHmax_calc=MAXVAL(MHcalc)
  MHmin=MIN(MHmin_exp,MHmin_calc)-0.08_wp*MAX(MHmax_exp,MHmax_calc)
  MHmax=MAX(MHmax_exp,MHmax_calc)+0.08_wp*MAX(MHmax_exp,MHmax_calc)


  IF (dbg) WRITE (StdOut,*) 'nH        = ',nH
  IF (dbg) WRITE (StdOut,*) 'hmin      = ',hmin
  IF (dbg) WRITE (StdOut,*) 'hmax      = ',hmax
  IF (dbg) WRITE (StdOut,*) 'MHmin_exp = ',MHmin_exp
  IF (dbg) WRITE (StdOut,*) 'MHmax_exp = ',MHmax_exp
  IF (dbg) WRITE (StdOut,*) 'MHmin_calc= ',MHmin_calc
  IF (dbg) WRITE (StdOut,*) 'MHmax_calc= ',MHmax_calc
  IF (dbg) WRITE (StdOut,*) 'MHmin     = ',MHmin
  IF (dbg) WRITE (StdOut,*) 'MHmax     = ',MHmax
  IF (dbg) WRITE (StdOut,*) 'zJ        = ',zJ
  IF (dbg) THEN
    DO iH=1,nH
      WRITE(StdOut,*) H(iH), ( MHcalc(iH,iTempMagn), MHexp(iH,iTempMagn), iTempMagn=1,nTempMagn)
    END DO
  END IF

  ! generate the GNUPLOT script in the $WorkDir
  lineOut=' '
  line1=' '
  line2=' '
  gnuplot_CMD=' '
  istat=0
  file_exist=.false.
  is_file_open=.false.
  file_size=0
  execute_gnuplot_cmd =.false.
  gnuplot_version=0.0_wp

  ! check if the file lineOUT exists
  INQUIRE(FILE="lineOUT",EXIST=file_exist,OPENED=is_file_open, NUMBER=file_number)

  IF (file_exist) then
     IF (dbg) WRITE (StdOut,'(A)') 'file  "lineOUT" exists in WorkDir'
     IF (is_file_open) then
       IF (dbg) WRITE (StdOut,'(A)') 'file  "lineOUT" is opened'
       ! close the file:
       CLOSE (UNIT=file_number,STATUS='DELETE')
     END IF
     ! delete the file
     IF (dbg) WRITE (StdOut,'(A)') 'deleting the file...'
     iErr=AixRm("lineOUT")
  ELSE
     IF (dbg) WRITE (StdOut,'(A)') 'file  "lineOUT" does not exist in WorkDir'
  END IF


  ! find the gnuplot
  IF (dbg) WRITE (StdOut,'(A)') 'inquire which GNUPLOT'

  CALL systemf ( "which gnuplot >> lineOUT", iErr )

  INQUIRE(FILE="lineOUT",EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number,SIZE=file_size)

  IF (dbg) WRITE (StdOut,*) 'File_number =',file_number
  IF (dbg) WRITE (StdOut,*) 'Is_file_open=',is_file_open
  IF (dbg) WRITE (StdOut,*) 'File_exist  =',file_exist
  IF (dbg) WRITE (StdOut,*) 'File_size   =',file_size

  IF (file_exist) then
    IF (file_size>0) then
      IF (dbg) WRITE (StdOut,'(A)') 'new file  "lineOUT"  exists in WorkDir'

      file_number=453
      Call molcas_open(file_number,'lineOUT')

      READ (file_number,'(A)') line1

      IF (dbg) WRITE (StdOut,*) 'line1=',line1
      IF (dbg) WRITE (StdOut,*) trim(line1)
      line2=trim(line1)
      IF (dbg) WRITE (StdOut,*) 'line2=',line2

      CLOSE(file_number)
      IF (dbg) WRITE (StdOut,*) 'Closing lineOUT file'
      execute_gnuplot_cmd =.true.
    ELSE
      ! file_size =0
      WRITE (StdOut,'(A)') 'file  "lineOUT" has a size=0. gnuplot was not found on the system.'
      WRITE (StdOut,'(A)') 'plots will not be created.'
    END IF
  ELSE
      WRITE (StdOut,'(A)') 'file  "lineOUT" does not exist in WorkDir'
  END IF
  ! remove file "lineOUT"
  iErr=AixRm("lineOUT")
!--------------------------------------------------------------------------------------------


! check the version of the gnuplot:
  IF ( execute_gnuplot_cmd ) Then
    IF (dbg) WRITE (StdOut,'(A)') 'inquire which version of GNUPLOT is installed'
    ! attempt to execute the script
    WRITE (gnuplot_CMD,'(2A)') trim(line2),' --version > lineOUT'
    IF (dbg) WRITE (StdOut,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
    CALL systemf ( gnuplot_CMD, iErr )
    file_number=452
    Call molcas_open(file_number,'lineOUT')
    READ (file_number,*) cdummy, gnuplot_version
    IF (dbg) WRITE (StdOut,'(A,F4.1)') 'gnuplot_version = ', gnuplot_version
    IF (abs(gnuplot_version)<0.1_wp) execute_gnuplot_cmd =.false.
    CLOSE (file_number)
    ! remove file "lineOUT"
    iErr=AixRm("lineOUT")
  END IF
!--------------------------------------------------------------------------------------------




  DO iTempMagn=1,nTempMagn
     ! generate the file "MH.dat":
     WRITE(filedat,'(A,I0,A)') 'MH_T_',iTempMagn,'.dat'
     INQUIRE(FILE=filedat,EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
     IF(file_exist)  iErr=AixRm(trim(filedat))
     LuData=554+iTempMagn
     Call molcas_open(LuData,filedat)
     IF (dbg) WRITE (StdOut,*) 'Opening'//trim(filedat)//' file'
     DO iH=1,nH
        WRITE (LuData,'(3ES24.14)') H(iH), MHexp(iH,iTempMagn), MHcalc(iH,iTempMagn)
     END DO
     IF (dbg) WRITE (StdOut,*) 'Writing into the '//trim(filedat)//' file. iTempMagn=',iTempMagn
     CLOSE (LuData)
     IF (dbg) WRITE (StdOut,*) 'Closing the '//trim(filedat)//' file. iTempMagn=',iTempMagn
     FLUSH (StdOut)



     ! generate the GNUPLOT script in the $WorkDir
     Write(fileplt,'(A,I0,A)') 'MH_T_',iTempMagn,'.plt'
     INQUIRE(FILE=fileplt,EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
     IF(file_exist)  iErr=AixRm(trim(fileplt))
     LuPlt=455+iTempMagn
     Call molcas_open(LuPlt,fileplt)
     IF (dbg) WRITE (StdOut,*) 'Opening '//trim(fileplt)//' file'


     IF ( gnuplot_version < 5.0_wp ) Then
     !===  GNUPLOT VERSION 4 and below ==>>  generate EPS
        WRITE (LuPlt,'(A)') 'set terminal postscript eps enhanced color  size 3.0, 2.0 font "arial, 10"'
        WRITE (LuPlt,'(A,i0,A)') 'set output "MH_T_',iTempMagn,'.eps"'
        WRITE (LuPlt,'(A)') 'set grid'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A)') '# Axes'
        WRITE (LuPlt,'(A)') '  set xlabel "Magnetic Field / Tesla" font "Arial,14"'
        WRITE (LuPlt,'(A)') '  set ylabel "M / {/Symbol m}_{B}" font "Arial,14"'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[', hmin,':', hmax,']'
        WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',MHmin,':',MHmax,']'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A)') '# Tics for axes'
        WRITE (LuPlt,'(A)') '  set xtics nomirror'
        WRITE (LuPlt,'(A)') '  set ytics nomirror'
        WRITE (LuPlt,'(A)') '  set mxtics'
        WRITE (LuPlt,'(A)') '  set mytics'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A)') '# Margins'
        WRITE (LuPlt,'(A)') '  set bmargin at screen 0.20'
        WRITE (LuPlt,'(A)') '  set lmargin at screen 0.15'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A)') '# legend'
        WRITE (LuPlt,'(A)') '  set key right bottom'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A)') '# actual plotting'

        WRITE (LuPlt,'(A,F7.3,A)') 'plot "'//trim(filedat)//'" using 1:2 with points  lt 1  lw 3 lc rgb "black"  title "Exp. T=',&
                                    TempMagn(iTempMagn),' K.", \'
        WRITE (LuPlt,'(A,F7.3,A)') '     "'//trim(filedat)//'" using 1:3 with lines   lt 1  lw 8 lc rgb "red"    title "Calc.T=',&
                                    TempMagn(iTempMagn),' K."'
        WRITE (LuPlt,'(A)')

     ELSE IF ( (gnuplot_version >= 5.0_wp) .AND.(gnuplot_version < 6.0_wp) ) Then
     !===  GNUPLOT VERSION 5 and above ==>>  generate PNG
        WRITE (LuPlt,'(A)') 'set terminal pngcairo transparent enhanced font "arial,10" fontscale 4.0 size 1800, 1200'
        WRITE (LuPlt,'(A,i0,A)') 'set output "MH_T_',iTempMagn,'.png"'
        WRITE (LuPlt,'(A)') 'set grid'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A)') '# Axes'
        WRITE (LuPlt,'(A)') '  set xlabel "Magnetic Field / Tesla" font "Arial,14"'
        WRITE (LuPlt,'(A)') '  set ylabel "M / {/Symbol m}_{B}" font "Arial,14"'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[', hmin,':', hmax,']'
        WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',MHmin,':',MHmax,']'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A)') '# Tics for axes'
        WRITE (LuPlt,'(A)') '  set xtics nomirror'
        WRITE (LuPlt,'(A)') '  set ytics nomirror'
        WRITE (LuPlt,'(A)') '  set mxtics'
        WRITE (LuPlt,'(A)') '  set mytics'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A)') '# Margins'
        WRITE (LuPlt,'(A)') '  set bmargin at screen 0.20'
        WRITE (LuPlt,'(A)') '  set lmargin at screen 0.15'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A)') '# legend'
        WRITE (LuPlt,'(A)') '  set key right bottom'
        WRITE (LuPlt,'(A)')
        WRITE (LuPlt,'(A)') '# actual plotting'

        WRITE (LuPlt,'(A,F7.3,A)') 'plot "'//trim(filedat)//'" using 1:2 with circles lt 1  lw 3 lc rgb "black"  title "Exp. T=',&
                                    TempMagn(iTempMagn),' K.", \'
        WRITE (LuPlt,'(A,F7.3,A)') '     "'//trim(filedat)//'" using 1:3 with lines   lt 1  lw 8 lc rgb "red"    title "Calc.T=',&
                                    TempMagn(iTempMagn),' K."'
        WRITE (LuPlt,'(A)')
     ELSE
       WRITE (StdOut,*) 'GNUPLOT has version: ', gnuplot_version
       WRITE (StdOut,*) 'This version of GNUPLOT is not known and thus, unsupported.'
       execute_gnuplot_cmd=.false.
     END IF
     CLOSE (LuPlt)

     IF (execute_gnuplot_cmd) Then
       ! attempt to execute the script
       WRITE (gnuplot_CMD,'(5A)') trim(line2),'  ',trim(fileplt)
       IF (dbg) WRITE (StdOut,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
       CALL systemf ( gnuplot_CMD, iErr )
       IF ( gnuplot_version < 5.0_wp ) Then
         WRITE (StdOut,'(A,i0,A)') 'File "MH_T_',iTempMagn,'.eps" was created in Working directory.'
       ELSE
         WRITE (StdOut,'(A,i0,A)') 'File "MH_T_',iTempMagn,'.png" was created in Working directory.'
       END IF
     END IF
  END DO ! iTempMagn

  RETURN
  END SUBROUTINE plot_MH_with_Exp





  Subroutine plot_MH_no_Exp( nH, H, nTempMagn, TempMagn, MHcalc, zJ )
#ifdef NAGFOR
      use f90_unix_env
      use f90_unix_proc
#endif

  IMPLICIT NONE

  Integer, parameter    :: wp=SELECTED_REAL_KIND(p=15,r=307)
  INTEGER, INTENT(in)   :: nH, nTempMagn
  REAL (wp), INTENT(in) :: H(nH), TempMagn(nTempMagn)
  REAL (wp), INTENT(in) :: MHcalc(nH,nTempMagn)
  REAL (wp), INTENT(in) :: zJ
  ! local variables
  REAL (wp) :: hmin, hmax, MHmin_calc, MHmax_calc, MHmin, MHmax, r
  REAL (wp) :: gnuplot_version
  INTEGER           :: file_number, istat, iH, iTempMagn, LuPlt, LuData, ik, ic, file_size, StdOut, iErr
  LOGICAL           :: file_exist, is_file_open, execute_gnuplot_cmd, dbg
  CHARACTER(LEN=100):: line1, line2, lineOut, fmtline, cdummy
  CHARACTER(LEN=100):: gnuplot_CMD
  CHARACTER(LEN=7) :: color(111)
  INTEGER, EXTERNAL:: AixRm

  color(  1)="#ffffff"; color(  2)="#000000"; color(  3)="#a0a0a0"; color(  4)="#ff0000"; color(  5)="#00c000"
  color(  6)="#0080ff"; color(  7)="#c000ff"; color(  8)="#00eeee"; color(  9)="#c04000"; color( 10)="#c8c800"
  color( 11)="#4169e1"; color( 12)="#ffc020"; color( 13)="#008040"; color( 14)="#c080ff"; color( 15)="#306080"
  color( 16)="#8b0000"; color( 17)="#408000"; color( 18)="#ff80ff"; color( 19)="#7fffd4"; color( 20)="#a52a2a"
  color( 21)="#ffff00"; color( 22)="#40e0d0"; color( 23)="#000000"; color( 24)="#1a1a1a"; color( 25)="#333333"
  color( 26)="#4d4d4d"; color( 27)="#666666"; color( 28)="#7f7f7f"; color( 29)="#999999"; color( 30)="#b3b3b3"
  color( 31)="#c0c0c0"; color( 32)="#cccccc"; color( 33)="#e5e5e5"; color( 34)="#ffffff"; color( 35)="#f03232"
  color( 36)="#90ee90"; color( 37)="#add8e6"; color( 38)="#f055f0"; color( 39)="#e0ffff"; color( 40)="#eedd82"
  color( 41)="#ffb6c1"; color( 42)="#afeeee"; color( 43)="#ffd700"; color( 44)="#00ff00"; color( 45)="#006400"
  color( 46)="#00ff7f"; color( 47)="#228b22"; color( 48)="#2e8b57"; color( 49)="#0000ff"; color( 50)="#00008b"
  color( 51)="#191970"; color( 52)="#000080"; color( 53)="#0000cd"; color( 54)="#87ceeb"; color( 55)="#00ffff"
  color( 56)="#ff00ff"; color( 57)="#00ced1"; color( 58)="#ff1493"; color( 59)="#ff7f50"; color( 60)="#f08080"
  color( 61)="#ff4500"; color( 62)="#fa8072"; color( 63)="#e9967a"; color( 64)="#f0e68c"; color( 65)="#bdb76b"
  color( 66)="#b8860b"; color( 67)="#f5f5dc"; color( 68)="#a08020"; color( 69)="#ffa500"; color( 70)="#ee82ee"
  color( 71)="#9400d3"; color( 72)="#dda0dd"; color( 73)="#905040"; color( 74)="#556b2f"; color( 75)="#801400"
  color( 76)="#801414"; color( 77)="#804014"; color( 78)="#804080"; color( 79)="#8060c0"; color( 80)="#8060ff"
  color( 81)="#808000"; color( 82)="#ff8040"; color( 83)="#ffa040"; color( 84)="#ffa060"; color( 85)="#ffa070"
  color( 86)="#ffc0c0"; color( 87)="#ffff80"; color( 88)="#ffffc0"; color( 89)="#cdb79e"; color( 90)="#f0fff0"
  color( 91)="#a0b6cd"; color( 92)="#c1ffc1"; color( 93)="#cdc0b0"; color( 94)="#7cff40"; color( 95)="#a0ff20"
  color( 96)="#bebebe"; color( 97)="#d3d3d3"; color( 98)="#d3d3d3"; color( 99)="#a0a0a0"; color(100)="#a0b6cd"
  color(101)="#000000"; color(102)="#1a1a1a"; color(103)="#333333"; color(104)="#4d4d4d"; color(105)="#666666"
  color(106)="#7f7f7f"; color(107)="#999999"; color(108)="#b3b3b3"; color(109)="#cccccc"; color(110)="#e5e5e5"
  color(111)="#ffffff"

  dbg=.false.
  StdOut=6
  iErr=0
  hmin=0.0_wp
  hmax=0.0_wp
  MHmin_calc=0.0_wp
  MHmax_calc=0.0_wp
  MHmin=0.0_wp
  MHmax=0.0_wp
  hmin=MINVAL(H)-0.02_wp*MAXVAL(H)
  hmax=MAXVAL(H)+0.02_wp*MAXVAL(H)
  MHmin_calc=MINVAL(MHcalc)
  MHmax_calc=MAXVAL(MHcalc)
  MHmin=MHmin_calc-0.01_wp*MHmax_calc
  MHmax=MHmax_calc+0.01_wp*MHmax_calc


  IF (dbg) WRITE (StdOut,*) 'nH        = ',nH
  IF (dbg) WRITE (StdOut,*) 'hmin      = ',hmin
  IF (dbg) WRITE (StdOut,*) 'hmax      = ',hmax
  IF (dbg) WRITE (StdOut,*) 'MHmin_calc= ',MHmin_calc
  IF (dbg) WRITE (StdOut,*) 'MHmax_calc= ',MHmax_calc
  IF (dbg) WRITE (StdOut,*) 'MHmin     = ',MHmin
  IF (dbg) WRITE (StdOut,*) 'MHmax     = ',MHmax
  IF (dbg) WRITE (StdOut,*) 'zJ        = ',zJ
  IF (dbg) THEN
    DO iH=1,nH
      WRITE(StdOut,*) H(iH), (MHcalc(iH,iTempMagn),iTempMagn=1,nTempMagn)
    END DO
  END IF

  ! generate the GNUPLOT script in the $WorkDir
  lineOut=' '
  line1=' '
  line2=' '
  gnuplot_CMD=' '
  istat=0
  file_exist=.false.
  is_file_open=.false.
  file_size=0
  execute_gnuplot_cmd =.false.
  gnuplot_version=0.0_wp

  ! check if the file lineOUT exists
  INQUIRE(FILE="lineOUT",EXIST=file_exist,OPENED=is_file_open, NUMBER=file_number)

  IF (file_exist) then
     IF (dbg) WRITE (StdOut,'(A)') 'file  "lineOUT" exists in WorkDir'
     IF (is_file_open) then
       IF (dbg) WRITE (StdOut,'(A)') 'file  "lineOUT" is opened'
       ! close the file:
       CLOSE (UNIT=file_number,STATUS='DELETE')
     END IF
     ! delete the file
     IF (dbg) WRITE (StdOut,'(A)') 'deleting the file...'
     iErr=AixRm("lineOUT")
  ELSE
     IF (dbg) WRITE (StdOut,'(A)') 'file  "lineOUT" does not exist in WorkDir'
  END IF


  ! find the gnuplot
  IF (dbg) WRITE (StdOut,'(A)') 'inquire which GNUPLOT'

  CALL systemf ( "which gnuplot >> lineOUT", iErr )
  INQUIRE(FILE="lineOUT",EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number,SIZE=file_size)

  IF (dbg) WRITE (StdOut,*) 'File_number =',file_number
  IF (dbg) WRITE (StdOut,*) 'Is_file_open=',is_file_open
  IF (dbg) WRITE (StdOut,*) 'File_exist  =',file_exist
  IF (dbg) WRITE (StdOut,*) 'File_size   =',file_size

  IF (file_exist) then
    IF (file_size>0) then
      IF (dbg) WRITE (StdOut,'(A)') 'new file  "lineOUT"  exists in WorkDir'

      file_number=453
      Call molcas_open(file_number,"lineOUT")

      READ (file_number,'(A)') line1

      IF (dbg) WRITE (StdOut,*) 'line1=',line1
      IF (dbg) WRITE (StdOut,*) trim(line1)
      line2=trim(line1)
      IF (dbg) WRITE (StdOut,*) 'line2=',line2

      CLOSE(file_number)
      IF (dbg) WRITE (StdOut,*) 'Closing lineOUT file'
      FLUSH (StdOut)
      execute_gnuplot_cmd=.true.
    ELSE
      ! file_size =0
      WRITE (StdOut,'(A)') 'file  "lineOUT" has a size=0. gnuplot was not found on the system.'
      WRITE (StdOut,'(A)') 'plots will not be created.'
    END IF
  ELSE
     WRITE (StdOut,'(A)') 'file  "lineOUT" does not exist in WorkDir'
  END IF
  ! remove file "lineOUT"
  iErr=AixRm("lineOUT")
!--------------------------------------------------------------------------------------------


! check the version of the gnuplot:
  IF ( execute_gnuplot_cmd ) Then
    IF (dbg) WRITE (StdOut,'(A)') 'inquire which version of GNUPLOT is installed'
    ! attempt to execute the script
    WRITE (gnuplot_CMD,'(2A)') trim(line2),' --version > lineOUT'
    IF (dbg) WRITE (StdOut,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
    CALL systemf ( gnuplot_CMD, iErr )
    file_number=452
    Call molcas_open(file_number,"lineOUT")
    READ (file_number,*) cdummy, gnuplot_version
    IF (dbg) WRITE (StdOut,'(A,F4.1)') 'gnuplot_version = ', gnuplot_version
    IF (abs(gnuplot_version)<0.1_wp) execute_gnuplot_cmd =.false.
    CLOSE (file_number)
    ! remove file "lineOUT"
    iErr=AixRm("lineOUT")
  END IF


!----------------------------------------------------------------------------------------
! create the file "MH.dat"
  INQUIRE(FILE="MH.dat",EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
  IF(file_exist)  iErr=AixRm("MH.dat")
  LuData=454
  Call molcas_open(LuData,"MH.dat")
  IF (dbg) WRITE (StdOut,*) 'Opening "MH.dat" file'
  write(fmtline,'(A,i0,A)') '(',nTempMagn+1,'ES24.14)'
  DO iH=1,nH
     WRITE (LuData,fmtline) H(iH), ( MHcalc(iH,iTempMagn), iTempMagn=1,nTempMagn )
  END DO
  IF (dbg) WRITE (StdOut,*) 'Writing into the "MH.dat" file'
  CLOSE (LuData)
  IF (dbg) WRITE (StdOut,*) 'Closing the "MH.dat" file'
  FLUSH (StdOut)

  ! generate the GNUPLOT script in the $WorkDir
  INQUIRE(FILE="MH.plt",EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
  IF(file_exist) iErr=AixRm("MH.plt")
  LuPlt=455
  Call molcas_open(LuPlt,"MH.plt")
  IF (dbg) WRITE (StdOut,*) 'Opening "MH.plt" file'


  IF ( gnuplot_version < 5.0_wp ) Then
  !===  GNUPLOT VERSION 4 and below ==>>  generate EPS
     WRITE (LuPlt,'(A)') 'set terminal postscript eps enhanced color  size 3.0, 2.0 font "arial, 10"'
     WRITE (LuPlt,'(A)') 'set output "MH.eps"'
     WRITE (LuPlt,'(A)') 'set grid'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A)') '# Axes'
     WRITE (LuPlt,'(A)') '  set xlabel "Magnetic Field / Tesla" font "Arial,14"'
     WRITE (LuPlt,'(A)') '  set ylabel "M / {/Symbol m}_{B}" font "Arial,14"'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[', hmin,':', hmax,']'
     WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',MHmin,':',MHmax,']'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A)') '# Tics for axes'
     WRITE (LuPlt,'(A)') '  set xtics nomirror'
     WRITE (LuPlt,'(A)') '  set ytics nomirror'
     WRITE (LuPlt,'(A)') '  set mxtics'
     WRITE (LuPlt,'(A)') '  set mytics'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A)') '# Margins'
     WRITE (LuPlt,'(A)') '  set bmargin at screen 0.20'
     WRITE (LuPlt,'(A)') '  set lmargin at screen 0.15'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A)') '# legend'
     WRITE (LuPlt,'(A)') '  set key right bottom'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A)') '# actual plotting'
     DO iTempMagn=1,nTempMagn
       ik=iTempMagn+1
       Call RANDOM_NUMBER(r)
       ic=INT(111.0_wp*r)
       IF((iTempMagn.eq.1).AND.(nTempMagn.gt.1)) THEN
         WRITE (LuPlt,'(A,i0,A,F7.3,A)') 'plot "MH.dat" using 1:',ik,' with lines lt 1  lw  8 lc rgb "'//color(ic)//&
                                          '"  title "Calc. M at T=',TempMagn(iTempMagn),'K.", \'
       ELSE IF ((iTempMagn.eq.1).AND.(nTempMagn.eq.1)) THEN
         WRITE (LuPlt,'(A,i0,A,F7.3,A)') 'plot "MH.dat" using 1:',ik,' with lines lt 1  lw  8 lc rgb "'//color(ic)//&
                                         '"  title "Calc. M at T=',TempMagn(iTempMagn),'K."'
       END IF
       IF ( (iTempMagn.lt.nTempMagn).AND.(iTempMagn.ne.1) ) THEN
         WRITE (LuPlt,'(A,i0,A,F7.3,A)') '     "MH.dat" using 1:',ik,' with lines lt 1  lw  8 lc rgb "'//color(ic)//&
                                         '"  title "Calc. M at T=',TempMagn(iTempMagn),'K.", \'
       END IF
       IF ( (iTempMagn.eq.nTempMagn).AND.(nTempMagn.gt.1) ) THEN
         WRITE (LuPlt,'(A,i0,A,F7.3,A)') '     "MH.dat" using 1:',ik,' with lines lt 1  lw  8 lc rgb "'//color(ic)//&
                                         '"  title "Calc. M at T=',TempMagn(iTempMagn),'K."'
       END IF
     END DO


  ELSE IF ( (gnuplot_version >= 5.0_wp) .AND.(gnuplot_version < 6.0_wp) ) Then
  !===  GNUPLOT VERSION 5 and above ==>>  generate PNG
     WRITE (LuPlt,'(A)') 'set terminal pngcairo transparent enhanced font "arial,10" fontscale 4.0 size 1800, 1200'
     WRITE (LuPlt,'(A)') 'set output "MH.png"'
     WRITE (LuPlt,'(A)') 'set grid'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A)') '# Axes'
     WRITE (LuPlt,'(A)') '  set xlabel "Magnetic Field / Tesla" font "Arial,14"'
     WRITE (LuPlt,'(A)') '  set ylabel "M / {/Symbol m}_{B}" font "Arial,14"'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[', hmin,':', hmax,']'
     WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',MHmin,':',MHmax,']'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A)') '# Tics for axes'
     WRITE (LuPlt,'(A)') '  set xtics nomirror'
     WRITE (LuPlt,'(A)') '  set ytics nomirror'
     WRITE (LuPlt,'(A)') '  set mxtics'
     WRITE (LuPlt,'(A)') '  set mytics'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A)') '# Margins'
     WRITE (LuPlt,'(A)') '  set bmargin at screen 0.20'
     WRITE (LuPlt,'(A)') '  set lmargin at screen 0.15'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A)') '# legend'
     WRITE (LuPlt,'(A)') '  set key right bottom'
     WRITE (LuPlt,'(A)')
     WRITE (LuPlt,'(A)') '# actual plotting'
     DO iTempMagn=1,nTempMagn
       ik=iTempMagn+1
       Call RANDOM_NUMBER(r)
       ic=INT(110.0_wp*r)
       IF((iTempMagn.eq.1).AND.(nTempMagn.gt.1)) THEN
         WRITE (LuPlt,'(A,i0,A,F7.3,A)') 'plot "MH.dat" using 1:',ik,' with lines lt 1  lw  8 lc rgb "'//color(ic)//&
                                          '"  title "Calc. M at T=',TempMagn(iTempMagn),'K.", \'
       ELSE IF ((iTempMagn.eq.1).AND.(nTempMagn.eq.1)) THEN
         WRITE (LuPlt,'(A,i0,A,F7.3,A)') 'plot "MH.dat" using 1:',ik,' with lines lt 1  lw  8 lc rgb "'//color(ic)//&
                                         '"  title "Calc. M at T=',TempMagn(iTempMagn),'K."'
       END IF
       IF ( (iTempMagn.lt.nTempMagn).AND.(iTempMagn.ne.1) ) THEN
         WRITE (LuPlt,'(A,i0,A,F7.3,A)') '     "MH.dat" using 1:',ik,' with lines lt 1  lw  8 lc rgb "'//color(ic)//&
                                         '"  title "Calc. M at T=',TempMagn(iTempMagn),'K.", \'
       END IF
       IF ( (iTempMagn.eq.nTempMagn).AND.(nTempMagn.gt.1) ) THEN
         WRITE (LuPlt,'(A,i0,A,F7.3,A)') '     "MH.dat" using 1:',ik,' with lines lt 1  lw  8 lc rgb "'//color(ic)//&
                                         '"  title "Calc. M at T=',TempMagn(iTempMagn),'K."'
       END IF
     END DO
  ELSE
    WRITE (StdOut,*) 'GNUPLOT has version: ', gnuplot_version
    WRITE (StdOut,*) 'This version of GNUPLOT is not known and thus, unsupported.'
    execute_gnuplot_cmd=.false.
  END IF
  WRITE (LuPlt,'(A)')
  CLOSE (LuPlt)


  WRITE (StdOut,'(A)') 'File "MH.png" was created in Working directory.'

  IF (execute_gnuplot_cmd) Then
    ! attempt to execute the script
    WRITE (gnuplot_CMD,'(A,A)') trim(line2),'  MH.plt'
    IF (dbg) WRITE (StdOut,'(A,A)') 'gnuplot_CMD=',gnuplot_CMD
    CALL systemf ( gnuplot_CMD, iErr )
    IF ( gnuplot_version < 5.0_wp ) Then
      WRITE (StdOut,'(A,i0,A)') 'File "MH.eps" was created in Working directory.'
    ELSE
      WRITE (StdOut,'(A,i0,A)') 'File "MH.png" was created in Working directory.'
    END IF
  END IF

  RETURN
  END SUBROUTINE plot_MH_no_Exp



