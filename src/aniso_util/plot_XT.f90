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
  Subroutine plot_XT_with_Exp(label, nT, T, XTcalc, XTexp, zJ )

  IMPLICIT NONE

  Integer, parameter    :: wp=SELECTED_REAL_KIND(p=15,r=307)
  INTEGER, INTENT(in)   :: nT
  REAL (wp), INTENT(in) :: T(nT)
  REAL (wp), INTENT(in) :: XTexp(nT)
  REAL (wp), INTENT(in) :: XTcalc(nT)
  REAL (wp), INTENT(in) :: zJ
  CHARACTER(LEN=50), intent(in) :: label
  ! local variables
  REAL (wp) :: tmin, tmax, XTmin_exp, XTmax_exp, XTmin_calc, XTmax_calc, XTmin, XTmax
  REAL (wp) :: gnuplot_version
  INTEGER               :: file_number, iT, LuPlt, LuData, file_size, StdOut, iErr
  LOGICAL               :: file_exist, is_file_open, execute_gnuplot_cmd, dbg
  CHARACTER(LEN=100)    :: line1, line2, cdummy
  CHARACTER(LEN=100)    :: datafile, plotfile, imagefile, epsfile
  INTEGER, EXTERNAL     :: AixRm
  INTEGER               :: Length
  CHARACTER(LEN=1023)   :: realname_plt, realname_dat, realname_png, realname_eps, gnuplot_CMD

  !INTEGER               :: file_number, istat, iT, LuPlt, LuData, file_size, StdOut, iErr
  !LOGICAL               :: file_exist, is_file_open, execute_gnuplot_cmd, dbg
  !CHARACTER(LEN=100)    :: line1, line2, lineOut, cdummy

  StdOut=6
  iErr=0
  dbg=.false.
  tmin=0.0_wp
  tmax=0.0_wp
  XTmin_exp=0.0_wp
  XTmax_exp=0.0_wp
  XTmin_calc=0.0_wp
  XTmax_calc=0.0_wp
  XTmin=0.0_wp
  XTmax=0.0_wp
  tmin=MINVAL(T)-0.05_wp*MAXVAL(T)
  tmax=MAXVAL(T)+0.05_wp*MAXVAL(T)
  XTmin_exp=MINVAL(XTexp)
  XTmax_exp=MAXVAL(XTexp)
  XTmin_calc=MINVAL(XTcalc)
  XTmax_calc=MAXVAL(XTcalc)
  XTmin=MIN(XTmin_exp,XTmin_calc)-0.08_wp*MAX(XTmax_exp,XTmax_calc)
  XTmax=MAX(XTmax_exp,XTmax_calc)+0.08_wp*MAX(XTmax_exp,XTmax_calc)


  IF (dbg) WRITE (StdOut,*) 'nT        = ',nT
  IF (dbg) WRITE (StdOut,*) 'tmin      = ',tmin
  IF (dbg) WRITE (StdOut,*) 'tmax      = ',tmax
  IF (dbg) WRITE (StdOut,*) 'XTmin_exp = ',XTmin_exp
  IF (dbg) WRITE (StdOut,*) 'XTmax_exp = ',XTmax_exp
  IF (dbg) WRITE (StdOut,*) 'XTmin_calc= ',XTmin_calc
  IF (dbg) WRITE (StdOut,*) 'XTmax_calc= ',XTmax_calc
  IF (dbg) WRITE (StdOut,*) 'XTmin     = ',XTmin
  IF (dbg) WRITE (StdOut,*) 'XTmax     = ',XTmax
  IF (dbg) WRITE (StdOut,*) 'zJ        = ',zJ
  IF (dbg) WRITE (StdOut,*) 'label     = ',label
  IF (dbg) THEN
    DO iT=1,nT
      WRITE(StdOut,*) T(iT), XTcalc(iT), XTexp(iT)
    END DO
  END IF

  ! generate the GNUPLOT script in the $WorkDir
  line1=' '
  line2=' '
  gnuplot_CMD=' '
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
    IF(file_size>0) then
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
    Call molcas_open(file_number,"lineOUT")
    READ (file_number,*) cdummy, gnuplot_version
    IF (dbg) WRITE (StdOut,'(A,F4.1)') 'gnuplot_version = ', gnuplot_version
    IF (abs(gnuplot_version)<0.1_wp) execute_gnuplot_cmd =.false.
    CLOSE (file_number)
    ! remove file "lineOUT"
    iErr=AixRm("lineOUT")
  END IF



!--------------------------------------------------------------------------------------------
! get the true real names of the files on disk:
  WRITE(datafile ,'(3A)')  'XT_'//trim(label)//'.dat'
  WRITE(imagefile,'(3A)')  'XT_'//trim(label)//'.png'
  WRITE(epsfile  ,'(3A)')  'XT_'//trim(label)//'.eps'
  WRITE(plotfile ,'(3A)')  'XT_'//trim(label)//'.plt'
  Call prgmtranslate(datafile ,realname_dat,Length)
  Call prgmtranslate(imagefile,realname_png,Length)
  Call prgmtranslate(epsfile  ,realname_eps,Length)
  Call prgmtranslate(plotfile ,realname_plt,Length)
  IF (dbg) THEN
    WRITE(StdOut,'(3A)') 'realname_dat=',trim(realname_dat)
    WRITE(StdOut,'(3A)') 'realname_png=',trim(realname_png)
    WRITE(StdOut,'(3A)') 'realname_eps=',trim(realname_eps)
    WRITE(StdOut,'(3A)') 'realname_plt=',trim(realname_plt)
  END IF




!  GENERATE FILES:
!--------------------------------------------------------------------------------------------
! generate the file "XT.dat":
  INQUIRE(FILE=datafile,EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
  IF(file_exist) iErr=AixRm( trim(datafile) )
  LuData=454
  Call molcas_open(LuData,datafile)
  IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(datafile)//'" file'
  IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(realname_dat)//'" file'
  DO iT=1,nT
     WRITE (LuData,'(3ES24.14)') T(iT), XTexp(iT), XTcalc(iT)
  END DO
  IF (dbg) WRITE (StdOut,*) 'Writing into the "'//trim(datafile)//'" file'
  IF (dbg) WRITE (StdOut,*) 'Writing into the "'//trim(realname_dat)//'" file'
  CLOSE (LuData)
  IF (dbg) WRITE (StdOut,*) 'Closing the "'//trim(datafile)//'" file'
  IF (dbg) WRITE (StdOut,*) 'Closing the "'//trim(realname_dat)//'" file'
  FLUSH (StdOut)



!--------------------------------------------------------------------------------------------
  ! generate the GNUPLOT script in the $WorkDir
  INQUIRE(FILE=plotfile,EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
  IF(file_exist)  iErr=AixRm( trim(plotfile) )
  LuPlt=455
  Call molcas_open(LuPlt,plotfile)
  IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(plotfile)//'" file'
  IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(realname_plt)//'" file'

  IF ( gnuplot_version < 5.0_wp ) Then
  !===  GNUPLOT VERSION 4 and below ==>>  generate EPS
    WRITE (LuPlt,'(A)') 'set terminal postscript eps enhanced color  size 3.0, 2.0 font "arial, 10"'
    WRITE (LuPlt,'(A)') 'set output "'//trim(realname_eps)//'" '
    WRITE (LuPlt,'(A)') 'set grid'
    WRITE (LuPlt,'(A)')
    WRITE (LuPlt,'(A)') '# Axes'
    WRITE (LuPlt,'(A)') '  set xlabel "Temperature / Kelvin" font "Arial,14"'
    WRITE (LuPlt,'(A)') '  set ylabel "{/Symbol c}T / cm^3 K mol^{-1}" font "Arial,14"'
    WRITE (LuPlt,'(A)')
    WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[', tmin,':', tmax,']'
    WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',XTmin,':',XTmax,']'
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
    WRITE (LuPlt,'(A)') 'plot "'//trim(realname_dat)//'" using 1:2 with points  lt 1  lw  3 lc rgb "black"  title "Experiment", \'
    WRITE (LuPlt,'(A)') '     "'//trim(realname_dat)//'" using 1:3 with lines   lt 1  lw 10 lc rgb "red"    title "Calculation"'
    WRITE (LuPlt,'(A)')

  ELSE IF ( (gnuplot_version >= 5.0_wp) .AND.(gnuplot_version < 6.0_wp) ) Then
  !===  GNUPLOT VERSION 5 and above ==>>  generate PNG
    WRITE (LuPlt,'(A)') 'set terminal pngcairo transparent enhanced font "arial,10" fontscale 4.0 size 1800, 1200'
    WRITE (LuPlt,'(A)') 'set output "'//trim(realname_png)//'" '
    WRITE (LuPlt,'(A)') 'set grid'
    WRITE (LuPlt,'(A)')
    WRITE (LuPlt,'(A)') '# Axes'
    WRITE (LuPlt,'(A)') '  set xlabel "Temperature / Kelvin" font "Arial,14"'
    WRITE (LuPlt,'(A)') '  set ylabel "{/Symbol c}T / cm^3 K mol^{-1}" font "Arial,14"'
    WRITE (LuPlt,'(A)')
    WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[', tmin,':', tmax,']'
    WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',XTmin,':',XTmax,']'
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
    WRITE (LuPlt,'(A)') 'plot "'//trim(realname_dat)//'" using 1:2 with points  lt 1  lw  3 lc rgb "black"  title "Experiment", \'
    WRITE (LuPlt,'(A)') '     "'//trim(realname_dat)//'" using 1:3 with lines   lt 1  lw 10 lc rgb "red"    title "Calculation"'
    WRITE (LuPlt,'(A)')
    WRITE (LuPlt,'(A)')
  ELSE
    WRITE (StdOut,*) 'GNUPLOT has version: ', gnuplot_version
    WRITE (StdOut,*) 'This version of GNUPLOT is not known and thus, unsupported.'
    execute_gnuplot_cmd=.false.
  END IF
  CLOSE (LuPlt)

  IF (execute_gnuplot_cmd) Then
    ! attempt to execute the script
    IF (dbg) THEN
      WRITE (StdOut,*) trim(realname_plt)
      INQUIRE(FILE=trim(realname_plt),EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
      IF(file_exist) THEN
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_plt)//'" exists.'
      ELSE
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_plt)//'" does not exist.'
      END IF
    END IF

    WRITE (gnuplot_CMD,'(5A)') trim(line2),' ',trim(realname_plt)
    IF (dbg) WRITE (StdOut,'(A,A)') 'gnuplot_CMD=',trim(gnuplot_CMD)

    CALL execute_command_line ( trim(gnuplot_CMD) )

    IF ( gnuplot_version < 5.0_wp ) Then
      INQUIRE(FILE=trim(realname_eps),EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
      IF(file_exist) THEN
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_eps)//'" was created in Working directory.'
      ELSE
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_eps)//'" was NOT created in Working directory.'
      END IF
    ELSE
      INQUIRE(FILE=trim(realname_png),EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
      IF(file_exist) THEN
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_png)//'" was created in Working directory.'
      ELSE
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_png)//'" was NOT created in Working directory.'
      END IF
    END IF
  END IF

  dbg=.false.
  RETURN
#ifdef _WARNING_WORKAROUND_
  IF (.FALSE.) CALL UNUSED_CHARACTER(cdummy)
#endif
  END SUBROUTINE plot_XT_with_Exp







  Subroutine plot_XT_no_Exp(label, nT, T, XTcalc, zJ )

  IMPLICIT NONE

  Integer, parameter    :: wp=SELECTED_REAL_KIND(p=15,r=307)
  INTEGER, INTENT(in)   :: nT
  REAL (wp), INTENT(in) :: T(nT)
  REAL (wp), INTENT(in) :: XTcalc(nT)
  REAL (wp), INTENT(in) :: zJ
  CHARACTER(LEN=50), intent(in) :: label
  ! local variables
  REAL (wp) :: tmin, tmax, XTmin_calc, XTmax_calc, XTmin, XTmax
  REAL (wp) :: gnuplot_version
!<<<<<<< HEAD
  INTEGER               :: file_number, istat, iT, LuPlt, LuData, file_size, StdOut, iErr
  LOGICAL               :: file_exist, is_file_open, execute_gnuplot_cmd, dbg
  CHARACTER(LEN=100)    :: line1, line2, lineOut, cdummy
  CHARACTER(LEN=100)    :: datafile, plotfile, imagefile, epsfile
  INTEGER, EXTERNAL     :: AixRm
  INTEGER               :: Length
  CHARACTER(LEN=1024)   :: realname_plt, realname_dat, realname_png, realname_eps, gnuplot_CMD
!=======
!  INTEGER           :: file_number, iT, LuPlt, LuData, file_size, StdOut, iErr
!  LOGICAL           :: file_exist, is_file_open, execute_gnuplot_cmd, dbg
!  CHARACTER(LEN=100):: line1, line2, cdummy
!  CHARACTER(LEN=100):: gnuplot_CMD, datafile, plotfile
!  INTEGER, EXTERNAL :: AixRm
!>>>>>>> upstream/master

  StdOut=6
  iErr=0
  dbg=.false.
  tmin=0.0_wp
  tmax=0.0_wp
  XTmin_calc=0.0_wp
  XTmax_calc=0.0_wp
  XTmin=0.0_wp
  XTmax=0.0_wp
  tmin=MINVAL(T)-0.02_wp*MAXVAL(T)
  tmax=MAXVAL(T)+0.02_wp*MAXVAL(T)
  XTmin_calc=MINVAL(XTcalc)
  XTmax_calc=MAXVAL(XTcalc)
  XTmin=XTmin_calc-0.01_wp*XTmax_calc
  XTmax=XTmax_calc+0.01_wp*XTmax_calc


  IF (dbg) WRITE (StdOut,*) 'nT        = ',nT
  IF (dbg) WRITE (StdOut,*) 'tmin      = ',tmin
  IF (dbg) WRITE (StdOut,*) 'tmax      = ',tmax
  IF (dbg) WRITE (StdOut,*) 'XTmin_calc= ',XTmin_calc
  IF (dbg) WRITE (StdOut,*) 'XTmax_calc= ',XTmax_calc
  IF (dbg) WRITE (StdOut,*) 'XTmin     = ',XTmin
  IF (dbg) WRITE (StdOut,*) 'XTmax     = ',XTmax
  IF (dbg) WRITE (StdOut,*) 'zJ        = ',zJ
  IF (dbg) WRITE (StdOut,*) 'label     = ',trim(label)
  IF (dbg) THEN
    DO iT=1,nT
      WRITE(StdOut,*) T(iT), XTcalc(iT)
    END DO
  END IF

  ! generate the GNUPLOT script in the $WorkDir
  line1=' '
  line2=' '
  gnuplot_CMD=' '
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
    IF(file_size>0) then
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
      execute_gnuplot_cmd =.true.
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



!--------------------------------------------------------------------------------------------
! get the true real names of the files on disk:
  WRITE(datafile ,'(3A)')  'XT_'//trim(label)//'.dat'
  WRITE(imagefile,'(3A)')  'XT_'//trim(label)//'.png'
  WRITE(epsfile  ,'(3A)')  'XT_'//trim(label)//'.eps'
  WRITE(plotfile ,'(3A)')  'XT_'//trim(label)//'.plt'
  Call prgmtranslate(datafile ,realname_dat,Length)
  Call prgmtranslate(imagefile,realname_png,Length)
  Call prgmtranslate(epsfile  ,realname_eps,Length)
  Call prgmtranslate(plotfile ,realname_plt,Length)
  IF (dbg) THEN
    WRITE(StdOut,'(3A)') 'realname_dat=',trim(realname_dat)
    WRITE(StdOut,'(3A)') 'realname_png=',trim(realname_png)
    WRITE(StdOut,'(3A)') 'realname_eps=',trim(realname_eps)
    WRITE(StdOut,'(3A)') 'realname_plt=',trim(realname_plt)
  END IF



!  GENERATE FILES:
!--------------------------------------------------------------------------------------------
  INQUIRE(FILE=datafile,EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
  IF(file_exist) iErr=AixRm(trim(datafile))
  LuData=454
  Call molcas_open(LuData,datafile)
  IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(datafile)//'" file'
  IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(realname_dat)//'" file'
  DO iT=1,nT
     WRITE (LuData,'(3ES24.14)') T(iT), XTcalc(iT)
  END DO
  IF (dbg) WRITE (StdOut,*) 'Writing into the "'//trim(datafile)//'" file'
  IF (dbg) WRITE (StdOut,*) 'Writing into the "'//trim(realname_dat)//'" file'
  CLOSE (LuData)
  IF (dbg) WRITE (StdOut,*) 'Closing the "'//trim(datafile)//'" file'
  IF (dbg) WRITE (StdOut,*) 'Closing the "'//trim(realname_dat)//'" file'
  FLUSH (StdOut)



!--------------------------------------------------------------------------------------------
  ! generate the GNUPLOT script in the $WorkDir
  INQUIRE(FILE=plotfile,EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
  IF(file_exist)  iErr=AixRm( trim(plotfile) )
  LuPlt=455
  Call molcas_open(LuPlt,plotfile)
  IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(plotfile)//'" file'
  IF (dbg) WRITE (StdOut,*) 'Opening "'//trim(realname_plt)//'" file'

  IF ( gnuplot_version < 5.0_wp ) Then
  !===  GNUPLOT VERSION 4 and below ==>>  generate EPS
    WRITE (LuPlt,'(A)') 'set terminal postscript eps enhanced color  size 3.0, 2.0 font "arial, 10"'
    WRITE (LuPlt,'(A)') 'set output "'//trim(realname_eps)//'" '
    WRITE (LuPlt,'(A)') 'set grid'
    WRITE (LuPlt,'(A)')
    WRITE (LuPlt,'(A)') '# Axes'
    WRITE (LuPlt,'(A)') '  set xlabel "Temperature / Kelvin" font "Arial,14"'
    WRITE (LuPlt,'(A)') '  set ylabel "{/Symbol c}T / cm^3 K mol^{-1}" font "Arial,14"'
    WRITE (LuPlt,'(A)')
    WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[', tmin,':', tmax,']'
    WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',XTmin,':',XTmax,']'
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
    WRITE (LuPlt,'(A)') 'plot "'//trim(realname_dat)//'" using 1:2 with lines lt 1  lw  8 lc rgb "red"  title "Calculation"'
    WRITE (LuPlt,'(A)')

  ELSE IF ( (gnuplot_version >= 5.0_wp) .AND.(gnuplot_version < 6.0_wp) ) Then
    !===  GNUPLOT VERSION 5 and above ==>>  generate PNG
    WRITE (LuPlt,'(A)') 'set terminal pngcairo transparent enhanced font "arial,10" fontscale 4.0 size 1800, 1200'
    WRITE (LuPlt,'(A)') 'set output "'//trim(realname_png)//'" '
    WRITE (LuPlt,'(A)') 'set grid'
    WRITE (LuPlt,'(A)')
    WRITE (LuPlt,'(A)') '# Axes'
    WRITE (LuPlt,'(A)') '  set xlabel "Temperature / Kelvin" font "Arial,14"'
    WRITE (LuPlt,'(A)') '  set ylabel "{/Symbol c}T / cm^3 K mol^{-1}" font "Arial,14"'
    WRITE (LuPlt,'(A)')
    WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set xrange[', tmin,':', tmax,']'
    WRITE (LuPlt,'(A,F8.4,A,F10.4,A)') '  set yrange[',XTmin,':',XTmax,']'
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
    WRITE (LuPlt,'(A)') 'plot "'//trim(realname_dat)//'" using 1:2 with lines lt 1  lw  8 lc rgb "red"  title "Calculation"'
    WRITE (LuPlt,'(A)')
  ELSE
    WRITE (StdOut,*) 'GNUPLOT has version: ', gnuplot_version
    WRITE (StdOut,*) 'This version of GNUPLOT is not known and thus, unsupported.'
    execute_gnuplot_cmd=.false.
  END IF

  IF (execute_gnuplot_cmd) Then
    ! attempt to execute the script
    IF (dbg) THEN
      WRITE (StdOut,*) trim(realname_plt)
      INQUIRE(FILE=trim(realname_plt),EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
      IF(file_exist) THEN
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_plt)//'" exists.'
      ELSE
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_plt)//'" does not exist.'
      END IF
    END IF

    WRITE (gnuplot_CMD,'(5A)') trim(line2),' ',trim(realname_plt)
    IF (dbg) WRITE (StdOut,'(A,A)') 'gnuplot_CMD=',trim(gnuplot_CMD)

    CALL execute_command_line ( trim(gnuplot_CMD) )

    IF ( gnuplot_version < 5.0_wp ) Then
      INQUIRE(FILE=trim(realname_eps),EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
      IF(file_exist) THEN
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_eps)//'" was created in Working directory.'
      ELSE
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_eps)//'" was NOT created in Working directory.'
      END IF
    ELSE
      INQUIRE(FILE=trim(realname_png),EXIST=file_exist,OPENED=is_file_open,NUMBER=file_number)
      IF(file_exist) THEN
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_png)//'" was created in Working directory.'
      ELSE
        WRITE (StdOut,'(A,i0,A)') 'File "'//trim(realname_png)//'" was NOT created in Working directory.'
      END IF
    END IF
  END IF

  CLOSE (LuPlt)
  dbg=.false.
  RETURN
#ifdef _WARNING_WORKAROUND_
  IF (.FALSE.) CALL UNUSED_CHARACTER(cdummy)
#endif
  END SUBROUTINE plot_XT_no_Exp
