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
! Copyright (C) 1991, Bjorn O. Roos                                    *
!***********************************************************************

subroutine Vibinp(ncase,ngrid,nvib,Umin,Umax,Rout,PotR,E0,dE0,Redm,Teas,Req,sc,temp)
!***********************************************************************
!  Object: Read input and construct default input parameters           *
!  written by Bjoern Roos in March 1991                                *
!***********************************************************************

use Vibrot_globals, only: Atom1, Atom2, dRo, EoutO, iad12, iad13, iadvib, iallrot, IfPrWf, iobs, iplot, iscale, ispc, J1A, J1B, &
                          J2A, J2B, lambda, n0, n02, nop, npin, npobs, npoint, nRot_Max, nvib1, nvib21, Obsin, R0o, R1o, RinO, &
                          Titobs, Vibwvs, Vibwvs1, Vibwvs2
use Constants, only: Zero, One, Five, UTOAU
use Definitions, only: wp, iwp, r8, u6

implicit none
integer(kind=iwp), intent(out) :: ncase, ngrid, nvib
real(kind=wp), intent(out) :: Umin, Umax, E0, dE0, Redm, Teas, Req, sc, temp
real(kind=wp), intent(out) :: Rout(npoint+4), PotR(npoint+4)
integer(kind=iwp) :: LuIn, LuIn1, ntit1, ntit2, i, j, k, ii, kk, iadvi1, iadvi2, ich1, ich2, IERR, iOpt, iplotp, ipot, isn1, &
                     isn1x, isn2, isn2x, ist, lambdx, ngridx, nobsi
logical(kind=iwp) :: skip, exists
real(kind=wp) :: del, dRp, O0, Oeq, R0p, R1p, Redm1, Redmx, Reqx, Rmax, Rmin, Rn1, Rout0, U, Umaxx, Uminx, xMass1, xMass2, xxx, &
                 Rin(npin), Ein(npin)
character(len=2) :: At1x, At2x
character(len=4) :: word, Diatom, Diatomx
! For storing character data using gather/scatter DAFILE operations, it
! is imperative that the strings are aligned on integers.
! Assume an 8-byte character string is always dividible into integers.
! Else, IntCh and (maybe) Title1 and Title2 must be changed.
character(len=8) :: IntCh
character(len=80) :: Title1(10), Title2(10)
character(len=180) :: Line, l84, l84x
integer(kind=iwp), parameter :: ntab = 19
character(len=4), parameter :: tabinp(ntab) = ['TITL','ATOM','GRID','RANG','VIBR','ROTA','ORBI','NOSP','OBSE','STEP', &
                                               'POTE','ROVI','TRAN','ASYM','PRWF','SCAL','TEMP','ALLR','END ']
integer(kind=iwp), external :: IsFreeUnit, iNuclearChargeFromSymbol, iMostAbundantIsotope
real(kind=r8), external :: dNuclearMass
character(len=180), external :: Get_Ln, Get_Ln_EOF

LuIn = IsFreeUnit(11)
call SpoolInp(LuIn)

! Set default values to input variables

Atom1 = ''
Atom2 = ''
ipot = 0        ! Indicator for potential input
ngrid = 199     ! Maximum number of grid points
Rmin = One
Rmax = Five     ! Integration range
Umin = log(Rmin)
Umax = log(Rmax)
nvib = 3
nvib1 = nvib-1  ! Number of vibrational states
J1A = 0
J2A = 5         ! Range of rotational quantum numbers
lambda = 0      ! Orbital angular momentum quantum number
dE0 = 4.0e-4_wp ! Start value for step size in eigenvalue search
ispc = 1        ! Spectroscopic constants will be computed
iobs = 0        ! No calculation of matrix elements for operators
Teas = Zero     ! Asymtotic energy difference between two potentials
iplotp = 0      ! No plot file for the potential
iscale = 0      ! No scaling of input potential (BOR in June 2001)
temp = 300.0_wp ! Temperature for Boltzmann weighting
iallrot = 0     ! No calculation of transition between all rot. levels
IfPrWf = 0      ! PAM97: New key word PRWF, flagged by IfPrWf > 0, default=0:
Title1(:) = ''

! Position input file

rewind(LuIn)
call RdNLst(LuIn,'VibRot')

! Read input data from input file

ntit1 = 0
skip = .false.
input: do
  if (skip) then
    skip = .false.
  else
    read(LuIn,'(a)') line
    call Upcase(line)
    if (line(1:1) == '*') cycle input
    word = line(1:4)
    if (Word == '') Word = 'END'
  end if

  select case (word)

    case (tabinp(1))
      ! Read title lines. Maximum 10 is allowed.
      do
        read(LuIn,'(a)') line
        word = line(1:4)
        call Upcase(word)
        do i=1,ntab
          if (word == tabinp(i)) then
            skip = .true.
            cycle input
          end if
        end do
        if (ntit1 < 10) then
          ntit1 = ntit1+1
          read(line,'(a)') Title1(ntit1)
        end if
      end do

    case (tabinp(2))
      ! Read isotope numbers
      ! isn1 and isn2 are the isotope numbers for atoms 1 and 2,
      ! Atom1 and Atom2 the corresponding chemical symbols
      ! Masses are read from the next one/two lines if one/two
      ! of the isotope numbers are negative.
      ! Else the masses are obtained from the tables in ISOTOPE.
      ! if isn=0 the mass for the most abundant isotope will be used
      Line = Get_Ln(LuIn)
      call Upcase(line)
      ii = 0
      ist = 0
      do i=1,80
        if (line(i:i) == ' ') cycle
        if (i < ist) cycle
        ii = ii+1
        kk = 0
        do k=i,80
          if (line(k:k) == ' ') exit
          kk = kk+1
        end do
        ist = i+kk
        if (ii == 1) then
          call Get_I1(1,isn1)
        else if (ii == 2) then
          call Get_S(2,Atom1(1:kk),1)
        else if (ii == 3) then
          call Get_I1(3,isn2)
        else if (ii == 4) then
          call Get_S(4,Atom2(1:kk),1)
          exit
        end if
      end do

      if (Atom1 == 'D') then
        ich1 = 1
      else if (Atom1 == 'T') then
        ich1 = 1
      else
        ich1 = iNuclearChargeFromSymbol(Atom1)
      end if
      if (isn1 == 0) then
        if (Atom1 == 'D') then
          isn1 = 2
        else if (Atom1 == 'T') then
          isn1 = 3
        else
          isn1 = iMostAbundantIsotope(ich1)
        end if
        xMass1 = dNuclearMass(ich1,isn1)
      else if (isn1 > 0) then
        xMass1 = dNuclearMass(ich1,isn1)
      else
        read(LuIn,*) xMass1
        xMass1 = UTOAU*xMass1
      end if

      if (Atom2 == 'D') then
        ich2 = 1
      else if (Atom2 == 'T') then
        ich2 = 1
      else
        ich2 = iNuclearChargeFromSymbol(Atom2)
      end if
      if (isn2 == 0) then
        if (Atom2 == 'D') then
          isn2 = 2
        else if (Atom2 == 'T') then
          isn2 = 3
        else
          isn2 = iMostAbundantIsotope(ich2)
        end if
        xMass2 = dNuclearMass(ich2,isn2)
      else if (isn2 > 0) then
        xMass2 = dNuclearMass(ich2,isn2)
      else
        read(LuIn,*) xMass2
        xMass2 = UTOAU*xMass2
      end if

      Redm = xMass1*xMass2/(xMass1+xMass2)

    case (tabinp(3))
      ! Read number of grid points for numerical integration (max: npoint)
      Line = Get_Ln(LuIn)
      call Get_I1(1,ngrid)
      if (mod(ngrid,2) == 0) ngrid = ngrid-1 ! ngrid should be odd
      if (ngrid >= npoint) ngrid = npoint-1

    case (tabinp(4))
      ! Read upper and lower integration range in atomic units
      Line = Get_Ln(LuIn)
      call Get_F1(1,Rmin)
      call Get_F1(2,RMax)
      Umin = log(Rmin)
      Umax = log(Rmax)

    case (tabinp(5))
      ! Read number of vibrational quantum numbers
      Line = Get_Ln(LuIn)
      call Get_I1(1,nvib)
      n0 = 0
      nvib1 = nvib-1

    case (tabinp(6))
      ! Read range for rotational quantum numbers
      Line = Get_Ln(LuIn)
      call Get_I1(1,J1A)
      call Get_I1(2,J2A)
      if (J2A >= nRot_Max) then
        write(u6,*)
        write(u6,*) '********************************'
        write(u6,*) ' VIBINP Error: J2A > nRot_Max. '
        write(u6,'(1x,a,2i6)') 'J2A:',J2A
        write(u6,*) '********************************'
        call Quit_OnUserError()
      end if

    case (tabinp(7))
      ! Read orbital angular momentum quantum number
      Line = Get_Ln(LuIn)
      call Get_I1(1,lambda)

    case (tabinp(8))
      ! Read flag for spectroscopic constants
      ispc = 0

    case (tabinp(9))
      ! Read input for calculation of matrix elements of observables
      ! like the dipole operator, etc.
      iobs = iobs+1
      iplot(iobs) = 0
      if (iobs > 10) then
        write(u6,*)
        write(u6,*) '***************************'
        write(u6,*) ' VIBINP Error: IOBS > 10. '
        write(u6,'(1x,a,2i6)') 'IOBS:',IOBS
        write(u6,*) '***************************'
        call Quit_OnUserError()
      end if

      !VV: make it really nasty. if it looks like a file - read from the file
      !    else-  this is a title
      read(LuIn,'(a)') Line
      call f_inquire(Line(1:index(Line,' ')-1),exists)
      if (exists) then
        LuIn1 = IsFreeUnit(15)
        call molcas_open(LuIn1,Line(1:index(Line,' ')-1))
        Line = Get_Ln(LuIn1)
      else
        LuIn1 = LuIn
      end if
      Titobs(iobs) = trim(Line)
      !write(u6,*) ' In VIBINP. IOBS=',IOBS
      !write(u6,*) ' TITOBS just read:'
      !write(u6,'(a80)') TITOBS(IOBS)
      nobsi = 0
      do
        Line = Get_Ln_EOF(LuIn1)
        if ((Line(1:3) == 'EOF') .and. (LuIn1 /= LuIn)) then
          close(LuIn1)
          Line = Get_Ln(LuIn)
        end if
        word = line(1:4)
        call Upcase(word)
        do i=1,ntab
          if (word == tabinp(i)) then
            skip = .true.
            cycle input
          end if
        end do
        if (word == 'PLOT') exit
        nobsi = nobsi+1
        if (NOBSI > NPIN) then
          write(u6,*)
          write(u6,*) '**********************************'
          write(u6,*) ' VIBINP Error: NOBSI > NPIN      '
          write(u6,'(1x,a,2i6)') 'NOBSI,NPIN:',NOBSI,NPIN
          write(u6,*) '**********************************'
          call Quit_OnUserError()
        end if
        npobs(iobs) = nobsi
        call Get_F1(1,RinO(nobsi,iobs))
        call Get_F1(2,Obsin(nobsi,iobs))
      end do
      do
        Line = Get_Ln(LuIn)
        if (line(1:1) /= '*') exit
      end do
      call Get_F1(1,R0o(iobs))
      call Get_F1(2,R1o(iobs))
      call Get_F1(3,dRo(iobs))
      iplot(iobs) = 1

    case (tabinp(10))
      ! Read starting value for step size in eigenvalue search
      Line = Get_Ln(LuIn)
      call Get_F1(1,dE0)

    case (tabinp(11))
      ! Read potential
      ipot = 1
      nop = 0
      Line = Get_Ln(LuIn)
      xxx = 999.97_wp
      read(line,*,iostat=IERR) xxx
      if ((xxx == 999.97_wp) .or. (IERR /= 0)) then
        call f_inquire(line,exists)
        write(u6,*) line
        if (.not. exists) then
          write(u6,*) 'File with potential non-existent'
          write(u6,*) 'File=',trim(line)
          call Abend()
        end if
        LuIn1 = IsFreeUnit(15)
        call molcas_open(LuIn1,line)
      else
        backspace(LuIn)
        Luin1 = Luin
      end if
      do
        Line = Get_Ln_EOF(LuIn1)
        if ((Line(1:3) == 'EOF') .and. (LuIn1 /= LuIn)) then
          close(LuIn1)
          Line = Get_Ln(LuIn)
        end if
        word = line(1:4)
        call Upcase(word)
        do i=1,ntab
          if (word == tabinp(i)) then
            skip = .true.
            cycle input
          end if
        end do
        if (word == 'PLOT') exit
        nop = nop+1
        if (NOP > NPIN) then
          write(u6,*)
          write(u6,*) '**********************************'
          write(u6,*) ' VIBINP Error: NOP > NPIN'
          write(u6,'(1x,a,2i6)') 'NOP,NPIN:',NOP,NPIN
          write(u6,*) '**********************************'
          call Quit_OnUserError()
        end if
        call Get_F(1,Rin(nop),1)
        call Get_F(2,Ein(nop),1)
      end do
      Line = Get_Ln(LuIn)
      call Get_F1(1,R0p)
      call Get_F1(2,R1p)
      call Get_F1(3,dRp)
      iplotp = 1
      if (LuIn1 /= LuIn) close(LuIn1)

    case (tabinp(12))
      ! Calculation of ro-vibrational wave functions (ncase=1)
      ncase = 1

    case (tabinp(13))
      ! Calculation of transition moments (ncase=2)
      ncase = 2

    case (tabinp(14))
      ! Asymptotic energy difference between two potentials
      Line = Get_Ln(LuIn)
      call Get_F1(1,Teas)

    case (tabinp(15))
      ! Flag for printing the wave function.
      IfPrWf = 1

    case (tabinp(16))
      ! Scaling of input potential such that the binding energy is 0.1 au.
      iscale = 1

    case (tabinp(17))
      ! Temperature for vibrational averaging
      Line = Get_Ln(LuIn)
      call Get_F1(1,Temp)

    case (tabinp(18))
      ! ALLRotational
      iallrot = 1

    case (tabinp(19))
      exit input

    case default
      write(u6,*)
      write(u6,*) '******************************************'
      write(u6,*) ' VIBINP Error: Input line not recognized. '
      write(u6,*) ' Input line, in upper case:               '
      write(u6,'(a)') line
      write(u6,*) ' Extracted keyword: ',word
      write(u6,*) '******************************************'
      call Quit_OnUserError()
  end select
end do input

if (J1A < lambda) then
  write(u6,*)
  write(u6,*) '********************************'
  write(u6,*) ' VIBINP Warning: J1A < Lambda. '
  write(u6,'(1x,a,2i6)') 'J1A,Lambda:',J1A,Lambda
  write(u6,*) ' J1A is now reset=Lambda.       '
  write(u6,*) '********************************'
end if

if (J2A < J1A) then
  write(u6,*)
  write(u6,*) '***************************'
  write(u6,*) ' VIBINP Error: J2A < J1A. '
  write(u6,'(1x,a,2i6)') 'J2A,J1A:',J2A,J1A
  write(u6,*) '***************************'
  call Quit_OnUserError()
end if

! Check for input error

if (ncase == 1) then
  if ((Atom1 == '  ') .or. (Atom2 == '  ')) then
    write(u6,*)
    write(u6,*) '**********************************'
    write(u6,*) ' VIBINP Error: No atoms in input.'
    write(u6,*) '**********************************'
    call Quit_OnUserError()
  end if
end if

! Retrieve data from Vibwvs1 and Vibwvs2 (ncase=2)

if (ncase == 2) then
  iadvi1 = 0
  call iDafile(Vibwvs1,2,iad12,100,iadvi1)
  iOpt = 2
  call WR_VibRot_Info1(Vibwvs1,iOpt,iadvi1,ntit1,J1A,J2A,lambda,n0,nvib1,Redm,Umax,Umin,ngrid,isn1,isn2,Req,xMass1,xMass2)
  call cDaFile(Vibwvs1,iOpt,Title1,10*80,iadvi1)
  call cDaFile(Vibwvs1,iOpt,IntCh,8,iadvi1)

  Atom1 = IntCh(1:2)
  Atom2 = IntCh(5:6)
  Rmin = exp(Umin)
  Rmax = exp(Umax)
  iadvi2 = 0
  call iDafile(Vibwvs2,2,iad13,100,iadvi2)
  iOpt = 2
  call WR_VibRot_Info1(Vibwvs2,iOpt,iadvi2,ntit2,J1B,J2B,lambdx,n02,nvib21,Redmx,Umaxx,Uminx,ngridx,isn1x,isn2x,Reqx,xMass1,xMass2)
  call cDaFile(Vibwvs2,iOpt,Title2,10*80,iadvi2)
  call cDaFile(Vibwvs2,iOpt,IntCh,8,iadvi2)
  At1x = IntCh(1:2)
  At2x = IntCh(5:6)

  ! Check for consistency of data on files

  IERR = 0
  if ((Atom1 /= At1x) .or. (Atom2 /= At2x)) then
    write(u6,*)
    write(u6,*) '***************************************'
    write(u6,*) ' VIBINP Error: Inconsistent data.'
    write(u6,'(1x,a,1x,a,1x,a)') 'ATOM1,AT1X:',ATOM1,AT1X
    IERR = 1
  end if
  if (abs(Redm-Redmx) > 1.d-06) then
    write(u6,*)
    write(u6,*) '***************************************'
    write(u6,*) ' VIBINP Error: REDM /= REDMX'
    write(u6,'(1x,a,2f16.6)') 'REDM,REDMX:',REDM,REDMX
    IERR = 1
  end if
  if (Umax /= Umaxx) then
    write(u6,*)
    write(u6,*) '***************************************'
    write(u6,*) ' VIBINP Error: UMAX /= UMAXX'
    write(u6,'(1x,a,2f16.6)') 'UMAX,UMAXX:',UMAX,UMAXX
    IERR = 1
  end if
  if (Umin /= Uminx) then
    write(u6,*)
    write(u6,*) '***************************************'
    write(u6,*) ' VIBINP Error: UMIN /= UMINX'
    write(u6,'(1x,a,2f16.6)') 'UMIN,UMINX:',UMIN,UMINX
    IERR = 1
  end if
  if (NGrid /= NGridx) then
    write(u6,*)
    write(u6,*) '***************************************'
    write(u6,*) ' VIBINP Error: NGrid /= NGridX'
    write(u6,'(1x,a,2i8)') 'NGrid,NGridX:',NGrid,NGridX
    IERR = 1
  end if
  if (Isn1 /= Isn1x) then
    write(u6,*)
    write(u6,*) '***************************************'
    write(u6,*) ' VIBINP Error: Isn1 /= Isn1X'
    write(u6,'(1x,a,2i8)') 'Isn1,Isn1X:',Isn1,Isn1X
    IERR = 1
  end if
  if (Isn2 /= Isn2x) then
    write(u6,*)
    write(u6,*) '***************************************'
    write(u6,*) ' VIBINP Error: Isn2 /= Isn2X'
    write(u6,'(1x,a,2i8)') 'Isn2,Isn2X:',Isn2,Isn2X
    IERR = 1
  end if
  if (IERR /= 0) then
    write(u6,*) '*****************************************'
    write(u6,*) ' VIBINP: Irrecoverable errors.'
    write(u6,*) ' Transition calculation, but data on the '
    write(u6,*) ' two VibWvs files are not compatible.'
    write(u6,*) '*****************************************'
    call Quit_OnUserError()
  end if
end if

! Section for print output of input information

ii = 0
Diatom = ''
Diatomx(1:2) = Atom1
Diatomx(3:4) = Atom2
do i=1,4
  if (Diatomx(i:i) /= ' ') then
    ii = ii+1
    Diatom(ii:ii) = Diatomx(i:i)
  end if
end do
if (ncase == 1) write(u6,1100) Diatom
if (ncase == 2) write(u6,1101) Diatom

if (ncase == 2) write(u6,1190)
do i=1,ntit1
  write(u6,*) Title1(i)
end do
write(u6,1200) J1A,J2A,lambda,n0,nvib1,isn1,Atom1,xMass1,isn2,Atom2,xMass2,Redm
if (ncase == 2) then
  write(u6,1191)
  do i=1,ntit2
    write(u6,*) Title2(i)
  end do
  write(u6,1201) J1B,J2B,lambdx,n02,nvib21,isn1,At1x,xMass1,isn2,At2x,xMass2,Redmx
end if

! Numerical integration data

del = (Umax-Umin)/(ngrid-1)
write(u6,1300) ngrid,del,Rmin,Rmax,Umin,Umax

if ((ispc == 0) .and. (ncase == 1)) write(u6,1400)
if ((ispc /= 0) .and. (ncase == 1)) write(u6,1500)
if (iobs == 0) write(u6,1600)
if (iobs /= 0) write(u6,1700) iobs
if ((ipot == 0) .and. (ncase == 1)) then
  write(u6,*)
  write(u6,*) '**********************************'
  write(u6,*) ' VIBINP Error: IPOT=0 and NCASE=1 '
  write(u6,*) '**********************************'
  call Quit_OnUserError()
end if

! Print input potential

if (ipot /= 0) then
  write(u6,1800)
  do i=1,nop
    write(u6,1810) Rin(i),Ein(i)
  end do
  if ((Rin(1) > Rmax) .or. (Rin(nop) < Rmin)) then
    write(u6,*)
    write(u6,*) '****************************'
    write(u6,*) ' VIBINP Error: Potential is '
    write(u6,*) ' not within defined range.  '
    write(u6,*) '****************************'
    call Quit_OnUserError()
  end if
  if (iplotp /= 0) write(u6,1820)
end if

! Compute radial coordinates

U = Umin-del
Rout0 = exp(U)
do i=1,ngrid
  U = U+del
  Rout(i) = exp(U)
end do
Rn1 = exp(U+del)
Rout(ngrid+1) = Rout0
Rout(ngrid+2) = Rn1
if (ipot /= 0) call POT(Rin,Ein,Rout,PotR,ngrid+2,1,E0,Req,R0p,R1p,dRp,nop,Title1(1),iplotp,Redm,sc,0)
iadvib = 0
if (ncase == 1) then
  ! Store data on Vibwvs (ncase=1)
  call iDafile(Vibwvs,1,iad12,100,iadvib)
  iOpt = 1
  Redm1 = Redm*sc
  isn1 = abs(isn1)
  isn2 = abs(isn2)
  call WR_VibRot_Info1(Vibwvs,iOpt,iadvib,ntit1,J1A,J2A,lambda,n0,nvib1,Redm1,Umax,Umin,ngrid,isn1,isn2,Req,xMass1,xMass2)
  IntCh(1:4) = Atom1
  IntCh(5:8) = Atom2
  call cDaFile(Vibwvs,iOpt,Title1,10*80,iadvib)
  call cDaFile(Vibwvs,iOpt,IntCh,8,iadvib)
end if

! Fit observable input
Rout0 = Rout(ngrid+1)
do i=1,iobs
  Rout(ngrid+1) = Req
  call POT(RinO(1,i),Obsin(1,i),Rout,EoutO(1,i),ngrid+1,2,O0,Oeq,R0o(i),R1o(i),dRo(i),npobs(i),Titobs(i),iplot(i),Redm,sc,iobs)
end do
Rout(ngrid+1) = Rout0

! Print observable input data

if (iobs /= 0) then
  write(u6,1900)
  do i=1,iobs,4
    write(u6,1910) (Titobs(k)(1:18),k=i,min(i+3,iobs))
    l84 = ' '
    do k=i,min(i+3,iobs)
      l84(24*(k-i)+1:24*(k-i)+10) = '     R(au)'
      l84(24*(k-i)+11:24*(k-i+1)) = '        Value '
    end do
    write(u6,'(a)') trim(l84)
    do j=1,npin
      l84 = ' '
      do k=i,min(i+3,iobs)
        if (j <= npobs(k)) write(l84(24*(k-i)+1:24*(k-i+1)),'(1x,F9.4,F14.6)') RinO(j,k),Obsin(j,k)
      end do
      if (l84 /= ' ') write(u6,'(a)') trim(l84)
    end do
    l84 = ' '
    l84x = ' '
    do k=i,min(i+3,iobs)
      if (iplot(i) /= 0) l84(24*(k-i)+1:24*(k-i)+14) = '     Plot file'
      if (iplot(i) == 0) l84(24*(k-i)+1:24*(k-i)+17) = '     No plot file'
      write(l84x(24*(k-i)+1:24*(k-i+1)),'(F10.4,F14.6)') Req,EoutO(ngrid+1,k)
    end do
    write(u6,'(a)') trim(l84)
    write(u6,*) ' Interpolated value at equilibrium for state 1:'
    write(u6,'(a)') trim(l84x)
  end do
end if

return

1100 format(/1x,'Vibration-Rotation spectrum for the ',a4,' molecule.')
1101 format(/1x,'Transition moments for the ',a4,' molecule.')
1190 format(/1x,'State number 1')
1191 format(/1x,'State number 2')
1200 format(/1x,'Rotational quantum number range ',2i3 &
            /1x,'Electronic angular momentum     ',i3 &
            /1x,'Vibrational quantum number range',2i3 &
            /1x,'Mass of atom ',i3,1x,a2,f15.6,' au' &
            /1x,'Mass of atom ',i3,1x,a2,f15.6,' au' &
            /1x,'Reduced mass       ',f15.6,' au')
1201 format(/1x,'Rotational quantum number range ',2i3 &
            /1x,'Electronic angular momentum     ',i3 &
            /1x,'Vibrational quantum number range',2i3 &
            /1x,'Mass of atom ',i3,1x,a2,f15.6,' au' &
            /1x,'Mass of atom ',i3,1x,a2,f15.6,' au' &
            /1x,'Reduced mass       ',f15.6,' au')
1300 format(/1x,'Statistics for numerical integration' &
            /1x,'Number of steps         ',i5 &
            /1x,'Step length             ',f10.6 &
            /1x,'Radial integration range',2f10.6,' au' &
            /1x,'logarithmic range       ',2f10.6)
1400 format(/1x,'Spectroscopic constants will not be computed')
1500 format(/1x,'Spectroscopic constants will be computed')
1600 format(1x,'Matrix elements of operators will not be computed')
1700 format(1x,'Matrix elements of',i3,' operators will be computed')
1800 format(/1x,'Potential read from input'/6x,' R(au)     E(au)')
1810 format(1x,2f12.6)
1820 format(1x,'A plot file of the potential will be generated')
1900 format(/1x,'Input data for observables')
1910 format(/1x,4a24)

end subroutine Vibinp
