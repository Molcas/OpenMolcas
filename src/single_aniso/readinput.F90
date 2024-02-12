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

subroutine readin_single(iprint,nmult,ndim,ldim,ndimcf,ldimcf,nlanth,axisoption,poly_file,Ifrestart,input_to_read,nk,mg,zmagn, &
                         Do_structure_abc,cryst,coord,encut_definition,compute_g_tensors,compute_CF,nDirTot,nss,nstate, &
                         compute_magnetization,compute_torque,smagn,tinput,hinput,compute_Mdir_vector,zeeman_energy,LUZee,doplot, &
                         encut_rate,ncut,nTempMagn,TempMagn,m_paranoid,compute_barrier,nBlock,AngPoints,input_file_name,nT,nH, &
                         texp,chit_exp,zJ,hexp,magn_exp,hmin,hmax,nDir,nDirZee,dirX,dirY,dirZ,dir_weight,xfield,tmin,tmax,thrs, &
                         H_torq,T_torq)
! THIS ROUTINE READS THE FILE "SINGLE_ANISO.INPUT".

implicit none
integer, parameter :: wp = kind(0.d0)
#include "warnings.h"
#include "mgrid.fh"
!-----------------------------------------------------------------------
! magnetization vectors:
integer :: nDir, nDirZee
real(kind=8) :: dirX(nDir), dirY(nDir), dirZ(nDir)
real(kind=8) :: dir_weight(nDirZee,3)
logical :: compute_Mdir_vector, zeeman_energy
!common/MVL/ compute_Mdir_vector
!common/MZEL/ zeeman_energy
!-----------------------------------------------------------------------
integer :: nss, nstate
integer :: iprint, nt, nh, nk, mg, l, jEnd
integer :: nlanth, ndimcf, ldimcf, axisoption, i_OxStat
integer :: input_to_read, encut_definition, ncut, ntempmagn
integer :: ndirtot
integer :: nBlock
integer :: nmult, ndim(nMult), ldim
integer :: AngPoints
integer :: LUZee(nDirZee)
real(kind=8) :: tmin, tmax, hmin, hmax, t1, t2, zj, tempmagn(nTempMagn), encut_rate
real(kind=8) :: texp(nT), chit_exp(nT)
real(kind=8) :: hexp(nH), magn_exp(nH,ntempmagn)
real(kind=8) :: zmagn(3,3), sum, tmp
real(kind=8) :: cryst(6), coord(3)
real(kind=8) :: column_check(3,3), row_check(3,3)
real(kind=8) :: check_dir_weight(nDirZee)
real(kind=8) :: zr(3,3), det_zmagn
real(kind=8) :: FindDetR
real(kind=8) :: Xfield
real(kind=8), intent(out) :: thrs
real(kind=8), intent(out) :: H_torq, T_torq
logical :: hcheck, tcheck, poly_file
logical :: Ifrestart, Do_structure_abc
logical :: compute_magnetization, encut_check, compute_cf
logical :: compute_g_tensors
logical :: checktmag
logical :: compute_barrier
logical :: compute_torque
logical :: smagn
logical :: tinput, hinput
logical, intent(out) :: m_paranoid
logical :: doplot
character(len=2) :: cME, clanth(37)
character(len=21) :: namefile_energy
character(len=180) :: input_file_name, tmpline, err_msg
external :: FindDetR
integer :: IsFreeUnit
external :: IsFreeUnit
!=======================================================================
!common/CHISUBR/ TMIN,TMAX,T1,T2
!common/CHISUBL/ TINPUT
!common/MAGNSUBI/ NK,MG
!common/MAGNSUBR/ HMIN,HMAX
!common/MAGNSUBL/ HINPUT
integer :: I, LINENR, j
character(len=280) :: LINE
logical :: DBG

DBG = .false.
!============ Some default settings=====================================
! variables in "mgrid.fh"
do i=1,32
  do j=1,3
    get_nP(j,i) = 0
  end do
end do
nsymm = 1
ngrid = 15
get_nP(1,1) = 5
get_nP(1,2) = 9
get_nP(1,3) = 17
get_nP(1,4) = 25
get_nP(1,5) = 29
get_nP(1,6) = 45
get_nP(1,7) = 49
get_nP(1,8) = 61
get_nP(1,9) = 77
get_nP(1,10) = 93
get_nP(1,11) = 105
get_nP(1,12) = 125
get_nP(1,13) = 141
get_nP(1,14) = 161
get_nP(1,15) = 185
get_nP(1,16) = 229
get_nP(1,17) = 309
get_nP(1,18) = 401
get_nP(1,19) = 505
get_nP(1,20) = 621
get_nP(1,21) = 749
get_nP(1,22) = 889
get_nP(1,23) = 1041
get_nP(1,24) = 1205
get_nP(1,25) = 1381
get_nP(1,26) = 1569
get_nP(1,27) = 1769
get_nP(1,28) = 1981
get_nP(1,29) = 2205
get_nP(1,30) = 2441
get_nP(1,31) = 2689
get_nP(1,32) = 2949
get_nP(2,1) = 4
get_nP(2,2) = 6
get_nP(2,3) = 11
get_nP(2,4) = 16
get_nP(2,5) = 17
get_nP(2,6) = 27
get_nP(2,7) = 28
get_nP(2,8) = 34
get_nP(2,9) = 41
get_nP(2,10) = 51
get_nP(2,11) = 57
get_nP(2,12) = 68
get_nP(2,13) = 75
get_nP(2,14) = 86
get_nP(2,15) = 98
get_nP(2,16) = 121
get_nP(2,17) = 162
get_nP(2,18) = 209
get_nP(2,19) = 262
get_nP(2,20) = 321
get_nP(2,21) = 386
get_nP(2,22) = 457
get_nP(2,23) = 534
get_nP(2,24) = 617
get_nP(2,25) = 706
get_nP(2,26) = 801
get_nP(2,27) = 902
get_nP(2,28) = 1009
get_nP(2,29) = 1122
get_nP(2,30) = 1241
get_nP(2,31) = 1366
get_nP(2,32) = 1497
get_nP(3,1) = 3
get_nP(3,2) = 4
get_nP(3,3) = 7
get_nP(3,4) = 10
get_nP(3,5) = 10
get_nP(3,6) = 16
get_nP(3,7) = 16
get_nP(3,8) = 19
get_nP(3,9) = 22
get_nP(3,10) = 28
get_nP(3,11) = 31
get_nP(3,12) = 37
get_nP(3,13) = 40
get_nP(3,14) = 46
get_nP(3,15) = 52
get_nP(3,16) = 64
get_nP(3,17) = 85
get_nP(3,18) = 109
get_nP(3,19) = 136
get_nP(3,20) = 166
get_nP(3,21) = 199
get_nP(3,22) = 235
get_nP(3,23) = 274
get_nP(3,24) = 316
get_nP(3,25) = 361
get_nP(3,26) = 409
get_nP(3,27) = 460
get_nP(3,28) = 514
get_nP(3,29) = 571
get_nP(3,30) = 631
get_nP(3,31) = 694
get_nP(3,32) = 760
compute_Mdir_vector = .false.
zeeman_energy = .false.
do i=1,nDir
  DirX(i) = 0.0_wp
  DirY(i) = 0.0_wp
  DirZ(i) = 0.0_wp
end do
do i=1,nDirZee
  dir_weight(i,1) = 0.0_wp
  dir_weight(i,2) = 0.0_wp
  dir_weight(i,3) = 0.0_wp
end do
!========== Initializations of arrays ==================================
do I=1,nTempMagn
  TempMagn(i) = 0.0_wp
end do
!============ Initializations of constants =============================

thrs = 1.0D-10
ldim = 1
ncut = 1
TMIN = 0.0_wp
TMAX = 300.0_wp
XFIELD = 0.0_wp
!NT = 301
HMIN = 0.0_wp
HMAX = 10.0_wp
!NH = 21
NK = 200
MG = 200
!NDIR = 0
!nDirZee = 0
!nTempMagn = 1
if (nTempMagn > 0) TempMagn(1) = 2.0_wp
T1 = 5.0_wp
T2 = 6.0_wp
ZJ = 0.0_wp
m_paranoid = .true.
checkTMAG = .false.
compute_g_tensors = .false.
compute_magnetization = .false.
compute_Mdir_vector = .false.
compute_CF = .false.
compute_barrier = .false.
compute_torque = .false.
smagn = .false.
TINPUT = .false.
HINPUT = .false.
TCHECK = .false.
HCHECK = .false.
POLY_FILE = .false.
Do_structure_abc = .false.
zeeman_energy = .false.
ENCUT_check = .false.
doplot = .false.
nlanth = 0
nDIMcf = 0
cME = '  '
! -- lanthanides
clanth(1) = 'CE'
clanth(2) = 'PR'
clanth(3) = 'ND'
clanth(4) = 'PM'
clanth(5) = 'SM'
clanth(6) = 'EU'
clanth(7) = 'GD'
clanth(8) = 'TB'
clanth(9) = 'DY'
clanth(10) = 'HO'
clanth(11) = 'ER'
clanth(12) = 'TM'
clanth(13) = 'YB'
clanth(14) = 'LU'
! -- actinides
clanth(15) = 'TH'
clanth(16) = 'PA'
clanth(17) = 'U'
clanth(18) = 'NP'
clanth(19) = 'PU'
clanth(20) = 'AM'
clanth(21) = 'CM'
clanth(22) = 'BK'
clanth(23) = 'CF'
clanth(24) = 'ES'
clanth(25) = 'FM'
clanth(26) = 'MD'
clanth(27) = 'NO'
clanth(28) = 'LR'
! -- transition metals
clanth(29) = 'SC'
clanth(30) = 'TI'
clanth(31) = 'V'
clanth(32) = 'CR'
clanth(33) = 'MN'
clanth(34) = 'FE'
clanth(35) = 'CO'
clanth(36) = 'NI'
clanth(37) = 'CU'

do i=1,6
  cryst(i) = 0.0_wp
end do
do i=1,3
  coord(i) = 0.0_wp
end do
axisoption = 1
input_to_read = 0
encut_definition = 2
encut_rate = 1
do i=1,3
  coord(i) = 0.0_wp
  do j=1,3
    zmagn(i,j) = 0.0_wp
  end do
end do
AngPoints = 46

!=========== End of default settings====================================
rewind(5)
50 read(5,'(A280)',end=998) LINE
if (DBG) write(6,'(A)') trim(LINE)
call NORMAL(LINE)
if (LINE(1:7) /= '&SINGLE') go to 50
LINENR = 0
100 read(5,'(A280)',end=998) LINE
if (DBG) write(6,'(A)') trim(LINE)
LINENR = LINENR+1
call NORMAL(LINE)
if (LINE(1:1) == '*') go to 100
if (LINE == ' ') go to 100
if ((LINE(1:4) == 'End ') .or. (LINE(1:4) == '    ')) go to 200

!------------------------------------------
!if (LINE(1:4) == 'TYPE') then
!  read(5,*,err=997) ICALC
!  if (icalc == 1) then
!    compute_g_tensors = .true.
!  else if (icalc == 2) then
!    compute_chiT = .true.
!  else if (icalc == 3) then
!    compute_magnetization = .true.
!  else if (icalc == 4) then
!    compute_g_tensors = .true.
!    compute_chiT = .true.
!  else if (icalc == 5) then
!    compute_g_tensors = .true.
!    compute_magnetization = .true.
!  else if (icalc == 6) then
!    compute_chiT = .true.
!    compute_magnetization = .true.
!  else if (icalc == 7) then
!    compute_g_tensors = .true.
!    compute_chiT = .true.
!    compute_magnetization = .true.
!  else
!    write(6,'(A)') 'ICALC: the maximum value is 7. However, the calculation will continue by computing the magnetism.'
!    compute_g_tensors = .true.
!    compute_chiT = .true.
!    compute_magnetization = .true.
!  end if
!  LINENR = LINENR+1
!  go to 100
!end if
!------------------------------------------
if (LINE(1:4) == 'MLTP') then
  read(5,*,err=997) NMULT
  if (DBG) write(6,*) 'MLTP:  NMULT=',NMULT
  compute_g_tensors = .true.
  read(5,*,err=997) (NDIM(i),i=1,NMULT)
  if (DBG) write(6,*) 'MLTP: NDIM()=',(NDIM(i),i=1,NMULT)
  LINENR = LINENR+2
  go to 100
end if
!------------------------------------------
if (LINE(1:4) == 'REST') then
  Ifrestart = .true.
  read(5,*,err=997) input_to_read
  if (DBG) write(6,*) 'REST: input_to_read=',input_to_read
  if ((input_to_read == 2) .or. (input_to_read == 3) .or. (input_to_read == 4)) then
    backspace(5)
    read(5,*) input_to_read,tmpline
    input_file_name = trim(tmpline)
  end if
  if (input_to_read == 1) then
    write(6,*) 'RESTART: -- The SINGLE_ANISO will take all ab initio information from the binary "$Project.aniso" file.'
  else if (input_to_read == 2) then
    write(6,*) 'RESTART: -- The SINGLE_ANISO will take all ab initio information from the ASCII '//trim(input_file_name)//' file.'
  else if (input_to_read == 3) then
    write(6,*) 'RESTART: -- The SINGLE_ANISO will take all ab initio information from the RASSI-HDF5 binary file.'
  else if (input_to_read == 4) then
    write(6,*) 'RESTART: -- The SINGLE_ANISO will take all ab initio information from the ASCII '//trim(input_file_name)// &
               ' file -- molcas-8.0 format.'
  else
    call WarningMessage(2,'SINGLE_ANISO:: RESTART  option is not known.')
    call Quit_OnUserError()
  end if
  go to 100
end if

!-------------------------------------------

if (line(1:4) == 'DATA') then
  Ifrestart = .true.
  read(5,*) tmpline
  input_file_name = trim(tmpline)
  input_to_read = 6
  if (DBG) write(6,*) 'restart_check: DATA, input_file_name='
  if (DBG) write(6,*) input_file_name
  LINENR = LINENR+1
  go to 100
end if

!-------------------------------------------
if (LINE(1:4) == 'TINT') then
  if (.not. TINPUT) then
    TCHECK = .true.

    t1 = 0.0_wp
    t2 = 0.0_wp

    read(5,*,err=997) t1,t2,nT

    if ((t1 < 0) .or. (t2 < 0)) then
      call WarningMessage(2,'TINT: negative temperature requested! ')
      call Quit_OnUserError()
    end if
    if ((t1-t2) > 0.0_wp) then
      Tmin = t2
      Tmax = t1
    else if ((t1-t2) < 0.0_wp) then
      Tmin = t1
      Tmax = t2
    else ! t1==t2
      call WarningMessage(2,'TINT: temperature interval == 0! ')
      call Quit_OnUserError()
    end if

    if (DBG) write(6,*) 'TINT: Tmin, Tmax, nT=',Tmin,Tmax,nT
  else
    goto 590
  end if
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'XFIE') then
  read(5,*,err=997) Xfield
  if (DBG) write(6,*) 'XFIE: Xfield=',Xfield
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'HINT') then
  if (.not. HINPUT) then
    HCHECK = .true.
    compute_magnetization = .true.

    t1 = 0.0_wp
    t2 = 0.0_wp

    read(5,*,err=997) t1,t2,nH

    if ((t1 < 0) .or. (t2 < 0)) then
      call WarningMessage(2,'HINT: negative field requested! ')
      call Quit_OnUserError()
    end if

    if ((t1-t2) > 0.0_wp) then
      Hmin = t2
      Hmax = t1
    else if ((t1-t2) < 0.0_wp) then
      Hmin = t1
      Hmax = t2
    else ! t1 == t2
      call WarningMessage(2,'HINT: temperature interval == 0! ')
      call Quit_OnUserError()
    end if

    if (DBG) write(6,*) 'HINT: Hmin, Hmax, nH=',Hmin,Hmax,nH
  else
    go to 591
  end if
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'NCUT') then
  if (ENCUT_check) then
    go to 595
  else
    ENCUT_check = .true.
    compute_magnetization = .true. ! request for computation of M(H)
    encut_definition = 1

    read(5,*,err=997) NCUT  ! E_cut=ESO(Ncut)

    if (NCUT < 0) then
      call WarningMessage(2,'NCUT: negative NCUT requested! ')
      call Quit_OnUserError()
    else if (NCUT == 0) then
      call WarningMessage(2,'NCUT: zero NCUT requested! ')
      call Quit_OnUserError()
    end if

    if (DBG) write(6,*) 'NCUT: NCUT=',NCUT
    LINENR = LINENR+1
    go to 100
  end if
end if
!-------------------------------------------
if (LINE(1:4) == 'ENCU') then
  if (ENCUT_check) then
    go to 595
  else
    ENCUT_check = .true.
    compute_magnetization = .true.
    encut_definition = 2

    read(5,*,err=997) NK,MG

    if ((NK <= 0) .or. (MG <= 0)) then
      call WarningMessage(2,'ENCU: zero or negative NK,MG requested! ')
      call Quit_OnUserError()
    end if

    if (DBG) write(6,*) 'ENCU: NK, MG=',NK,MG
    LINENR = LINENR+1
    go to 100
  end if
end if
!-------------------------------------------
if (LINE(1:4) == 'ERAT') then
  if (ENCUT_check) then
    go to 595
  else
    ENCUT_check = .true.
    compute_magnetization = .true.
    encut_definition = 3
    ! Ncut = INT(nss*encut_rate)
    ! E_cut = E(Ncut)

    read(5,*,err=997) encut_rate

    if (encut_rate <= 0.0_wp) then
      call WarningMessage(2,'ERAT: zero or negative encut rate requested! ')
      call Quit_OnUserError()
    end if

    if (DBG) write(6,*) 'ERAT: encut_rate=',encut_rate
    LINENR = LINENR+1
    go to 100
  end if
end if
!-------------------------------------------
if (LINE(1:4) == 'MVEC') then
  compute_magnetization = .true.   ! request for computation of M(H)
  compute_Mdir_vector = .true.
  read(5,*,err=997) nDir
  if (DBG) write(6,*) 'MVEC: nDir=',nDir
  do i=1,nDir
    read(5,*,err=997) DirX(i),DirY(i),DirZ(i)
    if (DBG) write(6,*) 'MVEC: DirX,DirY,DirZ=',DirX(i),DirY(i),DirZ(i)
  end do
  ! some processing:
  do i=1,nDir
    sum = 0.0_wp
    sum = DirX(i)*DirX(i)+DirY(i)*DirY(i)+DirZ(i)*DirZ(i)
    if (sum == 0.0_wp) then
      write(err_msg,'(a,i3,a)') 'error: MVEC  vector ',i,'has the modulus = 0.0_wp.'
      call WarningMessage(2,err_msg)
      call Quit_OnUserError()
    end if
    if (sum /= 1.0_wp) then
      write(6,'(a,i3,a)') 'the vector',i,'was re-normalized.'
      tmp = dirX(i)/sqrt(sum)
      dirX(i) = tmp
      tmp = dirY(i)/sqrt(sum)
      dirY(i) = tmp
      tmp = dirZ(i)/sqrt(sum)
      dirZ(i) = tmp
    end if
  end do

  LINENR = LINENR+NDIR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'MAVE') then
  compute_magnetization = .true.

  read(5,*,err=997) nsymm,ngrid

  if (DBG) write(6,*) 'MAVE: nsymm, ngrid=',nsymm,ngrid
  if ((nsymm < 1) .or. (nsymm > 3)) then
    write(6,'(A)') '"nsymm" must take Integer values 1, 2 or 3.'
    write(6,'(A,i5)') '"nsymm" = ',nsymm
    call Quit_OnUserError()
  end if
  if ((ngrid < 1) .or. (ngrid > 32)) then
    write(6,'(A)') '"ngrid" must take Integer values 1, 2, ... 32.'
    write(6,'(A,i5)') '"ngrid" = ',ngrid
    call Quit_OnUserError()
  end if
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
!if (LINE(1:4) == 'TLIN') then
!  read(5,*,err=997) T1,T2
!  LINENR = LINENR+1
!  go to 100
!end if
!-------------------------------------------
if (LINE(1:4) == 'SMAG') then
  smagn = .true.
  if (DBG) write(6,*) 'SMAG: =',smagn
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'PLOT') then
  doplot = .true.
  if (DBG) write(6,*) 'PLOT: =',doplot
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'TEXP') then
  if (.not. TCHECK) then
    TINPUT = .true.
    read(5,*,err=997) NT
    if (DBG) write(6,*) 'TEXP: nT=',nT
    do i=1,NT
      texp(i) = 0.0_wp
      chit_exp(i) = 0.0_wp
      read(5,*,err=997) texp(i),chit_exp(i)
      if (DBG) write(6,*) 'TEXP: texp(i), chit_exp(i)=',texp(i),chit_exp(i)
      ! check and clean negative values:
      if (texp(i) < 0.0_wp) texp(i) = abs(texp(i))
      if (chit_exp(i) < 0.0_wp) chit_exp(i) = abs(chit_exp(i))
    end do
    tmin = texp(1)
    tmax = texp(nT)
  else
    go to 590
  end if
  LINENR = LINENR+NT+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'HEXP') then
  compute_magnetization = .true.
  if (checkTMAG) write(6,'(A)') 'The data provided in TMAG will be ignored.'
  if (.not. HCHECK) then
    HINPUT = .true.
    read(5,*) nTempMagn,(TempMagn(i),i=1,nTempMagn)
    if (DBG) write(6,*) 'HEXP: nTempMagn =',nTempMagn
    if (DBG) write(6,*) 'HEXP: TempMagn()=',(TempMagn(i),i=1,nTempMagn)
    read(5,*) nH
    if (DBG) write(6,*) 'HEXP: nH =',nH
    if (nH < 0) nH = abs(nH)
    if (nH == 0) call Quit_OnUserError()
    do i=1,nH
      hexp(i) = 0.0_wp
      do j=1,nTempMagn
        magn_exp(i,j) = 0.0_wp
      end do
    end do
    do i=1,nH
      read(5,*,err=997) Hexp(i),(magn_exp(i,j),j=1,nTempMagn)
      if (DBG) write(6,*) 'HEXP: Hexp(i),  magn_exp(i,j)=',Hexp(i),(magn_exp(i,j),j=1,nTempMagn)
      ! check and clean negative values:
      if (hexp(i) < 0.0_wp) hexp(i) = abs(hexp(i))
      do j=1,nTempMagn
        if (magn_exp(i,j) < 0.0_wp) magn_exp(i,j) = abs(magn_exp(i,j))
      end do
    end do
    hmin = hexp(1)
    hmax = hexp(nH)
  else
    go to 591
  end if
  LINENR = LINENR+NH+2
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'ZJPR') then
  read(5,*,err=997) ZJ
  if (DBG) write(6,*) 'ZJPR: zJ =',zJ
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'TORQ') then
  compute_torque = .true.
  read(5,*,err=997) AngPoints,H_torq,T_torq
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'TMAG') then
  if (.not. HINPUT) then
    compute_magnetization = .true.
    checkTMAG = .true.

    read(5,*,err=997) nTempMagn,(TempMagn(i),i=1,nTempMagn)

    do i=1,nTempMagn
      if (TempMagn(i) <= 0.0_wp) then
        call WarningMessage(2,'TMAG: zero or negative temperature requested! ')
        if (TempMagn(i) < 0.0_wp) TempMagn(i) = abs(TempMagn(i))
        if (TempMagn(i) == 0.0_wp) TempMagn(i) = 0.0001_wp
      end if
    end do

    if (DBG) write(6,*) 'TMAG: nTempMagn =',nTempMagn
    if (DBG) write(6,*) 'TMAG: TempMagn()=',(TempMagn(i),i=1,nTempMagn)
    ! check and clean negative values:
  else
    write(6,'(A)') 'TMAG data is taken from HEXP.'
  end if
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'PRLV') then
  read(5,*,err=997) IPRINT
  if (DBG) write(6,*) 'PRLV: IPRINT =',iPrint
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'POLY') then
  if (DBG) write(6,*) 'POLY:'
  POLY_FILE = .true.
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'CRYS') then
  compute_CF = .true.
  read(5,*,err=997) cME
  if (DBG) write(6,*) 'CRYS: cME =',cME

  if ((cME == 'ce') .or. (cME == 'Ce') .or. (cME == 'cE') .or. (cME == 'CE')) then
    nlanth = 1
    nDIMcf = 6  ! f1; multiplet J=L-S=3-1/2=5/2  =>  J = 2F_5/2
    lDIMCF = 7  ! (L=3)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'pr') .or. (cME == 'Pr') .or. (cME == 'pR') .or. (cME == 'PR')) then
    nlanth = 2
    nDIMcf = 9  ! f2; multiplet J=L-S=5-1=4  => J = 3H_4
    lDIMCF = 11 ! (L=5)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'nd') .or. (cME == 'Nd') .or. (cME == 'nD') .or. (cME == 'ND')) then
    nlanth = 3
    nDIMcf = 10 ! f3; multiplet J=L-S=6-3/2=9/2  => J = 4I_9/2
    lDIMCF = 13 ! (L=6)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'pm') .or. (cME == 'Pm') .or. (cME == 'pM') .or. (cME == 'PM')) then
    nlanth = 4
    nDIMcf = 9  ! f4; multiplet J=L-S=6-2=4  => J = 5I_4
    lDIMCF = 13 ! (L=6)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'sm') .or. (cME == 'Sm') .or. (cME == 'sM') .or. (cME == 'SM')) then
    nlanth = 5
    nDIMcf = 6  ! f5; multiplet J=L-S=5-5/2=5/2  => J = 6H_5/2
    lDIMCF = 11 ! (L=5)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'eu') .or. (cME == 'Eu') .or. (cME == 'eU') .or. (cME == 'EU')) then
    nlanth = 6
    nDIMcf = 1  ! f6; multiplet J=L-S=3-3=0  => J = 3F_0
    lDIMCF = 7  ! (L=3)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'gd') .or. (cME == 'Gd') .or. (cME == 'gD') .or. (cME == 'GD')) then
    nlanth = 7
    nDIMcf = 8  ! f7; multiplet J=L+S=0+7/2=0  => J = 8S_7/2
    lDIMCF = 1  ! (L=0)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'tb') .or. (cME == 'Tb') .or. (cME == 'tB') .or. (cME == 'TB')) then
    nlanth = 8
    nDIMcf = 13 ! f8; multiplet J=L+S=3+3=0  => J = 7F_6
    lDIMCF = 7  ! (L=3)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'dy') .or. (cME == 'Dy') .or. (cME == 'dY') .or. (cME == 'DY')) then
    nlanth = 9
    nDIMcf = 16 ! f9; multiplet J=L+S=5+5/2=15/2  => J = 6H_15/2
    lDIMCF = 11 ! (L=5)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'ho') .or. (cME == 'Ho') .or. (cME == 'hO') .or. (cME == 'HO')) then
    nlanth = 10
    nDIMcf = 17 ! f10; multiplet J=L+S=6+2=8  => J = 5I_8
    lDIMCF = 13 ! (L=6)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'er') .or. (cME == 'Er') .or. (cME == 'eR') .or. (cME == 'ER')) then
    nlanth = 11
    nDIMcf = 16 ! f11; multiplet J=L+S=6+3/2=15/2  => J = 4I_15/2
    lDIMCF = 13 ! (L=6)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'tm') .or. (cME == 'Tm') .or. (cME == 'tM') .or. (cME == 'TM')) then
    nlanth = 12
    nDIMcf = 13 ! f12; multiplet J=L+S=5+1=6  => J = 3H_6
    lDIMCF = 11 ! (L=5)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'yb') .or. (cME == 'Yb') .or. (cME == 'yB') .or. (cME == 'YB')) then
    nlanth = 13
    nDIMcf = 8  ! f13; multiplet J=L+S=3+1/2=7/2  => J = 2F_7/2
    lDIMCF = 7  ! (L=3)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'lu') .or. (cME == 'Lu') .or. (cME == 'lU') .or. (cME == 'LU')) then
    nlanth = 14
    nDIMcf = 1  ! f14; multiplet J=L+S=0+0=0  => J = 1S_0
    lDIMCF = 1  ! (L=0)

    !- - - - - - - - - - - - - - - - - - - -
    ! ACTINIDES
  else if ((cME == 'th') .or. (cME == 'Th') .or. (cME == 'tH') .or. (cME == 'TH')) then
    nlanth = 15
    nDIMcf = 6  ! f1; multiplet J=L-S=3-1/2=5/2  =>  J = 2F_5/2
    lDIMCF = 7 ! (L=3)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'pa') .or. (cME == 'Pa') .or. (cME == 'pA') .or. (cME == 'PA')) then
    nlanth = 16
    nDIMcf = 9  ! f2; multiplet J=L-S=5-1=4  => J = 3H_4
    lDIMCF = 11 ! (L=5)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'U') .or. (cME == 'u') .or. (cME == 'u ') .or. (cME == 'U ')) then
    nlanth = 17
    nDIMcf = 10  ! f3; multiplet J=L-S=6-3/2=9/2  => J = 4I_9/2
    lDIMCF = 13 ! (L=6)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'np') .or. (cME == 'Np') .or. (cME == 'nP') .or. (cME == 'NP')) then
    nlanth = 18
    nDIMcf = 9  ! f4; multiplet J=L-S=6-2=4  => J = 5I_4
    lDIMCF = 13 ! (L=6)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'pu') .or. (cME == 'Pu') .or. (cME == 'pU') .or. (cME == 'PU')) then
    nlanth = 19
    nDIMcf = 6  ! f5; multiplet J=L-S=5-5/2=5/2  => J = 6H_5/2
    lDIMCF = 11 ! (L=5)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'am') .or. (cME == 'Am') .or. (cME == 'aM') .or. (cME == 'AM')) then
    nlanth = 20
    nDIMcf = 1  ! f6; multiplet J=L-S=3-3=0  => J = 3F_0
    lDIMCF = 7 ! (L=3)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'cm') .or. (cME == 'Cm') .or. (cME == 'cM') .or. (cME == 'CM')) then
    nlanth = 21
    nDIMcf = 8  ! f7; multiplet J=L+S=0+7/2=0  => J = 8S_7/2
    lDIMCF = 1 ! (L=0)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'bk') .or. (cME == 'Bk') .or. (cME == 'bK') .or. (cME == 'BK')) then
    nlanth = 22
    nDIMcf = 13  ! f8; multiplet J=L+S=3+3=0  => J = 7F_6
    lDIMCF = 7 ! (L=3)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'cf') .or. (cME == 'Cf') .or. (cME == 'cF') .or. (cME == 'CF')) then
    nlanth = 23
    nDIMcf = 16  ! f9; multiplet J=L+S=5+5/2=15/2  => J = 6H_15/2
    lDIMCF = 11 ! (L=5)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'es') .or. (cME == 'Es') .or. (cME == 'eS') .or. (cME == 'ES')) then
    nlanth = 24
    nDIMcf = 17  ! f10; multiplet J=L+S=6+2=8  => J = 5I_8
    lDIMCF = 13 ! (L=6)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'fm') .or. (cME == 'Fm') .or. (cME == 'fM') .or. (cME == 'FM')) then
    nlanth = 25
    nDIMcf = 16  ! f11; multiplet J=L+S=6+3/2=15/2  => J = 4I_15/2
    lDIMCF = 13 ! (L=6)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'md') .or. (cME == 'Md') .or. (cME == 'mD') .or. (cME == 'MD')) then
    nlanth = 26
    nDIMcf = 13  ! f12; multiplet J=L+S=5+1=6  => J = 3H_6
    lDIMCF = 11 ! (L=5)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'no') .or. (cME == 'No') .or. (cME == 'nO') .or. (cME == 'NO')) then
    nlanth = 27
    nDIMcf = 8  ! f13; multiplet J=L+S=3+1/2=7/2  => J = 2F_7/2
    lDIMCF = 7 ! (L=3)
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'lr') .or. (cME == 'Lr') .or. (cME == 'lR') .or. (cME == 'LR')) then
    nlanth = 28
    nDIMcf = 1  ! f14; multiplet J=L+S=0+0=0  => J = 1S_0
    lDIMCF = 1 ! (L=0)

    !--------------------- transition metals --------------------------!

  else if ((cME == 'Sc') .or. (cME == 'Sc') .or. (cME == 'sC') .or. (cME == 'SC')) then

    nlanth = 29
    ! Sc2+ -- d^1
    read(5,*,err=997) i_OxStat
    if (DBG) write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

    if (i_OxStat < 0) then
      write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
      write(6,'(A)') 'It was re-set to positive.'
      i_OxStat = abs(i_OxStat)
    end if
    if (i_OxStat < 2) then
      lDIMCF = 1 ! (L=0)
      write(6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    else if (i_OxStat == 2) then
      lDIMCF = 5 ! (L=2) d^1
    else
      lDIMCF = 1 ! (L=0)
      write(6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    end if

    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'Ti') .or. (cME == 'Ti') .or. (cME == 'tI') .or. (cME == 'TI')) then

    nlanth = 30
    ! Ti2+ -- d^2
    ! Ti3+ -- d^1
    ! Ti4+ -- d^0
    read(5,*,err=997) i_OxStat
    if (DBG) write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

    if (i_OxStat < 0) then
      write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
      write(6,'(A)') 'It was re-set to positive.'
      i_OxStat = abs(i_OxStat)
    end if
    if (i_OxStat < 2) then
      lDIMCF = 1 ! (L=0)
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    else if (i_OxStat == 2) then
      lDIMCF = 7 ! (L=3) d^2  3F
    else if (i_OxStat == 3) then
      lDIMCF = 5 ! (L=2) d^1  2D
    else
      lDIMCF = 1 ! (L=0) d^4
      write(6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    end if

    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'v') .or. (cME == 'V') .or. (cME == 'V ') .or. (cME == 'v ')) then

    nlanth = 31
    ! V2+ -- d^3
    ! V3+ -- d^2
    ! V4+ -- d^1
    read(5,*,err=997) i_OxStat
    if (DBG) write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

    if (i_OxStat < 0) then
      write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
      write(6,'(A)') 'It was re-set to positive.'
      i_OxStat = abs(i_OxStat)
    end if
    if (i_OxStat < 2) then
      lDIMCF = 1 ! (L=0)
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    else if (i_OxStat == 2) then
      lDIMCF = 7 ! (L=3) d^3
    else if (i_OxStat == 3) then
      lDIMCF = 7 ! (L=3) d^2
    else if (i_OxStat == 4) then
      lDIMCF = 5 ! (L=2) d^1
    else
      lDIMCF = 1 ! (L=0) d^4
      write(6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    end if

    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'cr') .or. (cME == 'cR') .or. (cME == 'Cr') .or. (cME == 'CR')) then

    nlanth = 32
    ! Cr3+ -- d^4
    read(5,*,err=997) i_OxStat

    if (i_OxStat < 0) then
      write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
      write(6,'(A)') 'It was re-set to positive.'
      i_OxStat = abs(i_OxStat)
    end if
    if (i_OxStat == 2) then
      lDIMCF = 5 ! (L=2) d^4
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
    else if (i_OxStat == 3) then
      lDIMCF = 7 ! (L=3) d^3
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
    else if (i_OxStat == 4) then
      lDIMCF = 7 ! (L=3) d^2
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
    else if (i_OxStat == 5) then
      lDIMCF = 5 ! (L=2) d^1
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
    else
      lDIMCF = 1 ! (L=0)
      write(6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    end if
    write(6,'(A)') 'Crystal field will not be computed'
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'mn') .or. (cME == 'mN') .or. (cME == 'Mn') .or. (cME == 'MN')) then

    nlanth = 33
    ! Mn3+ -- d^4
    read(5,*,err=997) i_OxStat
    if (DBG) write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

    if (i_OxStat < 0) then
      write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
      write(6,'(A)') 'It was re-set to positive.'
      i_OxStat = abs(i_OxStat)
    end if
    if (i_OxStat < 3) then
      lDIMCF = 1 ! (L=0)
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    else if (i_OxStat == 3) then
      lDIMCF = 5 ! (L=2) d^4
    else
      lDIMCF = 1 ! (L=0)
      write(6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    end if
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'fe') .or. (cME == 'fE') .or. (cME == 'Fe') .or. (cME == 'FE')) then

    nlanth = 34
    ! Co2+ -- d^6 or d^4
    read(5,*,err=997) i_OxStat
    if (DBG) write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

    if (i_OxStat < 0) then
      write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
      write(6,'(A)') 'It was re-set to positive.'
      i_OxStat = abs(i_OxStat)
    end if
    if (i_OxStat < 2) then
      lDIMCF = 1 ! (L=0)
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    else if (i_OxStat == 2) then
      lDIMCF = 5 ! (L=2)  d^6  or  d^4
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
    else if (i_OxStat == 3) then
      lDIMCF = 1 ! (L=2)  d^5
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    else if (i_OxStat == 4) then
      lDIMCF = 1 ! (L=2)
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    else
      lDIMCF = 1 ! (L=0)
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    end if
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'co') .or. (cME == 'cO') .or. (cME == 'Co') .or. (cME == 'CO')) then

    nlanth = 35
    ! Co2+ -- d^7
    read(5,*,err=997) i_OxStat
    if (DBG) write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

    if (i_OxStat < 0) then
      write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
      write(6,'(A)') 'It was re-set to positive.'
      i_OxStat = abs(i_OxStat)
    end if
    if (i_OxStat < 2) then
      lDIMCF = 1 ! (L=0)
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    else if (i_OxStat == 2) then
      lDIMCF = 7 ! (L=3) d^7
    else
      lDIMCF = 1 ! (L=0)
      write(6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    end if
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'ni') .or. (cME == 'nI') .or. (cME == 'Ni') .or. (cME == 'NI')) then

    nlanth = 36
    ! Ni2+ -- d^8
    read(5,*,err=997) i_OxStat
    if (DBG) write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

    if (i_OxStat < 0) then
      write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
      write(6,'(A)') 'It was re-set to positive.'
      i_OxStat = abs(i_OxStat)
    end if
    if (i_OxStat < 2) then
      lDIMCF = 1 ! (L=0)
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    else if (i_OxStat == 2) then
      lDIMCF = 7 ! (L=2) d^8
    else
      lDIMCF = 1 ! (L=0)
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    end if
    !- - - - - - - - - - - - - - - - - - - -
  else if ((cME == 'cu') .or. (cME == 'cU') .or. (cME == 'Cu') .or. (cME == 'CU')) then

    nlanth = 37
    ! Cu2+ -- d^9
    read(5,*,err=997) i_OxStat
    if (DBG) write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

    if (i_OxStat < 0) then
      write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
      write(6,'(A)') 'It was re-set to positive.'
      i_OxStat = abs(i_OxStat)
    end if
    if (i_OxStat < 2) then
      lDIMCF = 1 ! (L=0)
      write(6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    else if (i_OxStat == 2) then
      lDIMCF = 5 ! (L=2) d^9
    else
      lDIMCF = 1 ! (L=0) d^4
      write(6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
      write(6,'(A)') 'Crystal field will not be computed'
    end if
    !- - - - - - - - - - - - - - - - - - - -
  else
    write(6,'(A)') 'Label of the metal is not understood.'
    write(6,'(A)') 'Crystal field will not be computed'
  end if

  if (IPRINT > 2) then
    write(6,'(5x,3A)') 'SINGLE_ANISO will calculate the parameters of the crystal field for Ln = ',clanth(nlanth),','
    write(6,'(5x,A,I2,a)') 'for the ground multiplet J. Multiplicity of J = ',nDIMcf,' and'
    write(6,'(5x,A,I2)') 'for the ground LS term. Multiplicity of L = ',lDIMcf
  end if
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'QUAX') then
  !if (check_CRYS) then
  read(5,*,err=997) axisoption
  if (DBG) write(6,*) 'QUAX: axisoption =',axisoption
  LINENR = LINENR+1

  if ((axisoption < 1) .or. (axisoption > 3)) &
    call WarningMessage(2,'QUAX: axisoption out of range! Calculation will continue by employing the default option.')
  if (axisoption == 3) then
    do j=1,3
      read(5,*,err=997) (zmagn(i,j),i=1,3)
      if (DBG) write(6,*) 'QUAX: zmagn(i,j) =',(zmagn(i,j),i=1,3)
    end do
    LINENR = LINENR+3
  end if
  !else
  !  write(6,'(A)') 'The CRYS keyword must be declared above QUAX in the input!'
  !end if
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'PREX') then
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'UBAR') then
  compute_barrier = .true.
  if (DBG) write(6,*) 'UBAR:'
  LINENR = LINENR+1
  go to 100
end if
!-------------------------------------------
if (LINE(1:4) == 'ABCC') then
  Do_structure_abc = .true.
  read(5,*,err=997) (cryst(i),i=1,6)

  do i=1,6
    if (cryst(i) <= 0) then
      call WarningMessage(2,'ABCC: zero or negative crystallographic parameters requested! ')
      call Quit_OnUserError()
    end if
  end do

  if (DBG) write(6,*) 'ABCC: (cryst(i),i=1,6)=',(cryst(i),i=1,6)
  read(5,*,err=997) (coord(i),i=1,3)
  if (DBG) write(6,*) 'ABCC: (coord(i),i=1,3)=',(coord(i),i=1,3)
  LINENR = LINENR+2
  go to 100
end if
! array "cryst" collects the crystallographic data:
!  cryst(1)= a
!  cryst(2)= b
!  cryst(3)= c
!  cryst(4)= alpha
!  cryst(5)= beta
!  cryst(6)= gamma
!  coord(i) =the coordinates of the magnetic center in "abc" axes
!  logical variable 'Do_structure_abc' will make the program compute
!  the magnetic and anisotropy axes in the "abc" coordinate system
!-------------------------------------------
if (LINE(1:4) == 'ZEEM') then
  zeeman_energy = .true.
  compute_magnetization = .true.

  read(5,*,err=997) nDirZee
  if (DBG) write(6,*) 'ZEEM: nDirZee=',nDirZee

  do i=1,nDirZee
    ! open the zeeman_energy_xxx.txt file where Zeeman eigenstates will
    ! be further written in mangetization() subroutine
    write(namefile_energy,'(5A)') 'zeeman_energy_',char(48+mod(int((i)/100),10)),char(48+mod(int((i)/10),10)), &
                                  char(48+mod(int(i),10)),'.txt'
    !print *, 'namefile_energy: ', namefile_energy
    LUZee(i) = IsFreeUnit(30+i)
    call molcas_open(LUZee(i),namefile_energy)
    !open(30+i,file=namefile_energy)

    read(5,*,err=997) (dir_weight(i,l),l=1,3)
    if (DBG) write(6,*) 'ZEEM: (dir_weight(i,l),l=1,3)=',(dir_weight(i,l),l=1,3)

    check_dir_weight(i) = 0.0_wp
    check_dir_weight(i) = sqrt(dir_weight(i,1)**2+dir_weight(i,2)**2+dir_weight(i,3)**2)

    if ((check_dir_weight(i) < 0.995_wp) .or. (check_dir_weight(i) > 1.005_wp)) then
      write(6,'(A)') 'The directions for the magnetic field for the computation of the Zeeman splitting are wrong.'
      write(6,'(A)') '( px^2 + py^2 + pz^2 ) must give 1.!'
      write(6,'(A,I3,2x,A,F9.5)') 'In the present case for direction Nr.',i,' the dir_weight = px^2 + py^2 + pz^2 = ', &
                                  check_dir_weight(i)**2
      LINENR = LINENR+2+i
      go to 997
    end if

  end do
  LINENR = LINENR+nDirZee+1
  go to 100
end if
!-------------------------------------------

200 continue
if (IPRINT > 2) write(6,'(5X,A)') 'NO ERROR WAS LOCATED WHILE READING INPUT'

if (compute_CF) then
  if (axisoption == 3) then
    ! check the determinant of the ZMAGN
    Det_zmagn = 0.0_wp
    do I=1,3
      do J=1,3
        ZR(I,J) = 0.0_wp
        ZR(I,J) = zmagn(I,J)
      end do
    end do
    Det_zmagn = FindDetR(ZR,3)
    if (Det_zmagn < 0.0_wp) then
      write(6,'(A)') 'QUAX: The determinant of the rotation matrix provided in the input is NEGATIVE.'
      write(6,'(A,F22.14)') 'Determinant = ',Det_zmagn
      write(6,'(A)') 'This means that the matrix you have provided can be decomposed in a product of two '
      write(6,'(A)') 'matrices: Rotation*Inversion'
      write(6,'(A)') 'The determinant of the Rotation matrix must be POSITIVE.'
      write(6,'(A)') 'The program will stop.'
      return
    end if

    ! check the orthogonality of the ZMAGN:
    do i=1,3
      do j=1,3
        column_check(i,j) = 0.0_wp
        row_check(i,j) = 0.0_wp
        do L=1,3
          column_check(i,j) = column_check(i,j)+zmagn(i,L)*zmagn(j,L)
          row_check(i,j) = row_check(i,j)+zmagn(L,i)*zmagn(L,j)
        end do
      end do
    end do

    do i=1,3
      do j=i+1,3
        if (i == j) go to 112
        if ((abs(column_check(1,2)) > 0.0001_wp) .or. (abs(column_check(1,3)) > 0.0001_wp) .or. &
            (abs(column_check(2,3)) > 0.0001_wp) .or. (abs(row_check(1,2)) > 0.0001_wp) .or. &
            (abs(row_check(1,3)) > 0.0001_wp) .or. (abs(row_check(2,3)) > 0.0001_wp)) then
          write(6,'(A)') 'QUAX: The rotation matrix is not UNITARY.'
          write(6,'(A,F19.14)') 'column_check(1,2) = ',column_check(1,2)
          write(6,'(A,F19.14)') 'column_check(1,3) = ',column_check(1,3)
          write(6,'(A,F19.14)') 'column_check(2,3) = ',column_check(2,3)
          write(6,'(A,F19.14)') '   row_check(1,2) = ',row_check(1,2)
          write(6,'(A,F19.14)') '   row_check(1,3) = ',row_check(1,3)
          write(6,'(A,F19.14)') '   row_check(2,3) = ',row_check(2,3)
          write(6,'(A)') 'All above values must be exact 0.0.'
          write(6,'(A)') 'Or at least less than than 0.0001.'
          write(6,'(A)') 'Did you employ enough digits for the rotation matrix?'
          write(6,'(A)') 'The program will stop.'
          return
        end if
112     continue
      end do
    end do
  end if ! axisoption
end if ! compute_CF

! preparing the info for computation of molar magnetization
if (compute_magnetization) then
  ! calculate the total number of directions for the average procedure
  nDirTot = 0
  if (zeeman_energy) nDirTot = nDirZee
  if (compute_Mdir_vector) nDirTot = nDirTot+nDir
  nDirTot = nDirTot+get_nP(nsymm,ngrid)
end if

!------ CHECK the data from INPUT ------------------------------
!if (iprint > 10) then

if (dbg) write(6,'(A,  F9.5)') 'ZJPR :         = ',zJ
if (dbg) write(6,'(A,  I3  )') 'PRLV :         = ',iprint

!if (.not. compute_g_tensors) then
!  !generate an array of 10 low-lying groups of states
!
!  ndim(:) = 0
!  nmult_try = 10
!  j = 0
!  ndim(i) = 1
!  do i=1,nss
!    etmp = eso(i)
!    do k=i+1,nss
!      if (abs(eso(k)-etmp) < 0.01_wp) ndim(i) = ndim(i)+1
!    end do
!  end do
!end if

if (compute_g_tensors) then
  if (NMULT > 0) then
    write(6,'(A,I3)') 'MLTP :         = ',NMULT
    if (NMULT <= 20) then
      write(6,'(A,20I3)') '               = ',(NDIM(i),i=1,NMULT)
    else
      write(6,'(A,20I3)') '               = ',(NDIM(i),i=1,20)
      do j=21,NMULT,20
        jEnd = min(NMULT,J+19)
        write(6,'(A,20I3)') '                 ',(NDIM(i),i=j,jEnd)
      end do
    end if
  else
    write(6,'(A)') 'MLTP :         =  No pseudospin Hamiltonians will be computed. Is MLTP defined?'
  end if
end if

if (compute_CF .and. (nDIMcf <= nss) .and. (lDIMcf <= nstate)) then
  if (nlanth < 15) then
    write(6,'(3A)') 'The Crystal-Field acting on the ground atomic multiplet of Ln = ',clanth(nlanth),' is computed.'
  else if ((nlanth >= 15) .and. (nlanth < 29)) then
    write(6,'(3A)') 'The Crystal-Field acting on the ground atomic multiplet of Ac = ',clanth(nlanth),' is computed.'
  else if (nlanth >= 29) then
    write(6,'(3A)') 'The Crystal-Field acting on the ground atomic |L,ML> multiplet of TM = ',clanth(nlanth),' is computed.'
  end if

  write(6,'(A,A )') 'CHIT :         = ',' molar magnetic susceptibility is computed'
  if (TINPUT) write(6,'(A)') 'TEXP :         = the experimental temperature interval is read from the file "chitexp.input"'
  write(6,'(A, I3)') 'TINT :      nT = ',nT
  write(6,'(A,F7.3)') '          Tmin = ',Tmin
  write(6,'(A,F7.3)') '          Tmax = ',Tmax

end if
!--------------------------------------------------------------------
if (compute_magnetization) then

  write(6,'(A,A )') 'MAGN :         = ',' molar magnetization is computed'
  write(6,'(A, I3)') 'NDIRTOT        = ',nDirTot
  write(6,'(A, I3)') 'TMAG :         = ',nTempMagn
  write(6,'(6x,A,20F7.3)') 'TempMagn = ',(TempMagn(i),i=1,nTempMagn)
  write(6,'(A, I3)') 'HINT :      nH = ',nH
  write(6,'(A,F7.3)') '          Hmin = ',Hmin
  write(6,'(A,F7.3)') '          Hmax = ',Hmax
  write(6,'(A, I3)') 'MAVE :   nDir = ',get_nP(nsymm,ngrid)

  if (HINPUT) write(6,'(A)') 'HEXP :         = the experimental field interval is read from the file "mexp.input"'
  if (encut_definition == 1) then
    write(6,'(A, I3)') 'NCUT :         = ',ncut
  else if (encut_definition == 2) then
    write(6,'(A,I4,a,i4)') 'ECUT :         = ',nk,', ',mg
  else if (encut_definition == 3) then
    write(6,'(A,F7.3)') 'ERAT :         = ',encut_rate
  end if

  if (compute_Mdir_vector) then
    write(6,'(A,20I3)') 'MVEC :         = ',nDir
    if (nDir > 0) then
      do i=1,nDir
        write(6,'(A,I2,A,3F11.6)') '   Dir :',i,' : ',dirX(i),dirY(i),dirZ(i)
      end do
    end if
  end if
  if (zeeman_energy) then
    if (nDirZee == 1) then
      write(6,'(2A,I2,1x,A)') 'ZEEM :         = ',' Zeeman splitting for the following direction of the '
      write(6,'(18x,A)') 'applied magnetic field is given in the "zeeman_energy_xxx.txt" file in $WorkDir/'
    else if (nDirZee > 1) then
      write(6,'(2A,I2,1x,A)') 'ZEEM :         = ',' Zeeman splitting for the following',nDirZee,' directions of the '
      write(6,'(18x,A)') 'applied magnetic field are given in the "zeeman_energy_xxx.txt" files in $WorkDir/.'
    else
      write(6,'(A)') 'Error in input processing. nDirZee<0!'
      call Quit_OnUserError()
    end if
    do i=1,nDirZee
      write(6,'(17x,3F11.6)') (dir_weight(i,l),l=1,3)
    end do
  end if
end if ! magnetization
!--------------------------------------------------------------------

if (compute_torque)  write(6,'(A,A )') 'TORQ :         = ',' torque magnetization is computed'

if (doplot) write(6,'(A,A )') 'PLOT :         = ',' GNUPLOT scripts and corresponding XT, M and UBAR plots will be generated'

if (Do_structure_abc) then
  write(6,'(2A)') 'ABCC :         = ','the main magnetic axes for the computed pseudospins are written also in the '
  write(6,'( A)') 'crystallographic "abc" axes'
  write(6,'(10x,A,F9.4)') 'a       = ',cryst(1)
  write(6,'(10x,A,F9.4)') 'b       = ',cryst(2)
  write(6,'(10x,A,F9.4)') 'c       = ',cryst(3)
  write(6,'(10x,A,F9.4)') 'alpha   = ',cryst(4)
  write(6,'(10x,A,F9.4)') 'beta    = ',cryst(5)
  write(6,'(10x,A,F9.4)') 'gamma   = ',cryst(6)
  write(6,'(10x,a,3F9.4)') 'coords: = ',(coord(i),i=1,3)
end if

if (compute_g_tensors) then
  if (compute_barrier) then
    Nblock = 0
    do i=1,nmult
      Nblock = Nblock+ndim(i)
    end do
    write(6,'(A,i4)') 'nBlock = ',nBlock
  end if
end if

go to 190
!------ errors ------------------------------
write(6,*) ' The following input line was not understood:'
write(6,'(A)') LINE
go to 999

997 continue
write(6,*) ' READIN: Error reading standard input.'
write(6,*) ' SINGLE_ANISO input near line nr.',LINENR+1
go to 999

998 continue
write(6,*) ' READIN: Unexpected End of input file.'

999 continue
call XFLUSH(6)
call ABEnd()

590 continue
write(6,*) 'READIN: the TINT command is incompatible with TEXP'
call ABEnd()

591 continue
write(6,*) 'READIN: the HINT command is incompatible with HEXP'
call ABEnd()

595 continue
write(6,*) 'READIN: NCUT, ERAT and ENCU are mutually exclusive.'
call ABEnd()

190 continue

return

end subroutine readin_single
