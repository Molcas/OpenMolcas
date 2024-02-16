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

subroutine Readin_poly(nneq,neq,neqv,exch,nCenter,nT,nH,nTempMagn,nDir,nDirZee,nMult,nPair,nexch,nDim,i_pair,lant,multLn,iPrint, &
                       keopt,encut_definition,nK,mG,iopt,nP,AngPoints,ncut,LUZee,MxRank1,MxRank2,imaxrank,TempMagn,R_LG,R_ROT,Jex, &
                       JAex,JAex9,JDMex,JITOexR,JITOexI,tpar,upar,cryst,coord,Xfield,gtens_input,D_fact,EoverD_fact,riso, &
                       MagnCoords,thrs,tmin,tmax,hmin,hmax,Texp,chit_exp,Hexp,Mexp,encut_rate,zJ,dirX,dirY,dirZ,dir_weight,Title, &
                       itype,ifHDF,compute_g_tensors,compute_magnetization,TINPUT,HINPUT,Do_structure_abc,DoPlot, &
                       compute_Mdir_vector,zeeman_energy,m_paranoid,m_accurate,smagn,compute_susceptibility,decompose_exchange,KE, &
                       fitCHI,fitM,compute_torque,compute_barrier,Dipol,check_title,AnisoLines1,AnisoLines3,AnisoLines9, &
                       DM_exchange,JITO_exchange)
! THIS ROUTINE READS THE standard input.

use Constants, only: Zero, One, Three
use Definitions, only: wp, u5, u6

implicit none
#include "mgrid.fh"
#include "warnings.h"
! definition of the cluster:
integer :: nneq, neqv, neq(nneq), nCenter
logical :: ifHDF
! definition of the local metal sites
real(kind=8) :: R_LG(nneq,neqv,3,3)
real(kind=8) :: R_ROT(nneq,neqv,3,3)
real(kind=8) :: gtens_input(3,nneq)
real(kind=8) :: D_fact(nneq)
real(kind=8) :: EoverD_fact(nneq)
real(kind=8), intent(out) :: riso(nneq,3,3)
character(len=1) :: itype(nneq)
! definition of the exchange:
integer :: exch ! total number of exchange states
integer :: nPair ! number of metal pairs (number of interactions)
integer :: nexch(nneq) ! exchange basis, nmax= MAX(nexch(:))
integer :: i_pair(nPair,2) ! index of the metal site in a given interacting pair
logical :: AnisoLines1, AnisoLines3, AnisoLines9
logical :: Dipol, DM_exchange
real(kind=8) :: Jex(nPair) ! Lines exchange    ( 1 parameter / interacting pair)
real(kind=8) :: JAex(nPair,3) ! Anisotropic Lines ( 3 parameter / interacting pair)
real(kind=8) :: JAex9(nPair,3,3) ! Anisotropic Lines full ( 9 parameters / interacting pair)
real(kind=8) :: JDMex(nPair,3)
logical :: JITO_exchange ! options used in connection with ITO exchange:
integer, intent(in) :: MxRank1, MxRank2
integer, intent(out) :: imaxrank(npair,2)
real(kind=8), intent(out) :: JITOexR(nPair,MxRank1,-MxRank1:MxRank1,MxRank2,-MxRank2:MxRank2)
real(kind=8), intent(out) :: JITOexI(nPair,MxRank1,-MxRank1:MxRank1,MxRank2,-MxRank2:MxRank2)
! options used in connection with KE
integer :: lant, KEOPT, multLn
logical :: KE
real(kind=8) :: tpar, upar
! options used in connection with Dipol-Dipol interaction
real(kind=8) :: MagnCoords(nneq,3)
! definition of g and D tensors
integer :: nMult
integer :: nDim(nMult)
logical :: compute_g_tensors
! definition of data for susceptibility
integer :: nT
logical :: tinput, compute_susceptibility
real(kind=8) :: tmin, tmax
real(kind=8) :: chit_exp(nT), Texp(nT)
! options related to XT_MoverH
real(kind=8) :: Xfield
! definition of data for magnetization:
integer :: nH
integer :: nTempMagn
integer :: iopt
real(kind=8) :: TempMagn(nTempMagn)
real(kind=8) :: Hexp(nH), Mexp(nH,nTempMagn)
real(kind=8) :: thrs
real(kind=8) :: hmin, hmax
logical :: hinput
logical :: compute_magnetization
logical :: compute_Mdir_vector
logical :: zeeman_energy
logical :: m_paranoid
logical :: m_accurate
logical :: smagn
! options used to set up nM and EM
integer :: encut_definition
integer :: nK, mG ! encut_definition=1;
integer :: ncut   ! encut_definition=2;
real(kind=8) :: encut_rate ! encut_definition=3;
! decompose exchange
logical :: decompose_exchange
! magnetization torque
integer :: nP
integer :: AngPoints
logical :: compute_torque
! Zeeman energy and M vector
integer :: nDir, nDirZee
integer :: LUZee(nDirZee)
real(kind=8) :: dirX(nDir), dirY(nDir), dirZ(nDir)
real(kind=8) :: dir_weight(nDirZee,3)
! definition of mean field parameter
real(kind=8) :: zJ
! definition of the crystal axes:
logical :: Do_structure_abc
real(kind=8) :: cryst(6) ! a, b, c, alpha, beta, gamma
! Cartesian coordinates of the main metal site, or center
real(kind=8) :: coord(3)
! definitions for blocking barrier
logical :: compute_barrier
! options for automatic fitting of parameters:
logical :: fitCHI !-- not used so far
logical :: fitM !-- not used so far
! definition of print level
integer :: iPrint
logical :: check_title
character(len=180) :: Title
logical :: DoPlot
!--------- LOCAL VARIABLES --------------------
real(kind=8) :: check_dir_weight(3*nDirZee), tmp, sum
integer :: ll, i, j, l, m, n, icount_b_sites, lp, ic, jc
integer :: linenr, inneq, irank1, irank2, iproj1, iproj2
integer :: duplicate_check(nPair), nind(exch,2)
integer :: nst, ASUM, jrank1, jrank2, jproj1, jproj2
integer :: i1, i2, lb1, lb2
! index of the metal site in a given interacting pair
!integer :: nind(nPair,2)
!integer :: ind_exch(nneq)
!real(kind=wp) :: magncoords(2*maxanisofiles,3)
real(kind=8) :: finddetr, detR
real(kind=8) :: tmpR(3,3)
logical :: nosym
external :: finddetr
! variables connected to computation of g and d tensors
logical :: ab_initio_all
logical :: tcheck, hcheck, encut_check
logical :: check_symm_presence
logical :: checktmag
real(kind=8) :: t2, t1
character(len=2) :: lanth
!character(len=14) :: namefile_energy(nDirZee)
character(len=21) :: namefile_energy
character(len=288) :: Line, ctmp, string
integer :: IsFreeUnit
external :: IsFreeUnit
logical :: DBG

DBG = .false.

check_title = .false.
icount_B_sites = 0
i_pair = 0
nosym = .true.
ENCUT_check = .false.
TINPUT = .false.
HINPUT = .false.
TCHECK = .false.
HCHECK = .false.

do i=1,nneq
  do j=1,Neq(i)
    R_rot(i,j,1,1) = One
    R_rot(i,j,2,2) = One
    R_rot(i,j,3,3) = One
    R_lg(i,j,1,1) = One
    R_lg(i,j,2,2) = One
    R_lg(i,j,3,3) = One
  end do
end do
check_symm_presence = .true.
if (neqv > 1) check_symm_presence = .false.

if (DBG) write(u6,'(A,  i6)') 'RDIN:      nneq=',nneq
if (DBG) write(u6,'(A,  i6)') 'RDIN:      neqv=',neqv
if (DBG) write(u6,'(A,  i6)') 'RDIN:     nPair=',nPair
!If (DBG) write(u6,'(A,99I6)') 'RDIN:  i_pair(:,1)=',(i_pair(i,1),i=1,nPair)
!if (DBG) write(u6,'(A,99I6)') 'RDIN:  i_pair(:,2)=',(i_pair(i,2),i=1,nPair)
if (DBG) write(u6,'(A,99i6)') 'RDIN:     neq()=',(neq(i),i=1,nneq)
if (DBG) write(u6,'(A,99i6)') 'RDIN:   nexch()=',(nexch(i),i=1,nneq)

!=========== End of default settings====================================
rewind(u5)
50 read(u5,'(A72)',end=998) LINE
if (DBG) write(u6,'(A)') LINE
call NORMAL(LINE)
if (LINE(1:5) /= '&POLY') Go To 50
LINENR = 0

100 call xFlush(u6)
read(u5,'(A72)',end=998) LINE
LINENR = LINENR+1
call NORMAL(LINE)
if (LINE(1:1) == '*') Go To 100
if (LINE == ' ') Go To 100
if (LINE(1:3) == 'END') Go To 210 !End

! ------------ TITL ---------------------------------------------------*
if (LINE(1:4) == 'TITL') then

  read(u5,*,err=997) ctmp

  if (DBG) write(u6,'(A)') ctmp
  check_title = .true.
  Title = trim(ctmp)
  LINENR = LINENR+1
  Go To 100
end if

! ------------ OLDA ---------------------------------------------------*
if (LINE(1:4) == 'OLDA') then
  LINENR = LINENR+1
  Go To 100
end if

!---  process MLTP command --------------------------------------------*
if (LINE(1:4) == 'MLTP') then

  read(u5,*,err=997) NMULT

  if (DBG) write(u6,'(A,i4)') 'NMULT =',NMULT

  read(u5,*,Err=997) (NDIM(i),i=1,NMULT)

  do i=1,NMULT
    if (NDIM(i) < 0) then
      ctmp = ''
      write(ctmp,'(A,I2,A)') 'MLTP: the dimension of the multiplet ',i,' is negative!!!'
      call WarningMessage(2,ctmp)
      call quit(_RC_INPUT_ERROR_)
    else if (NDIM(i) == 0) then
      ctmp = ''
      write(ctmp,'(A,I2,A)') 'MLTP: the dimension of the multiplet ',i,' is zero!!!'
      call WarningMessage(2,ctmp)
      call quit(_RC_INPUT_ERROR_)
    end if
  end do

  if (DBG) write(u6,'(A,100i4)') 'NDIM: ',(NDIM(i),i=1,NMULT)
  compute_g_tensors = .true.
  LINENR = LINENR+2
  Go To 100
end if

!---  process TINT command --------------------------------------------*
if (LINE(1:4) == 'TINT') then
  compute_susceptibility = .true.
  if (.not. TINPUT) then
    TCHECK = .true.

    read(u5,*,err=997) t1,t2,nT

    if ((t1 < 0) .or. (t2 < 0)) then
      call WarningMessage(2,'TINT: negative temperature requested! ')
      call quit(_RC_INPUT_ERROR_)
    end if
    if ((t1-t2) > Zero) then
      Tmin = t2
      Tmax = t1
    else if ((t1-t2) < Zero) then
      Tmin = t1
      Tmax = t2
    else ! t1==t2
      call WarningMessage(2,'TINT: temperature interval == 0! ')
      call quit(_RC_INPUT_ERROR_)
    end if

    if (DBG) write(u6,'(A,2ES15.7,i6)') 'Tmin, Tmax, nT: ',Tmin,Tmax,nT
  else
    goto 590
  end if
  LINENR = LINENR+1
  Go To 100
end if

!---  process HINT command --------------------------------------------*
if (LINE(1:4) == 'HINT') then
  if (.not. HINPUT) then
    HCHECK = .true.
    compute_magnetization = .true.

    read(u5,*,err=997) t1,t2,nH

    if ((t1 < 0) .or. (t2 < 0)) then
      call WarningMessage(2,'HINT: negative field requested! ')
      call quit(_RC_INPUT_ERROR_)
    end if
    if ((t1-t2) > Zero) then
      Hmin = t2
      Hmax = t1
    else if ((t1-t2) < Zero) then
      Hmin = t1
      Hmax = t2
    else ! t1==t2
      call WarningMessage(2,'HINT: temperature interval == 0! ')
      call quit(_RC_INPUT_ERROR_)
    end if

    if (DBG) write(u6,'(A,2ES15.7,i6)') 'Hmin, Hmax, nH: ',Hmin,Hmax,nH
  else
    Go To 591
  end if
  LINENR = LINENR+1
  Go To 100
end if

!---  process THRS command --------------------------------------------*
if (LINE(1:4) == 'THRS') then

  read(u5,*,err=997) THRS

  if (thrs < Zero) then
    call WarningMessage(2,'THRS: negative threshold for  average M!!! ')
    write(u6,'(A)') 'Set to default thrs=1.0e-10'
    thrs = 1.0e-10_wp
  end if

  if (DBG) write(u6,'(A,ES15.7)') 'THRS: ',THRS
  LINENR = LINENR+1
  Go To 100
end if

!---  process XFIE command --------------------------------------------*
if (LINE(1:4) == 'XFIE') then
  compute_susceptibility = .true.

  read(u5,*,err=997) tmp

  if (tmp < Zero) then
    call WarningMessage(2,'XFIE: negative value of the applied field !')
    write(u6,'(A)') 'Set to positive !'
    Xfield = abs(tmp)
  else if (tmp == Zero) then
    call WarningMessage(2,'XFIE: zero value of the applied field !')
    write(u6,'(A)') 'Field-applied XT will not be computed!'
  else
    Xfield = tmp
  end if
  if (DBG) write(u6,'(A,ES15.7)') 'XFIE: ',xField

  LINENR = LINENR+1
  Go To 100
end if

!---  process PLOT command --------------------------------------------*
if (LINE(1:4) == 'PLOT') then
  DoPlot = .true.

  if (DBG) write(u6,'(A,L2)') 'PLOT: ',DoPlot

  LINENR = LINENR+1
  Go To 100
end if

!---  process IOPT command --------------------------------------------*
if (LINE(1:4) == 'IOPT') then

  read(u5,*,err=997) I  !option for computing MSUM and XTSUM

  if ((i < 0) .or. (i > 3)) then
    call WarningMessage(2,'IOPT: value out of range!!!')
    write(u6,'(A)') 'Set to default !'
    iopt = 1
  else
    iopt = i
  end if

  if (DBG) write(u6,'(A,I6)') 'IOPT: ',IOPT
  LINENR = LINENR+1
  Go To 100
end if

!---  process SMAG command --------------------------------------------*
if (LINE(1:4) == 'SMAG') then
  smagn = .true.
  compute_magnetization = .true.
  if (DBG) write(u6,'(A,L2)') 'SMAG: ',smagn
  LINENR = LINENR+1
  Go To 100
end if

!---  process MACC command --------------------------------------------*
if (LINE(1:4) == 'MACC') then
  compute_magnetization = .true.  ! request for computation of M(H)
  m_accurate = .true.             ! request for computation of M(H)
  if (DBG) write(u6,'(A,L2)') 'MACC: ',m_accurate
  LINENR = LINENR+1
  Go To 100
end if

!---  process FITX command --------------------------------------------*
if (LINE(1:4) == 'FITX') then
  fitCHI = .true.
  LINENR = LINENR+1
  Go To 100
end if

!---  process FITM command --------------------------------------------*
if (LINE(1:4) == 'FITM') then
  fitM = .true.
  LINENR = LINENR+1
  Go To 100
end if

!---  process MPAR command --------------------------------------------*
if (LINE(1:4) == 'MPAR') then
  m_paranoid = .true.             ! request for computation of M(H)
  compute_magnetization = .true.
  if (DBG) write(u6,'(A,L2)') 'MPAR: ',m_paranoid
  LINENR = LINENR+1
  Go To 100
end if

!---  process TORQ command --------------------------------------------*
if (LINE(1:4) == 'TORQ') then

  compute_torque = .true.         ! request for computation of M(H)
  read(u5,*,err=997) i        ! number of angular points

  if (i <= 0) then
    call WarningMessage(2,'TORQ: nP value out of range!!!')
    write(u6,'(A)') 'Set to default !'
    nP = 45
  else
    nP = i
  end if

  if (DBG) write(u6,'(A)') 'TORQ: nP=',nP
  AngPoints = nP+1
  LINENR = LINENR+1
  Go To 100
end if

!---  process MAVE command --------------------------------------------*
if (LINE(1:4) == 'MAVE') then
  compute_magnetization = .true.  ! request for computation of M(H)

  read(u5,*,err=997) i,j  !nsymm, ngrid

  if ((i <= 0) .or. (i >= 4)) then
    call WarningMessage(2,'MAVE: nSYMM value out of range!!!')
    write(u6,'(A)') 'Set to default !'
    nsymm = 1
  else
    nsymm = i
  end if

  if ((j <= 0) .or. (j >= 33)) then
    call WarningMessage(2,'MAVE: nGRID value out of range!!!')
    write(u6,'(A)') 'Set to default !'
    ngrid = 15
  else
    ngrid = i
  end if

  if (DBG) write(u6,'(2(A,i6))') ' nsymm:  ',nsymm,' ngrid:  ',ngrid
  LINENR = LINENR+1
  Go To 100
end if

!---  process NCUT command --------------------------------------------*
if (LINE(1:4) == 'NCUT') then
  if (ENCUT_check) then
    Go To 595
  else
    ENCUT_check = .true.
    ! request for computation of M(H)
    compute_magnetization = .true.
    encut_definition = 1

    read(u5,*,err=997) i ! NCUT
    !E_cut = exchange_energy(Ncut)

    if ((i <= 0) .or. (i > exch)) then
      call WarningMessage(2,'NCUT: value out of range!!!')
      write(u6,'(A)') 'Set to full exchange basis.'
      nCUT = exch
    else
      nCUT = i
    end if
    if (DBG) write(u6,'(A,2i6)') 'ncut:  ',nCut

    LINENR = LINENR+1
    Go To 100
  end if
end if

!---  process ENCU command --------------------------------------------*
if (LINE(1:4) == 'ENCU') then
  if (ENCUT_check) then
    Go To 595
  else
    ENCUT_check = .true.
    ! request for computation of M(H)
    compute_magnetization = .true.
    encut_definition = 2

    read(u5,*,err=997) NK,MG
    !E_cut = NK*K_Boltz+MG*mu_Bohr
    if (DBG) write(u6,'(A,2i6)') 'encu:  nK, mG=',NK,MG

    LINENR = LINENR+1
    Go To 100
  end if
end if

!---  process ERAT command --------------------------------------------*
if (LINE(1:4) == 'ERAT') then
  if (ENCUT_check) then
    Go To 595
  else
    ENCUT_check = .true.
    ! request for computation of M(H)
    compute_magnetization = .true.
    encut_definition = 3

    read(u5,*,err=997) encut_rate
    !Ncut = int(nexch*encut_rate)
    !E_cut = E(Ncut)
    if (DBG) write(u6,'(A,i6)') 'encut_rate=',encut_rate

    LINENR = LINENR+1
    Go To 100
  end if
end if

!---  process ZJPR command --------------------------------------------*
if (LINE(1:4) == 'ZJPR') then
  read(u5,*,err=997) ZJ
  if (DBG) write(u6,'(A,ES18.10)') 'zJ    =',zJ
  LINENR = LINENR+1
  Go To 100
end if

!---  process PRLV command --------------------------------------------*
if (LINE(1:4) == 'PRLV') then
  read(u5,*,err=997) iPrint
  if (DBG) write(u6,'(A,i6)') 'iPrint=',iPrint
  LINENR = LINENR+1
  Go To 100
end if

!---  process COOR command --------------------------------------------*
if (LINE(1:4) == 'COOR') then
  Dipol = .true.
  if (DBG) write(u6,'(A)') 'isite   MagnCoords:'
  do i=1,nneq
    read(u5,*,err=997) (MagnCoords(i,l),l=1,3)
    if (DBG) write(u6,'(i3,5x,3ES18.10)') i,(MagnCoords(i,l),l=1,3)
  end do
  LINENR = LINENR+nneq
  Go To 100
end if

!---  process MVEC command --------------------------------------------*
if (LINE(1:4) == 'MVEC') then
  compute_magnetization = .true.  ! request for computation of M(H)
  compute_Mdir_vector = .true.

  read(u5,*,err=997) nDir
  if (DBG) write(u6,'(A,i3)') 'nDir = ',nDir

  do i=1,nDir
    read(u5,*,err=997) DirX(i),DirY(i),DirZ(i)
    if (DBG) write(u6,'(i3,5x,3ES18.10)') i,DirX(i),DirY(i),DirZ(i)
  end do
  ! some processing:
  do i=1,nDir
    sum = DirX(i)*DirX(i)+DirY(i)*DirY(i)+DirZ(i)*DirZ(i)
    if (sum == Zero) then
      write(u6,'(a,i3,a)') 'error: MVEC  vector ',i,'has the modulus = 0.0.'
      write(u6,'(a     )') 'the program will stop now.'
      call quit(_RC_INPUT_ERROR_)
    end if
    if (abs(sum-One) > 0.5e-13_wp) then
      write(u6,'(a,i3,a)') 'the vector ',i,'was re-normalized.'
      tmp = dirX(i)/sqrt(sum)
      dirX(i) = tmp
      tmp = dirY(i)/sqrt(sum)
      dirY(i) = tmp
      tmp = dirZ(i)/sqrt(sum)
      dirZ(i) = tmp
    end if
  end do

  LINENR = LINENR+NDIR+1
  Go To 100
end if

!---  process TEXP command --------------------------------------------*
if (LINE(1:4) == 'TEXP') then
  compute_susceptibility = .true.
  if (.not. TCHECK) then
    TINPUT = .true.

    read(u5,*,err=997) NT
    if (DBG) write(u6,'(A,i3)') 'nT = ',nT

    do i=1,NT

      read(u5,*,err=997) texp(i),chit_exp(i)

      ! check and clean negative values:
      if (texp(i) < Zero) texp(i) = abs(texp(i))
      if (chit_exp(i) < Zero) chit_exp(i) = abs(chit_exp(i))
    end do
    tmin = texp(1)
    tmax = texp(nT)
  else
    Go To 590
  end if
  LINENR = LINENR+NT+1
  Go To 100
end if

!---  process HEXP command --------------------------------------------*
if (LINE(1:4) == 'HEXP') then
  compute_magnetization = .true.
  if (checkTMAG) then
    write(u6,'(A)') 'The data provided in TMAG will be ignored.'
  end if

  if (.not. HCHECK) then
    HINPUT = .true.

    read(u5,*) nTempMagn,(TempMagn(i),i=1,nTempMagn)
    read(u5,*) nH
    if (DBG) write(u6,*) 'HEXP: nTempMagn =',nTempMagn
    if (DBG) write(u6,*) 'HEXP: TempMagn()=',(TempMagn(i),i=1,nTempMagn)
    if (DBG) write(u6,*) 'HEXP: nH        =',nH

    if (nH < 0) nH = abs(nH)
    if (nH == 0) call Quit_OnUserError()

    do i=1,nH
      read(u5,*,err=997) Hexp(i),(Mexp(i,j),j=1,nTempMagn)
      ! check and clean negative values:
      if (hexp(i) < Zero) hexp(i) = abs(hexp(i))
      do j=1,nTempMagn
        if (Mexp(i,j) < Zero) Mexp(i,j) = abs(Mexp(i,j))
      end do
    end do
    hmin = hexp(1)
    hmax = hexp(nH)
  else
    Go To 591
  end if
  LINENR = LINENR+NH+2
  Go To 100
end if
!---  process TMAG command --------------------------------------------*
if (LINE(1:4) == 'TMAG') then
  if (.not. HINPUT) then
    compute_magnetization = .true.
    checkTMAG = .true.

    read(u5,*,err=997) nTempMagn,(TempMagn(i),i=1,nTempMagn)
    if (DBG) write(u6,*) 'TMAG: nTempMagn =',nTempMagn
    if (DBG) write(u6,*) 'TMAG: TempMagn()=',(TempMagn(i),i=1,nTempMagn)

    ! check and clean for zero / negative values:
    do i=1,nTempMagn
      if (TempMagn(i) < Zero) then
        call WarningMessage(2,'TMAG: negative temperature requested! ')
        write(u6,'(A)') 'Set to positive.'
        TempMagn(i) = abs(TempMagn(i))
      else if (TempMagn(i) == Zero) then
        call WarningMessage(2,'TMAG: zero temperature requested! ')
        write(u6,'(A)') 'Set to 0.0001 K.'
        TempMagn(i) = 0.0001_wp
      end if
    end do
  else
    write(u6,'(A)') 'TMAG data is taken from HEXP.'
  end if
  LINENR = LINENR+1
  Go To 100
end if
!---  process NNEQ command --------------------------------------------*
! this is the most important keyword for Poly_Aniso
if (LINE(1:4) == 'NNEQ') then
  ! number of non-equivalent centers; type of all centers
  read(u5,*,err=997) NNEQ,ab_initio_all,ifHDF
  if (DBG) write(u6,'(A,i4,A,L2,A,L2)') 'NNEQ=',NNEQ,' ab_initio_all=',ab_initio_all,' ifHDF=',ifHDF
  ! number of equivalent centers of type "i"
  read(u5,*,err=997) (NEQ(i),i=1,Nneq)
  if (DBG) write(u6,'(A,100I4)') 'NEQ(I)=',(NEQ(i),i=1,nneq)
  ! number of RASSI wf for exchange
  read(u5,*,err=997) (Nexch(i),i=1,Nneq)
  if (DBG) write(u6,'(A,100I4)') 'NExch(I)=',(NExch(i),i=1,nneq)

  do i=1,nneq
    if (Nexch(i) < 0) then
      write(u6,'(A,i1,a)') 'The number of functions taken in the exchange interaction from center of type ',i,' is negative!'
      write(u6,'(A     )') 'The program has to stop, since the input is not reasonable.'
      call quit(_RC_INPUT_ERROR_)
    else if (Nexch(i) == 0) then
      write(u6,'(A,i1,a)') 'The number of functions taken in the exchange interaction from center of type ',i,' is 0 (zero).'
      write(u6,'(A)') 'The program has to stop, since the input is not reasonable.'
      call quit(_RC_INPUT_ERROR_)
    end if
  end do
  do i=1,nneq
    if (neq(i) > 1) then
      nosym = .false.
    end if
  end do
  !write(u6,'(A,i5)') 'exch = ',exch
  if (exch == 1) then
    write(u6,'(3/)')
    write(u6,'(100A)') ('#',i=1,100)
    write(u6,'(3/)')
    write(u6,'(A)') 'The size of the exchange matrix is 1. Is this really what you intended to compute?'
    write(u6,'(A)') 'The program will continue...'
    write(u6,'(3/)')
    write(u6,'(100A)') ('#',i=1,100)
    write(u6,'(3/)')
  end if

  !do i=1,nCenter
  !  ind_exch(i) = i
  !end do

  ! If the EXCH is above the limit => exit with an error
  if (exch > 15000) then
    write(u6,'(A)') 'The number of exchange states is very large'
    write(u6,'(A)') 'EXCH=',exch
    write(u6,'(A)') 'The calculation will continue, but might take a LOT of time'
    write(u6,'(A)') 'We recomend to switch OFF the computation of powder magnetization'
  end if

  if (.not. ab_initio_all) then

    !type of the center:   A -- the information is read from aniso_ion.input
    !                      B -- the center is isotropic with g factor read from the input
    !                      C -- the center is anisotropic with
    read(u5,*,err=997) (itype(i),i=1,Nneq)
    if (DBG) write(u6,'(A,100A3)') 'itype: ',(itype(i),i=1,nneq)
    if (DBG) call xFlush(u6)
    icount_B_sites = 0
    do i=1,nneq
      if ((itype(i) == 'B') .or. (itype(i) == 'C')) then
        icount_B_sites = icount_B_sites+1
        read(u5,*,err=997) (gtens_input(l,i),l=1,3),D_fact(i),EoverD_fact(i)
        if (DBG) write(u6,'(A,i4,A,3ES20.10, 2(A,ES20.10) )') 'gtens_input(',i,')=',(gtens_input(l,i),l=1,3),' D = ',D_fact(i), &
                                                              ' E/D =',EoverD_fact(i)
        if ((itype(i) == 'C') .and. ((gtens_input(1,i) /= gtens_input(2,i)) .or. (gtens_input(1,i) /= gtens_input(3,i)) .or. &
            (gtens_input(2,i) /= gtens_input(3,i)))) then
          do ic=1,3
            read(u5,*,err=997) (riso(i,jc,ic),jc=1,3)
          end do

        else
          call unitmat(riso(i,:,:),3)
        end if
      end if

      if (abs(D_fact(i)) > Zero) then
        if (abs(EoverD_fact(i)/D_fact(i)) > One/Three) then
          write(string,'(A,i3,A)') 'NNEQ: |E/D| for center',i,' > 1/3'
          call WarningMessage(2,trim(string))
          write(u6,'(A,i3,A)') '|E/D| for center ',i,' is set to 1/3'
        end if
      end if
    end do
  else
    do i=1,nneq
      itype(i) = 'A'
    end do
  end if !ab_initio_all
  LINENR = LINENR+3+icount_B_sites
  if (DBG) call xFlush(u6)
  Go To 100
end if

!---  process LIN9 command --------------------------------------------*
if (LINE(1:4) == 'LIN9') then
  AnisoLines9 = .true.

  read(u5,*,err=997) npair
  if (DBG) write(u6,'(A,i6)') 'LIN9:  nPair=',nPair

  do i=1,npair
    ! the convention for 9 exchange interactions is
    ! Jxx, Jxy, Jxz, Jyx, Jyy, Jyz, Jzx, Jzy, Jzz
    read(u5,*,err=997) i_pair(i,1),i_pair(i,2),(JAex9(i,1,j),j=1,3),(JAex9(i,2,j),j=1,3),(JAex9(i,3,j),j=1,3)
    if (DBG) write(u6,'(A,2I3,9F14.8)') 'LIN9: ',i_pair(i,1),i_pair(i,2),(JAex9(i,1,j),j=1,3),(JAex9(i,2,j),j=1,3), &
                                        (JAex9(i,3,j),j=1,3)
  end do
  LINENR = LINENR+npair+1
  Go To 100
end if

!---  process LIN3 command --------------------------------------------*
if ((LINE(1:4) == 'ALIN') .or. (LINE(1:4) == 'LIN3')) then
  AnisoLines3 = .true.

  read(u5,*,err=997) npair
  if (DBG) write(u6,'(A,i6)') 'nPair=',nPair

  do i=1,npair
    ! the convention for 3 exchange interactions is
    ! Jxx, Jyy, Jzz
    read(u5,*,err=997) i_pair(i,1),i_pair(i,2),(JAex(i,j),j=1,3)
    if (DBG) write(u6,'(A,2i3,3F14.8)') 'ALIN/LIN3: ',i_pair(i,1),i_pair(i,2),(JAex(i,j),j=1,3)
  end do
  LINENR = LINENR+npair+1
  Go To 100
end if

!---  process JITO command --------------------------------------------*
if (LINE(1:4) == 'ITOJ') then
  JITO_exchange = .true.

  read(u5,*,err=997) npair
  if (DBG) write(u6,'(A,i6)') 'nPair=',nPair
  JITOexR(:,:,:,:,:) = Zero
  JITOexI(:,:,:,:,:) = Zero

  do i=1,npair
    ! the convention for exchange interactions is

    read(u5,*,err=997) i_pair(i,1),i_pair(i,2),imaxrank(i,1),imaxrank(i,2)
    do irank1=1,imaxrank(i,1),2
      do iproj1=-irank1,irank1
        do irank2=1,imaxrank(i,2),2
          do iproj2=-irank2,irank2
            read(u5,*,err=997) jrank1,jproj1,jrank2,jproj2,JITOexR(i,jrank1,jproj1,jrank2,jproj2),JITOexI(i,jrank1,jproj1,jrank2, &
                               jproj2)
          end do
        end do
      end do
    end do

    if (DBG) then
      write(u6,'(A,I3)') 'ITO Exchange parameters for pair:',i
      do jrank1=1,imaxrank(i,1),2
        do jproj1=-jrank1,jrank1
          do jrank2=1,imaxrank(i,2),2
            do jproj2=-jrank2,jrank2
              write(u6,'(4I3,2x,2ES21.14)') jrank1,jproj1,jrank2,jproj2,JITOexR(i,jrank1,jproj1,jrank2,jproj2),JITOexI(i,jrank1, &
                                            jproj1,jrank2,jproj2)
            end do
          end do
        end do
      end do
    end if
  end do
  LINENR = LINENR+npair+1
  Go To 100
end if

!---  process PAIR = LIN1 command -------------------------------------*
if ((LINE(1:4) == 'PAIR') .or. (LINE(1:4) == 'LIN1')) then
  AnisoLines1 = .true.

  read(u5,*,err=997) npair
  if (DBG) write(u6,'(A,i6)') 'nPair=',nPair

  do i=1,npair

    read(u5,*,err=997) i_pair(i,1),i_pair(i,2),Jex(i)
    if (DBG) write(u6,'(i4,2x,2I4,2x,ES18.10)') i,i_pair(i,1),i_pair(i,2),Jex(i)

  end do
  LINENR = LINENR+npair+1
  Go To 100
end if

!---  process DMEX command --------------------------------------------*
if (LINE(1:4) == 'DMEX') then
  DM_exchange = .true.
  if (DBG) write(u6,'(A,L2)') 'DMEX::  DM_exchange=',DM_exchange

  read(u5,*,err=997) npair
  if (DBG) write(u6,'(A,i6)') 'nPair=',nPair

  do i=1,npair
    ! the convention for 3 exchange interactions is
    ! JDMex(x, y, z)
    read(u5,*,err=997) i_pair(i,1),i_pair(i,2),(JDMex(i,j),j=1,3)
  end do
  LINENR = LINENR+npair+1
  Go To 100
  LINENR = LINENR+1
  Go To 100
end if

!---  process ZEEM command --------------------------------------------*
if (LINE(1:4) == 'ZEEM') then
  zeeman_energy = .true.
  compute_magnetization = .true.
  LUZEE = 0

  read(u5,*,err=997) nDirZee

  do i=1,nDirZee
    ! open the zeeman_energy_xxx.txt file where Zeeman eigenstates will
    ! be further written in mangetization() subroutine
    write(namefile_energy,'(5A)') 'zeeman_energy_',char(48+mod(int((i)/100),10)),char(48+mod(int((i)/10),10)), &
                                  char(48+mod(int(i),10)),'.txt'
    if (DBG) write(u6,'(2A)') 'namefile_energy: ',namefile_energy
    LUZee(i) = IsFreeUnit(30+i)
    call molcas_open(LUZee(i),namefile_energy)

    read(u5,*,err=997) (dir_weight(i,l),l=1,3)

    check_dir_weight(i) = sqrt(dir_weight(i,1)**2+dir_weight(i,2)**2+dir_weight(i,3)**2)

    if (abs(check_dir_weight(i)-One) > 0.005_wp) then
      write(u6,'(A)') 'The directions for the magnetic field for the computation of the Zeeman splitting are wrong.'
      write(u6,'(A)') '( px^2 + py^2 + pz^2 ) must give 1.!'
      write(u6,'(A,I3,2x,A,F9.5)') 'In the present case for direction Nr.',i,' the dir_weight = px^2 + py^2 + pz^2 = ', &
                                   check_dir_weight(i)**2
      LINENR = LINENR+2+i
      Go To 997
    end if

  end do
  LINENR = LINENR+nDirZee+1
  Go To 100
end if

!---  process MAGN command --------------------------------------------*
!if (LINE(1:4) == 'MAGN') then
!  compute_magnetization = .true.
!  go to 100
!end if

!---  process SYMM command --------------------------------------------*
if (LINE(1:4) == 'SYMM') then
  nosym = .false.
  check_symm_presence = .true.
  R_lg(:,:,:,:) = Zero
  R_rot(:,:,:,:) = Zero
  if (DBG) write(u6,'(A,i6)') 'SYMM - at init'
  ll = 0
  do i=1,nneq
    if (DBG) write(u6,'(A,i6)') 'SYMM:  i=',i

    read(u5,*,err=997) inneq
    if (DBG) write(u6,'(A,i6)') 'inneq=',inneq

    do j=1,Neq(i)
      if (DBG) write(u6,'(A,i6)') 'SYMM:  j=',j
      do m=1,3
        ll = ll+1
        read(u5,*,err=997) (R_lg(i,j,m,n),n=1,3)
        if (DBG) write(u6,'(3ES20.12)') (R_lg(i,j,m,n),n=1,3)
      end do
    end do

    do j=1,neq(i)
      do m=1,3
        do n=1,3
          tmpR(m,n) = R_lg(i,j,m,n)
        end do
      end do

      detR = FindDetR(tmpR(1:3,1:3),3)
      if (DBG) write(u6,'(A,3ES20.12)') 'SYMM:  detR=',detR

      if (abs(abs(detR)-One) > 0.001_wp) then
        write(u6,'(A)') 'The rotation matrices must be UNITARY.'
        write(u6,'(A)') 'and of RIGHT hand system'
        write(u6,'(A,F11.6)') 'DET = ',detR
        LINENR = LINENR+ll+i+1
        Go To 997
      end if
      if (detR < Zero) then
        do m=1,3
          do n=1,3
            R_rot(i,j,m,n) = -R_lg(i,j,m,n)
          end do !n
        end do !m
      else if (detR > Zero) then
        do m=1,3
          do n=1,3
            R_rot(i,j,m,n) = R_lg(i,j,m,n)
          end do !n
        end do !m
      end if
    end do !neq(i)
  end do !nneq
  LINENR = LINENR+nneq+ll
  Go To 100
end if

!---  process EXCH command --------------------------------------------*
if (LINE(1:4) == 'EXCH') then
  decompose_exchange = .true.
  Go To 100
end if

!---  process OLDA command --------------------------------------------*
!if (LINE(1:4) == 'OLDA') then
!  old_aniso_format = .true.
!  go to 100
!end if

!---  process EXCH command --------------------------------------------*
!if (LINE(1:4) == 'END') then
!  go to 100
!end if

!---  process UBAR command --------------------------------------------*
if (LINE(1:4) == 'UBAR') then
  compute_barrier = .true.
  !read(u5,*,err=997) icase
  !if (icase == 1) then
  !! icase =1  --> magnetic field is applied along the main magnetic
  !!               axis of each doublet (multiplet)
  !  continue
  !else if (icase == 2) then
  !! icase =2  --> magnetic field is applied along the main magnetic
  !!               axis of the specified doublet (number) number (NDim)
  !  read(u5,*,err=997) NmagMult
  !else if (icase == 3) then
  !! icase =3  --> a new coordination system is defined by the user,
  !!               and the magnetic field is applied along gZ (third axis)
  !  do i=1,3
  !    read(u5,*,err=997) (uBar_Rot(i,j),j=1,3)
  !  end do
  !  LINENR = LINENR+3
  !else
  !  write(u6,'(A)') 'Is the UBAR keyword used correctly?'
  !  write(u6,'(A)') 'The ICASE parameter is not understood.'
  !end if
  LINENR = LINENR+1
  Go To 100
end if

!----------------------------------------------------------------------*
if (LINE(1:4) == 'ABCC') then
  Do_structure_abc = .true.
  read(u5,*,err=997) (cryst(i),i=1,6)
  coord(:) = Zero
  !read(u5,*,err=997) (coord(i),i=1,3)
  LINENR = LINENR+2
  Go To 100
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

!---  process LONG command --------------------------------------------*
if (LINE(1:4) == 'LONG') then
  KE = .true.

  read(u5,*,err=997) lanth,tpar,upar,KEOPT

  if ((lanth == 'GD') .or. (lanth == 'gd') .or. (lanth == 'gD') .or. (lanth == 'Gd')) then
    lant = 1
    multLn = 8
  else if ((lanth == 'TB') .or. (lanth == 'tb') .or. (lanth == 'tB') .or. (lanth == 'Tb')) then
    lant = 2
    multLn = 13
  else if ((lanth == 'DY') .or. (lanth == 'dy') .or. (lanth == 'dY') .or. (lanth == 'Dy')) then
    lant = 3
    multLn = 16
  else if ((lanth == 'HO') .or. (lanth == 'ho') .or. (lanth == 'hO') .or. (lanth == 'Ho')) then
    lant = 4
    multLn = 17
  else if ((lanth == 'ER') .or. (lanth == 'er') .or. (lanth == 'eR') .or. (lanth == 'Er')) then
    lant = 5
    multLn = 16
  else if ((lanth == 'TM') .or. (lanth == 'tm') .or. (lanth == 'tM') .or. (lanth == 'Tm')) then
    lant = 6
    multLn = 13
  else if ((lanth == 'YB') .or. (lanth == 'yb') .or. (lanth == 'yB') .or. (lanth == 'Yb')) then
    lant = 7
    multLn = 8
  else
    write(u6,'( A)') 'Error in getting the type of the lanthanide!'
    write(u6,'(2A)') 'The program has this info: lanth =',lanth
    write(u6,'(26x,A,i2)') 'multLn =',multLn
  end if
  LINENR = LINENR+2
  Go To 100
end if

! end of reading input keywords
210 continue
!-----------------------------------------------------------------------

!=====  Perform some PROCESSING of the data ============================

if (.not. nosym) then
  if (.not. check_symm_presence) then
    write(u6,'(A)') 'The SYMM keyword is mandatory if the system '
    write(u6,'(A)') 'contains more than one center of the same type!'
    call quit(_RC_INPUT_ERROR_)
  end if !check_symm_presence
end if !nosym

! preparing the info for computation of molar magnetization
! calculate the total number of directions for the average procedure

!--------  definition g and D tensors ----------------------------------
if (compute_g_tensors) then
  nst = 0
  do i=1,nMult
    nst = nst+ndim(i)
    if (nst > exch) then
      write(u6,'(A)') 'You have requested the computation of properties of more states'
      write(u6,'(A)') 'than the total number of exchange states.'
      write(u6,'(A)') 'The program does not know how to compute'
      write(u6,'(A)') 'properties of inexistent states'
      write(u6,'(A)') 'New settings:'
      write(u6,'(A,i6)') 'NMULT = ',i-1

      nmult = i-1
      Go To 14
    end if
  end do
14 continue
end if
if (DBG) write(u6,*) 'READIN_POLY:  after proc g and D'
if (DBG) call xFlush(u6)

!--------  definition of exchange --------------------------------------
if (npair > 0) then
  ASUM = 0
  do lp=1,nPair
    ASUM = ASUM+i_pair(lp,1)+i_pair(lp,2)
  end do
  if (Dipol .and. (ASUM == 0)) then
    lp = 0
    do i=1,nCenter-1
      do j=i+1,nCenter
        if (i >= j) goto 17
        lp = lp+1
        i_pair(lp,1) = i
        i_pair(lp,2) = j
        if (DBG) write(u6,'(A,i3,A,2I3)') 'lp=',lp,' i_pair(lp,1:2)=',i_pair(lp,1),i_pair(lp,2)
17      continue
      end do
    end do
  end if

  Duplicate_check(1:nPair) = 0
  do i=1,npair
    Duplicate_check(i) = 1000*i_pair(i,1)+i_pair(i,2)
    if (DBG) write(u6,*) 'Duplicate_check: ',Duplicate_check(i)

    if (i_pair(i,1) == i_pair(i,2)) then
      write(u6,'(A,i2,a)') 'The center ',i_pair(i,1),' interacts with itself.'
      write(u6,'(A)') 'This is not possible. The program has to stop.'
      call quit(_RC_INPUT_ERROR_)

    else if (i_pair(i,1) > i_pair(i,2)) then
      write(u6,'(A)') 'The convention of this program enforces the label of the first magnetic center of a given '
      write(u6,'(A)') 'interaction to be smaller than the label of the second magnetic center.'
      write(u6,'(A)') 'In order to avoid confusions, please respect this convention.'
      write(u6,'(A)') 'Avoid duplicate interactions!'
      write(u6,'(A)') 'The program will stop now.'
      call quit(_RC_INPUT_ERROR_)
    end if
  end do
  if (DBG) call xFlush(u6)

  ! check on the duplicate
  do i=1,npair
    do j=i+1,npair
      if (j == i) Go To 197
      if (Duplicate_check(i) == Duplicate_check(j)) then
        write(u6,'(A)') 'Some interactions are declared twice in the input. Please declare all interactions only once!'
        write(u6,'(A)') 'The program has to stop.'
        call quit(_RC_INPUT_ERROR_)
      end if
197   continue
    end do
  end do

  ! check on the indices of the exchange couplings:
  do i=1,nPair
    if ((i_pair(i,1) > nCenter) .or. (i_pair(i,2) > nCenter)) then
      write(u6,'(A)') 'The numbering of the magnetic centers within NPAIR keyword is wrong.'
      write(u6,'(A,i2,a)') 'For the interaction Nr. = ',i,' you numerate individual centers with numbers larger than the total '// &
                           'number of centers in this molecule.'
      write(u6,'(A)') 'NPAIR keyword is wrong.'
      write(u6,'(A,I4)') 'i_pair(i,1) =',i_pair(i,1)
      write(u6,'(A,I4)') 'i_pair(i,2) =',i_pair(i,2)
      write(u6,'(A)') 'The program has to stop.'
      call quit(_RC_INPUT_ERROR_)
    end if
  end do

  ! check on the size of MxRank1 and MxRank2 wrt nexch(1) and nexch(2)
  if (JITO_exchange) then
    l = 0
    do i=1,nneq
      do j=1,neq(i)
        l = l+1
        nind(l,1) = i
        nind(l,2) = j
      end do
    end do

    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      if (nExch(i1) < (imaxrank(lp,1)+1)) then
        write(u6,'(100A)') ('#',i=1,100)
        write(u6,'(A)') 'interacting pair = ',lp
        write(u6,'(A)') 'type of site 1 = ',i1
        write(u6,'(A,i2,A,i2)') 'nExch(',i1,') = ',nExch(i1)
        write(u6,'(A,i2,A,i2)') 'Rank2(',lp,',1) = ',imaxrank(lp,1)
        write(u6,'(A)') 'nExch < Rank+1 !!!'
        write(u6,'(A)') 'Rank of ITO operators for site 1 is larger than the number of defined exchange states for this site'
        write(u6,'(A)') 'The program will use only parameters which bring non-zero contribution to exchange.'
        write(u6,'(100A)') ('#',i=1,100)
      end if
      if (nExch(i2) < (imaxrank(lp,2)+1)) then
        write(u6,'(100A)') ('#',i=1,100)
        write(u6,'(A)') 'interacting pair = ',lp
        write(u6,'(A)') 'type of site 2 = ',i2
        write(u6,'(A,i2,A,i2)') 'nExch(',i2,') = ',nExch(i2)
        write(u6,'(A,i2,A,i2)') 'Rank2(',lp,',2) = ',imaxrank(lp,2)
        write(u6,'(A)') 'nExch < Rank+1 !!!'
        write(u6,'(A)') 'Rank of ITO operators for site 2 is larger than the number of defined exchange states for this site'
        write(u6,'(A)') 'The program will use only parameters which bring non-zero contribution to exchange.'
        write(u6,'(100A)') ('#',i=1,100)
      end if
    end do
  end if ! JITO_exchange
end if ! nPair

if (DBG) write(u6,*) 'READIN_POLY:  before 200 '
if (DBG) call xFlush(u6)
!--------  definition of exchange --------------------------------------

Go To 200
!-----------------------------------------------------------------------
write(u6,*) ' The following input line was not understood:'
write(u6,'(A)') LINE
Go To 999

997 continue
write(u6,*) ' READIN_POLY: Error reading "poly_aniso.input" '
write(u6,*) ' near line nr.',LINENR+1
Go To 999
998 continue
write(u6,*) ' READIN_POLY: Unexpected End of input file.'
999 continue
call quit(_RC_INPUT_ERROR_)
590 continue
write(u6,*) 'READIN_POLY: THE TINT command is incompatible with TEXP'
call quit(_RC_INPUT_ERROR_)
591 continue
write(u6,*) 'READIN_POLY: THE HINT command is incompatible with HEXP'
call quit(_RC_INPUT_ERROR_)
595 continue
write(u6,*) 'READIN_POLY: THE NCUT, ENCU, and ERAT are mutually exclusive. You cannot use more than one keyword at the same time.'
call quit(_RC_INPUT_ERROR_)
! ===============   NORMAL EndING  =====================================
200 continue

if (IPRINT > 2) write(u6,'(5X,A)') 'NO ERROR WAS LOCATED WHILE READING INPUT'

return

end subroutine Readin_poly
