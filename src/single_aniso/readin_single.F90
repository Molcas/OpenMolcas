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

subroutine readin_single(iprint,nmult,ndim,ndimcf,ldimcf,nlanth,axisoption,poly_file,Ifrestart,input_to_read,nk,mg,zmagn, &
                         Do_structure_abc,cryst,coord,encut_definition,compute_g_tensors,compute_CF,nDirTot,nss,nstate, &
                         compute_magnetization,compute_torque,smagn,tinput,hinput,compute_Mdir_vector,zeeman_energy,LUZee,doplot, &
                         encut_rate,ncut,nTempMagn,TempMagn,m_paranoid,compute_barrier,nBlock,AngPoints,input_file_name,nT,nH, &
                         texp,chit_exp,zJ,hexp,magn_exp,hmin,hmax,nDir,nDirZee,dirX,dirY,dirZ,dir_weight,xfield,tmin,tmax,thrs, &
                         H_torq,T_torq,nsymm,ngrid)
! THIS ROUTINE READS THE FILE "SINGLE_ANISO.INPUT".

use Lebedev_quadrature, only: order_table
use Constants, only: Zero, One, Two, Five, Six, Ten
use Definitions, only: wp, iwp, u5, u6

implicit none
integer(kind=iwp), intent(inout) :: nmult, ntempmagn, nt, nh, nDir, nDirZee
integer(kind=iwp), intent(out) :: iprint, ndim(nMult), ndimcf, ldimcf, nlanth, axisoption, input_to_read, nk, mg, &
                                  encut_definition, ndirtot, LUZee(nDirZee), ncut, nBlock, AngPoints, nsymm, ngrid
logical(kind=iwp), intent(out) :: poly_file, Ifrestart, Do_structure_abc, compute_g_tensors, compute_cf, compute_magnetization, &
                                  compute_torque, smagn, tinput, hinput, compute_Mdir_vector, zeeman_energy, doplot, m_paranoid, &
                                  compute_barrier
real(kind=wp), intent(out) :: zmagn(3,3), cryst(6), coord(3), encut_rate, tempmagn(nTempMagn), texp(nT), chit_exp(nT), zj, &
                              hexp(nH), magn_exp(nH,ntempmagn), hmin, hmax, dirX(nDir), dirY(nDir), dirZ(nDir), &
                              dir_weight(nDirZee,3), Xfield, tmin, tmax, thrs, H_torq, T_torq
integer(kind=iwp), intent(in) :: nss, nstate
character(len=180), intent(inout) :: input_file_name
integer(kind=iwp) :: I, i_OxStat, istatus, j, jEnd, l, LINENR
real(kind=wp) :: check_dir_weight, column_check(3,3), det_zmagn, row_check(3,3), rsum, t1, t2, tmp, zr(3,3)
logical(kind=iwp) :: checktmag, encut_check, hcheck, tcheck
character(len=280) :: LINE
character(len=180) :: err_msg, tmpline
character(len=21) :: namefile_energy
character(len=2) :: cME, uME
integer(kind=iwp), parameter :: ngrid_map(32) = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59, &
                                                 62,65]
character(len=*), parameter :: clanth(37) = ['CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU', & ! lanthanides
                                             'TH','PA','U ','NP','PU','AM','CM','BK','CF','ES','FM','MD','NO','LR', & ! actinides
                                             'SC','TI','V ','CR','MN','FE','CO','NI','CU'] ! transition metals
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: FindDetR

!========== Initializations of arrays ==================================
DirX(:) = Zero
DirY(:) = Zero
DirZ(:) = Zero
dir_weight(:,:) = Zero
TempMagn(:) = Zero
!============ Initializations of constants =============================

nsymm = 1
ngrid = 15
IPRINT = 2
ndim(:) = 0
Ifrestart = .false.
nDirTot = 0
LuZee(:) = 0
nBlock = 0
texp(:) = Zero
chit_exp(:) = Zero
hexp(:) = Zero
magn_exp(:,:) = Zero
H_torq = 0.1_wp ! in tesla
T_torq = Two    ! in K

thrs = 1.0e-10_wp
ncut = 1
TMIN = Zero
TMAX = 300.0_wp
XFIELD = Zero
!NT = 301
HMIN = Zero
HMAX = Ten
!NH = 21
NK = 200
MG = 200
!NDIR = 0
!nDirZee = 0
!nTempMagn = 1
if (nTempMagn > 0) TempMagn(1) = Two
T1 = Five
T2 = Six
ZJ = Zero
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
lDIMcf = 0
cME = '  '

cryst(:) = Zero
coord(:) = Zero
axisoption = 1
input_to_read = 0
encut_definition = 2
encut_rate = 1
zmagn(:,:) = Zero
AngPoints = 46

!=========== End of default settings====================================
rewind(u5)
do
  read(u5,'(A280)',iostat=istatus) LINE
  if (istatus < 0) call Error(1)
# ifdef _DEBUGPRINT_
  write(u6,'(A)') trim(LINE)
# endif
  call NORMAL(LINE)
  if (LINE(1:7) == '&SINGLE') exit
end do
LINENR = 0
do
  read(u5,'(A280)',iostat=istatus) LINE
  if (istatus < 0) call Error(1)
# ifdef _DEBUGPRINT_
  write(u6,'(A)') trim(LINE)
# endif
  LINENR = LINENR+1
  call NORMAL(LINE)
  if ((LINE(1:1) == '*') .or. (LINE == ' ')) cycle

  select case (LINE(1:4))
    case ('END ','    ')
      exit
    !------------------------------------------
    !case ('TYPE')
    !  read(u5,*,iostat=istatus) ICALC
    !  if (istatus < 0) call Error(2)
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
    !    write(u6,'(A)') 'ICALC: the maximum value is 7. However, the calculation will continue by computing the magnetism.'
    !    compute_g_tensors = .true.
    !    compute_chiT = .true.
    !    compute_magnetization = .true.
    !  end if
    !  LINENR = LINENR+1
    !------------------------------------------
    case ('MLTP')
      read(u5,*,iostat=istatus) NMULT
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'MLTP:  NMULT=',NMULT
#     endif
      compute_g_tensors = .true.
      read(u5,*,iostat=istatus) (NDIM(i),i=1,NMULT)
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'MLTP: NDIM()=',(NDIM(i),i=1,NMULT)
#     endif
      LINENR = LINENR+2
    !------------------------------------------
    case ('REST')
      Ifrestart = .true.
      read(u5,*,iostat=istatus) input_to_read
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'REST: input_to_read=',input_to_read
#     endif
      if ((input_to_read == 2) .or. (input_to_read == 3) .or. (input_to_read == 4)) then
        backspace(u5)
        read(u5,*) input_to_read,tmpline
        input_file_name = trim(tmpline)
      end if
      if (input_to_read == 1) then
        write(u6,*) 'RESTART: -- The SINGLE_ANISO will take all ab initio information from the binary "$Project.aniso" file.'
      else if (input_to_read == 2) then
        write(u6,*) 'RESTART: -- The SINGLE_ANISO will take all ab initio information from the ASCII '//trim(input_file_name)// &
                    ' file.'
      else if (input_to_read == 3) then
        write(u6,*) 'RESTART: -- The SINGLE_ANISO will take all ab initio information from the RASSI-HDF5 binary file.'
      else if (input_to_read == 4) then
        write(u6,*) 'RESTART: -- The SINGLE_ANISO will take all ab initio information from the ASCII '//trim(input_file_name)// &
                    ' file -- molcas-8.0 format.'
      else
        call WarningMessage(2,'SINGLE_ANISO:: RESTART  option is not known.')
        call Quit_OnUserError()
      end if
    !-------------------------------------------
    case ('DATA')
      Ifrestart = .true.
      read(u5,*) tmpline
      input_file_name = trim(tmpline)
      input_to_read = 6
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: DATA, input_file_name='
      write(u6,*) input_file_name
#     endif
      LINENR = LINENR+1
    !-------------------------------------------
    case ('TINT')
      if (.not. TINPUT) then
        TCHECK = .true.

        read(u5,*,iostat=istatus) t1,t2,nT
        if (istatus < 0) call Error(2)

        if ((t1 < 0) .or. (t2 < 0)) then
          call WarningMessage(2,'TINT: negative temperature requested! ')
          call Quit_OnUserError()
        end if
        if ((t1-t2) > Zero) then
          Tmin = t2
          Tmax = t1
        else if ((t1-t2) < Zero) then
          Tmin = t1
          Tmax = t2
        else ! t1==t2
          call WarningMessage(2,'TINT: temperature interval == 0! ')
          call Quit_OnUserError()
        end if

#       ifdef _DEBUGPRINT_
        write(u6,*) 'TINT: Tmin, Tmax, nT=',Tmin,Tmax,nT
#       endif
      else
        write(u6,*) 'READIN_SINGLE: the TINT command is incompatible with TEXP'
        call ABEnd()
      end if
      LINENR = LINENR+1
    !-------------------------------------------
    case ('XFIE')
      read(u5,*,iostat=istatus) Xfield
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'XFIE: Xfield=',Xfield
#     endif
      LINENR = LINENR+1
    !-------------------------------------------
    case ('HINT')
      if (.not. HINPUT) then
        HCHECK = .true.
        compute_magnetization = .true.

        read(u5,*,iostat=istatus) t1,t2,nH
        if (istatus < 0) call Error(2)

        if ((t1 < 0) .or. (t2 < 0)) then
          call WarningMessage(2,'HINT: negative field requested! ')
          call Quit_OnUserError()
        end if

        if ((t1-t2) > Zero) then
          Hmin = t2
          Hmax = t1
        else if ((t1-t2) < Zero) then
          Hmin = t1
          Hmax = t2
        else ! t1 == t2
          call WarningMessage(2,'HINT: temperature interval == 0! ')
          call Quit_OnUserError()
        end if

#       ifdef _DEBUGPRINT_
        write(u6,*) 'HINT: Hmin, Hmax, nH=',Hmin,Hmax,nH
#       endif
      else
        write(u6,*) 'READIN_SINGLE: the HINT command is incompatible with HEXP'
        call ABEnd()
      end if
      LINENR = LINENR+1
    !-------------------------------------------
    case ('NCUT')
      if (ENCUT_check) then
        write(u6,*) 'READIN_SINGLE: NCUT, ERAT and ENCU are mutually exclusive.'
        call ABEnd()
      end if
      ENCUT_check = .true.
      compute_magnetization = .true. ! request for computation of M(H)
      encut_definition = 1

      read(u5,*,iostat=istatus) NCUT  ! E_cut=ESO(Ncut)
      if (istatus < 0) call Error(2)

      if (NCUT < 0) then
        call WarningMessage(2,'NCUT: negative NCUT requested! ')
        call Quit_OnUserError()
      else if (NCUT == 0) then
        call WarningMessage(2,'NCUT: zero NCUT requested! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'NCUT: NCUT=',NCUT
#     endif
      LINENR = LINENR+1
    !-------------------------------------------
    case ('ENCU')
      if (ENCUT_check) then
        write(u6,*) 'READIN_SINGLE: NCUT, ERAT and ENCU are mutually exclusive.'
        call ABEnd()
      end if
      ENCUT_check = .true.
      compute_magnetization = .true.
      encut_definition = 2

      read(u5,*,iostat=istatus) NK,MG
      if (istatus < 0) call Error(2)

      if ((NK <= 0) .or. (MG <= 0)) then
        call WarningMessage(2,'ENCU: zero or negative NK,MG requested! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'ENCU: NK, MG=',NK,MG
#     endif
      LINENR = LINENR+1
    !-------------------------------------------
    case ('ERAT')
      if (ENCUT_check) then
        write(u6,*) 'READIN_SINGLE: NCUT, ERAT and ENCU are mutually exclusive.'
        call ABEnd()
      end if
      ENCUT_check = .true.
      compute_magnetization = .true.
      encut_definition = 3
      ! Ncut = INT(nss*encut_rate)
      ! E_cut = E(Ncut)

      read(u5,*,iostat=istatus) encut_rate
      if (istatus < 0) call Error(2)

      if (encut_rate <= Zero) then
        call WarningMessage(2,'ERAT: zero or negative encut rate requested! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'ERAT: encut_rate=',encut_rate
#     endif
      LINENR = LINENR+1
    !-------------------------------------------
    case ('MVEC')
      compute_magnetization = .true.   ! request for computation of M(H)
      compute_Mdir_vector = .true.
      read(u5,*,iostat=istatus) nDir
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'MVEC: nDir=',nDir
#     endif
      do i=1,nDir
        read(u5,*,iostat=istatus) DirX(i),DirY(i),DirZ(i)
        if (istatus < 0) call Error(2)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'MVEC: DirX,DirY,DirZ=',DirX(i),DirY(i),DirZ(i)
#       endif
      end do
      ! some processing:
      do i=1,nDir
        rsum = DirX(i)*DirX(i)+DirY(i)*DirY(i)+DirZ(i)*DirZ(i)
        if (rsum == Zero) then
          write(err_msg,'(a,i3,a)') 'error: MVEC  vector ',i,'has the modulus = 0.0.'
          call WarningMessage(2,err_msg)
          call Quit_OnUserError()
        end if
        if (rsum /= One) then
          write(u6,'(a,i3,a)') 'the vector',i,'was re-normalized.'
          tmp = dirX(i)/sqrt(rsum)
          dirX(i) = tmp
          tmp = dirY(i)/sqrt(rsum)
          dirY(i) = tmp
          tmp = dirZ(i)/sqrt(rsum)
          dirZ(i) = tmp
        end if
      end do

      LINENR = LINENR+NDIR+1
    !-------------------------------------------
    case ('MAVE')
      compute_magnetization = .true.

      read(u5,*,iostat=istatus) nsymm,ngrid
      if (istatus < 0) call Error(2)

#     ifdef _DEBUGPRINT_
      write(u6,*) 'MAVE: nsymm, ngrid=',nsymm,ngrid
#     endif
      if ((nsymm < 1) .or. (nsymm > 3)) then
        write(u6,'(A)') '"nsymm" must take Integer values 1, 2 or 3.'
        write(u6,'(A,i5)') '"nsymm" = ',nsymm
        call Quit_OnUserError()
      end if
      if ((ngrid < 1) .or. (ngrid > 32)) then
        write(u6,'(A)') '"ngrid" must take Integer values 1, 2, ... 32.'
        write(u6,'(A,i5)') '"ngrid" = ',ngrid
        call Quit_OnUserError()
      end if
      LINENR = LINENR+1
    !-------------------------------------------
    !case ('TLIN')
    !  read(u5,*,iostat=istatus) T1,T2
    !  if (istatus < 0) call Error(2)
    !  LINENR = LINENR+1
    !-------------------------------------------
    case ('SMAG')
      smagn = .true.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'SMAG: =',smagn
#     endif
      LINENR = LINENR+1
    !-------------------------------------------
    case ('PLOT')
      doplot = .true.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'PLOT: =',doplot
#     endif
      LINENR = LINENR+1
    !-------------------------------------------
    case ('TEXP')
      if (.not. TCHECK) then
        TINPUT = .true.
        read(u5,*,iostat=istatus) NT
        if (istatus < 0) call Error(2)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'TEXP: nT=',nT
#       endif
        do i=1,NT
          read(5,*,iostat=istatus) texp(i),chit_exp(i)
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'TEXP: texp(i), chit_exp(i)=',texp(i),chit_exp(i)
#         endif
          ! check and clean negative values:
          if (texp(i) < Zero) texp(i) = abs(texp(i))
          if (chit_exp(i) < Zero) chit_exp(i) = abs(chit_exp(i))
        end do
        tmin = texp(1)
        tmax = texp(nT)
      else
        write(u6,*) 'READIN_SINGLE: the TINT command is incompatible with TEXP'
        call ABEnd()
      end if
      LINENR = LINENR+NT+1
    !-------------------------------------------
    case ('HEXP')
      compute_magnetization = .true.
      if (checkTMAG) write(u6,'(A)') 'The data provided in TMAG will be ignored.'
      if (.not. HCHECK) then
        HINPUT = .true.
        read(u5,*) nTempMagn,(TempMagn(i),i=1,nTempMagn)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'HEXP: nTempMagn =',nTempMagn
        write(u6,*) 'HEXP: TempMagn()=',(TempMagn(i),i=1,nTempMagn)
#       endif
        read(u5,*) nH
#       ifdef _DEBUGPRINT_
        write(u6,*) 'HEXP: nH =',nH
#       endif
        if (nH < 0) nH = abs(nH)
        if (nH == 0) call Quit_OnUserError()
        do i=1,nH
          read(u5,*,iostat=istatus) Hexp(i),(magn_exp(i,j),j=1,nTempMagn)
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'HEXP: Hexp(i),  magn_exp(i,j)=',Hexp(i),(magn_exp(i,j),j=1,nTempMagn)
#         endif
          ! check and clean negative values:
          if (hexp(i) < Zero) hexp(i) = abs(hexp(i))
          do j=1,nTempMagn
            if (magn_exp(i,j) < Zero) magn_exp(i,j) = abs(magn_exp(i,j))
          end do
        end do
        hmin = hexp(1)
        hmax = hexp(nH)
      else
        write(u6,*) 'READIN_SINGLE: the HINT command is incompatible with HEXP'
        call ABEnd()
      end if
      LINENR = LINENR+NH+2
    !-------------------------------------------
    case ('ZJPR')
      read(u5,*,iostat=istatus) ZJ
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'ZJPR: zJ =',zJ
#     endif
      LINENR = LINENR+1
    !-------------------------------------------
    case ('TORQ')
      compute_torque = .true.
      read(u5,*,iostat=istatus) AngPoints,H_torq,T_torq
      if (istatus < 0) call Error(2)
      LINENR = LINENR+1
    !-------------------------------------------
    case ('TMAG')
      if (.not. HINPUT) then
        compute_magnetization = .true.
        checkTMAG = .true.

        read(u5,*,iostat=istatus) nTempMagn,(TempMagn(i),i=1,nTempMagn)
        if (istatus < 0) call Error(2)

        do i=1,nTempMagn
          if (TempMagn(i) <= Zero) then
            call WarningMessage(2,'TMAG: zero or negative temperature requested! ')
            if (TempMagn(i) < Zero) TempMagn(i) = abs(TempMagn(i))
            if (TempMagn(i) == Zero) TempMagn(i) = 0.0001_wp
          end if
        end do

#       ifdef _DEBUGPRINT_
        write(u6,*) 'TMAG: nTempMagn =',nTempMagn
        write(u6,*) 'TMAG: TempMagn()=',(TempMagn(i),i=1,nTempMagn)
#       endif
        ! check and clean negative values:
      else
        write(u6,'(A)') 'TMAG data is taken from HEXP.'
      end if
      LINENR = LINENR+1
    !-------------------------------------------
    case ('PRLV')
      read(u5,*,iostat=istatus) IPRINT
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'PRLV: IPRINT =',iPrint
#     endif
      LINENR = LINENR+1
    !-------------------------------------------
    case ('POLY')
#     ifdef _DEBUGPRINT_
      write(u6,*) 'POLY:'
#     endif
      POLY_FILE = .true.
    !-------------------------------------------
    case ('CRYS')
      compute_CF = .true.
      read(u5,*,iostat=istatus) cME
      uME = cME
      call UpCase(uME)
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'CRYS: cME =',cME
#     endif

      select case (uME)
        ! LANTHANIDES
        case (clanth(1)) ! Ce
          nlanth = 1
          nDIMcf = 6  ! f1; multiplet J=L-S=3-1/2=5/2  =>  J = 2F_5/2
          lDIMCF = 7  ! (L=3)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(2)) ! Pr
          nlanth = 2
          nDIMcf = 9  ! f2; multiplet J=L-S=5-1=4  => J = 3H_4
          lDIMCF = 11 ! (L=5)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(3)) ! Nd
          nlanth = 3
          nDIMcf = 10 ! f3; multiplet J=L-S=6-3/2=9/2  => J = 4I_9/2
          lDIMCF = 13 ! (L=6)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(4)) ! Pm
          nlanth = 4
          nDIMcf = 9  ! f4; multiplet J=L-S=6-2=4  => J = 5I_4
          lDIMCF = 13 ! (L=6)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(5)) ! Sm
          nlanth = 5
          nDIMcf = 6  ! f5; multiplet J=L-S=5-5/2=5/2  => J = 6H_5/2
          lDIMCF = 11 ! (L=5)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(6)) ! Eu
          nlanth = 6
          nDIMcf = 1  ! f6; multiplet J=L-S=3-3=0  => J = 3F_0
          lDIMCF = 7  ! (L=3)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(7)) !Gd
          nlanth = 7
          nDIMcf = 8  ! f7; multiplet J=L+S=0+7/2=0  => J = 8S_7/2
          lDIMCF = 1  ! (L=0)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(8)) !Tb
          nlanth = 8
          nDIMcf = 13 ! f8; multiplet J=L+S=3+3=0  => J = 7F_6
          lDIMCF = 7  ! (L=3)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(9)) ! Dy
          nlanth = 9
          nDIMcf = 16 ! f9; multiplet J=L+S=5+5/2=15/2  => J = 6H_15/2
          lDIMCF = 11 ! (L=5)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(10)) ! Ho
          nlanth = 10
          nDIMcf = 17 ! f10; multiplet J=L+S=6+2=8  => J = 5I_8
          lDIMCF = 13 ! (L=6)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(11)) ! Er
          nlanth = 11
          nDIMcf = 16 ! f11; multiplet J=L+S=6+3/2=15/2  => J = 4I_15/2
          lDIMCF = 13 ! (L=6)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(12)) ! Tm
          nlanth = 12
          nDIMcf = 13 ! f12; multiplet J=L+S=5+1=6  => J = 3H_6
          lDIMCF = 11 ! (L=5)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(13)) ! Yb
          nlanth = 13
          nDIMcf = 8  ! f13; multiplet J=L+S=3+1/2=7/2  => J = 2F_7/2
          lDIMCF = 7  ! (L=3)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(14)) !Lu
          nlanth = 14
          nDIMcf = 1  ! f14; multiplet J=L+S=0+0=0  => J = 1S_0
          lDIMCF = 1  ! (L=0)
        !- - - - - - - - - - - - - - - - - - - -
        ! ACTINIDES
        case (clanth(15)) ! Th
          nlanth = 15
          nDIMcf = 6  ! f1; multiplet J=L-S=3-1/2=5/2  =>  J = 2F_5/2
          lDIMCF = 7  ! (L=3)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(16)) ! Pa
          nlanth = 16
          nDIMcf = 9  ! f2; multiplet J=L-S=5-1=4  => J = 3H_4
          lDIMCF = 11 ! (L=5)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(17)) ! U
          nlanth = 17
          nDIMcf = 10 ! f3; multiplet J=L-S=6-3/2=9/2  => J = 4I_9/2
          lDIMCF = 13 ! (L=6)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(18)) ! Np
          nlanth = 18
          nDIMcf = 9  ! f4; multiplet J=L-S=6-2=4  => J = 5I_4
          lDIMCF = 13 ! (L=6)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(19)) ! Pu
          nlanth = 19
          nDIMcf = 6  ! f5; multiplet J=L-S=5-5/2=5/2  => J = 6H_5/2
          lDIMCF = 11 ! (L=5)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(20)) ! An
          nlanth = 20
          nDIMcf = 1  ! f6; multiplet J=L-S=3-3=0  => J = 3F_0
          lDIMCF = 7  ! (L=3)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(21)) ! Cm
          nlanth = 21
          nDIMcf = 8  ! f7; multiplet J=L+S=0+7/2=0  => J = 8S_7/2
          lDIMCF = 1  ! (L=0)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(22)) ! Bk
          nlanth = 22
          nDIMcf = 13 ! f8; multiplet J=L+S=3+3=0  => J = 7F_6
          lDIMCF = 7  ! (L=3)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(23)) ! Cf
          nlanth = 23
          nDIMcf = 16 ! f9; multiplet J=L+S=5+5/2=15/2  => J = 6H_15/2
          lDIMCF = 11 ! (L=5)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(24)) ! Es
          nlanth = 24
          nDIMcf = 17 ! f10; multiplet J=L+S=6+2=8  => J = 5I_8
          lDIMCF = 13 ! (L=6)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(25)) ! Fm
          nlanth = 25
          nDIMcf = 16 ! f11; multiplet J=L+S=6+3/2=15/2  => J = 4I_15/2
          lDIMCF = 13 ! (L=6)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(26)) ! Md
          nlanth = 26
          nDIMcf = 13 ! f12; multiplet J=L+S=5+1=6  => J = 3H_6
          lDIMCF = 11 ! (L=5)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(27)) ! No
          nlanth = 27
          nDIMcf = 8  ! f13; multiplet J=L+S=3+1/2=7/2  => J = 2F_7/2
          lDIMCF = 7  ! (L=3)
        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(28)) ! Lr
          nlanth = 28
          nDIMcf = 1  ! f14; multiplet J=L+S=0+0=0  => J = 1S_0
          lDIMCF = 1  ! (L=0)

        !------------------- transition metals ------------------------!
        case (clanth(29)) ! Sc

          nlanth = 29
          ! Sc2+ -- d^1
          read(u5,*,iostat=istatus) i_OxStat
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth
#         endif

          if (i_OxStat < 0) then
            write(u6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
            write(u6,'(A)') 'It was re-set to positive.'
            i_OxStat = abs(i_OxStat)
          end if
          if (i_OxStat < 2) then
            lDIMCF = 1 ! (L=0)
            write(u6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          else if (i_OxStat == 2) then
            lDIMCF = 5 ! (L=2) d^1
          else
            lDIMCF = 1 ! (L=0)
            write(u6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          end if

        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(30)) ! Ti

          nlanth = 30
          ! Ti2+ -- d^2
          ! Ti3+ -- d^1
          ! Ti4+ -- d^0
          read(u5,*,iostat=istatus) i_OxStat
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth
#         endif

          if (i_OxStat < 0) then
            write(u6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
            write(u6,'(A)') 'It was re-set to positive.'
            i_OxStat = abs(i_OxStat)
          end if
          if (i_OxStat < 2) then
            lDIMCF = 1 ! (L=0)
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          else if (i_OxStat == 2) then
            lDIMCF = 7 ! (L=3) d^2  3F
          else if (i_OxStat == 3) then
            lDIMCF = 5 ! (L=2) d^1  2D
          else
            lDIMCF = 1 ! (L=0) d^4
            write(u6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          end if

        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(31)) ! V

          nlanth = 31
          ! V2+ -- d^3
          ! V3+ -- d^2
          ! V4+ -- d^1
          read(u5,*,iostat=istatus) i_OxStat
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth
#         endif

          if (i_OxStat < 0) then
            write(u6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
            write(u6,'(A)') 'It was re-set to positive.'
            i_OxStat = abs(i_OxStat)
          end if
          if (i_OxStat < 2) then
            lDIMCF = 1 ! (L=0)
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          else if (i_OxStat == 2) then
            lDIMCF = 7 ! (L=3) d^3
          else if (i_OxStat == 3) then
            lDIMCF = 7 ! (L=3) d^2
          else if (i_OxStat == 4) then
            lDIMCF = 5 ! (L=2) d^1
          else
            lDIMCF = 1 ! (L=0) d^4
            write(u6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          end if

        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(32)) ! Cr

          nlanth = 32
          ! Cr3+ -- d^4
          read(u5,*,iostat=istatus) i_OxStat
          if (istatus < 0) call Error(2)

          if (i_OxStat < 0) then
            write(u6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
            write(u6,'(A)') 'It was re-set to positive.'
            i_OxStat = abs(i_OxStat)
          end if
          if (i_OxStat == 2) then
            lDIMCF = 5 ! (L=2) d^4
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
          else if (i_OxStat == 3) then
            lDIMCF = 7 ! (L=3) d^3
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
          else if (i_OxStat == 4) then
            lDIMCF = 7 ! (L=3) d^2
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
          else if (i_OxStat == 5) then
            lDIMCF = 5 ! (L=2) d^1
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
          else
            lDIMCF = 1 ! (L=0)
            write(u6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          end if
          write(u6,'(A)') 'Crystal field will not be computed'

        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(33)) ! Mn

          nlanth = 33
          ! Mn3+ -- d^4
          read(u5,*,iostat=istatus) i_OxStat
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth
#         endif

          if (i_OxStat < 0) then
            write(u6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
            write(u6,'(A)') 'It was re-set to positive.'
            i_OxStat = abs(i_OxStat)
          end if
          if (i_OxStat < 3) then
            lDIMCF = 1 ! (L=0)
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          else if (i_OxStat == 3) then
            lDIMCF = 5 ! (L=2) d^4
          else
            lDIMCF = 1 ! (L=0)
            write(u6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          end if

        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(34)) ! Fe

          nlanth = 34
          ! Co2+ -- d^6 or d^4
          read(u5,*,iostat=istatus) i_OxStat
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth
#         endif

          if (i_OxStat < 0) then
            write(u6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
            write(u6,'(A)') 'It was re-set to positive.'
            i_OxStat = abs(i_OxStat)
          end if
          if (i_OxStat < 2) then
            lDIMCF = 1 ! (L=0)
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          else if (i_OxStat == 2) then
            lDIMCF = 5 ! (L=2)  d^6  or  d^4
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
          else if (i_OxStat == 3) then
            lDIMCF = 1 ! (L=2)  d^5
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          else if (i_OxStat == 4) then
            lDIMCF = 1 ! (L=2)
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          else
            lDIMCF = 1 ! (L=0)
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          end if

        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(35)) ! Co

          nlanth = 35
          ! Co2+ -- d^7
          read(u5,*,iostat=istatus) i_OxStat
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth
#         endif

          if (i_OxStat < 0) then
            write(u6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
            write(u6,'(A)') 'It was re-set to positive.'
            i_OxStat = abs(i_OxStat)
          end if
          if (i_OxStat < 2) then
            lDIMCF = 1 ! (L=0)
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          else if (i_OxStat == 2) then
            lDIMCF = 7 ! (L=3) d^7
          else
            lDIMCF = 1 ! (L=0)
            write(u6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          end if

        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(36)) ! Ni

          nlanth = 36
          ! Ni2+ -- d^8
          read(u5,*,iostat=istatus) i_OxStat
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth
#         endif

          if (i_OxStat < 0) then
            write(u6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
            write(u6,'(A)') 'It was re-set to positive.'
            i_OxStat = abs(i_OxStat)
          end if
          if (i_OxStat < 2) then
            lDIMCF = 1 ! (L=0)
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          else if (i_OxStat == 2) then
            lDIMCF = 7 ! (L=2) d^8
          else
            lDIMCF = 1 ! (L=0)
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          end if

        !- - - - - - - - - - - - - - - - - - - -
        case (clanth(37)) ! Cu

          nlanth = 37
          ! Cu2+ -- d^9
          read(u5,*,iostat=istatus) i_OxStat
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth
#         endif

          if (i_OxStat < 0) then
            write(u6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',i_OxStat
            write(u6,'(A)') 'It was re-set to positive.'
            i_OxStat = abs(i_OxStat)
          end if
          if (i_OxStat < 2) then
            lDIMCF = 1 ! (L=0)
            write(u6,'(3A,i5)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          else if (i_OxStat == 2) then
            lDIMCF = 5 ! (L=2) d^9
          else
            lDIMCF = 1 ! (L=0) d^4
            write(u6,'(A)') 'Oxidation state of ',cME,' is:',i_OxStat
            write(u6,'(A)') 'Crystal field will not be computed'
          end if

        !- - - - - - - - - - - - - - - - - - - -
        case default
          write(u6,'(A)') 'Label of the metal is not understood.'
          write(u6,'(A)') 'Crystal field will not be computed'
      end select

      if (IPRINT > 2) then
        write(u6,'(5x,3A)') 'SINGLE_ANISO will calculate the parameters of the crystal field for Ln = ',clanth(nlanth),','
        write(u6,'(5x,A,I2,a)') 'for the ground multiplet J. Multiplicity of J = ',nDIMcf,' and'
        write(u6,'(5x,A,I2)') 'for the ground LS term. Multiplicity of L = ',lDIMcf
      end if
      LINENR = LINENR+1
    !-------------------------------------------
    case ('QUAX')
      !if (check_CRYS) then
      read(u5,*,iostat=istatus) axisoption
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'QUAX: axisoption =',axisoption
#     endif
      LINENR = LINENR+1

      if ((axisoption < 1) .or. (axisoption > 3)) &
        call WarningMessage(2,'QUAX: axisoption out of range! Calculation will continue by employing the default option.')
      if (axisoption == 3) then
        do j=1,3
          read(u5,*,iostat=istatus) (zmagn(i,j),i=1,3)
          if (istatus < 0) call Error(2)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'QUAX: zmagn(i,j) =',(zmagn(i,j),i=1,3)
#         endif
        end do
        LINENR = LINENR+3
      end if
      !else
      !  write(u6,'(A)') 'The CRYS keyword must be declared above QUAX in the input!'
      !end if
    !-------------------------------------------
    case ('PREX')
      LINENR = LINENR+1
    !-------------------------------------------
    case ('UBAR')
      compute_barrier = .true.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'UBAR:'
#     endif
      LINENR = LINENR+1
    !-------------------------------------------
    case ('ABCC')
      Do_structure_abc = .true.
      read(u5,*,iostat=istatus) (cryst(i),i=1,6)
      if (istatus < 0) call Error(2)

      do i=1,6
        if (cryst(i) <= 0) then
          call WarningMessage(2,'ABCC: zero or negative crystallographic parameters requested! ')
          call Quit_OnUserError()
        end if
      end do

#     ifdef _DEBUGPRINT_
      write(u6,*) 'ABCC: (cryst(i),i=1,6)=',(cryst(i),i=1,6)
#     endif
      read(u5,*,iostat=istatus) (coord(i),i=1,3)
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'ABCC: (coord(i),i=1,3)=',(coord(i),i=1,3)
#     endif
      LINENR = LINENR+2
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
    case ('ZEEM')
      zeeman_energy = .true.
      compute_magnetization = .true.

      read(u5,*,iostat=istatus) nDirZee
      if (istatus < 0) call Error(2)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'ZEEM: nDirZee=',nDirZee
#     endif

      do i=1,nDirZee
        ! open the zeeman_energy_xxx.txt file where Zeeman eigenstates will
        ! be further written in mangetization() subroutine
        write(namefile_energy,'(5A)') 'zeeman_energy_',char(48+mod(i/100,10)),char(48+mod(i/10,10)),char(48+mod(i,10)),'.txt'
        !print *, 'namefile_energy: ', namefile_energy
        LUZee(i) = IsFreeUnit(30+i)
        call molcas_open(LUZee(i),namefile_energy)
        !open(30+i,file=namefile_energy)

        read(u5,*,iostat=istatus) (dir_weight(i,l),l=1,3)
        if (istatus < 0) call Error(2)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'ZEEM: (dir_weight(i,l),l=1,3)=',(dir_weight(i,l),l=1,3)
#       endif

        check_dir_weight = sqrt(dir_weight(i,1)**2+dir_weight(i,2)**2+dir_weight(i,3)**2)

        if (abs(check_dir_weight-One) > 0.005_wp) then
          write(u6,'(A)') 'The directions for the magnetic field for the computation of the Zeeman splitting are wrong.'
          write(u6,'(A)') '( px^2 + py^2 + pz^2 ) must give 1.!'
          write(u6,'(A,I3,2x,A,F9.5)') 'In the present case for direction Nr.',i,' the dir_weight = px^2 + py^2 + pz^2 = ', &
                                       check_dir_weight**2
          LINENR = LINENR+2+i
          write(u6,*) ' READIN_SINGLE: Error reading standard input.'
          write(u6,*) ' SINGLE_ANISO input near line nr.',LINENR+1
          call ABEnd()
        end if

      end do
      LINENR = LINENR+nDirZee+1
    !-------------------------------------------
  end select
end do

if (IPRINT > 2) write(u6,'(5X,A)') 'NO ERROR WAS LOCATED WHILE READING INPUT'

if (compute_CF) then
  if (axisoption == 3) then
    ! check the determinant of the ZMAGN
    ZR(:,:) = zmagn(:,:)
    Det_zmagn = FindDetR(ZR,3)
    if (Det_zmagn < Zero) then
      write(u6,'(A)') 'QUAX: The determinant of the rotation matrix provided in the input is NEGATIVE.'
      write(u6,'(A,F22.14)') 'Determinant = ',Det_zmagn
      write(u6,'(A)') 'This means that the matrix you have provided can be decomposed in a product of two '
      write(u6,'(A)') 'matrices: Rotation*Inversion'
      write(u6,'(A)') 'The determinant of the Rotation matrix must be POSITIVE.'
      write(u6,'(A)') 'The program will stop.'
      return
    end if

    ! check the orthogonality of the ZMAGN:
    column_check(:,:) = Zero
    row_check(:,:) = Zero
    do j=1,3
      do L=1,3
        column_check(:,j) = column_check(:,j)+zmagn(:,L)*zmagn(j,L)
        row_check(:,j) = row_check(:,j)+zmagn(L,:)*zmagn(L,j)
      end do
    end do

    do i=1,3
      do j=i+1,3
        if (i == j) cycle
        if ((abs(column_check(1,2)) > 0.0001_wp) .or. (abs(column_check(1,3)) > 0.0001_wp) .or. &
            (abs(column_check(2,3)) > 0.0001_wp) .or. (abs(row_check(1,2)) > 0.0001_wp) .or. &
            (abs(row_check(1,3)) > 0.0001_wp) .or. (abs(row_check(2,3)) > 0.0001_wp)) then
          write(u6,'(A)') 'QUAX: The rotation matrix is not UNITARY.'
          write(u6,'(A,F19.14)') 'column_check(1,2) = ',column_check(1,2)
          write(u6,'(A,F19.14)') 'column_check(1,3) = ',column_check(1,3)
          write(u6,'(A,F19.14)') 'column_check(2,3) = ',column_check(2,3)
          write(u6,'(A,F19.14)') '   row_check(1,2) = ',row_check(1,2)
          write(u6,'(A,F19.14)') '   row_check(1,3) = ',row_check(1,3)
          write(u6,'(A,F19.14)') '   row_check(2,3) = ',row_check(2,3)
          write(u6,'(A)') 'All above values must be exact 0.0.'
          write(u6,'(A)') 'Or at least less than than 0.0001.'
          write(u6,'(A)') 'Did you employ enough digits for the rotation matrix?'
          write(u6,'(A)') 'The program will stop.'
          return
        end if
      end do
    end do
  end if ! axisoption
end if ! compute_CF

! preparing the info for computation of molar magnetization
ngrid = ngrid_map(ngrid)
if (compute_magnetization) then
  ! calculate the total number of directions for the average procedure
  nDirTot = order_table(nsymm,ngrid)
  if (zeeman_energy) nDirTot = nDirTot+nDirZee
  if (compute_Mdir_vector) nDirTot = nDirTot+nDir
end if

!------ CHECK the data from INPUT ------------------------------
!if (iprint > 10) then

# ifdef _DEBUGPRINT_
write(u6,'(A,  F9.5)') 'ZJPR :         = ',zJ
write(u6,'(A,  I3  )') 'PRLV :         = ',iprint
# endif

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
    write(u6,'(A,I3)') 'MLTP :         = ',NMULT
    if (NMULT <= 20) then
      write(u6,'(A,20I3)') '               = ',(NDIM(i),i=1,NMULT)
    else
      write(u6,'(A,20I3)') '               = ',(NDIM(i),i=1,20)
      do j=21,NMULT,20
        jEnd = min(NMULT,J+19)
        write(u6,'(A,20I3)') '                 ',(NDIM(i),i=j,jEnd)
      end do
    end if
  else
    write(u6,'(A)') 'MLTP :         =  No pseudospin Hamiltonians will be computed. Is MLTP defined?'
  end if
end if

if (compute_CF .and. (nDIMcf <= nss) .and. (lDIMcf <= nstate)) then
  if (nlanth < 15) then
    write(u6,'(3A)') 'The Crystal-Field acting on the ground atomic multiplet of Ln = ',clanth(nlanth),' is computed.'
  else if ((nlanth >= 15) .and. (nlanth < 29)) then
    write(u6,'(3A)') 'The Crystal-Field acting on the ground atomic multiplet of Ac = ',clanth(nlanth),' is computed.'
  else if (nlanth >= 29) then
    write(u6,'(3A)') 'The Crystal-Field acting on the ground atomic |L,ML> multiplet of TM = ',clanth(nlanth),' is computed.'
  end if

  write(u6,'(A,A )') 'CHIT :         = ',' molar magnetic susceptibility is computed'
  if (TINPUT) write(u6,'(A)') 'TEXP :         = the experimental temperature interval is read from the file "chitexp.input"'
  write(u6,'(A, I3)') 'TINT :      nT = ',nT
  write(u6,'(A,F7.3)') '          Tmin = ',Tmin
  write(u6,'(A,F7.3)') '          Tmax = ',Tmax

end if
!--------------------------------------------------------------------
if (compute_magnetization) then

  write(u6,'(A,A )') 'MAGN :         = ',' molar magnetization is computed'
  write(u6,'(A, I3)') 'NDIRTOT        = ',nDirTot
  write(u6,'(A, I3)') 'TMAG :         = ',nTempMagn
  write(u6,'(6x,A,20F7.3)') 'TempMagn = ',(TempMagn(i),i=1,nTempMagn)
  write(u6,'(A, I3)') 'HINT :      nH = ',nH
  write(u6,'(A,F7.3)') '          Hmin = ',Hmin
  write(u6,'(A,F7.3)') '          Hmax = ',Hmax
  write(u6,'(A, I3)') 'MAVE :   nDir = ',order_table(nsymm,ngrid)

  if (HINPUT) write(u6,'(A)') 'HEXP :         = the experimental field interval is read from the file "mexp.input"'
  if (encut_definition == 1) then
    write(u6,'(A, I3)') 'NCUT :         = ',ncut
  else if (encut_definition == 2) then
    write(u6,'(A,I4,a,i4)') 'ECUT :         = ',nk,', ',mg
  else if (encut_definition == 3) then
    write(u6,'(A,F7.3)') 'ERAT :         = ',encut_rate
  end if

  if (compute_Mdir_vector) then
    write(u6,'(A,20I3)') 'MVEC :         = ',nDir
    if (nDir > 0) then
      do i=1,nDir
        write(u6,'(A,I2,A,3F11.6)') '   Dir :',i,' : ',dirX(i),dirY(i),dirZ(i)
      end do
    end if
  end if
  if (zeeman_energy) then
    if (nDirZee == 1) then
      write(u6,'(2A,I2,1x,A)') 'ZEEM :         = ',' Zeeman splitting for the following direction of the '
      write(u6,'(18x,A)') 'applied magnetic field is given in the "zeeman_energy_xxx.txt" file in $WorkDir/'
    else if (nDirZee > 1) then
      write(u6,'(2A,I2,1x,A)') 'ZEEM :         = ',' Zeeman splitting for the following',nDirZee,' directions of the '
      write(u6,'(18x,A)') 'applied magnetic field are given in the "zeeman_energy_xxx.txt" files in $WorkDir/.'
    else
      write(u6,'(A)') 'Error in input processing. nDirZee<0!'
      call Quit_OnUserError()
    end if
    do i=1,nDirZee
      write(u6,'(17x,3F11.6)') (dir_weight(i,l),l=1,3)
    end do
  end if
end if ! magnetization
!--------------------------------------------------------------------

if (compute_torque) write(u6,'(A,A )') 'TORQ :         = ',' torque magnetization is computed'

if (doplot) write(u6,'(A,A )') 'PLOT :         = ',' GNUPLOT scripts and corresponding XT, M and UBAR plots will be generated'

if (Do_structure_abc) then
  write(u6,'(2A)') 'ABCC :         = ','the main magnetic axes for the computed pseudospins are written also in the '
  write(u6,'( A)') 'crystallographic "abc" axes'
  write(u6,'(10x,A,F9.4)') 'a       = ',cryst(1)
  write(u6,'(10x,A,F9.4)') 'b       = ',cryst(2)
  write(u6,'(10x,A,F9.4)') 'c       = ',cryst(3)
  write(u6,'(10x,A,F9.4)') 'alpha   = ',cryst(4)
  write(u6,'(10x,A,F9.4)') 'beta    = ',cryst(5)
  write(u6,'(10x,A,F9.4)') 'gamma   = ',cryst(6)
  write(u6,'(10x,a,3F9.4)') 'coords: = ',(coord(i),i=1,3)
end if

if (compute_g_tensors .and. compute_barrier) then
  Nblock = sum(ndim(1:nMult))
  write(u6,'(A,i4)') 'nBlock = ',nBlock
end if

return

contains

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  select case (code)
    case (1)
      write(u6,*) ' READIN_SINGLE: Unexpected End of input file.'
    case (2)
      write(u6,*) ' READIN_SINGLE: Error reading standard input.'
      write(u6,*) ' SINGLE_ANISO input near line nr.',LINENR+1
  end select
  call ABEnd()

end subroutine Error

end subroutine readin_single
