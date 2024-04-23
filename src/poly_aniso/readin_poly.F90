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
                       keopt,encut_definition,nK,mG,iopt,nP,AngPoints,ncut,LUZee,MxRank1,MxRank2,imaxrank,nsymm,ngrid,TempMagn, &
                       R_LG,R_ROT,Jex,JAex,JAex9,JDMex,JITOexR,JITOexI,tpar,upar,cryst,coord,Xfield,gtens_input,D_fact, &
                       EoverD_fact,riso,MagnCoords,thrs,tmin,tmax,hmin,hmax,Texp,chit_exp,Hexp,Mexp,encut_rate,zJ,dirX,dirY,dirZ, &
                       dir_weight,Title,itype,ifHDF,compute_g_tensors,compute_magnetization,TINPUT,HINPUT,Do_structure_abc,DoPlot, &
                       compute_Mdir_vector,zeeman_energy,m_paranoid,m_accurate,smagn,compute_susceptibility,decompose_exchange,KE, &
                       fitCHI,fitM,compute_torque,compute_barrier,Dipol,check_title,AnisoLines1,AnisoLines3,AnisoLines9, &
                       DM_exchange,JITO_exchange)
! THIS ROUTINE READS THE standard input.
! definition of the cluster:
!  nneq, neqv, neq, nCenter, ifHDF
! definition of the local metal sites
!  R_LG, R_ROT, gtens_input, D_fact, EoverD_fact, riso, itype
! definition of the exchange:
!  exch          : total number of exchange states
!  nPair         : number of metal pairs (number of interactions)
!  nexch         : exchange basis, nmax= MAX(nexch(:))
!  i_pair        : index of the metal site in a given interacting pair
!  Jex           : Lines exchange    ( 1 parameter / interacting pair)
!  JAex          : Anisotropic Lines ( 3 parameter / interacting pair)
!  JAex9         : Anisotropic Lines full ( 9 parameters / interacting pair)
!  JITO_exchange : options used in connection with ITO exchange:
!  JDMex, AnisoLines1, AnisoLines3, AnisoLines9, Dipol, DM_exchange, MxRank1, MxRank2, imaxrank, JITOexR, JITOexI
! options used in connection with KE
!  lant, KEOPT, multLn, KE, tpar, upar
! options used in connection with Dipol-Dipol interaction
!  MagnCoords
! definition of g and D tensors
!  nMult, nDim, compute_g_tensors
! definition of data for susceptibility
!  nT, tinput, compute_susceptibility, tmin, tmax, chit_exp, Texp
! options related to XT_MoverH
!  Xfield
! definition of data for magnetization:
!  nH, nTempMagn, iopt, TempMagn, Hexp, Mexp, thrs, hmin, hmax, hinput, compute_magnetization, compute_Mdir_vector, zeeman_energy,
!  m_paranoid, m_accurate, smagn
! options used to set up nM and EM
!  nK, mG     : encut_definition=1;
!  ncut       : encut_definition=2;
!  encut_rate : encut_definition=3;
!  encut_definition
! decompose exchange
!  decompose_exchange
! magnetization torque
!  nP, AngPoints, compute_torque
! Zeeman energy and M vector
!  nDir, nDirZee, LUZee, dirX, dirY, dirZ, dir_weight
! definition of mean field parameter
!  zJ
! definition of the crystal axes:
!  cryst ! a, b, c, alpha, beta, gamma
!  Do_structure_abc
! Cartesian coordinates of the main metal site, or center
!  coord
! definitions for blocking barrier
!  compute_barrier
! options for automatic fitting of parameters:
!  fitCHI : not used so far
!  fitM   : not used so far
! definition of print level
!  iPrint, check_title, Title, DoPlot

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three
use Definitions, only: wp, iwp, u5, u6

implicit none
integer(kind=iwp), intent(inout) :: nneq, neq(nneq), nT, nH, nTempMagn, nDir, nDirZee, nMult, nPair, nexch(nneq), nDim(nMult), &
                                    lant, multLn, iPrint, KEOPT, encut_definition, nK, mG, iopt, nP, AngPoints, ncut, &
                                    LUZee(nDirZee), imaxrank(npair,2), nsymm, ngrid
integer(kind=iwp), intent(in) :: neqv, exch, nCenter, MxRank1, MxRank2
integer(kind=iwp), intent(out) :: i_pair(nPair,2)
real(kind=wp), intent(inout) :: TempMagn(nTempMagn), R_LG(nneq,neqv,3,3), R_ROT(nneq,neqv,3,3), Jex(nPair), JAex(nPair,3), &
                                JAex9(nPair,3,3), JDMex(nPair,3), &
                                JITOexR(nPair,MxRank1,-MxRank1:MxRank1,MxRank2,-MxRank2:MxRank2), &
                                JITOexI(nPair,MxRank1,-MxRank1:MxRank1,MxRank2,-MxRank2:MxRank2), tpar, upar, cryst(6), coord(3), &
                                Xfield, gtens_input(3,nneq), D_fact(nneq), EoverD_fact(nneq), riso(3,3,nneq), MagnCoords(nneq,3), &
                                thrs, tmin, tmax, hmin, hmax, Texp(nT), chit_exp(nT), Hexp(nH), Mexp(nH,nTempMagn), encut_rate, &
                                zJ, dirX(nDir), dirY(nDir), dirZ(nDir), dir_weight(nDirZee,3)
character(len=180), intent(inout) :: Title
character, intent(inout) :: itype(nneq)
logical(kind=iwp), intent(inout) :: ifHDF, compute_g_tensors, compute_magnetization, Do_structure_abc, DoPlot, &
                                    compute_Mdir_vector, zeeman_energy, m_paranoid, m_accurate, smagn, compute_susceptibility, &
                                    decompose_exchange, KE, fitCHI, fitM, compute_torque, compute_barrier, Dipol, AnisoLines1, &
                                    AnisoLines3, AnisoLines9, DM_exchange, JITO_exchange
logical(kind=iwp), intent(out) :: tinput, hinput, check_title
#include "warnings.h"
integer(kind=iwp) :: ASUM, i, i1, i2, ic, icount_b_sites, inneq, iproj1, iproj2, irank1, irank2, istatus, j, jc, jproj1, jproj2, &
                     jrank1, jrank2, l, lb1, lb2, linenr, ll, lp, m, n, nst
real(kind=wp) :: check_dir_weight, detR, rsum, t1, t2, tmp, tmpR(3,3)
logical(kind=iwp) :: ab_initio_all, check_symm_presence, checktmag, encut_check, hcheck, nosym, tcheck
character(len=2) :: lanth
character(len=21) :: namefile_energy
character(len=288) :: ctmp, Line, string
integer(kind=iwp), allocatable :: duplicate_check(:), nind(:,:)
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: finddetr

#include "macros.fh"

check_title = .false.
icount_B_sites = 0
i_pair(:,:) = 0
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

#ifdef _DEBUGPRINT_
write(u6,'(A,  i6)') 'RDIN:      nneq=',nneq
write(u6,'(A,  i6)') 'RDIN:      neqv=',neqv
write(u6,'(A,  i6)') 'RDIN:     nPair=',nPair
!write(u6,'(A,99I6)') 'RDIN:  i_pair(:,1)=',(i_pair(i,1),i=1,nPair)
!write(u6,'(A,99I6)') 'RDIN:  i_pair(:,2)=',(i_pair(i,2),i=1,nPair)
write(u6,'(A,99i6)') 'RDIN:     neq()=',(neq(i),i=1,nneq)
write(u6,'(A,99i6)') 'RDIN:   nexch()=',(nexch(i),i=1,nneq)
#endif

!=========== End of default settings====================================
rewind(u5)
do
  read(u5,'(A72)',iostat=istatus) LINE
  if (istatus < 0) call Error(5)
# ifdef _DEBUGPRINT_
  write(u6,'(A)') LINE
# endif
  call NORMAL(LINE)
  if (LINE(1:5) == '&POLY') exit
end do
LINENR = 0

do
  call xFlush(u6)
  read(u5,'(A72)',iostat=istatus) LINE
  if (istatus < 0) call Error(5)
  LINENR = LINENR+1
  call NORMAL(LINE)
  if ((LINE(1:1) == '*') .or. (LINE == ' ')) cycle
  if (LINE(1:3) == 'END') exit
  select case (LINE(1:4))

    ! ------------ TITL -----------------------------------------------*
    case ('TITL')

      read(u5,*,iostat=istatus) ctmp
      if (istatus /= 0) call Error(4)

#     ifdef _DEBUGPRINT_
      write(u6,'(A)') ctmp
#     endif
      check_title = .true.
      Title = trim(ctmp)
      LINENR = LINENR+1

    ! ------------ OLDA -----------------------------------------------*
    case ('OLDA')
      LINENR = LINENR+1

    !---  process MLTP command ----------------------------------------*
    case ('MLTP')

      read(u5,*,iostat=istatus) NMULT
      if (istatus /= 0) call Error(4)

#     ifdef _DEBUGPRINT_
      write(u6,'(A,i4)') 'NMULT =',NMULT
#     endif

      read(u5,*,iostat=istatus) (NDIM(i),i=1,NMULT)
      if (istatus /= 0) call Error(4)

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

#     ifdef _DEBUGPRINT_
      write(u6,'(A,100i4)') 'NDIM: ',(NDIM(i),i=1,NMULT)
#     endif
      compute_g_tensors = .true.
      LINENR = LINENR+2

    !---  process TINT command ----------------------------------------*
    case ('TINT')
      compute_susceptibility = .true.
      if (TINPUT) then
        call Error(1)
      else
        TCHECK = .true.

        read(u5,*,iostat=istatus) t1,t2,nT
        if (istatus /= 0) call Error(4)

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

#       ifdef _DEBUGPRINT_
        write(u6,'(A,2ES15.7,i6)') 'Tmin, Tmax, nT: ',Tmin,Tmax,nT
#       endif
      end if
      LINENR = LINENR+1

    !---  process HINT command ----------------------------------------*
    case ('HINT')
      if (HINPUT) then
        call Error(2)
      else
        HCHECK = .true.
        compute_magnetization = .true.

        read(u5,*,iostat=istatus) t1,t2,nH
        if (istatus /= 0) call Error(4)

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

#       ifdef _DEBUGPRINT_
        write(u6,'(A,2ES15.7,i6)') 'Hmin, Hmax, nH: ',Hmin,Hmax,nH
#       endif
      end if
      LINENR = LINENR+1

    !---  process THRS command ----------------------------------------*
    case ('THRS')

      read(u5,*,iostat=istatus) THRS
      if (istatus /= 0) call Error(4)

      if (thrs < Zero) then
        call WarningMessage(2,'THRS: negative threshold for  average M!!! ')
        write(u6,'(A)') 'Set to default thrs=1.0e-10'
        thrs = 1.0e-10_wp
      end if

#     ifdef _DEBUGPRINT_
      write(u6,'(A,ES15.7)') 'THRS: ',THRS
#     endif
      LINENR = LINENR+1

    !---  process XFIE command ----------------------------------------*
    case ('XFIE')
      compute_susceptibility = .true.

      read(u5,*,iostat=istatus) tmp
      if (istatus /= 0) call Error(4)

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
#     ifdef _DEBUGPRINT_
      write(u6,'(A,ES15.7)') 'XFIE: ',xField
#     endif

      LINENR = LINENR+1

    !---  process PLOT command ----------------------------------------*
    case ('PLOT')
      DoPlot = .true.

#     ifdef _DEBUGPRINT_
      write(u6,'(A,L2)') 'PLOT: ',DoPlot
#     endif

      LINENR = LINENR+1

    !---  process IOPT command ----------------------------------------*
    case ('IOPT')

      read(u5,*,iostat=istatus) I  !option for computing MSUM and XTSUM
      if (istatus /= 0) call Error(4)

      if ((i < 0) .or. (i > 3)) then
        call WarningMessage(2,'IOPT: value out of range!!!')
        write(u6,'(A)') 'Set to default !'
        iopt = 1
      else
        iopt = i
      end if

#     ifdef _DEBUGPRINT_
      write(u6,'(A,I6)') 'IOPT: ',IOPT
#     endif
      LINENR = LINENR+1

    !---  process SMAG command ----------------------------------------*
    case ('SMAG')
      smagn = .true.
      compute_magnetization = .true.
#     ifdef _DEBUGPRINT_
      write(u6,'(A,L2)') 'SMAG: ',smagn
#     endif
      LINENR = LINENR+1

    !---  process MACC command ----------------------------------------*
    case ('MACC')
      compute_magnetization = .true.  ! request for computation of M(H)
      m_accurate = .true.             ! request for computation of M(H)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,L2)') 'MACC: ',m_accurate
#     endif
      LINENR = LINENR+1

    !---  process FITX command ----------------------------------------*
    case ('FITX')
      fitCHI = .true.
      LINENR = LINENR+1

    !---  process FITM command ----------------------------------------*
    case ('FITM')
      fitM = .true.
      LINENR = LINENR+1

    !---  process MPAR command ----------------------------------------*
    case ('MPAR')
      m_paranoid = .true.             ! request for computation of M(H)
      compute_magnetization = .true.
#     ifdef _DEBUGPRINT_
      write(u6,'(A,L2)') 'MPAR: ',m_paranoid
#     endif
      LINENR = LINENR+1

    !---  process TORQ command ----------------------------------------*
    case ('TORQ')

      compute_torque = .true.         ! request for computation of M(H)
      read(u5,*,iostat=istatus) i        ! number of angular points
      if (istatus /= 0) call Error(4)

      if (i <= 0) then
        call WarningMessage(2,'TORQ: nP value out of range!!!')
        write(u6,'(A)') 'Set to default !'
        nP = 45
      else
        nP = i
      end if

#     ifdef _DEBUGPRINT_
      write(u6,'(A)') 'TORQ: nP=',nP
#     endif
      AngPoints = nP+1
      LINENR = LINENR+1

    !---  process MAVE command ----------------------------------------*
    case ('MAVE')
      compute_magnetization = .true.  ! request for computation of M(H)

      read(u5,*,iostat=istatus) i,j  !nsymm, ngrid
      if (istatus /= 0) call Error(4)

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

#     ifdef _DEBUGPRINT_
      write(u6,'(2(A,i6))') ' nsymm:  ',nsymm,' ngrid:  ',ngrid
#     endif
      LINENR = LINENR+1

    !---  process NCUT command ----------------------------------------*
    case ('NCUT')
      if (ENCUT_check) then
        call Error(3)
      else
        ENCUT_check = .true.
        ! request for computation of M(H)
        compute_magnetization = .true.
        encut_definition = 1

        read(u5,*,iostat=istatus) i ! NCUT
        if (istatus /= 0) call Error(4)
        !E_cut = exchange_energy(Ncut)

        if ((i <= 0) .or. (i > exch)) then
          call WarningMessage(2,'NCUT: value out of range!!!')
          write(u6,'(A)') 'Set to full exchange basis.'
          nCUT = exch
        else
          nCUT = i
        end if
#       ifdef _DEBUGPRINT_
        write(u6,'(A,2i6)') 'ncut:  ',nCut
#       endif

        LINENR = LINENR+1
      end if

    !---  process ENCU command ----------------------------------------*
    case ('ENCU')
      if (ENCUT_check) then
        call Error(3)
      else
        ENCUT_check = .true.
        ! request for computation of M(H)
        compute_magnetization = .true.
        encut_definition = 2

        read(u5,*,iostat=istatus) NK,MG
        if (istatus /= 0) call Error(4)
        !E_cut = NK*K_Boltz+MG*mu_Bohr
#       ifdef _DEBUGPRINT_
        write(u6,'(A,2i6)') 'encu:  nK, mG=',NK,MG
#       endif

        LINENR = LINENR+1
      end if

    !---  process ERAT command ----------------------------------------*
    case ('ERAT')
      if (ENCUT_check) then
        call Error(3)
      else
        ENCUT_check = .true.
        ! request for computation of M(H)
        compute_magnetization = .true.
        encut_definition = 3

        read(u5,*,iostat=istatus) encut_rate
        if (istatus /= 0) call Error(4)
        !Ncut = int(nexch*encut_rate)
        !E_cut = E(Ncut)
#       ifdef _DEBUGPRINT_
        write(u6,'(A,i6)') 'encut_rate=',encut_rate
#       endif

        LINENR = LINENR+1
      end if

    !---  process ZJPR command ----------------------------------------*
    case ('ZJPR')
      read(u5,*,iostat=istatus) ZJ
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,ES18.10)') 'zJ    =',zJ
#     endif
      LINENR = LINENR+1

    !---  process PRLV command ----------------------------------------*
    case ('PRLV')
      read(u5,*,iostat=istatus) iPrint
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i6)') 'iPrint=',iPrint
#     endif
      LINENR = LINENR+1

    !---  process COOR command ----------------------------------------*
    case ('COOR')
      Dipol = .true.
#     ifdef _DEBUGPRINT_
      write(u6,'(A)') 'isite   MagnCoords:'
#     endif
      do i=1,nneq
        read(u5,*,iostat=istatus) (MagnCoords(i,l),l=1,3)
        if (istatus /= 0) call Error(4)
#       ifdef _DEBUGPRINT_
        write(u6,'(i3,5x,3ES18.10)') i,(MagnCoords(i,l),l=1,3)
#       endif
      end do
      LINENR = LINENR+nneq

    !---  process MVEC command ----------------------------------------*
    case ('MVEC')
      compute_magnetization = .true.  ! request for computation of M(H)
      compute_Mdir_vector = .true.

      read(u5,*,iostat=istatus) nDir
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i3)') 'nDir = ',nDir
#     endif

      do i=1,nDir
        read(u5,*,iostat=istatus) DirX(i),DirY(i),DirZ(i)
        if (istatus /= 0) call Error(4)
#       ifdef _DEBUGPRINT_
        write(u6,'(i3,5x,3ES18.10)') i,DirX(i),DirY(i),DirZ(i)
#       endif
      end do
      ! some processing:
      do i=1,nDir
        rsum = DirX(i)*DirX(i)+DirY(i)*DirY(i)+DirZ(i)*DirZ(i)
        if (rsum == Zero) then
          write(u6,'(a,i3,a)') 'error: MVEC  vector ',i,'has the modulus = 0.0.'
          write(u6,'(a     )') 'the program will stop now.'
          call quit(_RC_INPUT_ERROR_)
        end if
        if (abs(rsum-One) > 0.5e-13_wp) then
          write(u6,'(a,i3,a)') 'the vector ',i,'was re-normalized.'
          tmp = dirX(i)/sqrt(rsum)
          dirX(i) = tmp
          tmp = dirY(i)/sqrt(rsum)
          dirY(i) = tmp
          tmp = dirZ(i)/sqrt(rsum)
          dirZ(i) = tmp
        end if
      end do

      LINENR = LINENR+NDIR+1

    !---  process TEXP command ----------------------------------------*
    case ('TEXP')
      compute_susceptibility = .true.
      if (TCHECK) then
        call Error(1)
      else
        TINPUT = .true.

        read(u5,*,iostat=istatus) NT
        if (istatus /= 0) call Error(4)
#       ifdef _DEBUGPRINT_
        write(u6,'(A,i3)') 'nT = ',nT
#       endif

        do i=1,NT

          read(u5,*,iostat=istatus) texp(i),chit_exp(i)
          if (istatus /= 0) call Error(4)

          ! check and clean negative values:
          if (texp(i) < Zero) texp(i) = abs(texp(i))
          if (chit_exp(i) < Zero) chit_exp(i) = abs(chit_exp(i))
        end do
        tmin = texp(1)
        tmax = texp(nT)
      end if
      LINENR = LINENR+NT+1

    !---  process HEXP command ----------------------------------------*
    case ('HEXP')
      compute_magnetization = .true.
      if (checkTMAG) write(u6,'(A)') 'The data provided in TMAG will be ignored.'

      if (HCHECK) then
        call Error(2)
      else
        HINPUT = .true.

        read(u5,*) nTempMagn,(TempMagn(i),i=1,nTempMagn)
        read(u5,*) nH
#       ifdef _DEBUGPRINT_
        write(u6,*) 'HEXP: nTempMagn =',nTempMagn
        write(u6,*) 'HEXP: TempMagn()=',(TempMagn(i),i=1,nTempMagn)
        write(u6,*) 'HEXP: nH        =',nH
#       endif

        if (nH < 0) nH = abs(nH)
        if (nH == 0) call Quit_OnUserError()

        do i=1,nH
          read(u5,*,iostat=istatus) Hexp(i),(Mexp(i,j),j=1,nTempMagn)
          if (istatus /= 0) call Error(4)
          ! check and clean negative values:
          if (hexp(i) < Zero) hexp(i) = abs(hexp(i))
          do j=1,nTempMagn
            if (Mexp(i,j) < Zero) Mexp(i,j) = abs(Mexp(i,j))
          end do
        end do
        hmin = hexp(1)
        hmax = hexp(nH)
      end if
      LINENR = LINENR+NH+2

    !---  process TMAG command ----------------------------------------*
    case ('TMAG')
      if (.not. HINPUT) then
        compute_magnetization = .true.
        checkTMAG = .true.

        read(u5,*,iostat=istatus) nTempMagn,(TempMagn(i),i=1,nTempMagn)
        if (istatus /= 0) call Error(4)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'TMAG: nTempMagn =',nTempMagn
        write(u6,*) 'TMAG: TempMagn()=',(TempMagn(i),i=1,nTempMagn)
#       endif

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

    !---  process NNEQ command ----------------------------------------*
    ! this is the most important keyword for Poly_Aniso
    case ('NNEQ')
      ! number of non-equivalent centers; type of all centers
      read(u5,*,iostat=istatus) NNEQ,ab_initio_all,ifHDF
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i4,A,L2,A,L2)') 'NNEQ=',NNEQ,' ab_initio_all=',ab_initio_all,' ifHDF=',ifHDF
#     endif
      ! number of equivalent centers of type "i"
      read(u5,*,iostat=istatus) (NEQ(i),i=1,Nneq)
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,100I4)') 'NEQ(I)=',(NEQ(i),i=1,nneq)
#     endif
      ! number of RASSI wf for exchange
      read(u5,*,iostat=istatus) (Nexch(i),i=1,Nneq)
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,100I4)') 'NExch(I)=',(NExch(i),i=1,nneq)
#     endif

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
        if (neq(i) > 1) nosym = .false.
      end do
      !write(u6,'(A,i5)') 'exch = ',exch
      if (exch == 1) then
        write(u6,'(3/)')
        write(u6,'(A)') repeat('#',100)
        write(u6,'(3/)')
        write(u6,'(A)') 'The size of the exchange matrix is 1. Is this really what you intended to compute?'
        write(u6,'(A)') 'The program will continue...'
        write(u6,'(3/)')
        write(u6,'(A)') repeat('#',100)
        write(u6,'(3/)')
      end if

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
        read(u5,*,iostat=istatus) (itype(i),i=1,Nneq)
        if (istatus /= 0) call Error(4)
#       ifdef _DEBUGPRINT_
        write(u6,'(A,100A3)') 'itype: ',(itype(i),i=1,nneq)
        call xFlush(u6)
#       endif
        icount_B_sites = 0
        do i=1,nneq
          if ((itype(i) == 'B') .or. (itype(i) == 'C')) then
            icount_B_sites = icount_B_sites+1
            read(u5,*,iostat=istatus) (gtens_input(l,i),l=1,3),D_fact(i),EoverD_fact(i)
            if (istatus /= 0) call Error(4)
#           ifdef _DEBUGPRINT_
            write(u6,'(A,i4,A,3ES20.10, 2(A,ES20.10) )') 'gtens_input(',i,')=',(gtens_input(l,i),l=1,3),' D = ',D_fact(i), &
                                                         ' E/D =',EoverD_fact(i)
#           endif
            if ((itype(i) == 'C') .and. &
                ((gtens_input(1,i) /= gtens_input(2,i)) .or. (gtens_input(1,i) /= gtens_input(3,i)) .or. &
                 (gtens_input(2,i) /= gtens_input(3,i)))) then
              do ic=1,3
                read(u5,*,iostat=istatus) (riso(jc,ic,i),jc=1,3)
                if (istatus /= 0) call Error(4)
              end do

            else
              call unitmat(riso(:,:,i),3)
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
        itype(:) = 'A'
      end if !ab_initio_all
      LINENR = LINENR+3+icount_B_sites
#     ifdef _DEBUGPRINT_
      call xFlush(u6)
#     endif

    !---  process LIN9 command ----------------------------------------*
    case ('LIN9')
      AnisoLines9 = .true.

      read(u5,*,iostat=istatus) npair
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i6)') 'LIN9:  nPair=',nPair
#     endif

      do i=1,npair
        ! the convention for 9 exchange interactions is
        ! Jxx, Jxy, Jxz, Jyx, Jyy, Jyz, Jzx, Jzy, Jzz
        read(u5,*,iostat=istatus) i_pair(i,1),i_pair(i,2),(JAex9(i,1,j),j=1,3),(JAex9(i,2,j),j=1,3),(JAex9(i,3,j),j=1,3)
        if (istatus /= 0) call Error(4)
#       ifdef _DEBUGPRINT_
        write(u6,'(A,2I3,9F14.8)') 'LIN9: ',i_pair(i,1),i_pair(i,2),(JAex9(i,1,j),j=1,3),(JAex9(i,2,j),j=1,3),(JAex9(i,3,j),j=1,3)
#       endif
      end do
      LINENR = LINENR+npair+1

    !---  process LIN3 command ----------------------------------------*
    case ('ALIN','LIN3')
      AnisoLines3 = .true.

      read(u5,*,iostat=istatus) npair
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i6)') 'nPair=',nPair
#     endif

      do i=1,npair
        ! the convention for 3 exchange interactions is
        ! Jxx, Jyy, Jzz
        read(u5,*,iostat=istatus) i_pair(i,1),i_pair(i,2),(JAex(i,j),j=1,3)
        if (istatus /= 0) call Error(4)
#       ifdef _DEBUGPRINT_
        write(u6,'(A,2i3,3F14.8)') 'ALIN/LIN3: ',i_pair(i,1),i_pair(i,2),(JAex(i,j),j=1,3)
#       endif
      end do
      LINENR = LINENR+npair+1

    !---  process JITO command ----------------------------------------*
    case ('ITOJ')
      JITO_exchange = .true.

      read(u5,*,iostat=istatus) npair
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i6)') 'nPair=',nPair
#     endif
      JITOexR(:,:,:,:,:) = Zero
      JITOexI(:,:,:,:,:) = Zero

      do i=1,npair
        ! the convention for exchange interactions is

        read(u5,*,iostat=istatus) i_pair(i,1),i_pair(i,2),imaxrank(i,1),imaxrank(i,2)
        if (istatus /= 0) call Error(4)
        do irank1=1,imaxrank(i,1),2
          do iproj1=-irank1,irank1
            do irank2=1,imaxrank(i,2),2
              do iproj2=-irank2,irank2
                read(u5,*,iostat=istatus) jrank1,jproj1,jrank2,jproj2,JITOexR(i,jrank1,jproj1,jrank2,jproj2), &
                                          JITOexI(i,jrank1,jproj1,jrank2,jproj2)
                if (istatus /= 0) call Error(4)
              end do
            end do
          end do
        end do

#       ifdef _DEBUGPRINT_
        write(u6,'(A,I3)') 'ITO Exchange parameters for pair:',i
        do jrank1=1,imaxrank(i,1),2
          do jproj1=-jrank1,jrank1
            do jrank2=1,imaxrank(i,2),2
              do jproj2=-jrank2,jrank2
                write(u6,'(4I3,2x,2ES21.14)') jrank1,jproj1,jrank2,jproj2,JITOexR(i,jrank1,jproj1,jrank2,jproj2), &
                                              JITOexI(i,jrank1,jproj1,jrank2,jproj2)
              end do
            end do
          end do
        end do
#       endif
      end do
      LINENR = LINENR+npair+1

    !---  process PAIR = LIN1 command ---------------------------------*
    case ('PAIR','LIN1')
      AnisoLines1 = .true.

      read(u5,*,iostat=istatus) npair
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i6)') 'nPair=',nPair
#     endif

      do i=1,npair

        read(u5,*,iostat=istatus) i_pair(i,1),i_pair(i,2),Jex(i)
        if (istatus /= 0) call Error(4)
#       ifdef _DEBUGPRINT_
        write(u6,'(i4,2x,2I4,2x,ES18.10)') i,i_pair(i,1),i_pair(i,2),Jex(i)
#       endif

      end do
      LINENR = LINENR+npair+1

    !---  process DMEX command ----------------------------------------*
    case ('DMEX')
      DM_exchange = .true.
#     ifdef _DEBUGPRINT_
      write(u6,'(A,L2)') 'DMEX::  DM_exchange=',DM_exchange
#     endif

      read(u5,*,iostat=istatus) npair
      if (istatus /= 0) call Error(4)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i6)') 'nPair=',nPair
#     endif

      do i=1,npair
        ! the convention for 3 exchange interactions is
        ! JDMex(x, y, z)
        read(u5,*,iostat=istatus) i_pair(i,1),i_pair(i,2),(JDMex(i,j),j=1,3)
        if (istatus /= 0) call Error(4)
      end do
      LINENR = LINENR+npair+1
      !LINENR = LINENR+1

    !---  process ZEEM command ----------------------------------------*
    case ('ZEEM')
      zeeman_energy = .true.
      compute_magnetization = .true.
      LUZEE(:) = 0

      read(u5,*,iostat=istatus) nDirZee
      if (istatus /= 0) call Error(4)

      do i=1,nDirZee
        ! open the zeeman_energy_xxx.txt file where Zeeman eigenstates will
        ! be further written in mangetization() subroutine
        write(namefile_energy,'(5A)') 'zeeman_energy_',char(48+mod(i/100,10)),char(48+mod(i/10,10)),char(48+mod(i,10)),'.txt'
#       ifdef _DEBUGPRINT_
        write(u6,'(2A)') 'namefile_energy: ',namefile_energy
#       endif
        LUZee(i) = IsFreeUnit(30+i)
        call molcas_open(LUZee(i),namefile_energy)

        read(u5,*,iostat=istatus) (dir_weight(i,l),l=1,3)
        if (istatus /= 0) call Error(4)

        check_dir_weight = sqrt(dir_weight(i,1)**2+dir_weight(i,2)**2+dir_weight(i,3)**2)

        if (abs(check_dir_weight-One) > 0.005_wp) then
          write(u6,'(A)') 'The directions for the magnetic field for the computation of the Zeeman splitting are wrong.'
          write(u6,'(A)') '( px^2 + py^2 + pz^2 ) must give 1.!'
          write(u6,'(A,I3,2x,A,F9.5)') 'In the present case for direction Nr.',i,' the dir_weight = px^2 + py^2 + pz^2 = ', &
                                       check_dir_weight**2
          LINENR = LINENR+2+i
          call Error(4)
        end if

      end do
      LINENR = LINENR+nDirZee+1

    !---  process MAGN command ----------------------------------------*
    !case ('MAGN')
    !  compute_magnetization = .true.

    !---  process SYMM command ----------------------------------------*
    case ('SYMM')
      nosym = .false.
      check_symm_presence = .true.
      R_lg(:,:,:,:) = Zero
      R_rot(:,:,:,:) = Zero
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i6)') 'SYMM - at init'
#     endif
      ll = 0
      do i=1,nneq
#       ifdef _DEBUGPRINT_
        write(u6,'(A,i6)') 'SYMM:  i=',i
#       endif

        read(u5,*,iostat=istatus) inneq
        if (istatus /= 0) call Error(4)
        unused_var(inneq)
#       ifdef _DEBUGPRINT_
        write(u6,'(A,i6)') 'inneq=',inneq
#       endif

        do j=1,Neq(i)
#         ifdef _DEBUGPRINT_
          write(u6,'(A,i6)') 'SYMM:  j=',j
#         endif
          do m=1,3
            ll = ll+1
            read(u5,*,iostat=istatus) (R_lg(i,j,m,n),n=1,3)
            if (istatus /= 0) call Error(4)
#           ifdef _DEBUGPRINT_
            write(u6,'(3ES20.12)') (R_lg(i,j,m,n),n=1,3)
#           endif
          end do
        end do

        do j=1,neq(i)
          tmpR(:,:) = R_lg(i,j,:,:)

          detR = FindDetR(tmpR,3)
#         ifdef _DEBUGPRINT_
          write(u6,'(A,3ES20.12)') 'SYMM:  detR=',detR
#         endif

          if (abs(abs(detR)-One) > 0.001_wp) then
            write(u6,'(A)') 'The rotation matrices must be UNITARY.'
            write(u6,'(A)') 'and of RIGHT hand system'
            write(u6,'(A,F11.6)') 'DET = ',detR
            LINENR = LINENR+ll+i+1
            call Error(4)
          end if
          if (detR < Zero) then
            R_rot(i,j,:,:) = -R_lg(i,j,:,:)
          else if (detR > Zero) then
            R_rot(i,j,:,:) = R_lg(i,j,:,:)
          end if
        end do !neq(i)
      end do !nneq
      LINENR = LINENR+nneq+ll

    !---  process EXCH command ----------------------------------------*
    case ('EXCH')
      decompose_exchange = .true.

    !---  process OLDA command ----------------------------------------*
    !case ('OLDA')
    !  old_aniso_format = .true.

    !---  process END command -----------------------------------------*
    !case ('END ')

    !---  process UBAR command ----------------------------------------*
    case ('UBAR')
      compute_barrier = .true.
      !read(u5,*,iostat=istatus) icase
      !if (istatus /= 0) call Error(4)
      !if (icase == 1) then
      !! icase =1  --> magnetic field is applied along the main magnetic
      !!               axis of each doublet (multiplet)
      !  continue
      !else if (icase == 2) then
      !! icase =2  --> magnetic field is applied along the main magnetic
      !!               axis of the specified doublet (number) number (NDim)
      !  read(u5,*,iostat=istatus) NmagMult
      !  if (istatus /= 0) call Error(4)
      !else if (icase == 3) then
      !! icase =3  --> a new coordination system is defined by the user,
      !!               and the magnetic field is applied along gZ (third axis)
      !  do i=1,3
      !    read(u5,*,iostat=istatus) (uBar_Rot(i,j),j=1,3)
      !    if (istatus /= 0) call Error(4)
      !  end do
      !  LINENR = LINENR+3
      !else
      !  write(u6,'(A)') 'Is the UBAR keyword used correctly?'
      !  write(u6,'(A)') 'The ICASE parameter is not understood.'
      !end if
      LINENR = LINENR+1

    !------------------------------------------------------------------*
    case ('ABCC')
      Do_structure_abc = .true.
      read(u5,*,iostat=istatus) (cryst(i),i=1,6)
      if (istatus /= 0) call Error(4)
      coord(:) = Zero
      !read(u5,*,iostat=istatus) (coord(i),i=1,3)
      !if (istatus /= 0) call Error(4)
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

    !---  process LONG command ----------------------------------------*
    case ('LONG')
      KE = .true.

      read(u5,*,iostat=istatus) lanth,tpar,upar,KEOPT
      if (istatus /= 0) call Error(4)
      call UpCase(lanth)

      select case (lanth)
        case ('GD')
          lant = 1
          multLn = 8
        case ('TB')
          lant = 2
          multLn = 13
        case ('DY')
          lant = 3
          multLn = 16
        case ('HO')
          lant = 4
          multLn = 17
        case ('ER')
          lant = 5
          multLn = 16
        case ('TM')
          lant = 6
          multLn = 13
        case ('YB')
          lant = 7
          multLn = 8
        case default
          write(u6,'( A)') 'Error in getting the type of the lanthanide!'
          write(u6,'(2A)') 'The program has this info: lanth =',lanth
          write(u6,'(26x,A,i2)') 'multLn =',multLn
      end select
      LINENR = LINENR+2

  end select
end do

! end of reading input keywords

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
      exit
    end if
  end do
end if
#ifdef _DEBUGPRINT_
write(u6,*) 'READIN_POLY:  after proc g and D'
call xFlush(u6)
#endif

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
        lp = lp+1
        i_pair(lp,1) = i
        i_pair(lp,2) = j
#       ifdef _DEBUGPRINT_
        write(u6,'(A,i3,A,2I3)') 'lp=',lp,' i_pair(lp,1:2)=',i_pair(lp,1),i_pair(lp,2)
#       endif
      end do
    end do
  end if

  call mma_allocate(Duplicate_check,nPair,label='Duplicate_check')
  do i=1,npair
    Duplicate_check(i) = 1000*i_pair(i,1)+i_pair(i,2)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Duplicate_check: ',Duplicate_check(i)
#   endif

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
# ifdef _DEBUGPRINT_
  call xFlush(u6)
# endif

  ! check on the duplicate
  do i=1,npair
    do j=i+1,npair
      if (Duplicate_check(i) == Duplicate_check(j)) then
        write(u6,'(A)') 'Some interactions are declared twice in the input. Please declare all interactions only once!'
        write(u6,'(A)') 'The program has to stop.'
        call quit(_RC_INPUT_ERROR_)
      end if
    end do
  end do
  call mma_deallocate(Duplicate_check)

  ! check on the indices of the exchange couplings:
  do i=1,nPair
    if ((i_pair(i,1) > nCenter) .or. (i_pair(i,2) > nCenter)) then
      write(u6,'(A)') 'The numbering of the magnetic centers within NPAIR keyword is wrong.'
      write(u6,'(A,i2,a)') 'For the interaction Nr. = ',i, &
                           ' you numerate individual centers with numbers larger than the total number of centers in this molecule.'
      write(u6,'(A)') 'NPAIR keyword is wrong.'
      write(u6,'(A,I4)') 'i_pair(i,1) =',i_pair(i,1)
      write(u6,'(A,I4)') 'i_pair(i,2) =',i_pair(i,2)
      write(u6,'(A)') 'The program has to stop.'
      call quit(_RC_INPUT_ERROR_)
    end if
  end do

  ! check on the size of MxRank1 and MxRank2 wrt nexch(1) and nexch(2)
  if (JITO_exchange) then
    call mma_allocate(nind,exch,2,label='nind')
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
        write(u6,'(A)') repeat('#',100)
        write(u6,'(A)') 'interacting pair = ',lp
        write(u6,'(A)') 'type of site 1 = ',i1
        write(u6,'(A,i2,A,i2)') 'nExch(',i1,') = ',nExch(i1)
        write(u6,'(A,i2,A,i2)') 'Rank2(',lp,',1) = ',imaxrank(lp,1)
        write(u6,'(A)') 'nExch < Rank+1 !!!'
        write(u6,'(A)') 'Rank of ITO operators for site 1 is larger than the number of defined exchange states for this site'
        write(u6,'(A)') 'The program will use only parameters which bring non-zero contribution to exchange.'
        write(u6,'(A)') repeat('#',100)
      end if
      if (nExch(i2) < (imaxrank(lp,2)+1)) then
        write(u6,'(A)') repeat('#',100)
        write(u6,'(A)') 'interacting pair = ',lp
        write(u6,'(A)') 'type of site 2 = ',i2
        write(u6,'(A,i2,A,i2)') 'nExch(',i2,') = ',nExch(i2)
        write(u6,'(A,i2,A,i2)') 'Rank2(',lp,',2) = ',imaxrank(lp,2)
        write(u6,'(A)') 'nExch < Rank+1 !!!'
        write(u6,'(A)') 'Rank of ITO operators for site 2 is larger than the number of defined exchange states for this site'
        write(u6,'(A)') 'The program will use only parameters which bring non-zero contribution to exchange.'
        write(u6,'(A)') repeat('#',100)
      end if
    end do
    call mma_deallocate(nind)
  end if ! JITO_exchange
end if ! nPair

#ifdef _DEBUGPRINT_
write(u6,*) 'READIN_POLY:  before 200 '
call xFlush(u6)
#endif
!--------  definition of exchange --------------------------------------

!-----------------------------------------------------------------------
!write(u6,*) ' The following input line was not understood:'
!write(u6,'(A)') LINE

! ===============   NORMAL ENDING  =====================================

if (IPRINT > 2) write(u6,'(5X,A)') 'NO ERROR WAS LOCATED WHILE READING INPUT'

return

contains

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  select case (code)
    case (1)
      write(u6,*) 'READIN_POLY: THE TINT command is incompatible with TEXP'
    case (2)
      write(u6,*) 'READIN_POLY: THE HINT command is incompatible with HEXP'
    case (3)
      write(u6,*) 'READIN_POLY: THE NCUT, ENCU, and ERAT are mutually exclusive. '// &
                  'You cannot use more than one keyword at the same time.'
    case (4)
      write(u6,*) ' READIN_POLY: Error reading "poly_aniso.input" '
      write(u6,*) ' near line nr.',LINENR+1
    case (5)
      write(u6,*) ' READIN_POLY: Unexpected End of input file.'
  end select
  call quit(_RC_INPUT_ERROR_)

end subroutine Error

end subroutine Readin_poly
