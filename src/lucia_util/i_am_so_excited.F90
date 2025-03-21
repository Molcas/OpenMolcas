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
! Copyright (C) 2015, Lasse Kragh Soerensen                            *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine I_AM_SO_EXCITED(NBATCH,IBATCH,LBATCH,I1BATCH)
! Subroutine by Lasse from October 2015
!
! Updated March 2018 for doubly excited states
!
! Will give single excited states in from the desired GAS (or GAS's)
! And now also doubly excited states
!
! The difference between the HEXS and DEXS is controlled by
! I_ELIMINATE_GAS
! Notice that for the doubly excited states (DEXS) all
! singly excited states (HEXS) are in effect.

use lucia_data, only: I2ELIMINATED_IN_GAS, I_AM_OUT, I_ELIMINATE_GAS, IBSPGPFTP, IELIMINATED_IN_GAS, ISPGPFTP, MXPSTT, &
                      N_2ELIMINATED_GAS, N_ELIMINATED_BATCHES, N_ELIMINATED_GAS, NELFGP, NGAS, NSPGPFTP
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NBATCH, IBATCH(8,*), LBATCH(NBATCH), I1BATCH(NBATCH)
integer(kind=iwp) :: I, IBLOCK, IEL, IGAS, IGAS_ELIM, IITYPE, IMATCH_ALPHA, IMATCH_ALPHAM1, IMATCH_BETA, IMATCH_BETAM1, &
                     IMATCH_BLOCK, IMAX_OCC(2,NGAS,2), IOFF, ISPGP, ITYPE_A, ITYPE_B, J, JBATCH, MAX_E_GAS_ALPHA(2,MXPSTT), &
                     MAX_E_GAS_BETA(2,MXPSTT), MAXM1_E_GAS_ALPHA(2,MXPSTT), MAXM1_E_GAS_BETA(2,MXPSTT), NALPHA, NALPHAM1, NBETA, &
                     NBETAM1

#ifdef _DEBUGPRINT_
write(u6,*) ' Oh I am so excited'
write(u6,*)
write(u6,*) ' Number of GAS without max (max-1) occupation = ',N_ELIMINATED_GAS+N_2ELIMINATED_GAS
write(u6,*)
if ((I_ELIMINATE_GAS == 1) .or. (I_ELIMINATE_GAS == 3)) then
  write(u6,*) ' GAS without maximum occupation (HEXS)'
  write(u6,*)
  do I=1,N_ELIMINATED_GAS
    write(u6,*) IELIMINATED_IN_GAS(I)
  end do
end if
if (I_ELIMINATE_GAS > 1) then
  write(u6,*) ' GAS without maximum-1 occupation (DEXS)'
  write(u6,*)
  do I=1,N_2ELIMINATED_GAS
    write(u6,*) I2ELIMINATED_IN_GAS(I)
  end do
end if
#endif

! First we need to find the GAS spaces for which we will eliminate
! the maximum occupation.

IMAX_OCC = 0

do JBATCH=1,min(NBATCH,2) ! only alpha and beta
  do IBLOCK=I1BATCH(JBATCH),I1BATCH(JBATCH)+LBATCH(JBATCH)-1
    do ISPGP=1,NSPGPFTP(JBATCH)
      IOFF = IBSPGPFTP(JBATCH)
      do IGAS=1,NGAS
        IITYPE = ISPGPFTP(IGAS,IOFF-1+ISPGP)
        IEL = NELFGP(IITYPE)
        if (IEL > IMAX_OCC(JBATCH,IGAS,1)) then
          IMAX_OCC(JBATCH,IGAS,1) = IEL
          IMAX_OCC(JBATCH,IGAS,2) = IITYPE
        end if
      end do
    end do
  end do
end do

#ifdef _DEBUGPRINT_
do JBATCH=1,min(NBATCH,2)
  if (JBATCH == 1) then
    write(u6,*) ' Maximum number of alpha electrons in each GAS'
  else
    write(u6,*) ' Maximum number of beta electrons in each GAS'
  end if
  write(u6,*)
  write(u6,*) ' GAS, Electrons, Group'
  do IGAS=1,NGAS
    write(u6,*) IGAS,IMAX_OCC(JBATCH,IGAS,1),IMAX_OCC(JBATCH,IGAS,2)
  end do
  write(u6,*)
end do
#endif

! Find which types contains the groups with a maximum number
! of alpha or beta electrons in a GAS

NALPHA = 0
NBETA = 0
NALPHAM1 = 0
NBETAM1 = 0
do JBATCH=1,min(NBATCH,2) ! only alpha and beta
  do IBLOCK=I1BATCH(JBATCH),I1BATCH(JBATCH)+LBATCH(JBATCH)-1
    do ISPGP=1,NSPGPFTP(JBATCH)
      IOFF = IBSPGPFTP(JBATCH)
      do IGAS=1,NGAS
        IITYPE = ISPGPFTP(IGAS,IOFF-1+ISPGP)
        IEL = NELFGP(IITYPE)
        if (IEL == IMAX_OCC(JBATCH,IGAS,1)) then
          if (JBATCH == 1) then
            NALPHA = NALPHA+1
            MAX_E_GAS_ALPHA(1,NALPHA) = IGAS
            MAX_E_GAS_ALPHA(2,NALPHA) = ISPGP
          else
            NBETA = NBETA+1
            MAX_E_GAS_BETA(1,NBETA) = IGAS
            MAX_E_GAS_BETA(2,NBETA) = ISPGP
          end if
        end if
        if (I_ELIMINATE_GAS > 1) then ! DEXS
          if (IEL == IMAX_OCC(JBATCH,IGAS,1)-1) then
            if (JBATCH == 1) then
              NALPHAM1 = NALPHAM1+1
              MAXM1_E_GAS_ALPHA(1,NALPHAM1) = IGAS
              MAXM1_E_GAS_ALPHA(2,NALPHAM1) = ISPGP
            else
              NBETAM1 = NBETAM1+1
              MAXM1_E_GAS_BETA(1,NBETAM1) = IGAS
              MAXM1_E_GAS_BETA(2,NBETAM1) = ISPGP
            end if
          end if
        end if
      end do
    end do
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) 'Maximum number of alpha supergroups that can be eliminated',NALPHA+NALPHAM1
write(u6,*)
write(u6,*) ' GAS Supergroup for HEXS'
write(u6,*)
do IGAS=1,NGAS
  do I=1,NALPHA
    if (MAX_E_GAS_ALPHA(1,I) == IGAS) write(u6,*) MAX_E_GAS_ALPHA(1,I),MAX_E_GAS_ALPHA(2,I)
  end do
end do
if (I_ELIMINATE_GAS > 1) then ! DEXS
  write(u6,*)
  write(u6,*) ' GAS Supergroup for DEXS'
  write(u6,*)
  do IGAS=1,NGAS
    do I=1,NALPHAM1
      if (MAXM1_E_GAS_ALPHA(1,I) == IGAS) write(u6,*) MAXM1_E_GAS_ALPHA(1,I),MAXM1_E_GAS_ALPHA(2,I)
    end do
  end do
end if
write(u6,*)
write(u6,*) 'Maximum number of beta supergroups that can be eliminated',NBETA+NBETAM1
write(u6,*)
write(u6,*) ' GAS Supergroup for HEXS'
write(u6,*)
do IGAS=1,NGAS
  do I=1,NBETA
    if (MAX_E_GAS_BETA(1,I) == IGAS) write(u6,*) MAX_E_GAS_BETA(1,I),MAX_E_GAS_BETA(2,I)
  end do
end do
if (I_ELIMINATE_GAS > 1) then ! DEXS
  write(u6,*)
  write(u6,*) ' GAS Supergroup for DEXS'
  write(u6,*)
  do IGAS=1,NGAS
    do I=1,NBETAM1
      if (MAXM1_E_GAS_BETA(1,I) == IGAS) write(u6,*) MAXM1_E_GAS_BETA(1,I),MAXM1_E_GAS_BETA(2,I)
    end do
  end do
end if
write(u6,*)
#endif

! Now find the batches to possibly eliminate

N_ELIMINATED_BATCHES = 0

do JBATCH=1,NBATCH
  do IBLOCK=I1BATCH(JBATCH),I1BATCH(JBATCH)+LBATCH(JBATCH)-1
    ITYPE_A = IBATCH(1,IBLOCK)
    ITYPE_B = IBATCH(2,IBLOCK)
    IMATCH_BLOCK = 0
    if ((I_ELIMINATE_GAS == 1) .or. (I_ELIMINATE_GAS == 3)) then
      ! HEXS
      do I=1,N_ELIMINATED_GAS
        IGAS_ELIM = IELIMINATED_IN_GAS(I)
        ! Will first check if it matches a beta type
        IMATCH_BETA = 0
        do J=1,NBETA
          if ((ITYPE_B == MAX_E_GAS_BETA(2,J)) .and. (IGAS_ELIM == MAX_E_GAS_BETA(1,J))) then
            IMATCH_BETA = 1
            exit
          end if
        end do
        ! Now check it also matches an alpha type
        IMATCH_ALPHA = 0
        do J=1,NALPHA
          if ((ITYPE_A == MAX_E_GAS_ALPHA(2,J)) .and. (IGAS_ELIM == MAX_E_GAS_ALPHA(1,J))) then
            IMATCH_ALPHA = 1
            exit
          end if
        end do
        if ((IMATCH_BETA == 1) .and. (IMATCH_ALPHA == 1)) IMATCH_BLOCK = 1
      end do
      if (IMATCH_BLOCK == 1) then
        N_ELIMINATED_BATCHES = N_ELIMINATED_BATCHES+1
        I_AM_OUT(N_ELIMINATED_BATCHES) = iblock
        cycle
      end if
    end if
    if (I_ELIMINATE_GAS > 1) then
      ! DEXS
      do I=1,N_2ELIMINATED_GAS
        IGAS_ELIM = I2ELIMINATED_IN_GAS(I)
        ! Will first check if it matches a beta type
        IMATCH_BETA = 0
        do J=1,NBETA
          if ((ITYPE_B == MAX_E_GAS_BETA(2,J)) .and. (IGAS_ELIM == MAX_E_GAS_BETA(1,J))) then
            IMATCH_BETA = 1
            exit
          end if
        end do
        ! Now check it also matches an alpha type
        IMATCH_ALPHA = 0
        do J=1,NALPHA
          if ((ITYPE_A == MAX_E_GAS_ALPHA(2,J)) .and. (IGAS_ELIM == MAX_E_GAS_ALPHA(1,J))) then
            IMATCH_ALPHA = 1
            exit
          end if
        end do
        ! Check if OCC matches a DEXS beta type
        IMATCH_BETAM1 = 0
        do J=1,NBETAM1
          if ((ITYPE_B == MAXM1_E_GAS_BETA(2,J)) .and. (IGAS_ELIM == MAXM1_E_GAS_BETA(1,J))) then
            IMATCH_BETAM1 = 1
            exit
          end if
        end do
        IMATCH_ALPHAM1 = 0
        do J=1,NALPHAM1
          if ((ITYPE_A == MAXM1_E_GAS_ALPHA(2,J)) .and. (IGAS_ELIM == MAXM1_E_GAS_ALPHA(1,J))) then
            IMATCH_ALPHAM1 = 1
            exit
          end if
        end do
        if (((IMATCH_BETA == 1) .and. (IMATCH_ALPHA == 1)) .or. ((IMATCH_BETAM1 == 1) .and. (IMATCH_ALPHA == 1)) .or. &
            ((IMATCH_BETA == 1) .and. (IMATCH_ALPHAM1 == 1))) IMATCH_BLOCK = 1
      end do
    end if
    if (IMATCH_BLOCK == 1) then
      N_ELIMINATED_BATCHES = N_ELIMINATED_BATCHES+1
      I_AM_OUT(N_ELIMINATED_BATCHES) = iblock
    end if
  end do
end do

if (N_ELIMINATED_BATCHES > MXPSTT) then
  write(u6,*) ' Increase MXPSTT to ',N_ELIMINATED_BATCHES
  call SYSABENDMSG('lucia_util/i_am_so_excited','Dimension of I_AM_OUT is too small','Increase MXPSTT')
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Number of eliminated blocks ',N_ELIMINATED_BATCHES
write(u6,*)
write(u6,*) ' The blocks eliminated'
write(u6,*)
do I=1,N_ELIMINATED_BATCHES
  write(u6,*) I_AM_OUT(I)
end do
write(u6,*)
#endif

end subroutine I_AM_SO_EXCITED
