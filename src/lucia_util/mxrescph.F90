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

!#define _DEBUGPRINT_
subroutine MXRESCPH(IAB,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NSMST,NSTFSMSPGP,MXPNSMST,NSMOB,MXPTOB,NTPOB,NTSOB,MXPKA,NEL1234,MXCJ,MXCIJA, &
                    MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
! Find largest dimension of matrix C(Ka,Ib,J)
! Find largest dimension of matrix C(ij,Ka,Ib)
! Find largest dimension of matrix C(ij,Ia,Kb)
! Find largest dimension of matrix C(ij,Ka,Kb)
! Find largest dimension of matrix S(P,Ia,I,K) for a single K-string
!
! Particle hole version : hole electrons added, particle elec removed
!
! Largest block of single excitations MXSXBL
!
! Input
! IAB    : allowed combination of alpha and beta supergroups
! IOCPTA : Number of first active alpha supergroup
! IOCPTB : Number of first active beta  supergroup
! NOCTPA : Number of active alpha supergroups
! NOCTPB : Number of active alpha supergroups
!
! Version of Jan 98 : IPHGAS added

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NOCTPA, NOCTPB, IAB(NOCTPA,NOCTPB), IOCTPA, IOCTPB, NSMST, MXPNSMST, NSTFSMSPGP(MXPNSMST,*), &
                                 NSMOB, MXPTOB, NTPOB, NTSOB(MXPTOB,NSMOB), MXPKA, NEL1234(MXPTOB,*), IPHGAS(*), NHLFSPGP(*), &
                                 MNHL, IADVICE
integer(kind=iwp), intent(out) :: MXCJ, MXCIJA, MXCIJB, MXCIJAB, MXSXBL, MXADKBLK, MXCJ_ALLSYM, MXADKBLK_AS, MX_NSPII
integer(kind=iwp) :: IAORC, IATP, IATPABS, IBTP, IBTPABS, IOBTP, ISM, ISMOB, ITOTA, ITOTB, JAORC, JATP, JBTP, JOBTP, JSMOB, K1ATP, &
                     K1BTP, KATP, KBTP, KSM, LBLK, LCJBLK, LCJBLK_ALLSYM, MXA, MXB, MXIA, MXIB, MXIOB, MXISOB, MXJOB, MXJSOB, &
                     MXKA, MXKAB, MXKACTEL, MXKAO, MXKB, MXKBO, MXSOB, NSOB_AS
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IATP_MX, IBTP_MX, IOBTP_MX, KATP_MX

write(u6,*) ' MXRESCPH : MXPKA ',MXPKA
#endif

! matrix C(j,Ka,Ib)

MXKAB = 0
MXCJ = 0
MXCJ_ALLSYM = 0
MXADKBLK = 0
MXADKBLK_AS = 0
MX_NSPII = 0
do IAORC=1,2
  do IATP=1,NOCTPA
    IATPABS = IATP+IOCTPA-1
    do IBTP=1,NOCTPB
      IBTPABS = IBTP+IOCTPB-1

      if (IAB(IATP,IBTP) /= 0) then
#       ifdef _DEBUGPRINT_
        write(u6,*) ' allowed IATP,IBTP',IATP,IBTP
#       endif
        MXB = 0
        ITOTB = 0
        do ISM=1,NSMST
          MXB = max(MXB,NSTFSMSPGP(ISM,IBTPABS))
          ITOTB = ITOTB+NSTFSMSPGP(ISM,IBTPABS)
        end do
#       ifdef _DEBUGPRINT_
        write(u6,*) ' MXB,ITOTB = ',MXB,ITOTB
#       endif
        do IOBTP=1,NTPOB
          ! No K strings obtained from creation in particle space
          if ((IAORC == 2) .and. (IPHGAS(IOBTP) == 1)) cycle
          ! type of K string obtained
          call NEWTYP(IATPABS,IAORC,IOBTP,KATP)
#         ifdef _DEBUGPRINT_
          write(u6,*) ' IOBTP KATP ',IOBTP,KATP
#         endif
          ! addi constraint to avoid calc with long columns and few rows
          ! Works only in connection with active advice routine !
          if (KATP > 0) then
            if ((IAORC == 1) .and. (IADVICE == 1) .and. (NHLFSPGP(IBTPABS)+NHLFSPGP(KATP) < MNHL) .and. &
                (NHLFSPGP(IATPABS) > (NHLFSPGP(IBTPABS)+1))) then
              !write(u6,*) ' N-1 hole space eliminated'
              !write(u6,*) ' IOBTP,IBTPABS,KATP',IOBTP,IBTPABS,KATP
              KATP = 0
            end if
          end if

          if (KATP > 0) then
            !     DIM_SPII(IASPGRP,IBSPGRP,IOBTP,IAB,IAC,NSPII)
            !call DIM_SPII(IATPABS,IBTPABS,IOBTP,1,IAORC,NSPII)
            !MX_NSPII = MAX(MX_NSPII,NSPII)
            MX_NSPII = 0
          end if

          if (KATP > 0) then
            MXKA = max(0,maxval(NSTFSMSPGP(1:NSMST,KATP)))
#           ifdef _DEBUGPRINT_
            write(u6,*) ' MXKA = ',MXKA
#           endif
            MXKAO = MXKA
            if ((MXPKA > 0) .and. (MXKA > MXPKA)) MXKA = MXPKA
            MXSOB = 0
            NSOB_AS = 0
            do ISMOB=1,NSMOB
              MXSOB = max(MXSOB,NTSOB(IOBTP,ISMOB))
              NSOB_AS = NSOB_AS+NTSOB(IOBTP,ISMOB)
            end do
#           ifdef _DEBUGPRINT_
            write(u6,*) ' MXSOB = ',MXSOB
#           endif

            MXADKBLK = max(MXADKBLK,MXSOB*MXKAO)
            MXADKBLK_AS = max(MXADKBLK_AS,NSOB_AS*MXKAO)
            LCJBLK = MXSOB*MXKA*MXB
            LCJBLK_ALLSYM = NSOB_AS*MXKA*ITOTB

            if (LCJBLK > MXCJ) then
              MXCJ = LCJBLK
#             ifdef _DEBUGPRINT_
              IATP_MX = IATP
              IBTP_MX = IBTP
              KATP_MX = KATP
              IOBTP_MX = IOBTP
#             endif
            end if
            MXCJ_ALLSYM = max(MXCJ_ALLSYM,LCJBLK_ALLSYM)
            MXKAB = max(MXKAB,MXKA)

          end if
        end do
      end if
    end do
  end do
end do
! End of anni/crea map

! matrix C(j,Ia,Kb)

do IAORC=1,2
  do IATP=1,NOCTPA
    IATPABS = IATP+IOCTPA-1
    do IBTP=1,NOCTPB
      IBTPABS = IBTP+IOCTPB-1

      if (IAB(IATP,IBTP) /= 0) then
#       ifdef _DEBUGPRINT_
        write(u6,*) ' allowed IATP,IBTP',IATP,IBTP
#       endif
        MXA = 0
        ITOTA = 0
        do ISM=1,NSMST
          MXA = max(MXA,NSTFSMSPGP(ISM,IATPABS))
          ITOTA = ITOTA+NSTFSMSPGP(ISM,IATPABS)
        end do
#       ifdef _DEBUGPRINT_
        write(u6,*) ' MXA = ',MXA
#       endif
        do IOBTP=1,NTPOB
          ! type of K string obtained by removing one elec of type IOPBTP from IATP
          if ((IAORC == 2) .and. (IPHGAS(IOBTP) == 1)) cycle
          call NEWTYP(IBTPABS,IAORC,IOBTP,KBTP)
#         ifdef _DEBUGPRINT_
          write(u6,*) ' IOBTP KBTP ',IOBTP,KBTP
#         endif
          if (KBTP > 0) then
            if ((IAORC == 1) .and. (IADVICE == 1) .and. (NHLFSPGP(IATPABS)+NHLFSPGP(KBTP) < MNHL) .and. &
                (NHLFSPGP(IBTPABS) > NHLFSPGP(IATPABS)+1)) then
              !write(u6,*) ' N-1 hole space eliminated'
              !write(u6,*) ' IOBTP,IATPABS,KBTP',IOBTP,IATPABS,KBTP
              KBTP = 0
            end if
          end if
          if (KBTP > 0) then
            !call DIM_SPII(IATPABS,IBTPABS,IOBTP,2,IAORC,NSPII)
            !MX_NSPII = MAX(MX_NSPII,NSPII)
            MX_NSPII = 0
          end if

          if (KBTP > 0) then
            MXKB = 0
            do KSM=1,NSMST
              MXKB = max(MXKB,NSTFSMSPGP(KSM,KBTP))
            end do
#           ifdef _DEBUGPRINT_
            write(u6,*) ' MXKB = ',MXKB
#           endif
            MXKBO = MXKB
            if ((MXPKA > 0) .and. (MXKB > MXPKA)) MXKB = MXPKA
            MXSOB = 0
            NSOB_AS = 0
            do ISMOB=1,NSMOB
              MXSOB = max(MXSOB,NTSOB(IOBTP,ISMOB))
              NSOB_AS = NSOB_AS+NTSOB(IOBTP,ISMOB)
            end do
#           ifdef _DEBUGPRINT_
            write(u6,*) ' MXSOB = ',MXSOB
#           endif

            MXADKBLK = max(MXADKBLK,MXSOB*MXKBO)
            MXADKBLK_AS = max(MXADKBLK_AS,NSOB_AS*MXKBO)
            !JULY29 LCJBLK = MXSOB*MXKB*MXB
            LCJBLK = MXSOB*MXKB*MXA
            LCJBLK_ALLSYM = NSOB_AS*MXKB*ITOTA
            MXCJ = max(MXCJ,LCJBLK)
            MXCJ_ALLSYM = max(MXCJ_ALLSYM,LCJBLK_ALLSYM)
            MXKAB = max(MXKAB,MXKB)

          end if
        end do
      end if
    end do
  end do
end do
! End of loop over creation/annihilation
#ifdef _DEBUGPRINT_
write(u6,*) 'MXRESCPH : MXADKBLK,MXCJ ',MXADKBLK,MXCJ
write(u6,*) ' MXCJ_ALLSYM = ',MXCJ_ALLSYM
#endif

! matrix C(ij,Ka,Ib)
! both Ka and Ib blocked

MXCIJA = 0
do IATP=1,NOCTPA
  IATPABS = IATP+IOCTPA-1
  do IBTP=1,NOCTPB
    IBTPABS = IBTP+IOCTPB-1

    if (IAB(IATP,IBTP) /= 0) then
      MXIB = 0
      do ISM=1,NSMST
        MXIB = max(MXIB,NSTFSMSPGP(ISM,IBTPABS))
      end do
      if (MXIB > MXPKA) MXIB = MXPKA
#     ifdef _DEBUGPRINT_
      write(u6,*) ' MXIB = ',MXIB
#     endif
      do IAORC=1,2
        do IOBTP=1,NTPOB
          ! type of K string obtained by removing one elec of type IOPBTP from IATP
          call NEWTYP(IATPABS,IAORC,IOBTP,K1ATP)
          ! No N+1 mappings for particle spaces
          if ((IAORC == 2) .and. (IPHGAS(IOBTP) == 1)) K1ATP = 0
#         ifdef _DEBUGPRINT_
          write(u6,*) ' IOBTP K1ATP ',IOBTP,K1ATP
#         endif
          if (K1ATP > 0) then
            MXISOB = 0
            do ISMOB=1,NSMOB
              MXISOB = max(MXISOB,NTSOB(IOBTP,ISMOB))
            end do
#           ifdef _DEBUGPRINT_
            write(u6,*) ' MXISOB = ',MXISOB
#           endif
            do JAORC=1,2
              do JOBTP=1,NTPOB
                ! type of K string obtained by removing one elec of type JOPBTP from K1ATP
                call NEWTYP(K1ATP,JAORC,JOBTP,KATP)
                if ((JAORC == 2) .and. (IPHGAS(JOBTP) == 1)) KATP = 0
                if (KATP > 0) then
                  MXKA = 0
                  do KSM=1,NSMST
                    MXKA = max(MXKA,NSTFSMSPGP(KSM,KATP))
                  end do
#                 ifdef _DEBUGPRINT_
                  write(u6,*) ' MXKA = ',MXKA
#                 endif
                  if ((MXPKA > 0) .and. (MXKA > MXPKA)) MXKA = MXPKA
                  MXJSOB = 0
                  do JSMOB=1,NSMOB
                    MXJSOB = max(MXJSOB,NTSOB(JOBTP,JSMOB))
                  end do
#                 ifdef _DEBUGPRINT_
                  write(u6,*) ' MXJSOB = ',MXJSOB
#                 endif

                  LBLK = MXISOB*MXJSOB*MXKA*MXIB
                  MXCIJA = max(MXCIJA,LBLK)
                end if
              end do
            end do
            ! End of loop over JOBTP, JAORC
          end if
        end do
      end do
      ! End of loop over IOBTP, IAORC
    end if
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) 'MXRESCPH : MXCIJA ',MXCIJA
#endif

! matrix C(ij,Ia,kb)
! both Ka and Ib blocked

MXCIJB = 0
do IATP=1,NOCTPA
  IATPABS = IATP+IOCTPA-1
  do IBTP=1,NOCTPB
    IBTPABS = IBTP+IOCTPB-1
    if (IAB(IATP,IBTP) /= 0) then
      MXIA = max(0,maxval(NSTFSMSPGP(1:NSMST,IATPABS)))
      if (MXIA > MXPKA) MXIA = MXPKA
#     ifdef _DEBUGPRINT_
      write(u6,*) ' MXIA = ',MXIA
#     endif
      do IAORC=1,2
        do IOBTP=1,NTPOB
          ! type of K string obtained by removing one elec of type IOPBTP from IBTP
          call NEWTYP(IBTPABS,IAORC,IOBTP,K1BTP)
#         ifdef _DEBUGPRINT_
          write(u6,*) ' IOBTP K1BTP ',IOBTP,K1BTP
#         endif
          if ((IAORC == 2) .and. (IPHGAS(IOBTP) == 1)) K1BTP = 0
          if (K1BTP > 0) then
            MXISOB = 0
            do ISMOB=1,NSMOB
              MXISOB = max(MXISOB,NTSOB(IOBTP,ISMOB))
            end do
#           ifdef _DEBUGPRINT_
            write(u6,*) ' MXISOB = ',MXISOB
#           endif
            do JAORC=1,2
              do JOBTP=1,NTPOB
                ! type of K string obtained by removing one elec of type JOPBTP from K1ATP
                call NEWTYP(K1BTP,JAORC,JOBTP,KBTP)
                if ((JAORC == 2) .and. (IPHGAS(JOBTP) == 1)) KBTP = 0
                if (KBTP > 0) then
                  MXKB = 0
                  do KSM=1,NSMST
                    MXKB = max(MXKB,NSTFSMSPGP(KSM,KBTP))
                  end do
#                 ifdef _DEBUGPRINT_
                  write(u6,*) ' MXKB = ',MXKB
#                 endif
                  if ((MXPKA > 0) .and. (MXKB > MXPKA)) MXKB = MXPKA
                  MXJSOB = 0
                  do JSMOB=1,NSMOB
                    MXJSOB = max(MXJSOB,NTSOB(JOBTP,JSMOB))
                  end do
#                 ifdef _DEBUGPRINT_
                  write(u6,*) ' MXJSOB = ',MXJSOB
#                 endif

                  MXCIJB = max(MXCIJB,MXISOB*MXJSOB*MXKB*MXIA)
                end if
              end do
            end do
            ! End of loop over JOBTP,JAORC
          end if
        end do
      end do
      ! End of loop over IOBTP,IAORC
    end if
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) 'MXRESCPH : MXCIJB ',MXCIJB
#endif

! matrix C(ij,Ka,kb)
! both Ka and Kb blocked

! Modified : Only used if at most two elecs in i and j
!            No batching
!            Used for hardwired few electron code
MXCIJAB = 0
MXKACTEL = 1
do IATP=1,NOCTPA
  IATPABS = IATP+IOCTPA-1
  do IBTP=1,NOCTPB
    IBTPABS = IBTP+IOCTPB-1
    if (IAB(IATP,IBTP) /= 0) then
      do IOBTP=1,NTPOB
        ! type of Ka string obtained by removing one elec of type IOPBTP from IATP
        call NEWTYP(IATPABS,1,IOBTP,KATP)
#       ifdef _DEBUGPRINT_
        write(u6,*) ' IOBTP KATP ',IOBTP,KATP
#       endif
        if (KATP > 0) then
          !   NEL1234(JOBTP,IATPABS)
          if (NEL1234(IOBTP,KATP) > MXKACTEL) KATP = 0
        end if
        if (KATP > 0) then
          MXKA = 0
          do KSM=1,NSMST
            MXKA = max(MXKA,NSTFSMSPGP(KSM,KATP))
          end do
#         ifdef _DEBUGPRINT_
          write(u6,*) ' MXKA = ',MXKA
#         endif
          ! No partitioning
          !if ((MXPKA > 0) .and. (MXKA > MXPKA)) MXKA = MXPKA

          MXISOB = 0
          do ISMOB=1,NSMOB
            MXISOB = max(MXISOB,NTSOB(IOBTP,ISMOB))
          end do
#         ifdef _DEBUGPRINT_
          write(u6,*) ' MXISOB = ',MXISOB
#         endif
          do JOBTP=1,NTPOB
            ! type of K string obtained by removing one elec of type JOPBTP from IBTP
            call NEWTYP(IBTPABS,1,JOBTP,KBTP)
            if (KBTP > 0) then
              if (NEL1234(JOBTP,KBTP) > MXKACTEL) KBTP = 0
            end if
            if (KBTP > 0) then
              MXKB = 0
              do KSM=1,NSMST
                MXKB = max(MXKB,NSTFSMSPGP(KSM,KBTP))
              end do
#             ifdef _DEBUGPRINT_
              write(u6,*) ' MXKB = ',MXKB
#             endif
              ! No partitioning
              !if ((MXPKA > 0) .and. (MXKB > MXPKA)) MXKB = MXPKA
              MXJSOB = 0
              do JSMOB=1,NSMOB
                MXJSOB = max(MXJSOB,NTSOB(JOBTP,JSMOB))
              end do
#             ifdef _DEBUGPRINT_
              write(u6,*) ' MXJSOB = ',MXJSOB
#             endif

              LBLK = MXISOB*MXJSOB*MXKB*MXKA
              MXCIJAB = max(MXCIJAB,LBLK)
            end if
          end do
        end if
      end do
    end if
  end do
end do

! Largest block of single excitations :
! Strings of given type and sym, orbitals of given type and sym

! Largest block of creations : a+i !kstring> where K string is
! obtained as single annihilations
MXSXBL = 0
! For alpha strings :
do IATP=1,NOCTPA
  IATPABS = IATP+IOCTPA-1
  MXIA = max(0,maxval(NSTFSMSPGP(1:NSMST,IATPABS)))
# ifdef _DEBUGPRINT_
  write(u6,*) ' MXIA = ',MXIA
# endif
  ! Orbitals to be removed
  do JOBTP=1,NTPOB
    ! Is this removal allowed ??
    call NEWTYP(IATPABS,1,JOBTP,KATP)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' JOBTP KATP ',JOBTP,KATP
#   endif
    if (KATP > 0) then
      ! Number of possible choices of J orbitals
      MXJOB = max(0,maxval(NTSOB(JOBTP,:)))
      MXJOB = min(MXJOB,NEL1234(JOBTP,IATPABS))
#     ifdef _DEBUGPRINT_
      write(u6,*) ' MXJOB = ',MXJOB
#     endif
      ! Then  : add an electron
      do IOBTP=1,NTPOB
        ! Allowed ?
        call NEWTYP(KATP,2,IOBTP,JATP)
        if (JATP > 0) then
          MXIOB = max(0,maxval(NTSOB(IOBTP,:)))
          MXSXBL = max(MXSXBL,MXIOB*MXJOB*MXIA)
        end if
      end do
    end if
  end do
end do

! For beta  strings :
do IBTP=1,NOCTPB
  IBTPABS = IBTP+IOCTPB-1
  MXIB = max(0,maxval(NSTFSMSPGP(1:NSMST,IBTPABS)))
# ifdef _DEBUGPRINT_
  write(u6,*) ' MXIB = ',MXIB
# endif
  ! Orbitals to be removed
  do JOBTP=1,NTPOB
    ! Is this removal allowed ??
    call NEWTYP(IBTPABS,1,JOBTP,KBTP)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' JOBTP KBTP ',JOBTP,KBTP
#   endif
    if (KBTP > 0) then
      ! Number of possible choices of J orbitals
      MXJOB = max(0,maxval(NTSOB(JOBTP,:)))
      MXJOB = min(MXJOB,NEL1234(JOBTP,IBTP))
#     ifdef _DEBUGPRINT_
      write(u6,*) ' MXJOB = ',MXJOB
#     endif
      ! Then  : add an electron
      do IOBTP=1,NTPOB
        ! Allowed ?
        call NEWTYP(KBTP,2,IOBTP,JBTP)
        if (JATP > 0) then
          MXIOB = max(0,maxval(NTSOB(IOBTP,:)))
          MXSXBL = max(MXSXBL,MXIOB*MXJOB*MXIB)
        end if
      end do
    end if
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) 'MXRESCPH: MXSXBL : ',MXSXBL
write(u6,*) ' MXRESCPH : MXKAB = ',MXKAB
write(u6,*) ' Info on largest C(Ka,j,Jb) block'
write(u6,*) ' IATP_MX, IBTP_MX, KATP_MX, IOBTP_MX ',IATP_MX,IBTP_MX,KATP_MX,IOBTP_MX
write(u6,*) ' MX_NSPII = ',MX_NSPII
#endif

end subroutine MXRESCPH
