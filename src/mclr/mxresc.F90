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

subroutine MXRESC(IAB,IACLS,IBCLS,NOCTPA,NOCTPB,NSMST,NSSOA,NSSOB,KCLS,KSSOA,KOCTPA,KSSOB,KOCTPB,NSMOB,MXPTOB,NTPOB,NTSOB,MXPKA, &
                  K2SSOA,K2OCTPA,K2SSOB,K2OCTPB,NAEL123,NBEL123,MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXIJST,MXIJSTF)
! Find largest dimension of matrix C(Ka,Ib,J)
! Find largest dimension of matrix C(ij,Ka,Ib)
! Find largest dimension of matrix C(ij,Ia,Kb)
! Find largest dimension of matrix C(ij,Ka,Kb)
! Find largest dimension of a+ia+j!k> MXIJST, MXIJSTF
!
! Largest block of single excitations MXSXBL
!
! To mess it up for the enemy: kcls
!
! Input:
!    iab
!    iacls,ibcls
!    noctpa,noctpb
!    nsmst
!    nssoa,nssob
!    kssoa
!    koctps,koctpb
!    nsmob
!    mxtpob,ntpob,ntsob
!    mxpka
!    k2ssoa,k2ssob
!    k2octpb
!    nael123,nbel123
! Output
!    mxcj,mxcija,mxcijb,mxcijab,mxsxbl,MXIJST,MXIJSTF

implicit real*8(A-H,O-Z)
dimension IAB(NOCTPA,NOCTPB)
dimension NSSOA(NOCTPA,NSMST), NSSOB(NOCTPB,NSMST)
dimension KSSOA(KOCTPA,NSMST), KSSOB(KOCTPB,NSMST)
dimension K2SSOA(K2OCTPA,NSMST)
dimension K2SSOB(K2OCTPB,NSMST)
dimension NTSOB(MXPTOB,NSMOB)
dimension NAEL123(NTPOB,*), NBEL123(NTPOB,*)

! matrix C(j,Ka,Ib)

MXCJ = 0
do IATP=1,NOCTPA
  do IBTP=1,NOCTPB
    if (IAB(IATP,IBTP) /= 0) then
      MXB = 0
      do ISM=1,NSMST
        MXB = max(MXB,NSSOB(IBTP,ISM))
      end do
      do IOBTP=1,NTPOB
        ! type of K string obtained by removing one elec of type IOPBTP from IATP
        call NEWTYP_MCLR(IACLS,IATP,[1],[IOBTP],1,KACLS,KATP)
        if (KATP > 0) then
          MXKA = 0
          do KSM=1,NSMST
            MXKA = max(MXKA,KSSOA(KATP,KSM))
          end do
          if ((MXPKA > 0) .and. (MXKA > MXPKA)) MXKA = MXPKA
          MXSOB = 0
          do ISMOB=1,NSMOB
            MXSOB = max(MXSOB,NTSOB(IOBTP,ISMOB))
          end do

          LCJBLK = MXSOB*MXKA*MXB
          MXCJ = max(MXCJ,LCJBLK)

        end if
      end do
    end if
  end do
end do

! matrix C(j,Ia,Kb)

do IATP=1,NOCTPA
  do IBTP=1,NOCTPB
    if (IAB(IATP,IBTP) /= 0) then
      MXA = 0
      do ISM=1,NSMST
        MXA = max(MXA,NSSOA(IATP,ISM))
      end do
      do IOBTP=1,NTPOB
        ! type of K string obtained by removing one elec of type IOPBTP from IBTP
        call NEWTYP_MCLR(IBCLS,IBTP,[1],[IOBTP],1,KBCLS,KBTP)
        if (KBTP > 0) then
          MXKB = 0
          do KSM=1,NSMST
            MXKB = max(MXKB,KSSOB(KBTP,KSM))
          end do
          if ((MXPKA > 0) .and. (MXKB > MXPKA)) MXKB = MXPKA
          MXSOB = 0
          do ISMOB=1,NSMOB
            MXSOB = max(MXSOB,NTSOB(IOBTP,ISMOB))
          end do

          LCJBLK = MXSOB*MXKB*MXA
          MXCJ = max(MXCJ,LCJBLK)

        end if
      end do
    end if
  end do
end do

! matrix C(ij,Ka,Ib)
! both Ka and Ib blocked

MXCIJA = 0
MXIJSTA = 0
MXIJSTAF = 0
do IATP=1,NOCTPA
  do IBTP=1,NOCTPB

    if (IAB(IATP,IBTP) /= 0) then
      MXIB = 0
      do ISM=1,NSMST
        MXIB = max(MXIB,NSSOB(IBTP,ISM))
      end do
      if ((MXIB > MXPKA) .and. (MXPKA > 0)) MXIB = MXPKA
      do IOBTP=1,NTPOB
        ! type of K string obtained by removing one elec of type IOPBTP from IATP
        call NEWTYP_MCLR(IACLS,IATP,[1],[IOBTP],1,K1ACLS,K1ATP)
        if (K1ATP > 0) then
          MXISOB = 0
          do ISMOB=1,NSMOB
            MXISOB = max(MXISOB,NTSOB(IOBTP,ISMOB))
          end do
          do JOBTP=1,NTPOB
            ! type of K string obtained by removing one elec of type JOPBTP from K1ATP
            call NEWTYP_MCLR(K1ACLS,K1ATP,[1],[JOBTP],1,KACLS,KATP)
            if (KATP > 0) then
              MXKA = 0
              do KSM=1,NSMST
                MXKA = max(MXKA,K2SSOA(KATP,KSM))
              end do
              MXKAF = MXKA
              if ((MXPKA > 0) .and. (MXKA > MXPKA)) MXKA = MXPKA
              MXJSOB = 0
              do JSMOB=1,NSMOB
                MXJSOB = max(MXJSOB,NTSOB(JOBTP,JSMOB))
              end do

              MXIJSTA = max(MXIJSTA,MXISOB*MXJSOB*MXKA)
              MXIJSTAF = max(MXIJSTAF,MXISOB*MXJSOB*MXKAF)
              LBLK = MXISOB*MXJSOB*MXKA*MXIB
              MXCIJA = max(MXCIJA,LBLK)
            end if
          end do
        end if
      end do
    end if
  end do
end do

! matrix C(ij,Ia,kb)
! both Ka and Ib blocked

MXCIJB = 0
MXIJSTB = 0
MXIJSTBF = 0
MXIA = 0 ! dummy initialize
do IATP=1,NOCTPA
  do IBTP=1,NOCTPB
    if (IAB(IATP,IBTP) /= 0) then
      MXIA = 0
      do ISM=1,NSMST
        MXIA = max(MXIA,NSSOA(IATP,ISM))
      end do
      if ((MXIA > MXPKA) .and. (MXPKA > 0)) MXIA = MXPKA
      do IOBTP=1,NTPOB
        ! type of K string obtained by removing one elec of type IOPBTP from IBTP
        call NEWTYP_MCLR(IBCLS,IBTP,[1],[IOBTP],1,K1BCLS,K1BTP)
        if (K1BTP > 0) then
          MXISOB = 0
          do ISMOB=1,NSMOB
            MXISOB = max(MXISOB,NTSOB(IOBTP,ISMOB))
          end do
          do JOBTP=1,NTPOB
            ! type of K string obtained by removing one elec of type JOPBTP from K1ATP
            call NEWTYP_MCLR(K1BCLS,K1BTP,[1],[JOBTP],1,KBCLS,KBTP)
            if (KBTP > 0) then
              MXKB = 0
              do KSM=1,NSMST
                MXKB = max(MXKB,K2SSOB(KBTP,KSM))
              end do
              MXKBF = MXKB
              if ((MXPKA > 0) .and. (MXKB > MXPKA)) MXKB = MXPKA
              MXJSOB = 0
              do JSMOB=1,NSMOB
                MXJSOB = max(MXJSOB,NTSOB(JOBTP,JSMOB))
              end do

              MXIJSTB = max(MXIJSTB,MXISOB*MXJSOB*MXKB)
              MXIJSTBF = max(MXIJSTBF,MXISOB*MXJSOB*MXKBF)
              LBLK = MXISOB*MXJSOB*MXKB*MXIA
              MXCIJB = max(MXCIJB,LBLK)
            end if
          end do
        end if
      end do
    end if
  end do
end do

! matrix C(ij,Ka,kb)
! both Ka and Kb blocked

MXCIJAB = 0
do IATP=1,NOCTPA
  do IBTP=1,NOCTPB

    if (IAB(IATP,IBTP) /= 0) then
      do IOBTP=1,NTPOB
        ! type of Ka string obtained by removing one elec of type IOPBTP from IATP
        call NEWTYP_MCLR(IACLS,IATP,[1],[IOBTP],1,KACLS,KATP)
        if (KATP > 0) then
          MXKA = 0
          do KSM=1,NSMST
            MXKA = max(MXKA,KSSOA(KATP,KSM))
          end do
          if ((MXPKA > 0) .and. (MXKA > MXPKA)) MXKA = MXPKA
          MXISOB = 0
          do ISMOB=1,NSMOB
            MXISOB = max(MXISOB,NTSOB(IOBTP,ISMOB))
          end do
          do JOBTP=1,NTPOB
            ! type of K string obtained by removing one elec of type JOPBTP from IBTP
            call NEWTYP_MCLR(IBCLS,IBTP,[1],[JOBTP],1,KBCLS,KBTP)
            if (KBTP > 0) then
              MXKB = 0
              do KSM=1,NSMST
                MXKB = max(MXKB,KSSOB(KBTP,KSM))
              end do
              if ((MXPKA > 0) .and. (MXKB > MXPKA)) MXKB = MXPKA
              MXJSOB = 0
              do JSMOB=1,NSMOB
                MXJSOB = max(MXJSOB,NTSOB(JOBTP,JSMOB))
              end do

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

MXSXBL = 0
! For alpha strings:
do IATP=1,NOCTPA
  MXIA = 0
  do ISM=1,NSMST
    MXIA = max(MXIA,NSSOA(IATP,ISM))
  end do
  ! Orbitals to be removed
  do JOBTP=1,NTPOB
    ! Is this removal allowed ??
    call NEWTYP_MCLR(IACLS,IATP,[1],[JOBTP],1,KACLS,KATP)
    if (KATP > 0) then
      ! Number of possible choices of J orbitals
      MXJOB = 0
      do JSMOB=1,NSMOB
        MXJOB = max(MXJOB,NTSOB(JOBTP,JSMOB))
      end do
      MXJOB = min(MXJOB,NAEL123(JOBTP,IATP))
      ! Then  : add an electron
      do IOBTP=1,NTPOB
        ! Allowed ?
        call NEWTYP_MCLR(KACLS,KATP,[2],[IOBTP],1,JACLS,JATP)
        if (JATP > 0) then
          MXIOB = 0
          do ISMOB=1,NSMOB
            MXIOB = max(MXIOB,NTSOB(IOBTP,ISMOB))
          end do

          MXSXBL = max(MXSXBL,MXIOB*MXJOB*MXIA)
        end if
      end do
    end if
  end do
end do

! For beta strings:
do IBTP=1,NOCTPB
  MXIB = 0
  do ISM=1,NSMST
    MXIB = max(MXIB,NSSOB(IBTP,ISM))
  end do
  ! Orbitals to be removed
  do JOBTP=1,NTPOB
    ! Is this removal allowed ??
    call NEWTYP_MCLR(IBCLS,IBTP,[1],[JOBTP],1,KBCLS,KBTP)
    if (KBTP > 0) then
      ! Number of possible choices of J orbitals
      MXJOB = 0
      do JSMOB=1,NSMOB
        MXJOB = max(MXJOB,NTSOB(JOBTP,JSMOB))
      end do
      MXJOB = min(MXJOB,NBEL123(JOBTP,IBTP))
      ! Then  : add an electron
      do IOBTP=1,NTPOB
        ! Allowed ?
        call NEWTYP_MCLR(KBCLS,KBTP,[2],[IOBTP],1,JBCLS,JBTP)
        if (JATP > 0) then
          MXIOB = 0
          do ISMOB=1,NSMOB
            MXIOB = max(MXIOB,NTSOB(IOBTP,ISMOB))
          end do

          MXSXBL = max(MXSXBL,MXIOB*MXJOB*MXIA)
        end if
      end do
    end if
  end do
end do

MXIJST = max(MXIJSTA,MXIJSTB)
MXIJSTF = max(MXIJSTAF,MXIJSTBF)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(KCLS)

end subroutine MXRESC
