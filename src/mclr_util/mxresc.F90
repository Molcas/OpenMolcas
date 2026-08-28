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

subroutine MXRESC(IAB,IACLS,IBCLS,NOCTPA,NOCTPB,NSM,NSSOA,NSSOB,KSSOA,KOCTPA,KSSOB,KOCTPB,NTSOB,MXPKA,K2SSOA,K2OCTPA,K2SSOB, &
                  K2OCTPB,NAEL123,NBEL123,MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXIJST,MXIJSTF)
! Find largest dimension of matrix C(Ka,Ib,J)
! Find largest dimension of matrix C(ij,Ka,Ib)
! Find largest dimension of matrix C(ij,Ia,Kb)
! Find largest dimension of matrix C(ij,Ka,Kb)
! Find largest dimension of a+ia+j!k> MXIJST, MXIJSTF
!
! Largest block of single excitations MXSXBL
!
! Input:
!    iab
!    iacls,ibcls
!    noctpa,noctpb
!    nsm
!    nssoa,nssob
!    kssoa
!    koctps,koctpb
!    ntsob
!    mxpka
!    k2ssoa,k2ssob
!    k2octpb
!    nael123,nbel123
! Output
!    mxcj,mxcija,mxcijb,mxcijab,mxsxbl,MXIJST,MXIJSTF

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NOCTPA, NOCTPB, IAB(NOCTPA,NOCTPB), IACLS, IBCLS, NSM, NSSOA(NOCTPA,NSM), NSSOB(NOCTPB,NSM), &
                                 KOCTPA, KSSOA(KOCTPA,NSM), KOCTPB, KSSOB(KOCTPB,NSM), NTSOB(3,NSM), MXPKA, K2OCTPA, &
                                 K2SSOA(K2OCTPA,NSM), K2OCTPB, K2SSOB(K2OCTPB,NSM), NAEL123(3,*), NBEL123(3,*)
integer(kind=iwp), intent(out) :: MXCJ, MXCIJA, MXCIJB, MXCIJAB, MXSXBL, MXIJST, MXIJSTF
integer(kind=iwp) :: IATP, IBTP, IOBTP, JACLS, JATP, JBCLS, JBTP, JOBTP, K1ACLS, K1ATP, K1BCLS, K1BTP, KACLS, KATP, KBCLS, KBTP, &
                     LBLK, LCJBLK, MXA, MXB, MXIA, MXIB, MXIJSTA, MXIJSTAF, MXIJSTB, MXIJSTBF, MXIOB, MXISOB, MXJOB, MXJSOB, MXKA, &
                     MXKAF, MXKB, MXKBF, MXSOB

! matrix C(j,Ka,Ib)

MXCJ = 0
do IATP=1,NOCTPA
  do IBTP=1,NOCTPB
    if (IAB(IATP,IBTP) /= 0) then
      MXB = max(0,maxval(NSSOB(IBTP,:)))
      do IOBTP=1,3
        ! type of K string obtained by removing one elec of type IOPBTP from IATP
        call NEWTYP_MCLR(IACLS,IATP,[1],[IOBTP],1,KACLS,KATP)
        if (KATP > 0) then
          MXKA = max(0,maxval(KSSOA(KATP,:)))
          if ((MXPKA > 0) .and. (MXKA > MXPKA)) MXKA = MXPKA
          MXSOB = max(0,maxval(NTSOB(IOBTP,:)))

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
      MXA = max(0,maxval(NSSOA(IATP,:)))
      do IOBTP=1,3
        ! type of K string obtained by removing one elec of type IOPBTP from IBTP
        call NEWTYP_MCLR(IBCLS,IBTP,[1],[IOBTP],1,KBCLS,KBTP)
        if (KBTP > 0) then
          MXKB = max(0,maxval(KSSOB(KBTP,:)))
          if ((MXPKA > 0) .and. (MXKB > MXPKA)) MXKB = MXPKA
          MXSOB = max(0,maxval(NTSOB(IOBTP,:)))

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
      MXIB = max(0,maxval(NSSOB(IBTP,:)))
      if ((MXIB > MXPKA) .and. (MXPKA > 0)) MXIB = MXPKA
      do IOBTP=1,3
        ! type of K string obtained by removing one elec of type IOPBTP from IATP
        call NEWTYP_MCLR(IACLS,IATP,[1],[IOBTP],1,K1ACLS,K1ATP)
        if (K1ATP > 0) then
          MXISOB = max(0,maxval(NTSOB(IOBTP,:)))
          do JOBTP=1,3
            ! type of K string obtained by removing one elec of type JOPBTP from K1ATP
            call NEWTYP_MCLR(K1ACLS,K1ATP,[1],[JOBTP],1,KACLS,KATP)
            if (KATP > 0) then
              MXKA = max(0,maxval(K2SSOA(KATP,:)))
              MXKAF = MXKA
              if ((MXPKA > 0) .and. (MXKA > MXPKA)) MXKA = MXPKA
              MXJSOB = max(0,maxval(NTSOB(JOBTP,:)))

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
      MXIA = max(0,maxval(NSSOA(IATP,:)))
      if ((MXIA > MXPKA) .and. (MXPKA > 0)) MXIA = MXPKA
      do IOBTP=1,3
        ! type of K string obtained by removing one elec of type IOPBTP from IBTP
        call NEWTYP_MCLR(IBCLS,IBTP,[1],[IOBTP],1,K1BCLS,K1BTP)
        if (K1BTP > 0) then
          MXISOB = max(0,maxval(NTSOB(IOBTP,:)))
          do JOBTP=1,3
            ! type of K string obtained by removing one elec of type JOPBTP from K1ATP
            call NEWTYP_MCLR(K1BCLS,K1BTP,[1],[JOBTP],1,KBCLS,KBTP)
            if (KBTP > 0) then
              MXKB = max(0,maxval(K2SSOB(KBTP,:)))
              MXKBF = MXKB
              if ((MXPKA > 0) .and. (MXKB > MXPKA)) MXKB = MXPKA
              MXJSOB = max(0,maxval(NTSOB(JOBTP,:)))

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
      do IOBTP=1,3
        ! type of Ka string obtained by removing one elec of type IOPBTP from IATP
        call NEWTYP_MCLR(IACLS,IATP,[1],[IOBTP],1,KACLS,KATP)
        if (KATP > 0) then
          MXKA = max(0,maxval(KSSOA(KATP,:)))
          if ((MXPKA > 0) .and. (MXKA > MXPKA)) MXKA = MXPKA
          MXISOB = max(0,maxval(NTSOB(IOBTP,:)))
          do JOBTP=1,3
            ! type of K string obtained by removing one elec of type JOPBTP from IBTP
            call NEWTYP_MCLR(IBCLS,IBTP,[1],[JOBTP],1,KBCLS,KBTP)
            if (KBTP > 0) then
              MXKB = max(0,maxval(KSSOB(KBTP,:)))
              if ((MXPKA > 0) .and. (MXKB > MXPKA)) MXKB = MXPKA
              MXJSOB = max(0,maxval(NTSOB(JOBTP,:)))

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
  MXIA = max(0,maxval(NSSOA(IATP,:)))
  ! Orbitals to be removed
  do JOBTP=1,3
    ! Is this removal allowed ??
    call NEWTYP_MCLR(IACLS,IATP,[1],[JOBTP],1,KACLS,KATP)
    if (KATP > 0) then
      ! Number of possible choices of J orbitals
      MXJOB = max(0,maxval(NTSOB(JOBTP,:)))
      MXJOB = min(MXJOB,NAEL123(JOBTP,IATP))
      ! Then  : add an electron
      do IOBTP=1,3
        ! Allowed ?
        call NEWTYP_MCLR(KACLS,KATP,[2],[IOBTP],1,JACLS,JATP)
        if (JATP > 0) then
          MXIOB = max(0,maxval(NTSOB(IOBTP,:)))

          MXSXBL = max(MXSXBL,MXIOB*MXJOB*MXIA)
        end if
      end do
    end if
  end do
end do

! For beta strings:
do IBTP=1,NOCTPB
  MXIB = max(0,maxval(NSSOB(IBTP,:)))
  ! Orbitals to be removed
  do JOBTP=1,3
    ! Is this removal allowed ??
    call NEWTYP_MCLR(IBCLS,IBTP,[1],[JOBTP],1,KBCLS,KBTP)
    if (KBTP > 0) then
      ! Number of possible choices of J orbitals
      MXJOB = max(0,maxval(NTSOB(JOBTP,:)))
      MXJOB = min(MXJOB,NBEL123(JOBTP,IBTP))
      ! Then  : add an electron
      do IOBTP=1,3
        ! Allowed ?
        call NEWTYP_MCLR(KBCLS,KBTP,[2],[IOBTP],1,JBCLS,JBTP)
        if (JATP > 0) then
          MXIOB = max(0,maxval(NTSOB(IOBTP,:)))

          MXSXBL = max(MXSXBL,MXIOB*MXJOB*MXIA)
        end if
      end do
    end if
  end do
end do

MXIJST = max(MXIJSTA,MXIJSTB)
MXIJSTF = max(MXIJSTAF,MXIJSTBF)

end subroutine MXRESC
