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
! Copyright (C) 1989,2003, Jeppe Olsen                                 *
!***********************************************************************

subroutine DIHDJ_MOLCAS(IASTR,IBSTR,NIDET,JASTR,JBSTR,NJDET,NAEL,NBEL,IWORK,NORB,ONEBOD,HAMIL,ISYM,NINOB,ECORE,ICOMBI,PSIGN, &
                        IPRINT,TUVX,ExFac,IREOTS)
! A SET OF DETERMINANTS IA DEFINED BY ALPHA AND BETA STRINGS
! IASTR,IBSTR AND ANOTHER SET OF DETERMINATS DEFINED BY STRINGS
! JASTR AND JBSTR ARE GIVEN. OBTAIN CORRESPONDING HAMILTONIAN MATRIX
!
! IF ICOMBI /= 0 COMBINATIONS ARE USED FOR ALPHA AND BETA STRING
! THAT DIFFERS :
!   1/SQRT(2) * ( |I1A I2B| + PSIGN * |I2A I1B| )
!
! IF ISYM == 0 FULL HAMILTONIAN IS CONSTRUCTED
! IF ISYM /= 0 LOWER HALF OF HAMILTONIAN IS CONSTRUCTED
!
! JEPPE OLSEN JANUARY 1989
!             IREOTS added for new Molcas compatibility, August 2003

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NAEL, NBEL, NIDET, NJDET, JASTR(NAEL,NJDET), JBSTR(NBEL,NJDET), NORB, ISYM, NINOB, ICOMBI, &
                                 IPRINT, IREOTS(NORB)
integer(kind=iwp), intent(inout) :: IASTR(NAEL,NIDET), IBSTR(NBEL,NIDET)
integer(kind=iwp), intent(out) :: IWORK(4*NORB+NIDET)
real(kind=wp), intent(in) :: ONEBOD(NORB,NORB), ECORE, PSIGN, TUVX(*), ExFac
real(kind=wp), intent(_OUT_) :: HAMIL(*)
integer(kind=iwp) :: I1, I1_REO, I2, I2_REO, IA, IA_REO, IAB, IAEL, IAEQIB, IB, IB_REO, IBEL, IDET, IDIFF, IEL, IEL1, iii, ILOOP, &
                     IORB, IORB_REO, IPERM, J1, J1_REO, J2, J2_REO, JA, JA_REO, JAB, JAEL, JAEQJB, JB, JB_REO, JBEL, JDET, JDIFF, &
                     JEL, JEL1, jjj, JORB, JORB_REO, JPERM, KLFREE, KLIAB, KLIAE, KLIBE, KLJAE, KLJBE, LHAMIL, MINI, NACM, NADIF, &
                     NBCM, NBDIF, NDIF0, NDIF1, NDIF2, NIABEL, NJABEL, NLOOP, NTERMS, NTEST
real(kind=wp) :: CONST, SGN, SIGNA, SIGNB, XVAL
real(kind=wp), external :: GETH2A

NTEST = 0
! Initialization
KLIAB = 0
IAEQIB = 0
JAEQJB = 0
CONST = Zero
IEL1 = 0
JEL1 = 0
IPERM = 0
JPERM = 0
SGN = Zero
SIGNA = Zero
SIGNB = Zero
XVAL = Zero
IA = 0
IB = 0
JA = 0
JB = 0
I1 = 0
I2 = 0
J1 = 0
J2 = 0
! SCRATCH SPACE

! 1 : EXPANSION OF ALPHA AND BETA STRINGS OF TYPE I

KLFREE = 1
KLIAE = KLFREE
KLFREE = KLIAE+NORB
KLIBE = KLFREE
KLFREE = KLIBE+NORB

KLJAE = KLFREE
KLFREE = KLJAE+NORB
KLJBE = KLFREE
KLFREE = KLJBE+NORB
if (ICOMBI /= 0) then
  KLIAB = KLFREE
  KLFREE = KLFREE+NIDET
end if

if (ICOMBI /= 0) then
  ! SET UP ARRAY COMPARING ALPHA AND BETA STRINGS IN IDET LIST
  do IDET=1,NIDET
    IAEQIB = 1
    do IEL=1,NAEL
      if (IASTR(IEL,IDET) /= IBSTR(IEL,IDET)) IAEQIB = 0
    end do
    IWORK(KLIAB-1+IDET) = IAEQIB
  end do
end if

if (ISYM == 0) then
  LHAMIL = NIDET*NJDET
else
  LHAMIL = NIDET*(NIDET+1)/2
end if
HAMIL(1:LHAMIL) = Zero

NTERMS = 0
NDIF0 = 0
NDIF1 = 0
NDIF2 = 0
! LOOP OVER J DETERMINANTS

do JDET=1,NJDET

  ! EXPAND JDET
  call ICOPY(NORB,[0],0,IWORK(KLJAE),1)
  call ICOPY(NORB,[0],0,IWORK(KLJBE),1)

  if (ICOMBI /= 0) then
    JAEQJB = 1
    do IEL=1,NAEL
      if (JASTR(IEL,JDET) /= JBSTR(IEL,JDET)) JAEQJB = 0
    end do
  end if

  do IAEL=1,NAEL
    IWORK(KLJAE-1+JASTR(IAEL,JDET)) = 1
  end do

  do IBEL=1,NBEL
    IWORK(KLJBE-1+JBSTR(IBEL,JDET)) = 1
  end do

  if (NTEST >= 10) then
    write(u6,*) ' LOOP 1000 JDET =  ',JDET
    write(u6,*) ' JASTR AND JBSTR'
    call IWRTMA(JASTR(1,JDET),1,NAEL,1,NAEL)
    call IWRTMA(JBSTR(1,JDET),1,NBEL,1,NBEL)
    write(u6,*) ' EXPANDED ALPHA AND BETA STRING'
    call IWRTMA(IWORK(KLJAE),1,NORB,1,NORB)
    call IWRTMA(IWORK(KLJBE),1,NORB,1,NORB)
  end if

  if (ISYM == 0) then
    MINI = 1
  else
    MINI = JDET
  end if
  do IDET=MINI,NIDET
    if (ICOMBI == 0) then
      NLOOP = 1
    else
      IAEQIB = IWORK(KLIAB-1+IDET)
      if (IAEQIB+JAEQJB == 0) then
        NLOOP = 2
      else
        NLOOP = 1
      end if
    end if

    do ILOOP=1,NLOOP
      NTERMS = NTERMS+1

      ! COMPARE DETERMINANTS

      ! SWAP IA AND IB FOR SECOND PART OF COMBINATIONS
      if (ILOOP == 2) then
        do iii=1,NAEL
          jjj = IASTR(iii,IDET)
          IASTR(iii,IDET) = IBSTR(iii,IDET)
          IBSTR(iii,IDET) = jjj
        end do
      end if

      NACM = 0
      do IAEL=1,NAEL
        NACM = NACM+IWORK(KLJAE-1+IASTR(IAEL,IDET))
      end do

      NBCM = 0
      do IBEL=1,NBEL
        NBCM = NBCM+IWORK(KLJBE-1+IBSTR(IBEL,IDET))
      end do

      NADIF = NAEL-NACM
      NBDIF = NBEL-NBCM
      if (NTEST >= 10) then
        write(u6,*) '  LOOP 900 IDET ',IDET
        write(u6,*) ' COMPARISON, NADIF, NBDIF ',NADIF,NBDIF
      end if

      if (NADIF+NBDIF > 2) cycle

      ! FACTOR FOR COMBINATIONS
      if (ICOMBI == 0) then
        CONST = One
      else
        if ((JAEQJB+IAEQIB) == 2) then
          CONST = One
        else if ((JAEQJB+IAEQIB) == 1) then
          CONST = One/sqrt(Two)*(One+PSIGN)
        else if ((JAEQJB+IAEQIB) == 0) then
          if (ILOOP == 1) then
            CONST = One
          else
            CONST = PSIGN
          end if
        end if
      end if

      ! FIND DIFFERING ORBITALS AND SIGN FOR PERMUTATION

      ! EXPAND IDET
      call ICOPY(NORB,[0],0,IWORK(KLIAE),1)
      call ICOPY(NORB,[0],0,IWORK(KLIBE),1)

      do IAEL=1,NAEL
        IWORK(KLIAE-1+IASTR(IAEL,IDET)) = 1
      end do

      do IBEL=1,NBEL
        IWORK(KLIBE-1+IBSTR(IBEL,IDET)) = 1
      end do

      if (NADIF == 1) then
        do IAEL=1,NAEL
          if (IWORK(KLJAE-1+IASTR(IAEL,IDET)) == 0) then
            IA = IASTR(IAEL,IDET)
            IEL1 = IAEL
            exit
          end if
        end do

        do JAEL=1,NAEL
          if (IWORK(KLIAE-1+JASTR(JAEL,JDET)) == 0) then
            JA = JASTR(JAEL,JDET)
            JEL1 = JAEL
            exit
          end if
        end do
        SIGNA = real((-1)**(JEL1+IEL1),kind=wp)
      end if
      if (NBDIF == 1) then
        do IBEL=1,NBEL
          if (IWORK(KLJBE-1+IBSTR(IBEL,IDET)) == 0) then
            IB = IBSTR(IBEL,IDET)
            IEL1 = IBEL
            exit
          end if
        end do

        do JBEL=1,NBEL
          if (IWORK(KLIBE-1+JBSTR(JBEL,JDET)) == 0) then
            JB = JBSTR(JBEL,JDET)
            JEL1 = JBEL
            exit
          end if
        end do
        SIGNB = real((-1)**(JEL1+IEL1),kind=wp)
      end if
      if (NADIF == 2) then
        IDIFF = 0
        do IAEL=1,NAEL
          if (IWORK(KLJAE-1+IASTR(IAEL,IDET)) == 0) then
            if (IDIFF == 0) then
              IDIFF = 1
              I1 = IASTR(IAEL,IDET)
              IPERM = IAEL
            else
              I2 = IASTR(IAEL,IDET)
              IPERM = IAEL+IPERM
              exit
            end if
          end if
        end do

        JDIFF = 0
        do JAEL=1,NAEL
          if (IWORK(KLIAE-1+JASTR(JAEL,JDET)) == 0) then
            if (JDIFF == 0) then
              JDIFF = 1
              J1 = JASTR(JAEL,JDET)
              JPERM = JAEL
            else
              J2 = JASTR(JAEL,JDET)
              JPERM = JAEL+JPERM
              exit
            end if
          end if
        end do
        SGN = real((-1)**(IPERM+JPERM),kind=wp)
      end if

      if (NBDIF == 2) then
        IDIFF = 0
        do IBEL=1,NBEL
          if (IWORK(KLJBE-1+IBSTR(IBEL,IDET)) == 0) then
            if (IDIFF == 0) then
              IDIFF = 1
              I1 = IBSTR(IBEL,IDET)
              IPERM = IBEL
            else
              I2 = IBSTR(IBEL,IDET)
              IPERM = IBEL+IPERM
              exit
            end if
          end if
        end do

        JDIFF = 0
        do JBEL=1,NBEL
          if (IWORK(KLIBE-1+JBSTR(JBEL,JDET)) == 0) then
            if (JDIFF == 0) then
              JDIFF = 1
              J1 = JBSTR(JBEL,JDET)
              JPERM = JBEL
            else
              J2 = JBSTR(JBEL,JDET)
              JPERM = JBEL+JPERM
              exit
            end if
          end if
        end do
        SGN = real((-1)**(IPERM+JPERM),kind=wp)
      end if

      ! OBTAIN VALUE OF HAMILTONIAN ELEMENT

      if ((NADIF == 2) .or. (NBDIF == 2)) then
        NDIF2 = NDIF2+1
        ! SGN * (I1 J1 | I2 J2 ) - ( I1 J2 | I2 J1 )
        I1 = I1+NINOB
        I2 = I2+NINOB
        J1 = J1+NINOB
        J2 = J2+NINOB
        ! Well, there is no reordering in integrals any more. so
        I1_REO = IREOTS(I1)
        J1_REO = IREOTS(J1)
        I2_REO = IREOTS(I2)
        J2_REO = IREOTS(J2)
        XVAL = SGN*(GETH2A(I1_REO,J1_REO,I2_REO,J2_REO,TUVX)-GETH2A(I1_REO,J2_REO,I2_REO,J1_REO,TUVX))
      else if ((NADIF == 1) .and. (NBDIF == 1)) then
        NDIF2 = NDIF2+1
        ! SGN * (IA JA | IB JB )
        IA = IA+NINOB
        IB = IB+NINOB
        JA = JA+NINOB
        JB = JB+NINOB

        IA_REO = IREOTS(IA)
        IB_REO = IREOTS(IB)
        JA_REO = IREOTS(JA)
        JB_REO = IREOTS(JB)
        XVAL = SIGNA*SIGNB*GETH2A(IA_REO,JA_REO,IB_REO,JB_REO,TUVX)
      else if (((NADIF == 1) .and. (NBDIF == 0)) .or. ((NADIF == 0) .and. (NBDIF == 1))) then
        NDIF1 = NDIF1+1
        ! SGN * (  H(I1 J1 ) +
        !  (SUM OVER ORBITALS OF BOTH      SPIN TYPES  ( I1 J1 | JORB JORB )
        ! -(SUM OVER ORBITALS OF DIFFERING SPIN TYPE   ( I1 JORB | JORB J1 ) )
        if (NADIF == 1) then
          I1 = IA+NINOB
          J1 = JA+NINOB
          SGN = SIGNA
        else
          I1 = IB+NINOB
          J1 = JB+NINOB
          SGN = SIGNB
        end if

        I1_REO = IREOTS(I1)
        J1_REO = IREOTS(J1)
        XVAL = ONEBOD(IREOTS(I1-NINOB),IREOTS(J1-NINOB))
        do JAEL=1,NAEL
          JORB = JASTR(JAEL,JDET)+NINOB
          JORB_REO = IREOTS(JORB)
          XVAL = XVAL+GETH2A(I1_REO,J1_REO,JORB_REO,JORB_REO,TUVX)
        end do
        do JBEL=1,NBEL
          JORB = JBSTR(JBEL,JDET)+NINOB
          JORB_REO = IREOTS(JORB)
          XVAL = XVAL+GETH2A(I1_REO,J1_REO,JORB_REO,JORB_REO,TUVX)
        end do
        if (NADIF == 1) then
          do JAEL=1,NAEL
            JORB = JASTR(JAEL,JDET)+NINOB
            JORB_REO = IREOTS(JORB)
            XVAL = XVAL-GETH2A(I1_REO,JORB_REO,JORB_REO,J1_REO,TUVX)
          end do
        else
          do JBEL=1,NBEL
            JORB = JBSTR(JBEL,JDET)+NINOB
            JORB_REO = IREOTS(JORB)
            XVAL = XVAL-GETH2A(I1_REO,JORB_REO,JORB_REO,J1_REO,TUVX)
          end do
        end if
        XVAL = XVAL*SGN
      else if ((NADIF == 0) .and. (NBDIF == 0)) then
        NDIF0 = NDIF0+1
        ! SUM(I,J OF JDET) H(I,J) + (I I | J J ) - (I J | J I )

        XVAL = ECORE
        do IAB=1,2
          if (IAB == 1) then
            NIABEL = NAEL
          else
            NIABEL = NBEL
          end if
          do JAB=1,2
            if (JAB == 1) then
              NJABEL = NAEL
            else
              NJABEL = NBEL
            end if
            do IEL=1,NIABEL
              if (IAB == 1) then
                IORB = IASTR(IEL,IDET)
              else
                IORB = IBSTR(IEL,IDET)
              end if
              IORB_REO = IREOTS(IORB)
              if (IAB == JAB) XVAL = XVAL+ONEBOD(IORB_REO,IORB_REO)
              IORB = IORB+NINOB
              do JEL=1,NJABEL
                if (JAB == 1) then
                  JORB = IASTR(JEL,IDET)+NINOB
                else
                  JORB = IBSTR(JEL,IDET)+NINOB
                end if
                JORB_REO = IREOTS(JORB)
                XVAL = XVAL+Half*GETH2A(IORB_REO,IORB_REO,JORB_REO,JORB_REO,TUVX)
                if (IAB == JAB) XVAL = XVAL-ExFac*Half*GETH2A(IORB_REO,JORB_REO,JORB_REO,IORB_REO,TUVX)
              end do
            end do
          end do
        end do
      end if

      if (ISYM == 0) then
        HAMIL((JDET-1)*NIDET+IDET) = HAMIL((JDET-1)*NIDET+IDET)+CONST*XVAL
      else
        HAMIL((IDET-1)*IDET/2+JDET) = HAMIL((IDET-1)*IDET/2+JDET)+CONST*XVAL
      end if
      ! RESTORE ORDER
      if (ILOOP == 2) then
        do iii=1,NAEL
          jjj = IASTR(iii,IDET)
          IASTR(iii,IDET) = IBSTR(iii,IDET)
          IBSTR(iii,IDET) = jjj
        end do
      end if
    end do
  end do
end do

if (IPRINT >= 2) then
  write(u6,*) '  HAMILTONIAN MATRIX'
  if (ISYM == 0) then
    call WRTMAT(HAMIL,NIDET,NJDET,NIDET,NJDET)
  else
    call PRSYM(HAMIL,NIDET)
  end if
end if

return

end subroutine DIHDJ_MOLCAS
