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
! Copyright (C) 1989,1993, Jeppe Olsen                                 *
!***********************************************************************

subroutine DIHDJ2_MCLR(IASTR,IBSTR,NIDET,JASTR,JBSTR,NJDET,NAEL,NBEL,jWORK,LWORK,NORB,HAMIL,ISYM,NINOB,ECORE,ICOMBI,PSIGN,IASTRM, &
                       IBSTRM,JASTRM,JBSTRM,IGENSG,IASGN,IBSGN,JASGN,JBSGN,LIA,LIB,NDIF0,NDIF1,NDIF2)
! A set of left hand side determinants defined by string numbers
! IASTR and IBSTR and a set of right hand side determinants
! defined by JASTR and JBSTR are given.
!
! Obtain Hamiltonian matrix  < IA IB ! H ! JA JB >
!
! If Icombi /= 0 Spin combinations are assumed  for alpha and
! beta strings with different orbital configurations
!   1/SQRT(2) * (!I1A I2B! + PSIGN * !I2A I1B!)
!
! If ISYM == 0 FULL Hamiltonian is constructed
! If ISYM /= 0 LOWER half of hamiltonian is constructed
!
! JEPPE OLSEN JANUARY 1989
!
! Modifed to work with string numbers instead of strings
! March 93

use Constants, only: Zero, One, Half

implicit real*8(A-H,O-Z)
dimension IASTR(*), IBSTR(*)
dimension JASTR(*), JBSTR(*)
dimension IASTRM(NAEL,*), IBSTRM(NBEL,*)
dimension JASTRM(NAEL,*), JBSTRM(NBEL,*)
dimension IASGN(*), IBSGN(*), JASGN(*), JBSGN(*)
dimension jWORK(*), HAMIL(*)
dimension LIA(NAEL), LIB(NBEL)

! Scratch space: 4 vectors of length NORB
KLFREE = 1
KLIAE = KLFREE
KLFREE = KLIAE+NORB
KLIBE = KLFREE
KLFREE = KLIBE+NORB

KLJAE = KLFREE
KLFREE = KLJAE+NORB
KLJBE = KLFREE
KLFREE = KLJBE+NORB

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

! Loop over J determinants

! dummy initialize
JAEQJB = -1
CONST = Zero
IEL1 = -1
JEL1 = -1
SIGNA = Zero
SIGNB = Zero
IPERM = 0
JPERM = 0
SIGN = Zero
XVAL = Zero
do JDET=1,NJDET
  ! Expand JDET
  JASTAC = JASTR(JDET)
  JBSTAC = JBSTR(JDET)

  if (IGENSG > 0) then
    JXSGN = JASGN(JASTAC)*JBSGN(JBSTAC)
  else
    JXSGN = 1
  end if

  jWORK(KLJAE:KLJAE+NORB-1) = 0
  jWORK(KLJBE:KLJBE+NORB-1) = 0
  do IAEL=1,NAEL
    jWORK(KLJAE-1+JASTRM(IAEL,JASTAC)) = 1
  end do

  do IBEL=1,NBEL
    jWORK(KLJBE-1+JBSTRM(IBEL,JBSTAC)) = 1
  end do

  if (ICOMBI /= 0) then
    if (JASTAC == JBSTAC) then
      JAEQJB = 1
    else
      JAEQJB = 0
    end if
  end if

  if (ISYM == 0) then
    MINI = 1
  else
    MINI = JDET
  end if

  ! Loop over I determinants

  do IDET=MINI,NIDET
    IASTAC = IASTR(IDET)
    IBSTAC = IBSTR(IDET)

    if (IGENSG > 0) then
      IXSGN = IASGN(IASTAC)*IBSGN(IBSTAC)
    else
      IXSGN = 1
    end if

    if (IASTAC == IBSTAC) then
      IAEQIB = 1
    else
      IAEQIB = 0
    end if

    if ((ICOMBI == 1) .and. (IAEQIB+JAEQJB == 0)) then
      NLOOP = 2
    else
      NLOOP = 1
    end if
    do ILOOP=1,NLOOP
      NTERMS = NTERMS+1
      ! For second part of spin combinations strings should be swopped
      if (ILOOP == 1) then
        LIA(:) = IASTRM(:,IASTAC)
        LIB(:) = IBSTRM(:,IBSTAC)
      else if (ILOOP == 2) then
        LIB(1:NAEL) = IASTRM(:,IBSTAC)
        LIA(1:NBEL) = IBSTRM(:,IASTAC)
      end if

      ! =============================
      ! Number of orbital differences
      ! =============================

      NACM = 0
      do IAEL=1,NAEL
        NACM = NACM+jWORK(KLJAE-1+LIA(IAEL))
      end do
      NBCM = 0
      do IBEL=1,NBEL
        NBCM = NBCM+jWORK(KLJBE-1+LIB(IBEL))
      end do
      NADIF = NAEL-NACM
      NBDIF = NBEL-NBCM

      if (NADIF+NBDIF > 2) cycle
      ! Factor for combinations
      if (ICOMBI == 0) then
        CONST = One
      else
        if ((JAEQJB+IAEQIB) == 2) then
          CONST = One
        else if ((JAEQJB+IAEQIB) == 1) then
          CONST = sqrt(Half)*(One+PSIGN)
        else if ((JAEQJB+IAEQIB) == 0) then
          if (ILOOP == 1) then
            CONST = One
          else
            CONST = PSIGN
          end if
        end if
      end if
      ! External sign factor
      if (IXSGN*JXSGN == -1) CONST = -CONST

      ! ================================================
      ! Find differing orbitals and sign for permutation
      ! ================================================

      ! Expand idet
      jWORK(KLIAE:KLIAE+NORB-1) = 0
      jWORK(KLIBE:KLIBE+NORB-1) = 0

      do IAEL=1,NAEL
        jWORK(KLIAE-1+LIA(IAEL)) = 1
      end do

      do IBEL=1,NBEL
        jWORK(KLIBE-1+LIB(IBEL)) = 1
      end do

      ! One pair of differing alpha electrons

      if (NADIF == 1) then
        do IAEL=1,NAEL
          if (jWORK(KLJAE-1+LIA(IAEL)) == 0) then
            IA = LIA(IAEL)
            IEL1 = IAEL
            exit
          end if
        end do

        do JAEL=1,NAEL
          if (jWORK(KLIAE-1+JASTRM(JAEL,JASTAC)) == 0) then
            JA = JASTRM(JAEL,JASTAC)
            JEL1 = JAEL
            exit
          end if
        end do
        SIGNA = (-One)**(JEL1+IEL1)
      end if

      ! One pair of differing beta electrons

      if (NBDIF == 1) then
        do IBEL=1,NBEL
          if (jWORK(KLJBE-1+LIB(IBEL)) == 0) then
            IB = LIB(IBEL)
            IEL1 = IBEL
            exit
          end if
        end do
        do JBEL=1,NBEL
          if (jWORK(KLIBE-1+JBSTRM(JBEL,JBSTAC)) == 0) then
            JB = JBSTRM(JBEL,JBSTAC)
            JEL1 = JBEL
            exit
          end if
        end do
        SIGNB = (-One)**(JEL1+IEL1)
      end if

      ! Two pairs of differing alpha electrons

      if (NADIF == 2) then
        IDIFF = 0
        do IAEL=1,NAEL
          if (jWORK(KLJAE-1+LIA(IAEL)) == 0) then
            if (IDIFF == 0) then
              IDIFF = 1
              I1 = LIA(IAEL)
              IPERM = IAEL
            else
              I2 = LIA(IAEL)
              IPERM = IAEL+IPERM
              exit
            end if
          end if
        end do

        JDIFF = 0
        do JAEL=1,NAEL
          if (jWORK(KLIAE-1+JASTRM(JAEL,JASTAC)) == 0) then
            if (JDIFF == 0) then
              JDIFF = 1
              J1 = JASTRM(JAEL,JASTAC)
              JPERM = JAEL
            else
              J2 = JASTRM(JAEL,JASTAC)
              JPERM = JAEL+JPERM
              exit
            end if
          end if
        end do
        SIGN = (-One)**(IPERM+JPERM)
      end if

      ! Two pairs of differing beta electrons

      if (NBDIF == 2) then
        IDIFF = 0
        do IBEL=1,NBEL
          if (jWORK(KLJBE-1+LIB(IBEL)) == 0) then
            if (IDIFF == 0) then
              IDIFF = 1
              I1 = LIB(IBEL)
              IPERM = IBEL
            else
              I2 = LIB(IBEL)
              IPERM = IBEL+IPERM
              exit
            end if
          end if
        end do

        JDIFF = 0
        do JBEL=1,NBEL
          if (jWORK(KLIBE-1+JBSTRM(JBEL,JBSTAC)) == 0) then
            if (JDIFF == 0) then
              JDIFF = 1
              J1 = JBSTRM(JBEL,JBSTAC)
              JPERM = JBEL
            else
              J2 = JBSTRM(JBEL,JBSTAC)
              JPERM = JBEL+JPERM
              exit
            end if
          end if
        end do
        SIGN = (-One)**(IPERM+JPERM)
      end if

      ! =======================
      ! Value of matrix element
      ! =======================

      if ((NADIF == 2) .or. (NBDIF == 2)) then
        ! 2 differences in alpha or beta strings
        NDIF2 = NDIF2+1
        ! SIGN * (I1 J1 ! I2 J2) - (I1 J2 ! I2 J1)
        XVAL = SIGN*(GTIJKL_MCLR(I1,J1,I2,J2)-GTIJKL_MCLR(I1,J2,I2,J1))
      else if ((NADIF == 1) .and. (NBDIF == 1)) then
        ! 1 difference in alpha strings and one difference in beta string
        NDIF2 = NDIF2+1
        ! SIGN * (IA JA ! IB JB)
        XVAL = SIGNA*SIGNB*GTIJKL_MCLR(IA,JA,IB,JB)
        ! 1 differences in alpha or beta strings
      else if ((NADIF == 1) .and. (NBDIF == 0) .or. (NADIF == 0) .and. (NBDIF == 1)) then
        NDIF1 = NDIF1+1
        ! SIGN *(H(I1 J1) +
        !  (SUM OVER ORBITALS OF BOTH      SPIN TYPES  (I1 J1 ! JORB JORB)
        ! -(SUM OVER ORBITALS OF DIFFERING SPIN TYPE   (I1 JORB ! JORB J1)
        if (NADIF == 1) then
          I1 = IA
          J1 = JA
          SIGN = SIGNA
        else
          I1 = IB
          J1 = JB
          SIGN = SIGNB
        end if

        XVAL = GETH1I_MCLR(I1,J1)
        do JAEL=1,NAEL
          JORB = JASTRM(JAEL,JASTAC)
          XVAL = XVAL+GTIJKL_MCLR(I1,J1,JORB,JORB)
        end do
        do JBEL=1,NBEL
          JORB = JBSTRM(JBEL,JBSTAC)
          XVAL = XVAL+GTIJKL_MCLR(I1,J1,JORB,JORB)
        end do
        if (NADIF == 1) then
          do JAEL=1,NAEL
            JORB = JASTRM(JAEL,JASTAC)
            XVAL = XVAL-GTIJKL_MCLR(I1,JORB,JORB,J1)
          end do
        else
          do JBEL=1,NBEL
            JORB = JBSTRM(JBEL,JBSTAC)
            XVAL = XVAL-GTIJKL_MCLR(I1,JORB,JORB,J1)
          end do
        end if
        XVAL = XVAL*SIGN
      else if ((NADIF == 0) .and. (NBDIF == 0)) then
        ! Diagonal elements
        NDIF0 = NDIF0+1
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
                IORB = LIA(IEL)
              else
                IORB = LIB(IEL)
              end if
              if (IAB == JAB) XVAL = XVAL+GETH1I_MCLR(IORB,IORB)
              do JEL=1,NJABEL
                if (JAB == 1) then
                  JORB = LIA(JEL)
                else
                  JORB = LIB(JEL)
                end if
                XVAL = XVAL+Half*GTIJKL_MCLR(IORB,IORB,JORB,JORB)
                ! test

                if (IAB == JAB) XVAL = XVAL-Half*GTIJKL_MCLR(IORB,JORB,JORB,IORB)
                ! test
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
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(LWORK)
  call Unused_integer(NINOB)
end if

end subroutine DIHDJ2_MCLR
