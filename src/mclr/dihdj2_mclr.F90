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

subroutine DIHDJ2_MCLR(IASTR,IBSTR,NIDET,JASTR,JBSTR,NJDET,NAEL,NBEL,jWORK,NORB,HAMIL,ISYM,ECORE,ICOMBI,PSIGN,IASTRM,IBSTRM, &
                       JASTRM,JBSTRM,LIA,LIB,NDIF0,NDIF1,NDIF2)
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

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IASTR(*), IBSTR(*), NIDET, JASTR(*), JBSTR(*), NJDET, NAEL, NBEL, NORB, ISYM, ICOMBI, &
                                 IASTRM(NAEL,*), IBSTRM(NBEL,*), JASTRM(NAEL,*), JBSTRM(NBEL,*)
integer(kind=iwp), intent(out) :: jWORK(NORB,4), LIA(max(NAEL,NBEL)), LIB(max(NAEL,NBEL)), NDIF0, NDIF1, NDIF2
real(kind=wp), intent(out) :: HAMIL(*)
real(kind=wp), intent(in) :: ECORE, PSIGN
integer(kind=iwp) :: I1, I2, IA, IAB, IAEL, IAEQIB, IASTAC, IB, IBEL, IBSTAC, IDET, IDIFF, IEL, IEL1, ILOOP, IORB, IPERM, J1, J2, &
                     JA, JAB, JAEL, JAEQJB, JASTAC, JB, JBEL, JBSTAC, JDET, JDIFF, JEL, JEL1, JORB, JPERM, LHAMIL, MINI, NACM, &
                     NADIF, NBCM, NBDIF, NIABEL, NJABEL, NLOOP, NTERMS
real(kind=wp) :: CONST, SGN, SIGNA, SIGNB, XVAL
real(kind=wp), external :: GETH1I_MCLR, GTIJKL_MCLR

if (ISYM == 0) then
  LHAMIL = NIDET*NJDET
else
  LHAMIL = nTri_Elem(NIDET)
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
SGN = Zero
XVAL = Zero
do JDET=1,NJDET
  ! Expand JDET
  JASTAC = JASTR(JDET)
  JBSTAC = JBSTR(JDET)

  jWORK(:,3:4) = 0
  do IAEL=1,NAEL
    jWORK(JASTRM(IAEL,JASTAC),3) = 1
  end do

  do IBEL=1,NBEL
    jWORK(JBSTRM(IBEL,JBSTAC),4) = 1
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
      ! For second part of spin combinations strings should be swapped
      if (ILOOP == 1) then
        LIA(1:NAEL) = IASTRM(:,IASTAC)
        LIB(1:NBEL) = IBSTRM(:,IBSTAC)
      else if (ILOOP == 2) then
        LIB(1:NAEL) = IASTRM(:,IBSTAC)
        LIA(1:NBEL) = IBSTRM(:,IASTAC)
      end if

      ! =============================
      ! Number of orbital differences
      ! =============================

      NACM = 0
      do IAEL=1,NAEL
        NACM = NACM+jWORK(LIA(IAEL),3)
      end do
      NBCM = 0
      do IBEL=1,NBEL
        NBCM = NBCM+jWORK(LIB(IBEL),4)
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
      !if (IXSGN*JXSGN == -1) CONST = -CONST

      ! ================================================
      ! Find differing orbitals and sign for permutation
      ! ================================================

      ! Expand idet
      jWORK(:,1:2) = 0

      do IAEL=1,NAEL
        jWORK(LIA(IAEL),1) = 1
      end do

      do IBEL=1,NBEL
        jWORK(LIB(IBEL),2) = 1
      end do

      ! One pair of differing alpha electrons

      if (NADIF == 1) then
        do IAEL=1,NAEL
          if (jWORK(LIA(IAEL),3) == 0) then
            IA = LIA(IAEL)
            IEL1 = IAEL
            exit
          end if
        end do

        do JAEL=1,NAEL
          if (jWORK(JASTRM(JAEL,JASTAC),1) == 0) then
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
          if (jWORK(LIB(IBEL),4) == 0) then
            IB = LIB(IBEL)
            IEL1 = IBEL
            exit
          end if
        end do
        do JBEL=1,NBEL
          if (jWORK(JBSTRM(JBEL,JBSTAC),2) == 0) then
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
          if (jWORK(LIA(IAEL),3) == 0) then
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
          if (jWORK(JASTRM(JAEL,JASTAC),1) == 0) then
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
        SGN = (-One)**(IPERM+JPERM)
      end if

      ! Two pairs of differing beta electrons

      if (NBDIF == 2) then
        IDIFF = 0
        do IBEL=1,NBEL
          if (jWORK(LIB(IBEL),4) == 0) then
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
          if (jWORK(JBSTRM(JBEL,JBSTAC),2) == 0) then
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
        SGN = (-One)**(IPERM+JPERM)
      end if

      ! =======================
      ! Value of matrix element
      ! =======================

      if ((NADIF == 2) .or. (NBDIF == 2)) then
        ! 2 differences in alpha or beta strings
        NDIF2 = NDIF2+1
        ! SGN * (I1 J1 ! I2 J2) - (I1 J2 ! I2 J1)
        XVAL = SGN*(GTIJKL_MCLR(I1,J1,I2,J2)-GTIJKL_MCLR(I1,J2,I2,J1))
      else if ((NADIF == 1) .and. (NBDIF == 1)) then
        ! 1 difference in alpha strings and one difference in beta string
        NDIF2 = NDIF2+1
        ! SGN * (IA JA ! IB JB)
        XVAL = SIGNA*SIGNB*GTIJKL_MCLR(IA,JA,IB,JB)
        ! 1 differences in alpha or beta strings
      else if ((NADIF == 1) .and. (NBDIF == 0) .or. (NADIF == 0) .and. (NBDIF == 1)) then
        NDIF1 = NDIF1+1
        ! SGN *(H(I1 J1) +
        !  (SUM OVER ORBITALS OF BOTH      SPIN TYPES  (I1 J1 ! JORB JORB)
        ! -(SUM OVER ORBITALS OF DIFFERING SPIN TYPE   (I1 JORB ! JORB J1)
        if (NADIF == 1) then
          I1 = IA
          J1 = JA
          SGN = SIGNA
        else
          I1 = IB
          J1 = JB
          SGN = SIGNB
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
        XVAL = XVAL*SGN
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
        HAMIL(nTri_Elem(IDET-1)+JDET) = HAMIL(nTri_Elem(IDET-1)+JDET)+CONST*XVAL
      end if
    end do
  end do
end do

end subroutine DIHDJ2_MCLR
