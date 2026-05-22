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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

module SUPERINDEX

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), allocatable :: KAGEB(:,:), KAGTB(:,:), KIA(:,:), KIGEJ(:,:), KIGTJ(:,:), KTGEU(:,:), KTGTU(:,:), KTU(:,:), &
                                  KTUV(:,:,:), MAGEB(:,:), MAGTB(:,:), MAREL(:,:), MIA(:,:), MIGEJ(:,:), MIGTJ(:,:), MIREL(:,:), &
                                  MTGEU(:,:), MTGTU(:,:), MTREL(:,:), MTU(:,:), MTUV(:,:)

public :: KAGEB, KAGTB, KIGEJ, KIGTJ, KTGEU, KTGTU, KTU, KTUV, MAGEB, MAGTB, MAREL, MIA, MIGEJ, MIGTJ, MIREL, MTGEU, MTGTU, MTREL, &
          MTU, MTUV, SUPFREE, SUPINI

contains

subroutine SUPINI()

  use Symmetry_Info, only: Mul
  use caspt2_module, only: nAes, nAGEB, nAGEBES, nAGTB, nAGTBES, nAsh, nAshT, nAshT, nASUP, nCases, nIAES, nIES, nIGEJ, nIGEJES, &
                           nIGTJ, nIGTJES, nInDep, nIsh, nIshT, nISUP, nSES, nSES, nSsh, nSsh, nSshT, nSym, nTGEU, nTGEUES, nTGTU, &
                           nTGTUES, nTU, nTUES, nTUV, nTUVEs
  use stdalloc, only: mma_allocate

  integer(kind=iwp) :: IA, IAGEB, IAGTB, IAQ, IB, IBQ, ICASE, II, IIA, IIGEJ, IIGTJ, IIQ, IJ, IJQ, IS1, IS2, ISA, ISB, ISI, ISJ, &
                       IST, ISU, ISUV, ISV, ISYA, ISYI, ISYM, IT, ITGEU, ITGTU, ITQ, ITU, ITUV, IU, IUQ, IV, IVQ, JC0, JCM, &
                       JCOUNT, JCP, N, N10, N11, N5, N6, N7, NAT, NAU, NAV, NC0, NCM, NCOUNT, NCP, NII, NIJ, NMAGEB, NMAGTB, NMIA, &
                       NMIGEJ, NMIGTJ, NMTGEU, NMTGTU, NSA, NSB, NSUM

  call MMA_ALLOCATE(KTUV,NASHT,NASHT,NASHT,Label='KTUV')
  call MMA_ALLOCATE(MTUV,3,NASHT**3,Label='MTUV')
  ITUV = 0
  do ISYM=1,NSYM
    NCOUNT = 0
    NTUVES(ISYM) = ITUV
    do ISV=1,NSYM
      NAV = NASH(ISV)
      do ISU=1,NSYM
        NAU = NASH(ISU)
        ISUV = Mul(ISU,ISV)
        IST = Mul(ISUV,ISYM)
        NAT = NASH(IST)
        JCOUNT = NAV*NAU*NAT
        if (JCOUNT == 0) cycle
        NCOUNT = NCOUNT+JCOUNT
        do IV=1,NAV
          IVQ = NAES(ISV)+IV
          do IU=1,NAU
            IUQ = NAES(ISU)+IU
            do IT=1,NAT
              ITQ = NAES(IST)+IT
              ITUV = ITUV+1
              KTUV(ITQ,IUQ,IVQ) = ITUV
              MTUV(1,ITUV) = ITQ
              MTUV(2,ITUV) = IUQ
              MTUV(3,ITUV) = IVQ
            end do
          end do
        end do
      end do
    end do
    NTUV(ISYM) = NCOUNT
  end do

  call MMA_ALLOCATE(KTU,NASHT,NASHT,Label='KTU')
  call MMA_ALLOCATE(MTU,3,NASHT**2,Label='MTU')

  call MMA_ALLOCATE(KTGEU,NASHT,NASHT,Label='KTGEU')
  call MMA_ALLOCATE(KTGTU,NASHT,NASHT,Label='KTGTU')
  NMTGEU = (NASHT*(NASHT+1))/2
  NMTGTU = (NASHT*(NASHT-1))/2
  call MMA_ALLOCATE(MTGEU,2,NMTGEU,Label='MTGEU')
  call MMA_ALLOCATE(MTGTU,2,NMTGTU,Label='MTGTU')

  ITU = 0
  ITGEU = 0
  ITGTU = 0
  do ISYM=1,NSYM
    NC0 = 0
    NCP = 0
    NCM = 0
    NTUES(ISYM) = ITU
    NTGEUES(ISYM) = ITGEU
    NTGTUES(ISYM) = ITGTU
    do ISU=1,NSYM
      NAU = NASH(ISU)
      IST = Mul(ISU,ISYM)
      NAT = NASH(IST)
      JC0 = 0
      JCP = 0
      JCM = 0
      do IU=1,NAU
        IUQ = IU+NAES(ISU)
        do IT=1,NAT
          ITQ = IT+NAES(IST)
          JC0 = JC0+1
          ITU = ITU+1
          KTU(ITQ,IUQ) = ITU
          MTU(1,ITU) = ITQ
          MTU(2,ITU) = IUQ
          if (ITQ < IUQ) cycle
          JCP = JCP+1
          ITGEU = ITGEU+1
          KTGEU(ITQ,IUQ) = ITGEU
          MTGEU(1,ITGEU) = ITQ
          MTGEU(2,ITGEU) = IUQ
          if (ITQ <= IUQ) cycle
          JCM = JCM+1
          ITGTU = ITGTU+1
          KTGTU(ITQ,IUQ) = ITGTU
          MTGTU(1,ITGTU) = ITQ
          MTGTU(2,ITGTU) = IUQ
        end do
      end do
      NC0 = NC0+JC0
      NCP = NCP+JCP
      NCM = NCM+JCM
    end do
    NTU(ISYM) = NC0
    NTGEU(ISYM) = NCP
    NTGTU(ISYM) = NCM
  end do

  !PAM99 Use allocated workspace instead of MAGEB, MAGTB:
  NMAGEB = (NSSHT*(NSSHT+1))/2
  NMAGTB = (NSSHT*(NSSHT-1))/2
  call MMA_ALLOCATE(MAGEB,2,NMAGEB,Label='MAGEB')
  call MMA_ALLOCATE(MAGTB,2,NMAGTB,Label='MAGTB')
  NMIGEJ = (NISHT*(NISHT+1))/2
  NMIGTJ = (NISHT*(NISHT-1))/2
  call MMA_ALLOCATE(MIGEJ,2,NMIGEJ,Label='MIGEJ')
  call MMA_ALLOCATE(MIGTJ,2,NMIGTJ,Label='MIGTJ')

  call MMA_ALLOCATE(KIGEJ,NISHT,NISHT,Label='KIGEJ')
  call MMA_ALLOCATE(KIGTJ,NISHT,NISHT,Label='KIGTJ')
  call MMA_ALLOCATE(KAGEB,NSSHT,NSSHT,Label='KAGEB')
  call MMA_ALLOCATE(KAGTB,NSSHT,NSSHT,Label='KAGTB')
  KIGEJ(:,:) = 0
  KIGTJ(:,:) = 0

  call MMA_ALLOCATE(KIA,NISHT,NSSHT,Label='KIA')
  NMIA = NISHT*NSSHT
  call MMA_ALLOCATE(MIA,2,NMIA,Label='MIA')

  ! Construct tables for inactive and secondary pair indices:
  IIGEJ = 0
  IIGTJ = 0
  IAGEB = 0
  IAGTB = 0
  IIA = 0
  do ISYM=1,NSYM
    ! Inactive pair indices:
    NIGEJES(ISYM) = IIGEJ
    NIGTJES(ISYM) = IIGTJ
    NCM = 0
    NCP = 0
    do ISI=1,NSYM
      ISJ = Mul(ISI,ISYM)
      if (ISI < ISJ) cycle
      NII = NISH(ISI)
      NIJ = NISH(ISJ)
      do II=1,NII
        IIQ = NIES(ISI)+II
        do IJ=1,NIJ
          IJQ = NIES(ISJ)+IJ
          NCP = NCP+1
          IIGEJ = IIGEJ+1
          KIGEJ(IIQ,IJQ) = IIGEJ
          MIGEJ(1,IIGEJ) = IIQ
          MIGEJ(2,IIGEJ) = IJQ
          if (IIQ <= IJQ) exit
          NCM = NCM+1
          IIGTJ = IIGTJ+1
          KIGTJ(IIQ,IJQ) = IIGTJ
          MIGTJ(1,IIGTJ) = IIQ
          MIGTJ(2,IIGTJ) = IJQ
        end do
      end do
    end do

    NIGEJ(ISYM) = NCP
    NIGTJ(ISYM) = NCM

    ! Secondary pair indices:
    NAGEBES(ISYM) = IAGEB
    NAGTBES(ISYM) = IAGTB
    NCM = 0
    NCP = 0
    do ISA=1,NSYM
      ISB = Mul(ISA,ISYM)
      if (ISA < ISB) cycle
      NSA = NSSH(ISA)
      NSB = NSSH(ISB)
      do IA=1,NSA
        IAQ = NSES(ISA)+IA
        do IB=1,NSB
          IBQ = NSES(ISB)+IB
          NCP = NCP+1
          IAGEB = IAGEB+1
          KAGEB(IAQ,IBQ) = IAGEB
          MAGEB(1,IAGEB) = IAQ
          MAGEB(2,IAGEB) = IBQ
          if (IAQ <= IBQ) exit
          NCM = NCM+1
          IAGTB = IAGTB+1
          KAGTB(IAQ,IBQ) = IAGTB
          MAGTB(1,IAGTB) = IAQ
          MAGTB(2,IAGTB) = IBQ
        end do
      end do
    end do
    NAGEB(ISYM) = NCP
    NAGTB(ISYM) = NCM

    ! Inactive-Secondary pair indices:
    NIAES(ISYM) = IIA
    do ISYA=1,NSYM
      ISYI = Mul(ISYA,ISYM)
      do IA=1,NSSH(ISYA)
        IAQ = IA+NSES(ISYA)
        do II=1,NISH(ISYI)
          IIQ = II+NIES(ISYI)
          IIA = IIA+1
          KIA(IIQ,IAQ) = IIA
          MIA(1,IIA) = IIQ
          MIA(2,IIA) = IAQ
        end do
      end do
    end do

  end do
  ! End of loop over symmetries.

  do ISYM=1,NSYM
    NASUP(ISYM,1) = NTUV(ISYM)
    NASUP(ISYM,2) = NTGEU(ISYM)
    NASUP(ISYM,3) = NTGTU(ISYM)
    NASUP(ISYM,4) = NTUV(ISYM)
    NASUP(ISYM,5) = 2*NTU(ISYM)
    NASUP(ISYM,6) = NASH(ISYM)
    NASUP(ISYM,7) = NASH(ISYM)
    NASUP(ISYM,8) = NTGEU(ISYM)
    NASUP(ISYM,9) = NTGTU(ISYM)
    NASUP(ISYM,10) = NASH(ISYM)
    NASUP(ISYM,11) = NASH(ISYM)
    NASUP(ISYM,12) = NAGEB(ISYM)
    NASUP(ISYM,13) = NAGTB(ISYM)
    N5 = 0
    N6 = 0
    N7 = 0
    N10 = 0
    N11 = 0
    do IS1=1,NSYM
      IS2 = Mul(IS1,ISYM)
      N5 = N5+NSSH(IS1)*NISH(IS2)
      N6 = N6+NSSH(IS1)*NIGEJ(IS2)
      N7 = N7+NSSH(IS1)*NIGTJ(IS2)
      N10 = N10+NAGEB(IS1)*NISH(IS2)
      N11 = N11+NAGTB(IS1)*NISH(IS2)
    end do
    NISUP(ISYM,1) = NISH(ISYM)
    NISUP(ISYM,2) = NIGEJ(ISYM)
    NISUP(ISYM,3) = NIGTJ(ISYM)
    NISUP(ISYM,4) = NSSH(ISYM)
    NISUP(ISYM,5) = N5
    NISUP(ISYM,6) = N6
    NISUP(ISYM,7) = N7
    NISUP(ISYM,8) = NAGEB(ISYM)
    NISUP(ISYM,9) = NAGTB(ISYM)
    NISUP(ISYM,10) = N10
    NISUP(ISYM,11) = N11
    NISUP(ISYM,12) = NIGEJ(ISYM)
    NISUP(ISYM,13) = NIGTJ(ISYM)
  end do
  do ICASE=1,NCASES
    NSUM = 0
    do ISYM=1,NSYM
      N = NASUP(ISYM,ICASE)*NISUP(ISYM,ICASE)
      ! Preliminary value for NINDEP: Nr of independent active params:
      NINDEP(ISYM,ICASE) = NASUP(ISYM,ICASE)
      if (N == 0) NINDEP(ISYM,ICASE) = 0
      NSUM = NSUM+N
    end do
  end do

  !SVC: prepare tables to translate from absolute indices to
  !(index,symmetry) pairs.

  call MMA_ALLOCATE(MIREL,2,NISHT,Label='MIREL')
  call MMA_ALLOCATE(MTREL,2,NASHT,Label='MTREL')
  call MMA_ALLOCATE(MAREL,2,NSSHT,Label='MAREL')

  do ISYM=1,NSYM
    do II=1,NISH(ISYM)
      IIQ = II+NIES(ISYM)
      MIREL(1,IIQ) = II
      MIREL(2,IIQ) = ISYM
    end do
    do IT=1,NASH(ISYM)
      ITQ = IT+NAES(ISYM)
      MTREL(1,ITQ) = IT
      MTREL(2,ITQ) = ISYM
    end do
    do IA=1,NSSH(ISYM)
      IAQ = IA+NSES(ISYM)
      MAREL(1,IAQ) = IA
      MAREL(2,IAQ) = ISYM
    end do
  end do

end subroutine SUPINI

subroutine SUPFREE()

  use stdalloc, only: mma_deallocate

  ! deallocate the superindex tables
  call MMA_DEALLOCATE(KIGEJ)
  call MMA_DEALLOCATE(KIGTJ)
  call MMA_DEALLOCATE(MIGEJ)
  call MMA_DEALLOCATE(MIGTJ)
  call MMA_DEALLOCATE(MAGEB)
  call MMA_DEALLOCATE(MAGTB)
  call MMA_DEALLOCATE(KAGEB)
  call MMA_DEALLOCATE(KAGTB)
  call MMA_DEALLOCATE(MTGEU)
  call MMA_DEALLOCATE(MTGTU)
  call MMA_DEALLOCATE(KTGEU)
  call MMA_DEALLOCATE(KTGTU)
  call MMA_DEALLOCATE(KTU)
  call MMA_DEALLOCATE(MTU)
  call MMA_DEALLOCATE(KTUV)
  call MMA_DEALLOCATE(MTUV)
  call MMA_DEALLOCATE(KIA)
  call MMA_DEALLOCATE(MIA)
  call MMA_DEALLOCATE(MIREL)
  call MMA_DEALLOCATE(MTREL)
  call MMA_DEALLOCATE(MAREL)

end subroutine SUPFREE

end module SUPERINDEX
