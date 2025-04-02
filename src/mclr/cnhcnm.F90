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
! Copyright (C) 1990, Jeppe Olsen                                      *
!***********************************************************************

subroutine CNHCNM(HSUB,ISYM,ILCNF,NLCNF,IRCNF,NRCNF,NLCSF,SCR,ICONF,NEL,IREFSM,NAEL,NBEL,NACOB,IPRODT,DTOC,INTSPC,ICOMBI,PSSIGN)
! Calculate  Hamiltonian block defined by configuration
! lists ILCNF,IRCNF
! If ISYM /= 0 only the lower half of the matrix is constructed
!
! Jeppe Olsen April 1990
! ========================
!
! No modifications
! ================

use iso_c_binding, only: c_f_pointer, c_loc
use Index_Functions, only: nTri_Elem
use MCLR_Data, only: NTYP, NCPCNT

implicit none
integer ISYM, NLCNF, NRCNF, NLCSF, NEL, IREFSM, NAEL, NBEL, NACOB, INTSPC, ICOMBI
real*8 PSSIGN
! Specific input
integer ILCNF(*), IRCNF(*)
! General input
integer ICONF(*), IPRODT(*)
real*8 DTOC(*)
! Output
real*8 HSUB(*)
! Scratch
real*8 SCR(*)
! Length of scratch: 2 * NEL + MXCSFC                   (used in CNHCNM)
!                  + 6*MXDTFC+MXDTFC**2+MXDTFC+MXCSFC   (used in CNHCN2)
!                  + MAX(MXDTFC*NEL+2*NEL,4*NORB+2*NEL) (used in DIHDJ,CNFSTR)

call CNHCNM_INTERNAL(SCR)

return

! This is to allow type punning without an explicit interface
contains

subroutine CNHCNM_INTERNAL(SCR)

  real*8, target :: SCR(*)
  integer, pointer :: iSCRl(:), iSCRr(:)
  integer NDIF0, NDIF1, NDIF2, MXCSFC, ITYP, KLFREE, KLCONF, KRCONF, KLPHPS, IILB, ICNL, NCSFL, IIRB, MXR, ICNR, NCSFR, MDIF0, &
          MDIF1, MDIF2, IIL, IIRMAX, IIR, IIRACT, IILACT, ILRO, ILRI, ILTYP, IRTYP

  NDIF0 = 0
  NDIF1 = 0
  NDIF2 = 0
  ! Largest configuration block possible
  MXCSFC = 0
  do ITYP=1,NTYP
    MXCSFC = max(MXCSFC,NCPCNT(ITYP))
  end do

  KLFREE = 1

  KLCONF = KLFREE
  KLFREE = KLFREE+NEL

  KRCONF = KLFREE
  KLFREE = KLFREE+NEL

  KLPHPS = KLFREE
  KLFREE = KLFREE+MXCSFC**2

  ! LHR
  IILB = 1
  call c_f_pointer(c_loc(SCR(KLCONF)),iSCRl,[1])
  call c_f_pointer(c_loc(SCR(KRCONF)),iSCRr,[1])
  do ICNL=1,NLCNF
    call GETCNF_MCLR(iSCRl,ILTYP,ILCNF(ICNL),ICONF,IREFSM,NEL)
    NCSFL = NCPCNT(ILTYP)
    IIRB = 1
    if (ISYM == 0) then
      MXR = NRCNF
    else
      MXR = ICNL
    end if
    do ICNR=1,MXR
      call GETCNF_MCLR(iSCRr,IRTYP,IRCNF(ICNR),ICONF,IREFSM,NEL)
      NCSFR = NCPCNT(IRTYP)
      call CNHCN2(iSCRl,ILTYP,iSCRr,IRTYP,SCR(KLPHPS),SCR(KLFREE),NEL,NAEL,NBEL,INTSPC,IPRODT,DTOC,NACOB,ICOMBI,PSSIGN,MDIF0, &
                  MDIF1,MDIF2)
      NDIF0 = NDIF0+MDIF0
      NDIF1 = NDIF1+MDIF1
      NDIF2 = NDIF2+MDIF2

      ! Copy to HSUB matrix
      if (ISYM /= 0) then
        ! Copy to lower half format
        do IIL=1,NCSFL
          if (IILB == IIRB) then
            IIRMAX = IIL
          else
            IIRMAX = NCSFR
          end if
          do IIR=1,IIRMAX
            IIRACT = IIRB-1+IIR
            IILACT = IILB-1+IIL
            ILRO = nTri_Elem(IILACT-1)+IIRACT
            ILRI = (IIR-1)*NCSFL+IIL
            HSUB(ILRO) = SCR(KLPHPS-1+ILRI)
          end do
        end do
      else
        ! Pack to full format
        do IIL=1,NCSFL
          do IIR=1,NCSFR
            IIRACT = IIRB-1+IIR
            IILACT = IILB-1+IIL
            ILRO = (IIRACT-1)*NLCSF+IILACT
            ILRI = (IIR-1)*NCSFL+IIL
            HSUB(ILRO) = SCR(KLPHPS-1+ILRI)
          end do
        end do
      end if
      IIRB = IIRB+NCSFR
    end do
    IILB = IILB+NCSFL
  end do
  nullify(iSCRl,iSCRr)

end subroutine CNHCNM_INTERNAL

end subroutine CNHCNM
