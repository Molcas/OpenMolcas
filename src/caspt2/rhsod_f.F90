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
! Copyright (C) Steven Vancoillie                                      *
!***********************************************************************

subroutine RHSOD_F(IVEC)

use Symmetry_Info, only: Mul
use SUPERINDEX, only: MAGEB, MAGTB, MAREL, MAREL, MTGEU, MTGTU, MTREL
use CHOVEC_IO, only: ChoVec_Read, ChoVec_Size, NVTOT_CHOSYM
use PrintLevel, only: DEBUG
#ifndef _MOLCAS_MPP_
use fake_GA, only: GA_Arrays
#endif
use caspt2_global, only: iPrGlb
use caspt2_module, only: NAGEBES, NAGTBES, NASUP, NISUP, NSSH, NSYM, NTGEUES, NTGTUES
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVEC
integer(kind=iwp) :: IA, IAABS, IAEND, IAGEB, IAGEBTOT, IAGTB, IAGTBTOT, IASTA, IAT, IAV, IC, ICABS, iCASE, ICT, ICV, IDX, IIEND, &
                     IISTA, IOFFAT, IOFFAV, IOFFCT, IOFFCV, IOSYM(8,8), ISYA, ISYC, ISYM, ISYT, ISYV, IT, ITABS, ITGEU, ITGEUTOT, &
                     ITGTU, ITGTUTOT, IV, IVABS, lg_W, MW, NAS, NCHOBUF, NIS, NV, NW
real(kind=wp) :: ATCV, AVCT, FMTVAC, FPTVAC, SCL
real(kind=wp), allocatable :: CHOBUF(:)
real(kind=wp), parameter :: SQRTH = sqrt(Half)
real(kind=wp), external :: DDot_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

if (iPrGlb >= DEBUG) write(u6,*) 'RHS on demand: case F'

!***********************************************************************
! Case F (8,9):
! FP(tv,ac)=((at,cv)+(av,ct))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(a,c))
! FM(tv,ac)= -((at,cv)-(av,ct))/(2*SQRT(1+Kron(a,c))
!***********************************************************************

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
call CHOVEC_SIZE(3,NCHOBUF,IOSYM)

call mma_allocate(CHOBUF,NCHOBUF,LABEL='CHOBUF')

call CHOVEC_READ(3,CHOBUF,NCHOBUF)

iCASE = 8
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
do ISYM=1,NSYM

  NAS = NASUP(ISYM,ICASE)
  NIS = NISUP(ISYM,ICASE)
  NW = NAS*NIS

  if (NW == 0) cycle

  call RHS_ALLO(NAS,NIS,lg_W)
  call RHS_ACCESS(NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
  NW = NAS*(IIEND-IISTA+1)

  !*********************************************************************
  ! inner loop over RHS elements in symmetry ISYM
  !*********************************************************************
  do IAGEB=IISTA,IIEND
    IAGEBTOT = IAGEB+NAGEBES(ISYM)
    IAABS = MAGEB(1,IAGEBTOT)
    ICABS = MAGEB(2,IAGEBTOT)
    IA = MAREL(1,IAABS)
    ISYA = MAREL(2,IAABS)
    IC = MAREL(1,ICABS)
    ISYC = MAREL(2,ICABS)
    do ITGEU=IASTA,IAEND ! these are always all elements
      ITGEUTOT = ITGEU+NTGEUES(ISYM)
      ITABS = MTGEU(1,ITGEUTOT)
      IVABS = MTGEU(2,ITGEUTOT)
      IT = MTREL(1,ITABS)
      ISYT = MTREL(2,ITABS)
      IV = MTREL(1,IVABS)
      ISYV = MTREL(2,IVABS)
      ! compute integrals (ta,vc) and (tc,va)
      NV = NVTOT_CHOSYM(Mul(ISYA,ISYT)) ! JSYM=ISYA*ISYA=ISYC*ISYC
      IAT = IA-1+NSSH(ISYA)*(IT-1)
      ICV = IC-1+NSSH(ISYC)*(IV-1)
      IOFFAT = 1+IOSYM(ISYA,ISYT)+NV*IAT
      IOFFCV = 1+IOSYM(ISYC,ISYV)+NV*ICV
      ATCV = DDOT_(NV,CHOBUF(IOFFAT),1,CHOBUF(IOFFCV),1)

      NV = NVTOT_CHOSYM(Mul(ISYA,ISYV)) ! JSYM=ISYA*ISYA=ISYC*ISYC
      IAV = IA-1+NSSH(ISYA)*(IV-1)
      ICT = IC-1+NSSH(ISYC)*(IT-1)
      IOFFAV = 1+IOSYM(ISYA,ISYV)+NV*IAV
      IOFFCT = 1+IOSYM(ISYC,ISYT)+NV*ICT
      AVCT = DDOT_(NV,CHOBUF(IOFFAV),1,CHOBUF(IOFFCT),1)

      ! FP(tv,ac)=((at,cv)+(av,ct))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(a,c))
      SCL = Half
      if (ITABS == IVABS) SCL = SCL*Half
      if (IAABS == ICABS) SCL = SCL*SQRTH
      FPTVAC = SCL*(ATCV+AVCT)
      ! write element FP(tv,ac)
      IDX = ITGEU+NAS*(IAGEB-IISTA)
#     ifdef _MOLCAS_MPP_
      DBL_MB(MW+IDX-1) = FPTVAC
#     else
      GA_Arrays(lg_w)%A(IDX) = FPTVAC
#     endif
    end do
  end do
  !*********************************************************************

  call RHS_RELEASE_UPDATE(lg_W,IASTA,IAEND,IISTA,IIEND)
  call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_W)
end do
!***********************************************************************

iCASE = 9
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
do ISYM=1,NSYM

  NAS = NASUP(ISYM,ICASE)
  NIS = NISUP(ISYM,ICASE)
  NW = NAS*NIS

  if (NW == 0) cycle

  call RHS_ALLO(NAS,NIS,lg_W)
  call RHS_ACCESS(NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
  NW = NAS*(IIEND-IISTA+1)

  !*********************************************************************
  ! inner loop over RHS elements in symmetry ISYM
  !*********************************************************************
  do IAGTB=IISTA,IIEND
    IAGTBTOT = IAGTB+NAGTBES(ISYM)
    IAABS = MAGTB(1,IAGTBTOT)
    ICABS = MAGTB(2,IAGTBTOT)
    IA = MAREL(1,IAABS)
    ISYA = MAREL(2,IAABS)
    IC = MAREL(1,ICABS)
    ISYC = MAREL(2,ICABS)
    do ITGTU=IASTA,IAEND ! these are always all elements
      ITGTUTOT = ITGTU+NTGTUES(ISYM)
      ITABS = MTGTU(1,ITGTUTOT)
      IVABS = MTGTU(2,ITGTUTOT)
      IT = MTREL(1,ITABS)
      ISYT = MTREL(2,ITABS)
      IV = MTREL(1,IVABS)
      ISYV = MTREL(2,IVABS)
      ! compute integrals (at,cv) and (av,ct)
      NV = NVTOT_CHOSYM(Mul(ISYA,ISYT)) ! JSYM=ISYA*ISYA=ISYC*ISYC
      IAT = IA-1+NSSH(ISYA)*(IT-1)
      ICV = IC-1+NSSH(ISYC)*(IV-1)
      IOFFAT = 1+IOSYM(ISYA,ISYT)+NV*IAT
      IOFFCV = 1+IOSYM(ISYC,ISYV)+NV*ICV
      ATCV = DDOT_(NV,CHOBUF(IOFFAT),1,CHOBUF(IOFFCV),1)

      NV = NVTOT_CHOSYM(Mul(ISYA,ISYV)) ! JSYM=ISYA*ISYA=ISYC*ISYC
      IAV = IA-1+NSSH(ISYA)*(IV-1)
      ICT = IC-1+NSSH(ISYC)*(IT-1)
      IOFFAV = 1+IOSYM(ISYA,ISYV)+NV*IAV
      IOFFCT = 1+IOSYM(ISYC,ISYT)+NV*ICT
      AVCT = DDOT_(NV,CHOBUF(IOFFAV),1,CHOBUF(IOFFCT),1)

      ! FM(tv,ac)= -((at,cv)-(av,ct))/(2*SQRT(1+Kron(a,c))
      SCL = Half
      !if (ITABS == IVABS) SCL = SCL*Half
      !if (IAABS == ICABS) SCL = SCL*SQRTH
      FMTVAC = SCL*(AVCT-ATCV)
      ! write element FM(tv,ac)
      IDX = ITGTU+NAS*(IAGTB-IISTA)
#     ifdef _MOLCAS_MPP_
      DBL_MB(MW+IDX-1) = FMTVAC
#     else
      GA_Arrays(lg_w)%A(IDX) = FMTVAC
#     endif
    end do
  end do
  !*********************************************************************

  call RHS_RELEASE_UPDATE(lg_W,IASTA,IAEND,IISTA,IIEND)
  call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_W)
end do
!***********************************************************************

call mma_deallocate(CHOBUF)

end subroutine RHSOD_F
