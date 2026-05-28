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

subroutine RHSOD_H(IVEC)

use Symmetry_Info, only: Mul
use SUPERINDEX, only: MAGEB, MAGTB, MAREL, MIGEJ, MIGTJ, MIREL
use CHOVEC_IO, only: ChoVec_Read, ChoVec_Size, NVTOT_CHOSYM
use PrintLevel, only: DEBUG
#ifndef _MOLCAS_MPP_
use fake_GA, only: GA_Arrays
#endif
use caspt2_global, only: iPrGlb
use caspt2_module, only: NAGEB, NAGEBES, NAGTB, NAGTBES, NIGEJ, NIGEJES, NIGTJ, NIGTJES, NSSH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Three, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVEC
integer(kind=iwp) :: IA, IAABS, IAEND, IAGEB, IAGEBTOT, IAGTB, IAGTBTOT, IAJ, IAL, IASTA, IC, ICABS, iCASE, ICJ, ICL, IDX, IIEND, &
                     IISTA, IJ, IJABS, IJGEL, IJGELTOT, IJGTL, IJGTLTOT, IL, ILABS, IOFFAJ, IOFFAL, IOFFCJ, IOFFCL, IOSYM(8,8), &
                     ISYA, ISYC, ISYJ, ISYL, ISYM, lg_W, MW, NAS, NCHOBUF, NIS, NV, NW
real(kind=wp) :: AJCL, ALCJ, HMACJL, HPACJL, SCL
real(kind=wp), allocatable :: CHOBUF(:)
real(kind=wp), parameter :: SQRT3 = sqrt(Three), SQRTH = sqrt(Half)
real(kind=wp), external :: DDot_
!logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

if (iPrGlb >= DEBUG) write(u6,*) 'RHS on demand: case H'

!***********************************************************************
! Case H:
!   WP(jl,ac)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
!   WM(jl,ac)=((ajcl)-(alcj))*SQRT(Three)
!***********************************************************************

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
call CHOVEC_SIZE(4,NCHOBUF,IOSYM)

call mma_allocate(CHOBUF,NCHOBUF,LABEL='CHOBUF')

call CHOVEC_READ(4,CHOBUF,NCHOBUF)

iCASE = 12
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
do ISYM=1,NSYM

  NAS = NAGEB(ISYM)
  NIS = NIGEJ(ISYM)
  NW = NAS*NIS

  if (NW == 0) cycle

  call RHS_ALLO(NAS,NIS,lg_W)
  call RHS_ACCESS(NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
  NW = NAS*(IIEND-IISTA+1)

  !*********************************************************************
  ! inner loop over RHS elements in symmetry ISYM
  !*********************************************************************
  do IJGEL=IISTA,IIEND
    IJGELTOT = IJGEL+NIGEJES(ISYM)
    IJABS = MIGEJ(1,IJGELTOT)
    ILABS = MIGEJ(2,IJGELTOT)
    IJ = MIREL(1,IJABS)
    ISYJ = MIREL(2,IJABS)
    IL = MIREL(1,ILABS)
    ISYL = MIREL(2,ILABS)
    do IAGEB=IASTA,IAEND ! these are always all elements
      IAGEBTOT = IAGEB+NAGEBES(ISYM)
      IAABS = MAGEB(1,IAGEBTOT)
      ICABS = MAGEB(2,IAGEBTOT)
      IA = MAREL(1,IAABS)
      ISYA = MAREL(2,IAABS)
      IC = MAREL(1,ICABS)
      ISYC = MAREL(2,ICABS)
      ! compute integrals (ajcl) and (alcj)
      NV = NVTOT_CHOSYM(Mul(ISYA,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
      IAJ = IA-1+NSSH(ISYA)*(IJ-1)
      ICL = IC-1+NSSH(ISYC)*(IL-1)
      IOFFAJ = 1+IOSYM(ISYA,ISYJ)+NV*IAJ
      IOFFCL = 1+IOSYM(ISYC,ISYL)+NV*ICL
      AJCL = DDOT_(NV,CHOBUF(IOFFAJ),1,CHOBUF(IOFFCL),1)

      NV = NVTOT_CHOSYM(Mul(ISYA,ISYL))
      IAL = IA-1+NSSH(ISYA)*(IL-1)
      ICJ = IC-1+NSSH(ISYC)*(IJ-1)
      IOFFAL = 1+IOSYM(ISYA,ISYL)+NV*IAL
      IOFFCJ = 1+IOSYM(ISYC,ISYJ)+NV*ICJ
      ALCJ = DDOT_(NV,CHOBUF(IOFFAL),1,CHOBUF(IOFFCJ),1)

      ! HP(ac,jl)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
      SCL = One
      if (IAABS == ICABS) SCL = SCL*SQRTH
      if (ILABS == IJABS) SCL = SCL*SQRTH
      HPACJL = SCL*(AJCL+ALCJ)
      ! write element HP(ac,jl)
      IDX = IAGEB+NAS*(IJGEL-IISTA)
#     ifdef _MOLCAS_MPP_
      DBL_MB(MW+IDX-1) = HPACJL
#     else
      GA_Arrays(lg_w)%A(IDX) = HPACJL
#     endif
    end do
  end do
  !*********************************************************************

  call RHS_RELEASE_UPDATE(lg_W,IASTA,IAEND,IISTA,IIEND)
  call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_W)
end do
!***********************************************************************

iCASE = 13
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
do ISYM=1,NSYM

  NAS = NAGTB(ISYM)
  NIS = NIGTJ(ISYM)
  NW = NAS*NIS

  if (NW == 0) cycle

  call RHS_ALLO(NAS,NIS,lg_W)
  call RHS_ACCESS(NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
  NW = NAS*(IIEND-IISTA+1)

  !*********************************************************************
  ! inner loop over RHS elements in symmetry ISYM
  !*********************************************************************
  do IJGTL=IISTA,IIEND
    IJGTLTOT = IJGTL+NIGTJES(ISYM)
    IJABS = MIGTJ(1,IJGTLTOT)
    ILABS = MIGTJ(2,IJGTLTOT)
    IJ = MIREL(1,IJABS)
    ISYJ = MIREL(2,IJABS)
    IL = MIREL(1,ILABS)
    ISYL = MIREL(2,ILABS)
    do IAGTB=IASTA,IAEND ! these are always all elements
      IAGTBTOT = IAGTB+NAGTBES(ISYM)
      IAABS = MAGTB(1,IAGTBTOT)
      ICABS = MAGTB(2,IAGTBTOT)
      IA = MAREL(1,IAABS)
      ISYA = MAREL(2,IAABS)
      IC = MAREL(1,ICABS)
      ISYC = MAREL(2,ICABS)
      ! compute integrals (ajcl) and (alcj)
      NV = NVTOT_CHOSYM(Mul(ISYA,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
      IAJ = IA-1+NSSH(ISYA)*(IJ-1)
      ICL = IC-1+NSSH(ISYC)*(IL-1)
      IOFFAJ = 1+IOSYM(ISYA,ISYJ)+NV*IAJ
      IOFFCL = 1+IOSYM(ISYC,ISYL)+NV*ICL
      AJCL = DDOT_(NV,CHOBUF(IOFFAJ),1,CHOBUF(IOFFCL),1)

      NV = NVTOT_CHOSYM(Mul(ISYA,ISYL))
      IAL = IA-1+NSSH(ISYA)*(IL-1)
      ICJ = IC-1+NSSH(ISYC)*(IJ-1)
      IOFFAL = 1+IOSYM(ISYA,ISYL)+NV*IAL
      IOFFCJ = 1+IOSYM(ISYC,ISYJ)+NV*ICJ
      ALCJ = DDOT_(NV,CHOBUF(IOFFAL),1,CHOBUF(IOFFCJ),1)

      ! HP(ac,jl)=((ajcl)-(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
      SCL = SQRT3
      HMACJL = SCL*(AJCL-ALCJ)
      ! write element HP(ac,jl)
      IDX = IAGTB+NAS*(IJGTL-IISTA)
#     ifdef _MOLCAS_MPP_
      DBL_MB(MW+IDX-1) = HMACJL
#     else
      GA_Arrays(lg_W)%A(IDX) = HMACJL
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

end subroutine RHSOD_H
