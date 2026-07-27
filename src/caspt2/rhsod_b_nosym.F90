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

subroutine RHSOD_B_NOSYM(IVEC)

use Symmetry_Info, only: Mul
use SUPERINDEX, only: MIGEJ, MIGTJ, MIREL, MTGEU, MTGTU, MTREL
use CHOVEC_IO, only: ChoVec_Read, ChoVec_Size, NVTOT_ChoSym
use PrintLevel, only: DEBUG
#ifdef _MOLCAS_MPP_
use GA_Wrapper, only: DBL_MB
#else
use fake_GA, only: GA_Arrays
#endif
use caspt2_global, only: iPrGlb
use general_data, only:  nAsh
use caspt2_module, only: nASup, NIGEJES, NIGTJES, nISup, nSym, NTGEUES, NTGTUES
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iVec
integer(kind=iwp) :: IAEND, IASTA, iCASE, IDX, IIEND, IISTA, IJ, IJABS, IJGEL, IJGELTOT, IJGTL, IJGTLTOT, IL, ILABS, IOFFTJ, &
                     IOFFTL, IOFFVJ, IOFFVL, IOSYM(8,8), ISYJ, ISYL, ISYM, ISYT, ISYV, IT, ITABS, ITGEU, ITGEUTOT, ITGTU, &
                     ITGTUTOT, ITJ, ITL, IV, IVABS, IVJ, IVL, lg_W, MW, NAS, NCHOBUF, NIS, NV, NW
real(kind=wp) :: BMTVJL, BPTVJL, SCL, TJVL, TLVJ
real(kind=wp), allocatable :: CHOBUF(:)
real(kind=wp), parameter :: SQRTH = sqrt(Half)
real(kind=wp), external :: DDot_

if (iPrGlb >= DEBUG) write(u6,*) 'RHS on demand: case B'

!***********************************************************************
! Case B (2,3):
!   Let  W(tv,j,l)=(jt,lv):
!   BP(tv,jl)=((tj,vl)+(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
!   BM(tv,jl)=((tj,vl)-(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
!***********************************************************************

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
call CHOVEC_SIZE(1,NCHOBUF,IOSYM)

call mma_allocate(CHOBUF,NCHOBUF,LABEL='CHOBUF')

call CHOVEC_READ(1,CHOBUF,NCHOBUF)

iCASE = 2
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
  do IJGEL=IISTA,IIEND
    IJGELTOT = IJGEL+NIGEJES(ISYM)
    IJABS = MIGEJ(1,IJGELTOT)
    ILABS = MIGEJ(2,IJGELTOT)
    IJ = MIREL(1,IJABS)
    ISYJ = MIREL(2,IJABS)
    IL = MIREL(1,ILABS)
    ISYL = MIREL(2,ILABS)
    do ITGEU=IASTA,IAEND ! these are always all elements
      ITGEUTOT = ITGEU+NTGEUES(ISYM)
      ITABS = MTGEU(1,ITGEUTOT)
      IVABS = MTGEU(2,ITGEUTOT)
      IT = MTREL(1,ITABS)
      ISYT = MTREL(2,ITABS)
      IV = MTREL(1,IVABS)
      ISYV = MTREL(2,IVABS)
      ! compute integrals (ajcl) and (alcj)
      NV = NVTOT_CHOSYM(Mul(ISYT,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
      ITJ = IT-1+NASH(ISYT)*(IJ-1)
      IVL = IV-1+NASH(ISYV)*(IL-1)
      IOFFTJ = 1+IOSYM(ISYT,ISYJ)+NV*ITJ
      IOFFVL = 1+IOSYM(ISYV,ISYL)+NV*IVL
      TJVL = DDOT_(NV,CHOBUF(IOFFTJ),1,CHOBUF(IOFFVL),1)

      NV = NVTOT_CHOSYM(Mul(ISYT,ISYL))
      ITL = IT-1+NASH(ISYT)*(IL-1)
      IVJ = IV-1+NASH(ISYV)*(IJ-1)
      IOFFTL = 1+IOSYM(ISYT,ISYL)+NV*ITL
      IOFFVJ = 1+IOSYM(ISYV,ISYJ)+NV*IVJ
      TLVJ = DDOT_(NV,CHOBUF(IOFFTL),1,CHOBUF(IOFFVJ),1)

      ! BP(tv,jl)=((tj,vl)+(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
      SCL = Half
      if (ITABS == IVABS) SCL = SCL*Half
      if (ILABS == IJABS) SCL = SCL*SQRTH
      BPTVJL = SCL*(TJVL+TLVJ)
      ! write element HP(ac,jl)
      IDX = ITGEU+NAS*(IJGEL-IISTA)
#     ifdef _MOLCAS_MPP_
      DBL_MB(MW+IDX-1) = BPTVJL
#     else
      GA_Arrays(lg_w)%A(IDX) = BPTVJL
#     endif
    end do
  end do
  !*********************************************************************

  call RHS_RELEASE_UPDATE(lg_W,IASTA,IAEND,IISTA,IIEND)
  call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_W)
end do
!***********************************************************************

iCASE = 3
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
  do IJGTL=IISTA,IIEND
    IJGTLTOT = IJGTL+NIGTJES(ISYM)
    IJABS = MIGTJ(1,IJGTLTOT)
    ILABS = MIGTJ(2,IJGTLTOT)
    IJ = MIREL(1,IJABS)
    ISYJ = MIREL(2,IJABS)
    IL = MIREL(1,ILABS)
    ISYL = MIREL(2,ILABS)
    do ITGTU=IASTA,IAEND ! these are always all elements
      ITGTUTOT = ITGTU+NTGTUES(ISYM)
      ITABS = MTGTU(1,ITGTUTOT)
      IVABS = MTGTU(2,ITGTUTOT)
      IT = MTREL(1,ITABS)
      ISYT = MTREL(2,ITABS)
      IV = MTREL(1,IVABS)
      ISYV = MTREL(2,IVABS)
      ! compute integrals (tj,vl) and (tlvj)
      NV = NVTOT_CHOSYM(Mul(ISYT,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
      ITJ = IT-1+NASH(ISYT)*(IJ-1)
      IVL = IV-1+NASH(ISYV)*(IL-1)
      IOFFTJ = 1+IOSYM(ISYT,ISYJ)+NV*ITJ
      IOFFVL = 1+IOSYM(ISYV,ISYL)+NV*IVL
      TJVL = DDOT_(NV,CHOBUF(IOFFTJ),1,CHOBUF(IOFFVL),1)

      NV = NVTOT_CHOSYM(Mul(ISYT,ISYL))
      ITL = IT-1+NASH(ISYT)*(IL-1)
      IVJ = IV-1+NASH(ISYV)*(IJ-1)
      IOFFTL = 1+IOSYM(ISYT,ISYL)+NV*ITL
      IOFFVJ = 1+IOSYM(ISYV,ISYJ)+NV*IVJ
      TLVJ = DDOT_(NV,CHOBUF(IOFFTL),1,CHOBUF(IOFFVJ),1)

      ! BM(tv,jl)=((tj,vl)-(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
      SCL = Half
      !if (ITABS == IVABS) SCL = SCL*Half
      !if (ILABS == IJABS) SCL = SCL*SQRTH
      BMTVJL = SCL*(TJVL-TLVJ)
      ! write element BM(tv,jl)
      IDX = ITGTU+NAS*(IJGTL-IISTA)
#     ifdef _MOLCAS_MPP_
      DBL_MB(MW+IDX-1) = BMTVJL
#     else
      GA_Arrays(lg_w)%A(IDX) = BMTVJL
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

end subroutine RHSOD_B_NOSYM
