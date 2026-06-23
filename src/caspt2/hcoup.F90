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
! Copyright (C) 2014, Steven Vancoillie                                *
!***********************************************************************

subroutine HCOUP(IVEC,JVEC,OVL,TG1,TG2,NASHT,TG3,NTG3,HEL)
! Compute the coupling Hamiltonian element defined as
!     HEL = < ROOT1 | H * OMEGA | ROOT2 >
! assuming that IVEC contains a contravariant representation of
! H|ROOT1>, JVEC contains a contravariant representation of
! OMEGA|ROOT2>, and OVL, TG1, TG2, TG3 contain the overlap (normally
! expected to be 0 or 1) and active transition density matrices of ROOT1
! and ROOT2. See also subroutine TSVEC for explanations.

! SVC (March 2014): modification of original code to handle distributed
! RHS arrays. There is now a main HCOUP subroutine that loops over cases
! and irreps and gets access to the process-specific block of the RHS.
! The coupling for that block is computed by the subroutine HCOUP_BLK.

use PrintLevel, only: DEBUG
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use GA_Wrapper, only: DBL_MB
#endif
use fake_GA, only: GA_Arrays
use caspt2_global, only: iPrGlb
use caspt2_module, only: CASES, HZERO, NASUP, NINDEP, NISUP, NSYM
use SC_NEVPT2, only: Do_FIC, Do_SC, DVALUE_SC, SC_NEVPT2_amplitude
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC, NASHT, NTG3
real(kind=wp), intent(in) :: OVL, TG1(NASHT,NASHT), TG2(NASHT,NASHT,NASHT,NASHT), TG3(NTG3)
real(kind=wp), intent(out) :: HEL
integer(kind=iwp) :: IAEND1, IAEND2, IASTA1, IASTA2, IC, ICASE, iHi1, iHi2, IIEND1, IIEND2, IISTA1, IISTA2, iLo1, iLo2, IS, ISYM, &
                     jHi1, jHi2, jLo1, jLo2, lg_V1, lg_V2, MV1, MV2, NAS, NHECOMP, NIN, NIS, NV1, NV2
real(kind=wp) :: HEBLK, HEBLK_SC, HECOMP(14,9)

! The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

! Sketch of procedure:
!  HEL=Zero
!  Loop over every (case/symmetry)-block.
!     If (No such vector block) Skip to end of loop
!     Allocate two places for this block, VEC1 and VEC2
!     Read VEC1 as IVEC component from file.
!     Read VEC2 as JVEC component from file.
!     Loop nest, computing
!        HEL := HEL + VEC1*GOM*VEC2
!     End of loop nest
!     Deallocate VEC1 and VEC2
!  End of loop.

HEL = Zero
HECOMP(:,:) = Zero
if ((HZERO == 'DYALL') .and. Do_SC) DVALUE_SC = Zero
do ICASE=1,13
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIN = NINDEP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    HEBLK = Zero

    if (NAS*NIS*NIN /= 0) then

      call RHS_ALLO(NAS,NIS,lg_V1)
      call RHS_ALLO(NAS,NIS,lg_V2)
      call RHS_READ(NAS,NIS,lg_V1,ICASE,ISYM,IVEC)
      call RHS_ACCESS(NAS,NIS,lg_V1,iLo1,iHi1,jLo1,jHi1,MV1)
      NV1 = (iHi1-iLo1+1)*(jHi1-jLo1+1)

      if (Do_FIC) then
        call RHS_READ(NAS,NIS,lg_V2,ICASE,ISYM,JVEC)
        call RHS_ACCESS(NAS,NIS,lg_V2,iLo2,iHi2,jLo2,jHi2,MV2)
        NV2 = (iHi2-iLo2+1)*(jHi2-jLo2+1)

        if ((iLo1 /= iLo2) .or. (iHi1 /= iHi2) .or. (jLo1 /= jLo2) .or. (jHi1 /= jHi2)) then
          write(u6,'(1X,A)') 'HCOUP: Error: block mismatch, abort...'
          call ABEND()
        end if

#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          call HCOUP_BLK(ICASE,ISYM,NAS,jLo1,jHi1,DBL_MB(MV1),NV1,DBL_MB(MV2),NV2,OVL,HEBLK,TG1,TG2,NASHT,TG3,NTG3)
        else
#       endif
          call HCOUP_BLK(ICASE,ISYM,NAS,jLo1,jHi1,GA_Arrays(MV1)%A,NV1,GA_Arrays(MV2)%A,NV2,OVL,HEBLK,TG1,TG2,NASHT,TG3,NTG3)
#       ifdef _MOLCAS_MPP_
        end if
#       endif

        call RHS_RELEASE(lg_V2,IASTA2,IAEND2,IISTA2,IIEND2)
      end if

      if ((HZERO == 'DYALL') .and. Do_SC) then
        HEBLK_SC = Zero
        call RHS_READ(NAS,NIS,lg_V2,ICASE,ISYM,IVEC)
        call RHS_ACCESS(NAS,NIS,lg_V2,iLo2,iHi2,jLo2,jHi2,MV2)
        NV2 = (iHi2-iLo2+1)*(jHi2-jLo2+1)
        call SC_NEVPT2_amplitude(NAS,NIS,ICASE,ISYM,lg_V2)

#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          call HCOUP_BLK(ICASE,ISYM,NAS,jLo1,jHi1,DBL_MB(MV1),NV1,DBL_MB(MV2),NV2,OVL,HEBLK_SC,TG1,TG2,NASHT,TG3,NTG3)
        else
#       endif
          call HCOUP_BLK(ICASE,ISYM,NAS,jLo1,jHi1,GA_Arrays(MV1)%A,NV1,GA_Arrays(MV2)%A,NV2,OVL,HEBLK_SC,TG1,TG2,NASHT,TG3,NTG3)
#       ifdef _MOLCAS_MPP_
        end if
#       endif

        DVALUE_SC = DVALUE_SC+HEBLK_SC
        call RHS_RELEASE(lg_V2,IASTA2,IAEND2,IISTA2,IIEND2)
      end if

      call RHS_RELEASE(lg_V1,IASTA1,IAEND1,IISTA1,IIEND1)
      call RHS_FREE(lg_V1)
      call RHS_FREE(lg_V2)

    end if

    HECOMP(ICASE,ISYM) = HEBLK
    HEL = HEL+HEBLK
  end do
end do

! Sum-reduce the per-process contributions
call GADGOP_SCAL(HEL,'+')
if ((HZERO == 'DYALL') .and. Do_SC) call GADGOP_SCAL(DVALUE_SC,'+')
NHECOMP = 14*9
call GADGOP(HECOMP,NHECOMP,'+')

if (IPRGLB >= DEBUG) then
  do ICASE=1,13
    HECOMP(ICASE,NSYM+1) = sum(HECOMP(ICASE,1:NSYM))
  end do

  do ISYM=1,NSYM+1
    HECOMP(14,ISYM) = sum(HECOMP(1:13,ISYM))
  end do

  write(u6,'(A)') repeat('-',80)
  write(u6,*) 'HCOUP: The contributions to the Hamiltonian coupling'
  write(u6,*) ' elements, by case and by symmetry label.'
  do IC=1,13
    write(u6,'(1X,A8,9F12.8)') CASES(IC),(HECOMP(IC,IS),IS=1,NSYM+1)
  end do
  write(u6,'(1X,A8,9F12.8)') 'Summed: ',(HECOMP(14,IS),IS=1,NSYM+1)
  write(u6,*)
end if

end subroutine HCOUP
