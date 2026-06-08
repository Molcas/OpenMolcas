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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine MS_STRANS(IVEC,JVEC,NASHT,NTG3,OVL,TG1,TG2,TG3,HEL,SCAL_)
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

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use GA_Wrapper, only: DBL_MB
#endif
use PrintLevel, only: DEBUG
use fake_GA, only: GA_Arrays
use caspt2_global, only: iPrGlb
use caspt2_module, only: CASES, NASUP, NINDEP, NISUP, NSYM
use SC_NEVPT2, only: SC_prop
use Constants, only: Half, Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC, NASHT, NTG3
real(kind=wp), intent(inout) :: OVL, TG1(NASHT,NASHT), TG2(NASHT,NASHT,NASHT,NASHT), TG3(NTG3)
real(kind=wp), intent(out) :: HEL
real(kind=wp), intent(in) :: SCAL_
integer(kind=iwp) :: IC, ICASE, iHi1, iHi2, iLo1, iLo2, IS, ISYM, jHi1, jHi2, jLo1, jLo2, lg_V1, lg_V2, MV1, MV2, NAS, NHECOMP, &
                     NIN, NIS, nvlen
real(kind=wp) :: HEBLK, HECOMP(14,9)
character(len=1) :: CNT
real(kind=wp) :: SCAL

! The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

! Sketch of procedure:
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

CNT = 'N'
SCAL = SCAL_
if (SC_prop) then
  CNT = 'T'
  SCAL = Half*SCAL_
end if

HEL = Zero
HECOMP = Zero
do ICASE=1,13
  !if ((icase /= 12) .and. (icase /= 13)) cycle ! H
  !if ((icase /= 10) .and. (icase /= 11)) cycle ! G
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIN = NINDEP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    HEBLK = Zero

    if ((NAS*NIS /= 0) .and. (NIN /= 0)) then
      call RHS_ALLO(NAS,NIS,lg_V1)
      call RHS_ALLO(NAS,NIS,lg_V2)
      call RHS_READ(NAS,NIS,lg_V1,ICASE,ISYM,IVEC)
      call RHS_READ(NAS,NIS,lg_V2,ICASE,ISYM,JVEC)
      call RHS_ACCESS(NAS,NIS,lg_V1,iLo1,iHi1,jLo1,jHi1,MV1)
      call RHS_ACCESS(NAS,NIS,lg_V2,iLo2,iHi2,jLo2,jHi2,MV2)

      if ((iLo1 /= iLo2) .or. (iHi1 /= iHi2) .or. (jLo1 /= jLo2) .or. (jHi1 /= jHi2)) then
        write(u6,'(1X,A)') 'HCOUP: Error: block mismatch, abort...'
        call ABEND()
      end if

      nvlen = (iHi1-jLo1+1)*(jHi1-jLo1+1)

#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        call MS_STRANS_BLK(ICASE,ISYM,NAS,nvlen,NASHT,NTG3,jLo1,jHi1,DBL_MB(MV1),DBL_MB(MV2),OVL,TG1,TG2,TG3,SCAL,CNT)
      else
#     endif
        call MS_STRANS_BLK(ICASE,ISYM,NAS,nvlen,NASHT,NTG3,jLo1,jHi1,GA_Arrays(MV1)%A,GA_Arrays(MV2)%A,OVL,TG1,TG2,TG3,SCAL,CNT)
#     ifdef _MOLCAS_MPP_
      end if
#     endif
      !! Save T*S
      call RHS_SAVE(NAS,NIS,lg_V2,ICASE,ISYM,JVEC)

      call RHS_RELEASE(lg_V1,iLo1,iHi1,jLo1,jHi1)
      call RHS_RELEASE(lg_V2,iLo2,iHi2,jLo2,jHi2)
      call RHS_FREE(lg_V1)
      call RHS_FREE(lg_V2)
    end if
    HECOMP(ICASE,ISYM) = HEBLK
    HEL = HEL+HEBLK
  end do
end do

! Sum-reduce the per-process contributions
call GADGOP_SCAL(HEL,'+')
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
  call XFLUSH(u6)
  write(u6,'(1X,A8,9F12.8)') 'Summed: ',(HECOMP(14,IS),IS=1,NSYM+1)
  write(u6,*)
end if

end subroutine MS_STRANS
