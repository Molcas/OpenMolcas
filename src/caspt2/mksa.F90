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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

!***********************************************************************
! Case A (ICASE=1)
!***********************************************************************
subroutine MKSA(DREF,NDREF,PREF,NPREF,NG3,G3,idxG3)

use definitions, only: iwp, wp, u6, Byte
use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_GA, only: GA_Arrays
use caspt2_module, only: NSYM, NINDEP, NTUV

implicit none
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
integer(kind=iwp), intent(in) :: NDREF, NPREF, NG3
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
real(kind=wp), intent(in) :: G3(NG3)
integer(kind=Byte), intent(in) :: idxG3(6,NG3)
#ifdef _MOLCAS_MPP_
real(kind=wp) Dummy(1)
integer(kind=iwp) MYRANK, MA
#endif
integer(kind=iwp) ILO, IHI, JLO, JHI, LDA
integer(kind=iwp) ICASE, ISYM, lg_SA, NAS, NIN, NSA, MSA
real(kind=wp), external :: PSBMAT_FPRINT
real(kind=wp) DSA

ICASE = 1
! LONG loop over superindex symmetry.
do ISYM=1,NSYM
  NIN = NINDEP(ISYM,ICASE)
  if (NIN == 0) cycle
  NAS = NTUV(ISYM)
  NSA = (NAS*(NAS+1))/2
  if (NSA <= 0) cycle
  ! Set up the matrix SA(tuv,xyz) defined by the expression
  ! <ituv|kxyz> = dik SA(tuv,xyz)
  ! Formula used:
  !    SA(tuv,xyz) =  -Gvuxtyz -dyu Gvzxt - dyt Gvuxz - dxu Gvtyz - dxu dyt Gvz +2 dtx Gvuyz + 2 dtx dyu Gvz

  call PSBMAT_GETMEM('SA',lg_SA,NAS)

  ! fill in the 3-el parts
# ifdef _MOLCAS_MPP_
  if (IS_REAL_PAR()) then
    MYRANK = GA_NODEID()
    call GA_DISTRIBUTION(LG_SA,MYRANK,ILO,IHI,JLO,JHI)
    if ((JLO /= 0) .and. (JHI-JLO+1 /= NAS)) then
      write(u6,*) 'MKSA: MISMATCH IN RANGE OF THE SUPERINDICES'
      call ABEND()
    end if
    if ((ILO > 0) .and. (JLO > 0)) then
      call GA_ACCESS(LG_SA,ILO,IHI,JLO,JHI,MA,LDA)
      call MKSA_G3_MPP(ISYM,DBL_MB(MA),ILO,IHI,NAS,LDA,NG3,G3,IDXG3)
      MSA = LDA*(jHi-jLo+1)
      call MKSA_DP(DREF,NDREF,PREF,NPREF,ISYM,DBL_MB(MA),MSA,ILO,IHI,JLO,JHI,LDA)
      call GA_RELEASE_UPDATE(LG_SA,ILO,IHI,JLO,JHI)
    else
      call MKSA_G3_MPP(ISYM,DUMMY,ILO,IHI,NAS,LDA,NG3,G3,IDXG3)
    end if
  else
# endif
    iLo = 1
    iHi = NAS
    jLo = 1
    jHi = NAS
    LDA = 0
    MSA = NAS*(NAS+1)/2
    call MKSA_G3(ISYM,GA_Arrays(lg_SA)%A(:),MSA,NG3,G3,IDXG3)
    call MKSA_DP(DREF,NDREF,PREF,NPREF,ISYM,GA_Arrays(lg_SA)%A(:),MSA,ILO,IHI,JLO,JHI,LDA)
# ifdef _MOLCAS_MPP_
  end if
# endif

  call PSBMAT_WRITE('S',iCase,iSYM,lg_SA,NAS)

  if (IPRGLB >= DEBUG) then
    DSA = PSBMAT_FPRINT(lg_SA,NAS)
    write(u6,'("DEBUG> ",A4,1X,I3,1X,ES21.14)') 'A',ISYM,DSA
  end if

  call PSBMAT_FREEMEM(lg_SA)
end do

end subroutine MKSA
