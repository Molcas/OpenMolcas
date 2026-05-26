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
! Case C (ICASE=4)
!***********************************************************************
subroutine MKSC(DREF,NDREF,PREF,NPREF,NG3,G3,idxG3)

use Index_Functions, only: nTri_Elem
use PrintLevel, only: DEBUG
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_GA, only: GA_Arrays
use caspt2_global, only: iPrGlb
use caspt2_module, only: NINDEP, NSYM, NTUV
use Definitions, only: wp, iwp, u6, byte

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF, NG3
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
real(kind=wp), intent(inout) :: G3(NG3)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) :: ICASE, IHI, ILO, ISYM, JHI, JLO, LDC, lg_SC, NAS, NIN, NSC
real(kind=wp) :: DSC
real(kind=wp), external :: PSBMAT_FPRINT
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: MC, MSC, MYRANK
real(kind=wp) :: Dummy(1)
#include "global.fh"
#include "mafdecls.fh"
#endif

ICASE = 4
! LONG loop over superindex symmetry.
do ISYM=1,NSYM
  NIN = NINDEP(ISYM,ICASE)
  if (NIN == 0) cycle
  NAS = NTUV(ISYM)
  NSC = nTri_Elem(NAS)
  if (NSC <= 0) cycle

  ! Set up the matrix SC(tuv,xyz) defined by the expression
  ! <atuv|cxyz> = dac SC(tuv,xyz)
  ! Formula used:
  !    SC(tuv,xyz)
  !    = Gvutxyz +dyu Gvztx + dyx Gvutz + dtu Gvxyz + dtu dyx Gvz

  call PSBMAT_GETMEM('SC',lg_SC,NAS)

  ! fill in the 3-el parts
# ifdef _MOLCAS_MPP_
  if (IS_REAL_PAR()) then
    MYRANK = GA_NODEID()
    call GA_DISTRIBUTION(LG_SC,MYRANK,ILO,IHI,JLO,JHI)
    if ((JLO /= 0) .and. (JHI-JLO+1 /= NAS)) then
      write(u6,*) 'MKSC: MISMATCH IN RANGE OF THE SUPERINDICES'
      call ABEND()
    end if
    if ((ILO > 0) .and. (JLO > 0)) then
      call GA_ACCESS(LG_SC,ILO,IHI,JLO,JHI,MC,LDC)
      call MKSC_G3_MPP(ISYM,DBL_MB(MC),ILO,NAS,LDC,NG3,G3,IDXG3)
      MSC = LDC*(jHi-jLo+1)
      call MKSC_DP(DREF,NDREF,PREF,NPREF,ISYM,DBL_MB(MC),MSC,ILO,IHI,JLO,JHI,LDC)
      call GA_RELEASE_UPDATE(LG_SC,ILO,IHI,JLO,JHI)
    else
      call MKSC_G3_MPP(ISYM,DUMMY,ILO,NAS,LDC,NG3,G3,IDXG3)
    end if
  else
# endif
    iLo = 1
    iHi = NAS
    jLo = 1
    jHi = NAS
    LDC = 0
    call MKSC_G3(ISYM,GA_Arrays(lg_SC)%A(:),NSC,NG3,G3,IDXG3)
    call MKSC_DP(DREF,NDREF,PREF,NPREF,ISYM,GA_Arrays(lg_SC)%A(:),NSC,ILO,IHI,JLO,JHI,LDC)
# ifdef _MOLCAS_MPP_
  end if
# endif

  call PSBMAT_WRITE('S',iCase,iSYM,lg_SC,NAS)

  if (IPRGLB >= DEBUG) then
    DSC = PSBMAT_FPRINT(lg_SC,NAS)
    write(u6,'("DEBUG> ",A4,1X,I3,1X,ES21.14)') 'C',ISYM,DSC
  end if

  call PSBMAT_FREEMEM(lg_SC)
end do

end subroutine MKSC
