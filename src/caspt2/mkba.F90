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
subroutine MKBA(DREF,NDREF,PREF,NPREF,FD,FP,NG3,F3,idxG3)

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
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF), F3(NG3)
real(kind=wp), intent(in) :: FD(NDREF), FP(NPREF)
integer(kind=Byte), intent(in) :: idxG3(6,NG3)
#ifdef _MOLCAS_MPP_
real(kind=wp) Dummy(1)
integer(kind=iwp) MYRANK, MA
#endif
integer(kind=iwp) ILO, IHI, JLO, JHI, LDA
integer(kind=iwp) ICASE, ISYM, NIN, NAS, NBA, lg_BA, MBA
real(kind=wp) DBA
real(kind=wp), external :: PSBMAT_FPRINT

ICASE = 1
! LONG loop over superindex symmetry.
do ISYM=1,NSYM
  NIN = NINDEP(ISYM,ICASE)
  if (NIN == 0) cycle
  NAS = NTUV(ISYM)
  NBA = (NAS*(NAS+1))/2
  if (NBA <= 0) cycle

  ! Set up the matrix BA(tuv,xyz) defined by the expression
  ! <ituv|H0-E0|kxyz> = dik ( alpha(a) SA(tuv,xyz) + BA(tuv,xyz) )
  ! Formula used:
  ! BA(tuv,xyz) = (Ey+Eu+Ex+Et-EASUM)*SA(tuv,xyz) - Fvuxtyz
  ! - dyu ( Fvzxt - Eu Gvzxt ) - dyt ( Fvuxz - Et Gvuxz )
  ! - dxu ( Fvtyz - Eu Gvtyz ) - dxu dyt ( Fvz - (Et+Eu) Gvz )
  ! + 2dxt ( Fvuyz - Et Gvuyz ) + 2dxt dyu ( Fvz - (Et+Eu) Gvz )

  ! where dyu = Kronecker(y,u) etc. Gvutxyz=<Evutxyz>, etc.
  ! Similarly, Fvutxyz= Sum(w)(EPSA(w)<Evutxyzww>, etc.

  call PSBMAT_GETMEM('BA',lg_BA,NAS)
  call PSBMAT_READ('S',iCase,iSym,lg_BA,NAS)

  ! fill in the 3-el parts
# ifdef _MOLCAS_MPP_
  if (IS_REAL_PAR()) then
    MYRANK = GA_NODEID()
    call GA_DISTRIBUTION(LG_BA,MYRANK,ILO,IHI,JLO,JHI)
    if ((JLO /= 0) .and. (JHI-JLO+1 /= NAS)) then
      write(u6,*) 'MKBA: MISMATCH IN RANGE OF THE SUPERINDICES'
      call ABEND()
    end if
    if ((ILO > 0) .and. (JLO > 0)) then
      call GA_ACCESS(LG_BA,ILO,IHI,JLO,JHI,MA,LDA)
      MBA = LDA*(JHI-JLO+1)
      call MKBA_DP(DREF,NDREF,PREF,NPREF,FD,FP,iSYM,DBL_MB(MA),MBA,ILO,IHI,JLO,JHI,LDA)
      call MKBA_F3_MPP(ISYM,DBL_MB(MA),ILO,IHI,NAS,LDA,NG3,F3,IDXG3)
      call GA_RELEASE_UPDATE(LG_BA,ILO,IHI,JLO,JHI)
    else
      call MKBA_F3_MPP(ISYM,DUMMY,ILO,IHI,NAS,LDA,NG3,F3,IDXG3)
    end if
  else
# endif
    LDA = 0
    ILO = 1
    IHI = NAS
    JLO = 1
    JHI = NAS
    MBA = NAS*(NAS+1)/2
    call MKBA_DP(DREF,NDREF,PREF,NPREF,FD,FP,ISYM,GA_Arrays(lg_BA)%A(:),MBA,ILO,IHI,JLO,JHI,LDA)
    call MKBA_F3(ISYM,GA_Arrays(lg_BA)%A(:),MBA,NG3,F3,IDXG3)
# ifdef _MOLCAS_MPP_
  end if
# endif

  call PSBMAT_WRITE('B',iCase,iSYM,lg_BA,NAS)

  if (IPRGLB >= DEBUG) then
    DBA = PSBMAT_FPRINT(lg_BA,NAS)
    write(u6,'("DEBUG> ",A4,1X,I3,1X,ES21.14)') 'A',ISYM,DBA
  end if

  call PSBMAT_FREEMEM(lg_BA)
end do

end subroutine MKBA
