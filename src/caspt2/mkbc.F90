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
subroutine MKBC(DREF,NDREF,PREF,NPREF,FD,FP,NG3,F3,idxG3)

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
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF), FD(NDREF), FP(NPREF), F3(NG3)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) :: ICASE, IHI, ILO, ISYM, JHI, JLO, LDA, lg_BC, MBC, NAS, NBC, NIN
real(kind=wp) :: DBC
real(kind=wp), external :: PSBMAT_FPRINT
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: MA, MYRANK
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
  NBC = (NAS*(NAS+1))/2
  if (NBC <= 0) cycle

  ! Set up the matrix BC(tuv,xyz) defined by the expression
  ! <atuv|H0-E0|cxyz> = dac ( alpha(a) SC(tuv,xyz) + BC(tuv,xyz) )
  ! Formula used:
  !    BC(tuv,xyz)
  !    = Fvutxyz +dyu Fvztx + dyx Fvutz + dtu Fvxyz + dtu dyx Fvz
  !    +(Ey+Eu-EASUM)*SC(tuv,xyz)
  !    -Eu*( dyu Gvztx + dtu Gvxyz )
  !    -Ey dyx Gvutz
  !    -(Eu+Ey)*( dtu dyx Gvz )

  ! where dyu = Kronecker(y,u) etc. Gvutxyz=<Evutxyz>, etc.
  ! Similarly, Fvutxyz= Sum(w)(EPSA(w)<Evutxyzww>, etc.

  call PSBMAT_GETMEM('BC',lg_BC,NAS)
  call PSBMAT_READ('S',iCase,iSym,lg_BC,NAS)

  ! fill in the 3-el parts
# ifdef _MOLCAS_MPP_
  if (IS_REAL_PAR()) then
    MYRANK = GA_NODEID()
    call GA_DISTRIBUTION(LG_BC,MYRANK,ILO,IHI,JLO,JHI)
    if ((JLO /= 0) .and. (JHI-JLO+1 /= NAS)) then
      write(u6,*) 'MKBC: MISMATCH IN RANGE OF THE SUPERINDICES'
      call ABEND()
    end if
    if ((ILO > 0) .and. (JLO > 0)) then
      call GA_ACCESS(LG_BC,ILO,IHI,JLO,JHI,MA,LDA)
      MBC = LDA*(JHI-JLO+1)
      call MKBC_DP(DREF,NDREF,PREF,NPREF,FD,FP,iSYM,DBL_MB(MA),MBC,ILO,IHI,JLO,JHI,LDA)
      call MKBC_F3_MPP(ISYM,DBL_MB(MA),ILO,IHI,NAS,LDA,NG3,F3,IDXG3)
      call GA_RELEASE_UPDATE(LG_BC,ILO,IHI,JLO,JHI)
    else
      call MKBC_F3_MPP(ISYM,DUMMY,ILO,IHI,NAS,LDA,NG3,F3,IDXG3)
    end if
  else
# endif
    ILO = 1
    IHI = NAS
    JLO = 1
    JHI = NAS
    LDA = 0
    MBC = NAS*(NAS+1)/2
    call MKBC_DP(DREF,NDREF,PREF,NPREF,FD,FP,ISYM,GA_Arrays(lg_BC)%A(:),MBC,ILO,IHI,JLO,JHI,LDA)
    call MKBC_F3(ISYM,GA_Arrays(lg_BC)%A(:),MBC,NG3,F3,IDXG3)

# ifdef _MOLCAS_MPP_
  end if
# endif

  call PSBMAT_WRITE('B',iCase,iSYM,lg_BC,NAS)

  if (IPRGLB >= DEBUG) then
    DBC = PSBMAT_FPRINT(lg_BC,NAS)
    write(u6,'("DEBUG> ",A4,1X,I3,1X,ES21.14)') 'C',ISYM,DBC
  end if

  call PSBMAT_FREEMEM(lg_BC)
end do

end subroutine MKBC
