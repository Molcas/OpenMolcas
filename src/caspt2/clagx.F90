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

subroutine CLagX(IFF,nConf,nRoots,nState,nAshT,CLag,DEPSA,VECROT)

use PrintLevel, only: VERBOSE
use sguga, only: L2ACT, SGS
use caspt2_global, only: iPrGlb
use caspt2_module, only: EPSA, HZERO, ISCF, JSTATE, NG1, NG2, NG3, NG3TOT
use BDerNEV, only: BDerNEV_E4
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IFF, nConf, nRoots, nState, nAshT
real(kind=wp), intent(inout) :: CLag(nConf,nRoots), DEPSA(nAshT,nAshT)
real(kind=wp), intent(in) :: VECROT(nState)
integer(kind=iwp) :: ILEV, nLev
real(kind=wp) :: DEASUM
real(kind=wp) :: CPE, CPTF0, CPTF10, CPUT, TIOE, TIOTF0, TIOTF10, WALLT
real(kind=wp), allocatable :: G3(:)
real(kind=wp), allocatable, target :: DF1_H(:), DF2_H(:), DF3_H(:), DG1(:), DG2(:), DG3(:), F1_H(:), F2_H(:), G1(:), G2(:)
real(kind=wp), pointer :: DF1(:), DF2(:), DF3(:), F1(:), F2(:)

nLev = SGS%nLev

!! reduced density matrix and fock-weighted RDM
call mma_allocate(G1,NG1,Label='G1')
call mma_allocate(G2,NG2,Label='G2')
call mma_allocate(G3,NG3,Label='G3')

if (IFF == 1) then
  call mma_allocate(F1_H,NG1,Label='F1_H')
  call mma_allocate(F2_H,NG2,Label='F2_H')
  F1 => F1_H
  F2 => F2_H
else
  F1 => G1
  F2 => G2
end if

!! their derivative contributions
NG3tot = NG3
!! Use NG3tot (in caspt2_module) for the moment
#ifdef _MOLCAS_MPP_
if (is_real_par()) call gaigop_scal(ng3tot,'+')
#endif
call mma_allocate(DG1,NG1,Label='DG1')
call mma_allocate(DG2,NG2,Label='DG2')
call mma_allocate(DG3,NG3,Label='DG3')
call mma_allocate(DF1_H,NG1,Label='DF1_H')
call mma_allocate(DF2_H,NG2,Label='DF2_H')
if (IFF == 1) then
  call mma_allocate(DF3_H,NG3,Label='DF3_H')
else
  call mma_allocate(DF3_H,1,Label='DF3_H')
end if
DF1 => DF1_H
DF2 => DF2_H
DF3 => DF3_H

call PT2_GET(NG1,' GAMMA1',G1)
call PT2_GET(NG2,' GAMMA2',G2)
call PT2_GET(NG3,' GAMMA3',G3)

if (IFF == 1) then
  call PT2_GET(NG1,' DELTA1',F1)
  call PT2_GET(NG2,' DELTA2',F2)
end if

!! Initialize them
DG1(:) = Zero
DG2(:) = Zero
DG3(:) = Zero
DF1(:) = Zero
DF2(:) = Zero
DF3(:) = Zero
!! DEASUM is the derivative cont. of EASUM
DEASUM = Zero

call TIMING(CPTF0,CPE,TIOTF0,TIOE)
call CLagD(NASHT,NG3,NSTATE,G1,G2,G3,DG1,DG2,DG3,DF1,DF2,DF3,DEASUM,DEPSA,VECROT)
call TIMING(CPTF10,CPE,TIOTF10,TIOE)
if (IPRGLB >= VERBOSE) then
  CPUT = CPTF10-CPTF0
  WALLT = TIOTF10-TIOTF0
  write(u6,'(a,2f10.2)') ' CLagD   : CPU/WALL TIME=',cput,wallt
  !# ifdef _MOLCAS_MPP_
  !if (is_real_par()) call GADGOP_SCAL (deasum,'+')
  !# endif
  !write(u6,*) 'Deasum = ', deasum
  !# ifdef _MOLCAS_MPP_
  !if (is_real_par()) DEASUM = DEASUM/GA_NNODES()
  !# endif
end if

!! Some symmetrizations are likely required
call CLagSym(nAshT,DG1,DG2,DF1,DF2,0)

!! Do for the derivative of EASUM
!! EASUM=EASUM+EPSA(IT)*DREF(IT,IT)
if (IFF == 1) then
  do ILEV=1,nLev
    DG1(ILEV+nLev*(ILEV-1)) = DG1(ILEV+nLev*(ILEV-1))+DEASUM*EPSA(L2ACT(ILEV))
    if (ISCF == 0) then
      DEPSA(:,ILEV) = DEPSA(:,ILEV)+DEASUM*G1(nLev*(ILEV-1)+1:nLev*ILEV)
    else
      !! ?
    end if
  end do
end if

#ifdef _MOLCAS_MPP_
!! the master node does the job, so distribute to slave nodes
!! only for the G1 and G2 replicate arrays
if (is_real_par()) then
  call GADGOP(DG1,NG1,'+')
  call GADGOP(DG2,NG2,'+')
  if (IFF == 1) then
    call GADGOP(DF1,NG1,'+')
    call GADGOP(DF2,NG2,'+')
  end if
end if
#endif

call CnstCLag(IFF,nLev,NG3,NCONF,CLag(1,jState),DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G1,G2,G3)

call mma_deallocate(G1)
call mma_deallocate(G2)
call mma_deallocate(G3)
if (IFF == 1) then
  call mma_deallocate(F1_H)
  call mma_deallocate(F2_H)
end if

call mma_deallocate(DG1)
call mma_deallocate(DG2)
call mma_deallocate(DG3)
call mma_deallocate(DF1_H)
call mma_deallocate(DF2_H)
call mma_deallocate(DF3_H)

if (HZERO == 'DYALL') then
  if (IPRGLB >= VERBOSE) call TIMING(CPTF0,CPE,TIOTF0,TIOE)
  call BDerNEV_E4(nConf,nLev,CLag(1,jState))
  if (IPRGLB >= VERBOSE) then
    call TIMING(CPTF10,CPE,TIOTF10,TIOE)
    CPUT = CPTF10-CPTF0
    WALLT = TIOTF10-TIOTF0
    write(u6,'(a,2f10.2)') ' BDerNEV4: CPU/WALL TIME=',cput,wallt
  end if
end if

end subroutine CLagX
