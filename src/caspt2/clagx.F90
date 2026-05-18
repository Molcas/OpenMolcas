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
use caspt2_global, only: iPrGlb
use Constants, only: Zero
use definitions, only: wp, iwp, u6
use stdalloc, only: mma_allocate, mma_deallocate
use sguga, only: SGS
use caspt2_module, only: NASH, ISCF, JSTATE, EPSA
use caspt2_module, only: NG1, NG2, NG3, NG3TOT
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif

implicit none
#ifdef _MOLCAS_MPP_
#include "global.fh"
#endif
integer(kind=iwp), intent(in) :: IFF, nConf, nRoots, nState, nAshT
real(kind=wp), intent(inout) :: CLag(nConf,nRoots), DEPSA(nAshT,nAshT)
real(kind=wp), intent(in) :: VECROT(nState)
real(kind=wp), allocatable :: G1(:), G2(:), G3(:)
real(kind=wp), allocatable :: DG1(:), DG2(:), DG3(:), DF1(:), DF2(:), DF3(:)
integer(kind=iwp) :: nLev, iT, iU
real(kind=wp) :: DEASUM
real(kind=wp) :: CPUT, WALLT, CPE, CPTF0, CPTF10, TIOE, TIOTF0, TIOTF10

nLev = SGS%nLev

!! reduced density matrix and fock-weighted RDM
call mma_allocate(G1,NG1,Label='G1')
call mma_allocate(G2,NG2,Label='G2')
call mma_allocate(G3,NG3,Label='G3')

!! their derivative contributions
NG3tot = NG3
!! Use NG3tot (in caspt2_module) for the moment
#ifdef _MOLCAS_MPP_
if (is_real_par()) call gaigop_scal(ng3tot,'+')
#endif
call mma_allocate(DG1,NG1,Label='DG1')
call mma_allocate(DG2,NG2,Label='DG2')
call mma_allocate(DG3,NG3,Label='DG3')
call mma_allocate(DF1,NG1,Label='DF1')
call mma_allocate(DF2,NG2,Label='DF2')
call mma_allocate(DF3,NG3,Label='DF3')

call PT2_GET(NG1,' GAMMA1',G1)
call PT2_GET(NG2,' GAMMA2',G2)
call PT2_GET(NG3,' GAMMA3',G3)

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
do iT=1,nAsh(1)
  DG1(iT+nAsh(1)*(iT-1)) = DG1(iT+nAsh(1)*(iT-1))+DEASUM*EPSA(iT)
  if (ISCF == 0) then
    do iU=1,nAsh(1)
      DEPSA(iT,iU) = DEPSA(iT,iU)+DEASUM*G1(iT+nAsh(1)*(iU-1))
    end do
  else
    !! ?
  end if
end do

#ifdef _MOLCAS_MPP_
!! the master node does the job, so distribute to slave nodes
!! only for the G1 and G2 replicate arrays
if (is_real_par()) then
  call GADGOP(DG1,NG1,'+')
  call GADGOP(DG2,NG2,'+')
  call GADGOP(DF1,NG1,'+')
  call GADGOP(DF2,NG2,'+')
end if
#endif

call CnstCLag(IFF,nLev,NG3,NCONF,CLag(1,jState),DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G1,G2,G3)

call mma_deallocate(G1)
call mma_deallocate(G2)
call mma_deallocate(G3)

call mma_deallocate(DG1)
call mma_deallocate(DG2)
call mma_deallocate(DG3)
call mma_deallocate(DF1)
call mma_deallocate(DF2)
call mma_deallocate(DF3)

end subroutine CLagX
