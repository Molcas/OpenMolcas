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

subroutine CnstCLag(IFF,nLev,NG3,NCONF,CLag,DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G1,G2,G3)

use PrintLevel, only: VERBOSE
use caspt2_global, only: IDTCEX, iPrGlb, LUCIEX, LUSOLV, SGS
use caspt2_module, only: CITHR, EPSA, ETA, ISCF, JSTATE, MSTATE, NSTATE, STSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, byte, u6

implicit none
integer(kind=iwp), intent(in) :: IFF, nLev, NG3, NCONF
real(kind=wp), intent(inout) :: CLag(nConf), DG1(nLev**2), DG2(nLev**4), DG3(NG3), DF1(nLev**2), DF2(nLev**4), DF3(NG3), &
                                DEPSA(nLev**2)
real(kind=wp), intent(in) :: G1(nLev**2), G2(nLev**4), G3(NG3)
integer(kind=iwp) :: IDCI, ILEV, ILUID
real(kind=wp) :: CPE, CPTF0, CPTF10, CPUT, TIOE, TIOTF0, TIOTF10, WALLT
integer(kind=byte), allocatable :: idxG3(:,:)
real(kind=wp), allocatable :: CI1(:)

if (IFF == 1) then
  ! ORBITAL ENERGIES IN CI-COUPLING ORDER:
  do ILEV=1,NLEV
    ETA(ILEV) = EPSA(SGS%L2ACT(ILEV))
  end do
end if

!-SVC20100831: allocate local G3 matrices
call mma_allocate(idxG3,6,NG3,label='idxG3')
iLUID = 0
call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

call mma_allocate(CI1,NCONF,LABEL='CI')
if (ISCF == 0) then
  IDCI = IDTCEX(JSTATE)
  call DDAFILE(LUCIEX,2,CI1,NCONF,IDCI)
  if (IPRGLB >= VERBOSE) then
    write(u6,*)
    if (NSTATE > 1) then
      write(u6,'(A,I4)') ' With new orbitals, the CI array of state ',MSTATE(JSTATE)
    else
      write(u6,*) ' With new orbitals, the CI array is:'
    end if
    call PRWF_CP2(STSYM,NCONF,CI1,CITHR)
  end if
else
  CI1(1) = One
end if

call TIMING(CPTF0,CPE,TIOTF0,TIOE)
if (ISCF == 0) then
  call DERFG3(IFF,NCONF,NLEV,NG3,CI1,CLAG,DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G1,G2)
else
  call DERSPE(NLEV,NG3,DF1,DF2,DF3,idxG3,DEPSA,G1,G2,G3)
end if
call TIMING(CPTF10,CPE,TIOTF10,TIOE)
if (IPRGLB >= VERBOSE) then
  CPUT = CPTF10-CPTF0
  WALLT = TIOTF10-TIOTF0
  write(u6,*)
  write(u6,'(a,2f10.2)') ' DERFG3  : CPU/WALL TIME=',cput,wallt
end if

call mma_deallocate(CI1)
call mma_deallocate(idxG3)

return

end subroutine CnstCLag
