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
! Copyright (C) 2010,2012,2017, Francesco Aquilante                    *
!               2015,2017, Alexander Zech                              *
!***********************************************************************

#include "compiler_features.h"
#ifdef _NOT_USED_

subroutine VEMB_Exc_states(Vemb,nVemb,xKSDFT,Func_Bx)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nVemb
real(kind=wp), intent(inout) :: Vemb(nVemb)
character(len=*), intent(in) :: xKSDFT
real(kind=wp), intent(in) :: Func_Bx
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
integer(kind=iwp) :: IAD12, KROOT, nDummy
real(kind=wp) :: DFT_NAD, Dummy(1), Func_A, Func_AB, Vemb_Xstate
real(kind=wp), allocatable :: D1ao_b(:), DState(:), F_DFT(:), xxCMO(:), xxOCCN(:)
real(kind=wp), external :: ddot_

nDummy = 1
Dummy(1) = Zero

IAD12 = IADR15(12)

call mma_allocate(xxCMO,NTOT2,Label='xxCMO')
call mma_allocate(xxOCCN,NTOT,Label='xxOCCN')
call mma_allocate(DState,NTOT1,Label='DState')
call mma_allocate(F_DFT,nVemb,Label='F_DFT')
call mma_allocate(D1ao_b,nVemb,Label='D1ao_b')

do KROOT=1,LROOTS

  ! Read natural orbitals
  if (NAC > 0) then
    call DDAFILE(JOBIPH,2,xxCMO,NTOT2,IAD12)
    call DDAFILE(JOBIPH,2,xxOCCN,NTOT,IAD12)
  end if
  ! Get GS and excited state densities:
  ! Fill allocated mem with zeroes.
  DSTATE(:) = Zero

  call DONE_RASSCF(xxCMO,xxOCCN,DState) ! computes D=CnC'
  ! Nonelectr. Vemb with GS and excited state density
  Vemb_Xstate = ddot_(nVemb,Vemb,1,DState,1)
  !write(u6,*) 'Kroot, Vemb_K ',KROOT,Vemb_Xstate
  write(u6,'(A,F19.10,3X,A,I3)') 'Nonelectr. Vemb w. rhoA_emb =',Vemb_Xstate,'root = ',KROOT
  ! E_xc,T[rhoA]
  Func_A = Zero
  F_DFT(:) = Zero
  DState(1:nVemb) = Half*DState(1:nVemb)
  call wrap_DrvNQ(xKSDFT,F_DFT,1,Func_A,DState,nVemb,1,.false.,Dummy,nDummy,'SCF ')
  !write(u6,*) 'Kroot, Func_A ',KROOT,Func_A
  ! E_xc,T[rhoA+rhoB]
  call NameRun('AUXRFIL') ! switch RUNFILE name
  call Get_dArray_chk('D1ao',D1ao_b,nVemb)
  DState(1:nVemb) = DState(1:nVemb)+Half*D1ao_b(:)

  Func_AB = Zero
  F_DFT(:) = Zero
  call wrap_DrvNQ(xKSDFT,F_DFT,1,Func_AB,DState,nVemb,1,.false.,Dummy,nDummy,'SCF ')
  !write(u6,*) 'Kroot, Func_AB',KROOT,Func_AB
  !write(u6,*) 'Kroot, Func_Bx',KROOT,Func_Bx
  ! Calculate DFT NAD for all densities:
  DFT_NAD = Func_AB-Func_A-Func_Bx
  write(u6,'(A,F19.10,3X,A,I3)') 'DFT energy (NAD) =           ',DFT_NAD,'root = ',KROOT
  call NameRun('#Pop')    ! go back to previous name
end do
call mma_deallocate(D1ao_b)
call mma_deallocate(F_DFT)
call mma_deallocate(DState)
call mma_deallocate(xxCMO)
call mma_deallocate(xxOCCN)

return

end subroutine VEMB_Exc_states

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(VEMB_Exc_states)

#endif
