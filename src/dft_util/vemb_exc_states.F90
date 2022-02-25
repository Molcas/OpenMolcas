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

subroutine VEMB_Exc_states(Vemb,nVemb,xKSDFT,Func_Bx)

implicit real*8(a-h,o-z)
real*8 Vemb(nVemb)
real*8 Func_Bx
character*(*) xKSDFT
character*16 MyNamRfil
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "ciinfo.fh"
#include "rctfld.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
real*8, allocatable :: D1ao_b(:), F_DFT(:)
real*8, allocatable :: xxCMO(:), xxOCCN(:), DState(:)
real*8 :: Dummy(1) = [0.0d0]
integer :: nDummy = 1

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
  DSTATE(:) = 0.0d0

  call DONE_RASSCF(xxCMO,xxOCCN,DState) ! computes D=CnC'
  ! Nonelectr. Vemb with GS and excited state density
  Vemb_Xstate = ddot_(nVemb,Vemb,1,DState,1)
  !write(6,*) 'Kroot, Vemb_K ',KROOT,Vemb_Xstate
  write(6,'(A,F19.10,3X,A,I3)') 'Nonelectr. Vemb w. rhoA_emb =',Vemb_Xstate,'root = ',KROOT
  ! E_xc,T[rhoA]
  Func_A = 0.0d0
  F_DFT(:) = 0.0d0
  call dscal_(nVemb,0.5d0,DState,1)
  call wrap_DrvNQ(xKSDFT,F_DFT,1,Func_A,DState,nVemb,1,.false.,Dummy,nDummy,'SCF ')
  !write(6,*) 'Kroot, Func_A ',KROOT,Func_A
  ! E_xc,T[rhoA+rhoB]
  call Get_NameRun(MyNamRfil) ! save current Runfile name
  call NameRun('AUXRFIL') ! switch RUNFILE name
  call Get_D1ao(D1ao_b,nVemb)
  call daxpy_(nVemb,0.5d0,D1ao_b,1,DState,1)

  Func_AB = 0.0d0
  F_DFT(:) = 0.0d0
  call wrap_DrvNQ(xKSDFT,F_DFT,1,Func_AB,DState,nVemb,1,.false.,Dummy,nDummy,'SCF ')
  !write(6,*) 'Kroot, Func_AB',KROOT,Func_AB
  !write(6,*) 'Kroot, Func_Bx',KROOT,Func_Bx
  ! Calculate DFT NAD for all densities:
  DFT_NAD = Func_AB-Func_A-Func_Bx
  write(6,'(A,F19.10,3X,A,I3)') 'DFT energy (NAD) =           ',DFT_NAD,'root = ',KROOT
  call NameRun(MyNamRfil) ! go back to MyNamRfil
end do
call mma_deallocate(D1ao_b)
call mma_deallocate(F_DFT)
call mma_deallocate(DState)
call mma_deallocate(xxCMO)
call mma_deallocate(xxOCCN)

return

end subroutine VEMB_Exc_states
