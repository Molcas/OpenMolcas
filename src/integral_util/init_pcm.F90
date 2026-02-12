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
! Copyright (C) 1992,2000, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Init_PCM(NonEq,iCharg)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!                                                                      *
!             Modified for Langevin polarizabilities, March 2000 (RL)  *
!***********************************************************************

use PCM_arrays, only: Centr, DCntr, dPnt, dRad, dTes, IntSph, MxVert, NewSph, nVert, PCM_n, PCM_SQ, PCMDm, PCMiSph, PCMSph, &
                      PCMTess, SSPh, Vert
use Isotopes, only: MaxAtomNum, PTab
use UnixInfo, only: ProgName
use rctfld_module, only: DoDeriv, iCharge_Ref, iSlPar, NoNEQ_Ref, nPCM_Info, nS, nTs, PCM, PCM_Info_Dmp
use Molcas, only: MxAtom
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
logical(kind=iwp), intent(in) :: NonEq
integer(kind=iwp), intent(inout) :: iCharg
character(len=2) ::Elements(MxAtom*8)
integer(kind=iwp) :: i, iPrint, j, lCnAtm, nAtoms
real(kind=wp), allocatable :: Coor(:,:), LcCoor(:,:)
integer(kind=iwp), allocatable :: ANr(:), LcANr(:)

if (.not. PCM) return
!                                                                      *
!***********************************************************************
!                                                                      *
! Reinitialize always for gradient calculations

DoDeriv = .false.
!pcm_solvent
! added mckinley for pcm in second derivatives
if ((ProgName == 'alaska') .or. (ProgName == 'mckinley') .or. (ProgName == 'mclr')) DoDeriv = .true.
!pcm_solvent end
if (DoDeriv) then
  LcNAtm = ISlPar(42)
  call mma_allocate(dTes,nTs,lcNAtm,3,Label='dTes')
  call mma_allocate(dPnt,nTs,lcNAtm,3,3,Label='dPnt')
  call mma_allocate(dRad,nS,lcNAtm,3,Label='dRad')
  call mma_allocate(dCntr,nS,lcNAtm,3,3,Label='dCntr')
  call mma_allocate(PCM_SQ,2,nTs,Label='PCM_SQ')
  call Get_dArray('PCM Charges',PCM_SQ,2*nTs)
else
  !                                                                      *
  !***********************************************************************
  !                                                                      *
  ! Check if we have retrievable PCM data

  call Get_iScalar('PCM info length',nPCM_info)
  if (nPCM_info /= 0) then

    ! nonequilibrium model for the case of ionization:
    ! redo PCM initiation but with the fake charge corresponding to
    ! the one of the ground (reference) state
    if ((iCharge_ref < iCharg) .and. NonEq) then
      iCharg = iCharge_ref
    else if ((iCharge_ref == iCharg) .and. (NonEq_ref .eqv. NonEq)) then
      ! If charge or NonEq/Eq status in retrieved data not the same as
      ! for request redo the PCM initiation.

      ! Evolving the new code

      call mma_allocate(PCMSph,4,NS,Label='PCMSph')
      call mma_allocate(PCMTess,4,nTs,Label='PCMTess')
      call mma_allocate(Vert,3,MxVert,nTs,Label='Vert')
      call mma_allocate(Centr,3,MxVert,nTs,Label='Centr')
      call mma_allocate(SSph,NS,Label='SSph')
      call mma_allocate(PCMDM,nTs,nTs,Label='PCMDM')
      call mma_allocate(PCM_N,NS,Label='PCM_N')
      call mma_allocate(PCMiSph,nTs,Label='PCMiSph')
      call mma_allocate(NVert,nTs,Label='NVert')
      call mma_allocate(IntSph,MxVert,nTs,Label='IntSph')
      call mma_allocate(NewSph,2,NS,Label='NewSph')

      call Get_dArray('PCMSph',PCMSph,4*NS)
      call Get_dArray('PCMTess',PCMTess,4*nTs)
      call Get_dArray('Vert',Vert,3*MxVert*nTs)
      call Get_dArray('Centr',Centr,3*MxVert*nTs)
      call Get_dArray('SSph',SSph,NS)
      call Get_dArray('PCMDM',PCMDM,nTs**2)
      call Get_iArray('PCM_N',PCM_N,NS)
      call Get_iArray('PCMiSph',PCMiSph,nTs)
      call Get_iArray('NVert',NVert,nTs)
      call Get_iArray('IntSph',IntSph,MxVert*nTs)
      call Get_iArray('NewSph',NewSph,2*NS)

      return
    end if
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Initial processing for PCM

call Get_nAtoms_All(nAtoms)
call mma_allocate(Coor,3,nAtoms,Label='Coor')
call Get_Coord_All(Coor,nAtoms)
call Get_Name_All(Elements)
call mma_Allocate(ANr,nAtoms,Label='ANr')
do i=1,nAtoms
  do j=0,MaxAtomNum
    if (PTab(j) == Elements(i)) then
      ANr(i) = j
      exit
    end if
  end do
end do
call mma_allocate(LcCoor,3,nAtoms,Label='LcCoor')
call mma_allocate(LcANr,nAtoms,Label='LcANr')

! Initialize PCM model

! iPrint: Print level
! ICharg: Molecular charge
! nAtoms: total number of atoms
! Coor: Coordinates of atoms
! MxVert*nTs ANr: atomic numbers
! LcCoor: local array for atomic coordinates
! LcANr: local array for atomic numbers
! Solvent: string with explicit solvent name
! Conductor: logical flag to activate conductor approximation
! aArea: average area of a tessera
! r_min_sphere: minimum radius of smoothing sphere
! ip_Ts: pointer to tesserae
! nTs  : number of tesserae

#ifdef _DEBUGPRINT_
iPrint = 99
#else
iPrint = 5
#endif
call PCM_Init(iPrint,ICharg,nAtoms,Coor,ANr,LcCoor,LcANr,NonEq)
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*)
#endif

call mma_deallocate(LcANr)
call mma_deallocate(LcCoor)
call mma_deallocate(ANr)
call mma_deallocate(Coor)
!                                                                      *
!***********************************************************************
!                                                                      *
! Put the dynamic arrays on COMFILE

call Put_iScalar('PCM info length',nPCM_info)

call Put_dArray('PCMSph',PCMSph,4*NS)
call Put_dArray('PCMTess',PCMTess,4*nTs)
call Put_dArray('Vert',Vert,3*MxVert*nTs)
call Put_dArray('Centr',Centr,3*MxVert*nTs)
call Put_dArray('SSph',SSph,NS)
call Put_dArray('PCMDM',PCMDM,nTs**2)
call Put_iArray('PCM_N',PCM_N,NS)
call Put_iArray('PCMiSph',PCMiSph,nTs)
call Put_iArray('NVert',NVert,nTs)
call Put_iArray('IntSph',IntSph,MxVert*nTs)
call Put_iArray('NewSph',NewSph,2*NS)

iCharge_ref = iCharg
NonEq_ref = NonEq

! Put the reaction field common blocks on disk again!

call PCM_Info_Dmp()

end subroutine Init_PCM
