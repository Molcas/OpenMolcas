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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
!                                                                      *
! 2023, Matthew R. Hennefarth - modified lhrot -> si_pdft              *
! 2024, Matthew R. Henenfarth - upgraded to F90                        *
!***********************************************************************

subroutine MSPDFTGrad_Misc(si_pdft,states)
! This subroutine does miscellaneous things needed
! in MS-PDFT gradient calculation.

use mspdftgrad, only: D1aoMS, D1SAOMS, DIDA, F1MS, F2MS, FocMS, FxyMS, P2MOT
use rasscf_global, only: lroots, NACPR2, nTot4
use general_data, only: ispin, ntot1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: si_pdft(lroots,lroots)
integer(kind=iwp), intent(in) :: states(2)
real(kind=wp), allocatable :: r_active_dm1(:), r_casdm2(:), r_fock_occ(:), r_Fxy(:), u_rlx_sq(:)
! r_* : rotated (ie, in final MS-PDFT) eigenbasis quantities (i think...)

! Functions added by Paul Calio for MECI Opt *****
! Original calls are in slapaf_util/start_alasaks.f

call put_iArray('NACstatesOpt',states,2)
! End of stuff added by Paul

call Put_DArray('MS_FINAL_ROT',si_pdft,lRoots**2)
call Put_DArray('F1_PDFT',F1MS,nTot1*lRoots)
call Put_DArray('F2_PDFT',F2MS,NACPR2*lRoots)
call Put_DArray('D1AO_MS',D1AOMS,nTot1*lRoots)
if (ispin /= 1) call Put_DArray('D1SAO_MS',D1SAOMS,nTot1*lRoots)

! Takes the rlxroot column of the rotation matrix and squares it
call mma_allocate(u_rlx_sq,lroots,label='u_rlx_sq')
u_rlx_sq(:) = si_pdft(:,states(1))*si_pdft(:,states(2))

! Then add the matrix for each state to the ground state
! add put the ground state one in the runfile. Do not
! forget to multiply the (R_IK)^2, where K is "jRoot" below

! DIDA(:,lRoots+1) is currently DI
call Put_DArray('MSPDFTD5',DIDA(:,lRoots+1),nTot1)

call mma_allocate(r_fock_occ,ntot1,label='r_fock_occ')
call dgemv_('n',ntot1,lroots,One,FocMS,ntot1,u_rlx_sq,1,Zero,r_fock_occ,1)
call Put_dArray('FockOcc',r_fock_occ,ntot1)
call mma_deallocate(r_fock_occ)

! DIDA is currently DA over intermediate states
! DIDA for prepp
call mma_allocate(r_active_dm1,ntot1,label='r_active_dm1')
call dgemv_('n',ntot1,lroots,One,DIDA(:,:lroots),ntot1,-u_rlx_sq,1,Zero,r_active_dm1,1)
call Put_DArray('MSPDFTD6',r_active_dm1,nTot1)
call mma_deallocate(r_active_dm1)

! FT99 for bk
call mma_allocate(r_Fxy,ntot4,label='r_fxy')
call dgemv_('n',ntot4,lroots,One,FxyMS,ntot4,u_rlx_sq,1,Zero,r_Fxy,1)
call Put_DArray('FxyMS',r_Fxy,nTot4)
call mma_deallocate(r_Fxy)

! P2MOt for active 2RDM
call mma_allocate(r_casdm2,nacpr2,label='r_casdm2')
call dgemv_('n',nacpr2,lroots,One,P2MOt,nacpr2,u_rlx_sq,1,Zero,r_casdm2,1)
call Put_dArray('P2MOt',r_casdm2,NACPR2)
call mma_deallocate(r_casdm2)

call mma_deallocate(u_rlx_sq)

end subroutine MSPDFTGrad_Misc
