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
  use definitions,only:iwp,wp
  use constants,only:zero,one
  use mspdftgrad,only:F1MS,F2MS,FocMS,FxyMS,P2MOT,D1aoMS,DIDA,D1SAOMS
  use rasscf_global,only:lroots,NACPR2,nTot4
  use general_data,only:ispin,nsym,ntot1,nbas

  implicit none

  real(kind=wp),intent(in) :: si_pdft(lroots,lroots)
  integer(kind=iwp),intent(in) :: states(2)

  integer(kind=iwp) :: ij,iS,jRoot,iBas,jBas
  ! rotated (ie, in final MS-PDFT) eigenbasis quantities
  real(kind=wp) :: r_fock_occ(ntot1),r_active_dm1(ntot1),r_Fxy(ntot4),r_casdm2(nacpr2)
  real(kind=wp) :: u_rlx_sq(lroots)

! Functions added by Paul Calio for MECI Opt *****
! Original calls are in slapaf_util/start_alasaks.f

  Call put_iArray('NACstatesOpt    ',states,2)
! End of stuff added by Paul

  Call Put_DArray('MS_FINAL_ROT    ',si_pdft(:,:),lRoots**2)
  CALL Put_DArray('F1_PDFT         ',F1MS(:,:),nTot1*lRoots)
  CALL Put_DArray('F2_PDFT         ',F2MS(:,:),NACPR2*lRoots)
  CALL Put_DArray('D1AO_MS         ',D1AOMS(:,:),nTot1*lRoots)
  if(ispin /= 1) then
    CALL Put_DArray('D1SAO_MS        ',D1SAOMS(:,:),nTot1*lRoots)
  endif

  ! Takes the rlxroot column of the rotation matrix and squares it
  u_rlx_sq = si_pdft(:,states(1))*si_pdft(:,states(2))
  !u_rlx_sq = si_pdft(:,mcpdft_options%rlxroot)**2

  ! Now storing the density matrix needed for computing 1RDM
  ! First rescale the off-diagonal elements as done in
  ! integral_util/prep.f
  ij = 0
  Do iS = 1,nSym
    do iBas = 1,nBas(iS)
      do jBas = 1,iBas-1
        ij = ij+1
        do jRoot = 1,lRoots+1
          DIDA(ij,jRoot) = 0.5D0*DIDA(ij,jRoot)
        enddo
      enddo
      ij = ij+1
    enddo
  EndDo

  ! Then add the matrix for each state to the ground state
  ! add put the ground state one in the runfile. Do not
  ! forget to multiply the (R_IK)^2, where K is "jRoot" below

  ! DIDA(:,lRoots+1) is currently DI
  CALL Put_DArray('MSPDFTD5        ',DIDA(:,lRoots+1),nTot1)

  call dgemm_('n','n',ntot1,1,lroots,one,FocMS(:,:),ntot1,u_rlx_sq,lroots,zero,r_fock_occ,ntot1)
  Call Put_dArray('FockOcc',r_fock_occ,ntot1)

  ! DIDA is currently DA over intermediate states
  ! DIDA for prepp
  call dgemm_('n','n',ntot1,1,lroots,one,DIDA(:,:lroots),ntot1,-u_rlx_sq,lroots,zero,r_active_dm1,ntot1)
  CALL Put_DArray('MSPDFTD6        ',r_active_dm1,nTot1)

  ! FT99 for bk
  call dgemm_('n','n',ntot4,1,lroots,one,FxyMS(:,:),ntot4,u_rlx_sq,lroots,zero,r_Fxy,ntot4)
  CALL Put_DArray('FxyMS           ',r_Fxy,nTot4)

  ! P2MOt for active 2RDM
  call dgemm_('n','n',nacpr2,1,lroots,one,P2MOt(:,:),nacpr2,u_rlx_sq,lroots,zero,r_casdm2,nacpr2)
  Call Put_dArray('P2MOt',r_casdm2,NACPR2)

EndSubroutine

