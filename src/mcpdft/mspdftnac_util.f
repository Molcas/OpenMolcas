************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2023, Paul B Calio                                     *
* Based on MSPDFTGrad_Misc from Jie J. Bao                             *
************************************************************************
        Subroutine MSPDFTNAC_Misc(si_pdft)
********This subroutine does miscellaneous things needed
********in MS-PDFT NAC calculation.
      use definitions, only: wp
      use mspdft, only: F1MS, F2MS, FxyMS, FocMS, DIDA, P2MOt,
     &                  D1AOMS, D1SAOMS
      use wadr, only: FockOcc
      use mcpdft_input, only: mcpdft_options
#include "rasdim.fh"
#include "warnings.h"
#include "rasscf.fh"
#include "general.fh"

      real(kind=wp), dimension(lroots**2), intent(in) :: si_pdft
      INTEGER ij,iS,jRoot,iBas,jBas, bra_state, ket_state
      Real*8 RJKRIK

******* Functions added by Paul Calio for MECI Opt *****
******* Original calls are in slapaf_util/start_alasaks.f
      Logical :: CalcNAC_Opt = .False.
      Logical :: MECI_via_SLAPAF = .False.

      bra_state = mcpdft_options%nac_states(1)
      ket_state = mcpdft_options%nac_states(2)


      call put_lscalar('MECI_via_SLAPAF ', MECI_via_SLAPAF)
      Call put_iArray('NACstatesOpt    ', mcpdft_options%nac_states,2)
      Call Put_lscalar('CalcNAC_Opt     ', CalcNAC_Opt)
****** End of stuff added by Paul


      Call Put_DArray('MS_FINAL_ROT    ',si_pdft(:),lRoots**2)
      CALL Put_DArray('F1MS            ',F1MS(:,:),  nTot1*nRoots)
      CALL Put_DArray('F2MS            ',F2MS(:,:), NACPR2*nRoots)
      CALL Put_DArray('D1AO_MS         ',D1AOMS(:,:), nTot1*nRoots)
      if (ispin.ne.1)
     &CALL Put_DArray('D1SAO_MS        ',D1SAOMS(:,:),nTot1*nRoots)

**********Fock_Occ Part
      FockOcc(:)=0.0D0
      DO JRoot=1,lRoots

        RJKRIK = si_pdft((ket_state-1)*lroots+jroot) *
     &     si_pdft((bra_state-1)*lroots+jroot)

        CALL daXpY_(ntot1,RJKRIK,FocMS(:,JRoot),1,FockOcc,1)
      END DO
      Call Put_dArray('FockOcc',FockOcc,ntot1)
**********Now storing the density matrix needed for computing 1RDM
**********First rescale the off-diagonal elements as done in
**********integral_util/prep.f
      ij=0
      Do iS=1,nSym
       do iBas=1,nBas(iS)
       do jBas=1,iBas-1
        ij=ij+1
        do jRoot=1,nRoots+1
         DIDA(ij,jRoot)=0.5D0*DIDA(ij,jRoot)
        end do
       end do
       ij=ij+1
       end do
      End Do
**********Then add the matrix for each state to the ground state
**********add put the ground state one in the runfile. Do not
**********forget to multiply the (R_IK)^2, where K is "jRoot" below
      RJKRIK=si_pdft(1+(ket_state-1)*lRoots)*
     &    si_pdft(1+(bra_state-1)*lroots)
      CALL DScal_(nTot1,-RJKRIK,DIDA(:,1),1)
      CALL DScal_(nTot4,RJKRIK,FxyMS(:,1),1)
      CALL DScal_(NACPR2,RJKRIK,P2MOt(:,1),1)

      ij=0
      jRoot=1
      Do jRoot=2,lRoots
        ij=0
        RJKRIK = si_pdft((ket_state-1)*lroots+jroot) *
     &     si_pdft((bra_state-1)*lroots+jroot)
*******DIDA for prepp
       CALL DaXpY_(nTot1,-RJKRIK,DIDA(:,jRoot),1,DIDA(:,1),1)
*******FT99 for bk
       CALL DaXpY_(nTot4,RJKRIK,FxyMS(:,jRoot),1,FxyMS(:,1),1)
*******P2MOt for active 2RDM
       CALL DaXpY_(NACPR2,RJKRIK,P2MOt(:,jRoot),1,P2MOt(:,1),1)
      End Do

      CALL Put_DArray('MSPDFTD6        ',DIDA(:,1),nTot1)
********DIDA(:,lRoots+1) is currently DI
      CALL Put_DArray('FxyMS           ',FxyMS(:,1), nTot4)
      Call Put_dArray('P2MOt',P2MOt(:,1),NACPR2)

      ! Some other things that were initially in mcpdft.f
      Call Put_cArray('Relax Method','MSPDFT  ',8)
      Call Put_cArray('MCLR Root','****************',16)
      Call Put_iScalar('Relax CASSCF root',irlxroot)

      End Subroutine

