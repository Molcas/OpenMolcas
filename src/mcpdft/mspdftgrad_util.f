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
* Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
!                                                                      *
! 2023, Matthew R. Hennefarth - modified lhrot -> si_pdft              *
!***********************************************************************

        Subroutine MSPDFTGrad_Misc(si_pdft)
********This subroutine does miscellaneous things needed
********in MS-PDFT gradient calculation.
      use definitions, only: wp
      use mspdft, only: F1MS, F2MS, FxyMS, FocMS, DIDA, P2MOt,
     &                  D1AOMS, D1SAOMS
      use wadr, only: FockOcc
#include "rasdim.fh"
#include "warnings.h"
#include "rasscf.fh"
#include "general.fh"

      INTEGER ij,iS,jRoot,iBas,jBas
      Real*8 RIK2
      real(kind=wp), dimension(lroots,lroots), intent(in) :: si_pdft

******* Functions added by Paul Calio for MECI Opt *****
******* Original calls are in slapaf_util/start_alasaks.f
      Logical :: CalcNAC_Opt = .False.
      Logical :: MECI_via_SLAPAF = .False.
      INTEGER NACstatesOpt(2)

      NACstatesOpt(1)=irlxroot
      NACstatesOpt(2)=0

      Call put_iArray('NACstatesOpt    ', NACstatesOpt,2)
      Call Put_lscalar('CalcNAC_Opt     ', CalcNAC_Opt)
      call put_lscalar('MECI_via_SLAPAF ', MECI_via_SLAPAF)
****** End of stuff added by Paul


      Call Put_DArray('MS_FINAL_ROT    ',si_pdft(:,:),   lRoots**2)
      CALL Put_DArray('F1MS            ',F1MS(:,:),  nTot1*nRoots)
      CALL Put_DArray('F2MS            ',F2MS(:,:), NACPR2*nRoots)
      CALL Put_DArray('D1AO_MS         ',D1AOMS(:,:), nTot1*nRoots)
      if (ispin.ne.1)
     &CALL Put_DArray('D1SAO_MS        ',D1SAOMS(:,:),nTot1*nRoots)


**********Fock_Occ Part
      FockOcc(:)=0.0D0
      DO JRoot=1,lRoots
      call daXpY_(ntot1,si_pdft(jroot,irlxroot)**2,
     &           FocMS(:,JRoot),1,FockOcc,1)
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
      RIK2=si_pdft(1,irlxroot)**2
      CALL DScal_(nTot1,-RIK2,DIDA(:,1),1)
      CALL DScal_(nTot4,RIK2,FxyMS(:,1),1)
      CALL DScal_(NACPR2,RIK2,P2MOt(:,1),1)
      ij=0
      jRoot=1
      ! for the comment below, lhrot -> si_pdft
********LHRot(1,iRlxRoot) should give me R_I1
      Do jRoot=2,lRoots
      ij=0
       RIK2=si_pdft(jroot,irlxroot)**2
*******DIDA for prepp
       CALL DaXpY_(nTot1,-RIK2,DIDA(:,jRoot),1,DIDA(:,1),1)
*******FT99 for bk
       CALL DaXpY_(nTot4,RIK2,FxyMS(:,jRoot),1,FxyMS(:,1),1)
*******P2MOt for active 2RDM
       CALL DaXpY_(NACPR2,RIK2,P2MOt(:,jRoot),1,P2MOt(:,1),1)
      End Do
********DIDA is currently DA over intermediate states
      CALL Put_DArray('MSPDFTD6        ',DIDA(:,1),nTot1)
********DIDA(:,lRoots+1) is currently DI
      CALL Put_DArray('MSPDFTD5        ',DIDA(:,lRoots+1),nTot1)
      CALL Put_DArray('FxyMS           ',FxyMS(:,1), nTot4)
      Call Put_dArray('P2MOt',P2MOt(:,1),NACPR2)

      ! Some other things that were initially in mcpdft.f
      Call Put_cArray('Relax Method','MSPDFT  ',8)
      Call Put_cArray('MCLR Root','****************',16)
      Call Put_iScalar('Relax CASSCF root',irlxroot)
      RETURN
      End Subroutine

