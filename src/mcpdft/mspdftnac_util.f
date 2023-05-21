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
      use mspdft, only: iF1MS, iF2MS, iFxyMS, iFocMS, iDIDA, IP2MOt,
     &                  D1AOMS, D1SAOMS, cmsnacstates
#include "WrkSpc.fh"
#include "wadr.fh"
#include "rasdim.fh"
#include "warnings.h"
#include "input_ras_mcpdft.fh"
#include "rasscf.fh"
#include "general.fh"
#include "mspdft.fh"

      real(kind=wp), dimension(lroots**2), intent(in) :: si_pdft
      INTEGER ij,iS,jRoot,iBas,jBas
      Real*8 RJKRIK

******* Functions added by Paul Calio for MECI Opt *****
******* Original calls are in slapaf_util/start_alasaks.f
      Logical :: CalcNAC_Opt = .False.
      Logical :: MECI_via_SLAPAF = .False.
      INTEGER NACstatesOpt(2)

      NACstatesOpt(1)=cmsNACstates(1)
      NACstatesOpt(2)=cmsNACstates(2)

      call put_lscalar('MECI_via_SLAPAF ', MECI_via_SLAPAF)
      Call put_iArray('NACstatesOpt    ', NACstatesOpt,2)
      Call Put_lscalar('CalcNAC_Opt     ', CalcNAC_Opt)
****** End of stuff added by Paul


      Call Put_DArray('MS_FINAL_ROT    ',si_pdft,     lRoots**2)
      CALL Put_DArray('F1MS            ',Work(iF1MS),  nTot1*nRoots)
      CALL Put_DArray('F2MS            ',Work(iF2MS), NACPR2*nRoots)
      CALL Put_DArray('D1AO_MS         ',Work(D1AOMS), nTot1*nRoots)
      if (ispin.ne.1)
     &CALL Put_DArray('D1SAO_MS        ',Work(D1SAOMS),nTot1*nRoots)

**********Fock_Occ Part
      CALL FZero(Work(ipFocc),ntot1)
      DO JRoot=1,lRoots

        RJKRIK = si_pdft((cmsnacstates(2)-1)*lroots+jroot) *
     &     si_pdft((cmsnacstates(1)-1)*lroots+jroot)
!       RJKRIK = Work(LHRot+(cmsNACstates(2)-1)*lRoots+jRoot-1)*
!    &     Work(LHRot+(cmsNACstates(1)-1)*lRoots+jRoot-1)

        CALL daXpY_(ntot1,RJKRIK,
     &           Work(iFocMS+(JRoot-1)*nTot1),1,Work(ipFocc),1)
      END DO
      Call Put_dArray('FockOcc',Work(ipFocc),ntot1)
**********Now storing the density matrix needed for computing 1RDM
**********First rescale the off-diagonal elements as done in
**********integral_util/prep.f
      ij=-1
      Do iS=1,nSym
       do iBas=1,nBas(iS)
       do jBas=1,iBas-1
        ij=ij+1
        do jRoot=1,nRoots+1
         Work(iDIDA+ij+(jRoot-1)*nTot1)=0.5D0*
     &   Work(iDIDA+ij+(jRoot-1)*nTot1)
        end do
       end do
       ij=ij+1
       end do
      End Do
**********Then add the matrix for each state to the ground state
**********add put the ground state one in the runfile. Do not
**********forget to multiply the (R_IK)^2, where K is "jRoot" below
      RJKRIK=si_pdft(1+(cmsnacstates(2)-1)*lRoots)*
     &    si_pdft(1+(cmsnacstates(1)-1)*lroots)
!     RJKRIK=Work(LHRot+(cmsNACstates(2)-1)*lRoots)*
!    &    Work(LHRot+(cmsNACstates(1)-1)*lRoots)
      CALL DScal_(nTot1,-RJKRIK,Work(iDIDA),1)
      CALL DScal_(nTot4,RJKRIK,Work(iFxyMS),1)
      CALL DScal_(NACPR2,RJKRIK,Work(iP2MOt),1)

      ij=0
      jRoot=1
      Do jRoot=2,lRoots
        ij=0
        RJKRIK = si_pdft((cmsnacstates(2)-1)*lroots+jroot) *
     &     si_pdft((cmsnacstates(1)-1)*lroots+jroot)
!       RJKRIK = Work(LHRot+(cmsNACstates(2)-1)*lRoots+jRoot-1)*
!    &     Work(LHRot+(cmsNACstates(1)-1)*lRoots+jRoot-1)
*******DIDA for prepp
       CALL DaXpY_(nTot1,-RJKRIK,
     &            Work(iDIDA+(jRoot-1)*nTot1),1,Work(iDIDA),1)
*******FT99 for bk
       CALL DaXpY_(nTot4,RJKRIK,
     &            Work(iFxyMS+(jRoot-1)*nTot4),1,Work(iFxyMS),1)
*******P2MOt for active 2RDM
       CALL DaXpY_(NACPR2,RJKRIK,
     &            Work(IP2MOt+(jRoot-1)*NACPR2),1,Work(iP2MOt),1)
      End Do

      CALL Put_DArray('MSPDFTD6        ',Work(iDIDA),nTot1)
********Work(iDIDA+lRoots*nTot1) is currently DI
      CALL Put_DArray('FxyMS           ',Work(iFxyMS), nTot4)
      Call Put_dArray('P2MOt',Work(iP2MOt),NACPR2)

      ! Some other things that were initially in mcpdft.f
      Call Put_cArray('Relax Method','MSPDFT  ',8)
      Call Put_cArray('MCLR Root','****************',16)
      Call Put_iScalar('Relax CASSCF root',irlxroot)
      RETURN
      End Subroutine

