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
************************************************************************
        Subroutine MSPDFTGrad_Misc(LHRot)
********This subroutine does miscellaneous things needed
********in MS-PDFT gradient calculation.
#include "WrkSpc.fh"
#include "wadr.fh"
#include "rasdim.fh"
#include "warnings.h"
#include "input_ras_mcpdft.fh"
#include "rasscf.fh"
#include "general.fh"
#include "mspdft.fh"

      INTEGER LHRot,ij,iS,jRoot,iBas,jBas
      Real*8 RIK2

      Call Put_DArray('MS_FINAL_ROT    ',Work(LHRot),     lRoots**2)
      CALL Put_DArray('F1MS            ',Work(iF1MS),  nTot1*nRoots)
      CALL Put_DArray('F2MS            ',Work(iF2MS), NACPR2*nRoots)
      CALL Put_DArray('D1AO_MS         ',Work(D1AOMS), nTot1*nRoots)
      if (ispin.ne.1)
     &CALL Put_DArray('D1SAO_MS        ',Work(D1SAOMS),nTot1*nRoots)


**********Fock_Occ Part
      CALL FZero(Work(ipFocc),ntot1)
      DO JRoot=1,lRoots
      CALL daXpY_(ntot1,Work(LHRot+(irlxroot-1)*lRoots+jRoot-1)**2,
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
      RIK2=Work(LHRot+(irlxroot-1)*lRoots)**2
      CALL DScal_(nTot1,-RIK2,Work(iDIDA),1)
      CALL DScal_(nTot4,RIK2,Work(iFxyMS),1)
      CALL DScal_(NACPR2,RIK2,Work(iP2MOt),1)
      ij=0
      jRoot=1
********Work(LHRot+(iRlxRoot-1)*lRoots) should give me R_I1
      Do jRoot=2,lRoots
      ij=0
       RIK2=Work(LHRot+(irlxroot-1)*lRoots+jRoot-1)**2
*******DIDA for prepp
       CALL DaXpY_(nTot1,-RIK2,
     &            Work(iDIDA+(jRoot-1)*nTot1),1,Work(iDIDA),1)
*******FT99 for bk
       CALL DaXpY_(nTot4,RIK2,
     &            Work(iFxyMS+(jRoot-1)*nTot4),1,Work(iFxyMS),1)
*******P2MOt for active 2RDM
       CALL DaXpY_(NACPR2,RIK2,
     &            Work(IP2MOt+(jRoot-1)*NACPR2),1,Work(iP2MOt),1)
      End Do
********Work(iDIDA) is currently DA over intermediate states
      CALL Put_DArray('MSPDFTD6        ',Work(iDIDA),nTot1)
********Work(iDIDA+lRoots*nTot1) is currently DI
      CALL Put_DArray('MSPDFTD5        ',
     &                 Work(iDIDA+lRoots*nTot1),nTot1)
      CALL Put_DArray('FxyMS           ',Work(iFxyMS), nTot4)
      Call Put_dArray('P2MOt',Work(iP2MOt),NACPR2)
      RETURN
      End Subroutine

