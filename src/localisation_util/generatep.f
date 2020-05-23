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
* Copyright (C) Yannick Carissan                                       *
*               2005, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine GenerateP(Ovlp,cMO,Name,nBasis,nOrb2Loc,nAtoms,
     &                     nBas_per_Atom,nBas_Start,PA,Debug)
c
c     Author: Yannick Carissan.
c
c     Modifications:
c        - October 6, 2005 (Thomas Bondo Pedersen):
c          Reduce operation count and use BLAS.
c
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
#include "Molcas.fh"
      Real*8, Allocatable:: SBar(:,:)
      Integer nBas_per_Atom(*),nBas_Start(*)
      Real*8 cMO(nBasis,*),Ovlp(nBasis,nBasis)
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
      Logical Debug
      Character*(LENIN8) Name(*)

      Call mma_Allocate(SBar,nBasis,nOrb2Loc,Label='SBar')
      Call GenerateP_1(Ovlp,cMO,Sbar,Name,nBasis,nOrb2Loc,
     &                 nAtoms,nBas_per_Atom,nBas_Start,PA,Debug)
      Call mma_deallocate(SBar)

      End
      Subroutine GenerateP_1(Ovlp,cMO,Sbar,Name,nBasis,nOrb2Loc,nAtoms,
     &                       nBas_per_Atom,nBas_Start,PA,Debug)
c
c     Author: Yannick Carissan.
c
c     Modifications:
c        - October 6, 2005 (Thomas Bondo Pedersen):
c          Reduce operation count and use BLAS.
c
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
#include "Molcas.fh"
      Integer nBas_per_Atom(*),nBas_Start(*)
      Real*8 cMO(nBasis,*),Ovlp(nBasis,nBasis)
      Real*8 Sbar(nBasis,nOrb2Loc)
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
      Logical Debug
      Character*(LENIN8) Name(*),PALbl
c
c---- Compute Sbar(mu,s) = sum_{nu} Ovlp(mu,nu) * cMO(nu,s)
c
      Call DGEMM_('N','N',nBasis,nOrb2Loc,nBasis,
     &           One,Ovlp,nBasis,cMO,nBasis,
     &           Zero,Sbar,nBasis)
c
      Do iAt=1,nAtoms
c
c------ Compute MA(s,t) = sum_{mu_in_A} cMO(mu,s) * Sbar(mu,t)
c
        Call DGEMM_('T','N',
     &             nOrb2Loc,nOrb2Loc,nBas_per_Atom(iAt),
     &             One,cMO(nBas_Start(iAt),1),nBasis,
     &                 Sbar(nBas_Start(iAt),1),nBasis,
     &             Zero,PA(1,1,iAt),nOrb2Loc)
c
c------ Compute <s|PA|t> by symmetrization of MA.
c
        Do iMO_s=1,nOrb2Loc
          Do iMO_t=iMO_s+1,nOrb2Loc
            PAst = PA(iMO_s,iMO_t,iAt)
            PAts = PA(iMO_t,iMO_s,iAt)
            PA(iMO_s,iMO_t,iAt) = Half*(PAst+PAts)
            PA(iMO_t,iMO_s,iAt) = PA(iMO_s,iMO_t,iAt)
          End Do !iMO_t
        End Do !iMO_s
c
      End Do !iAt
c
      If (Debug) Then
        Write(6,*) 'In GenerateP'
        Write(6,*) '------------'
        Do iAt=1,nAtoms
          PALbl='PA__'//Name(nBas_Start(iAt))(1:LENIN)
          Call RecPrt(PALbl,' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
        End Do
      End If
c
      Return
      End
