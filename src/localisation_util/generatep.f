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
#define _Test4_
#ifdef _Test4_
      Subroutine GenerateP(Ovlp,cMO,Name,nBasis,nOrb2Loc,nAtoms,
     &                     iTab_ptr,nBas_per_Atom,nBas_Start,PA,Debug)
#else
      Subroutine GenerateP(Ovlp,cMO,Name,nBasis,nOrb2Loc,nAtoms,
     &                     iTab_ptr,nBas_per_Atom,nBas_Start,Debug)
#endif
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
      Integer nBas_per_Atom(*),nBas_Start(*),iTab_Ptr(*)
      Real*8 cMO(nBasis,*),Ovlp(nBasis,nBasis)
#ifdef _Test4_
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
#endif
      Logical Debug
      Character*(LENIN8) Name(*)

      Call mma_Allocate(SBar,nBasis,nOrb2Loc,Label='SBar')
#ifdef _Test4_
      Call GenerateP_1(Ovlp,cMO,Sbar,Name,nBasis,nOrb2Loc,
     &                 nAtoms,iTab_ptr,nBas_per_Atom,nBas_Start,PA,
     &                 Debug)
#else
      Call GenerateP_1(Ovlp,cMO,Sbar,Name,nBasis,nOrb2Loc,
     &                 nAtoms,iTab_ptr,nBas_per_Atom,nBas_Start,
     &                 Debug)
#endif
      Call mma_deallocate(SBar)

      End
#ifdef _Test4_
      Subroutine GenerateP_1(Ovlp,cMO,Sbar,Name,nBasis,nOrb2Loc,nAtoms,
     &                       iTab_ptr,nBas_per_Atom,nBas_Start,PA,Debug)
#else
      Subroutine GenerateP_1(Ovlp,cMO,Sbar,Name,nBasis,nOrb2Loc,nAtoms,
     &                       iTab_ptr,nBas_per_Atom,nBas_Start,Debug)
#endif
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
      Integer nBas_per_Atom(*),nBas_Start(*),iTab_Ptr(*)
      Real*8 cMO(nBasis,*),Ovlp(nBasis,nBasis)
      Real*8 Sbar(nBasis,nOrb2Loc)
#ifdef _Test4_
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
#endif
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
c------ The array iTab_ptr contains the value of the pointer to the
c       PA array for atom iAt
c
        ip  = iTab_ptr(iAt)
        ip0 = ip - 1
c
c------ Compute MA(s,t) = sum_{mu_in_A} cMO(mu,s) * Sbar(mu,t)
c
        Call DGEMM_('T','N',
     &             nOrb2Loc,nOrb2Loc,nBas_per_Atom(iAt),
     &             One,cMO(nBas_Start(iAt),1),nBasis,
     &                 Sbar(nBas_Start(iAt),1),nBasis,
     &             Zero,Work(ip),nOrb2Loc)
      Call RecPrt('Work(ip)',' ',Work(ip),nOrb2Loc,nOrb2Loc)
*
#ifdef _Test4_
        Call DGEMM_('T','N',
     &             nOrb2Loc,nOrb2Loc,nBas_per_Atom(iAt),
     &             One,cMO(nBas_Start(iAt),1),nBasis,
     &                 Sbar(nBas_Start(iAt),1),nBasis,
     &             Zero,PA(1,1,iAt),nOrb2Loc)
      Call RecPrt('PA(1,1,iAt)',' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
#endif
c
c------ Compute <s|PA|t> by symmetrization of MA.
c
        Do iMO_s=1,nOrb2Loc
          Do iMO_t=iMO_s+1,nOrb2Loc
            mAd_st = nOrb2Loc*(iMO_t-1)+iMO_s
            mAd_ts = nOrb2Loc*(iMO_s-1)+iMO_t
            PAst = Work(ip0+mAd_st)
            PAts = Work(ip0+mAd_ts)
            Work(ip0+mAd_st) = Half*(PAst+PAts)
            Work(ip0+mAd_ts) = Work(ip0+mAd_st)
      Call RecPrt('Work(iq)',' ',Work(ip),nOrb2Loc,nOrb2Loc)
#ifdef _Test4_
            PAst = PA(iMO_s,iMO_t,iAt)
            PAts = PA(iMO_t,iMO_s,iAt)
            PA(iMO_s,iMO_t,iAt) = Half*(PAst+PAts)
            PA(iMO_t,iMO_s,iAt) = PA(iMO_s,iMO_t,iAt)
      Call RecPrt('QA(1,1,iAt)',' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
#endif
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
          ip=iTab_ptr(iAt)
          Call RecPrt(PALbl,' ',Work(ip),nOrb2Loc,nOrb2Loc)
#ifdef _Test4_
          Call RecPrt(PALbl,' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
#endif
        End Do
      End If
c
      Return
      End
