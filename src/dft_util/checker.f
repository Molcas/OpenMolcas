************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Checker(mGrid,Rho,nRho,P2_ontop,
     &                   nP2_ontop,iSpin,F_xc,
     &                   dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_Z)
      Implicit Real*8 (A-H,O-Z)
      External LSDA, BLYP, BPBE, B3LYP, B2PLYP, HFS, HFB, HFO,
     &  XAlpha, LSDA5, B3LYP5,TLYP,NLYP, OLYP, O3LYP, OPBE,
     &  SSBSW, SSBD, PBE, PBESOL, PBE0, M06L, M06, M062X, M06HF, O2PLYP,
     & HFG, GLYP, GPBE, HFB86, B86LYP, B86PBE, BWIG, KT3,
     & KT2, RGE2, PTCA
      Integer Functional_type
#include "functional_types.fh"
#include "WrkSpc.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('F_xc1','Allo','Real',ip_F_xc1,mGrid)
      Call GetMem('F_xc2','Allo','Real',ip_F_xc2,mGrid)
      Call GetMem('F_xc3','Allo','Real',ip_F_xc3,mGrid)
      Call GetMem('F_xc4','Allo','Real',ip_F_xc4,mGrid)
      Call GetMem('dF_temp','Allo','Real',ip_dF_temp,ndF_dRho*mGrid)
      Call GetMem('Rho_temp','Allo','Real',ip_Rho_temp,nRho*mGrid)
      T_x=1.0D-8
*                                                                      *
************************************************************************
*                                                                      *
      Functional_type=LDA_type
*---- LSDA
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              LSDA,Functional_type,'LSDA',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- LSDA5
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              LSDA5,Functional_type,'LSDA5',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- HFS
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              HFS,Functional_type,'HFS',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- xalpha
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              XAlpha,Functional_type,'xAlpha',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*
      Functional_type=GGA_type
*---- HFB
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              HFB,Functional_type,'HFB',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- HFO
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              HFO,Functional_type,'HFO',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- HFB86
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              HFB86,Functional_type,'HFB86',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- HFG
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              HFG,Functional_type,'HFG',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- BWIG
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              BWIG,Functional_type,'BWIG',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- BLYP
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              BLYP,Functional_type,'BLYP',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- OLYP
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              OLYP,Functional_type,'OLYP',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- KT3
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              KT3,Functional_type,'KT3',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- KT2
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              KT2,Functional_type,'KT2',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- GLYP
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              GLYP,Functional_type,'GLYP',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- B86LYP
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              B86LYP,Functional_type,'B86LYP',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- BPBE
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              BPBE,Functional_type,'BPBE',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- OPBE
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              OPBE,Functional_type,'OPBE',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- GPBE
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              GPBE,Functional_type,'GPBE',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- B86PBE
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              B86PBE,Functional_type,'B86PBE',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- TLYP
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              TLYP,Functional_type,'TLYP',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
#ifdef _SKIP_
*---- NLYP
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              NLYP,Functional_type,'NLYP',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
#endif
*---- B3LYP
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              B3LYP,Functional_type,'B3LYP',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- O3LYP
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              O3LYP,Functional_type,'O3LYP',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- B2PLYP
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              B2PLYP,Functional_type,'B2PLYP',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- O2PLYP
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              O2PLYP,Functional_type,'O2PLYP',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- B3LYP5
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              B3LYP5,Functional_type,'B3LYP5',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- PBE0
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              PBE0,Functional_type,'PBE0',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- PBE
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              PBE,Functional_type,'PBE',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- SSBSW
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              SSBSW,Functional_type,'SSBSW',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- SSBD
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              SSBD,Functional_type,'SSBD',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- PBEsol
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              PBEsol,Functional_type,'PBESOL',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- RGE2
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              RGE2,Functional_type,'RGE2',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- PTCA
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              PTCA,Functional_type,'PTCA',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*
*
      Functional_type=meta_GGA_type1
*---- M06-L
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              M06L,Functional_type,'M06-L',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- M06
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              M06,Functional_type,'M06',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- M06-2X
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              M062X,Functional_type,'M06-2X',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*---- M06-HF
      Call Checker_(mGrid,Rho,nRho,P2_ontop,
     &              nP2_ontop,iSpin,F_xc,
     &              dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &              M06HF,Functional_type,'M06-HF',
     &              Work(ip_F_xc1),Work(ip_F_xc2),
     &              Work(ip_F_xc3),Work(ip_F_xc4),
     &              Work(ip_dF_temp),Work(ip_Rho_temp))
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('Rho_temp','Free','Real',ip_Rho_temp,nRho*mGrid)
      Call GetMem('dF_temp','Free','Real',ip_dF_temp,ndF_dRho*mGrid)
      Call GetMem('F_xc4','Free','Real',ip_F_xc4,mGrid)
      Call GetMem('F_xc3','Free','Real',ip_F_xc3,mGrid)
      Call GetMem('F_xc2','Free','Real',ip_F_xc2,mGrid)
      Call GetMem('F_xc1','Free','Real',ip_F_xc1,mGrid)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(T_Z)
      End
      Subroutine Checker_(mGrid,Rho,nRho,P2_ontop,
     &                   nP2_ontop,iSpin,F_xc,
     &                   dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X,
     &                   Kernel,Ftype,FName,F_xc1,F_xc2,F_xc3,F_xc4,
     &                   dF_temp,Rho_temp)
      Implicit Real*8 (A-H,O-Z)
      External Kernel
#include "functional_types.fh"
#include "nq_index.fh"
#include "real.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid), F_xc1(mGrid),F_xc2(mGrid),
     &       F_xc3(mGrid), F_xc4(mGrid), dF_temp(ndF_dRho,mGrid),
     &       Rho_temp(nRho,mGrid)
      Integer FType
      Character FName*(*)
*                                                                      *
************************************************************************
*                                                                      *
C     Write (*,*) 'Functional: ',FName
      Call FZero(F_xc,mGrid)
      Call FZero(dF_dRho,ndF_dRho*mGrid)
      Call Kernel(mGrid,Rho,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc,
     &            dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*                                                                      *
************************************************************************
*                                                                      *
      Delta=1.0D-3
      Error=1.0D-4
      Thr_V_a=1.0D-8
      Tf=(One-Two*Delta)
*                                                                      *
************************************************************************
*                                                                      *
*     Check for rho (alpha)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipR,iGrid) = (One+Delta)*Rho(ipR,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipRa,iGrid) = (One+Delta)*Rho(ipRa,iGrid)
         End Do
      End If
      Call FZero(F_xc1,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc1,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipR,iGrid) = (One-Delta)*Rho(ipR,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipRa,iGrid) = (One-Delta)*Rho(ipRa,iGrid)
         End Do
      End If
      Call FZero(F_xc2,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc2,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipR,iGrid) = (One+Two*Delta)*Rho(ipR,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipRa,iGrid) = (One+Two*Delta)*Rho(ipRa,iGrid)
         End Do
      End If
      Call FZero(F_xc3,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc3,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipR,iGrid) = (One-Two*Delta)*Rho(ipR,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipRa,iGrid) = (One-Two*Delta)*Rho(ipRa,iGrid)
         End Do
      End If
      Call FZero(F_xc4,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc4,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      Do iGrid = 1, mGrid
         If (iSpin.eq.1) Then
            If (Tf*Rho(ipR,iGrid).lt.T_x) Go To 99
            V_a=dF_dRho(ipR,iGrid)
            Delta0=Delta*Rho(ipR,iGrid)*Two
         Else
            If (Tf*Rho(ipRa,iGrid).lt.T_x) Go To 99
            V_a=dF_dRho(ipRa,iGrid)
            Delta0=Delta*Rho(ipRa,iGrid)
         End If
         Diff1=F_xc1(iGrid)-F_xc2(iGrid)
         Diff2=F_xc3(iGrid)-F_xc4(iGrid)
         V_n=(Two*Diff1-Diff2/Four)/(Three*Delta0)
         Rel_Diff=Abs((V_a-V_n)/V_n)
         If (V_a-V_n.eq.Zero) Rel_Diff=Zero
C        If (V_a.eq.Zero) Rel_Diff=Abs(V_n)
         If (Abs(V_a).lt.Thr_V_a) Rel_Diff=Abs(V_n)
         If (Rel_Diff.gt.Error) Then
            Call WarningMessage(2,'Rel_Diff.gt.Error')
            Write (6,'(A,A)') 'Functional:',FName
            Write (6,*) 'Error in dF_dRho(alpha)!'
            Write (6,*) 'V_a,V_n=',V_a,V_n
            Write (6,*) 'iGrid,mGrid=',iGrid,mGrid
            If (iSpin.eq.1) Then
               Write (6,*) 'Rho=',Rho(ipR,iGrid)
               Write (6,*) 'GradRho=',Rho(ipdRx,iGrid),
     &                                Rho(ipdRy,iGrid),
     &                                Rho(ipdRz,iGrid)
               Write (6,*) 'Tau=',    Rho(ipTau,iGrid)
            Else
               Write (6,*) 'Rho=',Rho(ipRa,iGrid),
     &                            Rho(ipRb,iGrid)
               Write (6,*) 'GradRho=',Rho(ipdRxa,iGrid),
     &                                Rho(ipdRya,iGrid),
     &                                Rho(ipdRza,iGrid)
               Write (6,*) '        ',Rho(ipdRxb,iGrid),
     &                                Rho(ipdRyb,iGrid),
     &                                Rho(ipdRzb,iGrid)
               Write (6,*) 'Tau=',    Rho(ipTaua,iGrid),
     &                                Rho(ipTaub,iGrid)
            End If
            Write (6,*) 'F_xc(iGrid),F_xc1(iGrid),F_xc2(iGrid)=',
     &                   F_xc(iGrid),F_xc1(iGrid),F_xc2(iGrid)
            Call Abend()
         End If
 99      Continue
      End Do
*
      If (iSpin.eq.1) Go To 100
*                                                                      *
************************************************************************
*                                                                      *
*     Check for rho (beta)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipRb,iGrid) = (One+Delta)*Rho(ipRb,iGrid)
      End Do
      Call FZero(F_xc1,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc1,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipRb,iGrid) = (One-Delta)*Rho(ipRb,iGrid)
      End Do
      Call FZero(F_xc2,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc2,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipRb,iGrid) = (One+Two*Delta)*Rho(ipRb,iGrid)
      End Do
      Call FZero(F_xc3,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc3,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipRb,iGrid) = (One-Two*Delta)*Rho(ipRb,iGrid)
      End Do
      Call FZero(F_xc4,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc4,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      Do iGrid = 1, mGrid
         If (Tf*Rho(ipRb,iGrid).lt.T_x) Go To 101
         V_a=dF_dRho(ipRb,iGrid)
         Delta0=Delta*Rho(ipRb,iGrid)
         Diff1=F_xc1(iGrid)-F_xc2(iGrid)
         Diff2=F_xc3(iGrid)-F_xc4(iGrid)
         V_n=(Two*Diff1-Diff2/Four)/(Three*Delta0)
         Rel_Diff=Abs((V_a-V_n)/V_n)
         If (V_a-V_n.eq.Zero) Rel_Diff=Zero
C        If (V_a.eq.Zero) Rel_Diff=Abs(V_n)
         If (Abs(V_a).lt.Thr_V_a) Rel_Diff=Abs(V_n)
         If (Rel_Diff.gt.Error) Then
            Call WarningMessage(2,'Rel_Diff.gt.Error')
            Write (6,'(A,A)') 'Functional:',FName
            Write (6,*) 'Error in dF_dRho(beta)!'
            Write (6,*) 'V_a,V_n=',V_a,V_n
            Write (6,*) 'iGrid,mGrid=',iGrid,mGrid
            Write (6,*) 'Rho=',Rho(ipRa,iGrid),
     &                         Rho(ipRb,iGrid)
            Write (6,*) 'GradRho=',Rho(ipdRxa,iGrid),
     &                             Rho(ipdRya,iGrid),
     &                             Rho(ipdRza,iGrid)
            Write (6,*) '        ',Rho(ipdRxb,iGrid),
     &                             Rho(ipdRyb,iGrid),
     &                             Rho(ipdRzb,iGrid)
            Write (6,*) 'Tau=',    Rho(ipTaua,iGrid),
     &                             Rho(ipTaub,iGrid)
            Write (6,*) 'F_xc(iGrid),F_xc1(iGrid),F_xc2(iGrid)=',
     &                   F_xc(iGrid),F_xc1(iGrid),F_xc2(iGrid)
            Call Abend()
         End If
 101     Continue
      End Do
*
 100  Continue
      If (FType.eq.LDA_type) Return
*                                                                      *
************************************************************************
*                                                                      *
*    Check for grad rho (alpha)
*    What we will do is actually check the derivative wrt the norm of
*    grad rho alpha.
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipdRx,iGrid) = (One+Delta)*Rho(ipdRx,iGrid)
            Rho_Temp(ipdRy,iGrid) = (One+Delta)*Rho(ipdRy,iGrid)
            Rho_Temp(ipdRz,iGrid) = (One+Delta)*Rho(ipdRz,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipdRxa,iGrid) = (One+Delta)*Rho(ipdRxa,iGrid)
            Rho_Temp(ipdRya,iGrid) = (One+Delta)*Rho(ipdRya,iGrid)
            Rho_Temp(ipdRza,iGrid) = (One+Delta)*Rho(ipdRza,iGrid)
         End Do
      End If
      Call FZero(F_xc1,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc1,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipdRx,iGrid) = (One-Delta)*Rho(ipdRx,iGrid)
            Rho_Temp(ipdRy,iGrid) = (One-Delta)*Rho(ipdRy,iGrid)
            Rho_Temp(ipdRz,iGrid) = (One-Delta)*Rho(ipdRz,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipdRxa,iGrid) = (One-Delta)*Rho(ipdRxa,iGrid)
            Rho_Temp(ipdRya,iGrid) = (One-Delta)*Rho(ipdRya,iGrid)
            Rho_Temp(ipdRza,iGrid) = (One-Delta)*Rho(ipdRza,iGrid)
         End Do
      End If
      Call FZero(F_xc2,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc2,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipdRx,iGrid) = (One+Two*Delta)*Rho(ipdRx,iGrid)
            Rho_Temp(ipdRy,iGrid) = (One+Two*Delta)*Rho(ipdRy,iGrid)
            Rho_Temp(ipdRz,iGrid) = (One+Two*Delta)*Rho(ipdRz,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipdRxa,iGrid) = (One+Two*Delta)*Rho(ipdRxa,iGrid)
            Rho_Temp(ipdRya,iGrid) = (One+Two*Delta)*Rho(ipdRya,iGrid)
            Rho_Temp(ipdRza,iGrid) = (One+Two*Delta)*Rho(ipdRza,iGrid)
         End Do
      End If
      Call FZero(F_xc3,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc3,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipdRx,iGrid) = (One-Two*Delta)*Rho(ipdRx,iGrid)
            Rho_Temp(ipdRy,iGrid) = (One-Two*Delta)*Rho(ipdRy,iGrid)
            Rho_Temp(ipdRz,iGrid) = (One-Two*Delta)*Rho(ipdRz,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipdRxa,iGrid) = (One-Two*Delta)*Rho(ipdRxa,iGrid)
            Rho_Temp(ipdRya,iGrid) = (One-Two*Delta)*Rho(ipdRya,iGrid)
            Rho_Temp(ipdRza,iGrid) = (One-Two*Delta)*Rho(ipdRza,iGrid)
         End Do
      End If
      Call FZero(F_xc4,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc4,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      Do iGrid = 1, mGrid
         If (iSpin.eq.1) Then
            If (Tf*Rho(ipR,iGrid).lt.T_x) Go To 199
            Gaa=Rho(ipdRx,iGrid)**2
     &         +Rho(ipdRy,iGrid)**2
     &         +Rho(ipdRz,iGrid)**2
            Gaa1o2=Sqrt(Gaa)
            Gab=Gaa
            V_a=dF_dRho(ipGxx,iGrid)*Two*Gaa1o2
     &         +dF_dRho(ipGxy,iGrid)*Gab/Gaa1o2
            Delta0=Delta*Sqrt(Gaa)*Two
         Else
            If (Tf*Rho(ipRa,iGrid).lt.T_x) Go To 199
            Gaa=Rho(ipdRxa,iGrid)**2
     &         +Rho(ipdRya,iGrid)**2
     &         +Rho(ipdRza,iGrid)**2
            Gab=Rho(ipdRxa,iGrid)*Rho(ipdRxb,iGrid)
     &         +Rho(ipdRya,iGrid)*Rho(ipdRyb,iGrid)
     &         +Rho(ipdRza,iGrid)*Rho(ipdRzb,iGrid)
            Gaa1o2=Sqrt(Gaa)
            V_a=dF_dRho(ipGaa,iGrid)*Two*Gaa1o2
     &         +dF_dRho(ipGab,iGrid)*Gab/Gaa1o2
            Delta0=Delta*Sqrt(Gaa)
         End If
         Diff1=F_xc1(iGrid)-F_xc2(iGrid)
         Diff2=F_xc3(iGrid)-F_xc4(iGrid)
         V_n=(Two*Diff1-Diff2/Four)/(Three*Delta0)
         Rel_Diff=Abs((V_a-V_n)/V_n)
         If (V_a-V_n.eq.Zero) Rel_Diff=Zero
C        If (V_a.eq.Zero) Rel_Diff=Abs(V_n)
         If (Abs(V_a).lt.Thr_V_a) Rel_Diff=Abs(V_n)
         If (Rel_Diff.gt.Error) Then
            Call WarningMessage(2,'Rel_Diff.gt.Error')
            Write (6,'(A,A)') 'Functional:',FName
            Write (6,*) 'Error in dF_dgradRho(alpha)!'
            Write (6,*) 'V_a,V_n=',V_a,V_n
            Write (6,*) 'iGrid,mGrid=',iGrid,mGrid
            Write (6,*) 'Diff1,Diff2=',Diff1,Diff2
            If (iSpin.eq.1) Then
               Write (6,*) 'Rho=',Rho(ipR,iGrid)
               Write (6,*) 'GradRho=',Rho(ipdRx,iGrid),
     &                                Rho(ipdRy,iGrid),
     &                                Rho(ipdRz,iGrid)
               Write (6,*) 'Tau=',    Rho(ipTau,iGrid)
            Else
               Write (6,*) 'Rho=',Rho(ipRa,iGrid),
     &                            Rho(ipRb,iGrid)
               Write (6,*) 'GradRho=',Rho(ipdRxa,iGrid),
     &                                Rho(ipdRya,iGrid),
     &                                Rho(ipdRza,iGrid)
               Write (6,*) '        ',Rho(ipdRxb,iGrid),
     &                                Rho(ipdRyb,iGrid),
     &                                Rho(ipdRzb,iGrid)
               Write (6,*) 'Tau=',    Rho(ipTaua,iGrid),
     &                                Rho(ipTaub,iGrid)
            End If
            Call Abend()
         End If
 199     Continue
      End Do
*
      If (iSpin.eq.1) Go To 200
*                                                                      *
************************************************************************
*                                                                      *
*     Check for grad rho (beta)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipdRxb,iGrid) = (One+Delta)*Rho(ipdRxb,iGrid)
         Rho_Temp(ipdRyb,iGrid) = (One+Delta)*Rho(ipdRyb,iGrid)
         Rho_Temp(ipdRzb,iGrid) = (One+Delta)*Rho(ipdRzb,iGrid)
      End Do
      Call FZero(F_xc1,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc1,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipdRxb,iGrid) = (One-Delta)*Rho(ipdRxb,iGrid)
         Rho_Temp(ipdRyb,iGrid) = (One-Delta)*Rho(ipdRyb,iGrid)
         Rho_Temp(ipdRzb,iGrid) = (One-Delta)*Rho(ipdRzb,iGrid)
      End Do
      Call FZero(F_xc2,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc2,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipdRxb,iGrid) = (One+Two*Delta)*Rho(ipdRxb,iGrid)
         Rho_Temp(ipdRyb,iGrid) = (One+Two*Delta)*Rho(ipdRyb,iGrid)
         Rho_Temp(ipdRzb,iGrid) = (One+Two*Delta)*Rho(ipdRzb,iGrid)
      End Do
      Call FZero(F_xc3,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc3,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipdRxb,iGrid) = (One-Two*Delta)*Rho(ipdRxb,iGrid)
         Rho_Temp(ipdRyb,iGrid) = (One-Two*Delta)*Rho(ipdRyb,iGrid)
         Rho_Temp(ipdRzb,iGrid) = (One-Two*Delta)*Rho(ipdRzb,iGrid)
      End Do
      Call FZero(F_xc4,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc4,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      Do iGrid = 1, mGrid
         If (Tf*Rho(ipRa,iGrid).lt.T_x) Go To 201
         Gbb=Rho(ipdRxb,iGrid)**2
     &      +Rho(ipdRyb,iGrid)**2
     &      +Rho(ipdRzb,iGrid)**2
         Gab=Rho(ipdRxa,iGrid)*Rho(ipdRxb,iGrid)
     &      +Rho(ipdRya,iGrid)*Rho(ipdRyb,iGrid)
     &      +Rho(ipdRza,iGrid)*Rho(ipdRzb,iGrid)
         Gbb1o2=Sqrt(Gbb)
         V_a=dF_dRho(ipGbb,iGrid)*Two*Gbb1o2
     &      +dF_dRho(ipGab,iGrid)*Gab/Gbb1o2
         Delta0=Delta*Sqrt(Gbb)
         Diff1=F_xc1(iGrid)-F_xc2(iGrid)
         Diff2=F_xc3(iGrid)-F_xc4(iGrid)
         V_n=(Two*Diff1-Diff2/Four)/(Three*Delta0)
         Rel_Diff=Abs((V_a-V_n)/V_n)
         If (V_a-V_n.eq.Zero) Rel_Diff=Zero
C        If (V_a.eq.Zero) Rel_Diff=Abs(V_n)
         If (Abs(V_a).lt.Thr_V_a) Rel_Diff=Abs(V_n)
         If (Rel_Diff.gt.Error) Then
            Call WarningMessage(2,'Rel_Diff.gt.Error')
            Write (6,'(A,A)') 'Functional:',FName
            Write (6,*) 'Error in dF_dgradRho(beta)!'
            Write (6,*) 'V_a,V_n=',V_a,V_n
            Write (6,*) 'iGrid,mGrid=',iGrid,mGrid
            Write (6,*) 'Rho=',Rho(ipRa,iGrid),
     &                         Rho(ipRb,iGrid)
            Write (6,*) 'GradRho=',Rho(ipdRxa,iGrid),
     &                             Rho(ipdRya,iGrid),
     &                             Rho(ipdRza,iGrid)
            Write (6,*) '        ',Rho(ipdRxb,iGrid),
     &                             Rho(ipdRyb,iGrid),
     &                             Rho(ipdRzb,iGrid)
            Write (6,*) 'Tau=',    Rho(ipTaua,iGrid),
     &                             Rho(ipTaub,iGrid)
            Call Abend()
         End If
 201     Continue
      End Do
*
 200  Continue
      If (FType.eq.GGA_type) Return
*                                                                      *
************************************************************************
*                                                                      *
*     Check for tau (alpha)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipTau,iGrid) = (One+Delta)*Rho(ipTau,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipTaua,iGrid) = (One+Delta)*Rho(ipTaua,iGrid)
         End Do
      End If
      Call FZero(F_xc1,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc1,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipTau,iGrid) = (One-Delta)*Rho(ipTau,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipTaua,iGrid) = (One-Delta)*Rho(ipTaua,iGrid)
         End Do
      End If
      Call FZero(F_xc2,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc2,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipTau,iGrid) = (One+Two*Delta)*Rho(ipTau,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipTaua,iGrid) = (One+Two*Delta)*Rho(ipTaua,iGrid)
         End Do
      End If
      Call FZero(F_xc3,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc3,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            Rho_Temp(ipTau,iGrid) = (One-Two*Delta)*Rho(ipTau,iGrid)
         End Do
      Else
         Do iGrid = 1, mGrid
            Rho_Temp(ipTaua,iGrid) = (One-Two*Delta)*Rho(ipTaua,iGrid)
         End Do
      End If
      Call FZero(F_xc4,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc4,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      Do iGrid = 1, mGrid
         If (iSpin.eq.1) Then
            If (Tf*Rho(ipR,iGrid).lt.T_x) Go To 299
            V_a=dF_dRho(ipT,iGrid)
            Delta0=Delta*Rho(ipTau,iGrid)*Two
         Else
            If (Tf*Rho(ipRa,iGrid).lt.T_x) Go To 299
            V_a=dF_dRho(ipTa,iGrid)
            Delta0=Delta*Rho(ipTaua,iGrid)
         End If
         Diff1=F_xc1(iGrid)-F_xc2(iGrid)
         Diff2=F_xc3(iGrid)-F_xc4(iGrid)
         V_n=(Two*Diff1-Diff2/Four)/(Three*Delta0)
         Rel_Diff=Abs((V_a-V_n)/V_n)
         If (V_a-V_n.eq.Zero) Rel_Diff=Zero
C        If (V_a.eq.Zero) Rel_Diff=Abs(V_n)
         If (Abs(V_a).lt.Thr_V_a) Rel_Diff=Abs(V_n)
         If (Rel_Diff.gt.Error) Then
            Call WarningMessage(2,'Rel_Diff.gt.Error')
            Write (6,'(A,A)') 'Functional:',FName
            Write (6,*) 'Error in tau(alpha)!'
            Write (6,*) 'V_a,V_n=',V_a,V_n
            If (iSpin.eq.1) Then
               Write (6,*) 'Rho=',Rho(ipR,iGrid)
               Write (6,*) 'GradRho=',Rho(ipdRx,iGrid),
     &                                Rho(ipdRy,iGrid),
     &                                Rho(ipdRz,iGrid)
               Write (6,*) 'Tau=',    Rho(ipTau,iGrid)
            Else
               Write (6,*) 'Rho=',Rho(ipRa,iGrid),
     &                            Rho(ipRb,iGrid)
               Write (6,*) 'GradRho=',Rho(ipdRxa,iGrid),
     &                                Rho(ipdRya,iGrid),
     &                                Rho(ipdRza,iGrid)
               Write (6,*) '        ',Rho(ipdRxb,iGrid),
     &                                Rho(ipdRyb,iGrid),
     &                                Rho(ipdRzb,iGrid)
               Write (6,*) 'Tau=',    Rho(ipTaua,iGrid),
     &                                Rho(ipTaub,iGrid)
            End If
            Call Abend()
         End If
 299     Continue
      End Do
*
      If (iSpin.eq.1) Go To 300
*                                                                      *
************************************************************************
*                                                                      *
*     Check for tau (beta)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipTaub,iGrid) = (One+Delta)*Rho(ipTaub,iGrid)
      End Do
      Call FZero(F_xc1,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc1,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipTaub,iGrid) = (One-Delta)*Rho(ipTaub,iGrid)
      End Do
      Call FZero(F_xc2,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc2,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipTaub,iGrid) = (One+Two*Delta)*Rho(ipTaub,iGrid)
      End Do
      Call FZero(F_xc3,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc3,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      call dcopy_(nRho*mGrid,Rho,1,Rho_Temp,1)
      Do iGrid = 1, mGrid
         Rho_Temp(ipTaub,iGrid) = (One-Two*Delta)*Rho(ipTaub,iGrid)
      End Do
      Call FZero(F_xc4,mGrid)
      Call Kernel(mGrid,Rho_Temp,nRho,P2_ontop,
     &            nP2_ontop,iSpin,F_xc4,
     &            dF_temp,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
*
      Do iGrid = 1, mGrid
         If (Tf*Rho(ipRa,iGrid).lt.T_x) Go To 301
         V_a=dF_dRho(ipTb,iGrid)
         Delta0=Delta*Rho(ipTaub,iGrid)
         Diff1=F_xc1(iGrid)-F_xc2(iGrid)
         Diff2=F_xc3(iGrid)-F_xc4(iGrid)
         V_n=(Two*Diff1-Diff2/Four)/(Three*Delta0)
         Rel_Diff=Abs((V_a-V_n)/V_n)
         If (V_a-V_n.eq.Zero) Rel_Diff=Zero
C        If (V_a.eq.Zero) Rel_Diff=Abs(V_n)
         If (Abs(V_a).lt.Thr_V_a) Rel_Diff=Abs(V_n)
         If (Rel_Diff.gt.Error) Then
            Call WarningMessage(2,'Rel_Diff.gt.Error')
            Write (6,'(A,A)') 'Functional:',FName
            Write (6,*) 'Error in tau(beta)!'
            Write (6,*) 'V_a,V_n=',V_a,V_n
            Write (6,*) 'Rho=',Rho(ipRa,iGrid),
     &                         Rho(ipRb,iGrid)
            Write (6,*) 'GradRho=',Rho(ipdRxa,iGrid),
     &                             Rho(ipdRya,iGrid),
     &                             Rho(ipdRza,iGrid)
            Write (6,*) '        ',Rho(ipdRxb,iGrid),
     &                             Rho(ipdRyb,iGrid),
     &                             Rho(ipdRzb,iGrid)
            Write (6,*) 'Tau=',    Rho(ipTaua,iGrid),
     &                             Rho(ipTaub,iGrid)
            Call Abend()
         End If
 301     Continue
      End Do
 300  Continue
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
