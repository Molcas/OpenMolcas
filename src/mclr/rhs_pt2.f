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
* Copyright (C) 1998, Anders Bernhardsson                              *
************************************************************************
#define NOCODE
      Subroutine RHS_PT2(rkappa,iprci)
#ifdef NOCODE
      Real*8 rKappa(*)
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(rkappa)
         Call Unused_integer(iprci)
      End If
#else
*
*    Calculate and reads in from CASPT2 the RHS needed
*    for the calculation of the Lagrangian multipliers
*    for CASPT2.
*
*    From the same Fock matrix the effective Fock matrix
*    needed for the renormalization term is  calculated.

*    Here is the structure, it is not debugged and one
*    needs to check the detail, but all "input" is there.
*
      Implicit Real*8(a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "glbbas_mclr.fh"
#include "Files_mclr.fh"
#ifndef _DEBUGED_
#include "detdim.fh"
#include "csfbas_mclr.fh"
#endif
      Real*8 rKappa(*)
      Real*8, Allocatable:: DCAS(:)
      Half=0.5d0
*
*     Read in a and b part of effective gradient from CASPT2
*
      If (imethod.eq.2)Then
      i=0
      Call Getmem('TEMPCI','ALLO','REAL',ipT,nconf1)
      Call Getmem('TEMPCI','ALLO','REAL',ipT2,nconf1)
      Call Getmem('TEMPKAP','ALLO','REAL',ipK,ndens2)
      Call dDaFile(LuPT2,2,Work(ipK),ndens2,i)
#ifdef _DEBUGED_
*
*    CASPT2 mode
*
      Call dDaFile(LuPT2,2,Work(ipT),nconf1,i)
#else
*
* lucia mode
*
      n=nint(xispsm(State_SYM,1))
      Call dDaFile(LuPT2,2,Work(ipT),n,i)
#endif
*
*---- Transform from split GUGA to symmetric group
*
#ifdef _DEBUGED_
*
*    CASPT2 mode (split graph)
*
      Call Gugactl_MCLR(ipT,1)
#else
*
* lucia mode (Symmetric group)
*
      iprdia=0
      Call INCSFSD(STATE_SYM,STATE_SYM,.false.)
      CALL CSDTVC_MCLR(Work(ipT2),Work(ipT),2,
     &                 WORK(KDTOC),iWORK(KICTS(1)),
     &                 State_SYM,1,IPRDIA)
#endif
*
      Else
*
*     MP2
*
        i=0
        Call Getmem('TEMPKAP','ALLO','REAL',ipK,ndens2)
        Call dDaFile(LuPT2,2,Work(ipK),ndens2,i)
*       Call ThreeP(Kappa)
*       MP2
      End if
* ----
*
      Call GetMem('Temp1','ALLO','REAL',ipT1,ndens2)
      Call GetMem('Temp2','ALLO','REAL',ipT2,ndens2)
      Call GetMem('FockAO1','ALLO','REAL',ipFAO1,ndens2)
      Call GetMem('FockAO2','ALLO','REAL',ipFAO2,ndens2)
      Call GetMem('FockMO1','ALLO','REAL',ipFMO1,ndens2)
      Call GetMem('FockMO2','ALLO','REAL',ipFMO2,ndens2)
      Call GetMem('dens1','ALLO','REAL',ipDP,ndens2)
      Call GetMem('Dens2','ALLO','REAL',ipDP2,ndens2)
      Call GetMem('Dens4','ALLO','REAL',ipDCAS2,ndens2)
      Call GetMem('Kappa1','ALLO','REAL',ipK1,ndens2)
      Call GetMem('Kappa2','ALLO','REAL',ipK2,ndens2)
      Call GetMem('Dens5','ALLO','REAL',ipD,ndens2)
*
*---  Read in necessary densities.
*
      Call Qpg_dArray('D1ao',Found,nDens)
      If (Found .and. nDens/=0) Then
         Call mma_allocate(DCAS,nDens,Label='DCAS*)
      Else
         Write (6,*) 'RHS_PT2: Density not found'
         Call Abend()
      End If
      Call Get_D1ao(DCAS,nDens)
      irc=-1
      iopt=0
      Call RdRlx(irc,iopt,'D1PT22',Work(ipDP))
      If (irc.ne.0) Goto 100
      irc=-1
      iopt=0
      Call RdRlx(irc,iopt,'OVLP',rovlp)
      If (irc.ne.0) Goto 100
*
*--- Squared density
*
      Call UnFold_MCLR(Work(ipDP),Work(ipDP2))
*
      Call UnFold_MCLR(DCAS,Work(ipDCAS2))
*
*==============================================================================*
*
*
*---  Make Fock matrixes
*
      ExFac=1.0D0
*
*  1) P(PT2)
*
      nFlt=0
      nBMX=0
      Do iSym = 1, nSym
         nFlt=nFlt+nBas(iSym)*(nBas(iSym)+1)/2
         nBMX=Max(nBMX,nBas(iSym))
      End Do
*
      Call FZero(Work(ipFAO1),nDens2)
      Call FockTwo_Drv(nSym,nBas,nBas,nSkip,
     &                 Work(ipDP),Work(ipP2),Work(ipFAO1),nFlt,
     &                 ExFac,nDens2,nBMX)
      Call AO2MO(Work(ipFAO1),Work(ipFMO1))
*
*  2) P(CAS)
*
      Call FZero(Work(ipFAO2),nDens2)
      Call FockTwo_Drv(nSym,nBas,nBas,nSkip,
     &                 Work(ipDAS),Work(ipDCAS2),Work(ipFAO2),nFlt,
     &                 ExFac,nDens2,nBMX)
      Call AO2MO(Work(ipFAO2),Work(ipFMO2))
*
*-- Add one particle hamiltonian here ???
*
  ??? Call DaXpY_(ndens2,-rovlp,Work(kint1),1,Work(ipFMO1),1)
*
*==============================================================================*
*
*---- CI Vector
*
*     <i|Sigma> = <i|F|0> - <0|F|0><i|0>+CI_a+CI_b
*
      Call CISigma(0,State_sym,State_sym,Work(ipFMO1),Work(0),0,ipci,iprci,'N')
      rE=ddot_(nconf1,Work(ipin(iprci)),1,Work(ipin(ipci)),1)
      Call Daxpy_(nconf1,1.0d0,Work(ipT),1,Work(ipin(iprci)),1)
      Call Daxpy_(nconf1,-rE,Work(ipin(ipCI)),1,Work(ipin(iprci)),1)
*==============================================================================*
*
*     D^CAS*P(PT2+CAS)
*
      Call DFock(Work(ipDCAS2),Work(ipFMO1),Work(ipK1))
*
*     D^(PT2+CAS)*P(CAS)
*
      Call DFock(Work(ipDTot2),Work(ipFMO2),Work(ipK2))
*
*==============================================================================*
*
*     Add all orbital terms together
*
      Call Daxpy_(nDens2,1.0d0,Work(ipK2),1,Work(ipK1),1)
      call daxpy_(ndens2,1.0d0,Work(ipK),1,Work(ipK1),1)
*
*==============================================================================*
*
*  OKIDOKI Two things are needed, first of all the symmetric fockmatrix for the
*  connection term:
*
*==============================================================================*
*
*---  Calculate efficent Fock matrix in AO basis (contravariant,folded)
*     and write to disk. Woba
*
      Do iS=1,nSym
        If (nbas(is).ne.0)
     *  Call DGEadd(Work(ipK1-1+ipMat(is,is)),nBas(is),'N',
     *              Work(ipK1-1+ipMat(is,is)),nBas(is),'T',
     *              Work(ipK2-1+ipMat(is,is)),nBas(is),
     *              nBas(is),nBas(is))
      End Do
      Do iS=1,nSym
        If (nBas(is).ne.0) Then
           Call DGEMM_('N','N',
     &                 nBas(iS),nBas(iS),nBas(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 Work(ipK2-1+ipMat(is,is)),nBas(iS),
     &                 0.0d0,Work(ipT2),nBas(iS))
           Call DGEMM_('N','T',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Work(ipT2),nBas(iS),
     &                 Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 0.0d0,Work(ipFAO1-1+ipMat(iS,iS)),nBas(is))
        End If
      End Do
      Call Fold2(nsym,nbas,Work(ipFAO1),Work(ipFAO2))
      Call Put_Fock_Occ(Work(ipFAO2),nDens2)
*
*     And then I want the unsymmetric part to the MCLR
*
*
*---  Keep gradient terms woba!!!!!!!
*
      Do iS=1,nSym
        If (nbas(is).ne.0)
     *  Call DGESUB(Work(ipK1-1+ipMat(is,is)),nBas(is),'N',
     *              Work(ipK1-1+ipMat(is,is)),nBas(is),'T',
     *              rKappa(ipMat(is,is)),nBas(is),
     *              nBas(is),nBas(is))
      End Do

*
*    OK Thats all folks!!
*
      Call GetMem('Temp1','FREE','REAL',ipT1,ndens2)
      Call GetMem('Temp2','FREE','REAL',ipT2,ndens2)
      Call GetMem('FockAO1','FREE','REAL',ipFAO1,ndens2)
      Call GetMem('FockAO2','FREE','REAL',ipFAO2,ndens2)
      Call GetMem('FockMO1','FREE','REAL',ipFMO1,ndens2)
      Call GetMem('FockMO2','FREE','REAL',ipFMO2,ndens2)
      Call GetMem('dens1','FREE','REAL',ipDTot,ndens2)
      Call GetMem('Dens2','FREE','REAL',ipDTOT2,ndens2)
      Call mma_deallocate(DCAS)
      Call GetMem('Dens4','FREE','REAL',ipDCAS2,ndens2)
      Call GetMem('Kappa1','FREE','REAL',ipK1,ndens2)
      Call GetMem('Kappa2','FREE','REAL',ipK2,ndens2)
      Call GetMem('Kappa','FREE','REAL',ipK,ndens2)
      Call GetMem('Dens5','FREE','REAL',ipD,ndens2)
      Call Getmem('TEMPCI','FREE','REAL',ipT,nconf1)
*
*............
      return
 100  Call SysHalt('rhs_pt2')
      end
      Subroutine DFock(DAO,FockMO,Fock)
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "WrkSpc.fh"
      Real*8 DAO(*),FockMO(*),Fock(*)
      Call Getmem('Temp1','ALLO','REAL',ipT1,ndens2)
      Call Getmem('Temp2','ALLO','REAL',ipT2,ndens2)
      Do iS=1,nSym
        If (nBas(is).ne.0) Then
           Call DGEMM_('T','N',
     &                 nBas(iS),nBas(iS),nBas(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 DAO(ipCM(is)),nBas(iS),
     &                 0.0d0,Work(ipT2),nBas(iS))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Work(ipT2),nBas(iS),
     &                 Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 0.0d0,Work(ipT1),nBas(is))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Work(ipT1),nBas(iS),
     &                 FockMO(ipCM(iS)),nBas(iS),
     &                 0.0d0,Fock(ipCM(is)),nBas(is))
        End If
      End Do
      Call Getmem('Temp1','FREE','REAL',ipT1,ndens2)
      Call Getmem('Temp2','FREE','REAL',ipT2,ndens2)
      Return
      End
      Subroutine AO2MO(FAO ,FMO)
      Implicit Real*8 (a-h,o-z)
      Real*8 FAO(*),FMO(*)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Call GetMem('Temp','ALLO','REAL',ipT1,ndens2)
      Call GetMem('Temp','ALLO','REAL',ipT2,ndens2)
      ip=1
      Do iS=1,nSym
        If (nBas(is).ne.0) Then
           Call Square(FAO(ip),
     *                   Work(ipT1),
     *                   1,nBas(is),nBas(is))
           Call DGEMM_('T','N',
     &                 nBas(iS),nBas(iS),nBas(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 Work(ipT1),nBas(iS),
     &                 0.0d0,Work(ipT2),nBas(iS))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Work(ipT2),nBas(iS),
     &                 Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 0.0d0,FMO(ipMat(iS,iS)),nBas(is))
           ip=ip+nBas(is)*(nBas(iS)+1)/2
        End If
      End Do
      Call GetMem('Temp','FREE','REAL',ipT1,ndens2)
      Call GetMem('Temp','FREE','REAL',ipT2,ndens2)
#endif
      Return
      End
