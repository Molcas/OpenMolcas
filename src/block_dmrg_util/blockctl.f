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
* Copyright (C) 2014, Naoki Nakatani                                   *
************************************************************************
      Subroutine BlockCtl(LW1,TUVX,IFINAL,IRST)
************************************************************************
*                                                                      *
*     Description                                                      *
*                                                                      *
*     This is an alternative of DavCtl for DMRG-CASSCF in DMRGCtl      *
*                                                                      *
************************************************************************
*                                                                      *
*     DMRG control section                                             *
*                                                                      *
*     calling arguments:                                               *
*     LW1     : Memory pointer to active Fock matrix                   *
*               array of real*8                                        *
*     TUVX    : array of real*8                                        *
*               two-electron integrals (tu!vx)                         *
*     IFINAL  : integer                                                *
*               termination flag                                       *
*     IRST    : integer                                                *
*               restart flag of DMRG                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     N. Nakatani, Hokkaido University, Japan, 2014                    *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Dimension LW1(*), TUVX(*)

      Integer iChMolpro(8)
      Character*3 Label

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='BLOCKCTL')
      Call qEnter(ROUTINE)

* Load symmetry info from RunFile
      Call Get_iScalar('NSYM',nIrrep)
      Call Get_iArray('Symmetry operations',iOper,nIrrep)
      Call Get_iScalar('Rotational Symmetry Number',iSigma)

* Get character table to convert MOLPRO symmetry format
      Call MOLPRO_ChTab(nSym,Label,iChMolpro)

* Convert orbital symmetry into MOLPOR format
      Call Getmem('OrbSym','Allo','Inte',lOrbSym,NAC)
      iOrb=1
      Do iSym=1,nSym
        Do jOrb=1,NASH(iSym)
          iWork(lOrbSym+iOrb-1)=iChMolpro(iSym)
          iOrb=iOrb+1
        End Do
      End Do
      lSymMolpro=iChMolpro(lSym)

      NRDM_ORDER=2
      If (NACTEL.EQ.1) NRDM_ORDER=1

* Default setting for the first iteraction
      JRST=IRST
      ThDMRG=THRE*1.0d1
      ThNoise=1.0d-4

      If (JRST.NE.0) Then
        If (IFINAL.EQ.0) Then
          If (ABS(ROTMAX).GE.1.0d-1) Then
* No restart
            JRST=0
            ThDMRG=THRE*1.0d1
            ThNoise=1.0d-4
          Else If (ABS(ROTMAX).GE.0.5d-1) Then
* Full-restart with noise
            JRST=2
            ThDMRG=THRE*5.0d0
            ThNoise=1.0d-5
          Else If (ABS(ROTMAX).GE.0.1d-1) Then
* Full-restart with smaller noise
            JRST=2
            ThDMRG=THRE
            ThNoise=1.0d-6
          Else
* Restart without noise
            JRST=1
            ThDMRG=THRE
            ThNoise=0.0d0
          End If
        Else If (IFINAL.EQ.1) Then
* Full-restart for CI-only
          JRST=2
          ThDMRG=THRE
          ThNoise=1.0d-6
        Else
* Full-restart for the final wavefunction
          If (iOrbTyp.EQ.2) Then
* with OutOrb = Canonical, this will be problematic for LMO-DMRG calc.
            JRST=2
            ThDMRG=THRE
            ThNoise=1.0d-6
          Else
* by default, just restarting
            JRST=1
            ThDMRG=THRE
            ThNoise=0.0d0
          End If
        End If
      End If

* In the following calls the extra array HFOCC has been passed
* from Molcas to Block to have a user customized reference det.
* This change must be done accordingly into the block code.
* In particular, in file molcas/block_calldmrg.C:
* line~23:  const char* hf_occ)
* line~24:  block_calldmrg(*Restart, *N_roots, *N_act, *N_elec, *M_s, Sym, *iSym, OrbSym, *E_core, h0, tuvx, *M_state, *N_pdm, *T_sweep, *T_noise, E_sweep, hf_occ);
* line~63:  const char* hf_occ)
* line~186: fcon << "hf_occ " << hf_occ << endl;
*
* In molcas/block_calldmrg.h:
* line~25:  const char* hf_occ);
* line~46:  const char* hf_occ);


* Compute DMRG
      Call block_calldmrg(JRST,lRoots,NAC,NACTEL,ISPIN-1,
     &                    Label,lSymMolpro,iWork(lOrbSym),
     &                    0.0d0,LW1,TUVX,MxDMRG,NRDM_ORDER,
     &                    ThDMRG,ThNoise,ENER(1,ITER),HFOCC,NRS2T)

      If (IFINAL.EQ.2 .AND. Do3RDM .AND. NACTEL.GT.2) Then
* Compute 3RDM for DMRG-cu4-CASPT2
        Call block_calldmrg(1,lRoots,NAC,NACTEL,ISPIN-1,
     &                      Label,lSymMolpro,iWork(lOrbSym),
     &                      0.0d0,LW1,TUVX,MxDMRG,3,
     &                      THRE,0.0d0,ENER(1,ITER),HFOCC,NRS2T)
      End If

      Call Getmem('OrbSym','Free','Inte',lOrbSym,NAC)

      Call qExit(ROUTINE)

      Return
      End
