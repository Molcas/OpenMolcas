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
* Copyright (C) 1991,1993,1995,1999,2000, Roland Lindh                 *
************************************************************************
      Subroutine Eval_Ints_New(iiS,jjS,kkS,llS,TInt,nTInt,
     &                         iTOffs,nShi,nShj,nShk,nShl,
     &                         Integ_Proc,
     &                         Dens,Fock,lDens,ExFac,nDens,
     &                         Ind,nInd,FckNoClmb,FckNoExch)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals, parallel region          *
*          contains memory partitioning and loops over uncontracted    *
*          functions...                                                *
*                                                                      *
*  Input:                                                              *
*          iiS,jjS,kkS,llS     : shell indices                         *
*          TInt                : Computed Integrals                    *
*          nTInt               : dimension of TInt                     *
*          iTOffs              : iTOffs holds symmetry block offsets   *
*                                                                      *
*     nShi,nShj,          Dimensions used for blocks in Tint (input)   *
*     nshk,nshl:          Symmetry block isym,jsym,ksym,lsym for       *
*                         shells iS,jS,kS,lS starts at                 *
*                         iTOffs(ksym,jsym,isym)+1 and is dimensioned  *
*                         [nshl(lsym),nshk(ksym),nshj(jsym,nshi(isym)] *
*                         Note that l runs fastest! The dimensions     *
*                         must be larger or equal to the number of     *
*                         SAOs in the specified shells and symmetries, *
*                         otherwise chaos!!                            *
*                                                                      *
*          Dens                : 1-particle density matrix             *
*          lDens               : length of density/Fock matrices       *
*          nDens               : # of density/Fock matrices            *
*          ExFac               : another scaling factor passed to      *
*                                Integ_Proc                            *
*          Ind,nInd            : auxiliary index list for Fock matrix  *
*                                construction                          *
*          FckNoClmb           : no Coulomb contributions to Fock mat. *
*          FckNoExch           : no exchange contributions to Fock mat.*
*          Thize               : int threshold for disk write (SD)     *
*          W2Disc,PreSch       : booleans 1st iter / prescreening      *
*          Disc_Mx             : # ints to write on disk (semidirect)  *
*          iDisk               : act. position in file TMPINT          *
*                                                                      *
*  Output: Fock                : 2el Hamiltonian                       *
*          Disc                : # ints written to disk (semidirect)   *
*                                                                      *
*                                                                      *
*  Local:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter,QExit                                            *
*              Int_Setup                                               *
*              Dens_Info                                               *
*              MemRys                                                  *
*              PSOAO0                                                  *
*              Picky_                                                  *
*              TwoEl_NoSym                                             *
*              TwoEl_Sym                                               *
*              Integ_Proc                                              *
*                                                                      *
* Author:     Roland Lindh                                             *
*             Dept. of Theoretical Chemistry, University of Lund,      *
*             SWEDEN.                                                  *
*                                                                      *
*             Modified for k2 loop. August 1991                        *
*             Modified for direct SCF. January 1993                    *
*             Modified to minimize overhead for calculations with      *
*             small basis sets and large molecules. Sept. 1993         *
*             parallel region split off in drvtwo.f, April 1995        *
*             Total rehack May 1999                                    *
*             Wrapper with old parameter list, Jan 2000                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      External Integ_Proc
#include "real.fh"
*     subroutine parameters
      Integer lDens
      Real*8  Thize,Fock(lDens,nDens),Dens(lDens,nDens),
     &        ExFac(nDens), Disc_Mx,Disc, TInt(nTInt)
      Integer iTOffs(8,8,8), nShi(0:7), nShj(0:7), nShk(0:7), nShl(0:7),
     &        Ind(nInd,nInd,2)
      Logical W2Disc,PreSch,FckNoClmb(nDens),FckNoExch(nDens),
     &        DoIntegrals,DoFock
*                                                                      *
************************************************************************
*                                                                      *
*     Set the additional parameters for integrals only option.
*
*     Disable semi-direct option
*
      W2Disc=.False.
      Disc=Zero
      Thize=Zero
      Disc_Mx=Zero
*
      PreSch=.True.
      DoFock=.False.
      DoIntegrals=.True.
*
      nij=max(iis,jjs)*(max(iis,jjs)+1)/2+min(iis,jjs)
      nkl=max(kks,lls)*(max(kks,lls)+1)/2+min(kks,lls)
      pmax=dble(max(nij,nkl))
      pmin=dble(min(nij,nkl))
      Quad_ijkl=pmax*(pmax+One)/Two+Pmin
*                                                                      *
************************************************************************
*                                                                      *
*     Call to subroutine with extended parameter list.
*
      Call Eval_Ints_New_Internal
     &               (iiS,jjS,kkS,llS,TInt,nTInt,
     &                iTOffs,nShi,nShj,nShk,nShl,
     &                Integ_Proc,
     &                Dens,Fock,lDens,ExFac,nDens,
     &                Ind,nInd,FckNoClmb,FckNoExch,
     &                Thize,W2Disc,PreSch,Disc_Mx,Disc, ! New arguments
     &                Quad_ijkl,DoIntegrals,DoFock)     ! New arguments
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
