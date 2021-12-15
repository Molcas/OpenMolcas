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
* Copyright (C) 1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine TwoDns(ianga,iCmp,shijij,ishll,ishell,iAO,
     &           nop,iBasi,jBasj,kBask,lBasl,
     &           Aux,nAux,Work2,nWork2,Work3,nWork3,work4,
     &           nWork4,PSO,nPSO,Fact)
*
      use Basis_Info
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "disp.fh"
#include "disp2.fh"
      Logical Shijij
      Integer nOp(4),iAnga(4),iCmp(4),iShell(4),iShll(4),iAO(4)
      Real*8  PSO(nPSO),Aux(nAux),Work2(nWork2),Work3(nWork3),
     &        Work4(nWork4)
*
*     Statement function
*
      nElem(i) = (i+1)*(i+2)/2
*
      nijkl=iBasi*jBasj*kBask*lBasl
      iShlla = iShll(1)
      jShllb = iShll(2)
      kShllc = iShll(3)
      lShlld = iShll(4)
      la=iAnga(1)
      lb=iAnga(2)
      lc=iAnga(3)
      ld=iAnga(4)
      iCmpa = iCmp(1)
      jCmpb = iCmp(2)
      kCmpc = iCmp(3)
      lCmpd = iCmp(4)
      mab = nElem(la)*nElem(lb)
      mcd = nElem(lc)*nElem(ld)

*
*----------------------------------------------------------------*
*
*              Fix the second order density matrix
*
*----------------------------------------------------------------*
*
*--------------Desymmetrize the second order density matrix
*
*--------------(faA fbR(B) | fcT(C) fdTS(D))ijkl
*     PSO->Work2
*
               Call DesymP(iAnga,iCmp(1),iCmp(2),
     &                     iCmp(3),iCmp(4),
     &                     Shijij,iShll,iShell,iAO,nOp,nijkl,
     &                     Aux,nAux,Work2,PSO,nPSO)
*
*
               If (Fact.ne.One) Call DScal_(nijkl*
     &             iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4),
     &             Fact,Work2,1)
*
*--------------Backtransform 2nd order density matrix from spherical
*              harmonic gaussians to cartesian gaussians.
*
               ijklab = nijkl * iCmp(1)*iCmp(2)
*  Work2->Work2  (Work3:Scratch)
               Call SphCr1(Work2,ijklab,
     &                     Work3,nWork3,
     &                     RSph(ipSph(lc)),nElem(lc),kCmpc,
     &                     Shells(kShllc)%Transf,
     &                     Shells(kShllc)%Prjct,
     &                     RSph(ipSph(ld)),nElem(ld),lCmpd,
     &                     Shells(lShlld)%Transf,
     &                     Shells(lShlld)%Prjct,
     &                     Work2,mcd)
*  Work2->Work4  (Work3:Scratch)
               Call SphCr2(Work2,nijkl,mcd,
     &                     Work3,nWork3,
     &                     RSph(ipSph(la)),nElem(la),iCmpa,
     &                     Shells(iShlla)%Transf,
     &                     Shells(iShlla)%Prjct,
     &                     RSph(ipSph(lb)),nElem(lb),jCmpb,
     &                     Shells(jShllb)%Transf,
     &                     Shells(jShllb)%Prjct,
     &                     Work4,mab)
*
*----------------------------------------------------------------*
*
*   P is now in cartisan AO base
*
      Return
      End
