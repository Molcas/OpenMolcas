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
          SubRoutine MkFck(iAnga,iCmpa,iCmp,
     &                        Shijij,
     &                        iShll,iShell,IndShl,
     &                        iBasi,jBasj,kBask,lBasl,
     &                        iAO,iAOst,nOp,jOp,
     &                        Dij,mDij,nDij,ij1,ij2,ij3,ij4,
     &                        Dkl,mDkl,nDkl,kl1,kl2,kl3,kl4,
     &                        Dik,mDik,nDik,ik1,ik2,ik3,ik4,
     &                        Dil,mDil,nDil,il1,il2,il3,il4,
     &                        Djk,mDjk,nDjk,jk1,jk2,jk3,jk4,
     &                        Djl,mDjl,nDjl,jl1,jl2,jl3,jl4,
     &                        AOInt,nAO,TwoHam,nFock,
     &                        Scrtch1,nS1,Scrtch2,nS2,
     &                        iDCRR,iDCRS,iDCRT,FckTmp,nFT,
     &                        pert,iuvwx,iCent,iCar,indgrd,ipDisp)
*
************************************************************************
*                                                                      *
* Object: Driver for the generation of the two electron contribution   *
*         to the Fock Matrix directly from the two electron integrals. *
*                                                                      *
* Called from: TwoEL                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy   (ESSL)                                          *
*              ICopy                                                   *
*              Trnsps                                                  *
*              Trns1                                                   *
*              Phase                                                   *
*              FckAcc                                                  *
*              QExit                                                   *
*                                                                      *
*              Anders Bernhardsson 1995                                *
************************************************************************
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
c#include "print.fh"
#include "disp.fh"
#include "disp2.fh"
*
      Real*8 Dij(mDij,nDij),Dkl(mDkl,nDkl),Dik(mDik,nDik),
     &       Dil(mDil,nDil),Djk(mDjk,nDjk),Djl(mDjl,nDjl),
     &       FckTmp(nFT),AOInt(nAO),TwoHam(nFock),
     &       Scrtch1(nS1),Scrtch2(nS2)
      Integer iCmp(4), nOp(4),iAnga(4), iShll(4),iShell(4),
     &        jOp(6),iCmpa(4) , iAO(4), iAOst(4), IndShl(4),
     &        indgrd(3,4,0:nirrep-1),ipdisp(*)
      Logical Shijij,pert(0:nIrrep-1)
*
*     Just the make a nice interface
*
c     iRout = 12
c     iPrint = nPrint(iRout)
      nijkl=iBasi*jBasj*kBask*lBasl
*
*--------------Accumulate contributions directly to the symmetry
*              adapted Fock matrix.
*
      Fact=DBLE(iuvwx)/DBLE(nIrrep)
*
      Call FckAcc_mck(iAnga,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &            Shijij,iShll,iShell,IndShl,nOp,nijkl,
     &            AOInt,TwoHam,nFock,Scrtch2,nS2,
     &            iAO,iAOst,
     &            iBasi,jBasj,kBask,lBasl,
     &            Dij(1,jOp(1)),ij1,ij2,ij3,ij4,
     &            Dkl(1,jOp(2)),kl1,kl2,kl3,kl4,
     &            Dik(1,jOp(3)),ik1,ik2,ik3,ik4,
     &            Dil(1,jOp(4)),il1,il2,il3,il4,
     &            Djk(1,jOp(5)),jk1,jk2,jk3,jk4,
     &            Djl(1,jOp(6)),jl1,jl2,jl3,jl4,
     &            FckTmp,nFT,fact,iCar,iCent,pert,indgrd,ipdisp)
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iCmpa)
         Call Unused_real_array(Scrtch1)
         Call Unused_integer(iDCRR)
         Call Unused_integer(iDCRS)
         Call Unused_integer(iDCRT)
      End If
      End
