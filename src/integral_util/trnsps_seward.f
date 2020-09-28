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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      Subroutine Trnsps_Seward(ijCmp, iCmp, jCmp, iAng, jAng, iShll,
     &                         jShll, kOp, ijkl, ij, AOInt, Scrtch)
************************************************************************
*                                                                      *
*  Object: to transpose the integrals in order to resolve the          *
*          redundancy (faA,fbB)=(fcC,fdD). In this case both sides will*
*          have the same DCR, i.e. (R)=(S). In this case we will only  *
*          need the unique combinations. For the off diagonal comb-    *
*          inations (R=/=S) we will pick up two terms. The terms are   *
*          (faA,fbR(B)|faT(A),fbTS(B)) and                             *
*          (faA,fbS(B)|faT(A),fbTR(B)). Since T and T-1 are the same   *
*          in D2h it is simple to see that after applying T on the     *
*          second integral we will end up with the first one.          *
*                                                                      *
*          However, since we compute the integrals in batches there    *
*          will not be a simple one to one correspondes between the    *
*          integrals in batch one and two. But after transposing the   *
*          pair arguments we will achive that one to one correspondens.*
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             May '90                                                  *
************************************************************************
      use Basis_Info
      use Real_Spherical, only: iSphCr
      use Symmetry_Info, only: iChBas
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 AOInt(ijkl,ijCmp,ijCmp), Scrtch(ijkl,ijCmp,ijCmp)
*
*     Statement Function
*
      iOff(ixyz)  = ixyz*(ixyz+1)*(ixyz+2)/6
*
      iRout = 67
      iPrint = nPrint(iRout)
*     Call RecPrt(' In Trnsps: AOInt ',' ',AOInt,ijkl,ijCmp*ijCmp)
*
*     Change phase factor. This is only nessecary if T=/=E.
*
      If (kOp.eq.0 .or. ijCmp.eq.0) Go To 14
      ii = iOff(iAng)
      jj = iOff(jAng)
      Do 10 i1 = 1, iCmp
       iChBs = iChBas(ii+i1)
       If (Shells(iShll)%Transf) iChBs = iChBas(iSphCr(ii+i1))
       pa1T = DBLE(iPrmt(kOp,iChBs))
       Do 11 i2 = 1, jCmp
        jChBs = iChBas(jj+i2)
        If (Shells(jShll)%Transf) jChBs = iChBas(iSphCr(jj+i2))
        pb1T = DBLE(iPrmt(kOp,jChBs))
        ij1 = iCmp*(i2-1)+i1
*
        Do 12 i3 = 1, iCmp
         kChBs = iChBas(ii+i3)
         If (Shells(iShll)%Transf) kChBs = iChBas(iSphCr(ii+i3))
         pa2T = DBLE(iPrmt(kOp,kChBs))
         Do 13 i4 = 1, jCmp
          lChBs = iChBas(jj+i4)
          If (Shells(jShll)%Transf) lChBs = iChBas(iSphCr(jj+i4))
          pb2T = DBLE(iPrmt(kOp,lChBs))
          ij2 = iCmp*(i4-1)+i3
          Factor=pa1T*pb1T*pa2T*pb2T
          If (Factor.ne.One) Call DScal_(ijkl,Factor,AOInt(1,ij1,ij2),1)
 13      Continue
 12     Continue
 11    Continue
 10   Continue
 14   Continue
*
*     Transpose ijkl,abcd to klij,cdab
*
      If (ijCmp.eq.1 .or. ij.eq.1) Then
         Call DGeTMI(AOInt,ijCmp*ij,ijCmp*ij)
      Else
         Do 100 i12 = 1, ijCmp
            Do 200 i34 = 1, ijCmp
               Call DGeTMO(AOInt(1,i12,i34),ij,ij,
     &                     ij,Scrtch(1,i34,i12),ij)
*
 200        Continue
 100     Continue
         call dcopy_(ijkl*ijCmp*ijCmp,Scrtch,1,AOInt,1)
      End If
*
*     Call RecPrt(' Exit Trnsps: AOInt ',' ',AOInt,ijkl,ijCmp*ijCmp)
*     Call GetMem(' Exit Trnsps','CHECK','REAL',iDum,iDum)
      Return
      End
