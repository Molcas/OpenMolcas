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
* Copyright (C) 2007,2008, Roland Lindh                                *
************************************************************************
      Subroutine Mk_RICD_Shells(Info,nInfo,DInf,nDInf)
************************************************************************
*                                                                      *
*    Objective: To generate aCD auxiliary basis sets on-the-fly.       *
*                                                                      *
* Called from: RdCtl_Seward                                            *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chem. Phys., Lund Univ., Sweden.  *
*                                                                      *
*             Final implementation for aCD and acCD auxiliary          *
*             basis sets developed while visiting N. Ferre' at the     *
*             Univ. of Provance (champus Univ. Paul Cezanne) in        *
*             Marseille, France, 20 March - 19 April.                  *
*                                                                      *
*             Modified to transform the auxiliary basis to a true      *
*             Cholesky basis set while on TACC 2008 conference in      *
*             Songjiang District, Shanghai, China, 23-27 Sept. 2008.   *
*                                                                      *
************************************************************************
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "SysDef.fh"
#include "real.fh"
#include "print.fh"
#include "status.fh"
#include "stdalloc.fh"
      Real*8 DInf(nDInf)
      Logical DoRys, Save_Logical, W2L
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
       iPrint=49
C      iPrint=99
#else
       iPrint=5
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call Gen_RelPointers(-(Info-1))
      Call qEnter('Mk_aCD_Shells')
*
      Call StatusLine('Gateway:',
     &                ' Generating aCD or acCD auxiliary basis set')
*                                                                      *
************************************************************************
*                                                                      *
*     Preamble: Compute kOffAO  and lOffAO
*
      Call Setup_OffAO()
*
*     Set up transformation matrix from Cartesian to real spherical
*     harmonics.
*
      Call Sphere(iAngMx)
*
*     Setup of tables for coefficients for the Rys roots and weights.
*
      nDiff=0
      If (iAngMx.eq.0) nDiff=2
      DoRys=.True.
      Call SetUp_RW(DoRys,nDiff)
*
      iShll=Mx_Shll - 1
      iStrt=ipExp(iShll)
      mCnttp=nCnttp
*                                                                      *
************************************************************************
*                                                                      *
*     Add the DUMMY SHELL!
*
      Call Mk_Dummy_Shell(Info,nInfo,DInf,nDInf)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Loop now over all unique valence basis sets and generate the
*     corresponding aCD auxiliary basis sets. Note that there are two
*     different types of aCD auxiliary basis sets, aCD and acCD.
*
      Call Gen_RelPointers(Info-1)
      Do 1100 iCnttp = 1, mCnttp
         If (FragCnttp(iCnttp).or.nVal_Shells(iCnttp).eq.0) goto 1100
#ifdef _DEBUG_
         If (iPrint.ge.99)
     &   Write (6,*) 'Generating auxiliary basis set for valence basis'
     &             //':',iCnttp
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        Procrastinate the printing of the RICD basis set to library
*        until the last unique valence basis set is processed.
*
         W2L=.True.
         Do jCnttp = iCnttp+1, mCnttp
            If (Bsl_Old(iCnttp).eq.Bsl_Old(jCnttp)) Then
               W2L=.False.
               Exit
            End If
         End Do
*                                                                      *
************************************************************************
*                                                                      *
         If (Do_nacCD_Basis) Then
            Do_acCD_Basis=.False.
*                                                                      *
************************************************************************
*                                                                      *
*           nacCD section
*
*           Creat first a virgin aCD auxiliary basis set
*
            Thrshld_CD_Save = Thrshld_CD
            Thrshld_CD = Zero
            Save_Logical = Skip_High_AC
            Skip_High_AC = .False.
*
            kCnttp = nCnttp
            Call Mk_aCD_acCD_Shells(Info,nInfo,iCnttp,W2L)
            lCnttp = nCnttp
*
*           Now let us use the aCD auxiliary basis set to generate the
*           nacCD auxiliary basis set.
*
            Thrshld_CD = Thrshld_CD_Save
            Skip_High_AC = Save_Logical
            Call Mk_nacCD_Shells(Info,nInfo,kCnttp,lCnttp)
*
*           Remove the temporary aCD auxiliary basis set
*
            Do jCnttp = kCnttp+1, lCnttp
               Call rm_AuxShell(Info,nInfo,jCnttp)
            End Do
*                                                                      *
************************************************************************
*                                                                      *
         Else
*                                                                      *
************************************************************************
*                                                                      *
*        aCD and acCD section
*
*
            Call Mk_aCD_acCD_Shells(Info,nInfo,iCnttp,W2L)
*
         End If
*                                                                      *
************************************************************************
*                                                                      *
 1100 Continue ! iCnttp
      Call Gen_RelPointers(-(Info-1))
*                                                                      *
************************************************************************
*                                                                      *
*     Let us now Gram-Schmidt orthonormalize the auxiliary basis for
*     better numerics and balance.
*
*     Call ReNorm()
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
*                                                                      *
************************************************************************
*                                                                      *
*     Cleanup the mess!
*
      Call CloseR()
*
      If (Allocated(RSph)) Call mma_deallocate(RSph)
      If (Allocated(ipSph)) Call mma_deallocate(ipSph)
*                                                                      *
************************************************************************
*                                                                      *
      Call Gen_RelPointers(Info-1)
      Call qExit('Mk_aCD_Shells')
      Return
      End
      Subroutine Remove_High_Exponents(iD,nD,List2,mData,nTheta_All)
      Implicit Real*8 (a-h,o-z)
************************************************************************
*                                                                      *
*     Experimental code to be used with care.                          *
*                                                                      *
************************************************************************
#include "itmax.fh"
#include "info.fh"
      Integer iD(nD), List2(mData,nTheta_All)
      Logical Skip
*
      Call iVcPrt('Remove_High_Exponents: iD',' ',iD,nD)
      mD = nD
      i = 1
 100  Continue
         iTheta_All=iD(i)
         Skip=.False.
         kAng  = List2(1,iTheta_All)
         lAng  = List2(2,iTheta_All)
         k     = List2(5,iTheta_All)
         l     = List2(6,iTheta_All)
         kShll = List2(7,iTheta_All)
         lShll = List2(8,iTheta_All)
         If (kAng.eq.lAng) Then
            l     = List2(6,iTheta_All)
            Skip = (k.eq.1.and.l.eq.1).and.nExp(kShll).ne.1
         Else
            Skip=l.eq.1.and.nExp(lShll).ne.1
         End If
         If (Skip) Then
            If (mD.eq.i) Then
               mD = mD -1
               Go To 200
            End If
            Do j = i+1, mD
               iD(j-1) = iD(j)
            End Do
            mD = mD -1
            Go To 100
         End If
         i = i + 1
         If (i.le.mD) Go To 100
 200  Continue
      nD = mD
      Call iVcPrt('Remove_High_Exponents: iD',' ',iD,nD)
*
      Return
      End
      Subroutine Mk_AngList(iAL,nCompA,nCompB,
     &                      iD_c,nD_c,
     &                      List2,nList2,mData,
     &                      iAng,jAng)
      Integer iAL(nCompA,nCompB), iD_c(nD_c),
     &        List2(mData,nList2)
*
      Call IZero(iAL,nCompA*nCompB)
      Do jD_c = 1, nD_c
         ijSO=iD_c(jD_c)
         If (List2(1,ijSO).eq.iAng .and.
     &       List2(2,ijSO).eq.jAng ) Then
            iA = List2(3,ijSO)
            iB = List2(4,ijSO)
            iAL(iA,iB) = 1
         End If
      End Do
*
      Return
      End
