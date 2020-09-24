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
      SubRoutine Cho_GnVc_GetInt(xInt,lInt,nVecRS,iVecRS,ListSp,
     &                           mSym,mPass,mmShl,iPass1,NumPass,NumSP)
C
C     Purpose: compute integrals for NumPass integral passes starting at
C              pass iPass1.
C
#include "implicit.fh"
      Real*8  xInt(lInt)
      Integer nVecRS(mSym,mPass), iVecRS(mSym,mPass), ListSP(mmShl)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Character*15 SecNam
      Parameter (SecNam = 'Cho_GnVc_GetInt')

      Parameter (N2 = InfVec_N2)

      Integer  Cho_F2SP
      External Cho_F2SP

      InfVec(i,j,k)=iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
      IndRSh(i)=iWork(ip_IndRSh-1+i)

#if defined (_DEBUG_)
#endif

C     Initialization and input check.
C     -------------------------------

      If (NumPass .lt. 1) Then
         NumSP = 0
         Go To 1 ! exit
      End If

      If (mSym .ne. nSym) Then
         Call Cho_Quit('Input error [1] in '//SecNam,103)
      End If

      If (iPass1 .lt. 1) Then
         Call Cho_Quit('Input error [2] in '//SecNam,103)
      End If

      iPass2 = iPass1 + NumPass - 1
      If (iPass2 .gt. mPass) Then
         Call Cho_Quit('Input error [3] in '//SecNam,103)
      End If

      If (mmShl .lt. nnShl) Then
         Call Cho_Quit('Input error [4] in '//SecNam,103)
      End If

C     Set up list of shell pairs to compute.
C     --------------------------------------

      l_SPTmp = nnShl
      Call Cho_Mem('SPTmp','Allo','Inte',ip_SPTmp,l_SPTmp)
      Call Cho_iZero(iWork(ip_SPTmp),l_SPTmp)
      kOff0 = ip_SPTmp - 1
      NumSP = 0
      Do iPass = iPass1,iPass2
         Do iSym = 1,nSym
            iV1 = iVecRS(iSym,iPass)
            iV2 = iV1 + nVecRS(iSym,iPass) - 1
            Do iV = iV1,iV2
               iAB   = InfVec(iV,1,iSym) ! addr in 1st reduced set
               jShAB = IndRsh(iAB) ! shell pair (full)
               iShAB = Cho_F2SP(jShAB) ! reduced shell pair
               If (iShAB .gt. 0) Then
                  If (iWork(kOff0+iShAB) .eq. 0) Then ! register SP
                     iWork(kOff0+iShAB) = 1
                     NumSP = NumSP + 1
                     ListSP(NumSP) = iShAB
                  End If
               Else
                  Call Cho_Quit('SP not found in reduced list!',103)
               End If
            End Do
         End Do
      End Do
      Call Cho_Mem('SPTmp','Free','Inte',ip_SPTmp,l_SPTmp)

C     Set memory used by Seward.
C     --------------------------

      Call GetMem('Int.Max','Max ','Real',kSewInt,lSewInt)
      Call xSetMem_Ints(lSewInt)

C     Loop through shell pair list.
C     -----------------------------

      Do iSP = 1,NumSP

C        Get shell pair index.
C        ---------------------

         iShAB = ListSP(iSP)

C        Compute integral distribution (**|A B).
C        ---------------------------------------

         Call Cho_MCA_CalcInt_3(xInt,lInt,iShAB)

      End Do

C     Deallocation.
C     -------------

      Call xRlsMem_Ints
    1 Continue
#if defined (_DEBUG_)
#endif

      End
