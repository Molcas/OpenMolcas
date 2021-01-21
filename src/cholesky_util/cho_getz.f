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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_GetZ(irc,
     &                    NVT,l_NVT,
     &                    nBlock,l_nBlock,
     &                    nV,l_nV1,l_nV2,
     &                    iV1,l_iV11,l_iV12,
     &                    ip_Z,l_Z1,l_Z2)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: get the Z vectors in core.
C
C     Input:
C        NVT(i): total no. of vectors in symmetry i
C        nBlock(i): no. of vector blocks in symmetry i
C        nV(j,i): number of vectors in block j of symmetry i
C        iV1(j,i): first vector in block j of symmetry i
C        ip_Z(j,i): pointer to triangular block j of symmetry i
C                   (here, j is a compound index j=iTri(k,l) for
C                    block k,l)
C
C     On exit, the Z vector blocks are stored in memory according
C     to ip_Z.
C
      use ChoSwp, only: InfVec
      Implicit None
      Integer irc
      Integer l_NVT
      Integer l_nBlock
      Integer l_nV1, l_nV2
      Integer l_iV11, l_iV12
      Integer l_Z1, l_Z2
      Integer NVT(l_NVT)
      Integer nBlock(l_nBlock)
      Integer nV(l_NV1,l_NV2)
      Integer iV1(l_IV11,l_iV12)
      Integer ip_Z(l_Z1,l_Z2)
#include "cholesky.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#if defined(_DEBUGPRINT_)
#include "choprint.fh"
#endif

      Integer iSym
      Integer iLoc, iRedC, iRed
      Integer ip_iRS2RS, l_iRS2RS
      Integer ip_Wrk, l_Wrk
      Integer idRS2RS, KK1, nVRead, mUsed, kOffV
      Integer KK, KKK, iJ
      Integer jBlock, kBlock
      Integer kOffZ
      Integer J_inBlock, K_inBlock

      Integer  Cho_iRange
      External Cho_iRange

      Character(LEN=8), Parameter:: SecNam='Cho_GetZ'

      Real*8 C0, C1, W0, W1

#if defined (_DEBUGPRINT_)
      Integer nBlock_Max, nnB, n
      Integer, Parameter:: myDebugInfo=100
      Real*8, Parameter:: Tol=1.0d-14
#endif

      Integer, Pointer:: InfVct(:,:,:)

      Integer i, j, k, iRS2RS, iTri

      iRS2RS(i)=iWork(ip_iRS2RS-1+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      Subroutine Cho_X_GetIP_InfVec(InfVcT)
      Integer, Pointer:: InfVct(:,:,:)
      End Subroutine Cho_X_GetIP_InfVec
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
C     Set return code.
C     ----------------

      irc = 0

#if defined (_DEBUGPRINT_)
      ! Check input variables
      If (l_NVT.lt.nSym .or. l_nBlock.lt.nSym .or.
     &    l_nV2.lt.nSym .or. l_iV12.lt.nSym .or.
     &    l_Z2.lt.nSym) Then
         irc=-1
         Return
      End If
      nBlock_Max=nBlock(1)
      Do iSym=2,nSym
         nBlock_Max=max(nBlock_Max,nBlock(iSym))
      End Do
      nnB=nBlock_Max*(nBlock_Max+1)/2
      If (l_nV1.lt.nBlock_Max .or. l_iV11.lt.nBlock_Max .or.
     &    l_Z1.lt.nnB) Then
         irc=-2
         Return
      End If
      If (iPrint.ge.myDebugInfo) Then
         Call Cho_Head('Entering '//SecNam//':','*',80,LuPri)
         Write(LuPri,'(A,8I6)') 'NVT   :',(NVT(iSym),iSym=1,nSym)
         Write(LuPri,'(A,8I6)') 'nBlock:',(nBlock(iSym),iSym=1,nSym)
      End If
      n=0
      Do iSym=1,nSym
         If (iPrint.ge.myDebugInfo) Then
            Do i=1,nBlock(iSym)
               Write(LuPri,'(2X,A,I6,A,I2,A,I6,A,I6)')
     &         'Block',i,' of sym.',iSym,': nV=',nV(i,iSym),
     &         ' iV1=',iV1(i,iSym)
            End Do
         End If
         Do j=1,nBlock(iSym)
            Do i=j,nBlock(iSym)
               If (i.eq.j) Then
                  k=nV(i,iSym)*(nV(i,iSym)+1)/2
               Else
                  k=nV(i,iSym)*nV(j,iSym)
               End If
               n=n+k
               If (iPrint.ge.myDebugInfo) Then
                  Write(Lupri,'(5X,A,I6,A,I6,A,I2,A,I9)')
     &            'iBlock',i,' jBlock',j,' of Sym.',iSym,': ip_Z=',
     &            ip_Z(iTri(i,j),iSym)
                  Write(LuPri,'(5X,A,I6,A,I9)')
     &            '--> Block dimension:',k,'  Block ends at:',
     &            ip_Z(iTri(i,j),iSym)+k-1
               End If
            End Do
         End Do
      End Do
      k=NVT(1)*(NVT(1)+1)/2
      Do iSym=2,nSym
         k=k+NVT(iSym)*(NVT(iSym)+1)/2
      End Do
      If (iPrint.ge.myDebugInfo) Then
         Write(Lupri,'(A,I8)') 'Total dimension of Z (from blocks):',n
         Write(Lupri,'(A,I8)') 'Total dimension of Z (from NVT)   :',k
         Write(Lupri,'(A,I8)') '                        Difference:',n-k
      End If
      If (n.ne.k) Then
         irc=-3
         If (iPrint.ge.myDebugInfo) Then
            Write(LuPri,'(A)') 'Ooops! They disagree....'
            Write(LuPri,'(A,I4)')
     &      'Returning with code:',irc
         End If
         Return
      End If
#endif

C     Zero result array.
C     ------------------

      Do iSym=1,nSym
         Do kBlock=1,nBlock(iSym)
            Call Cho_dZero(Work(ip_Z(iTri(kBlock,kBlock),iSym)),
     &                     nV(kBlock,iSym)*(nV(kBlock,iSym)+1)/2)
            Do jBlock=kBlock+1,nBlock(iSym)
               Call Cho_dZero(Work(ip_Z(iTri(jBlock,kBlock),iSym)),
     &                        nV(jBlock,iSym)*nV(kBlock,iSym))
            End Do
         End Do
      End Do

C     Scratch location in index arrays.
C     ---------------------------------

      iLoc=3 ! do NOT change (used implicitly by reading routine)

C     Get pointer to InfVec array for all vectors.
C     Needed for parallel runs.
C     --------------------------------------------

      Call Cho_X_GetIP_InfVec(InfVcT)

C     Copy rs1 to location 2.
C     -----------------------

      ! Note: location 2 must contain rs1 throughout this routine! It is
      ! an assumption which is never checked, so do NOT change this call
      ! or modify contents of location 2.
      Call Cho_X_RSCopy(irc,1,2)
      If (irc .ne. 0) Then
         Write(LuPri,'(A,A,I5)') SecNam,': Cho_X_RSCopy returned code',
     &                           irc
         irc = 1
         Go To 1 ! exit after deallocation
      End If

C     Get Z vectors.
C     --------------

      iRedC = -1
      Do iSym = 1,nSym

         l_iRS2RS = nnBstR(iSym,1)
         Call GetMem('RS-TO-RS','Allo','Inte',ip_iRS2RS,l_iRS2RS)
         Call mma_maxDBLE(l_Wrk)
         Call GetMem('Wrk','Allo','Real',ip_Wrk,l_Wrk)
         Call iZero(iWork(ip_iRS2RS),l_iRS2RS)
         idRS2RS = -2
         KK1 = 1
         Do While (KK1 .le. NumCho(iSym))
            nVRead = 0
            mUsed = 0
            Call Cho_Timer(C0,W0)
            Call Cho_X_VecRd(Work(ip_Wrk),l_Wrk,KK1,NumCho(iSym),iSym,
     &                       nVRead,iRedC,mUsed)
            If (nVRead .lt. 1) Then
               irc = 2
               Go To 1 ! exit after deallocation
            End If
            Call Cho_Timer(C1,W1)
            tDecom(1,2)=tDecom(1,2)+(C1-C0)
            tDecom(2,2)=tDecom(2,2)+(W1-W0)
            nSys_Call=nSys_Call+1
            kOffV = ip_Wrk - 1
            Do KKK = 0,nVRead-1
               KK = KK1 + KKK
               iRed = InfVec(KK,2,iSym)
               If (iRedC .ne. iRed) Then
                  Call Cho_X_SetRed(irc,iLoc,iRed)
                  If (irc .ne. 0) Then
                     irc = 3
                     Go To 1 ! exit after deallocation
                  End If
                  iRedC = iRed
               End If
               If (idRS2RS .ne. iRedC) Then
                  Call Cho_RS2RS(iWork(ip_iRS2RS),l_iRS2RS,
     &                           2,iLoc,iRedC,iSym)
                  idRS2RS = iRedC
               End If
               K = InfVec(KK,5,iSym)
               kBlock=Cho_iRange(K+1,iV1(1,iSym),nBlock(iSym),.True.)
#if defined (_DEBUGPRINT_)
               If (kBlock.lt.1 .or. kBlock.gt.nBlock(iSym)) Then
                  Call Cho_Quit('[BLOCK] Error detected in '//SecNam,
     &                          104)
               End If
#endif
               K_inBlock=K-iV1(kBlock,iSym)+1
               kOffZ=ip_Z(iTri(kBlock,kBlock),iSym)-1
               Do J_inBlock=K_inBlock,nV(kBlock,iSym)
                  J=iV1(kBlock,iSym)+J_inBlock-1
                  iJ=iRS2RS(InfVcT(J,1,iSym)-iiBstR(iSym,1))
                  Work(kOffZ+iTri(J_inBlock,K_inBlock))=Work(kOffV+iJ)
               End Do
               Do jBlock=kBlock+1,nBlock(iSym)
                  kOffZ=ip_Z(iTri(jBlock,kBlock),iSym)-1
     &                 +nV(jBlock,iSym)*(K_inBlock-1)
                  Do J_inBlock=1,nV(jBlock,iSym)
                     J=iV1(jBlock,iSym)+J_inBlock-1
                     iJ=iRS2RS(InfVcT(J,1,iSym)-iiBstR(iSym,1))
                     Work(kOffZ+J_inBlock)=Work(kOffV+iJ)
                  End Do
               End Do
               kOffV=kOffV+nnBstR(iSym,iLoc)
            End Do
            KK1=KK1+nVRead
         End Do
         Call GetMem('Wrk','Free','Real',ip_Wrk,l_Wrk)
         Call GetMem('RS-TO-RS','Free','Inte',ip_iRS2RS,l_iRS2RS)

      End Do

C     Synchronize result array across nodes.
C     --------------------------------------

      Do iSym=1,nSym
         Do kBlock=1,nBlock(iSym)
            Call Cho_GAdGOp(Work(ip_Z(iTri(kBlock,kBlock),iSym)),
     &                      nV(kBlock,iSym)*(nV(kBlock,iSym)+1)/2,'+')
            Do jBlock=kBlock+1,nBlock(iSym)
               Call Cho_GAdGOp(Work(ip_Z(iTri(jBlock,kBlock),iSym)),
     &                         nV(jBlock,iSym)*nV(kBlock,iSym),'+')
            End Do
         End Do
      End Do

#if defined (_DEBUGPRINT_)
C     Check that diagonal elements of Z are strictly positive.
C     --------------------------------------------------------

      If (iPrint.ge.myDebugInfo) Then
         Call Cho_Head(SecNam//': Diagonal of Z Vector Matrix','=',80,
     &                 LuPri)
         Write(Lupri,'(A)') ' '
      End If
      n=0
      Do iSym=1,nSym
         Do jBlock=1,nBlock(iSym)
            kOffZ=ip_Z(iTri(jBlock,jBlock),iSym)-1
            Do J_inBlock=1,nV(jBlock,iSym)
               If (iPrint.ge.myDebugInfo) Then
                  Write(LuPri,'(A,I2,A,I6,A,1P,D15.6,A,D15.6)')
     &            'Sym=',iSym,
     &            '  J=',iV1(jBlock,iSym)+J_inBlock-1,
     &            '  Z(J,J)=',
     &            Work(kOffZ+iTri(J_inBlock,J_inBlock)),
     &            '  Squared=',
     &            Work(kOffZ+iTri(J_inBlock,J_inBlock))**2
               End If
               If (abs(Work(kOffZ+iTri(J_inBlock,J_inBlock))).lt.Tol
     &             .or. Work(kOffZ+iTri(J_inBlock,J_inBlock)).lt.-Tol)
     &         Then
                  n=n+1
                  If (iPrint.ge.myDebugInfo) Then
                     Write(LuPri,'(A)')
     &               '  --> Small or negative Z diagonal!'
                  End If
               End If
            End Do
         End Do
      End Do
      If (n.ne.0) Then
         irc=20
         Go To 1 ! return
      End If
      ! Check diagonal elements
      Call Cho_CheckDiagFromZ(irc,NVT,l_NVT,nBlock,l_nBlock,
     &                        nV,l_nV1,l_nV2,
     &                        iV1,l_iV11,l_iV12,
     &                        ip_Z,l_Z1,l_Z2,
     &                        iPrint.ge.myDebugInfo)
      If (irc.ne.0) Then
         Go To 1 ! return
      End If
#endif

C     Exit. If error termination.
C     -----------------------------------------------

    1 Continue

#ifndef _DEBUGPRINT_
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(NVT)
#endif
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine for checking integral diagonal diagonal
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SubRoutine Cho_CheckDiagFromZ(irc,
     &                              NVT,l_NVT,
     &                              nBlock,l_nBlock,
     &                              nV,l_nV1,l_nV2,
     &                              iV1,l_iV11,l_iV12,
     &                              ip_Z,l_Z1,l_Z2,
     &                              Report)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: Check integral diagonal from Z vectors:
C
C              (J|J) = sum[K=1,J] Z(J,K)*Z(J,K)
C
C              Return codes
C              irc=0: all is fine
C              irc<0: too negative diagonals encountered, but otherwise
C                     calculation seems converged
C              irc>0: calculation not converged
C
      Implicit None
      Integer irc
      Integer l_NVT
      Integer l_nBlock
      Integer l_nV1, l_nV2
      Integer l_iV11, l_iV12
      Integer l_Z1, l_Z2
      Integer NVT(l_NVT)
      Integer nBlock(l_nBlock)
      Integer nV(l_NV1,l_NV2)
      Integer iV1(l_IV11,l_iV12)
      Integer ip_Z(l_Z1,l_Z2)
      Logical Report
#include "cholesky.fh"
#include "WrkSpc.fh"

      Character(LEN=18), Parameter:: SecNam='Cho_CheckDiagFromZ'

      Integer ip_D, l_D
      Integer iSym
      Integer jBlock, kblock
      Integer J_inBlock, K_inBlock
      Integer kOffD, kOffZ
      Integer iD
      Integer n1, n2, n3, n4, n5
      Integer nTot

      Real*8 Dmax, Damax, Dmin, Damin
      Integer, Pointer:: InfVct(:,:,:)

      Integer i, j
      Integer iTri
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      Subroutine Cho_X_GetIP_InfVec(InfVcT)
      Integer, Pointer:: InfVct(:,:,:)
      End Subroutine Cho_X_GetIP_InfVec
      End Interface
*                                                                      *
************************************************************************
*                                                                      *

      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      ! Get pointer to global InfVec array
      Call Cho_X_getIP_InfVec(InfVcT)

      ! Allocate memory for exact integral diagonal
      l_D=nnBstRT(1)
      Call GetMem('IntDia','Allo','Real',ip_D,l_D)

      ! Read diagonal
      Call Cho_IODiag(Work(ip_D),2)

      ! Subtract Z vector contributions
      kOffD=ip_D-1
      Do iSym=1,nSym
         Do kBlock=1,nBlock(iSym)
            Do K_inBlock=1,nV(kBlock,iSym)
               kOffZ=ip_Z(iTri(kBlock,kBlock),iSym)-1
               Do J_inBlock=K_inBlock,nV(kBlock,iSym)
                  J=iV1(kBlock,iSym)+J_inBlock-1
                  iD=InfVcT(J,1,iSym)
                  Work(kOffD+iD)=Work(kOffD+iD)
     &            -Work(kOffZ+iTri(J_inBlock,K_inBlock))**2
               End Do
            End Do
            Do jBlock=kBlock+1,nBlock(iSym)
               Do K_inBlock=1,nV(kBlock,iSym)
                  kOffZ=ip_Z(iTri(jBlock,kBlock),iSym)-1
     &                 +nV(jBlock,iSym)*(K_inBlock-1)
                  Do J_inBlock=1,nV(jBlock,iSym)
                     J=iV1(jBlock,iSym)+J_inBlock-1
                     iD=InfVcT(J,1,iSym)
                     Work(kOffD+iD)=Work(kOffD+iD)
     &               -Work(kOffZ+J_inBlock)**2
                  End Do
               End Do
            End Do
         End Do
      End Do

      ! Total number of vectors
      ! ...should equal number of converged diagonals
      nTot=NVT(1)
      Do iSym=2,nSym
         nTot=nTot+NVT(iSym)
      End Do

      ! Count diagonals smaller than threshold
      n1=0
      n2=0
      n3=0
      n4=0
      n5=0
      Dmax=-9.0d9
      Damax=0.0d0
      Dmin=9.0d9
      Damin=9.0d9
      Do iSym=1,nSym
         Do J=1,NVT(iSym)
            iD=InfVcT(J,1,iSym)-1
            Dmax=max(Dmax,Work(ip_D+iD))
            Damax=max(Damax,abs(Work(ip_D+iD)))
            Dmin=min(Dmin,Work(ip_D+iD))
            Damin=min(Damin,abs(Work(ip_D+iD)))
            If (Work(ip_D+iD).le.ThrCom) n1=n1+1
            If (Work(ip_D+iD).lt.0.0d0)  n2=n2+1
            If (Work(ip_D+iD).lt.ThrNeg) n3=n3+1
            If (Work(ip_D+iD).lt.WarNeg) n4=n4+1
            If (Work(ip_D+iD).lt.TooNeg) n5=n5+1
         End Do
      End Do

      ! Write a report if requested
      If (Report) Then
         Call Cho_Head(SecNam//': Report on (J|J) Diagonal from Z','=',
     &                 80,LuPri)
         Write(LuPri,'(/,A,I8)')
     &   'Total dimension of diagonal............',nnBstRT(1)
         Write(LuPri,'(A,I8)')
     &   'Number of Cholesky vectors.............',nTot
         Write(LuPri,'(A,I8)')
     &   'Converged diagonals....................',n1
         Write(LuPri,'(A,I8)')
     &   'Unconverged diagonals..................',nTot-n1
         Write(LuPri,'(A,I8)')
     &   'Negative diagonals.....................',n2
         Write(LuPri,'(A,I8)')
     &   'Neg. diag. that would be zeroed........',n3
         Write(LuPri,'(A,I8)')
     &   'Neg. diag. that would cause warning....',n4
         Write(LuPri,'(A,I8)')
     &   'Neg. diag. that would cause crash......',n5
         Write(LuPri,'(A,1P,D15.6)')
     &   'Max diagonal...........................',Dmax
         Write(LuPri,'(A,1P,D15.6)')
     &   'Min diagonal...........................',Dmin
         Write(LuPri,'(A,1P,D15.6)')
     &   'Max abs diagonal.......................',Damax
         Write(LuPri,'(A,1P,D15.6)')
     &   'Min abs diagonal.......................',Damin
      End If

      ! Deallocation
      Call GetMem('IntDia','Free','Real',ip_D,l_D)

      ! Set return code and return
      If (n1.eq.nTot) Then
         If (n5.ne.0) Then
            irc=-10
         Else
            irc=0
         End If
      Else
         irc=10
      End If

      End
