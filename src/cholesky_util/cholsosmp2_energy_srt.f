!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************
      Subroutine ChoLSOSMP2_Energy_Srt(N,w,t,EOcc,EVir,Delete,EMP2,irc)
!
!     Thomas Bondo Pedersen, December 2012.
!
!     Compute Laplace-SOS-MP2 energy correction from sorted Cholesky
!     vectors (i.e., occupied orbitals processed in batches).
!
      use ChoMP2, only: iFirstS, LnOcc, LnT1am, LiT1am, lUnit
      use Constants
      Implicit None
      Integer N
      Real*8  w(N)
      Real*8  t(N)
      Real*8  EOcc(*)
      Real*8  EVir(*)
      Logical Delete
      Real*8  EMP2
      Integer irc

#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "stdalloc.fh"

      Character(LEN=21), Parameter:: SecNam='ChoLSOSMP2_Energy_Srt'

      Real*8, External:: dDot_

#if !defined (_I8_) || defined (_DEBUGPRINT_)
      Character*2 Unt
      Real*8  Byte
#endif

      Integer nEnrVec(8)

      Integer iSym
      Integer l_X
      Integer Nai
      Integer nBlock
      Integer iOpt, iAddr, l_Tot
      Integer q
      Integer iBlock, jBlock
      Integer iSymi, iSyma
      Integer i, a, ii
      Integer ipi, ipj
      Integer nVeci, nVecj
      Integer iVec
      Integer ip0, ip1
      Integer ipX, lenX
      Integer iBatch
      Integer bsize
      Integer blast

      Real*8  lX, xM, xn, xb, xbp
      Real*8  tq, Eq, wq

      Real*8, Allocatable:: X(:), V(:)
      Integer j, k
      Real*8  epsi, epsa
      Integer MulD2h

      epsi(j,k)=EOcc(iOcc(k)+j)
      epsa(j,k)=EVir(iVir(k)+j)
      MulD2h(j,k)=iEOr(j-1,k-1)+1

      ! init return code
      irc=0

      ! init energy
      EMP2=Zero

      ! check input (incl. common block variables)
      If (nBatch.lt.2) Then
         irc=-1
         Return
      End If
      If (N.ne.Laplace_nGridPoints) Then
         irc=-2
         Return
      End If
      If (Laplace_BlockSize.lt.1) Then
         irc=-3
         Return
      End If

      ! set number of vectors
      If (DecoMP2) Then
         Call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
      Else
         Call iCopy(nSym,NumCho,1,nEnrVec,1)
      End If

      ! allocate X
      lX=Zero
      Do iSym=1,nSym
         If (nT1am(iSym).gt.0 .and. nEnrVec(iSym).gt.0) Then
            bsize=min(Laplace_BlockSize,nEnrVec(iSym))
            nBlock=(nEnrVec(iSym)-1)/bsize+1
            blast=nEnrVec(iSym)-bsize*(nBlock-1)
            xM=dble(nEnrVec(iSym))
            xn=dble(nBlock)
            xb=dble(bsize)
            xbp=dble(blast)
            lX=max(lX,0.5d0*(xM*(xM+1.0d0)
     &                      +(xn-1.0d0)*xb*(xb-1.0d0)
     &                      +xbp*(xbp-1.0d0)))
         End If
      End Do
      l_X=int(lX)
#if !defined (_I8_) || defined (_DEBUGPRINT_)
      If (l_X .lt. 0) Then
         Write(Lupri,'(A,A)')
     &   SecNam,': dimension of X matrix is negative!'
         Write(Lupri,'(A,I15)') 'l_X=',l_X
         If (lX .gt. 0.0d0) Then
            Write(LuPri,'(A)') 'This seems to be an integer overflow!'
            Call Cho_RWord2Byte(lX,Byte,Unt)
            Write(LuPri,'(A,1P,D15.6,A,D15.6,1X,A,A)')
     &      'In double precision, lX=',lX,
     &      ' words (',Byte,Unt,')'
         End If
         irc=1
         Return
      End If
#endif
      Call mma_allocate(X,l_X,Label='X')

      ! allocate vector array
      Nai=0
      Do iBatch=1,nBatch
         Do iSym=1,nSym
            Nai=max(Nai,LnT1am(iSym,iBatch)*nEnrVec(iSym))
         End Do
      End Do
      Call mma_allocate(V,Nai,Label='V')

      ! compute energy correction
      Do q=1,N
         ! init energy for this q
         Eq=0.0d0
         ! scale grid point by 1/2
         tq=0.5d0*t(q)
         ! scale weight by 2 (only lower blocks of X computed)
         wq=2.0d0*w(q)
         Do iSym=1,nSym
            If (nEnrVec(iSym).gt.0) Then
               ! init X for this symmetry
               bsize=min(Laplace_BlockSize,nEnrVec(iSym))
               nBlock=(nEnrVec(iSym)-1)/bsize+1
               blast=nEnrVec(iSym)-bsize*(nBlock-1)
               lenX=nEnrVec(iSym)*(nEnrVec(iSym)+1)/2
     &             +(nBlock-1)*(bsize*(bsize-1)/2)
     &             +blast*(blast-1)/2
#if defined (_DEBUGPRINT_)
               If (lenX.gt.l_X) Then
                  Call WarningMessage(2,
     &                            SecNam//': insufficient X allocation')
                  Write(6,'(A,2(1X,I10))') 'lenX,l_X=',lenX,l_X
                  Call Abend()
               End If
#endif
               X(1:lenX)=Zero
               Do iBatch=1,nBatch
                  Nai=LnT1am(iSym,iBatch)
                  If (Nai.gt.0) Then
                     ! open file, read vectors, close file
                     ! do not delete file here - it may be needed later
                     ! (because of the loop over q)
                     Call ChoMP2_OpenB(1,iSym,iBatch)
                     iOpt=2
                     l_Tot=Nai*nEnrVec(iSym)
                     iAddr=1
#if defined (_DEBUGPRINT_)
                     If (l_Tot.gt.SIZE(V)) Then
                        Call WarningMessage(2,
     &                            SecNam//': insufficient V allocation')
                        Write(6,'(A,2(1X,I10))') 'l_Tot,SIZE(V)=',
     &                                            l_Tot,SIZE(V)
                        Call Abend()
                     End If
#endif
                     Call dDAFile(lUnit(iSym,iBatch),iOpt,V,l_Tot,iAddr)
                     Call ChoMP2_OpenB(2,iSym,iBatch)
                     ! scale vectors
                     Do iVec=1,nEnrVec(iSym)
                        ip0=Nai*(iVec-1)
                        Do iSymi=1,nSym
                           If (LnOcc(iSymi,iBatch).gt.0) Then
                              iSyma=MulD2h(iSym,iSymi)
                              ip1=ip0+LiT1am(iSyma,iSymi,iBatch)
                              Do i=0,LnOcc(iSymi,iBatch)-1
                                 ii=iFirstS(iSymi,iBatch)+i
                                 Call dScal_(nVir(iSyma),
     &                                      exp(epsi(ii,iSymi)*tq),
     &                                      V(ip1+nVir(iSyma)*i+1),1)
                              End Do
                              Do a=1,nVir(iSyma)
                                 Call dScal_(LnOcc(iSymi,iBatch),
     &                                      exp(-epsa(a,iSyma)*tq),
     &                                      V(ip1+a),nVir(iSyma))
                              End Do
                           End If
                        End Do
                     End Do
                     ! loop over vector blocks to compute
                     ! X(J,K) +=
                     ! sum_ai L(ai,J)*L(ai,K)*exp(-(e(a)-e(i))*t(q)/2)
                     ipX=1
                     Do jBlock=1,nBlock
                        ipj=1+Nai*Laplace_BlockSize*(jBlock-1)
                        If (jBlock.eq.nBlock) Then
                           nVecj=nEnrVec(iSym)
     &                          -Laplace_BlockSize*(nBlock-1)
                        Else
                           nVecj=Laplace_BlockSize
                        End If
                        Do iBlock=jBlock,nBlock
                           ipi=1+Nai*Laplace_BlockSize*(iBlock-1)
                           If (iBlock.eq.nBlock) Then
                              nVeci=nEnrVec(iSym)
     &                             -Laplace_BlockSize*(nBlock-1)
                           Else
                              nVeci=Laplace_BlockSize
                           End If
                           Call dGEMM_('T','N',nVeci,nVecj,Nai,
     &                                1.0d0,V(ipi),Nai,
     &                                      V(ipj),Nai,
     &                                1.0d0,X(ipX),nVeci)
                           ipX=ipX+nVeci*nVecj
                        End Do
                     End Do
#if defined (_DEBUGPRINT_)
                     If (lenX.ne.(ipX-1)) Then
                        Call WarningMessage(2,
     &                                SecNam//': dimension problem [1]')
                        Write(6,'(A,I10,A,I10)')
     &                  'lenX=',lenX,' ipX-1=',ipX-1
                        Call Abend()
                     End If
#endif
                  End If
               End Do
               ! compute energy contribution
               ! Eq += sum_JK [X(J,K)]**2
               ipX=1
               Do jBlock=1,nBlock
                  If (jBlock.eq.nBlock) Then
                     nVecj=nEnrVec(iSym)
     &                    -Laplace_BlockSize*(nBlock-1)
                  Else
                     nVecj=Laplace_BlockSize
                  End If
                  Do iBlock=jBlock,nBlock
                     If (iBlock.eq.nBlock) Then
                        nVeci=nEnrVec(iSym)
     &                       -Laplace_BlockSize*(nBlock-1)
                     Else
                        nVeci=Laplace_BlockSize
                     End If
                     If (iBlock.eq.jBlock) Then
                        Eq=Eq+0.5d0*dDot_(nVeci*nVecj,X(ipX),1,X(ipX),1)
                     Else
                        Eq=Eq+dDot_(nVeci*nVecj,X(ipX),1,X(ipX),1)
                     End If
                     ipX=ipX+nVeci*nVecj
                  End Do
               End Do
#if defined (_DEBUGPRINT_)
               If (lenX.ne.(ipX-1)) Then
                  Call WarningMessage(2,
     &                             SecNam//': dimension problem [2]')
                  Write(6,'(A,I10,A,I10)')
     &            'lenX=',lenX,' ipX-1=',ipX-1
                  Call Abend()
               End If
#endif
            End If
         End Do
         ! scale Eq
         Eq=-wq*Eq
         ! Accumulate in EMP2
         EMP2=EMP2+Eq
      End Do

      ! deallocations
      Call mma_deallocate(V)
      Call mma_deallocate(X)

      ! delete files if requested
      If (Delete) Then
         Do iBatch = 1,nBatch
            Do iSym = 1,nSym
               Call ChoMP2_OpenB(1,iSym,iBatch)
               Call ChoMP2_OpenB(3,iSym,iBatch)
            End Do
         End Do
      End If

      End
