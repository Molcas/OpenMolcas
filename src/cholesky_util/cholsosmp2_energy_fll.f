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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine ChoLSOSMP2_Energy_Fll(N,w,t,EOcc,EVir,Delete,EMP2,irc)
C
C     Thomas Bondo Pedersen, December 2012.
C
C     Compute Laplace-SOS-MP2 energy correction from full Cholesky
C     vectors (i.e., not batched).
C
      Implicit None
      Integer N
      Real*8  w(N)
      Real*8  t(N)
      Real*8  EOcc(*)
      Real*8  EVir(*)
      Logical Delete
      Real*8  EMP2
      Integer irc
#include "WrkSpc.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "cholesky.fh"

      Character*21 SecNam
      Parameter (SecNam='ChoLSOSMP2_Energy_Fll')

      Integer nEnrVec(8)

      Integer l_X
      Integer l_V
      Integer iSym
      Integer Nai
      Integer need
      Integer ip, l

      ! check if there is enough memory to read through vector files
      ! only once
      If (DecoMP2) Then
         Call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
      Else
         Call iCopy(nSym,NumCho,1,nEnrVec,1)
      End If
      l_X=0
      l_V=0
      Do iSym=1,nSym
         Nai=nT1am(iSym)
         If (Nai.gt.0 .and. nEnrVec(iSym).gt.0) Then
            l_X=max(l_X,min(Laplace_BlockSize,nEnrVec(iSym)))
            l_V=max(l_V,Nai*nEnrVec(iSym))
         End If
      End Do
      need=l_X+2*l_V
      Call GetMem('LSMTst','Max ','Real',ip,l)
      l=l-100
      If (l.lt.1 .or. need.ge.l) Then ! not enough memory for one read
         Call ChoLSOSMP2_Energy_Fll2(N,w,t,EOcc,EVir,Delete,EMP2,irc)
         If (irc.ne.0) Then
            Write(6,'(A,A,I10)')
     &      SecNam,': Cho_LSOSMP2_Energy_Fll2 returned',irc
         End If
      Else ! enough memory for one read through vector files
         Call ChoLSOSMP2_Energy_Fll1(N,w,t,EOcc,EVir,Delete,EMP2,irc)
         If (irc.ne.0) Then
            Write(6,'(A,A,I10)')
     &      SecNam,': Cho_LSOSMP2_Energy_Fll1 returned',irc
         End If
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine ChoLSOSMP2_Energy_Fll1(N,w,t,EOcc,EVir,Delete,EMP2,irc)
C
C     Thomas Bondo Pedersen, December 2012.
C
C     Compute Laplace-SOS-MP2 energy correction from full Cholesky
C     vectors (i.e., not batched), reading the vectors only once
C     at the expense of memory.
C
      Implicit None
      Integer N
      Real*8  w(N)
      Real*8  t(N)
      Real*8  EOcc(*)
      Real*8  EVir(*)
      Logical Delete
      Real*8  EMP2
      Integer irc
#include "WrkSpc.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "cholesky.fh"

      Character*22 SecNam
      Parameter (SecNam='ChoLSOSMP2_Energy_Fll1')

      Real*8   dDot_
      external ddot_

      Integer nEnrVec(8)

      Integer iClos
      Integer iTyp
      Integer iSym
      Integer ip_X, l_X
      Integer ip_V, l_V, ip_SV
      Integer Nai
      Integer nBlock
      Integer iOpt, iAddr, l_Tot
      Integer q
      Integer iBlock, jBlock
      Integer iSymi, iSyma
      Integer i, a
      Integer ipi, ipj
      Integer nVeci, nVecj
      Integer iVec
      Integer ip0, ip1

      Real*8  tq, Eq

      Integer j, k
      Real*8  epsi, epsa
      Integer MulD2h
      epsi(j,k)=EOcc(iOcc(k)+j)
      epsa(j,k)=EVir(iVir(k)+j)
      MulD2h(j,k)=iEOr(j-1,k-1)+1

      ! init return code
      irc=0

      ! init energy
      EMP2=0.0d0

      ! check input (incl. common block variables)
      If (nBatch.ne.1) Then
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

      ! determine if files are to be deleted after use
      If (Delete) Then
         iClos=3
      Else
         iClos=2
      End If

      ! set number and type of vectors
      If (DecoMP2) Then
         iTyp=2
         Call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
      Else
         iTyp=1
         Call iCopy(nSym,NumCho,1,nEnrVec,1)
      End If

      ! compute energy correction
      Do iSym=1,nSym
         Nai=nT1am(iSym)
         If (Nai.gt.0 .and. nEnrVec(iSym).gt.0) Then
            ! set vector block size and allocate X matrix
            l_X=min(Laplace_BlockSize,nEnrVec(iSym))
            l_X=l_X**2
            Call GetMem('LSMX','Allo','Real',ip_X,l_X)
            ! compute number of vector blocks
            nBlock=(nEnrVec(iSym)-1)/Laplace_BlockSize+1
            ! open vector file
            Call ChoMP2_OpenF(1,iTyp,iSym)
            ! allocate memory for vectors
            l_Tot=Nai*nEnrVec(iSym)
            l_V=2*l_Tot
            Call GetMem('LSMV','Allo','Real',ip_V,l_V)
            ip_SV=ip_V+l_Tot
            ! Read all vectors
            iOpt=2
            iAddr=1
            Call dDAFile(lUnit_F(iSym,iTyp),iOpt,Work(ip_V),l_Tot,iAddr)
            ! loop over Laplace grid
            Do q=1,N
               ! init energy for this q
               Eq=0.0d0
               ! scale grid point by 1/2
               tq=0.5d0*t(q)
               ! scale vectors
               Call dCopy_(l_Tot,Work(ip_V),1,Work(ip_SV),1)
               Do iVec=1,nEnrVec(iSym)
                  ip0=ip_SV-1+Nai*(iVec-1)
                  Do iSymi=1,nSym
                     iSyma=MulD2h(iSym,iSymi)
                     ip1=ip0+iT1am(iSyma,iSymi)
                     Do i=1,nOcc(iSymi)
                        Call dScal_(nVir(iSyma),exp(epsi(i,iSymi)*tq),
     &                             Work(ip1+nVir(iSyma)*(i-1)+1),1)
                     End Do
                     Do a=1,nVir(iSyma)
                        Call dScal_(nOcc(iSymi),exp(-epsa(a,iSyma)*tq),
     &                             Work(ip1+a),nVir(iSyma))
                     End Do
                  End Do
               End Do
               ! loop over vector blocks to compute
               ! X(J,K) = sum_ai L(ai,J)*L(ai,K)*exp(-(e(a)-a(i))*t(q)/2)
               ! Eq += w(q)*sum_JK [X(J,K)]**2
               Do jBlock=1,nBlock
                  ipj=ip_SV+Nai*Laplace_BlockSize*(jBlock-1)
                  If (jBlock.eq.nBlock) Then
                     nVecj=nEnrVec(iSym)-Laplace_BlockSize*(nBlock-1)
                  Else
                     nVecj=Laplace_BlockSize
                  End If
                  Do iBlock=jBlock,nBlock
                     ipi=ip_SV+Nai*Laplace_BlockSize*(iBlock-1)
                     If (iBlock.eq.nBlock) Then
                        nVeci=nEnrVec(iSym)-Laplace_BlockSize*(nBlock-1)
                     Else
                        nVeci=Laplace_BlockSize
                     End If
                     Call dGEMM_('T','N',nVeci,nVecj,Nai,
     &                           1.0d0,Work(ipi),Nai,Work(ipj),Nai,
     &                           0.0d0,Work(ip_X),nVeci)
                     If (iBlock.eq.jBlock) Then
                        Eq=Eq
     &               +0.5d0*dDot_(nVeci*nVecj,Work(ip_X),1,Work(ip_X),1)
                     Else
                        Eq=Eq
     &                     +dDot_(nVeci*nVecj,Work(ip_X),1,Work(ip_X),1)
                     End If
                  End Do
               End Do
               ! accumulate in EMP2
               EMP2=EMP2-w(q)*Eq
            End Do
            ! deallocate memory
            Call GetMem('LSMV','Free','Real',ip_V,l_V)
            Call GetMem('LSMX','Free','Real',ip_X,l_X)
            ! close (and possibly delete) vector file
            Call ChoMP2_OpenF(iClos,iTyp,iSym)
         End If
      End Do

      ! Scale energy (only half of X computed)
      EMP2=2.0d0*EMP2

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine ChoLSOSMP2_Energy_Fll2(N,w,t,EOcc,EVir,Delete,EMP2,irc)
C
C     Thomas Bondo Pedersen, December 2012.
C
C     Compute Laplace-SOS-MP2 energy correction from unsorted Cholesky
C     vectors (i.e., no batching). This uses the exact same loop
C     ordering as the sorted algorithm and thus reads through the
C     vector files once for each grid point.
C
      Implicit None
      Integer N
      Real*8  w(N)
      Real*8  t(N)
      Real*8  EOcc(*)
      Real*8  EVir(*)
      Logical Delete
      Real*8  EMP2
      Integer irc
#include "WrkSpc.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "cholesky.fh"

      Character*22 SecNam
      Parameter (SecNam='ChoLSOSMP2_Energy_Fll2')

      Real*8   dDot_
      external ddot_

#if !defined (_I8_) || defined (_DEBUG_)
      Character*2 Unt
      Real*8  Byte
#endif

      Integer nEnrVec(8)

      Integer iTyp
      Integer iSym
      Integer ip_X, l_X
      Integer ip_V, l_V
      Integer Nai
      Integer nBlock
      Integer iOpt, iAddr, l_Tot
      Integer q
      Integer iBlock, jBlock
      Integer iSymi, iSyma
      Integer i, a
      Integer ipi, ipj
      Integer nVeci, nVecj
      Integer iVec
      Integer ip0, ip1
      Integer ipX, lenX
      Integer bsize
      Integer blast

      Real*8  lX, xM, xn, xb, xbp
      Real*8  tq, Eq, wq

      Integer j, k
      Real*8  epsi, epsa
      Integer MulD2h
      epsi(j,k)=EOcc(iOcc(k)+j)
      epsa(j,k)=EVir(iVir(k)+j)
      MulD2h(j,k)=iEOr(j-1,k-1)+1

      ! init return code
      irc=0

      ! init energy
      EMP2=0.0d0

      ! check input (incl. common block variables)
      If (N.ne.Laplace_nGridPoints) Then
         irc=-2
         Return
      End If
      If (Laplace_BlockSize.lt.1) Then
         irc=-3
         Return
      End If

      ! set number and type of vectors
      If (DecoMP2) Then
         iTyp=2
         Call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
      Else
         iTyp=1
         Call iCopy(nSym,NumCho,1,nEnrVec,1)
      End If

      ! allocate X
      lX=0.0d0
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
#if !defined (_I8_) || defined (_DEBUG_)
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
      Call GetMem('LSMX','Allo','Real',ip_X,l_X)

      ! allocate vector array
      Nai=nT1am(1)*nEnrVec(1)
      Do iSym=2,nSym
         Nai=max(Nai,nT1am(iSym)*nEnrVec(iSym))
      End Do
      l_V=Nai
      Call GetMem('LSMV','Allo','Real',ip_V,l_V)

      ! compute energy correction
      Do q=1,N
         ! init energy for this q
         Eq=0.0d0
         ! scale grid point by 1/2
         tq=0.5d0*t(q)
         ! scale weight by 2 (only lower blocks of X computed)
         wq=2.0d0*w(q)
         Do iSym=1,nSym
            Nai=nT1am(iSym)
            If (Nai.gt.0 .and. nEnrVec(iSym).gt.0) Then
               ! init X for this symmetry
               bsize=min(Laplace_BlockSize,nEnrVec(iSym))
               nBlock=(nEnrVec(iSym)-1)/bsize+1
               blast=nEnrVec(iSym)-bsize*(nBlock-1)
               lenX=nEnrVec(iSym)*(nEnrVec(iSym)+1)/2
     &             +(nBlock-1)*(bsize*(bsize-1)/2)
     &             +blast*(blast-1)/2
               Call fZero(Work(ip_X),lenX)
               ! open file, read vectors, close file
               ! do not delete file here - it may be needed later
               ! (because of the loop over q)
               Call ChoMP2_OpenF(1,iTyp,iSym)
               iOpt=2
               l_Tot=Nai*nEnrVec(iSym)
               iAddr=1
               Call dDAFile(lUnit_F(iSym,iTyp),iOpt,Work(ip_V),
     &                      l_Tot,iAddr)
               Call ChoMP2_OpenF(2,iTyp,iSym)
               ! scale vectors
               Do iVec=1,nEnrVec(iSym)
                  ip0=ip_V-1+Nai*(iVec-1)
                  Do iSymi=1,nSym
                     If (nOcc(iSymi).gt.0) Then
                        iSyma=MulD2h(iSym,iSymi)
                        ip1=ip0+iT1am(iSyma,iSymi)
                        Do i=1,nOcc(iSymi)
                           Call dScal_(nVir(iSyma),
     &                                exp(epsi(i,iSymi)*tq),
     &                                Work(ip1+nVir(iSyma)*(i-1)+1),1)
                        End Do
                        Do a=1,nVir(iSyma)
                           Call dScal_(nOcc(iSymi),
     &                                exp(-epsa(a,iSyma)*tq),
     &                                Work(ip1+a),nVir(iSyma))
                        End Do
                     End If
                  End Do
               End Do
               ! loop over vector blocks to compute
               ! X(J,K) +=
               ! sum_ai L(ai,J)*L(ai,K)*exp(-(e(a)-e(i))*t(q)/2)
               ipX=ip_X
               Do jBlock=1,nBlock
                  ipj=ip_V+Nai*Laplace_BlockSize*(jBlock-1)
                  If (jBlock.eq.nBlock) Then
                     nVecj=nEnrVec(iSym)
     &                    -Laplace_BlockSize*(nBlock-1)
                  Else
                     nVecj=Laplace_BlockSize
                  End If
                  Do iBlock=jBlock,nBlock
                     ipi=ip_V+Nai*Laplace_BlockSize*(iBlock-1)
                     If (iBlock.eq.nBlock) Then
                        nVeci=nEnrVec(iSym)
     &                       -Laplace_BlockSize*(nBlock-1)
                     Else
                        nVeci=Laplace_BlockSize
                     End If
                     Call dGEMM_('T','N',nVeci,nVecj,Nai,
     &                          1.0d0,Work(ipi),Nai,Work(ipj),Nai,
     &                          1.0d0,Work(ipX),nVeci)
                     ipX=ipX+nVeci*nVecj
                  End Do
               End Do
#if defined (_DEBUG_)
               If (lenX.ne.(ipX-ip_X)) Then
                  Call WarningMessage(2,
     &                          SecNam//': dimension problem [1]')
                  Write(6,'(A,I10,A,I10)')
     &            'lenX=',lenX,' ipX-ip_X=',ipX-ip_X
                  Call Abend()
               End If
#endif
               ! compute energy contribution
               ! Eq += sum_JK [X(J,K)]**2
               ipX=ip_X
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
                        Eq=Eq+0.5d0*dDot_(nVeci*nVecj,
     &                                   Work(ipX),1,Work(ipX),1)
                     Else
                        Eq=Eq+dDot_(nVeci*nVecj,
     &                             Work(ipX),1,Work(ipX),1)
                     End If
                     ipX=ipX+nVeci*nVecj
                  End Do
               End Do
#if defined (_DEBUG_)
               If (lenX.ne.(ipX-ip_X)) Then
                  Call WarningMessage(2,
     &                             SecNam//': dimension problem [2]')
                  Write(6,'(A,I10,A,I10)')
     &            'lenX=',lenX,' ipX-ip_X=',ipX-ip_X
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
      Call GetMem('LSMV','Free','Real',ip_V,l_V)
      Call GetMem('LSMX','Free','Real',ip_X,l_X)

      ! delete files if requested
      If (Delete) Then
         Do iSym=1,nSym
            Call ChoMP2_OpenF(1,iTyp,iSym)
            Call ChoMP2_OpenF(3,iTyp,iSym)
         End Do
      End If

      End
