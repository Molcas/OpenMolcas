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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************
#ifdef _OLD_CODE_
      Subroutine NTMul(A,B,C,nRowA,nColA,nRowB)
      use Constants, only: Zero
      Implicit None
      Integer nRowA,nColA,nRowB
      Real*8 A(nRowA,nColA), B(nRowB,nColA),                            &
     &       C(nRowA,nRowB)
!
      Integer nCache, mCache, Indj, jj, njVec, i, j, k
      nCache=(64/8)*1024
      mCache=(nCache*3)/4 - nRowA*nColA
      Incj=mCache/(nRowA+nColA)
!
!-----Sectioning of long index
!
      Do jj = 1, nRowB, Incj
         njVec=Min(Incj,nRowB-jj+1)
!
         Do i = 1, nRowA
!-----------Set target to zero
            Do j = jj, jj+njVec-1
               C(i,j) = Zero
            End Do
            Do k = 1, nColA
               If (A(i,k).ne.Zero) Then
                  Do j = jj, jj+njVec-1
                     C(i,j) = C(i,j) + A(i,k)*B(j,k)
                  End Do
               End If
            End Do
         End Do
!
      End Do    ! End of sectioning
!
      Return
      End Subroutine NTMul
#else
!--------------------------------------------------------------------
      subroutine ntmul(a,b,r,ncol,nlink,nrow)
!--------------------------------------------------------------------
      Implicit None
      Integer nCol, nRow, nLink
      Real*8 r(ncol,*),a(ncol,*),b(nrow,*)

      Integer, parameter :: mxind=2000
      Integer ind(mxind)
      Integer i, nnot, k, j, nr1
      Real*8 S1, S2, S3, S4, S5, S6, S7, S8,                            &
     &       T1, T2, T3, T4, T5, T6, T7, T8
!
      do 100 i=1,ncol
!
      nnot=0
      do k=1,min(nlink,mxind)
        if ( a(i,k) .ne. 0.0d0 ) then
          nnot = nnot + 1
          ind(nnot) = k
        end if
      end do
!
      do 20 j=1,nrow-15,16
      s1=0.0D0
      s2=0.0D0
      s3=0.0D0
      s4=0.0D0
      s5=0.0D0
      s6=0.0D0
      s7=0.0D0
      s8=0.0D0
      t1=0.0D0
      t2=0.0D0
      t3=0.0D0
      t4=0.0D0
      t5=0.0D0
      t6=0.0D0
      t7=0.0D0
      t8=0.0D0
      do 15 k=1,nnot
      s1=s1+a(i,ind(k))*b(j,  ind(k))
      s2=s2+a(i,ind(k))*b(j+1,ind(k))
      s3=s3+a(i,ind(k))*b(j+2,ind(k))
      s4=s4+a(i,ind(k))*b(j+3,ind(k))
      s5=s5+a(i,ind(k))*b(j+4,ind(k))
      s6=s6+a(i,ind(k))*b(j+5,ind(k))
      s7=s7+a(i,ind(k))*b(j+6,ind(k))
      s8=s8+a(i,ind(k))*b(j+7,ind(k))
      t1=t1+a(i,ind(k))*b(j+8,ind(k))
      t2=t2+a(i,ind(k))*b(j+9,ind(k))
      t3=t3+a(i,ind(k))*b(j+10,ind(k))
      t4=t4+a(i,ind(k))*b(j+11,ind(k))
      t5=t5+a(i,ind(k))*b(j+12,ind(k))
      t6=t6+a(i,ind(k))*b(j+13,ind(k))
      t7=t7+a(i,ind(k))*b(j+14,ind(k))
      t8=t8+a(i,ind(k))*b(j+15,ind(k))
15    continue
      r(i,j  )=s1
      r(i,j+1)=s2
      r(i,j+2)=s3
      r(i,j+3)=s4
      r(i,j+4)=s5
      r(i,j+5)=s6
      r(i,j+6)=s7
      r(i,j+7)=s8
      r(i,j+8)=t1
      r(i,j+9)=t2
      r(i,j+10)=t3
      r(i,j+11)=t4
      r(i,j+12)=t5
      r(i,j+13)=t6
      r(i,j+14)=t7
      r(i,j+15)=t8
20    continue
!
      nr1=mod(nrow,16)
      if(nr1.eq.0) goto 100
      j=nrow-nr1+1
!
      if(nr1.ge.8) then
      s1=0.0D0
      s2=0.0D0
      s3=0.0D0
      s4=0.0D0
      s5=0.0D0
      s6=0.0D0
      s7=0.0D0
      s8=0.0D0
      do 25 k=1,nnot
      s1=s1+a(i,ind(k))*b(j,  ind(k))
      s2=s2+a(i,ind(k))*b(j+1,ind(k))
      s3=s3+a(i,ind(k))*b(j+2,ind(k))
      s4=s4+a(i,ind(k))*b(j+3,ind(k))
      s5=s5+a(i,ind(k))*b(j+4,ind(k))
      s6=s6+a(i,ind(k))*b(j+5,ind(k))
      s7=s7+a(i,ind(k))*b(j+6,ind(k))
      s8=s8+a(i,ind(k))*b(j+7,ind(k))
25    continue
      r(i,j  )=s1
      r(i,j+1)=s2
      r(i,j+2)=s3
      r(i,j+3)=s4
      r(i,j+4)=s5
      r(i,j+5)=s6
      r(i,j+6)=s7
      r(i,j+7)=s8
      nr1=nr1-8
      j=j+8
      end if
!
      if(nr1.ge.4) then
      s1=0.0D0
      s2=0.0D0
      s3=0.0D0
      s4=0.0D0
      do 35 k=1,nnot
      s1=s1+a(i,ind(k))*b(j,  ind(k))
      s2=s2+a(i,ind(k))*b(j+1,ind(k))
      s3=s3+a(i,ind(k))*b(j+2,ind(k))
      s4=s4+a(i,ind(k))*b(j+3,ind(k))
35    continue
      r(i,j  )=s1
      r(i,j+1)=s2
      r(i,j+2)=s3
      r(i,j+3)=s4
      nr1=nr1-4
      j=j+4
      end if
!
      If (nr1.eq.1) Then
         s1=0.0D0
         do 45 k=1,nnot
         s1=s1+a(i,ind(k))*b(j,ind(k))
45       continue
         r(i,j)=s1
      Else If (nr1.eq.2) Then
         s1=0.0D0
         s2=0.0D0
         do 50 k=1,nnot
         s1=s1+a(i,ind(k))*b(j,  ind(k))
         s2=s2+a(i,ind(k))*b(j+1,ind(k))
50       continue
         r(i,j  )=s1
         r(i,j+1)=s2
      Else If (nr1.eq.3) Then
         s1=0.0D0
         s2=0.0D0
         s3=0.0D0
         do 55 k=1,nnot
         s1=s1+a(i,ind(k))*b(j,  ind(k))
         s2=s2+a(i,ind(k))*b(j+1,ind(k))
         s3=s3+a(i,ind(k))*b(j+2,ind(k))
55       continue
         r(i,j  )=s1
         r(i,j+1)=s2
         r(i,j+2)=s3
      Else If (nr1.gt.3) Then
         Call WarningMessage(2,'nr1.gt.3')
         Call Abend()
      End If
 100  Continue
!
      if(mxind.ge.nlink) return
      Call WarningMessage(2,'MxInd.lt.nLink')
      Write (6,*) 'mxind,nlink=',mxind,nlink
      Call Abend()
      Return
      End subroutine ntmul
#endif
