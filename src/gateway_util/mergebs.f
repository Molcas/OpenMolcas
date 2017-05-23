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
      SubRoutine MergeBS (z1,n1,z2,n2,z,n,RatioThres,iDominantSet)
      Implicit Real*8 (a-h,o-z)
      Real*8 z1(n1), z2(n2), z(*)
      Parameter (mPrim=60)
      Integer ix1(mPrim), ix2(mPrim)
      Logical IfTest
      Data IfTest/.False./
*
#ifdef _DEBUG_
      Call QEnter('MergeBS')
      IfTest=.True.
#endif
      If (n1.gt.mPrim .or. n2.gt.mPrim) Then
        Call WarningMessage(2,'Error in MergeBS')
        Write (6,*) ' MergeBS: n1,n2 .gt.mPrim',n1,n2,mPrim
        Write (6,*) ' MergeBS: rise mPrim and recompile'
        Call Abend()
      End If
*
      iSetPrev     = 0
*
      Do i = 1, mPrim
        ix1(i) = i
        ix2(i) = i
      End Do
*
      Do i = 1, n1 - 1
        Do j = i + 1, n1
          If (z1(ix1(i)).lt.z1(ix1(j)))  then
            idum=ix1(i)
            ix1(i)=ix1(j)
            ix1(j)=idum
          End If
        End Do
      End Do
*
      Do i = 1, n2 - 1
        Do j = i + 1, n2
          If (z2(ix2(i)).lt.z2(ix2(j)))  then
            idum=ix2(i)
            ix2(i)=ix2(j)
            ix2(j)=idum
          End If
        End Do
      End Do
*
      If (IfTest) Write (6,'(A)')
      If (IfTest) Write (6,'(4f20.4)') (z1(ix1(i)),i=1,n1)
      If (IfTest) Write (6,'(A)')
      If (IfTest) Write (6,'(4f20.4)') (z2(ix2(i)),i=1,n2)
*
      i = 0
      i1 = 1
      i2 = 1
      Do while (i1.le.n1 .or. i2.le.n2)
         i = i + 1
        If (i.gt.mPrim) Then
           Call WarningMessage(2,'Error in MergeBS')
           Write (6,*) ' MergeBS: i.gt.mPrim',i,mPrim
           Write (6,*) ' MergeBS: rise mPrim and recompile'
           Call Abend()
        End If
        If ( i1.gt.n1 ) Then
          z(i) = z2(ix2(i2))
          iSetThis = 2
          i2 = i2 + 1
        Else If ( i2.gt.n2 ) Then
          z(i) = z1(ix1(i1))
          iSetThis = 1
          i1 = i1 + 1
        Else If ( z1(ix1(i1)) .gt. z2(ix2(i2)) ) then
          z(i) = z1(ix1(i1))
          iSetThis = 1
          i1 = i1 + 1
        Else
          z(i) = z2(ix2(i2))
          iSetThis = 2
          i2 = i2 + 1
        End If
*
        If (i.eq.1) Then
          iSetPrev = iSetThis
        Else
          ratio = z(i-1)/z(i)
          If (ratio.ge.RatioThres) Then
            iSetPrev = iSetThis
          Else
            If (iSetThis.ne.iDominantSet) Then
              i = i - 1
            Else If (iSetThis.eq.iDominantSet) Then
              If (iSetPrev.ne.iDominantSet) Then
                z(i-1) = z(i)
                i = i - 1
              Else
                continue
              End If
              iSetPrev = iSetThis
            End If
          End If
        End If

      End Do
*
      n = i
      If (IfTest) Write (6,'(I4)') n
      If (IfTest) Write (6,'(4f20.4)') (z(i),i=1,n)
#ifdef _DEBUG_
      Call QExit('MergeBS')
#endif
*
      Return
      End
