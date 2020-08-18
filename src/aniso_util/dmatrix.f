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
      Subroutine DMATRIX(E,F,Z,IPRINT)

C   THIS ROUTINE COMPUTES THE USUAL D-TENSOR ON THA BASIS OF COEFFICINETS
C   OF THE STEWENS OPERATORS OF ORDER 2 (ES AND FS) AND DIAGONALIZE IT
C   TO OBTAIN THE MAIN ANISOTROPY AXES
C

      Implicit None
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)        :: iprint
      Real(kind=8),intent(in)   :: Z(3,3)
      Complex(kind=8),intent(in):: E(0:2), F(0:2)
      Integer                    :: I, INFO, J
      Real(kind=8)              :: WD(3), CF,Unity(3,3), SMAT(3,3),
     &                              DMATR(3,3), ZD(3,3),
     &                              dtens(3), daxes(3,3), D_factor,
     &                              E_factor, diff12, diff23, ZD2(3,3)
      Complex(kind=8)           :: DMAT(3,3)

      Call qEnter('dmatrix')

      CF=sqrt(3.0_wp/2.0_wp)

      DMAT(1,1) =  CF*E(2)-E(0)
      DMAT(2,2) = -CF*E(2)-E(0)
      DMAT(3,3) =  2.0_wp*E(0)
      DMAT(1,2) =  CF*F(2)
      DMAT(1,3) =  CF*E(1)
      DMAT(2,3) =  CF*F(1)
      DMAT(2,1) =  DMAT(1,2)
      DMAT(3,1) =  DMAT(1,3)
      DMAT(3,2) =  DMAT(2,3)

      Do i=1,3
        Do j=1,3
          DMATR(i,j)=DBLE(DMAT(i,j))
        End Do
      End Do

      Call DIAG_R2(DMATR,3,INFO,WD,ZD)

c  calculate the rotation matrix:
      Call dcopy_(3*3,[0.0_wp],0,Unity,1)
      Do i=1,3
        Unity(i,i)=1.0_wp
      End Do

      Call DGEMM_('N','N',3,3,3,1.0_wp,ZD,3,Unity,3,0.0_wp,SMAT,3)

C Set the dtens and daxes with respect to the gtens and maxes.
cccccccccccccccccccccccccccccccc
      dtens(:)=0.0_wp
      If ( (abs(SMAT(1,1)).gt.abs(SMAT(1,2))) .AND.
     &     (abs(SMAT(1,1)).gt.abs(SMAT(1,3))) ) Then
c  the WD(1) and ZD(i,1) correspond to gtens(1) and maxes(i,1)
        dtens(1)=WD(1)
        If (SMAT(1,1).gt.0.0_wp) Then
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(1,1) is larger than SMAT(1,2) and ' //
     &                     'SMAT(1,3) and is positive'
          End If
          Do i=1,3
            daxes(i,1)=ZD(i,1)
          End Do
        Else
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(1,1) is larger than SMAT(1,2) and ' //
     &                     'SMAT(1,3) and is negative'
          End If
          Do i=1,3
            daxes(i,1)=-ZD(i,1)
          End Do
        End If

      Else If ( (abs(SMAT(1,2)).gt.abs(SMAT(1,1))) .AND.
     &          (abs(SMAT(1,2)).gt.abs(SMAT(1,3))) ) Then
c the WD(1) and ZD(i,1) correspond to gtens(2) and maxes(i,2)
        dtens(1)=WD(2)
        If (SMAT(1,2).gt.0.0_wp) Then
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(1,2) is larger than SMAT(1,1) and ' //
     &                     'SMAT(1,3) and is positive'
          End If
          Do i=1,3
            daxes(i,1)=ZD(i,2)
          End Do
        Else
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(1,2) is larger than SMAT(1,1) and ' //
     &                     'SMAT(1,3) and is negative'
          End If
          Do i=1,3
            daxes(i,1)=-ZD(i,2)
          End Do
        End If

      Else If ( (abs(SMAT(1,3)).gt.abs(SMAT(1,1))) .AND.
     &          (abs(SMAT(1,3)).gt.abs(SMAT(1,2))) ) Then
c the WD(1) and ZD(i,1) correspond to gtens(3) and maxes(i,3)
        dtens(1)=WD(3)
        If (SMAT(1,3).gt.0.0_wp) Then
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(1,3) is larger than SMAT(1,1) and ' //
     &                     'SMAT(1,2) and is positive'
          End If
          Do i=1,3
            daxes(i,1)=ZD(i,3)
          End Do
        Else
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(1,3) is larger than SMAT(1,1) and ' //
     &                     'SMAT(1,2) and is negative'
          End If
          Do i=1,3
            daxes(i,1)=-ZD(i,3)
          End Do
        End If
      End If
cccccccccccccccccccccccccccccccc

      If ( (abs(SMAT(2,1)).gt.abs(SMAT(2,2))) .AND.
     &     (abs(SMAT(2,1)).gt.abs(SMAT(2,3))) ) Then
c  the WD(2) and ZD(i,2) correspond to gtens(1) and maxes(i,1)
        dtens(2)=WD(1)
        If (SMAT(2,1).gt.0.0_wp) Then
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(2,1) is larger than SMAT(2,2) and ' //
     &                     'SMAT(2,3) and is positive'
          End If
          Do i=1,3
            daxes(i,2)=ZD(i,1)
          End Do
        Else
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(2,1) is larger than SMAT(2,2) and ' //
     &                     'SMAT(2,3) and is negative'
          End If
          Do i=1,3
            daxes(i,2)=-ZD(i,1)
          End Do
        End If

      Else If ( (abs(SMAT(2,2)).gt.abs(SMAT(2,1))) .AND.
     &          (abs(SMAT(2,2)).gt.abs(SMAT(2,3))) ) Then
c the WD(2) and ZD(i,2) correspond to gtens(2) and maxes(i,2)
        dtens(2)=WD(2)
        If (SMAT(2,2).gt.0.0_wp) Then
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(2,2) is larger than SMAT(2,1) and ' //
     &                     'SMAT(2,3) and is positive'
          End If
          Do i=1,3
            daxes(i,2)=ZD(i,2)
          End Do
        Else
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(2,2) is larger than SMAT(2,1) and ' //
     &                     'SMAT(2,3) and is negative'
          End If
          Do i=1,3
            daxes(i,2)=-ZD(i,2)
          End Do
        End If

      Else If ( (abs(SMAT(2,3)).gt.abs(SMAT(2,1))) .AND.
     &          (abs(SMAT(2,3)).gt.abs(SMAT(2,2))) ) Then
c the WD(2) and ZD(i,2) correspond to gtens(3) and maxes(i,3)
        dtens(2)=WD(3)
        If (SMAT(2,3).gt.0.0_wp) Then
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(2,3) is larger than SMAT(2,1) and ' //
     &                     'SMAT(2,2) and is positive'
          End If
          Do i=1,3
            daxes(i,2)=ZD(i,3)
          End Do
        Else
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(2,3) is larger than SMAT(2,1) and ' //
     &                     'SMAT(2,2) and is negative'
          End If
          Do i=1,3
            daxes(i,2)=-ZD(i,3)
          End Do
        End If
      End If

cccccccccccccccccccccccccccccccc
      If ( (abs(SMAT(3,1)).gt.abs(SMAT(3,2))) .AND.
     &     (abs(SMAT(3,1)).gt.abs(SMAT(3,3))) ) Then
c  the WD(3) and ZD(i,3) correspond to gtens(1) and maxes(i,1)
        dtens(3)=WD(1)
        If (SMAT(3,1).gt.0.0_wp) Then
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(3,1) is larger than SMAT(3,2) and ' //
     &                     'SMAT(3,3) and is positive'
          End If
          Do i=1,3
            daxes(i,3) = ZD(i,1)
          End Do
        Else
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(3,1) is larger than SMAT(3,2) and ' //
     &                     'SMAT(3,3) and is negative'
          End If
          Do i=1,3
            daxes(i,3) =-ZD(i,1)
          End Do
        End If

      Else If ( (abs(SMAT(3,2)).gt.abs(SMAT(3,1))) .AND.
     &          (abs(SMAT(3,2)).gt.abs(SMAT(3,3))) ) Then
c the WD(3) and ZD(i,3) correspond to gtens(2) and maxes(i,2)
        dtens(3)=WD(2)
        If (SMAT(3,2).gt.0.0_wp) Then
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(3,2) is larger than SMAT(3,1) and ' //
     &                     'SMAT(3,3) and is positive'
          End If
          Do i=1,3
            daxes(i,3) = ZD(i,2)
          End Do
        Else
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(3,2) is larger than SMAT(3,1) and ' //
     &                     'SMAT(3,3) and is negative'
          End If
          Do i=1,3
            daxes(i,3) =-ZD(i,2)
          End Do
        End If

      Else If ( (abs(SMAT(3,3)).gt.abs(SMAT(3,1))) .AND.
     &          (abs(SMAT(3,3)).gt.abs(SMAT(3,2))) ) Then
c the WD(3) and ZD(i,3) correspond to gtens(3) and maxes(i,3)
        dtens(3)=WD(3)
        If (SMAT(3,3).gt.0.0_wp) Then
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(3,3) is larger than SMAT(3,1) and ' //
     &                     'SMAT(3,2) and is positive'
          End If
          Do i=1,3
            daxes(i,3)=ZD(i,3)
          End Do
        Else
          If(IPRINT.gt.2) Then
            Write(6,'(a)') 'SMAT(3,3) is larger than SMAT(3,1) and ' //
     &                     'SMAT(3,2) and is negative'
          End If
          Do i=1,3
            daxes(i,3)=-ZD(i,3)
          End Do
        End If
      End If

      Call DGEMM_('N','N',3,3,3,1.0_wp,Z,3,daxes,3,0.0_wp,ZD2,3)

      diff12=abs(dtens(1)-dtens(2))
      diff23=abs(dtens(2)-dtens(3))

      If(diff12.gt.diff23) Then
        D_factor=1.5_wp*dtens(1)
        E_factor=(dtens(2)-dtens(3))/2.0_wp
      Else
        D_factor=1.5_wp*dtens(3)
        E_factor=(dtens(1)-dtens(2))/2.0_wp
      End If



      If(iprint.gt.2) Then
        Write(6,'(20X,A)') 'D-TENSOR:'
        Write(6,*)
        Write(6,'(10X,A,10X,3(F9.5,2X))') '|  xx    xy    xz  |',
     &   (DMATR(1,J),J=1,3)
        Write(6,'(10X,A,10X,3(F9.5,2X))') '|  yx    yy    yz  |',
     &   (DMATR(2,J),J=1,3)
        Write(6,'(10X,A,10X,3(F9.5,2X))') '|  zx    zy    zz  |',
     &   (DMATR(3,J),J=1,3)
        Write(6,*)
      End If

      If(iprint.GE.2) Then
        Write(6,*)
        Write(6,'(A)') 'D TENSOR:'
        Write(6,'(90a)') ('-',i=1,84),'|'
        Write(6,'(A,4x,A,27x,A,21x,a,3x,a)') 'MAIN VALUES','|',
     &           'MAIN ANISOTROPY AXES','|',
     &           'x , y , z  -- initial Cartesian axes'
        Write(6,'(85a,3x,a)') ('-',i=1,15),'|',('-',i=1,36),'|',
     &           ('-',i=1,31), '|', 'Xm, Ym, Zm -- main magnetic axes'
        Write(6,
     &    '(15x,a,4x,a,5x,a,8x,a,8x,a,4x,a,5x,a,9x,a,9x,a,5x,a,3x,a)')
     &           '|','|','Xm','Ym','Zm','|','x','y','z','|',
     &           'Xa, Ya, Za -- main anisotropy axes'
        Write(6,'(90a)') ('-',i=1,15),'|',('-',i=1,4),'|',('-',i=1,31),
     &           '|',('-',i=1,31),'|'
        Write(6,'(A,F9.3,A,3F10.6,1x,A,3F10.6,1x,A)') ' Dx =',dtens(1),
     &           ' | Xa |',(daxes(j,1),j=1,3),'|',(ZD2(j,1),j=1,3),'|'
        Write(6,'(A,F9.3,A,3F10.6,1x,A,3F10.6,1x,A)') ' Dy =',dtens(2),
     &           ' | Ya |',(daxes(j,2),j=1,3),'|',(ZD2(j,2),j=1,3),'|'
        Write(6,'(A,F9.3,A,3F10.6,1x,A,3F10.6,1x,A)') ' Dz =',dtens(3),
     &           ' | Za |',(daxes(j,3),j=1,3),'|',(ZD2(j,3),j=1,3),'|'
        Write(6,'(90a)') ('-',i=1,84),'|'
        Write(6,*)
        Write(6,'(A)') '2-nd order ZFS Hamiltonian:'
        Write(6,*)
        Write(6,'(A)') 'H_zfs^{2}= D * [S_{Za}^2 - S*(S+1)/3] + E * '//
     &           '[S_{Xa}^2 - S_{Ya}^2]'
        Write(6,*)
        Write(6,'(A)') 'Anisotropy parameters: D = 3/2 * Dz;  E = '//
     &           '(Dx-Dy)/2;'
        Write(6,'(a,F9.4)') 'D = ', D_factor
        Write(6,'(a,F9.4)') 'E = ', E_factor
      End If

      Call qExit('dmatrix')
      Return
      End
