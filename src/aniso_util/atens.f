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
      Subroutine atens(moment, dim, gtens, maxes, iprint)

      Implicit None
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
c Calling variables:
      Integer,          intent(in)  :: dim, iprint
      Complex (kind=wp),intent(in)  :: moment(3,dim,dim)
      Real (kind=wp),   intent(out) :: gtens(3)
      Real (kind=wp),   intent(out) :: maxes(3,3)
c----------------------------------------------
c  dim    -- size of the magnetic moment
c            dim = muliplicity of the pseuDospin ( 2*S+1, where S is the pseuDospin);
c  moment -- matrix of size (3,dim,dim) of the moment (magnetic, spin or angular)
c  gtens  -- array of size (3) keeping the main values of the A tensor ( sqrt(main_values) )
c  maxes  -- array of size (3,3) keeping the main axes of the A tensor Writen in
c            the right coordinate system (Determinant = +1)
c  iprint -- the print level of the Subroutine
c            iprint = 1 => no output
c            iprint = 2 => standard
c            iprint = 3 => print for debug
c----------------------------------------------
c local variables
      Integer           :: ic1, ic2, i, j, info
      Real (kind=wp)    :: A_TENS_TERM(3,3), W(3), MAIN(3), Z(3,3),
     &                     factor, Det_gtens, diff12, diff23,
     &                     ZR(3,3), dnorm
      Real (kind=wp)    :: dznrm2, FindDetR
      Complex (kind=wp) :: AC_TENS(3,3), trace
      External          :: dznrm2, FindDetR, trace

      Call qEnter('atens')

      dnorm=0.0_wp
      dnorm =  dznrm2(3*dim*dim, moment, 1 )

      If ( dnorm.eq.0._wp ) Then
         Write(6,'(A)') 'Norm of the magnetic moment is zero.'
         Write(6,'(A)') 'Returning the default (dummy) values'
         gtens=0.0_wp
         maxes=0.0_wp
         Do i=1,3
            maxes(i,i)=1.0_wp
         End Do
         Return
      End If

c initialization:
      Do I=1,3
         Do J=1,3
            AC_TENS(I,J)=(0.0_wp,0.0_wp)
            A_TENS_TERM(I,J)=0.0_wp
         End Do
      End Do

      Do ic1=1,3
         Do ic2=1,3
            Ac_tens(ic1,ic2)=trace(dim,moment(ic1,1:dim,1:dim),
     &                                 moment(ic2,1:dim,1:dim))
         End Do
      End Do

      Do ic1=1,3
         Do ic2=1,3
            A_TENS_TERM(ic1,ic2)=DBLE( Ac_tens( ic1, ic2 ) +
     &                                 Ac_tens( ic2, ic1 ) ) /2.0_wp
         End Do
      End Do

      factor=0.0_wp
      factor=12.0_wp/dble( dim**3-dim )
      Do ic1=1,3
         Do ic2=1,3
            A_TENS_TERM(ic1,ic2)=factor*A_TENS_TERM(ic1,ic2)
         End Do
      End Do

      If (iprint.gt.2) Then
         Write(6,'(/)')
         Write(6,'(5X,A)') 'A_TENS_TERM(ic1,ic2):'
         Write(6,*)
         Do ic1=1,3
            Write(6,'(5X,3(2F14.7,3x))') (A_TENS_TERM(ic1,ic2),ic2=1,3)
         End Do
      End If
c
C   Diagonalization of A_tens - g tensors
C
      Do I=1,3
      main(I)=0.0_wp
      w(I)=0.0_wp
         Do J=1,3
         z(I,J)=0.0_wp
         End Do
      End Do
      info=0

      Call DIAG_R2(A_TENS_TERM(1:3,1:3), 3, info, W(1:3),  Z(1:3,1:3))

      If (INFO.ne.0) Go To 199
      If ((w(1).lt.0._wp).AND.(w(2).lt.0._wp).AND.(w(3).lt.0._wp)) Then
      Write(6,'(2x,A)') 'ALL EIGENVALUES OF THE A-TENSOR ARE NEGATIVE'
      Write(6,'(2X,A)') 'THIS IS A VERY UNUSUAL SITUATION. PLEASE'//
     & 'CHECK MANUALLY '
      Write(6,'(2x,A)') 'THE FOLLOWING PART OF THE PSEUDoSPIN SECTION'
      Write(6,'(2x,A)') 'MUST BE DISREGARDED. THE RESULTS ARE NOT' //
     & 'TRUSTABLE.'
      Go To 199
      End If
c
      If (iprint.gt.2) Then
      Write(6,*)
      Write(6,'(4x,A)') 'A_TENS_TERM TENSOR:'
      Write(6,'(65a)') ('-',i=1,56),'|'
      Write(6,'(4x,A,4x,A,13x,A,5x,a,3x,a)') 'MAIN VALUES','|',
     & 'MAIN MAGNETIC AXES','|', 'x , y , z  -- initial Cartesian axes'
      Write(6,'(57a,3x,a)') ('-',i=1,19),'|',('-',i=1,36),'|',
     & 'Xm, Ym, Zm -- main magnetic axes'
      Write(6,'(19x,a,4x,a,5x,a,9x,a,9x,a,5x,a)') '|','|','x','y','z',
     & '|'
      Write(6,'(65a)') ('-',i=1,19),'|',('-',i=1,4),'|',('-',i=1,31),
     & '|'
      Write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gX = ',w(1),' | Xm |',
     & (z(j,1),j=1,3),'|'
      Write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',w(2),' | Ym |',
     & (z(j,2),j=1,3),'|'
      Write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',w(3),' | Zm |',
     & (z(j,3),j=1,3),'|'
      Write(6,'(65a)') ('-',i=1,56),'|'
      End If

      Do I=1,3
      If (W(I).lt.0.0_wp) Then
         W(I)= 1.0e-24_wp
      End If
         MAIN(i)=sqrt(W(i))
      End Do

C  Check the sign of the coordinate system. If CS is Left-handed,
C  Then change it to RIGHT-handed
      Det_gtens=0.0_wp
      Do I=1,3
         Do J=1,3
            ZR(I,J)=0.0_wp
            ZR(I,J)=Z(I,J)
         End Do
      End Do

      Det_gtens=FindDetR(ZR,3)

      If (Det_gtens.lt.0.0_wp) Then
         Do i=1,3
            Z(i,1)=-Z(i,1)
         End Do
         If (iprint.gt.2) Then
            Write(6,'(a)') 'The coordinate system is LEFT-handed.'
            Write(6,'(a)') 'It has been changed to RIGHT-handed'
         End If
      End If

      diff12=0.0_wp
      diff23=0.0_wp
      diff12=MAIN(2)-MAIN(1)
      diff23=MAIN(3)-MAIN(2)

      If (iprint.gt.2) Then
         Write(6,'(5x,a,3F19.15)') 'diff12 = ', diff12
         Write(6,'(5x,a,3F19.15)') 'diff23 = ', diff23
      End If

      Do i=1,3
         gtens(i)=0.0_wp
         Do j=1,3
            maxes(i,j)=0.0_wp
         End Do
      End Do
c set the main Z axis:
      If (diff12 > diff23) Then
         gtens(3)=MAIN(1)
         gtens(2)=MAIN(2)
         gtens(1)=MAIN(3)

         If (Z(3,1).ge.0._wp) Then
            Do i=1,3
            maxes(i,3)=Z(i,1)
            maxes(i,1)=Z(i,3)
            End Do
         Else If (Z(3,1).lt.0._wp) Then
            Do i=1,3
            maxes(i,3)=-1._wp*Z(i,1)
            maxes(i,1)=Z(i,3)
            End Do
         End If

          maxes(1,2)=maxes(2,3)*maxes(3,1)-maxes(2,1)*maxes(3,3)
          maxes(2,2)=maxes(1,1)*maxes(3,3)-maxes(1,3)*maxes(3,1)
          maxes(3,2)=maxes(1,3)*maxes(2,1)-maxes(1,1)*maxes(2,3)

      Else If (diff23 > diff12) Then
         gtens(3)=MAIN(3)
         gtens(2)=MAIN(2)
         gtens(1)=MAIN(1)

         If (Z(3,3).ge.0._wp) Then
            Do i=1,3
            maxes(i,3)=Z(i,3)
            maxes(i,1)=Z(i,1)
            End Do
         Else If (Z(3,3).lt.0.0_wp) Then
            Do i=1,3
            maxes(i,3)=-1.0_wp*Z(i,3)
            maxes(i,1)=Z(i,1)
            End Do
         End If

         maxes(1,2)=maxes(2,3)*maxes(3,1)-maxes(2,1)*maxes(3,3)
         maxes(2,2)=maxes(1,1)*maxes(3,3)-maxes(1,3)*maxes(3,1)
         maxes(3,2)=maxes(1,3)*maxes(2,1)-maxes(1,1)*maxes(2,3)
      Else !( diff23==diff12)
         ! this special case is isotropic:
         ! therefore assign main axes as close as to be to the cartesian xyz
         Do i=1,3
            gtens(i)=MAIN(i)
            maxes(i,i)=1.0_wp
         End Do

      End If
      If (iprint.gt.2) Then
         Write(6,*)
         Write(6,'(20X,A)') 'A-TENSOR:'
         Write(6,*)
         Write(6,'(10X,A,10X,3(F11.5,2X))') '|  xx    xy    xz  |',
     &                                      (A_TENS_TERM(1,i),i=1,3)
         Write(6,'(10X,A,10X,3(F11.5,2X))') '|  yx    yy    yz  |',
     &                                      (A_TENS_TERM(2,i),i=1,3)
         Write(6,'(10X,A,10X,3(F11.5,2X))') '|  zx    zy    zz  |',
     &                                      (A_TENS_TERM(3,i),i=1,3)
      End If

      If (iprint.GE.2) Then
         Write(6,*)
         Write(6,'(4x,A)') 'g TENSOR:'
         Write(6,'(65a)') ('-',i=1,56),'|'
         Write(6,'(4x,A,4x,A,13x,A,5x,a,3x,a)') 'MAIN VALUES','|',
     &                      'MAIN MAGNETIC AXES','|',
     &                      'x , y , z  -- initial Cartesian axes'
         Write(6,'(26a,3x,a)') ('-',i=1,19),'|',('-',i=1,4),'|',
     &                      '----- x ------- y ------- z ---|',
     &                      'Xm, Ym, Zm -- main magnetic axes'
c         Write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gX = ',gtens(1),' | Xm |',
c     &                                      (maxes(j,1),j=1,3),'|'
c         Write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',gtens(2),' | Ym |',
c     &                                      (maxes(j,2),j=1,3),'|'
c         Write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',gtens(3),' | Zm |',
c     &                                      (maxes(j,3),j=1,3),'|'
         Write(6,'(A,F18.14,A,3F18.14,1x,A)')
     &       ' gX = ',gtens(1),' | Xm |',(maxes(j,1),j=1,3),'|'
         Write(6,'(A,F18.14,A,3F18.14,1x,A)')
     &       ' gY = ',gtens(2),' | Ym |',(maxes(j,2),j=1,3),'|'
         Write(6,'(A,F18.14,A,3F18.14,1x,A)')
     &       ' gZ = ',gtens(3),' | Zm |',(maxes(j,3),j=1,3),'|'
         Write(6,'(65a)') ('-',i=1,56),'|'

         Call Add_Info('GTENS_MAIN',gtens,3,4)
      End If

 199  Continue
      Call qExit('atens')

      Return
      End Subroutine atens
