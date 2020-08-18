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
      Subroutine pa_pseudo(M,S,dim,iopt,zfin,MF,SF,coord,iprint)
c dim - dimension of the pseuDospin
c M(3,dim,dim) -- magnetic moment in the initial basis
c S(3,dim,dim) -- spin moment in the initial basis
c iopt - option for choosing the quantization axis
c       = 1 : quantization axis is the main magnetic axis of the entire manIfold (dim)
c       = 2 : quantization axis is the main magnetic axis of the ground Doublet (low-lying two states)
c       = 3 : quantization axis is provided by the user by means of coord
c       = 4 : quantization axis is the unit matrix, i.e. the original Z axis
c  coord(3,3): matrix specIfying the rotation of initial coordinate system to the
c              coordinate system of used for determination of the quantization axis
c              by diagonalizing the Zeeman hamiltonian
c zfin - pseuDospin eigenfunctions
c MF(3,dim,dim) -- magnetic moment in the pseuDospin basis
c SF(3,dim,dim) -- spin moment in the pseuDospin basis


      Implicit None

      Integer, parameter :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer            :: dim,info,i,j,k,l,i1,i2,iopt
      Integer            :: iprint
      Real (kind=8) :: gtens(3),w(dim),maxes(3,3),det,FindDetR
      Real (kind=8) :: coord(3,3),coord2(3,3)
      Complex (kind=8) :: M( 3,dim,dim),S( 3,dim,dim),z(dim,dim)
      Complex (kind=8) :: MF(3,dim,dim),SF(3,dim,dim),zfin(dim,dim)
      Complex (kind=8) :: dipso2(3,dim,dim)
      Complex (kind=8) :: s_so2(3,dim,dim)
      Complex (kind=8) :: hzee(dim,dim)
      external FindDetR
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c set to zero important variables:
c  choose the quantization axis:
      If (    iopt .eq. 1 ) Then
        gtens=0.0_wp
        maxes=0.0_wp
        Call atens(M,dim,gtens,maxes,2)

      Else If( iopt .eq. 2 ) Then
        gtens=0.0_wp
        maxes=0.0_wp
        Call atens(M(1:3,1:2,1:2),2,gtens,maxes,2)

      Else If( iopt .eq. 3 ) Then
        Write(6,'(A)') 'User provided the COORD to the PA_PSEUDo:'

      Else If( iopt .eq. 4 ) Then
        Write(6,'(A)') 'COORD is set to unity'
        coord=0.0_wp
        Do i=1,3
        coord(i,i)=1.0_wp
        End Do
      Else
        Write(6,'(A)') 'check the iopt parameter to PA_PSEUDo!!!'
c        Call Abort()
      End If
cccccccccccccccccccccccccccccccccccc
      Write(6,'(A)') 'input "coord" matrix to PA_PSEUDo'
      Do i=1,3
      Write(6,'(3F20.14)') (coord(j,i),j=1,3)
      End Do
c check If "coord" is empty:
      Call rZeroMatrix(coord2,3)
      Do i=1,3
         Do j=1,3
         coord2(i,j)=coord(i,j)
         End Do
      End Do
      det=0.0_wp
      det=FindDetR(coord2,3)
      Write(6,'(A, f20.13)') 'det = ', det
      If ( abs(det-1.0_wp) .lt. 0.0001_wp) Then  ! 'coord is not empty'
      Call rZeroMatrix(maxes,3)
        Do i=1,3
           Do j=1,3
           maxes(i,j)=coord(i,j)
           End Do
        End Do
      Else                                     ! 'coord is empty'
      Call rZeroMatrix(coord,3)
        Do i=1,3
           Do j=1,3
           coord(i,j)=maxes(i,j)
           End Do
        End Do
      End If
      Write(6,'(A)') 'Employed  axes for diagonalization of the '//
     &               'Zeeman Hamiltonian:'
      Do i=1,3
      Write(6,'(3F20.14)') (maxes(j,i),j=1,3)
      End Do

      Call cZeroMoment( s_so2,dim)
      Call cZeroMoment(dipso2,dim)
      Do l=1,3
         Do i=1,dim
            Do j=1,dim
               Do k=1,3
      dipso2(l,i,j) = dipso2(l,i,j) + M(k,i,j) *
     &                          cmplx(maxes(k,l),0.0_wp,wp)
       s_so2(l,i,j) =  s_so2(l,i,j) + S(k,i,j) *
     &                          cmplx(maxes(k,l),0.0_wp,wp)
               End Do
            End Do
         End Do
      End Do

      Call cZeroMatrix(hzee,dim)
      Do i=1,dim
        Do j=1,dim
          hzee(i,j)=hzee(i,j) + (-1.0_wp,0.0_wp)*dipso2(3,i,j)
        End Do
      End Do
      info=0
      Call rZeroVector(w,dim)
      Call cZeroMatrix(z,dim)
      Call diag_c2(hzee,dim,info,w,z)
      If (info.ne.0) Then
      Write(6,'(5x,a)') 'diagonalization of the zeeman hamiltonian'//
     & ' failed.'
      Go To 199
      End If

      Call cZeroMatrix(zfin,dim)
      Call spin_phase(dipso2,dim,z,zfin)

      If(iprint.gt.2) Then
      Write(6,'(5X,A)') 'MAIN VALUES OF THE ZEEMAN HAMILTONIAN:'
      Write(6,*)
         If(MOD(dim,2).eq.1) Then
      Do I=1,dim
      Write(6,'(3X,A,I3,A,F17.3)') '|',(dim-1)/2+(1-I),'> = ',W(i)
      End Do
      Else
      Do I=1,dim
      Write(6,'(3X,A,I3,A,F17.3)') '|',(dim-1)-2*(I-1),'/2 > = ',W(i)
      End Do
         End If
      Write(6,*)
      Write(6,'(1X,A)') 'EIGENFUNCTIONS OF THE EFFECTIVE SPIN:'
      Write(6,*)
         If(mod(dim,2).eq.1) Then
      Do i=1,dim
      Write(6,'(10x,a,i2,a,10x,20(2f16.12,2x))') 'eigenvector of |',
     & (dim-1)/2+(1-i),' > :',(zfin(j,i),j=1,dim)
      End Do
         Else
      Do i=1,dim
      Write(6,'(10x,a,i2,a,10x,20(2f16.12,2x))') 'eigenvector of |',
     & (dim-1)-2*(i-1),'/2 > :',(zfin(j,i),j=1,dim)
      End Do
         End If
      End If ! printing of eigenfunctions with identical phase.

      Call cZeroMoment(MF,dim)
      Call cZeroMoment(SF,dim)
      Do l=1,3
        Do i=1,dim
          Do j=1,dim
            Do i1=1,dim
              Do i2=1,dim
      MF(l,i,j)=MF(l,i,j)+dipso2(l,i1,i2)*conjg(zfin(i1,i))*zfin(i2,j)
      SF(l,i,j)=SF(l,i,j)+ s_so2(l,i1,i2)*conjg(zfin(i1,i))*zfin(i2,j)
              End Do
            End Do
          End Do
        End Do
      End Do

      If(IPRINT.gt.2) Then
      Call prMom('PseudoSpin basis:  MAGNETIC MOMENT: MF :',MF,dim)
      Call prMom('PseudoSpin basis:      SPIN MOMENT: SF :',SF,dim)
      End If

 199  continue
      Return
      End
