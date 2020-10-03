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
      Subroutine barrier( nBlock, dipIn, W, imanifold, NMULT, NDIM,
     &                    doPLOT, iprint)
c the present Subroutine computes the matrix elements of the transitions from states forming the blocking barrier;
c the states of opposite magnetization are given in the input, under keyword MLTP:
c the magnetic field is applied to each group of states delared at MLTP in order to form states of a definite
c projection of M on the quantization axis.
c  N --  dimension of the barrier

      Implicit None
#include "stdalloc.fh"
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)          :: nBlock, nMult, iprint, imanifold
      Integer, intent(in)          :: nDim(nMult)
      Real(kind=8), intent(in)    :: W(nBlock)
      Complex(kind=8), intent(in) :: dipIn(3,nBlock,nBlock)
      Logical, intent(in)          :: doPLOT
!-----------------------------------------------------------------------
      Integer          :: k,l,i,j,i1,i2,il,idim,Ifunct,j1,j2,iMult,ipar
      Integer          :: nb,maxmult
      Real(kind=8)    :: mave
      Character(len=90):: string2
      Character(len=5) :: s1,s2

      Integer, allocatable          :: ibas(:,:)
!                                          (nmult,nBlock)
      Real(kind=8), allocatable    :: gtens(:)
!                                           (3)
      Real(kind=8), allocatable    :: maxes(:,:), E(:)
!                                           (3,3), E(nmult)
      Real(kind=8), allocatable    :: wz(:,:)
!                                        (nmult,nBlock)
      Complex(kind=8), allocatable :: dipN(:,:,:)
!                                          (3,nBlock,nBlock)
      Complex(kind=8), allocatable :: dipN3(:,:,:)
!                                           (3,ndim(imanifold),ndim(imanifold))
      Complex(kind=8), allocatable :: CZ(:,:,:)
!                                        (nmult,nBlock,nBlock)
      Complex(kind=8), allocatable :: Ztr(:,:)
!                                         (nBlock,nBlock)
      Complex(kind=8), allocatable :: ML(:,:,:)
!                                        (3,nBlock,nBlock)
      Complex(kind=8), allocatable :: tmp(:,:)
!                                         (nBlock,nBlock)
      Complex(kind=8), allocatable :: MM(:)
!                                        (3)
      Complex(kind=8), allocatable :: dipso5(:,:,:,:,:)
!                                            (3,nmult,10,nmult,10)
!-----------------------------------------------------------------------

      If((nmult>0).and.(nBlock>0)) Then
         Call mma_allocate(ibas,nmult,nBlock,'ibas')
         Call mma_allocate(wz,nmult,nBlock,'wz')
         Call mma_allocate(cz,nmult,nBlock,nBlock,'cz')
         Call icopy( nmult*nBlock,[0],0,ibas,1)
         Call dcopy_(nmult*nBlock,[0.0_wp],0,wz,1)
         Call zcopy_(nmult*nBlock*nBlock,[(0.0_wp,0.0_wp)],0,CZ,1)
      End If

      If(nmult>0) Then
         Call mma_allocate(E,nmult,'E')
         Call mma_allocate(dipso5,3,nmult,10,nmult,10,'dipso5')
         Call dcopy_(nmult,[0.0_wp],0,E,1)
         Call zcopy_(3*nmult*10*nmult*10,[(0.0_wp,0.0_wp)],0,dipso5,1)
      End If

      If(nBlock>0) Then
         Call mma_allocate(dipN,3,nBlock,nBlock,'dipN')
         Call mma_allocate(ML,3,nBlock,nBlock,'ML')
         Call mma_allocate(Ztr,nBlock,nBlock,'Ztr')
         Call mma_allocate(tmp,nBlock,nBlock,'tmp')
         Call zcopy_(3*nBlock*nBlock,[(0.0_wp,0.0_wp)],0,dipN,1)
         Call zcopy_(3*nBlock*nBlock,[(0.0_wp,0.0_wp)],0,ML,1)
         Call zcopy_(nBlock*nBlock,[(0.0_wp,0.0_wp)],0,Ztr,1)
         Call zcopy_(nBlock*nBlock,[(0.0_wp,0.0_wp)],0,tmp,1)
      End If

      k=nDim(imanifold)
      If(k>0) Then
         Call mma_allocate(dipN3,3,k,k,'dipN3')
         Call zcopy_(3*k*k,[(0.0_wp,0.0_wp)],0, dipN3,1)
      End If

      Call mma_allocate(gtens,3,'gtens')
      Call mma_allocate(MM,3,'MM')
      Call mma_allocate(maxes,3,3,'maxes')
      Call dcopy_(3,[0.0_wp],0,gtens,1)
      Call dcopy_(3*3,[0.0_wp],0,maxes,1)
      Call zcopy_(3,[(0.0_wp,0.0_wp)],0,MM,1)


!-----------------------------------------------------------------------
      Do i=1,3
        maxes(i,i)=1.0_wp
      End Do
c  rotate the magnetic moment to the magnetic axes of the ground multiplet ( NDIM(1) )
      If(ndim(imanifold).le.1) Then
        Write(6,'(a,i2,a)') 'The manifold ',imanifold,
     &                      'was chosen for determination of the '//
     &                      'quantization axis.'
        Write(6,'(a     )') 'However, the size of this manifold is:'
        Write(6,'(a,i1  )') 'size = ',ndim(imanifold)
        Write(6,'(a     )') 'in this case, quantization axis will '//
     &                      'remain the original Z axis'
      Else
        Do i=1,ndim(imanifold)
          Do j=1,ndim(imanifold)
            Do l=1,3
              dipN3(l,i,j)=dipIn(l,i,j)
            End Do
          End Do
        End Do
        Call atens(dipN3, ndim(imanifold), gtens, maxes, 1 )
      End If
      Call rotmom2( dipIn(1:3,1:nBlock,1:nBlock), nBlock,
     &              maxes(1:3,1:3), dipN(1:3,1:nBlock,1:nBlock) )
      If(iprint.gt.2) Then
        Write(6,*)
        Write(6,'(10X,A)') 'Magnetic moment (dipIN) in the '//
     &                     'original coordinate system'
        Write(6,*)
        Write(6,'(A,11X,A,20X,A,20X,A,15x,A)') '< I | moment | J >',
     &                                         'projection = X',
     &                                         'projection = Y',
     &                                         'projection = Z',
     &                                   'ABS(< i | moment | j >)/3'
        Write(6,*)
        Do I=1,nBlock
          Do J=1,nBlock
            MAVE=0.0_wp
            MAVE=( abs(dipIN(1,I,J)) +
     &             abs(dipIN(2,I,J)) +
     &             abs(dipIN(3,I,J)) )/3.0_wp
            Write(6,'(A,i2,1x,A,i2,A, 3(2F16.11,2x),8X,F17.11)')
     &              '<',i,'| moment |',j,' >', (dipIN(L,I,J),L=1,3),MAVE
          End Do
        End Do
        Write(6,'(///)')
        Write(6,'(10X,A)') 'Magnetic moment (dipN) in the '//
     &                     'coordinate system of the magnetic axes'
        Write(6,*)
        Write(6,'(A,11X,A,20X,A,20X,A,15x,A)') '< I | moment | J >',
     &                                         'projection = X',
     &                                         'projection = Y',
     &                                         'projection = Z',
     &                                   'ABS(< i | moment | j >)/3'
        Write(6,*)
        Do I=1,nBlock
          Do J=1,nBlock
            MAVE=0.0_wp
            MAVE=( abs(dipN(1,I,J)) +
     &             abs(dipN(2,I,J)) +
     &             abs(dipN(3,I,J)) )/3.0_wp
            Write(6,'(A,i2,1x,A,i2,A, 3(2F16.11,2x),8X,F17.11)')
     &              '<',i,'| moment |',j,' >', (dipN(L,I,J),L=1,3),MAVE
          End Do
        End Do
      End If
c------------------------------------------------------------------------
c  determine the "CZ matrix":
      Ifunct=1
      Do il=1,nmult
        idim=ndim(il)
        If(idim.eq.0) Then
          Return
        Else If(idim.eq.1) Then
          CZ(il,1,1)=(1.0_wp,0.0_wp)
          WZ(il,1)=W(Ifunct)
        Else ! idim > 1

       Call pseudospin(  dipN(1:3,Ifunct:(Ifunct+ndim(il)-1),
     &                            Ifunct:(Ifunct+ndim(il)-1) ),
     &                      ndim(il),
     &                  CZ( il, 1:ndim(il), 1:ndim(il)),
     &                   3, 1, iprint)

       If (iPrint.gt.2) Call pa_prMat('barrier:  CZ',
     &      CZ(il,1:nDim(il),1:nDim(il)),nDim(il))
        End If
        Ifunct=Ifunct+ndim(il)
      End Do
c------------------------------------------------------------------------
c  compute the matrix elements between multiplets
      l=0
      ibas=0
      Do i=1,nmult
        Do j=1,ndim(i)
          l=l+1
          ibas(i,j)=l
        End Do
      End Do

      maxmult=0
      Do il=1,nmult
        If(ndim(il).gt.maxmult) Then
          maxmult=ndim(il)
        End If
      End Do
c build the transformation matrix Z(nBlock,nBlock)
      Call zcopy_(nBlock*nBlock,[(0.0_wp,0.0_wp)],0,ZTR,1)
      Do iMult=1,nmult
        Do j1=1,ndim(iMult)
          Do j2=1,ndim(iMult)
            i=ibas(iMult,j1)
            j=ibas(iMult,j2)
            Ztr( i, j )=CZ(iMult,j1,j2)
          End Do
        End Do
      End Do

      Call zcopy_(3*nBlock*nBlock,[(0.0_wp,0.0_wp)],0,ML,1)
      Do L=1,3
         Call zcopy_(nBlock*nBlock,[(0.0_wp,0.0_wp)],0,TMP,1)
         Call ZGEMM_('C','N',nBlock,nBlock,nBlock,(1.0_wp,0.0_wp),
     &                  Ztr,nBlock,
     &          dipN(L,:,:),nBlock,              (0.0_wp,0.0_wp),
     &                  TMP,nBlock )
         Call ZGEMM_('N','N',nBlock,nBlock,nBlock,(1.0_wp,0.0_wp),
     &                  TMP,nBlock,
     &                  Ztr,nBlock,              (0.0_wp,0.0_wp),
     &            ML(L,:,:),nBlock)
      End Do !L

      If(iprint.gt.2) Then
        Write(6,*)
        Write(6,'(10X,A)') 'Magnetic moment (ML) in the '//
     &                     'coordinate system of the magnetic axes'
        Write(6,*)
        Write(6,'(A,11X,A,20X,A,20X,A,15x,A)') '< I | moment | J >',
     &                                         'projection = X',
     &                                         'projection = Y',
     &                                         'projection = Z',
     &                                   'ABS(< i | moment | j >)/3'
        Write(6,*)
        Do I=1,nBlock
          Do J=1,nBlock
            MAVE=0.0_wp
            MAVE=( abs(ML(1,I,J)) +
     &             abs(ML(2,I,J)) +
     &             abs(ML(3,I,J)) )/3.0_wp
            Write(6,'(A,i2,1x,A,i2,A, 3(2F16.11,2x),8X,F17.11)')
     &              '<',i,'| moment |',j,' >', (ML(L,I,J),L=1,3),MAVE
          End Do
        End Do
      End If

      Call zcopy_(3*nmult*10*nmult*10,[(0.0_wp,0.0_wp)],0,DIPSO5,1)
      Do l=1,3
        Do i1=1,nmult
          Do j1=1,ndim(i1)
            Do i2=1,nmult
              Do j2=1,ndim(i2)
                i=ibas(i1,j1)
                j=ibas(i2,j2)
                dipso5(l,i1,j1,i2,j2)=ML(l,i,j)

c                write(6,'(A,5i3,A,2F20.12,3x,2I5)')
c     &                  'DIPSO5(',l,i1,j1,i2,j2,')=',
c     &                   dipso5(  l,i1,j1,i2,j2), i,j
              End Do
            End Do
          End Do
        End Do
      End Do

      If(doPLOT) then
         CALL plot_barrier(nBlock,nMult,nDIM,W,dipso5)
         Write(6,'(A)') 'The following files '
         Write(6,'(A)') '#-->  $WorkDir/BARRIER.plt'
         Write(6,'(A)') '#-->  $WorkDir/BARRIER_ENE.dat'
         Write(6,'(A)') '#-->  $WorkDir/BARREIR_TME.dat'
         Write(6,'(A)') '#-->  $WorkDir/BARREIR.plt'
         Write(6,'(A)') 'Have been generated successfully.'
      End If


ccccccccccccccccccccccccccccccccccccc


c   new print-out code :
      Write(6,*)
      Write(6,'(100A)') ('%',i=1,95)
      Write(6,'(15X,A)') 'AB INITIO BLOCKING BARRIER'
      Write(6,'(100A)') ('%',i=1,95)
      Write(6,'(A)') 'please, acknowledge the fact that the '//
     &               'information printed below provides'
      Write(6,'(A)') 'only a qualitative relaxation path of a '//
     &               'single-molecule magnet'
      Write(6,*)
c check the parity of all manIfolds:
      ipar=0
      i=1
      Do il=1,nmult
        ipar=ipar + mod(ndim(il),2)
        If(mod(ndim(il),2).eq.1)  i= i + 1
!        count the number of odd manifolds
c        Write(6,'(3(A,i2,2x))') 'il=',il,'ipar=',ipar,'i=',i
      End Do
      If (i.gt.1) Then
        ipar=ipar/(i-1)
      End If

      Ifunct=0
      Do il=1,nmult
        E(il)=0.0_wp
        Do i=1,ndim(il)
          E(il)=E(il)+W(Ifunct+i)/ndim(il)
        End Do
        Ifunct=Ifunct+ndim(il)
      End Do
      Write(6,'(A)') 'Zeeman eigenstates:'
c
c the convention to label states in the blocing barrier is the following:
c   Size of the             labelling scheme
c   manifold
c      1    ---------------                     0
c      2    ---------------                 1+,    1-
c      3    ---------------                 1+, 0, 1-
c      4    ---------------             2+, 1+,    1-, 2-
c      5    ---------------             2+, 1+, 0, 1-, 2-
c      6    ---------------         3+, 2+, 1+,    1-, 2-, 3-
c      7    ---------------         3+, 2+, 1+, 0, 1-, 2-, 3-
c      8    ---------------     4+, 3+, 2+, 1+,    1-, 2-, 3-, 4-
c      9    ---------------     4+, 3+, 2+, 1+, 0, 1-, 2-, 3-, 4-
c      9    --------------- 5+, 4+, 3+, 2+, 1+,    1-, 2-, 3-, 4-, 5-
c      9    --------------- 5+, 4+, 3+, 2+, 1+, 0, 1-, 2-, 3-, 4-, 5-
c    label(il,istate) is a Character*5 array of dimension (nmult,10)
c
c  notation:   multiplet . Lbl+
c  two Characters are assigned for multiplet
c

cccccccccccccccccccccccccccccccccccccccccccccccccc      cccc

      If(ipar.eq.0) Then
c  all multiplets have the same parity
        Write(string2, '(a, i2, a)') '(A,', maxmult, 'A,A)'
        Write(6,string2) '-------',
     &    ('--------------',i=1,maxmult),'----------------'
        If(mod(maxmult,2).eq.0) Then
           Write(string2, '(a, i2, a,i2,a)')
     &        '(A,',maxmult/2,'(A,i2,a),',
     &              maxmult/2,'(A,i2,a),a)'
           Write(6,string2) ' Mult.|',
     &       ('     ',i,'+     |',i=maxmult/2,1,-1),
     &       ('     ',i,'-     |',i=1,maxmult/2,1),
     &        '    E (cm-1)   |'
        Else If(mod(maxmult,2).eq.1) Then
           Write(string2, '(a,i2,a,i2,a)') '(A,',
     &                                      (maxmult-1)/2,'(A,i2,a),a,',
     &                                      (maxmult-1)/2,'(A,i2,a),a)'
           Write(6,string2)
     &        ' Mult.|',
     &       ('     ',i,'+     |',i=int((maxmult-1)/2),1,-1),
     &        '      0      |',
     &       ('     ',i,'-     |',i=1,int((maxmult-1)/2),1),
     &        '    E (cm-1)   |'
        End If
      Write(string2, '(a, i2, a)') '(A,', maxmult, 'A,A)'
      Write(6,string2) '------|',
     & ('-------------|', i=1,maxmult),
     &  '---------------|'

        Do il=1,nmult
      If(ndim(il).lt.maxmult) Then
        nb=int((maxmult-ndim(il))/2)
           Write(string2, '(3(a,i2),a)') '(2x,i2,a,',nb,'A,',
     &            ndim(il),'(F11.7,a),',nb,'A,F13.7,a)'
           Write(6,string2)
     &          il,'. | ',
     &          ('            | ',i=1,nb),
     & ( dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)),
     &          ('            | ',i=1,nb),
     &          E(il),' |'
      Else
         Write(string2, '(a, i2, a)') '(2x,i2,a,',maxmult,
     &                             '(F11.7,a),F13.7,a)'
         Write(6,string2) il,'. | ',
     &   ( dble(dipso5(3,il,i,il,i)),
     &    ' | ',i=1,ndim(il)), E(il),' |'
      End If
      End Do !il
      Write(string2, '(a, i2, a)') '(A,', maxmult, 'A,A)'
      Write(6,string2) '------|',
     & ('-------------|', i=1,maxmult),
     &  '---------------|'

      Else !ipar
c  multiplets have different parity (even and odd)
        Write(string2, '(a, i2, a, i2, a)')
     &  '(A,', maxmult/2, 'A,A,', maxmult/2,'A,A)'
        Write(6,string2) '-------',
     &    ('--------------',i=1,maxmult/2),
     &     '--------------',
     &    ('--------------',i=1,maxmult/2),
     &     '----------------'
         If(mod(maxmult,2).eq.0) Then  !maxmult = even
            Write(string2, '(a, i2, a, i2, a)')
     &           '(A,',maxmult/2,'(A,i2,a),a,'
     &                ,maxmult/2,'(A,i2,a),a)'
            Write(6,string2) ' Mult.|',
     &        ('     ',i,'+     |',i=maxmult/2,1,-1),
     &         '      0      |',
     &        ('     ',i,'-     |',i=1,maxmult/2,1),
     &         '    E (cm-1)   |'
            Write(string2, '(a, i2, a, i2, a)')
     &         '(A,', maxmult/2, 'A,A,',maxmult/2,'A,A)'
            Write(6,string2) '------|',
     &        ('-------------|', i=1,maxmult/2),
     &         '-------------|',
     &        ('-------------|', i=1,maxmult/2),
     &         '---------------|'

           Do il=1,nmult
             If(mod(ndim(il),2).eq.0) Then
!               il = even > the same parity as maxmult
                nb=(maxmult-ndim(il))/2
                If(nb.eq.0) Then !ndim(il) => maxmult
                Write(string2, '(2(a,i2),a)') '(2x,i2,a,',
     &                  ndim(il)/2,'(F11.7,a),    A,',
     &                  ndim(il)/2,'(F11.7,a),F13.7,a)'
                Write(6,string2) il,'. | ',
     &   ( dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),
     &     '            | ' ,
     &   ( dble(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,
     &    ndim(il)), E(il),' |'
                Else !nb>0, Integer
                Write(string2, '(4(a,i2),a)') '(2x,i2,a,',nb,'A,',
     &                            ndim(il)/2,'(F11.7,a),    A,',
     &                            ndim(il)/2,'(F11.7,a),',
     &                                      nb,'A,F13.7,a)'
                Write(6,string2) il,'. | ',
     &          ('            | ',i=1,nb),
     &          ( dble(dipso5(3,il,i,il,i)),
     &                    ' | ',i=1,ndim(il)/2),
     &           '            | ' ,
     &          ( dble(dipso5(3,il,i,il,i)),
     &                    ' | ',i=1+ndim(il)/2,ndim(il)),
     &          ('            | ',i=1,nb),
     &          E(il),' |'
                End If
             Else ! il = odd  < maxmult
               nb=(maxmult+1-ndim(il))/2
                Write(string2, '(4(a,i2),a)') '(2x,i2,a,',nb,'A,',
     &              ndim(il),'(F11.7,a),',nb,'A,F13.7,a)'
                Write(6,string2) il,'. | ',
     &          ('            | ',i=1,nb),
     &   ( dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)),
     &          ('            | ',i=1,nb),
     &          E(il),' |'
             End If
           End Do


         Else ! maxmult = odd
           Write(string2, '(a, i2, a, i2, a)') '(A,',
     &                           (maxmult-1)/2,'(a,i2,a),a,',
     &                           (maxmult-1)/2,'(a,i2,a),A)'
           Write(6,string2)
     &        ' Mult.|',
     &       ('     ',i,'+     |',i=  int((maxmult-1)/2),1,-1),
     &        '      0      |',
     &       ('     ',i,'-     |',i=1,int((maxmult-1)/2),1),
     &        '    E (cm-1)   |'
           Write(string2, '(a, i2, a, i2, a)')
     &        '(A,', maxmult/2, 'A,A,',maxmult/2,'A,A)'
           Write(6,string2) '------|',
     &       ('-------------|', i=1,maxmult/2),
     &        '-------------|',
     &       ('-------------|', i=1,maxmult/2),
     &        '---------------|'

            Do il=1,nmult
             If(mod(ndim(il),2).eq.1) Then
!               il = odd,  the parity of il = maxmult
                nb=(maxmult-ndim(il))/2
                If(nb.eq.0) Then !ndim(il) => maxmult
                Write(string2, '(2(a,i2),a)') '(2x,i2,a,',
     &                  ndim(il),'(F11.7,a),F13.7,a)'
                Write(6,string2) il,'. | ',
     &   ( dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)),
     &          E(il),' |'
                Else !nb>0, Integer
                Write(string2, '(4(a,i2),a)') '(2x,i2,a,',
     &          nb,'A,', ndim(il),'(F11.7,a),', nb,'A,F13.7,a)'
                Write(6,string2) il,'. | ',
     &          ('            | ',i=1,nb),
     &   ( dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)),
     &          ('            | ',i=1,nb),
     &          E(il),' |'
                End If

             Else ! il = even  < maxmult
               nb=(maxmult-1-ndim(il))/2
                If(nb.eq.0) Then !ndim(il) => maxmult
                Write(string2, '(2(a,i2),a)') '(2x,i2,a,',
     &                  ndim(il)/2,'(F11.7,a),    A,',
     &                  ndim(il)/2,'(F11.7,a),F13.7,a)'
                Write(6,string2) il,'. | ',
     &   ( dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),
     &     '            | ' ,
     &   ( dble(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,
     &    ndim(il)), E(il),' |'
                Else !nb>0, Integer
                Write(string2, '(4(a,i2),a)') '(2x,i2,a,',nb,'A,',
     &                            ndim(il)/2,'(F11.7,a),    A,',
     &                            ndim(il)/2,'(F11.7,a),',
     &                                      nb,'A,F13.7,a)'
                Write(6,string2) il,'. | ',
     &          ('            | ',i=1,nb),
     &   ( dble(dipso5(3,il,i,il,i)),' | ',i=1,ndim(il)/2),
     &           '            | ' ,
     &   ( dble(dipso5(3,il,i,il,i)),' | ',i=1+ndim(il)/2,
     &                                                 ndim(il)),
     &          ('            | ',i=1,nb),E(il),' |'
                End If
             End If
            End Do

         End If
      Write(string2, '(a, i2, a, i2, a)')
     & '(A,', maxmult/2, 'A,A,',maxmult/2,'A,A)'
      Write(6,string2) '------|',
     & ('-------------|', i=1,maxmult/2),
     &  '-------------|',
     & ('-------------|', i=1,maxmult/2),
     &  '---------------|'
      End If ! ipar , line 271
      Write(6,*)
c  printing of all off-diagonal matrix elements
      Write(6,'(A)') 'Matrix elements of the magnetic moment '//
     &               'connecting Zeeman eigenstates'
      Write(6,'(A)') 'The average is done according to the formula:'
      Write(6,'(A)') '<i|m|j> = ( ABS(<i|m_X|j>) + ABS(<i|m_Y|j>) + '//
     &                           'ABS(<i|m_Z|j>) ) / 3'
      Write(6,*)
      Write(6,'(A)') 'MATRIX ELEMENTS BETWEEN STATES WITH '//
     & 'OPPOSITE MAGNETIZATION'
      Write(6,'(A)') 'in cases with even number of electrons, '//
     &                'these values cannot be used'
      Write(6,'(A)') 'check the tunnelling splitting instead'
      Write(6,'(4A)') '-------','-------------------------',
     & '----------------------------------------',
     & '------------------------'
      Write(6,'(A,5x,A,5x,A,13x,A,13x,A,8x,a,8x,a)') ' Mult.|',
     & 'Matrix Element','|','Complex VALUE','|','AVERAGE','|'
      Write(6,'(4A)') '------|','------------------------|',
     & '---------Real------------imaginary-----|',
     & '-----------------------|'

      Do il=1,nmult
      If(ndim(il).eq.1) Then
               write(s1,'(i2,A1,A1,A1)') il,'.','0',' '
               write(s2,'(i2,A1,A1,A1)') il,'.','0',' '
               Call prbar(il,s1,s2,dipso5(1:3,il,1,il,1))
      Else
c         If( mod(ndim(il),2)==0 ) Then   ! even multiplicity

        Do i=1,int(ndim(il)/2)
      If(i.gt.1) Write(6,'(4A)') '      |',
     & '------------------------|',
     & '---------------------------------------|',
     & '-----------------------|'
            write(s1,'(i2,A1,i1,A1)') il,'.',i,'+'
            Do j=i,int(ndim(il)/2)
              write(s2,'(i2,A1,i1,A1)') il,'.',j,'+'
c            !Call prbar(il,s1,s2,dipso5(1:3,il,i,il,ndim(il)-i+1))
              Call prbar(il,s1,s2,dipso5(1:3,il,i,il, j ) )
            End Do

            Do j=i,int(ndim(il)/2)
              write(s2,'(i2,A1,i1,A1)') il,'.',j,'-'
c            !Call prbar(il,s1,s2,dipso5(1:3,il,i,il,ndim(il)-i+1))
              Call prbar(il,s1,s2,dipso5(1:3,il,i,il, ndim(il)-j+1 ) )
            End Do
        End Do

c        !Else ! mod(ndim(il),2)==1 , i.e. odd multiplicity
c        End If

      End If
      Write(6,'(4A)') '------|','------------------------|',
     & '---------------------------------------|',
     & '-----------------------|'
      End Do !il
      Write(6,*)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Write(6,'(10A)') ('##########',i=1,10)
      Do k=1,nmult-1
      Write(6,*)
      Write(6,'(A,i2,a)') 'MATRIX ELEMENTS BETWEEN STATES ARISING '//
     & 'FROM NEIGHBORING MULTIPLETS: I -> I+',k,' :'
      Write(6,'(4A)') '-------','-------------------------',
     & '----------------------------------------',
     & '------------------------'
      Write(6,'(A,5x,A,5x,A,13x,A,13x,A,8x,a,8x,a)') ' Mult.|',
     & 'Matrix Element','|','Complex VALUE','|','AVERAGE','|'
      Write(6,'(4A)') '------|','------------------------|',
     & '---------Real------------imaginary-----|',
     & '-----------------------|'
      Do il=1,nmult-k
      If( mod(ndim(il),2) .eq. 0 )  Then
       Do i=ndim(il)/2,1,-1
          If(mod(ndim(il+k),2) .eq. 0 ) Then
                Do j=ndim(il+k)/2,1,-1

               write(s1,'(i2,A1,i1,A1)') il  ,'.',i,'+'
               write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
               Call prbar(il,s1,s2,dipso5(1:3,il,i,il+k,j))

                End Do

                Do j=1,ndim(il+k)/2

               write(s1,'(i2,A1,i1,A1)') il  ,'.',i,'+'
               write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
               Call prbar(il,s1,s2,dipso5(1:3,il,i,il+k,ndim(il+k)-j+1))

                End Do

          Else !ndim(il+k) is odd
             If(ndim(il+k).eq.1) Then

               write(s1,'(i2,A1,i1,A1)') il  ,'.', i ,'+'
               write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
               Call prbar(il,s1,s2,dipso5(1:3,il,i,il+k,1))

             Else

                Do j=(ndim(il+k)-1)/2,1,-1
                   write(s1,'(i2,A1,i1,A1)') il  ,'.',i,'+'
                   write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
                   Call prbar(il,s1,s2,dipso5(1:3,il,i,il+k,j))
                End Do

                write(s1,'(i2,A1,i1,A1)') il  ,'.', i ,'+'
                write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
                MM(1:3)=dipso5(1:3,il,i,il+k,(ndim(il+k)+1)/2)
                Call prbar(il,s1,s2,MM(1:3))

                Do j=1,(ndim(il+k)-1)/2
                  write(s1,'(i2,A1,i1,A1)') il  ,'.', i ,'+'
                  write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
                  MM(1:3)=dipso5(1:3,il,i,il+k,ndim(il+k)-j+1)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do

             End If ! ndim(il+k) = 1
          End If ! ndim(il+k)
      If(i.gt.1) Write(6,'(4A)') '      |','------------------------|',
     & '---------------------------------------|',
     & '-----------------------|'
       End Do !i

      Else ! ndim(il) = odd

       If(ndim(il).eq.1) Then
          If(mod(ndim(il+k),2) .eq. 0 ) Then
                Do j=ndim(il+k)/2,1,-1
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,i1,A1)') il+k,'.', j ,'+'
                  MM(1:3)=dipso5(1:3,il,1,il+k,j)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do

                Do j=1,ndim(il+k)/2
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,i1,A1)') il+k,'.', j ,'-'
                  MM(1:3)=dipso5(1:3,il,1,il+k,ndim(il+k)-j+1)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do

          Else !ndim(il+k) is odd
             If(ndim(il+k).eq.1) Then
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
                  MM(1:3)=dipso5(1:3,il,1,il+k,1)
                  Call prbar(il,s1,s2,MM(1:3))
             Else
                Do j=(ndim(il+k)-1)/2,1,-1
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,i1,A1)') il+k,'.', j ,'+'
                  MM(1:3)=dipso5(1:3,il,1,il+k,j)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
                  MM(1:3)=dipso5(1:3,il,1,il+k,(ndim(il+k)+1)/2)
                  Call prbar(il,s1,s2,MM(1:3))

                Do j=1,(ndim(il+k)-1)/2
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,i1,A1)') il+k,'.', j ,'-'
                  MM(1:3)=dipso5(1:3,il,1,il+k,ndim(il+k)-j+1)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do
             End If ! ndim(il+k) = 1
          End If ! ndim(il+k), parity

       Else !ndim(il) > 1, odd

       Do i=(ndim(il)-1)/2,1,-1
          If(mod(ndim(il+k),2) .eq. 0 ) Then
                Do j=ndim(il+k)/2,1,-1
                  write(s1,'(i2,A1,i1,A1)') il  ,'.',i,'+'
                  write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
                  MM(1:3)=dipso5(1:3,il,i,il+k,j)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do
                Do j=1,ndim(il+k)/2
                  write(s1,'(i2,A1,i1,A1)') il  ,'.',i,'+'
                  write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'-'
                  MM(1:3)=dipso5(1:3,il,i,il+k,ndim(il+k)-j+1)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do
          Else !ndim(il+k) is odd
             If(ndim(il+k).eq.1) Then
                  write(s1,'(i2,A1,i1,A1)') il  ,'.', i ,'+'
                  write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
                  MM(1:3)=dipso5(1:3,il,i,il+k,1)
                  Call prbar(il,s1,s2,MM(1:3))
             Else
                Do j=(ndim(il+k)-1)/2,1,-1
                  write(s1,'(i2,A1,i1,A1)') il  ,'.',i,'+'
                  write(s2,'(i2,A1,i1,A1)') il+k,'.',j,'+'
                  MM(1:3)=dipso5(1:3,il,i,il+k,j)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do
                  write(s1,'(i2,A1,i1,A1)') il  ,'.', i ,'+'
                  write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
                  MM(1:3)=dipso5(1:3,il,i,il+k,(ndim(il+k)+1)/2)
                  Call prbar(il,s1,s2,MM(1:3))
                Do j=1,(ndim(il+k)-1)/2
                  write(s1,'(i2,A1,i1,A1)') il  ,'.', i ,'+'
                  write(s2,'(i2,A1,i1,A1)') il+k,'.', j ,'-'
                  MM(1:3)=dipso5(1:3,il,i,il+k,ndim(il+k)-j+1)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do
             End If ! ndim(il+k) = 1
          End If ! ndim(il+k), parity
      Write(6,'(4A)') '      |','------------------------|',
     & '---------------------------------------|',
     & '-----------------------|'
       End Do ! i


          If(mod(ndim(il+k),2) .eq. 0 ) Then
                Do j=ndim(il+k)/2,1,-1
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,i1,A1)') il+k,'.', j ,'+'
                  MM(1:3)=dipso5(1:3,il,(ndim(il)+1)/2,il+k,j)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do
                Do j=1,ndim(il+k)/2
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,i1,A1)') il+k,'.', j ,'-'
                  MM(1:3)=dipso5(1:3,il  ,(ndim(il)+1)/2,
     &                               il+k, ndim(il+k)-j+1)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do
          Else !ndim(il+k) is odd
             If(ndim(il+k).eq.1) Then
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
                  MM(1:3)=dipso5(1:3,il,(ndim(il)+1)/2,il+k,1)
                  Call prbar(il,s1,s2,MM(1:3))
             Else
                Do j=(ndim(il+k)-1)/2,1,-1
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,i1,A1)') il+k,'.', j ,'+'
                  MM(1:3)=dipso5(1:3,il,(ndim(il)+1)/2,il+k,j)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,A1,A1)') il+k,'.','0',' '
                  MM(1:3)=dipso5(1:3,il,(ndim(il)+1)/2,
     &                             il+k,(ndim(il+k)+1)/2)
                  Call prbar(il,s1,s2,MM(1:3))
                Do j=1,(ndim(il+k)-1)/2
                  write(s1,'(i2,A1,A1,A1)') il  ,'.','0',' '
                  write(s2,'(i2,A1,i1,A1)') il+k,'.', j ,'-'
                  MM(1:3)=dipso5(1:3,il,(ndim(il)+1)/2,
     &                             il+k,ndim(il+k)-j+1)
                  Call prbar(il,s1,s2,MM(1:3))
                End Do
             End If ! ndim(il+k) = 1
          End If ! ndim(il+k), parity

       End If ! ndim(il), size

      End If !ndim(il), parity
      Write(6,'(4A)') '------|','------------------------|',
     & '---------------------------------------|',
     & '-----------------------|'
      End Do ! il
      End Do ! k

!-----------------------------------------------------------------------
      If((nmult>0).and.(nBlock>0)) Then
         Call mma_deallocate(ibas)
         Call mma_deallocate(wz)
         Call mma_deallocate(cz)
      End If

      If(nmult>0) Then
         Call mma_deallocate(E)
         Call mma_deallocate(dipso5)
      End If

      If(nBlock>0) Then
         Call mma_deallocate(dipN)
         Call mma_deallocate(ML)
         Call mma_deallocate(Ztr)
         Call mma_deallocate(tmp)
      End If

      k=nDim(imanifold)
      If(k>0) Then
         Call mma_deallocate(dipN3)
      End If

      Call mma_deallocate(gtens)
      Call mma_deallocate(MM)
      Call mma_deallocate(maxes)


      Return
      End



      Subroutine prbar(ist,s1,s2,M)

      Implicit none
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)          :: ist
      Character(len=5), intent(in) :: s1, s2
      Complex(kind=8), intent(in) :: M(3)
!     local
      Character(len=30) :: fx, fy, fz
      Character(len=40) :: f1, f2
      Real(kind=8)     :: R

      write(fx,'(i2,5a)') ist,'. | <',s1,' | mu_X |',s2,' > |'
      write(fy,'(i2,5a)') ist,'. | <',s1,' | mu_Y |',s2,' > |'
      write(fz,'(i2,5a)') ist,'. | <',s1,' | mu_Z |',s2,' > |'
      R=0.0_wp
      R=( abs(M(1))+abs(M(2))+abs(M(3)) )/3.0_wp

      f1='(2x,a,2E19.11,1x,A,      23x,A)'
      f2='(2x,a,2E19.11,1x,A,E22.12,1x,A)'
      Write(6,f1) fx,M(1),'|',  '|'
      Write(6,f2) fy,M(2),'|',R,'|'
      Write(6,f1) fz,M(3),'|',  '|'
      Return
      End subroutine prbar

