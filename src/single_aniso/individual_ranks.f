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
      Subroutine individual_ranks(nDIMCF,BC,BS,Hinit,LJ,iprint)
      Implicit None
      Integer, Parameter          :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)         :: nDIMCF,iprint
      Real(kind=8),intent(in)    :: BC(nDIMcf,0:nDIMcf)
      Real(kind=8),intent(in)    :: BS(nDIMcf,0:nDIMcf)
      Complex(kind=8),intent(in) :: Hinit(nDIMcf,nDIMcf)
      Character(len=1),intent(in) :: LJ
      ! local variables:
      Integer                     :: i,j,k,n,m,ik,jEnd,nfields,info,
     &                               ikmax,ip,ir,iq
      Integer, allocatable        :: rankKQ(:) ! nDIMCF*(2*nDIMCF+1) )
      Integer, allocatable        :: projKQ(:) ! nDIMCF*(2*nDIMCF+1) )
      Real(kind=8)               :: Tnrm, TnrmKQ, wt
      Real(kind=8)               :: RnrmKQ( nDIMcf,-nDIMcf:nDIMcf )
      Real(kind=8), allocatable  :: ListKQ(:) ! nDIMCF*(2*nDIMCF+1) )
      Real(kind=8), allocatable  :: Rnrm(:) !nDIMCF)
      Real(kind=8), allocatable  :: Snrm(:) !nDIMCF)
      Real(kind=8), allocatable  :: Wk(:,:) !nDIMCF,nDIMCF)
      Real(kind=8), allocatable  :: Ws(:,:) !nDIMCF,nDIMCF)
      Real(kind=8), allocatable  :: Winit(:) !nDIMCF)
      Real(kind=8), external     :: dznrm2_

      Complex(kind=8)              :: zf, redME
      Complex(kind=8), allocatable :: O(:,:)  !nDIMCF,nDIMCF)! real ITO
      Complex(kind=8), allocatable :: W(:,:)  !nDIMCF,nDIMCF)! imag ITO
      Complex(kind=8), allocatable :: HCF(:,:,:) !nDIMCF,nDIMCF,nDIMCF)
      Complex(kind=8), allocatable :: HCFS(:,:,:)!nDIMCF,nDIMCF,nDIMCF)
      Complex(kind=8), allocatable :: Zk(:,:,:)  !nDIMCF,nDIMCF,nDIMCF)
      Complex(kind=8), allocatable :: Zs(:,:,:)  !nDIMCF,nDIMCF,nDIMCF)
      Complex(kind=8), allocatable :: Zinit(:,:) !nDIMCF,nDIMCF)
      Complex(kind=8), allocatable :: HKQ(:,:),Zkq(:,:)

      Real(kind=8)             :: dznrm2
      External                  :: dznrm2
      Character(len=16)         :: field(8)
      Character(len=6)          :: iprog
!-----------------------------------------------------------------------
      Call mma_allocate(Rnrm,nDIMCF,'Rnrm')
      Call mma_allocate(Snrm,nDIMCF,'Snrm')
      Call mma_allocate(Wk,nDIMCF,nDIMCF,'Wk')
      Call mma_allocate(Ws,nDIMCF,nDIMCF,'Ws')
      Call mma_allocate(Winit,nDIMCF,'Winit')
      Call mma_allocate(ListKQ,nDIMCF*(2*nDIMCF+1),'ListKQ')
      Call mma_allocate(rankKQ,nDIMCF*(2*nDIMCF+1),'rankKQ')
      Call mma_allocate(projKQ,nDIMCF*(2*nDIMCF+1),'projKQ')

      Call mma_allocate(O,nDIMCF,nDIMCF,'O')
      Call mma_allocate(W,nDIMCF,nDIMCF,'W')
      Call mma_allocate(HKQ,nDIMCF,nDIMCF,'HKQ')
      Call mma_allocate(ZKQ,nDIMCF,nDIMCF,'ZKQ')
      Call mma_allocate(HCF,nDIMCF,nDIMCF,nDIMCF,'HCF')
      Call mma_allocate(HCFS,nDIMCF,nDIMCF,nDIMCF,'HCFS')
      Call mma_allocate(Zk,nDIMCF,nDIMCF,nDIMCF,'Zk')
      Call mma_allocate(Zs,nDIMCF,nDIMCF,nDIMCF,'Zs')
      Call mma_allocate(Zinit,nDIMCF,nDIMCF,'Zinit')

      Call dcopy_(nDIMCF,[0.0_wp],0,Rnrm,1)
      Call dcopy_(nDIMCF,[0.0_wp],0,Snrm,1)
      Tnrm= 0.0_wp
      Call dcopy_(nDIMCF*nDIMCF,[0.0_wp],0,Wk,1)
      Call dcopy_(nDIMCF*nDIMCF,[0.0_wp],0,Ws,1)
      Call dcopy_(nDIMCF,[0.0_wp],0,Winit,1)

      Call zcopy_(nDIMCF*nDIMCF,[(0.0_wp,0.0_wp)],0,Zinit,1)
      Call zcopy_(nDIMCF*nDIMCF*nDIMCF,[(0.0_wp,0.0_wp)],0,HCF,1)
      Call zcopy_(nDIMCF*nDIMCF*nDIMCF,[(0.0_wp,0.0_wp)],0,HCFS,1)
      Call zcopy_(nDIMCF*nDIMCF*nDIMCF,[(0.0_wp,0.0_wp)],0,Zk,1)
      Call zcopy_(nDIMCF*nDIMCF*nDIMCF,[(0.0_wp,0.0_wp)],0,Zs,1)
!-----------------------------------------------------------------------
      ! re-construct the  initial CF matrix:
      Do N=2,nDIMcf-1,2
        Do M=0,N
          Call Liviu_ESO(nDIMcf,N,M,O,W,redME)
          If(M==0) Then
            zf=cmplx(BC(N,0),0.0_wp,wp)
            Call zaxpy_(nDIMcf*nDIMcf,zf,O,1,HCF(N,1:nDIMcf,1:nDIMcf),1)
          Else
            zf=cmplx(BC(N,M),0.0_wp,wp)
            Call zaxpy_(nDIMcf*nDIMcf,zf,O,1,HCF(N,1:nDIMcf,1:nDIMcf),1)
            zf=cmplx(BS(N,M),0.0_wp,wp)
            Call zaxpy_(nDIMcf*nDIMcf,zf,W,1,HCF(N,1:nDIMcf,1:nDIMcf),1)
          End If
        End Do
      End Do

      Do k=2,nDIMCF-1,2
        Do ik=2,k
          ! add the ranks as follows:
          ! O2= O2
          ! O4= O2+O4
          ! O6= O2+O4+O6
          ! O8= O2+O4+O6+O8
          ! etc...
          ! compute the cumulative CF matrix
          Do i=1,nDIMcf
            Do j=1,nDIMcf
              HCFS(k,i,j) = HCFS(k,i,j) + HCF(ik,i,j)
            End Do
          End Do
        End Do
      End Do

      Do N=2,nDIMCF-1,2
        Rnrm(N)=dznrm2(nDIMcf*nDIMcf,HCF(N,:,:),1)
        Tnrm   =Tnrm+Rnrm(N)
        Do i=2,N,2
          Snrm(N)=Snrm(N)+Rnrm(i)
        End Do
      End Do
      ! compute the CF spinting of individual weight operators:
      Do k=2,nDIMcf-1,2
        CALL DIAG_C2(  HCF(k,:,:), nDIMcf, info, Wk(k,:), Zk(k,:,:) )
        CALL DIAG_C2( HCFS(k,:,:), nDIMcf, info, Ws(k,:), Zs(k,:,:) )
      End Do
      ! set the initial energies as to the sum of all contributions:
      CALL DIAG_C2( Hinit, nDIMcf, info, Winit, Zinit )

!-----------------------------------------------------------------------
! individual parameter contribution:
      TnrmKQ=0.0_wp
      RnrmKQ(:,:)=0.0_wp
      ListKQ(:)=0.0_wp
      projKQ(:)=0
      rankKQ(:)=0
      ik=0
      Do N=2,nDIMcf-1,2
        Do M=-N,N
          ik=ik+1
          ! generate ITO operators:
          Call Liviu_ITO(nDIMcf,N,ABS(M),O,W,redME)
          ! generate HCF for each parameter rank and projection:
          Call zcopy_(nDIMcf*nDIMcf,[(0.0_wp,0.0_wp)],0,HKQ,1)
          If (M<0) Then
            zf=cmplx(BS(N,ABS(M)),0.0_wp,wp)
            Call zaxpy_(nDIMcf*nDIMcf,zf,W,1,HKQ,1)
          Else If (M>=0) Then
            zf=cmplx(BC(N,ABS(M)),0.0_wp,wp)
            Call zaxpy_(nDIMcf*nDIMcf,zf,O,1,HKQ,1)
          End If
!         find the rank value of each NM operator
!         for further estimate the effect of the corresponding CF parameter
          RnrmKQ(N,M)=dznrm2_(nDIMcf*nDIMcf,HKQ(:,:),1)
          TnrmKQ=TnrmKQ+RnrmKQ(N,M)

          ! indexing lists
          ListKQ(ik)=RnrmKQ(N,M)
          rankKQ(ik)=N
          projKQ(ik)=M
        End Do
      End Do
      ikmax=ik
      ! make an ordered list of CF parameters
      ! following their descending effect:
      Call sort_KQ(ikmax,ListKQ,rankKQ,projKQ,2)
!--------- below we just print the data --------------------------------
      If(LJ=='L') then
          iprog='RASSCF'
      Else
          iprog=' RASSI'
      End If

      If ( iprint >= 4 ) Then
        Write(6,'(A,F20.12)') 'Tnrm=',Tnrm
        Do k=1,nDIMCF-1
          Write(6,'(2(A,i2,A,F20.12,2x))') 'Rnrm(',k,')=',Rnrm(k),
     &                                     'Snrm(',k,')=',Snrm(k)
        End Do
      End If


      ! print the cumulative and individual rank weight:
      Write(6,'(/)')
      Write(6,'(A)') 'CUMULATIVE WEIGHT OF INDIVIDUAL-RANK OPERATORS '//
     &               'ON THE CRYSTAL FIELD SPLITTING:'
      If((nDIMCF-1) >= 2) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O2 :------------------------------------------: ',
     &                             (Snrm(2)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 4) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O2 + O4 :-------------------------------------: ',
     &                             (Snrm(4)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 6) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O2 + O4 + O6 :--------------------------------: ',
     &                             (Snrm(6)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 8) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O2 + O4 + O6 + O8 :---------------------------: ',
     &                             (Snrm(8)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 10) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O2 + O4 + O6 + O8 + O10 :---------------------: ',
     &                             (Snrm(10)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 12) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O2 + O4 + O6 + O8 + O10 + O12 :---------------: ',
     &                             (Snrm(12)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 14) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O2 + O4 + O6 + O8 + O10 + O12 + O14 :---------: ',
     &                             (Snrm(14)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 16) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O2 + O4 + O6 + O8 + O10 + O12 + O14 + O16 :---: ',
     &                             (Snrm(16)/Tnrm)*100.0_wp ,' %.'
      End If

!-----------------------------------------------------------------------
      Write(6,*)
      Write(6,'(A)') 'ENERGY SPLITTING INDUCED BY CUMMULATIVE '//
     & 'INDIVIDUAL-RANK OPERATORS (in cm-1).'
      field(1)='      O2       |'
      field(2)='     O2+O4     |'
      field(3)='   O2+O4+O6    |'
      field(4)='  O2+O4+O6+O8  |'
      field(5)=' O2+O4+...+O10 |'
      field(6)=' O2+O4+...+O12 |'
      field(7)=' O2+O4+...+O14 |'
      field(8)=' O2+O4+...+O16 |'
      nfields=INT((nDIMCF-1)/2)
      Do J=1,nfields,4
        jEnd=MIN(nfields,J+3)
        Write(6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)

        If(MOD(nDIMcf,2).eq.1) Then
          Write(6,'(3A,I2,3x,A,5x,A,4x,A,10(5x,A,6x,A))')
     &     '  ',LJ,' =', (nDIMcf-1)/2,'|',
     &      iprog,'|',('ONLY','|',i=j,jEnd)
        Else
          Write(6,'(3A,I3,A,5x,A,4x,A,10(5x,A,6x,A))')
     &     ' ',LJ,' =',(nDIMcf-1),'/2 |',
     &      iprog,'|',('ONLY','|',i=j,jEnd)
        End If

        Write(6,'(10x,10A)') '|     INITIAL   |',(field(i),i=j,jEnd)
        Write(6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
        Do i=1,nDIMcf
          Write(6,'(1x,A,1x,i2,2x,a,12(f14.8,1x,a))') 'w.f.',i,'|',
     &            (Winit(i)-Winit(1)), '|',
     &           (( Ws(k,i)- Ws(k,1)), '|', k=2*j,2*jEnd,2)
        End Do
        Write(6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
      End Do ! j




      Write(6,'(/)')
      Write(6,'(A)') 'WEIGHT OF INDIVIDUAL-RANK OPERATORS ON THE '//
     &               'CRYSTAL FIELD SPLITTING:'
      If((nDIMCF-1) >= 2) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O2  :-----------------------------------------: ',
     &                             (Rnrm(2)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 4) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O4  :-----------------------------------------: ',
     &                             (Rnrm(4)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 6) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O6  :-----------------------------------------: ',
     &                             (Rnrm(6)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 8) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O8  :-----------------------------------------: ',
     &                             (Rnrm(8)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 10) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O10 :-----------------------------------------: ',
     &                             (Rnrm(10)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 12) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O12 :-----------------------------------------: ',
     &                             (Rnrm(12)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 14) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O14 :-----------------------------------------: ',
     &                             (Rnrm(14)/Tnrm)*100.0_wp ,' %.'
      End If
      If((nDIMCF-1) >= 16) Then
         Write(6,'(2x,A,F10.6,A)')
     &         'O16 :-----------------------------------------: ',
     &                             (Rnrm(16)/Tnrm)*100.0_wp ,' %.'
      End If

!-----------------------------------------------------------------------
      Write(6,'(/)')
      Write(6,'(A)') 'ENERGY SPLITTING INDUCED BY '//
     & 'INDIVIDUAL-RANK OPERATORS (in cm-1).'
      field(1)='      O2       |'
      field(2)='      O4       |'
      field(3)='      O6       |'
      field(4)='      O8       |'
      field(5)='      O10      |'
      field(6)='      O12      |'
      field(7)='      O14      |'
      field(8)='      O16      |'
      nfields=INT((nDIMCF-1)/2)

      Do J=1,nfields,4
        jEnd=MIN(nfields,J+3)
        Write(6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)

        If(MOD(nDIMcf,2).eq.1) Then
          Write(6,'(3A, I2,3x,A,5x,A,4x,A,10(5x,A,6x,A))')
     &     '  ',LJ,' =', (nDIMcf-1)/2,'|',
     &     iprog,'|',('ONLY','|',i=j,jEnd)
        Else
          Write(6,'(3A,I3,a,5x,A,4x,A,10(5x,A,6x,A))')
     &     ' ',LJ,' =',(nDIMcf-1),'/2 |',
     &        iprog,'|',('ONLY','|',i=j,jEnd)
        End If

        Write(6,'(10x,10A)') '|     INITIAL   |',(field(i),i=j,jEnd)
        Write(6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
        Do i=1,nDIMcf
          Write(6,'(1x,A,1x,i2,2x,a,12(f14.8,1x,a))') 'w.f.',i,'|',
     &            (Winit(i)-Winit(1)), '|',
     &           (( Wk(k,i)- Wk(k,1)), '|', k=2*j,2*jEnd,2)
        End Do
        Write(6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
      End Do ! j


!-----------------------------------------------------------------------
      Write(6,'(/)')
      Write(6,'(A)') 'WEIGHT OF INDIVIDUAL CRYSTAL FIELD PARAMETERS'//
     & ' ON THE CRYSTAL FIELD SPLITTING: (in descending order):'
      Write(6,'(A)') 'CFP are given in ITO used in J. Chem. Phys. '//
     &               '137, 064112 (2012).'

      Write(6,'(100A)') ('-',i=1,55),'|'
      Write(6,'(A)') '  k |  q  |         B(k,q)        |'//
     &                          '    Weight (in %)   |'
      Write(6,'(A)') '----|-----|-----------------------|'//
     &                          '--------------------|'
      Do ik=1,ikmax
         ip=projKQ(ik)
         iq=abs(projKQ(ik))
         ir=rankKQ(ik)
         wt=100_wp*ListKQ(ik)/TnrmKQ

         If(projKQ(ik)>=0) Then
            Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),(E22.14,1x,A)'//
     &                'F19.14,1x,A)')
     &               ir,'|',ip,'|',BC(ir,iq),'|',wt,'|'
         Else
            Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),(E22.14,1x,A)'//
     &                'F19.14,1x,A)')
     &               ir,'|',ip,'|',BS(ir,iq),'|',wt,'|'
         End If
      End Do
      Write(6,'(A)') '----|-----|-----------------------|'//
     &                          '--------------------|'




!-----------------------------------------------------------------------
      Call mma_deallocate(Rnrm)
      Call mma_deallocate(Snrm)
      Call mma_deallocate(Wk)
      Call mma_deallocate(Ws)
      Call mma_deallocate(Winit)

      Call mma_deallocate(O)
      Call mma_deallocate(W)
      Call mma_deallocate(HKQ)
      Call mma_deallocate(HCF)
      Call mma_deallocate(HCFS)
      Call mma_deallocate(Zk)
      Call mma_deallocate(Zkq)
      Call mma_deallocate(Zs)
      Call mma_deallocate(Zinit)
      Call mma_deallocate(ListKQ)
      Call mma_deallocate(rankKQ)
      Call mma_deallocate(projKQ)

      Return
      End Subroutine individual_ranks




      Subroutine sort_KQ(N,ARR,rank,proj,iopt)
      Implicit None
      Integer, Parameter          :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)         :: N, iopt
      Integer, intent(inout)      :: rank(N), proj(N)
      Real(kind=8), intent(inout):: ARR(N)
      ! local
      Integer       :: i,j,ir,ip
      Real(kind=8) :: a
!  iopt = 1   => sort in ascending order
!  iopt = 2   => sort in descending order

      If (iopt==1) Then
         Do j=2,N
            a=ARR(j)
            ir=rank(j)
            ip=proj(j)

            Do i = j-1,1,-1
               If ( ARR(i)<=a ) goto 10
               ARR(i+1) =ARR(i)
               rank(i+1)=rank(i)
               proj(i+1)=proj(i)
            End Do
            i=0
 10         Continue
            ARR(i+1) =a
            rank(i+1)=ir
            proj(i+1)=ip
         End Do
      Else If (iopt==2) Then
         Do j=2,N
            a=ARR(j)
            ir=rank(j)
            ip=proj(j)

            Do i = j-1,1,-1
               If ( ARR(i)>=a ) goto 11
               ARR(i+1) =ARR(i)
               rank(i+1)=rank(i)
               proj(i+1)=proj(i)
            End Do
            i=0
 11         Continue
            ARR(i+1) =a
            rank(i+1)=ir
            proj(i+1)=ip
         End Do
      Else
         Write(6,'(A)') 'sort_KQ error:  iopt parameter is wrong.'
         Write(6,    *) 'iopt = ', iopt
         Write(6,'(A)') 'iopt = 1, sort in ascending order'
         Write(6,'(A)') 'iopt = 2, sort in descending order'
         Write(6,'(A)') 'Return, wthout sorting'
         Call xFlush(6)
         Return
      End If

      Return
      End Subroutine sort_KQ
