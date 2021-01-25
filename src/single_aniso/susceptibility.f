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
      Subroutine SUSCEPTIBILITY( NSS, ESO, S_SO, DIPSO, nT, nTempMagn,
     &                           T, tmin, tmax, XTexp, zJ, tinput,
     &                           chiT_theta, doplot, iPrint, mem )

      Implicit None
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer          , intent(in) ::  nss, iprint, nT, nTempMagn, mem
      Real (kind=8)   , intent(in) :: eso(nss)
      Real (kind=8)   , intent(in) :: zJ, tmin, tmax
      Real (kind=8)   , intent(in) ::          T(nT+nTempMagn)
      Real (kind=8)   , intent(in) ::      XTexp(nT+nTempMagn)
      Real (kind=8)   , intent(out):: chit_theta(nT+nTempMagn )
      Complex (kind=8), intent(in) ::   s_so(3,nss,nss)
      Complex (kind=8), intent(in) ::  dipso(3,nss,nss)
      Logical          , intent(in) :: tinput
      Logical          , intent(in) :: doplot
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      this routine calculates the magnetic susceptibility and all related to it values.
c      the units are cgsemu: cm^3*k/mol
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "stdalloc.fh"
      Integer :: ic,jc,it,im,jm,j,i,info,jT,mem_local,RtoB
      Real (kind=8) :: det, xxm, Zst
      Real (kind=8) :: coeff_X, boltz_k, dev
      Logical        :: DBG
      External       :: dev
      Real(kind=8)  :: gtens(3), maxes(3,3)
      Real (kind=8), allocatable ::        chit_tens(:,:,:)
      Real (kind=8), allocatable ::  chit_theta_tens(:,:,:)
      Real (kind=8), allocatable ::           zstat1(:)
      Real (kind=8), allocatable ::             chit(:)
      Real (kind=8), allocatable ::      chi_theta_1(:)
      Real (kind=8), allocatable :: XMM(:,:)
      ! tensors for zJ /= 0
      Real (kind=8), allocatable :: XSM(:,:), XSS(:,:), XZJ(:,:),
     &                               unity(:,:), a_dir(:,:), a_inv(:,:)
      ! main values and axes of XT tensors:
      Real (kind=8), allocatable :: WT(:), ZT(:,:)
      Character(len=50) :: label
c constants used in this subrutine
      RtoB=8
      mem_local=0
      coeff_X=0.125048612_wp*3.0_wp
      boltz_k=0.6950356_wp !boltzmann constant

      DBG=.false.
      If(iPrint.gt.2) DBG=.true.

      Write(6,*)
      Write(6,'(100A)') (('%'),J=1,95)
      Write(6,'(35X,A)') 'CALCULATION OF THE MAGNETIC SUSCEPTIBILITY'
      Write(6,'(100A)') (('%'),J=1,95)

      Write(6,*)
      If (TINPUT) Then
        Write(6,'(5X,A)') 'Temperature dependence of the magnetic '//
     &                    'susceptibility calculated according '
        Write(6,'(5X,A)') 'to experimental values provided in the '//
     &                    'input'
      Else
        Write(6,'(5X,A)') 'Temperature dependence of the '//
     &                    'magnetic susceptibility calculated in'
        Write(6,'(5x,I3,a,f4.1,a,f6.1,a)') nT,' points, equally '//
     &                    'distributed in temperature range',
     &                     Tmin,' ---', Tmax,' K.'
      End If
        Write(6,'(5X,A)') 'The algorithm employed for XT=f(T) in '//
     &                    'this section is based on the zero '//
     &                    'magnetic field limit.'
      If(doplot) Then
        Write(6,'(5X,A)') 'The GNUPLOT script and correponding '//
     &                    'images are generated in $WorkDir'
      End If
      Write(6,*)
!-----------------------------------------------------------------------
      If(dbg) Write(6,*) 'Tmin       =',tmin
      If(dbg) Write(6,*) 'Tmax       =',tmax
      If(dbg) Write(6,*) '  nT       =',nT
      If(dbg) Write(6,*) '  nTempMagn=',nTempMagn
      If(dbg) Write(6,*) 'Temperature:', T(1:nT+nTempMagn)
      If(dbg) Write(6,*) 'ESO:', ESO(1:nss)
      If(dbg) Call prMom('SUSC: input  S_SO(l,i,j):', S_SO,nss)
      If(dbg) Call prMom('SUSC: input DIPSO(l,i,j):',DIPSO,nss)
      If(dbg) gtens=0.0_wp
      If(dbg) maxes=0.0_wp
      If(dbg) Call atens(DIPSO, nss, gtens, maxes, 2)
!-----------------------------------------------------------------------
C   ***********************************************************
C   *     Powder averaged high-T limit of susceptibility      *
C   *                   and g-tensor                          *
C   ***********************************************************
c      ee=0.0_wp
c      ee_2=0.0_wp
c      Do Iss=1,NSS
c        Do Jss=1,NSS
c          Do ic=1,3
c      ee=ee+DBLE(conjg(DipSO(ic,Iss,Jss))*DipSO(ic,Iss,Jss))
c       If(Iss.LE.4.AND.Jss.LE.4) Then
c      ee_2=ee_2+DBLE(conjg(DipSO(ic,Iss,Jss))*DipSO(ic,Iss,Jss))
c       End If
c          End Do ! ic
c        End Do ! Jss
c      End Do ! Iss
c      chiT_high=coeff_X*ee/DBLE(NSS)/3.0_wp
c      SS1=(DBLE(IGSM)-1.0_wp)**2.0_wp/4.0_wp+(DBLE(IGSM)-1.0_wp)/2.0_wp
c      chiT_high_2=coeff_X*ee_2/12.0_wp
c      g_high  =sqrt(chiT_high  /(coeff_X*SS1/3.0_wp))
c      g_high_2=sqrt(chiT_high_2/(coeff_X*SS1/3.0_wp))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Call mma_allocate(chiT,       (nT+nTempMagn),'chiT')
      Call mma_allocate(Zstat1,     (nT+nTempMagn),'Zstat1')
      Call mma_allocate(chi_theta_1,(nT+nTempMagn),'chi_theta_1')
      Call mma_allocate(chiT_tens,  (nT+nTempMagn),3,3,'chiT_tens')
      Call mma_allocate(wt,3,'wt')
      Call mma_allocate(zt,3,3,'zt')
      Call dcopy_(     (nT+nTempMagn), [0.0_wp], 0, chiT,1)
      Call dcopy_(     (nT+nTempMagn), [0.0_wp], 0, chiT_theta,1)
      Call dcopy_(     (nT+nTempMagn), [0.0_wp], 0, Zstat1,1)
      Call dcopy_(     (nT+nTempMagn), [0.0_wp], 0, chi_theta_1,1)
      Call dcopy_( 3*3*(nT+nTempMagn), [0.0_wp], 0, chiT_tens,1)
      mem_local=mem_local+4*(nT+nTempMagn)*RtoB+3*3*(nT+nTempMagn)*RtoB

      If(zJ.eq.0) Then
         If(dbg) Write(6,*) 'SUSC:  zJ = 0'
         Call mma_allocate(XMM,3,3,'XMM')
         mem_local=mem_local+9*RtoB

         If(dbg) Write(6,*) 'SUSC:  memory allocated (local):'
         If(dbg) Write(6,*) 'mem_local=', mem_local
         If(dbg) Write(6,*) 'SUSC:  memory allocated (total):'
         If(dbg) Write(6,*) 'mem_total=', mem+mem_local

         Do iT=1,nT+nTempMagn
           Zst=0.0_wp
           Call dcopy_(3*3,[0.0_wp], 0, XMM,1)
           ! compute XT tensor for this temperature:
           Call chi( DipSO, DipSO, Eso, Nss, T(iT), Zst, XMM)
           If(dbg) Write(6,'(A,9F12.6)') 'XMM:', XMM(1:3,1:3)
           If(dbg) Write(6,'(A,9F12.6)') 'chiT:',
     &              coeff_X*(XMM(1,1)+XMM(2,2)+XMM(3,3))/3.0_wp, Zst

           ! Call dscal_( 3*3, coeff_X, XMM, 1 )
           do i=1,3
             do j=1,3
               chiT_tens(iT,i,j)=coeff_X*XMM(i,j)
             end do
           end do

!           Call dcopy_( 3*3, XMM, 1, chiT_tens(iT,:,:),1)
           ! compute the powder XT for this temperature:
           chiT(iT)       =coeff_X*(XMM(1,1)+XMM(2,2)+XMM(3,3))/3.0_wp

           If( abs(chit(iT)) < 1.0d-20) Then
              chit(iT)=1.d-20
           End If

           chit_theta(iT) =chiT(iT)
           chi_theta_1(iT)=   T(iT)/chiT(iT)
           Zstat1(iT)     =Zst
         End Do
         Call mma_deallocate(XMM)
      Else  ! i.e. when zJ .ne. 0.0_wp
         If(dbg) Write(6,*) 'SUSC:  zJ \= 0'
         ! allocate matrices:
         Call mma_allocate(chiT_theta_tens,(nT+nTempMagn),3,3,'XTT')
         ! initialize:
         Call dcopy_( 3*3*(nT+nTempMagn),[0.0_wp], 0, chiT_theta_tens,1)
         Call mma_allocate(XMM,3,3,'XMM')
         Call mma_allocate(XSM,3,3,'XSM')
         Call mma_allocate(XSS,3,3,'XSS')
         Call mma_allocate(XZJ,3,3,'XZJ')
         Call mma_allocate(A_dir,3,3,'A_dir')
         Call mma_allocate(A_inv,3,3,'A_inv')
         Call mma_allocate(Unity,3,3,'Unity')
         mem_local=mem_local+3*3*3*(nT+nTempMagn)*RtoB+7*3*3*RtoB
         If(dbg) Write(6,*) 'SUSC:  memory allocated (local):'
         If(dbg) Write(6,*) 'mem_local=', mem_local
         If(dbg) Write(6,*) 'SUSC:  memory allocated (total):'
         If(dbg) Write(6,*) 'mem_total=', mem+mem_local

         Do iT=1,nT+nTempMagn
           ! initialize temporary matrices:
           Call dcopy_(3*3,[0.0_wp], 0, XMM,1)
           Call dcopy_(3*3,[0.0_wp], 0, XSM,1)
           Call dcopy_(3*3,[0.0_wp], 0, XSS,1)
           Call dcopy_(3*3,[0.0_wp], 0, XZJ,1)
           Call dcopy_(3*3,[0.0_wp], 0, A_dir,1)
           Call dcopy_(3*3,[0.0_wp], 0, A_inv,1)
           Call dcopy_(3*3,[0.0_wp], 0, Unity,1)
           Unity(1,1)=1.0_wp
           Unity(2,2)=1.0_wp
           Unity(3,3)=1.0_wp
           det=0.0_wp
           Zst=0.0_wp
           ! compute tensors:
           Call chi( DipSO, DipSO, Eso, Nss, T(it), Zst, XMM)
           Call chi(  S_SO, DipSO, Eso, Nss, T(it), Zst, XSM)
           Call chi(  S_SO,  S_SO, Eso, Nss, T(it), Zst, XSS)
          ! form the A_dir matrix:
          ! A_dir(:,:) = 1(:,:)* kB * T(iT) - zJ*XSS(:,:)
          Call daxpy_(3*3,Boltz_k*T(iT),Unity,1,A_dir,1)
          Call daxpy_(3*3,          -zJ,  XSS,1,A_dir,1)
          ! invert it:
          Call REVERSE(A_dir,A_inv,DET)

          Do ic=1,3
            Do jc=1,3
              ! compute the trace of the product XSM*A*XSM
              xxm=0.0_wp
              Do im=1,3
                Do jm=1,3
                  xxm = xxm + XSM(im,ic) * A_inv(im,jm) * XSM(jm,jc)
                End Do
              End Do
              ! add this contribution to the total susceptibility tensor
              XZJ(ic,jc) = XMM(ic,jc) + zJ*xxm
            End Do ! jc
          End Do ! ic

           ! Scale the tensors by coeff_X factor:
           Call dscal_( 3*3, coeff_X, XMM, 1 )
           Call dscal_( 3*3, coeff_X, XZJ, 1 )
!          place the tensors in the corresponding part of the "big" arrays:
           Call dcopy_( 3*3, XMM, 1, chiT_tens(iT,:,:),1)
           Call dcopy_( 3*3, XZJ, 1, chiT_theta_tens(iT,:,:),1)
           ! compute powder:
           chiT(iT)      =(XMM(1,1)+XMM(2,2)+XMM(3,3))/3.0_wp
           chiT_theta(iT)=(XZJ(1,1)+XZJ(2,2)+XZJ(3,3))/3.0_wp

           If( abs(chit(iT)) < 1.0d-20) Then
              chit(iT)=1.d-20
           End If
           If( abs(chiT_theta(iT)) < 1.0d-20) Then
              chiT_theta(iT)=1.d-20
           End If

           chi_theta_1(iT)=   T(iT)/chiT_theta(iT)
           Zstat1(iT)     =Zst
        End Do ! iT
        Call mma_deallocate(XMM)
        Call mma_deallocate(XSM)
        Call mma_deallocate(XSS)
        Call mma_deallocate(XZJ)
        Call mma_deallocate(A_dir)
        Call mma_deallocate(A_inv)
        Call mma_deallocate(Unity)
      End If !zJ

C
C  WRITING SOME OF THE OUTPUT...
C
      Write(6,'(/)')
      Write(6,'(A)') '     |     T      | Statistical |    CHI*T    '//
     &'|    CHI*T    |     CHI     |    1/CHI    |'
      Write(6,'(A)') '     |            |  Sum (Z)    |    (zJ=0)   '//
     &'|             |             |             |'
      Write(6,'(A)') '     |            |             |             '//
     &'|             |             |             |'
      Write(6,'(A)') '-----|----------------------------------------'//
     &'------------------------------------------|'
      Write(6,'(A)') 'Units|   Kelvin   |    ---      |  cm3*K/mol  '//
     &'|  cm3*K/mol  |   cm3/mol   |   mol/cm3   |'
      Write(6,'(A)') '-----|----------------------------------------'//
     &'------------------------------------------|'

      Do iT=1,nT
        jT=iT+nTempMagn
        Write(6,'(A,F11.6,A,E12.5,A,F12.8,A,F12.8,A,E12.5,A,E12.5,A)')
     &   '     |',         T(jT),       ' |',     zstat1(jT),
     &       ' |',      chiT(jT),       ' |', chiT_theta(jT),
     &       ' |',chiT_theta(jT)/T(jT), ' |',chi_theta_1(jT),' |'
      End Do
      Write(6,'(A)') '-----|----------------------------------------'//
     &'------------------------------------------|'



      If(TINPUT) Then
        Write(6,'(/)')
        Write(6,'(5X,A      )') 'STANDARD DEVIATION OF THE CALCULATED'//
     &                          ' MAGNETIC SUSCEPTIBILITY'
        Write(6,'(5X,A,F12.7)') 'FROM EXPERIMENTAL VALUES PROVIDED '//
     &                          'IN THE INPUT FILE IS:',
     &        dev( (nT-nTempMagn), chit_theta( (1+nTempMagn):(nT) ),
     &                                  XTexp( (1+nTempMagn):(nT) )  )

      End If

!-------------------------  PLOTs -------------------------------------!
      WRITE(label,'(A)') "no_field"
      IF ( DoPlot ) THEN
         IF ( tinput ) THEN
            Call plot_XT_with_Exp(label, nT-nTempMagn,
     &                                     T((1+nTempMagn):(nT) ),
     &                            chit_theta((1+nTempMagn):(nT) ),
     &                                 XTexp((1+nTempMagn):(nT) ), zJ )
         ELSE
            Call plot_XT_no_Exp( label, nT-nTempMagn,
     &                                    T((1+nTempMagn):(nT) ),
     &                           chit_theta((1+nTempMagn):(nT) ), zJ )
         END IF
      END IF
!------------------------- END PLOTs -------------------------------------!



c print out the main VAN VLECK SUSCEPTIBILITY TENSOR, its main values and main axes:
      If(zJ.eq.0.0_wp) Then

        Write(6,'(/)')
        Write(6,'(111A)') ('-',i=1,110),'|'
        Write(6,'(31X,A,22x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR'//
     &                           ' FOR zJ = 0,  in cm3*K/mol','|'
        Write(6,'(111A)') ('-',i=1,110),'|'
        Write(6,'(A   )') '     T(K)   |   |          '//
     &                    'Susceptibility Tensor'//
     &                    '      |    Main Values  |'//
     &                    '               Main Axes             |'
        Do iT=1,nT
          jT=iT+nTempMagn
          info=0
          Call dcopy_(  3, [0.0_wp], 0, wt,1)
          Call dcopy_(3*3, [0.0_wp], 0, zt,1)
          Call DIAG_R2( chiT_tens(jT,:,:) ,3,info,wt,zt)
          Write(6,'(A)') '------------|---'//
     &                   '|------- x --------- y --------- z ---'//
     &                   '|-----------------'//
     &                   '|------ x --------- y --------- z ----|'
          Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &             '            | x |',(chiT_tens(jT,1,j),j=1,3),
     &             ' |  X:',wt(1),'|',(zt(j,1),j=1,3),'|'
          Write(6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),
     &             '| y |',            (chiT_tens(jT,2,j),j=1,3),
     &             ' |  Y:',wt(2),'|',(zt(j,2),j=1,3),'|'
          Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &             '            | z |',(chiT_tens(jT,3,j),j=1,3),
     &             ' |  Z:',wt(3),'|',(zt(j,3),j=1,3),'|'
        End Do
        Write(6,'(111A)') ('-',i=1,110),'|'

      Else ! zJ .ne. 0.0_wp

        Write(6,'(/)')
        Write(6,'(111A)') ('-',i=1,110),'|'
        Write(6,'(31X,A,F7.4,15x,A)')
     &                'VAN VLECK SUSCEPTIBILITY TENSOR FOR zJ = ',zJ,
     &                ',  in cm3*K/mol','|'
        Write(6,'(111A)') ('-',i=1,110),'|'
        Write(6,'(A)')'     T(K)   |   |          Susceptibility '//
     &                'Tensor      |    Main Values  '//
     &                '|               Main Axes             |'
        Do iT=1,nT
          jT=iT+nTempMagn
          info=0
          Call dcopy_(  3, [0.0_wp], 0, wt,1)
          Call dcopy_(3*3, [0.0_wp], 0, zt,1)
          Call DIAG_R2( chiT_theta_tens(jT,:,:) ,3,info,wt,zt)
          Write(6,'(A)') '------------|---'//
     &                   '|------- x --------- y --------- z ---'//
     &                   '|-----------------'//
     &                   '|------ x --------- y --------- z ----|'
          Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &          '            | x |',(chiT_theta_tens(jT,1,j),j=1,3),
     &          ' |  X:',wt(1),'|',(zt(j,1),j=1,3),'|'
          Write(6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),
     &                      '| y |',(chiT_theta_tens(jT,2,j),j=1,3),
     &          ' |  Y:',wt(2),'|',(zt(j,2),j=1,3),'|'
          Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &          '            | z |',(chiT_theta_tens(jT,3,j),j=1,3),
     &          ' |  Z:',wt(3),'|',(zt(j,3),j=1,3),'|'
        End Do
        Write(6,'(111A)') ('-',i=1,110),'|'
         Call mma_deallocate(chiT_theta_tens)
      End If ! zJ

c      Write(6,'(/)')
c      Write(6,'(10A)') (('------------'), K=1,10)
c      Write(6,'(5X,A)') 'MAGNETIC SUSCEPTIBILITY IN THE DIRECTION '//
c     &'OF THE MAIN MAGNETIC AXES'
c      Write(6,'(10A)') (('------------'), K=1,10)
c      Write(6,*)
c      Write(6,'(7X,A)') 'T(K)         gx          gy          gz'
c      Write(6,*)
c      Do iT=1,nT
c      Write(6,'(4X,F6.1,6X,3(F8.4,4X))') T(iT),(ChiT_main(iT,ic),ic=1,3)
c      End Do
C  saving some information for tests:

      Call Add_Info('Temperature',T         ,nT+nTempMagn,5)
      Call Add_Info('CHIT'       ,chiT      ,nT+nTempMagn,5)
      Call Add_Info('CHIT_THETA' ,chiT_theta,nT+nTempMagn,5)

      Call mma_deallocate(Zstat1)
      Call mma_deallocate(chiT_tens)
      Call mma_deallocate(chiT)
      Call mma_deallocate(chi_theta_1)
      Call mma_deallocate(wt)
      Call mma_deallocate(zt)

      Return
      End

C Calculation of susceptibility in the direction of main magnetic axes:

c         Do i=1,3
c           Do j=1,3
c             Do k=1,3
c               ChiT_main(iT,i)=ChiT_main(iT,i)+ ZMAGN(k,i)*ZMAGN(j,i)
c       &                                       *ChiT_tens(iT,k,j)
c             End Do
c           End Do
c         End Do
