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
       Subroutine susceptibility_pa( exch, nLoc, nCenter, nneq,
     &                            neqv, neq, nss, nexch, nTempMagn,
     &                            nT, Tmin, Tmax, XTexp,
     &                            eso, dipso, s_so, w, dipexch,
     &                            s_exch, T, R_LG, zJ, tinput,
     &                            XLM, ZLM, XRM, ZRM, iopt,
     &                            chiT_theta, doplot, mem )

c       chi*t ----------- the units are cgsemu: [ cm^3*k/mol ]
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)

      Integer, intent(in) :: nLoc, nCenter, nTempMagn, nT, mem
      Integer, intent(in) :: exch, nneq, neqv, iopt
      Integer, intent(in) :: neq(nneq), nss(nneq), nexch(nneq)

      Real(kind=wp), intent(in) :: T(nT+nTempMagn)
      Real(kind=wp), intent(in) :: W(exch)
      Real(kind=wp), intent(in) :: eso(nneq,nLoc)
      Real(kind=wp), intent(in) :: zJ
      Real(kind=wp), intent(in) :: Tmin, Tmax
      Real(kind=wp), intent(in) :: XTexp(nT+nTempMagn)
      Real(kind=wp), intent(in) :: R_LG(nneq,neqv,3,3)
      Real(kind=wp), intent(out):: chit_theta(nT+nTempMagn)
c contributions from local excited states, computed in the XT section:
      Real(kind=wp), intent(out):: XLM( nCenter,nTempMagn,3,3)
      Real(kind=wp), intent(out):: ZLM( nCenter,nTempMagn)
      Real(kind=wp), intent(out):: XRM( nCenter,nTempMagn,3,3)
      Real(kind=wp), intent(out):: ZRM( nCenter,nTempMagn)
      Logical, intent(in)       :: tinput, doplot
      ! BIG matrices:
      Complex(kind=wp), intent(in) :: dipexch(3,exch,exch)
      Complex(kind=wp), intent(in) ::  s_exch(3,exch,exch)
      Complex(kind=wp), intent(in) :: dipso(nneq,3,nLoc,nLoc)
      Complex(kind=wp), intent(in) ::  s_so(nneq,3,nLoc,nLoc)
#include "stdalloc.fh"

c local variables
      Real(kind=wp), allocatable ::     chit_tens_l(:,:,:)
!                                       chit_tens_l( nneq,3,3)
      Real(kind=wp), allocatable :: smu_chit_tens_l(:,:,:)
!                                   smu_chit_tens_l( nneq,3,3)
      Real(kind=wp), allocatable ::  ss_chit_tens_l(:,:,:)
!                                    ss_chit_tens_l( nneq,3,3)
      Real(kind=wp), allocatable ::     chit_tens_lr(:,:,:)
!                                       chit_tens_lr(nneq,3,3)
      Real(kind=wp), allocatable :: smu_chit_tens_lr(:,:,:)
!                                   smu_chit_tens_lr(nneq,3,3)
      Real(kind=wp), allocatable ::  ss_chit_tens_lr(:,:,:)
!                                    ss_chit_tens_lr(nneq,3,3)
      Real(kind=wp), allocatable ::     chit_tens_ex(:,:)
!                                       chit_tens_ex(3,3)
      Real(kind=wp), allocatable :: smu_chit_tens_ex(:,:)
!                                   smu_chit_tens_ex(3,3)
      Real(kind=wp), allocatable ::  ss_chit_tens_ex(:,:)
!                                     ss_chit_tens_ex(3,3)
      Real(kind=wp), allocatable ::     chit_tens_tot(:,:,:)
!                                       chit_tens_tot(nT+nTempMagn,3,3)
      Real(kind=wp), allocatable :: smu_chit_tens_tot(:,:)
!                                   smu_chit_tens_tot(3,3)
      Real(kind=wp), allocatable ::  ss_chit_tens_tot(:,:)
!                                    ss_chit_tens_tot(3,3)
      Real(kind=wp), allocatable ::   chit_theta_tens(:,:,:)
!                                     chit_theta_tens(nT+nTempMagn,3,3)
      Real(kind=wp), allocatable :: zstat_l(:)       !zstat_l( nneq)
      Real(kind=wp), allocatable :: zstat_lr(:)      !zstat_lr(nneq)
      Real(kind=wp), allocatable :: zstat_tot(:)
!                                   zstat_tot(nT+nTempMagn)
      Real(kind=wp), allocatable :: chit(:)          !chit(nT+nTempMagn)
      Real(kind=wp), allocatable :: chi_theta_1(:)
!                                   chi_theta_1(nT+nTempMagn)
      Real(kind=wp), allocatable :: XL(:,:,:)        !XL(nCenter,3,3)
      Real(kind=wp), allocatable :: ZL(:)            !ZL(nCenter)
      Real(kind=wp), allocatable :: XR(:,:,:)        !XR(nCenter,3,3)
      Real(kind=wp), allocatable :: ZR(:)            !ZR(nCenter)
      Real(kind=wp), allocatable :: SMUR(:,:,:)      !SMUR(nCenter,3,3)
      Real(kind=wp), allocatable :: SMUL(:,:,:)      !SMUL(nCenter,3,3)
      Real(kind=wp), allocatable :: SSR(:,:,:)       !SSR( nCenter,3,3)
      Real(kind=wp), allocatable :: SSL(:,:,:)       !SSL( nCenter,3,3)
      Real(kind=wp), allocatable :: wt(:), zt(:,:)   !wt(3),zt(3,3)
      Real(kind=wp), allocatable :: A_dir(:,:)               !A_dir(3,3)
      Real(kind=wp), allocatable :: A_inv(:,:)               !A_inv(3,3)
      Real(kind=wp), allocatable :: unity(:,:)               !unity(3,3)
      Real(kind=wp) :: xxm
      Real(kind=wp) :: zstat_ex
      Real(kind=wp) :: boltz_k,coeff_chi
      Real(kind=wp) :: det
      Real(kind=wp) :: dev, Fa, Fb, Fc, Fd, Fe, Ff
      external dev
      Integer       :: i,iT,jT,ic,jc
      Integer       :: j,n1,n2,im,jm
      Integer       :: isite,info,mem_local,RtoB
      Logical       :: dbg
      Character(len=50) :: label
      Real(wp), external :: dnrm2_

      Call qEnter('PA_suscept')

      mem_local=0
      dbg=.false.
      RtoB=8
!     = n_a*mu_bohr^2/(k_boltz) in cm^3*k/mol
      coeff_chi=0.1250486120_wp*3.0_wp
      boltz_k=0.69503560_wp                    !   in cm^-1*k-1
!-----------------------------------------------------------------------
      If(dbg) Then
         Write(6,*) 'Verification of input data on entrance to PA-SUSC:'
         Write(6,*) 'exch:         ', exch
         Write(6,*) 'nLoc:         ', nLoc
         Write(6,*) 'nCenter:      ', nCenter
         Write(6,*) 'nneq:         ', nneq
         Write(6,*) 'neqv:         ', neqv
         Write(6,*) 'neq():        ', (neq(i),i=1,nneq)
         Write(6,*) 'nss():        ', (nss(i),i=1,nneq)
         Write(6,*) 'nexch():      ', (nexch(i),i=1,nneq)
         Write(6,*) 'nT:           ', nT
         Write(6,*) 'iopt:         ', iopt
         Write(6,*) 'mem:          ', mem
         Write(6,*) 'nTempMagn:    ', nTempMagn
         Write(6,*) 'Tmin:         ', Tmin
         Write(6,*) 'Tmax:         ', Tmax
         Write(6,*) 'XTexp():      ', (XTexp(i),i=1,nT+nTempMagn)
         Write(6,*) 'T():          ', (T(i),i=1,nT+nTempMagn)
         Write(6,*) 'chit_theta(): ', (chit_theta(i),i=1,nT+nTempMagn)
         Write(6,*) 'W()           ', (W(i),i=1,exch)
         Write(6,*) 'zJ:           ', zJ
         Write(6,*) 'tinput:       ', tinput
         Write(6,*) 'doplot:       ', doplot
!         Do i=1,nneq
!            Write(6,*) 'eso()         ', (eso(i,j),j=1,nss(i))
!            Call prMom('SUSC: input  s_so(i,:,:,:):', s_so(i,:,:,:),
!     &               nexch(i))
!            Call prMom('SUSC: input dipso(i,:,:,:):',dipso(i,:,:,:),
!     &               nexch(i))
!            gtens=0.0_wp
!            maxes=0.0_wp
!            Call atens(s_so(i,:,:,:), nexch(i), gtens, maxes, 2)
!            gtens=0.0_wp
!            maxes=0.0_wp
!            Call atens(dipso(i,:,:,:), nexch(i), gtens, maxes, 2)
!         End Do
!         Call prMom('SUSC: input  S_EXCH(l,i,j):', s_exch,exch)
!         Call prMom('SUSC: input DIPEXCH(l,i,j):',dipexch,exch)
!         gtens=0.0_wp
!         maxes=0.0_wp
!         Call atens(s_exch, exch, gtens, maxes, 2)
!         gtens=0.0_wp
!         maxes=0.0_wp
!         Call atens(dipexch, exch, gtens, maxes, 2)
!      change to pseudospin:
!        Call zcopy_(  exch*exch,[(0.0_wp,0.0_wp)],0,Z,1)
!        Call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,dipexch2,1)
!        Call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0, s_exch2,1)
!        Call zcopy_(3*exch*exch,dipexch,1,dipexch2,1)
!        Call zcopy_(3*exch*exch, s_exch,1, s_exch2,1)
!        Call pseudospin(dipexch2,exch,Z,3,1)

!        Call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,dipexch,1)
!        Call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,s_exch,1)
!        Call UTMU( exch, exch, Z, dipexch2, dipexch )
!        Call UTMU( exch, exch, Z, s_exch2, s_exch )

!        Call prMom('SUSC: input  S_EXCH2(l,i,j):', s_exch,exch)
!        Call prMom('SUSC: input DIPEXCH2(l,i,j):',dipexch,exch)
      End If ! dbg
!-----------------------------------------------------------------------
      Write(6,*)
      Write(6,'(100A)') (('%'),J=1,95)
      Write(6,'(35X,A)') 'CALCULATION OF THE MAGNETIC SUSCEPTIBILITY'
      Write(6,'(100A)') (('%'),J=1,95)
      Write(6,*)

      If(tinput) Then
         Write(6,'(2x,a)') 'Temperature dependence of the magnetic '//
     &                     'susceptibility and'
         Write(6,'(2x,a)') 'high-field magnetization will be '//
     &                     'calculated according to '
         Write(6,'(2x,a)') 'experimental values provided in the '//
     &                     'input.'
      Else
         Write(6,'(2x,a,i3,a)') 'Temperature dependence of the '//
     &                          'magnetic susceptibility will be '//
     &                          'calculated in',nT,' points, '
         Write(6,'(2x,a,f4.1,a,f6.1,a)') 'equally distributed in '//
     &                          'temperature range ',tmin,' ---',
     &                           tmax,' K.'
      End If
      Call XFlush(6)
!-----------------------------------------------------------------------
      Call mma_allocate(chit_tens_tot,(nT+nTempMagn),3,3,'XT_tens_tot')
      Call mma_allocate(chit_theta_tens,(nT+nTempMagn),3,3,'XT_theta_t')
      Call mma_allocate(chit_tens_l ,nneq,3,3,'XT_tens_l')
      Call mma_allocate(chit_tens_lr,nneq,3,3,'XT_tens_lr')
      Call mma_allocate(chit_tens_ex,3,3,'XT_tens_ex')
      Call mma_allocate(zstat_l,nneq,'Zstat_l')
      Call mma_allocate(zstat_lr,nneq,'Zstat_lr')
      mem_local=mem_local+(2+2*3*3)*nneq*RtoB+9*RtoB

      Call mma_allocate(ZR,nCenter,'ZR')
      Call mma_allocate(ZL,nCenter,'ZL')
      Call mma_allocate(XL,nCenter,3,3,'XL')
      Call mma_allocate(XR,nCenter,3,3,'XR')
      mem_local=mem_local+(5+2*3*3)*nCenter*RtoB

      Call mma_allocate(zstat_tot,(nT+nTempMagn),'Zstat_tot')
      Call mma_allocate(chiT,(nT+nTempMagn),'chiT')
      Call mma_allocate(chi_theta_1,(nT+nTempMagn),'chi_theta_1')
      mem_local=mem_local+(3+2*3*3)*(nT+nTempMagn)*RtoB


      Call dcopy_(3*3*(nT+nTempMagn),[0.0_wp],0,chit_tens_tot,1)
      Call dcopy_(3*3*(nT+nTempMagn),[0.0_wp],0,chit_theta_tens,1)
      Call dcopy_(    (nT+nTempMagn),[0.0_wp],0,zstat_tot,1)

      If (zJ == 0.0_wp) Then
         If(dbg) Write(6,*) 'SUSC:  memory allocated (local):'
         If(dbg) Write(6,*) 'mem_local=', mem_local
         If(dbg) Write(6,*) 'SUSC:  memory allocated (total):'
         If(dbg) Write(6,*) 'mem_total=', mem+mem_local

         Do iT=1,nT+nTempMagn
            ! initialize temporary variables:
            Call dcopy_(3*3,[0.0_wp],0,chit_tens_ex,1)
            Call dcopy_(nneq*3*3,[0.0_wp],0,chit_tens_l,1)
            Call dcopy_(nneq*3*3,[0.0_wp],0,chit_tens_lr,1)
            Call dcopy_(nneq,[0.0_wp],0,zstat_l,1)
            Call dcopy_(nneq,[0.0_wp],0,zstat_lr,1)
            Call dcopy_(nCenter,[0.0_wp],0,ZR,1)
            Call dcopy_(nCenter,[0.0_wp],0,ZL,1)
            Call dcopy_(3*3*nCenter,[0.0_wp],0,XL,1)
            Call dcopy_(3*3*nCenter,[0.0_wp],0,XR,1)
            zstat_ex=0.0_wp
c------------------------------------------------------------------------------------
cc  local susceptibility= total susceptibility coming from individual magnetic centers
            Do i=1, nneq
              If(dbg) Write(6,'(A,2I5)') 'nss(i)=',nss(i)
              If(dbg) Write(6,'(A,2I5)') 'nexch(i)=',nexch(i)
              If(dbg) Write(6,'(A,9F10.6)') 'eso(i,:)=',eso(i,1:nss(i))
              If(dbg) Write(6,'(A,9F10.6)') 'W(:)    =',W(1:exch)
              If(dbg) Write(6,'(A,9F10.6)') 'T(iT)=',T(iT)
              Call chi( dipso( i, 1:3, 1:nss(i), 1:nss(i) ),
     &                  dipso( i, 1:3, 1:nss(i), 1:nss(i) ),
     &                    eso( i, 1:nss(i)), nss(i), T(it),
     &                zstat_l( i), chit_tens_l(i,1:3,1:3) )

              Call chi( dipso( i, 1:3, 1:nexch(i), 1:nexch(i) ),
     &                  dipso( i, 1:3, 1:nexch(i), 1:nexch(i) ),
     &                    eso( i, 1:nexch(i)), nexch(i), T(iT),
     &               zstat_lr( i), chit_tens_lr(i,1:3,1:3) )

            End Do ! i (nneq)
            Call chi( dipexch, dipexch, W, exch, T(iT), zstat_ex,
     &                chit_tens_ex )

            Fa=0.0_wp; Fb=0.0_wp; Fc=0.0_wp; Fd=0.0_wp; Fe=0.0_wp;
            Fa=dnrm2_(9*nneq,chit_tens_l ,1)
            Fb=dnrm2_(9*nneq,chit_tens_lr,1)
            Fc=dnrm2_(9     ,chit_tens_ex,1)
            Fd=dnrm2_(nneq  ,zstat_l     ,1)
            Fe=dnrm2_(nneq  ,zstat_lr    ,1)
            Call Add_Info('XT:  chit_tens_l'   ,Fa,1,6)
            Call Add_Info('XT:  chit_tens_lr'  ,Fb,1,6)
            Call Add_Info('XT:  chit_tens_exch',Fc,1,6)
            Call Add_Info('XT:  zstat_exch'    ,zstat_ex,1,6)
            Call Add_Info('XT:  zstat_l'       ,Fd,1,6)
            Call Add_Info('XT:  zstat_lr'      ,Fe,1,6)
c expand the basis and rotate local tensors to the general
c coordinate system:
            isite=0
            Do i=1,nneq
               Do j=1,neq(i)
               isite=isite+1
               ZL(isite)=zstat_l(i)
               ZR(isite)=zstat_lr(i)
                  Do ic=1,3
                     Do jc=1,3
                        Do n1=1,3
                           Do n2=1,3
                           XL(isite,ic,jc)=XL(isite,ic,jc)+
     &                                     r_lg(i,j,ic,n1)*
     &                                     r_lg(i,j,jc,n2)*
     &                                     chit_tens_l(i,n1,n2)
                           XR(isite,ic,jc)=XR(isite,ic,jc)+
     &                                     r_lg(i,j,ic,n1)*
     &                                     r_lg(i,j,jc,n2)*
     &                                     chit_tens_lr(i,n1,n2)
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do

            Fa=0.0_wp; Fb=0.0_wp;
            Fa=dnrm2_(9*nCenter,XL,1)
            Fb=dnrm2_(9*nCenter,XR,1)
            Call Add_Info('XT:  XL',Fa,1,6)
            Call Add_Info('XT:  XR',Fb,1,6)
c save some data:
            If(it.le.nTempMagn) Then
               Call dscal_( 3*3*nCenter, coeff_chi, XR, 1 )
               Call dscal_( 3*3*nCenter, coeff_chi, XL, 1 )
               Call dcopy_( 3*3*nCenter, XR, 1, XRM(:,iT,:,:),1)
               Call dcopy_( 3*3*nCenter, XL, 1, XLM(:,iT,:,:),1)
               Call dcopy_(     nCenter, ZR, 1, ZRM(:,iT),1)
               Call dcopy_(     nCenter, ZL, 1, ZLM(:,iT),1)
            End If
            Call chi_sum( nCenter, chit_tens_ex, zstat_ex,
     &                    XL, ZL, XR, ZR, iopt,
     &                    chit_tens_tot(it,1:3,1:3), zstat_tot(it) )

            Fa=0.0_wp; Fb=0.0_wp;
            Fa=dnrm2_(nCenter,ZL,1)
            Fb=dnrm2_(nCenter,ZR,1)
            Call Add_Info('XT:  ZL',ZL,1,6)
            Call Add_Info('XT:  ZR',ZR,1,6)
            Call Add_Info('XT: ZEx',zstat_ex,1,6)

            chit(it)=coeff_chi * ( chit_tens_tot(iT,1,1)
     &                            +chit_tens_tot(iT,2,2)
     &                            +chit_tens_tot(iT,3,3) )/3.0_wp
            chit_theta(iT)=chit(iT)

            If( abs(chit(iT)) < 1.0e-20_wp) Then
               chit(iT)=1.0e-20_wp
               chit_theta(iT)=1.0e-20_wp
            End If
            chi_theta_1(iT)=T(iT)/chit(iT)
            ! add some verification data:
            Fa=0.0_wp
            Fa=dnrm2_(9,chit_theta_tens(iT,1:3,1:3),1)
            Call Add_Info('XT: chit_theta_tens',Fa,1,6)
         End Do ! iT
         Fb=0.0_wp
         Fb=dnrm2_(nT+nTempMagn,T,1)
         Call Add_Info('XT: T',T,nT+nTempMagn,6)

      Else ! i.e. when (zJ .ne. 0)

         ! allocate memory for temporary arrays:
         Call mma_allocate(smu_chit_tens_l,nneq,3,3,'smu_X_tens_l')
         Call mma_allocate(smu_chit_tens_lr,nneq,3,3,'smu_X_tens_lr')
         Call mma_allocate(smu_chit_tens_ex,3,3,'smu_chit_X_ex')
         Call mma_allocate(smu_chit_tens_tot,3,3,'SM_XT')
         Call mma_allocate(ss_chit_tens_l,nneq,3,3,'ss_chit_tens_l')
         Call mma_allocate(ss_chit_tens_lr,nneq,3,3,'ss_chit_tens_lr')
         Call mma_allocate(ss_chit_tens_ex,3,3,'ss_chit_tens_ex')
         Call mma_allocate(ss_chit_tens_tot,3,3,'SS_XT')
         Call mma_allocate(SMUR,nCenter,3,3,'SMUR')
         Call mma_allocate(SMUL,nCenter,3,3,'SMUL')
         Call mma_allocate(SSR,nCenter,3,3,'SSR')
         Call mma_allocate(SSL,nCenter,3,3,'SSL')
         Call mma_allocate(A_dir,3,3,'A_dir')
         Call mma_allocate(A_inv,3,3,'A_inv')
         Call mma_allocate(Unity,3,3,'Unity')
         mem_local=mem_local+4*3*3*(nCenter+nneq)*RtoB
         mem_local=mem_local+7*3*3*RtoB

         If(dbg) Write(6,*) 'SUSC:  memory allocated (local):'
         If(dbg) Write(6,*) 'mem_local=', mem_local
         If(dbg) Write(6,*) 'SUSC:  memory allocated (total):'
         If(dbg) Write(6,*) 'mem_total=', mem+mem_local

         Do iT=1,nT+nTempMagn
            ! initialization:
            Call dcopy_(     3*3   ,[0.0_wp],0,chit_tens_ex,1)
            Call dcopy_(     3*3   ,[0.0_wp],0,smu_chit_tens_ex,1)
            Call dcopy_(     3*3   ,[0.0_wp],0,ss_chit_tens_ex,1)
            Call dcopy_(nneq*3*3   ,[0.0_wp],0,chit_tens_l,1)
            Call dcopy_(nneq*3*3   ,[0.0_wp],0,smu_chit_tens_l,1)
            Call dcopy_(nneq*3*3   ,[0.0_wp],0,ss_chit_tens_l,1)
            Call dcopy_(nneq*3*3   ,[0.0_wp],0,chit_tens_lr,1)
            Call dcopy_(nneq*3*3   ,[0.0_wp],0,smu_chit_tens_lr,1)
            Call dcopy_(nneq*3*3   ,[0.0_wp],0,ss_chit_tens_lr,1)
            Call dcopy_(     3*3   ,[0.0_wp],0,smu_chit_tens_tot,1)
            Call dcopy_(     3*3   ,[0.0_wp],0,ss_chit_tens_tot,1)
            Call dcopy_(nneq       ,[0.0_wp],0,zstat_l,1)
            Call dcopy_(nneq       ,[0.0_wp],0,zstat_lr,1)
            Call dcopy_(nCenter    ,[0.0_wp],0,ZR,1)
            Call dcopy_(nCenter    ,[0.0_wp],0,ZL,1)
            Call dcopy_(nCenter*3*3,[0.0_wp],0,XL,1)
            Call dcopy_(nCenter*3*3,[0.0_wp],0,XR,1)
            Call dcopy_(nCenter*3*3,[0.0_wp],0,SMUR,1)
            Call dcopy_(nCenter*3*3,[0.0_wp],0,SMUL,1)
            Call dcopy_(nCenter*3*3,[0.0_wp],0,SSR,1)
            Call dcopy_(nCenter*3*3,[0.0_wp],0,SSL,1)
            Call dcopy_(        3*3,[0.0_wp],0,A_dir,1)
            Call dcopy_(        3*3,[0.0_wp],0,A_inv,1)
            Call dcopy_(        3*3,[0.0_wp],0,Unity,1)
            zstat_ex=0.0_wp
            det=0.0_wp
            Do ic=1,3
              unity(ic,ic)=1.0_wp
            End Do

            ! compute local tensors  L, and LR:
            Do i=1, nneq
               Call chi( dipso( i, 1:3, 1:nss(i), 1:nss(i) ),
     &                   dipso( i, 1:3, 1:nss(i), 1:nss(i) ),
     &                     eso( i, 1:nss(i)), nss(i), T(iT),
     &                 zstat_l( i), chit_tens_l(i,1:3,1:3) )

               Call chi(  s_so( i, 1:3, 1:nss(i), 1:nss(i) ),
     &                   dipso( i, 1:3, 1:nss(i), 1:nss(i) ),
     &                     eso( i, 1:nss(i)), nss(i), T(iT),
     &                 zstat_l( i), smu_chit_tens_l(i,1:3,1:3) )

               Call chi(  s_so( i, 1:3, 1:nss(i), 1:nss(i) ),
     &                    s_so( i, 1:3, 1:nss(i), 1:nss(i) ),
     &                     eso( i, 1:nss(i)), nss(i), T(iT),
     &                 zstat_l( i), ss_chit_tens_l(i,1:3,1:3) )

               Call chi( dipso( i, 1:3, 1:nexch(i), 1:nexch(i) ),
     &                   dipso( i, 1:3, 1:nexch(i), 1:nexch(i) ),
     &                     eso( i, 1:nexch(i)), nexch(i), T(iT),
     &                zstat_lr( i), chit_tens_lr(i,1:3,1:3) )

               Call chi(  s_so( i, 1:3, 1:nexch(i), 1:nexch(i) ),
     &                   dipso( i, 1:3, 1:nexch(i), 1:nexch(i) ),
     &                     eso( i, 1:nexch(i)), nexch(i), T(iT),
     &                zstat_lr( i), smu_chit_tens_lr(i,1:3,1:3) )

               Call chi(  s_so( i, 1:3, 1:nexch(i), 1:nexch(i) ),
     &                    s_so( i, 1:3, 1:nexch(i), 1:nexch(i) ),
     &                     eso( i, 1:nexch(i)), nexch(i), T(iT),
     &                zstat_lr( i), ss_chit_tens_lr(i,1:3,1:3) )
            End Do ! i (nneq)

            ! compute exchange tensors:
            Call chi( dipexch, dipexch, W, exch, T(it), zstat_ex,
     &                     chit_tens_ex )
            Call chi(  s_exch, dipexch, W, exch, T(iT), zstat_ex,
     &                 smu_chit_tens_ex )
            Call chi(  s_exch,  s_exch, W, exch, T(iT), zstat_ex,
     &                  ss_chit_tens_ex )
c expand the basis and rotate local tensors to the general
c coordinate system:
            isite=0
            Do i=1,nneq
               Do j=1,neq(i)
               isite=isite+1
               ZL(isite)=zstat_l(i)
               ZR(isite)=zstat_lr(i)
!              use R_lg matrices, which have arbitrary determinant:  +1  or -1;
!              reason:  X_ab is a bi-dimensional tensor.
!                       We need to rotate twice ==> the sign of R does not
!                       matter (+1 * +1) = (-1 * -1)
!               >> R_rot matrices have determinant strict +1, and are used to
!                       rotate vectors
                  Do ic=1,3
                     Do jc=1,3
                        Do n1=1,3
                           Do n2=1,3
                           XR(isite,ic,jc)=XR(isite,ic,jc)+
     &                                     r_lg(i,j,ic,n1)*
     &                                     r_lg(i,j,jc,n2)*
     &                                     chit_tens_lr(i,n1,n2)
                           XL(isite,ic,jc)=XL(isite,ic,jc)+
     &                                     r_lg(i,j,ic,n1)*
     &                                     r_lg(i,j,jc,n2)*
     &                                     chit_tens_l(i,n1,n2)

                         SMUL(isite,ic,jc)=SMUL(isite,ic,jc)+
     &                                     r_lg(i,j,ic,n1)*
     &                                     r_lg(i,j,jc,n2)*
     &                                     smu_chit_tens_l(i,n1,n2)
                         SMUR(isite,ic,jc)=SMUR(isite,ic,jc)+
     &                                     r_lg(i,j,ic,n1)*
     &                                     r_lg(i,j,jc,n2)*
     &                                     smu_chit_tens_lr(i,n1,n2)

                         SSL(isite,ic,jc)=SSL(isite,ic,jc)+
     &                                    r_lg(i,j,ic,n1)*
     &                                    r_lg(i,j,jc,n2)*
     &                                    ss_chit_tens_l(i,n1,n2)
                         SSR(isite,ic,jc)=SSR(isite,ic,jc)+
     &                                    r_lg(i,j,ic,n1)*
     &                                    r_lg(i,j,jc,n2)*
     &                                    ss_chit_tens_lr(i,n1,n2)
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do

            ! add verification data:
            Fa=0.0_wp; Fb=0.0_wp; Fc=0.0_wp; Fd=0.0_wp; Fe=0.0_wp;
            Ff=0.0_wp;
            Fa=dnrm2_(9*nCenter,XR  ,1)
            Fb=dnrm2_(9*nCenter,XL  ,1)
            Fc=dnrm2_(9*nCenter,SMUL,1)
            Fd=dnrm2_(9*nCenter,SMUR,1)
            Fe=dnrm2_(9*nCenter,SSL,1)
            Ff=dnrm2_(9*nCenter,SSR,1)
            Call Add_Info('XT:    XR',Fa,1,6)
            Call Add_Info('XT:    XL',Fb,1,6)
            Call Add_Info('XT:  SMUL',Fc,1,6)
            Call Add_Info('XT:  SMUR',Fd,1,6)
            Call Add_Info('XT:   SSL',Fe,1,6)
            Call Add_Info('XT:   SSR',Ff,1,6)

c save some data:
            If(iT.le.nTempMagn) Then
               Call dscal_( 3*3*nCenter, coeff_chi, XR, 1 )
               Call dscal_( 3*3*nCenter, coeff_chi, XL, 1 )
               Call dcopy_( 3*3*nCenter, XR, 1, XRM(:,iT,:,:),1)
               Call dcopy_( 3*3*nCenter, XL, 1, XLM(:,iT,:,:),1)
               Call dcopy_(     nCenter, ZR, 1, ZRM(:,iT),1)
               Call dcopy_(     nCenter, ZL, 1, ZLM(:,iT),1)
            End If
c
            ! compute total tensors:
            Call chi_sum( nCenter, chit_tens_ex, zstat_ex,
     &                    XL,   ZL,   XR, ZR, iopt,
     &                    chit_tens_tot(iT,1:3,1:3),     zstat_tot(iT) )

            Call chi_sum( nCenter, smu_chit_tens_ex, zstat_ex,
     &                    SMUL, ZL, SMUR, ZR, iopt,
     &                    smu_chit_tens_tot, zstat_tot(iT) )

            Call chi_sum( nCenter, ss_chit_tens_ex, zstat_ex,
     &                    SSL,  ZL,  SSR, ZR, iopt,
     &                    ss_chit_tens_tot,  zstat_tot(iT) )

            ! form the A_dir matrix:
            ! A_dir(:,:) = 1(:,:)* kB * T(iT) - zJ*XSS(:,:)
            Call daxpy_(3*3,Boltz_k*T(iT),Unity,1,A_dir,1)
            Call daxpy_(3*3,          -zJ, ss_chit_tens_tot,1,A_dir,1)
            ! invert it:
            Call REVERSE(A_dir,A_inv,DET)
            Do ic=1,3
              Do jc=1,3
                xxm=0.0_wp
                Do im=1,3
                  Do jm=1,3
                    xxm=xxm+smu_chit_tens_tot(im,ic)*a_inv(im,jm)*
     &                      smu_chit_tens_tot(jm,jc)
                  End Do
                End Do
                chit_theta_tens(iT,ic,jc)=chit_tens_tot(iT,ic,jc)+zj*xxm
              End Do ! jc
            End Do ! ic

            chit(iT)=coeff_chi * ( chit_tens_tot(iT,1,1)
     &                           + chit_tens_tot(iT,2,2)
     &                           + chit_tens_tot(iT,3,3) )/3.0_wp

            chit_theta(iT)=coeff_chi * ( chit_theta_tens(iT,1,1)
     &                                + chit_theta_tens(iT,2,2)
     &                                + chit_theta_tens(iT,3,3) )/3.0_wp
            If( abs(chit(iT)) < 1.0d-20) Then
               chit(iT)=1.d-20
               chit_theta(iT)=1.d-20
            End If
            If( abs(chit_theta(iT)) < 1.0d-20) Then
               chit_theta(iT)=1.d-20
            End If

            chi_theta_1(iT)=t(iT)/chit_theta(iT)

            ! add some verification data:
            Fa=0.0_wp; Fb=0.0_wp; Fc=0.0_wp; Fd=0.0_wp; Fe=0.0_wp;
            Ff=0.0_wp;
            Fa=dnrm2_(9,    chit_tens_tot(it,1:3,1:3),1)
            Fb=dnrm2_(9,  chit_theta_tens(it,1:3,1:3),1)
            Fc=dnrm2_(9,smu_chit_tens_tot(1:3,1:3),1)
            Fd=dnrm2_(9, ss_chit_tens_tot(1:3,1:3),1)

            Call Add_Info('XT: chit_tens_tot'    ,Fa,1,6)
            Call Add_Info('XT: chit_theta_tens'  ,Fb,1,6)
            Call Add_Info('XT: smu_chit_tens_tot',Fc,1,6)
            Call Add_Info('XT:  ss_chit_tens_tot',Fd,1,6)
         End Do ! it
         Fb=0.0_wp;
         Fb=dnrm2_(nT+nTempMagn,T,1)
         Call Add_Info('XT: T',Fb,1,6)

         Call mma_deallocate(smu_chit_tens_l)
         Call mma_deallocate(smu_chit_tens_lr)
         Call mma_deallocate(smu_chit_tens_ex)
         Call mma_deallocate(smu_chit_tens_tot)
         Call mma_deallocate(ss_chit_tens_l)
         Call mma_deallocate(ss_chit_tens_lr)
         Call mma_deallocate(ss_chit_tens_ex)
         Call mma_deallocate(ss_chit_tens_tot)
         Call mma_deallocate(SMUR)
         Call mma_deallocate(SMUL)
         Call mma_deallocate(SSR)
         Call mma_deallocate(SSL)
         Call mma_deallocate(A_dir)
         Call mma_deallocate(A_inv)
         Call mma_deallocate(Unity)

      End If  !zJ
c------------------------------------------------------------------------------------
c printing the results
      Do iT=1,nT+nTempMagn
        Do ic=1,3
          Do jc=1,3
             chit_tens_tot(iT,ic,jc)  =coeff_chi*
     &                                         chit_tens_tot(iT,ic,jc)
             If(zJ.ne.0.0_wp) Then
             chit_theta_tens(iT,ic,jc)=coeff_chi*
     &                                       chit_theta_tens(iT,ic,jc)
             End If
          End Do
        End Do
      End Do
      Write(6,*)
      Write(6,'(A)') '----------------------------------------------'//
     & '------------------------------------------|'
      Write(6,'(A)') '     |     T      | Statistical |    CHI*T    '//
     & '|    CHI*T    |     CHI     |    1/CHI    |'
      Write(6,'(A)') '     |            |  Sum (Z)    |    (zJ=0)   '//
     & '|             |             |             |'
      Write(6,'(A)') '-----|----------------------------------------'//
     & '------------------------------------------|'
      Write(6,'(A)') 'Units|   Kelvin   |    ---      |  cm3*K/mol  '//
     & '|  cm3*K/mol  |   cm3/mol   |   mol/cm3   |'
      Write(6,'(A)') '-----|----------------------------------------'//
     & '------------------------------------------|'

      Do iT=1,nT
      jT=iT+nTempMagn
      Write(6,'(A,F11.6,A,E12.5,A,F12.8,A,F12.8,A,E12.5,A,E12.5,A)')
     & '     |',         T(jT),       ' |',  zstat_tot(jT),
     &     ' |',      chiT(jT),       ' |', chiT_theta(jT),
     &     ' |',chit_theta(jT)/T(jT), ' |',chi_theta_1(jT),' |'
      End Do
      Write(6,'(A)') '-----|----------------------------------------'//
     & '------------------------------------------|'
      Fb=0.0_wp;
      Fb=dnrm2_(nT+nTempMagn,chiT,1)
      Call Add_Info('XT: T',Fb,1,6)
      Fa=0.0_wp;
      Fa=dnrm2_(nT+nTempMagn,chiT_theta,1)
      Call Add_Info('XT: CHIT_THETA',Fa,1,6)
      Fa=0.0_wp;
      Fa=dnrm2_(nT+nTempMagn,zstat_tot,1)
      Call Add_Info('XT: CHIT_THETA',Fa,1,6)
c  calcualtion of the standard deviation:
      If (tinput) Then
         Write(6,'(a,5x, f20.14)') 'ST.DEV.CHIT:',
     &        dev( (nT-nTempMagn), chit_theta( (1+nTempMagn):(nT) ),
     &                                  XTexp( (1+nTempMagn):(nT) )  )
      End If !tinput


!-------------------------  PLOTs -------------------------------------!
      WRITE(label,'(A)') "no_field"
      IF ( DoPlot ) THEN
         IF ( tinput ) THEN
            Call plot_XT_with_Exp(label, nT-nTempMagn,
     &                                 T((1+nTempMagn):(nT) ),
     &                        chit_theta((1+nTempMagn):(nT) ),
     &                             XTexp((1+nTempMagn):(nT)), zJ )
         ELSE
            Call plot_XT_no_Exp( label, nT-nTempMagn,
     &                                 T((1+nTempMagn):(nT) ),
     &                        chit_theta((1+nTempMagn):(nT) ), zJ )
         END IF
      END IF
! ------------------------- END PLOTs -------------------------------------!




c print out the main VAN VLECK SUSCEPTIBILITY TENSOR, its main values and main axes:
      Call mma_allocate(wt,3,'wt')
      Call mma_allocate(zt,3,3,'zt')
      If( zJ==0.0_wp ) Then
        Write(6,'(/)')
        Write(6,'(111A)') ('-',i=1,110),'|'
        Write(6,'(31X,A,22x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR '//
     &                           'FOR zJ = 0,  in cm3*K/mol','|'
        Write(6,'(111A)') ('-',i=1,110),'|'
        Write(6,'(A)') '     T(K)   |   |          Susceptibility '//
     &                 'Tensor      |    Main Values  |           '//
     &                 '    Main Axes             |'
        Do iT=1,nT
          jT=iT+nTempMagn
          info=0
          Call dcopy_(  3,[0.0_wp],0,wt,1)
          Call dcopy_(3*3,[0.0_wp],0,zt,1)
          Call DIAG_R2( chit_tens_tot(jT,:,:) ,3,info,wt,zt)
          Write(6,'(A)') '------------|---|'//
     &                   '------- x --------- y --------- z ---|'//
     &                   '-----------------|'//
     &                   '------ x --------- y --------- z ----|'
          Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &                   '            | x |',
     &                    (chit_tens_tot(jT,1,j),j=1,3),
     &                   ' |  X:',wt(1),'|',(zt(j,1),j=1,3),'|'
          Write(6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),
     &                               '| y |',
     &                    (chiT_tens_tot(jT,2,j),j=1,3),
     &                   ' |  Y:',wt(2),'|',(zt(j,2),j=1,3),'|'
          Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &                   '            | z |',
     &                    (chiT_tens_tot(jT,3,j),j=1,3),
     &                   ' |  Z:',wt(3),'|',(zt(j,3),j=1,3),'|'
        End Do
        Write(6,'(111A)') ('-',i=1,110),'|'

      Else ! zJ .ne. 0.0_wp

        Write(6,'(/)')
        Write(6,'(111A)') ('-',i=1,110),'|'

        Write(6,'(31X,A,F9.6,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY '//
     &                                'TENSOR FOR zJ =',zJ,
     &                                ',  in cm3*K/mol','|'
        Write(6,'(111A)') ('-',i=1,110),'|'
        Write(6,'(A)') '     T(K)   |   |          Susceptibility '//
     &                 'Tensor      |    Main Values  |           '//
     &                 '    Main Axes             |'
        Do iT=1,nT
          jT=iT+nTempMagn
          info=0
          Call dcopy_(  3,[0.0_wp],0,wt,1)
          Call dcopy_(3*3,[0.0_wp],0,zt,1)
          Call DIAG_R2( chit_theta_tens(jT,:,:) ,3,info,wt,zt)
          Write(6,'(A)') '------------|---|'//
     &                   '------- x --------- y --------- z ---|'//
     &                   '-----------------|'//
     &                   '------ x --------- y --------- z ----|'

          Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &                   '            | x |',
     &                    (chit_theta_tens(jT,1,j),j=1,3),
     &                   ' |  X:',wt(1),'|',(zt(j,1),j=1,3),'|'
          Write(6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),
     &                               '| y |',
     &                    (chiT_theta_tens(jT,2,j),j=1,3),
     &                   ' |  Y:',wt(2),'|',(zt(j,2),j=1,3),'|'
          Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &                   '            | z |',
     &                    (chiT_theta_tens(jT,3,j),j=1,3),
     &                   ' |  Z:',wt(3),'|',(zt(j,3),j=1,3),'|'
        End Do
        Write(6,'(111A)') ('-',i=1,110),'|'
      End If
      Call mma_deallocate(wt)
      Call mma_deallocate(zt)

      Call mma_deallocate(chit_tens_tot)
      Call mma_deallocate(chit_theta_tens)
      Call mma_deallocate(chit_tens_l)
      Call mma_deallocate(chit_tens_lr)
      Call mma_deallocate(chit_tens_ex)
      Call mma_deallocate(zstat_l)
      Call mma_deallocate(zstat_lr)
      Call mma_deallocate(zstat_tot)
      Call mma_deallocate(ZR)
      Call mma_deallocate(ZL)
      Call mma_deallocate(XL)
      Call mma_deallocate(XR)
      Call mma_deallocate(chiT)
      Call mma_deallocate(chi_theta_1)

      Go To 190
c------------------------------------------------------------------------------------
      Write(6,*)
      Write(6,'(5x,a)') 'on user request, the magnetic '//
     &                  'susceptibility and '
      Write(6,'(5x,a)') 'the magnetic susceptibility '//
     &                  'tensor were not calculated.'
      Write(6,*)
  190 continue
      Call qExit('PA_suscept')
      Return
      End

