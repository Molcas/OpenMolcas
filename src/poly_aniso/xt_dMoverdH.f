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
       Subroutine XT_dMoverdH( exch, nLoc, nCenter, nneq, neqv, neq,
     &                         nss, nexch, nTempMagn, nT, NM, iopt,
     &                         mem,
     &                         Tmin, Tmax, chit_exp, eso, w, T, RROT,
     &                         zJ, Xfield, EM, THRS, XT_no_field,
     &                         dipso, s_so, dipexch, s_exch,
     &                         tinput, smagn, m_paranoid, m_accurate )
c this Subroutine computes the XT as M/H as it is observed for most of experiments.
c The M is averaged over the grid for each temperature point.
c       chi*t ----------- the units are cgsemu: [ cm^3*k/mol ]

      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "mgrid.fh"
#include "stdalloc.fh"
      Integer, intent(in)          :: exch, nLoc, nCenter, nneq, neqv
      Integer, intent(in)          :: nTempMagn, nT, NM, iopt, mem
      Integer, intent(in)          :: neq(nneq), nss(nneq), nexch(nneq)


      Real(kind=wp), intent(in)    :: Tmin, Tmax
      Real(kind=wp), intent(in)    :: eso(nneq,nLoc)
      Real(kind=wp), intent(in)    :: W(exch)
      Real(kind=wp), intent(in)    :: zJ
      Real(kind=wp), intent(in)    :: RROT(nneq,neqv,3,3)
      Real(kind=wp), intent(in)    :: T(nT+nTempMagn)
      Real(kind=wp), intent(in)    :: chit_exp(nT)
      Real(kind=wp), intent(in)    :: XT_no_field(nT+nTempMagn)
      Real(kind=wp), intent(in)    :: Xfield
      Real(kind=wp), intent(in)    :: EM
      Real(kind=wp), intent(in)    :: THRS
      Complex(kind=wp), intent(in) :: dipso(nneq,3,nLoc,nLoc)
      Complex(kind=wp), intent(in) ::  s_so(nneq,3,nLoc,nLoc)
      Complex(kind=wp), intent(in) :: dipexch(3,exch,exch)
      Complex(kind=wp), intent(in) ::  s_exch(3,exch,exch)
      Logical, intent(in)          :: tinput, smagn, m_paranoid
cccc local variables ccc
      Integer                      :: nDirX
      Integer       :: iM,iT,jT,i,j,l,n,isite !,nP

      Real(kind=wp), allocatable :: WEX0(:)
      Real(kind=wp), allocatable :: WEX1(:)
      Real(kind=wp), allocatable :: WEX2(:)   ! Zeeman exchange energies
      Real(kind=wp), allocatable :: WL0(:,:)    ! Zeeman local energies
      Real(kind=wp), allocatable :: WL1(:,:)
      Real(kind=wp), allocatable :: WL2(:,:)    ! Zeeman local energies
!     Zeeman local reduced energies, using only NEXCH states;
      Real(kind=wp), allocatable :: WR0(:,:)
      Real(kind=wp), allocatable :: WR1(:,:)
!     Zeeman local reduced energies, using only NEXCH states;
      Real(kind=wp), allocatable :: WR2(:,:)

!     local statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ZL0(:,:)
!     local statistical sum, Boltzmann distribution, using only NEXCH states
      Real(kind=wp), allocatable :: ZR0(:,:)
!     local statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ZL1(:,:)
!     local statistical sum, Boltzmann distribution, using only NEXCH states
      Real(kind=wp), allocatable :: ZR1(:,:)
!     local statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ZL2(:,:)
!     local statistical sum, Boltzmann distribution, using only NEXCH states
      Real(kind=wp), allocatable :: ZR2(:,:)
!     spin magnetisation, from the local sites, using ALL states ;
      Real(kind=wp), allocatable :: SL0(:,:,:)
!     spin magnetisation, from the local sites, using only NEXCH states ;
      Real(kind=wp), allocatable :: SR0(:,:,:)
!     spin magnetisation, from the local sites, using ALL states ;
      Real(kind=wp), allocatable :: SL1(:,:,:)
!     spin magnetisation, from the local sites, using only NEXCH states ;
      Real(kind=wp), allocatable :: SR1(:,:,:)
!     spin magnetisation, from the local sites, using ALL states ;
      Real(kind=wp), allocatable :: SL2(:,:,:)
!     spin magnetisation, from the local sites, using only NEXCH states ;
      Real(kind=wp), allocatable :: SR2(:,:,:)
!     magnetisation, from local sites, using ALL states;
      Real(kind=wp), allocatable :: ML0(:,:,:)
!     magnetisation, from local sites, using only NEXCH states;
      Real(kind=wp), allocatable :: MR0(:,:,:)
!     magnetisation, from local sites, using ALL states;
      Real(kind=wp), allocatable :: ML1(:,:,:)
!     magnetisation, from local sites, using only NEXCH states;
      Real(kind=wp), allocatable :: MR1(:,:,:)
!     magnetisation, from local sites, using ALL states;
      Real(kind=wp), allocatable :: ML2(:,:,:)
!     magnetisation, from local sites, using only NEXCH states;
      Real(kind=wp), allocatable :: MR2(:,:,:)

!     spin magnetisation, from the exchange block;
      Real(kind=wp), allocatable :: SEX0(:,:)
!     spin magnetisation, from the exchange block;
      Real(kind=wp), allocatable :: SEX1(:,:)
!     spin magnetisation, from the exchange block;
      Real(kind=wp), allocatable :: SEX2(:,:)
!     magnetisation, form the exchange block
      Real(kind=wp), allocatable :: MEX0(:,:)
!     magnetisation, form the exchange block
      Real(kind=wp), allocatable :: MEX1(:,:)
!     magnetisation, form the exchange block
      Real(kind=wp), allocatable :: MEX2(:,:)
!     exchange statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ZEX0(:)
!     exchange statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ZEX1(:)
!     exchange statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ZEX2(:)
c total vectors in general coordinate system:
      Real(kind=wp), allocatable :: ZRT0(:,:)
      Real(kind=wp), allocatable :: ZLT0(:,:)
      Real(kind=wp), allocatable :: ZRT1(:,:)
      Real(kind=wp), allocatable :: ZLT1(:,:)
      Real(kind=wp), allocatable :: ZRT2(:,:)
      Real(kind=wp), allocatable :: ZLT2(:,:)
      Real(kind=wp), allocatable :: MRT0(:,:,:)
      Real(kind=wp), allocatable :: MLT0(:,:,:)
      Real(kind=wp), allocatable :: SRT0(:,:,:)
      Real(kind=wp), allocatable :: SLT0(:,:,:)
      Real(kind=wp), allocatable :: MRT1(:,:,:)
      Real(kind=wp), allocatable :: MLT1(:,:,:)
      Real(kind=wp), allocatable :: SRT1(:,:,:)
      Real(kind=wp), allocatable :: SLT1(:,:,:)
      Real(kind=wp), allocatable :: MRT2(:,:,:)
      Real(kind=wp), allocatable :: MLT2(:,:,:)
      Real(kind=wp), allocatable :: SRT2(:,:,:)
      Real(kind=wp), allocatable :: SLT2(:,:,:)
c data for total system:
! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ZT0(:)
! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ZT1(:)
! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ZT2(:)
      Real(kind=wp), allocatable :: MT0(:,:) ! total magnetisation
      Real(kind=wp), allocatable :: MT1(:,:) ! total magnetisation
      Real(kind=wp), allocatable :: ST0(:,:) ! total spin magnetisation,
      Real(kind=wp), allocatable :: ST1(:,:) ! total spin magnetisation,
      Real(kind=wp), allocatable :: ST2(:,:) ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT2(:,:) ! total magnetisation
c standard deviation data:
      Real(kind=wp)              :: dev
      Real(kind=wp), allocatable :: XTM_MH(:)
      Real(kind=wp), allocatable :: XTM_dMdH(:)
      Real(kind=wp), allocatable :: XTtens_MH(:,:,:)
      Real(kind=wp), allocatable :: XTtens_dMdH(:,:,:)

      Integer                    :: nTempTotal
      Real(kind=wp)              :: Xfield_1, Xfield_2, dltXF
      !Real(kind=wp)             :: XTS(nT+nTempMagn)
c Zeeman energy and M vector
c      Real(kind=wp)             :: dirX(nDir), dirY(nDir), dirZ(nDir)
c      Real(kind=wp)             :: dir_weight(nDirZee,3)

      Real(kind=wp)              :: dHX(3)
      Real(kind=wp)              :: dHY(3)
      Real(kind=wp)              :: dHZ(3)
      Real(kind=wp)              :: dHW(3)
      Integer                    :: RtoB, mem_local

      Integer                    :: info, ibuf, ibuf1, ibuf3
      Real(kind=wp)              :: wt(3), zt(3,3)

      Real(kind=wp)              :: dnrm2_
      External                   :: dev, dnrm2_
      Real(kind=wp)              :: cm3tomB
      logical                    :: DBG, m_accurate
      Call qEnter('PA_XTdMdH')
      DBG=.false.
      !m_paranoid=.true.!.false.
      cm3tomB=0.55849389040_wp   !   in cm3 * mol-1 * T
cccc-------------------------------------------------------cccc
      If(DBG) Then
         Write(6,'(A,   I5)') '      exch =',exch
         Write(6,'(A,   I5)') '      nLoc =',nLoc
         Write(6,'(A,   I5)') '   nCenter =',nCenter
         Write(6,'(A,   I5)') '      nneq =',nneq
         Write(6,'(A,   I5)') '      neqv =',neqv
         Write(6,'(A, 10I5)') '    neq(i) =',(neq(i),i=1,nneq)
         Write(6,'(A, 10I5)') '  nexch(i) =',(nexch(i),i=1,nneq)
         Write(6,'(A, 10I5)') '    nss(i) =',(nss(i),i=1,nneq)
         Write(6,'(A,   I5)') ' nTempMagn =',nTempMagn
         Write(6,'(A,   I5)') '        nT =',nT
         Write(6,'(A,   I5)') '        nM =',nM
         Write(6,'(A,   I5)') '      iopt =',iopt
         Write(6,'(A,F12.5)') '      Tmin =',Tmin
         Write(6,'(A,F12.5)') '      Tmax =',Tmax
         Write(6,'(A,F12.5)') '        zJ =',zJ
         Write(6,'(A,F12.5)') '    Xfield =',Xfield
         Write(6,'(A,F12.5)') '        EM =',EM
         Write(6,'(A,E12.5)') '      THRS =',THRS
         Write(6,          *) 'm_accurate =',m_accurate
         Write(6,          *) '     smagn =',smagn
         Write(6,          *) 'm_paranoid =',m_paranoid
         Write(6,          *) '    tinput =',tinput
         If (tinput) Then
            Do iT=1,nTempMagn
               Write(6,'(2(A,i3,A,F12.6,2x))') 'T(',iT,')=',T(iT)
            End Do
            Do iT=1,nT
               jT=iT+nTempMagn
               Write(6,'(2(A,i3,A,F12.6,2x))') 'T(',jT,')=',T(jT),
     &                          ' chiT_exp(',iT,')=',chit_exp(iT)
           End Do

         Else
           Do iT=1,nT+nTempMagn
             Write(6,'(2(A,i3,A,F12.6,2x))') 'T(',iT,')=',T(iT)
           End Do
         End If
         Write(6,'(A)') 'local ESO:'
         Do i=1,nneq
            Write(6,'(A,i3)') 'site:',i
            Write(6,'(10F12.5)') (eso(i,j),j=1,nss(i))
         End Do
         Write(6,'(A)') 'EXCHANGE ENERGY'
         Write(6,'(10F12.5)') (W(i),i=1,exch)
         Write(6,'(A)') 'Rotation Matrices:'
         Do i=1,nneq
            Write(6,'(A,i3)') 'site:',i
            Do j=1,neq(i)
               Do l=1,3
                  Write(6,'(10F12.5)') (RROT(i,j,l,n),n=1,3)
               End Do
            End Do
         End Do
         Call xFlush(6)
      End If !DBG
cccc-------------------------------------------------------cccc
      Write(6,*)
      Write(6,'(100A)') (('%'),J=1,95)
      Write(6,'(16X,A)') 'CALCULATION OF THE FIELD-DEPENDENT '//
     &                   'MAGNETIC SUSCEPTIBILITY'
      Write(6,'(18X,A)') 'within true (dM/dH) and "experimentalists"'//
     &                   ' (M/H) models'
      Write(6,'(100A)') (('%'),J=1,95)
      Write(6,*)
      Write(6,'(2x,A,F10.6,A)') 'Magnetic field strength:',Xfield,
     &                          ' Tesla.'
      If(tinput) Then
         Write(6,'(2x,a)') 'Temperature dependence of the magnetic '//
     &                     'susceptibility and'
         Write(6,'(2x,a)') 'high-field magnetization will be '//
     &                     'calculated according to '
         Write(6,'(2x,a)') 'experimental values provided by the '//
     &                     'user in file "chitexp.input".'
      Else
         Write(6,'(2x,a,i3,a)') 'Temperature dependence of the '//
     &                          'magnetic susceptibility will be '//
     &                          'calculated in',nT,' points, '
         Write(6,'(2x,a,f4.1,a,f6.1,a)') 'equally distributed in '//
     &                          'temperature range ',tmin,' ---',
     &                           tmax,' K.'
      End If

!=======================================================================
      mem_local=0
      RtoB=8
      If(nM>0) Then
        Call mma_allocate(WEX0,nM,'WEX0')
        Call mma_allocate(WEX1,nM,'WEX1')
        Call mma_allocate(WEX2,nM,'WEX2')
        Call dcopy_(nM,0.0_wp,0,WEX0,1)
        Call dcopy_(nM,0.0_wp,0,WEX1,1)
        Call dcopy_(nM,0.0_wp,0,WEX2,1)
        mem_local=mem_local+3*nM*RtoB
      End If
      If((nneq>0).and.(nLoc>0)) Then
        Call mma_allocate(WL0,nneq,nLoc,'WL0')
        Call mma_allocate(WL1,nneq,nLoc,'WL1')
        Call mma_allocate(WL2,nneq,nLoc,'WL2')
        Call mma_allocate(WR0,nneq,nLoc,'WR0')
        Call mma_allocate(WR1,nneq,nLoc,'WR1')
        Call mma_allocate(WR2,nneq,nLoc,'WR2')
        Call dcopy_(nneq*nLoc,0.0_wp,0,WL0,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WL1,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WL2,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WR0,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WR1,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WR2,1)
        mem_local=mem_local+6*nneq*nLoc*RtoB
      End If

      If((nneq>0).and.((nT+nTempMagn)>0)) Then
        Call mma_allocate(ZL0,nneq,(nT+nTempMagn),'ZL0')
        Call mma_allocate(ZR0,nneq,(nT+nTempMagn),'ZR0')
        Call mma_allocate(ZL1,nneq,(nT+nTempMagn),'ZL1')
        Call mma_allocate(ZR1,nneq,(nT+nTempMagn),'ZR1')
        Call mma_allocate(ZL2,nneq,(nT+nTempMagn),'ZL2')
        Call mma_allocate(ZR2,nneq,(nT+nTempMagn),'ZR2')
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZL0,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZR0,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZL1,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZR1,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZL2,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZR2,1)
        mem_local=mem_local+6*nneq*(nT+nTempMagn)*RtoB

        Call mma_allocate(SL0,nneq,3,(nT+nTempMagn),'SL0')
        Call mma_allocate(SR0,nneq,3,(nT+nTempMagn),'SR0')
        Call mma_allocate(SL1,nneq,3,(nT+nTempMagn),'SL1')
        Call mma_allocate(SR1,nneq,3,(nT+nTempMagn),'SR1')
        Call mma_allocate(SL2,nneq,3,(nT+nTempMagn),'SL2')
        Call mma_allocate(SR2,nneq,3,(nT+nTempMagn),'SR2')
        Call mma_allocate(ML0,nneq,3,(nT+nTempMagn),'ML0')
        Call mma_allocate(MR0,nneq,3,(nT+nTempMagn),'MR0')
        Call mma_allocate(ML1,nneq,3,(nT+nTempMagn),'ML1')
        Call mma_allocate(MR1,nneq,3,(nT+nTempMagn),'MR1')
        Call mma_allocate(ML2,nneq,3,(nT+nTempMagn),'ML2')
        Call mma_allocate(MR2,nneq,3,(nT+nTempMagn),'MR2')
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SL0,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SR0,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SL1,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SR1,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SL2,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SR2,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,ML0,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,MR0,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,ML1,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,MR1,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,ML2,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,MR2,1)
        mem_local=mem_local+12*3*nneq*(nT+nTempMagn)*RtoB
      End If

      If((nT+nTempMagn)>0) Then
        Call mma_allocate(ZEX0,(nT+nTempMagn),'ZEX0')
        Call mma_allocate(ZEX1,(nT+nTempMagn),'ZEX1')
        Call mma_allocate(ZEX2,(nT+nTempMagn),'ZEX2')
        Call dcopy_((nT+nTempMagn),0.0_wp,0,ZEX0,1)
        Call dcopy_((nT+nTempMagn),0.0_wp,0,ZEX1,1)
        Call dcopy_((nT+nTempMagn),0.0_wp,0,ZEX2,1)
        mem_local=mem_local+3*(nT+nTempMagn)*RtoB

        Call mma_allocate(ZT0,(nT+nTempMagn),'ZT0')
        Call mma_allocate(ZT1,(nT+nTempMagn),'ZT1')
        Call mma_allocate(ZT2,(nT+nTempMagn),'ZT2')
        Call dcopy_((nT+nTempMagn),0.0_wp,0,ZEX0,1)
        Call dcopy_((nT+nTempMagn),0.0_wp,0,ZEX1,1)
        Call dcopy_((nT+nTempMagn),0.0_wp,0,ZEX2,1)
        mem_local=mem_local+3*(nT+nTempMagn)*RtoB

        Call mma_allocate(SEX0,3,(nT+nTempMagn),'SEX0')
        Call mma_allocate(SEX1,3,(nT+nTempMagn),'SEX1')
        Call mma_allocate(SEX2,3,(nT+nTempMagn),'SEX2')
        Call mma_allocate(MEX0,3,(nT+nTempMagn),'MEX0')
        Call mma_allocate(MEX1,3,(nT+nTempMagn),'MEX1')
        Call mma_allocate(MEX2,3,(nT+nTempMagn),'MEX2')
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,SEX0,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,SEX1,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,SEX2,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MEX0,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MEX1,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MEX2,1)
        mem_local=mem_local+6*3*(nT+nTempMagn)*RtoB

        Call mma_allocate(ST0,3,(nT+nTempMagn),'ST0')
        Call mma_allocate(ST1,3,(nT+nTempMagn),'ST1')
        Call mma_allocate(ST2,3,(nT+nTempMagn),'ST2')
        Call mma_allocate(MT0,3,(nT+nTempMagn),'MT0')
        Call mma_allocate(MT1,3,(nT+nTempMagn),'MT1')
        Call mma_allocate(MT2,3,(nT+nTempMagn),'MT2')
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,ST0,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,ST1,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,ST2,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MT0,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MT1,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MT2,1)
        mem_local=mem_local+6*3*(nT+nTempMagn)*RtoB

        Call mma_allocate(XTM_MH  ,(nT+nTempMagn),'XTM_MH')
        Call mma_allocate(XTM_dMdH,(nT+nTempMagn),'XTM_dMdH')
        Call mma_allocate(XTtens_MH  ,3,3,(nT+nTempMagn),'XTtens_MH')
        Call mma_allocate(XTtens_dMdH,3,3,(nT+nTempMagn),'XTtens_dMdH')
        Call dcopy_((nT+nTempMagn),0.0_wp,0,XTM_MH,1)
        Call dcopy_((nT+nTempMagn),0.0_wp,0,XTM_dMdH,1)
        Call dcopy_(3*3*(nT+nTempMagn),0.0_wp,0,XTtens_MH,1)
        Call dcopy_(3*3*(nT+nTempMagn),0.0_wp,0,XTtens_dMdH,1)
        mem_local=mem_local+20*(nT+nTempMagn)*RtoB
      End If

      If((nCenter>0).and.((nT+nTempMagn)>0)) Then
        Call mma_allocate(ZRT0,nCenter,(nT+nTempMagn),'ZRT0')
        Call mma_allocate(ZLT0,nCenter,(nT+nTempMagn),'ZLT0')
        Call mma_allocate(ZRT1,nCenter,(nT+nTempMagn),'ZRT1')
        Call mma_allocate(ZLT1,nCenter,(nT+nTempMagn),'ZLT1')
        Call mma_allocate(ZRT2,nCenter,(nT+nTempMagn),'ZRT2')
        Call mma_allocate(ZLT2,nCenter,(nT+nTempMagn),'ZLT2')
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZRT0,1)
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZLT0,1)
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZRT1,1)
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZLT1,1)
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZRT2,1)
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZLT2,1)
        mem_local=mem_local+6*3*nCenter*(nT+nTempMagn)*RtoB

        Call mma_allocate(MRT0,nCenter,3,(nT+nTempMagn),'MRT0')
        Call mma_allocate(MLT0,nCenter,3,(nT+nTempMagn),'MLT0')
        Call mma_allocate(SRT0,nCenter,3,(nT+nTempMagn),'SRT0')
        Call mma_allocate(SLT0,nCenter,3,(nT+nTempMagn),'SLT0')
        Call mma_allocate(MRT1,nCenter,3,(nT+nTempMagn),'MRT1')
        Call mma_allocate(MLT1,nCenter,3,(nT+nTempMagn),'MLT1')
        Call mma_allocate(SRT1,nCenter,3,(nT+nTempMagn),'SRT1')
        Call mma_allocate(SLT1,nCenter,3,(nT+nTempMagn),'SLT1')
        Call mma_allocate(MRT2,nCenter,3,(nT+nTempMagn),'MRT2')
        Call mma_allocate(MLT2,nCenter,3,(nT+nTempMagn),'MLT2')
        Call mma_allocate(SRT2,nCenter,3,(nT+nTempMagn),'SRT2')
        Call mma_allocate(SLT2,nCenter,3,(nT+nTempMagn),'SLT2')
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MRT0,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MLT0,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SRT0,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SLT0,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MRT1,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MLT1,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SRT1,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SLT1,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MRT2,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MLT2,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SRT2,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SLT2,1)
        mem_local=mem_local+12*3*nCenter*(nT+nTempMagn)*RtoB
      End If

      If(dbg) Write(6,*) 'XT_MH:  memory allocated (local):'
      If(dbg) Write(6,*) 'mem_local=', mem_local
      If(dbg) Write(6,*) 'XT_MH:  memory allocated (total):'
      If(dbg) Write(6,*) 'mem_total=', mem+mem_local

!=======================================================================
      dHX=0.0_wp
      dHY=0.0_wp
      dHZ=0.0_wp
      dHW=0.0_wp

      nDirX=3

      dHX(1)=1.0_wp
      dHY(1)=0.0_wp
      dHZ(1)=0.0_wp

      dHX(2)=0.0_wp
      dHY(2)=1.0_wp
      dHZ(2)=0.0_wp

      dHX(3)=0.0_wp
      dHY(3)=0.0_wp
      dHZ(3)=1.0_wp

      Xfield_1=0.0_wp
      Xfield_2=0.0_wp
      dltXF=0.0_wp
      Xfield_1=Xfield*0.99999_wp
      Xfield_2=Xfield*1.00001_wp
      dltXF=Xfield_2-Xfield_1

      nTempTotal=nT+nTempMagn

c ///  opening the loop over different directions of the magnetic field
      Do iM=1,nDirX
        Call dcopy_(nM,0.0_wp,0,WEX0,1)
        Call dcopy_(nM,0.0_wp,0,WEX1,1)
        Call dcopy_(nM,0.0_wp,0,WEX2,1)
        Call dcopy_((nT+nTempMagn),0.0_wp,0,ZEX0,1)
        Call dcopy_((nT+nTempMagn),0.0_wp,0,ZEX1,1)
        Call dcopy_((nT+nTempMagn),0.0_wp,0,ZEX2,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,SEX0,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,SEX1,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,SEX2,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MEX0,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MEX1,1)
        Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MEX2,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WL0,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WL1,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WL2,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WR0,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WR1,1)
        Call dcopy_(nneq*nLoc,0.0_wp,0,WR2,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZL0,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZR0,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZL1,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZR1,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZL2,1)
        Call dcopy_(nneq*(nT+nTempMagn),0.0_wp,0,ZR2,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SL0,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SR0,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SL1,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SR1,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SL2,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,SR2,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,ML0,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,MR0,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,ML1,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,MR1,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,ML2,1)
        Call dcopy_(3*nneq*(nT+nTempMagn),0.0_wp,0,MR2,1)

        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZRT0,1)
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZLT0,1)
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZRT1,1)
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZLT1,1)
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZRT2,1)
        Call dcopy_(nCenter*(nT+nTempMagn),0.0_wp,0,ZLT2,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MRT0,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MLT0,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SRT0,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SLT0,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MRT1,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MLT1,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SRT1,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SLT1,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MRT2,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,MLT2,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SRT2,1)
        Call dcopy_(3*nCenter*(nT+nTempMagn),0.0_wp,0,SLT2,1)
c exchange magnetization:
         Call MAGN( EXCH, NM, dHX(iM), dHY(iM), dHZ(iM),
     &              XField, W, zJ, THRS,
     &              DIPEXCH(1:3,1:exch,1:exch),
     &               S_EXCH(1:3,1:exch,1:exch),
     &              nTempTotal, T(1:nTempTotal), smagn,
     &              Wex0(1:nM),
     &              Zex0(1:nTempTotal),
     &              Sex0(1:3,1:nTempTotal),
     &              Mex0(1:3,1:nTempTotal), m_paranoid, DBG )
         Call MAGN( EXCH, NM, dHX(iM), dHY(iM), dHZ(iM),
     &              XField_1, W, zJ, THRS,
     &              DIPEXCH(1:3,1:exch,1:exch),
     &               S_EXCH(1:3,1:exch,1:exch),
     &              nTempTotal, T(1:nTempTotal), smagn,
     &              Wex1(1:nM),
     &              Zex1(1:nTempTotal),
     &              Sex1(1:3,1:nTempTotal),
     &              Mex1(1:3,1:nTempTotal), m_paranoid, DBG )
         Call MAGN( EXCH, NM, dHX(iM), dHY(iM), dHZ(iM),
     &              XField_2, W, zJ, THRS,
     &              DIPEXCH(1:3,1:exch,1:exch),
     &               S_EXCH(1:3,1:exch,1:exch),
     &              nTempTotal, T(1:nTempTotal), smagn,
     &              Wex2(1:nM),
     &              Zex2(1:nTempTotal),
     &              Sex2(1:3,1:nTempTotal),
     &              Mex2(1:3,1:nTempTotal), m_paranoid, DBG )
      If(DBG) Then
        Do i=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex0(',i,')=',(Mex0(l,i),l=1,3)
        End Do
        Do i=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex1(',i,')=',(Mex1(l,i),l=1,3)
        End Do
        Do i=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex2(',i,')=',(Mex2(l,i),l=1,3)
        End Do
        Do i=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex 2-1 (',i,')=',
     &                                (Mex2(l,i)-Mex1(l,i),l=1,3)
        End Do
      End If !DBG

       Call Add_Info('dM/dH   Mex0',dnrm2_(3*nTempTotal,Mex0,1),1,8)
       Call Add_Info('dM/dH   Sex0',dnrm2_(3*nTempTotal,Sex0,1),1,8)
       Call Add_Info('dM/dH   Zex0',dnrm2_(  nTempTotal,Zex0,1),1,8)
       Call Add_Info('dM/dH   Wex0',dnrm2_(          nM,Wex0,1),1,8)

       Call Add_Info('dM/dH   Mex1',dnrm2_(3*nTempTotal,Mex1,1),1,8)
       Call Add_Info('dM/dH   Sex1',dnrm2_(3*nTempTotal,Sex1,1),1,8)
       Call Add_Info('dM/dH   Zex1',dnrm2_(  nTempTotal,Zex1,1),1,8)
       Call Add_Info('dM/dH   Wex1',dnrm2_(          nM,Wex1,1),1,8)

       Call Add_Info('dM/dH   Mex2',dnrm2_(3*nTempTotal,Mex2,1),1,8)
       Call Add_Info('dM/dH   Sex2',dnrm2_(3*nTempTotal,Sex2,1),1,8)
       Call Add_Info('dM/dH   Zex2',dnrm2_(  nTempTotal,Zex2,1),1,8)
       Call Add_Info('dM/dH   Wex2',dnrm2_(          nM,Wex2,1),1,8)

c compute local magnetizations:
         If(m_accurate) Then
            Do i=1,nneq
c all states:
               If ( NSS(i).gt.NEXCH(i) )  Then
c this check is to avoid the unnecessary computation, in cases when no local excited states are present
                  Call MAGN( NSS(i), NEXCH(i), dHX(iM),dHY(iM),dHZ(iM),
     &                       XField,
     &                       ESO(i,1:NSS(i)), zJ, THRS,
     &                       DIPSO(i,1:3,1:NSS(i),1:NSS(i)),
     &                        S_SO(i,1:3,1:NSS(i),1:NSS(i)),
     &                       nTempTotal, T(1:(nT+nTempMagn)), smagn,
     &                       WL0(i,1:NSS(i)),
     &                       ZL0(i,1:(nT+nTempMagn)),
     &                       SL0(i,1:3,1:(nT+nTempMagn)),
     &                       ML0(i,1:3,1:(nT+nTempMagn)),
     &                       m_paranoid, DBG )

                  Call MAGN( NSS(i), NEXCH(i), dHX(iM),dHY(iM),dHZ(iM),
     &                       XField_1,
     &                       ESO(i,1:NSS(i)), zJ, THRS,
     &                       DIPSO(i,1:3,1:NSS(i),1:NSS(i)),
     &                        S_SO(i,1:3,1:NSS(i),1:NSS(i)),
     &                       nTempTotal, T(1:(nT+nTempMagn)), smagn,
     &                       WL1(i,1:NSS(i)),
     &                       ZL1(i,1:(nT+nTempMagn)),
     &                       SL1(i,1:3,1:(nT+nTempMagn)),
     &                       ML1(i,1:3,1:(nT+nTempMagn)),
     &                       m_paranoid, DBG )

                  Call MAGN( NSS(i), NEXCH(i), dHX(iM),dHY(iM),dHZ(iM),
     &                       Xfield_2,
     &                       ESO(i,1:NSS(i)), zJ, THRS,
     &                       DIPSO(i,1:3,1:NSS(i),1:NSS(i)),
     &                        S_SO(i,1:3,1:NSS(i),1:NSS(i)),
     &                       nTempTotal, T(1:(nT+nTempMagn)), smagn,
     &                       WL2(i,1:NSS(i)),
     &                       ZL2(i,1:(nT+nTempMagn)),
     &                       SL2(i,1:3,1:(nT+nTempMagn)),
     &                       ML2(i,1:3,1:(nT+nTempMagn)),
     &                       m_paranoid, DBG )
c only local "exchange states":
                  Call MAGN( NEXCH(i), NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),
     &                       XField,
     &                       ESO(i,1:NEXCH(i)), zJ, THRS,
     &                       DIPSO( i,1:3,1:NEXCH(i),1:NEXCH(i) ),
     &                        S_SO( i,1:3,1:NEXCH(i),1:NEXCH(i) ),
     &                       nTempTotal, T(1:(nT+nTempMagn)), smagn,
     &                       WR0(i,1:Nexch(i)),
     &                       ZR0(i,1:(nT+nTempMagn)),
     &                       SR0(i,1:3,1:(nT+nTempMagn)),
     &                       MR0(i,1:3,1:(nT+nTempMagn)), m_paranoid )

                  Call MAGN( NEXCH(i), NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),
     &                       XField_1,
     &                       ESO(i,1:NEXCH(i)), zJ, THRS,
     &                       DIPSO( i,1:3,1:NEXCH(i),1:NEXCH(i) ),
     &                        S_SO( i,1:3,1:NEXCH(i),1:NEXCH(i) ),
     &                       nTempTotal, T(1:(nT+nTempMagn)), smagn,
     &                       WR1(i,1:Nexch(i)),
     &                       ZR1(i,1:(nT+nTempMagn)),
     &                       SR1(i,1:3,1:(nT+nTempMagn)),
     &                       MR1(i,1:3,1:(nT+nTempMagn)),
     &                       m_paranoid, DBG )

                  Call MAGN( NEXCH(i), NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),
     &                       XField_2,
     &                       ESO(i,1:NEXCH(i)), zJ, THRS,
     &                       DIPSO( i,1:3,1:NEXCH(i),1:NEXCH(i) ),
     &                        S_SO( i,1:3,1:NEXCH(i),1:NEXCH(i) ),
     &                       nTempTotal, T(1:(nT+nTempMagn)), smagn,
     &                       WR2(i,1:Nexch(i)),
     &                       ZR2(i,1:(nT+nTempMagn)),
     &                       SR2(i,1:3,1:(nT+nTempMagn)),
     &                       MR2(i,1:3,1:(nT+nTempMagn)),
     &                       m_paranoid, DBG )
                  If (DBG) Then
                     Do iT=1,nTempTotal
                        Write(6,'(A,i1,A,i3,A,3E22.14)')
     &                          'ML(',i,',L,',iT,')=',
     &                                   (ML2(i,l,iT)-ML1(i,l,iT),l=1,3)
                        Write(6,'(A,i1,A,i3,A,3E22.14)')
     &                          'MR(',i,',L,',iT,')=',
     &                                   (MR2(i,l,iT)-MR1(i,l,iT),l=1,3)
                     End Do
                  End If ! DBG
               End If ! NSS(i).gt.NEXCH(i)
            End Do ! i=1,nneq


       ibuf1=nTempTotal*nneq
       ibuf3=3*ibuf1
       Call Add_Info('dM/dH    ML0',dnrm2_(  ibuf3,ML0,1),1,8)
       Call Add_Info('dM/dH    SL0',dnrm2_(  ibuf3,SL0,1),1,8)
       Call Add_Info('dM/dH    ZL0',dnrm2_(  ibuf1,ZL0,1),1,8)
       Call Add_Info('dM/dH    WL0',dnrm2_(nM*nneq,WL0,1),1,8)

       Call Add_Info('dM/dH    MR0',dnrm2_(  ibuf3,MR0,1),1,8)
       Call Add_Info('dM/dH    SR0',dnrm2_(  ibuf3,SR0,1),1,8)
       Call Add_Info('dM/dH    ZR0',dnrm2_(  ibuf1,ZR0,1),1,8)
       Call Add_Info('dM/dH    WR0',dnrm2_(nM*nneq,WR0,1),1,8)

       Call Add_Info('dM/dH    ML1',dnrm2_(  ibuf3,ML1,1),1,8)
       Call Add_Info('dM/dH    SL1',dnrm2_(  ibuf3,SL1,1),1,8)
       Call Add_Info('dM/dH    ZL1',dnrm2_(  ibuf1,ZL1,1),1,8)
       Call Add_Info('dM/dH    WL1',dnrm2_(nM*nneq,WL1,1),1,8)

       Call Add_Info('dM/dH    MR1',dnrm2_(  ibuf3,MR1,1),1,8)
       Call Add_Info('dM/dH    SR1',dnrm2_(  ibuf3,SR1,1),1,8)
       Call Add_Info('dM/dH    ZR1',dnrm2_(  ibuf1,ZR1,1),1,8)
       Call Add_Info('dM/dH    WR1',dnrm2_(nM*nneq,WR1,1),1,8)

       Call Add_Info('dM/dH    ML2',dnrm2_(  ibuf3,ML2,1),1,8)
       Call Add_Info('dM/dH    SL2',dnrm2_(  ibuf3,SL2,1),1,8)
       Call Add_Info('dM/dH    ZL2',dnrm2_(  ibuf1,ZL2,1),1,8)
       Call Add_Info('dM/dH    WL2',dnrm2_(nM*nneq,WL2,1),1,8)

       Call Add_Info('dM/dH    MR2',dnrm2_(  ibuf3,MR2,1),1,8)
       Call Add_Info('dM/dH    SR2',dnrm2_(  ibuf3,SR2,1),1,8)
       Call Add_Info('dM/dH    ZR2',dnrm2_(  ibuf1,ZR2,1),1,8)
       Call Add_Info('dM/dH    WR2',dnrm2_(nM*nneq,WR2,1),1,8)

c expand the basis and rotate local vectors to the general
c coordinate system:

            isite=0
            Do i=1,NNEQ
               Do j=1,NEQ(i)
                  isite=isite+1
c statistical distributions
                  Do iT=1,nTempTotal
                     ZLT0(isite,iT)=ZL0(i,iT)
                     ZRT0(isite,iT)=ZR0(i,iT)
                     ZLT1(isite,iT)=ZL1(i,iT)
                     ZRT1(isite,iT)=ZR1(i,iT)
                     ZLT2(isite,iT)=ZL2(i,iT)
                     ZRT2(isite,iT)=ZR2(i,iT)
                  End Do
c magnetizations:
c    use R_rot matrices, which have determinant +1.
c  note that  R_lg matrices have arbitrary determinant.
                  Do iT=1,nTempTotal
                     Do l=1,3
                        Do n=1,3
        MLT0(isite,l,iT) = MLT0(isite,l,iT) + rrot(i,j,l,n)*ML0(i,n,iT)
        SLT0(isite,l,iT) = SLT0(isite,l,iT) + rrot(i,j,l,n)*SL0(i,n,iT)
        MRT0(isite,l,iT) = MRT0(isite,l,iT) + rrot(i,j,l,n)*MR0(i,n,iT)
        SRT0(isite,l,iT) = SRT0(isite,l,iT) + rrot(i,j,l,n)*SR0(i,n,iT)

        MLT1(isite,l,iT) = MLT1(isite,l,iT) + rrot(i,j,l,n)*ML1(i,n,iT)
        SLT1(isite,l,iT) = SLT1(isite,l,iT) + rrot(i,j,l,n)*SL1(i,n,iT)
        MRT1(isite,l,iT) = MRT1(isite,l,iT) + rrot(i,j,l,n)*MR1(i,n,iT)
        SRT1(isite,l,iT) = SRT1(isite,l,iT) + rrot(i,j,l,n)*SR1(i,n,iT)

        MLT2(isite,l,iT) = MLT2(isite,l,iT) + rrot(i,j,l,n)*ML2(i,n,iT)
        SLT2(isite,l,iT) = SLT2(isite,l,iT) + rrot(i,j,l,n)*SL2(i,n,iT)
        MRT2(isite,l,iT) = MRT2(isite,l,iT) + rrot(i,j,l,n)*MR2(i,n,iT)
        SRT2(isite,l,iT) = SRT2(isite,l,iT) + rrot(i,j,l,n)*SR2(i,n,iT)
                        End Do
                     End Do
                  End Do

               End Do ! j, neq(i)
            End Do ! i, nneq
            ibuf=0
            ibuf=3*(nT+nTempMagn)*nCenter
            Call Add_Info('dM/dH    MLT0',dnrm2_(ibuf,MLT0,1),1,8)
            Call Add_Info('dM/dH    SLT0',dnrm2_(ibuf,SLT0,1),1,8)
            Call Add_Info('dM/dH    MRT0',dnrm2_(ibuf,MRT0,1),1,8)
            Call Add_Info('dM/dH    SRT0',dnrm2_(ibuf,SRT0,1),1,8)
            !
            Call Add_Info('dM/dH    MLT1',dnrm2_(ibuf,MLT1,1),1,8)
            Call Add_Info('dM/dH    SLT1',dnrm2_(ibuf,SLT1,1),1,8)
            Call Add_Info('dM/dH    MRT1',dnrm2_(ibuf,MRT1,1),1,8)
            Call Add_Info('dM/dH    SRT1',dnrm2_(ibuf,SRT1,1),1,8)
            !
            Call Add_Info('dM/dH    MLT2',dnrm2_(ibuf,MLT2,1),1,8)
            Call Add_Info('dM/dH    SLT2',dnrm2_(ibuf,SLT2,1),1,8)
            Call Add_Info('dM/dH    MRT2',dnrm2_(ibuf,MRT2,1),1,8)
            Call Add_Info('dM/dH    SRT2',dnrm2_(ibuf,SRT2,1),1,8)
c compute the total magnetizations according to the derived formulas:
            If(dbg) Write(6,*) 'check point MLT0'
            Do iT=1,nTempTotal
               If (smagn) Then
                  Call MSUM( nCenter, Sex0(1:3,iT), Zex0(iT),
     &                       SLT0(1:nCenter,1:3,iT), ZLT0(1:nCenter,iT),
     &                       SRT0(1:nCenter,1:3,iT), ZRT0(1:nCenter,iT),
     &                       iopt, ST0(1:3,iT),  ZT0(iT)  )
                  Call MSUM( nCenter, Sex1(1:3,iT), Zex1(iT),
     &                       SLT1(1:nCenter,1:3,iT), ZLT1(1:nCenter,iT),
     &                       SRT1(1:nCenter,1:3,iT), ZRT1(1:nCenter,iT),
     &                       iopt, ST1(1:3,iT),  ZT1(iT)  )
                  Call MSUM( nCenter, Sex2(1:3,iT), Zex2(iT),
     &                       SLT2(1:nCenter,1:3,iT), ZLT2(1:nCenter,iT),
     &                       SRT2(1:nCenter,1:3,iT), ZRT2(1:nCenter,iT),
     &                       iopt, ST2(1:3,iT),  ZT2(iT)  )
               End If

               If(dbg) Write(6,*) 'check point MSUM'
                  Call MSUM( nCenter, Mex0(1:3,iT), Zex0(iT),
     &                       MLT0(1:nCenter,1:3,iT), ZLT0(1:nCenter,iT),
     &                       MRT0(1:nCenter,1:3,iT), ZRT0(1:nCenter,iT),
     &                       iopt, MT0(1:3,iT),  ZT0(iT)  )
                  Call MSUM( nCenter, Mex1(1:3,iT), Zex1(iT),
     &                       MLT1(1:nCenter,1:3,iT), ZLT1(1:nCenter,iT),
     &                       MRT1(1:nCenter,1:3,iT), ZRT1(1:nCenter,iT),
     &                       iopt, MT1(1:3,iT),  ZT1(iT)  )
                  Call MSUM( nCenter, Mex2(1:3,iT), Zex2(iT),
     &                       MLT2(1:nCenter,1:3,iT), ZLT2(1:nCenter,iT),
     &                       MRT2(1:nCenter,1:3,iT), ZRT2(1:nCenter,iT),
     &                       iopt, MT2(1:3,iT),  ZT2(iT)  )
            End Do

         End If ! m_accurate

         Call Add_Info('dM/dH    ST0',ST0,3*(nT+nTempMagn),6)
         Call Add_Info('dM/dH    ST1',ST1,3*(nT+nTempMagn),6)
         Call Add_Info('dM/dH    ST2',ST2,3*(nT+nTempMagn),6)
         Call Add_Info('dM/dH    MT0',MT0,3*(nT+nTempMagn),6)
         Call Add_Info('dM/dH    MT1',MT1,3*(nT+nTempMagn),6)
         Call Add_Info('dM/dH    MT2',MT2,3*(nT+nTempMagn),6)
c  computing the AVERAGE MOMENTS calculated at dIfferent temperatures (T(i))
         Do iT=1,nTempTotal

            ! dM/dH model
            XTtens_dMdH(iM,1,iT) = (MT2(1,iT)-MT1(1,iT))*T(iT)
     &                             * cm3tomB / dltXF
            XTtens_dMdH(iM,2,iT) = (MT2(2,iT)-MT1(2,iT))*T(iT)
     &                             * cm3tomB / dltXF
            XTtens_dMdH(iM,3,iT) = (MT2(3,iT)-MT1(3,iT))*T(iT)
     &                             * cm3tomB / dltXF

            ! M/H model
            XTtens_MH(iM,1,iT) = MT0(1,iT)*T(iT)*cm3tomB/Xfield
            XTtens_MH(iM,2,iT) = MT0(2,iT)*T(iT)*cm3tomB/Xfield
            XTtens_MH(iM,3,iT) = MT0(3,iT)*T(iT)*cm3tomB/Xfield
         End Do
c ///  closing the loops over field directions
      End Do ! iM (nDirX)
c computing the XT as tensor's average:
      XTM_dMdH=0.0_wp
      XTM_MH  =0.0_wp
      Do iT=1,nTempTotal
         XTM_dMdH(iT) = ( XTtens_dMdH(1,1,iT)
     &                  + XTtens_dMdH(2,2,iT)
     &                  + XTtens_dMdH(3,3,iT) ) / 3.0_wp

         XTM_MH(iT)   = ( XTtens_MH(1,1,iT)
     &                  + XTtens_MH(2,2,iT)
     &                  + XTtens_MH(3,3,iT)   ) / 3.0_wp
      End Do !iT

      Call Add_Info('dM/dH    XTtens_dMdH',XTtens_dMdH,
     &                                            9*(nT+nTempMagn),6)
      Call Add_Info('dM/dH    XTtens_MH',XTtens_MH,
     &                                            9*(nT+nTempMagn),6)
      Call Add_Info('dM/dH    XTM_dMdH',XTM_dMdH,(nT+nTempMagn),6)
      Call Add_Info('dM/dH    XTM_MH'  ,XTM_MH  ,(nT+nTempMagn),6)
C -------------------------------------------------------------------
C   WRITING SOME OF THE OUTPUT....
C -------------------------------------------------------------------
      Write(6,*)
      Write(6,'(A)') '----------------------------------------------'//
     &               '----------------------------|'
      Write(6,'(A)') '     |     T      | Statistical |   CHI*T     '//
     &               '|   CHI*T     |   CHI*T     |'
      Write(6,'(A,F6.3,A,F6.3,A)')
     &                '     |            |  Sum (Z)    | H = 0.000 T '//
     &               '| H =',Xfield,' T | H =',XField,' T |'
      Write(6,'(A)') '     |            |             |             '//
     &               '| X = dM/dH   | X = M/H     |'
      Write(6,'(A)') '-----|----------------------------------------'//
     &               '----------------------------|'
      Write(6,'(A)') 'Units|   Kelvin   |    ---      |  cm3*K/mol  '//
     &               '|  cm3*K/mol  |  cm3*K/mol  |'
      Write(6,'(A)') '-----|----------------------------------------'//
     &               '----------------------------|'

      Do iT=1,nT
         jT=iT+nTempMagn
         Write(6,'(A,F11.6,A,E12.5,A,F12.8,A,F12.8,A,F12.7,A)')
     & '     |',T(jT),' |', ZT0(jT),' |',XT_no_field(jT),' |',
     &          XTM_dMdH(jT),' |', XTM_MH(jT), ' |'
      End Do
      Write(6,'(A)') '-----|----------------------------------------'//
     &               '----------------------------|'

c  calcualtion of the standard deviation:
      If (tinput) Then
         Write(6,'(a,5x, f20.14)') 'ST.DEV: X= dM/dH:',
     &          dev( nT, XTM_dMdH((1+nTempMagn):(nT+nTempMagn)),
     &                   chit_exp((1+nTempMagn):(nT+nTempMagn))  )
         Write(6,'(a,5x, f20.14)') 'ST.DEV: X=  M/H :',
     &          dev( nT, XTM_MH((1+nTempMagn):(nT+nTempMagn)),
     &                   chit_exp((1+nTempMagn):(nT+nTempMagn))  )
      Write(6,'(A)') '-----|----------------------------------------'//
     &               '----------------------------|'
      End If !tinput






c print out the main VAN VLECK SUSCEPTIBILITY TENSOR, its main values and main axes:
      Write(6,'(/)')
      Write(6,'(111A)') ('-',i=1,110),'|'
      Write(6,'(25X,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR '//
     &                         'FOR the  X=dM/dH  model,  '//
     &                         'in cm3*K/mol','|'
      Write(6,'(111A)') ('-',i=1,110),'|'
      Write(6,'(A)') '     T(K)   |   |          Susceptibility '//
     &               'Tensor      |    Main Values  |           '//
     &               '    Main Axes             |'
      Do iT=1,nT
        jT=iT+nTempMagn
        info=0
        wt=0.0_wp
        zt=0.0_wp
        Call DIAG_R2( XTtens_dMdH(:,:,jT) ,3,info,wt,zt)
        Write(6,'(A)') '------------|---|'//
     &                 '------- x --------- y --------- z ---|'//
     &                 '-----------------|'//
     &                 '------ x --------- y --------- z ----|'
        Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &                 '            | x |',
     &                  (XTtens_dMdH(1,j,jT),j=1,3),
     &                 ' |  X:',wt(1),'|',(zt(j,1),j=1,3),'|'
        Write(6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),
     &                             '| y |',
     &                  (XTtens_dMdH(2,j,jT),j=1,3),
     &                 ' |  Y:',wt(2),'|',(zt(j,2),j=1,3),'|'
        Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &                 '            | z |',
     &                  (XTtens_dMdH(3,j,jT),j=1,3),
     &                 ' |  Z:',wt(3),'|',(zt(j,3),j=1,3),'|'
      End Do
      Write(6,'(111A)') ('-',i=1,110),'|'




      Write(6,'(/)')
      Write(6,'(111A)') ('-',i=1,110),'|'

      Write(6,'(25X,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR '//
     &                         'FOR the  X = M/H  model,  '//
     &                         'in cm3*K/mol','|'
      Write(6,'(111A)') ('-',i=1,110),'|'
      Write(6,'(A)') '     T(K)   |   |          Susceptibility '//
     &               'Tensor      |    Main Values  |           '//
     &               '    Main Axes             |'
      Do iT=1,nT
        jT=iT+nTempMagn
        info=0
        wt=0.0_wp
        zt=0.0_wp
        Call DIAG_R2( XTtens_MH(:,:,jT) ,3,info,wt,zt)
        Write(6,'(A)') '------------|---|'//
     &                 '------- x --------- y --------- z ---|'//
     &                 '-----------------|'//
     &                 '------ x --------- y --------- z ----|'
        Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &                 '            | x |',
     &                  (XTtens_MH(1,j,jT),j=1,3),
     &                 ' |  X:',wt(1),'|',(zt(j,1),j=1,3),'|'
        Write(6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),
     &                             '| y |',
     &                  (XTtens_MH(2,j,jT),j=1,3),
     &                 ' |  Y:',wt(2),'|',(zt(j,2),j=1,3),'|'
        Write(6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)')
     &                 '            | z |',
     &                  (XTtens_MH(3,j,jT),j=1,3),
     &                 ' |  Z:',wt(3),'|',(zt(j,3),j=1,3),'|'
      End Do
      Write(6,'(111A)') ('-',i=1,110),'|'




!=======================================================================
      If(nM>0) Then
        Call mma_deallocate(WEX0)
        Call mma_deallocate(WEX1)
        Call mma_deallocate(WEX2)
      End If
      If((nneq>0).and.(nLoc>0)) Then
        Call mma_deallocate(WL0)
        Call mma_deallocate(WL1)
        Call mma_deallocate(WL2)
        Call mma_deallocate(WR0)
        Call mma_deallocate(WR1)
        Call mma_deallocate(WR2)
      End If

      If((nneq>0).and.((nT+nTempMagn)>0)) Then
        Call mma_deallocate(ZL0)
        Call mma_deallocate(ZR0)
        Call mma_deallocate(ZL1)
        Call mma_deallocate(ZR1)
        Call mma_deallocate(ZL2)
        Call mma_deallocate(ZR2)
        Call mma_deallocate(SL0)
        Call mma_deallocate(SR0)
        Call mma_deallocate(SL1)
        Call mma_deallocate(SR1)
        Call mma_deallocate(SL2)
        Call mma_deallocate(SR2)
        Call mma_deallocate(ML0)
        Call mma_deallocate(MR0)
        Call mma_deallocate(ML1)
        Call mma_deallocate(MR1)
        Call mma_deallocate(ML2)
        Call mma_deallocate(MR2)
      End If

      If((nT+nTempMagn)>0) Then
        Call mma_deallocate(ZEX0)
        Call mma_deallocate(ZEX1)
        Call mma_deallocate(ZEX2)
        Call mma_deallocate(ZT0)
        Call mma_deallocate(ZT1)
        Call mma_deallocate(ZT2)
        Call mma_deallocate(SEX0)
        Call mma_deallocate(SEX1)
        Call mma_deallocate(SEX2)
        Call mma_deallocate(MEX0)
        Call mma_deallocate(MEX1)
        Call mma_deallocate(MEX2)
        Call mma_deallocate(ST0)
        Call mma_deallocate(ST1)
        Call mma_deallocate(ST2)
        Call mma_deallocate(MT0)
        Call mma_deallocate(MT1)
        Call mma_deallocate(MT2)
        Call mma_deallocate(XTM_MH)
        Call mma_deallocate(XTM_dMdH)
        Call mma_deallocate(XTtens_MH)
        Call mma_deallocate(XTtens_dMdH)
      End If

      If((nCenter>0).and.((nT+nTempMagn)>0)) Then
        Call mma_deallocate(ZRT0)
        Call mma_deallocate(ZLT0)
        Call mma_deallocate(ZRT1)
        Call mma_deallocate(ZLT1)
        Call mma_deallocate(ZRT2)
        Call mma_deallocate(ZLT2)
        Call mma_deallocate(MRT0)
        Call mma_deallocate(MLT0)
        Call mma_deallocate(SRT0)
        Call mma_deallocate(SLT0)
        Call mma_deallocate(MRT1)
        Call mma_deallocate(MLT1)
        Call mma_deallocate(SRT1)
        Call mma_deallocate(SLT1)
        Call mma_deallocate(MRT2)
        Call mma_deallocate(MLT2)
        Call mma_deallocate(SRT2)
        Call mma_deallocate(SLT2)
      End If

      Call qExit('PA_XTdMdH')
      Return
      End



