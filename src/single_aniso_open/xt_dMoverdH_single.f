       Subroutine XT_dMoverdH_single( nss, nTempMagn, nT, nM,
     &                                Tmin, Tmax, chit_exp, eso, T,
     &                                zJ, Xfield, EM, dM, sM,
     &                                XT_no_field, tinput, smagn, mem )
c       chi*t ----------- the units are cgsemu: [ cm^3*k/mol ]
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
!#include "mgrid.fh"
#include "stdalloc.fh"
      Integer, intent(in) :: nss, nTempMagn, nT, NM, mem
      Real(kind=wp), intent(in) :: Tmin, Tmax
      Real(kind=wp), intent(in) :: zJ
      Real(kind=wp), intent(in) :: Xfield
      Real(kind=wp), intent(in) :: EM
      Real(kind=wp), intent(in) :: eso(nss)
      Real(kind=wp), intent(in) :: T(nT+nTempMagn)
      Real(kind=wp), intent(in) :: chit_exp(nT)
      Real(kind=wp), intent(in) :: XT_no_field( nT+nTempMagn )
      Complex(kind=wp), intent(in) :: dM(3,nss,nss)
      Complex(kind=wp), intent(in) :: sM(3,nss,nss)
      Logical, intent(in)          :: tinput, smagn
cccc local variables ccc
      Integer       :: iM,iT,jT,i,j,l, nDirX
      Logical       :: m_paranoid
      Real(kind=wp) :: THRS
      Real(kind=wp), allocatable :: WM1(:), WM2(:), WM3(:), WM4(:),
     &                              WM5(:), WM6(:), WM7(:)
      !WM0(nm), WM1(nm), WM2(nm) ! Zeeman exchange energies
c data for total system:
      Real(kind=wp), allocatable :: ZT1(:)   ! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ST1(:,:) ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT1(:,:) ! total magnetisation
      Real(kind=wp), allocatable :: ZT2(:)   ! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ST2(:,:) ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT2(:,:) ! total magnetisation
      Real(kind=wp), allocatable :: ZT3(:)   ! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ST3(:,:) ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT3(:,:) ! total magnetisation
      Real(kind=wp), allocatable :: ZT4(:)   ! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ST4(:,:) ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT4(:,:) ! total magnetisation
      Real(kind=wp), allocatable :: ZT5(:)   ! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ST5(:,:) ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT5(:,:) ! total magnetisation
      Real(kind=wp), allocatable :: ZT6(:)   ! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ST6(:,:) ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT6(:,:) ! total magnetisation
      Real(kind=wp), allocatable :: ZT7(:)   ! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ST7(:,:) ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT7(:,:) ! total magnetisation
c standard deviation data:
      Real(kind=wp) :: dev
      external dev
      Real(kind=wp), allocatable :: XTM_MH(:) !XTM_MH(  nT+nTempMagn)
      Real(kind=wp), allocatable :: XTM_dMdH(:) !XTM_dMdH(nT+nTempMagn)
      Real(kind=wp), allocatable :: XTtens_MH(:,:,:) !XTtens_MH(  3,3,nT+nTempMagn)
      Real(kind=wp), allocatable :: XTtens_dMdH(:,:,:) !XTtens_dMdH(3,3,nT+nTempMagn)

      Integer       :: nTempTotal
      Real(kind=wp) :: Xfield_1, Xfield_2, Xfield_3, Xfield_4, Xfield_5,
     &                 Xfield_6, Xfield_7
      Real(kind=wp) :: hp

      Real(kind=wp) :: dHX(3)
      Real(kind=wp) :: dHY(3)
      Real(kind=wp) :: dHZ(3)
      Real(kind=wp) :: dHW(3)

      Integer       :: Info
      Real(kind=wp) :: WT(3), ZT(3,3)

      Integer       :: mem_local, RtoB
      Real(kind=wp) :: cm3tomB
      logical       :: DBG

      Call qEnter('XT_dMdH')
      DBG=.false.
      m_paranoid=.true.!.false.
      cm3tomB=0.5584938904_wp   !   in cm3 * mol-1 * T
      THRS=1.d-13 !threshold for convergence of average spin, in case (zJ .ne. 0)
      RtoB=8
      mem_local=0

cccc-------------------------------------------------------cccc
      If(DBG) Then
         Write(6,'(A, 10I5)') '       nss =',nss
         Write(6,'(A,   I5)') ' nTempMagn =',nTempMagn
         Write(6,'(A,   I5)') '        nT =',nT
         Write(6,'(A,   I5)') '        nM =',nM
         Write(6,'(A,F12.5)') '      Tmin =',Tmin
         Write(6,'(A,F12.5)') '      Tmax =',Tmax
         Write(6,'(A,F12.5)') '        zJ =',zJ
         Write(6,'(A,F12.5)') '    Xfield =',Xfield
         Write(6,'(A,F12.5)') '        EM =',EM
         Write(6,'(A,E12.5)') '      THRS =',THRS
         Write(6,          *) '     smagn =',smagn
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
         Write(6,'(A)') 'SPIN-ORBIT ENERGY'
         Write(6,'(10F12.5)') (ESO(i),i=1,NSS)
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
      Write(6,'(2x,A,F10.6,A)') 'dM/dH is computed numerically '//
     &                          'using 7 point stencil formula, with '//
     &                          'h=0.00001'

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
! allocate memory:
      If ( nM > 0 ) Then
         Call mma_allocate(WM1,nM,'WM1')
         Call mma_allocate(WM2,nM,'WM2')
         Call mma_allocate(WM3,nM,'WM3')
         Call mma_allocate(WM4,nM,'WM4')
         Call mma_allocate(WM5,nM,'WM5')
         Call mma_allocate(WM6,nM,'WM6')
         Call mma_allocate(WM7,nM,'WM7')
         mem_local=mem_local+7*nM*RtoB
      End If
      If ( (nT+nTempMagn) > 0 ) Then
         Call mma_allocate(ZT1,(nT+nTempMagn),'ZT1')
         Call mma_allocate(ZT2,(nT+nTempMagn),'ZT2')
         Call mma_allocate(ZT3,(nT+nTempMagn),'ZT3')
         Call mma_allocate(ZT4,(nT+nTempMagn),'ZT4')
         Call mma_allocate(ZT5,(nT+nTempMagn),'ZT5')
         Call mma_allocate(ZT6,(nT+nTempMagn),'ZT6')
         Call mma_allocate(ZT7,(nT+nTempMagn),'ZT7')
         mem_local=mem_local+7*(nT+nTempMagn)*RtoB

         Call mma_allocate(MT1,3,(nT+nTempMagn),'MT0')
         Call mma_allocate(MT2,3,(nT+nTempMagn),'MT1')
         Call mma_allocate(MT3,3,(nT+nTempMagn),'MT3')
         Call mma_allocate(MT4,3,(nT+nTempMagn),'MT4')
         Call mma_allocate(MT5,3,(nT+nTempMagn),'MT5')
         Call mma_allocate(MT6,3,(nT+nTempMagn),'MT6')
         Call mma_allocate(MT7,3,(nT+nTempMagn),'MT7')

         Call mma_allocate(ST1,3,(nT+nTempMagn),'ST1')
         Call mma_allocate(ST2,3,(nT+nTempMagn),'ST2')
         Call mma_allocate(ST3,3,(nT+nTempMagn),'ST3')
         Call mma_allocate(ST4,3,(nT+nTempMagn),'ST4')
         Call mma_allocate(ST5,3,(nT+nTempMagn),'ST5')
         Call mma_allocate(ST6,3,(nT+nTempMagn),'ST6')
         Call mma_allocate(ST7,3,(nT+nTempMagn),'ST7')
         mem_local=mem_local+14*3*(nT+nTempMagn)*RtoB

         Call mma_allocate(XTM_MH  ,(nT+nTempMagn),'XTM_MH')
         Call mma_allocate(XTM_dMdH,(nT+nTempMagn),'XTM_dMdH')
         Call mma_allocate(XTtens_MH  ,3,3,(nT+nTempMagn),'XTtens_MH')
         Call mma_allocate(XTtens_dMdH,3,3,(nT+nTempMagn),'XTtens_dMdH')
         mem_local=mem_local+(2+2*3*3)*(nT+nTempMagn)*RtoB
      End If
      If(dbg) Write(6,*) 'XTMG:  memory allocated (local):'
      If(dbg) Write(6,*) 'mem_local=', mem_local
      If(dbg) Write(6,*) 'XTMG:  memory allocated (total):'
      If(dbg) Write(6,*) 'mem_total=', mem+mem_local





      dHX=0.0_wp
      dHY=0.0_wp
      dHZ=0.0_wp
      dHW=0.0_wp

      nDirX=3

      ! field along X
      dHX(1)=1.0_wp
      dHY(1)=0.0_wp
      dHZ(1)=0.0_wp
      ! field along Y
      dHX(2)=0.0_wp
      dHY(2)=1.0_wp
      dHZ(2)=0.0_wp
      ! field along Z
      dHX(3)=0.0_wp
      dHY(3)=0.0_wp
      dHZ(3)=1.0_wp

      Call dcopy_(3*3*(nT+nTempMagn),0.0_wp,0,XTtens_MH,1)
      Call dcopy_(3*3*(nT+nTempMagn),0.0_wp,0,XTtens_dMdH,1)

      hp=0.0001_wp
      Xfield_1=Xfield-3.0_wp*hp
      Xfield_2=Xfield-2.0_wp*hp
      Xfield_3=Xfield-1.0_wp*hp
      Xfield_4=xField
      Xfield_5=Xfield+1.0_wp*hp
      Xfield_6=Xfield+2.0_wp*hp
      Xfield_7=Xfield+3.0_wp*hp
      If(dbg) Write(6,*) 'XTMG:  Xfield: ', Xfield_1, Xfield_2,
     &        Xfield_3, Xfield_4, Xfield_5, Xfield_6, Xfield_7
     &

      nTempTotal=nT+nTempMagn
      m_paranoid = .true.


c ///  opening the loop over different directions of the magnetic field
      Do iM=1,nDirX
         Call dcopy_(nM,0.0_wp,0,WM1,1)
         Call dcopy_(nM,0.0_wp,0,WM2,1)
         Call dcopy_(nM,0.0_wp,0,WM3,1)
         Call dcopy_(nM,0.0_wp,0,WM4,1)
         Call dcopy_(nM,0.0_wp,0,WM5,1)
         Call dcopy_(nM,0.0_wp,0,WM6,1)
         Call dcopy_(nM,0.0_wp,0,WM7,1)

         Call dcopy_(  (nT+nTempMagn),0.0_wp,0,ZT1,1)
         Call dcopy_(  (nT+nTempMagn),0.0_wp,0,ZT2,1)
         Call dcopy_(  (nT+nTempMagn),0.0_wp,0,ZT3,1)
         Call dcopy_(  (nT+nTempMagn),0.0_wp,0,ZT4,1)
         Call dcopy_(  (nT+nTempMagn),0.0_wp,0,ZT5,1)
         Call dcopy_(  (nT+nTempMagn),0.0_wp,0,ZT6,1)
         Call dcopy_(  (nT+nTempMagn),0.0_wp,0,ZT7,1)

         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MT1,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MT2,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MT3,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MT4,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MT5,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MT6,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,MT7,1)

         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,ST1,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,ST2,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,ST3,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,ST4,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,ST5,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,ST6,1)
         Call dcopy_(3*(nT+nTempMagn),0.0_wp,0,ST7,1)

        ! compute magnetization:
        ! seven points numerical field perturbation:
           !  Xfield =XF -3h
         Call MAGN( NSS, NM, dHX(iM), dHY(iM), dHZ(iM), XField_1,
     &              ESO, zJ, THRS,  DM, SM, nTempTotal, T, smagn,
     &              WM1, ZT1, ST1, MT1, m_paranoid, DBG )
           !  Xfield =XF -2h
         Call MAGN( NSS, NM, dHX(iM), dHY(iM), dHZ(iM), XField_2,
     &              ESO, zJ, THRS,  DM, SM, nTempTotal, T, smagn,
     &              WM2, ZT2, ST2, MT2, m_paranoid, DBG )
           !  Xfield =XF -h
         Call MAGN( NSS, NM, dHX(iM), dHY(iM), dHZ(iM), XField_3,
     &              ESO, zJ, THRS,  DM, SM, nTempTotal, T, smagn,
     &              WM3, ZT3, ST3, MT3, m_paranoid, DBG )
           !  Xfield =XF
         Call MAGN( NSS, NM, dHX(iM), dHY(iM), dHZ(iM), XField_4,
     &              ESO, zJ, THRS,  DM, SM, nTempTotal, T, smagn,
     &              WM4, ZT4, ST4, MT4, m_paranoid, DBG )
           !  Xfield =XF +h
         Call MAGN( NSS, NM, dHX(iM), dHY(iM), dHZ(iM), XField_5,
     &              ESO, zJ, THRS,  DM, SM, nTempTotal, T, smagn,
     &              WM5, ZT5, ST5, MT5, m_paranoid, DBG )
           !  Xfield =XF +2h
         Call MAGN( NSS, NM, dHX(iM), dHY(iM), dHZ(iM), XField_6,
     &              ESO, zJ, THRS,  DM, SM, nTempTotal, T, smagn,
     &              WM6, ZT6, ST6, MT6, m_paranoid, DBG )
           !  Xfield =XF +3h
         Call MAGN( NSS, NM, dHX(iM), dHY(iM), dHZ(iM), XField_7,
     &              ESO, zJ, THRS,  DM, SM, nTempTotal, T, smagn,
     &              WM7, ZT7, ST7, MT7, m_paranoid, DBG )

      If(DBG) Then
        Do iT=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex1(',iT,')=',(MT1(l,iT),l=1,3)
        End Do
        Do iT=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex2(',iT,')=',(MT2(l,iT),l=1,3)
        End Do
        Do iT=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex3(',iT,')=',(MT3(l,iT),l=1,3)
        End Do
        Do iT=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex4(',iT,')=',(MT4(l,iT),l=1,3)
        End Do
        Do iT=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex5(',iT,')=',(MT5(l,iT),l=1,3)
        End Do
        Do iT=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex6(',iT,')=',(MT6(l,iT),l=1,3)
        End Do
        Do iT=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex7(',iT,')=',(MT7(l,iT),l=1,3)
        End Do
        Do iT=1,nTempTotal
          Write(6,'(A,i3,A,9E22.14)') 'Mex 7-1 (',iT,')=',
     &                                (MT7(l,iT)-MT1(l,iT),l=1,3)
        End Do
      End If !DBG
         ! computing the AVERAGE MOMENTS calculated at different temperatures (T(i))
         Do iT=1,nTempTotal

            ! dM/dH model
            If(iM==1) Then
            XTtens_dMdH(iM,1,iT) = (-1.0_wp*MT1(1,iT)
     &                              +9.0_wp*MT2(1,iT)
     &                             -45.0_wp*MT3(1,iT)
     &                              +0.0_wp*MT4(1,iT)
     &                             +45.0_wp*MT5(1,iT)
     &                              -9.0_wp*MT6(1,iT)
     &                              +1.0_wp*MT7(1,iT) ) *T(iT)* cm3tomB
     &                             / (60.0_wp* hp)

            Else If (iM==2) Then
            XTtens_dMdH(iM,2,iT) = (-1.0_wp*MT1(2,iT)
     &                              +9.0_wp*MT2(2,iT)
     &                             -45.0_wp*MT3(2,iT)
     &                              +0.0_wp*MT4(2,iT)
     &                             +45.0_wp*MT5(2,iT)
     &                              -9.0_wp*MT6(2,iT)
     &                              +1.0_wp*MT7(2,iT) ) *T(iT)* cm3tomB
     &                             / (60.0_wp* hp)
            Else If (iM==3) Then
            XTtens_dMdH(iM,3,iT) = (-1.0_wp*MT1(3,iT)
     &                              +9.0_wp*MT2(3,iT)
     &                             -45.0_wp*MT3(3,iT)
     &                              +0.0_wp*MT4(3,iT)
     &                             +45.0_wp*MT5(3,iT)
     &                              -9.0_wp*MT6(3,iT)
     &                              +1.0_wp*MT7(3,iT) ) *T(iT)* cm3tomB
     &                             / (60.0_wp* hp)
            End If

            ! M/H model
            If(iM==1) Then
            XTtens_MH(iM,1,iT) = MT4(1,iT)*T(iT)*cm3tomB/Xfield
            Else If (iM==2) Then
            XTtens_MH(iM,2,iT) = MT4(2,iT)*T(iT)*cm3tomB/Xfield
            Else If (iM==3) Then
            XTtens_MH(iM,3,iT) = MT4(3,iT)*T(iT)*cm3tomB/Xfield
            End If
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
c      Call Add_Info('XTM_dMdH          ',XTM_dMdH,nTempTotal,5)
c      Call Add_Info('XTM_MH            ',XTM_MH  ,nTempTotal,5)

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
     & '     |',T(jT),' |', ZT3(jT),' |',XT_no_field(jT),' |',
     &          XTM_dMdH(jT),' |', XTM_MH(jT), ' |'
      End Do
      Write(6,'(A)') '-----|----------------------------------------'//
     &               '----------------------------|'

c  calcualtion of the standard deviation:
      If (tinput) Then
         Write(6,'(a,5x, f20.14)') 'ST.DEV: X= dM/dH:',
     &          dev(  (nT-nTempMagn),
     &                     XTM_dMdH((1+nTempMagn):(nT+nTempMagn)),
     &                     chit_exp((1+nTempMagn):(nT+nTempMagn))  )
         Write(6,'(a,5x, f20.14)') 'ST.DEV: X= M/H:',
     &          dev(  (nT-nTempMagn),
     &                     XTM_MH((1+nTempMagn):(nT+nTempMagn)),
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


      If ( nM > 0 ) Then
         Call mma_deallocate(WM1)
         Call mma_deallocate(WM2)
         Call mma_deallocate(WM3)
         Call mma_deallocate(WM4)
         Call mma_deallocate(WM5)
         Call mma_deallocate(WM6)
         Call mma_deallocate(WM7)
      End If
      If ( nT+nTempMagn > 0 ) Then
         Call mma_deallocate(ZT1)
         Call mma_deallocate(ZT2)
         Call mma_deallocate(ZT3)
         Call mma_deallocate(ZT4)
         Call mma_deallocate(ZT5)
         Call mma_deallocate(ZT6)
         Call mma_deallocate(ZT7)

         Call mma_deallocate(MT1)
         Call mma_deallocate(MT2)
         Call mma_deallocate(MT3)
         Call mma_deallocate(MT4)
         Call mma_deallocate(MT5)
         Call mma_deallocate(MT6)
         Call mma_deallocate(MT7)

         Call mma_deallocate(ST1)
         Call mma_deallocate(ST2)
         Call mma_deallocate(ST3)
         Call mma_deallocate(ST4)
         Call mma_deallocate(ST5)
         Call mma_deallocate(ST6)
         Call mma_deallocate(ST7)

         Call mma_deallocate(XTM_MH)
         Call mma_deallocate(XTM_dMdH)
         Call mma_deallocate(XTtens_MH)
         Call mma_deallocate(XTtens_dMdH)
      End If
      Call qExit('XT_dMdH')

      Return
      End



