!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine torque_pa(nneq,nCenter,neq,neqv,nLoc,exch,nTempMagn,nH,nM,AngPoints,nexch,iopt,nss,mem,smagn,m_paranoid,m_accurate, &
                     TempMagn,w,hmin,hmax,dltH0,EM,zJ,THRS,hexp,dipexch,s_exch,dipso,s_so,eso,hinput,r_rot,XLM,ZLM,XRM,ZRM)
! correction to M from the local excited states:
!  XLM, ZLM, XRM, ZRM
! rotation matrices for equivalent sites:
!  R_ROT
! exchange spectum:
! exchange energies printed out in the previous part
!  W
! local spin-orbit spectum:
! spin-orbit energies from ANISO files
!  ESO
! magnetic and spin moments (i.e. the BIG matrices):
!  DIPEXCH, S_EXCH, dipso, s_so
! exchange data:
!  WEX : Zeeman exchange energies
!  ZEX : exchange statistical sum, Boltzmann distribution
!  SEX : spin magnetisation, from the exchange block
!  MEX : magnetisation, form the exchange block
! data for individual sites (all states):
!  ZL : local statistical sum, Boltzmann distribution
!  WL : Zeeman local energies
!  SL : spin magnetisation, from the local sites, using ALL states
!  ML : magnetisation, from local sites, using ALL states
! data for individual sites (only states that enter exchange):
!  ZR : local statistical sum, Boltzmann distribution, using only NEXCH states
!  WR : Zeeman local reduced energies, using only NEXCH states
!  SR : spin magnetisation, from the local sites, using only NEXCH states
!  MR : magnetisation, from local sites, using only NEXCH states
! total vectors in general coordinate system:
!  ZRT, ZLT, MRT, MLT, SRT, SLT
! data for total system:
!  ZT : total statistical sum, Boltzmann distribution
!  ST : total spin magnetisation,
!  MT : total magnetisation
! magnetic field strength and orientation data:
!  nPlanes, dlth, H, dX, dY, dZ, Ang
! magnetic torque
!  tx : magnetization torque, X
!  ty : magnetization torque, Y
!  tz : magnetization torque, Z
!  sx : spin magnetization torque, X
!  sy : spin magnetization torque, Y
!  sz : spin magnetization torque, Z

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Ten, mBohr, rNAVO
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: nneq, nCenter, neq(nneq), neqv, nLoc, exch, nTempMagn, nH, nM, AngPoints, nexch(nneq), iopt, &
                                 nss(nneq), mem
logical(kind=iwp), intent(in) :: smagn, m_paranoid, m_accurate, hinput
real(kind=wp), intent(in) :: TempMagn(nTempMagn), W(exch), hmin, hmax, dltH0, EM, zJ, THRS, Hexp(nH), ESO(nneq,nLoc), &
                             R_ROT(nneq,neqv,3,3), XLM(nCenter,nTempMagn,3,3), ZLM(nCenter,nTempMagn), XRM(nCenter,nTempMagn,3,3), &
                             ZRM(nCenter,nTempMagn)
complex(kind=wp), intent(in) :: DIPEXCH(3,EXCH,EXCH), S_EXCH(3,EXCH,EXCH), dipso(nneq,3,nLoc,nLoc), s_so(nneq,3,nLoc,nLoc)
integer(kind=iwp) :: I, IH, IM, iPl, isite, it, J, k, mem_local, n
real(kind=wp) :: dlth
real(kind=wp), allocatable :: Ang(:), dX(:,:), dY(:,:), dZ(:,:), ESO_TMP(:), H(:), MEX(:,:), ML(:,:,:), MLT(:,:,:), MR(:,:,:), &
                              MRT(:,:,:), MT(:,:), SEX(:,:), SL(:,:,:), SLT(:,:,:), SR(:,:,:), SRT(:,:,:), ST(:,:), sx(:,:,:,:), &
                              sy(:,:,:,:), sz(:,:,:,:), tx(:,:,:,:), ty(:,:,:,:), tz(:,:,:,:), WEX(:), WL(:,:), WR(:,:), ZEX(:), &
                              ZL(:,:), ZLT(:,:), ZR(:,:), ZRT(:,:), ZT(:)
complex(kind=wp), allocatable :: dipso_tmp(:,:,:), dipso_tmp2(:,:,:), s_so_tmp(:,:,:), s_so_tmp2(:,:,:)
integer(kind=iwp), parameter :: nPlanes = 3
real(kind=wp), parameter :: cm3tomB = rNAVO*mBohr/Ten ! in cm3 * mol-1 * T
#ifdef _DEBUGPRINT_
#  define _DBG_ .true.
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: DBG = _DBG_

#ifndef _DEBUGPRINT_
#include "macros.fh"
unused_var(mem)
#endif

write(u6,*)
write(u6,'(A)') repeat('%',96)
write(u6,'(20X,A)') 'ANGULAR DEPENDENCE OF THE MAGNETIZATION TORQUE'
write(u6,'(A)') repeat('%',96)
write(u6,*)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
write(u6,'(2X,A,i3,A)') 'Magnetization torque is calculated for the ',NH,' field points, in the field Domain:'
write(u6,'(2X,F4.1,1x,a,1X,F4.1,a,30(F6.3,a))') HMIN,'--',HMAX,' T., at the following temperatures:'
do i=11,nTempMagn,10
  j = min(nTempMagn,i+9)
  write(u6,'(17x,10(F9.3,A))') (TempMagn(k),' K.;',k=i,j)
end do
write(u6,'(2x,A,i3,A)') 'Angular dependence of the magnetization torque is computed for ',AngPoints,' angular points distributed'
write(u6,'(2x,A)') 'in the domain 0-180 deg.'
write(u6,'(2X,A)') 'The cut-off energy for the exact diagonalization of the Zeeman Hamiltonian is:'
write(u6,'(2x,a,F15.9,A)') 'E = ',EM,' cm(-1).'
if (NM < 10) then
  write(u6,'(2X,A,i2,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 10) .and. (NM < 100)) then
  write(u6,'(2X,A,i3,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 100) .and. (NM < 1000)) then
  write(u6,'(2X,A,i4,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 1000) .and. (NM < 10000)) then
  write(u6,'(2X,A,i5,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
end if
if (m_accurate) then
  write(u6,'(2x,A)') 'The contribution of local excited states is computed exactly.'
else
  write(u6,'(2x,A)') 'The contribution of local excited states is computed approximately (by using the susceptibility data). '
end if
!-----------------------------------------------------------------------
#ifdef _DEBUGPRINT_
write(u6,*) 'nM       = ',nM
write(u6,*) 'nTempMagn= ',nTempMagn
write(u6,*) 'nneq     = ',nneq
write(u6,*) 'nCenter  = ',nCenter
write(u6,*) 'AngPoints= ',AngPoints
write(u6,*) 'nPlanes  = ',nPlanes
write(u6,*) 'nH       = ',nH
write(u6,*) 'nLoc     = ',nLoc
write(u6,*) 'nexch()  = ',(nexch(i),i=1,nneq)
write(u6,*) 'neq()    = ',(neq(i),i=1,nneq)
write(u6,*) 'exch     = ',exch
write(u6,*) 'iopt     = ',iopt
write(u6,*) 'neqv     = ',neqv
write(u6,*) 'nss()    = ',(nss(i),i=1,nneq)
write(u6,*) 'W()      = ',(W(i),i=1,exch)
write(u6,*) 'zJ       = ',zJ
write(u6,*) 'EM       = ',EM
write(u6,*) 'm_paranoi= ',m_paranoid
write(u6,*) 'm_accurat= ',m_accurate
write(u6,*) 'smagn    = ',smagn
write(u6,*) 'hinput   = ',hinput
#endif

! Allocate memory for this calculation:
mem_local = 0
! Zeeman exchange energy spectrum
call mma_allocate(Wex,nM,'Wex')
mem_local = mem_local+size(Wex)*RtoB

#ifdef _DEBUGPRINT_
write(u6,*) 'mem_local 1 = ',mem_local
#endif

! exchange statistical sum, Boltzmann distribution
call mma_allocate(Zex,nTempMagn,'Zex')
mem_local = mem_local+size(Zex)*RtoB
! spin magnetisation, from the exchange block
call mma_allocate(SEX,3,nTempMagn,'SEX')
mem_local = mem_local+size(SEX)*RtoB
! magnetisation, from the exchange block
call mma_allocate(MEX,3,nTempMagn,'MEX')
mem_local = mem_local+size(MEX)*RtoB

! total statistical sum, Boltzmann distribution
call mma_allocate(ZT,nTempMagn,'ZT')
mem_local = mem_local+size(ZT)*RtoB
! total spin magnetisation
call mma_allocate(ST,3,nTempMagn,'ST')
mem_local = mem_local+size(ST)*RtoB
! total magnetisation
call mma_allocate(MT,3,nTempMagn,'MT')
mem_local = mem_local+size(MT)*RtoB

#ifdef _DEBUGPRINT_
write(u6,*) 'mem_local 2 = ',mem_local
#endif
! local statistical sum, Boltzmann distribution
call mma_allocate(ZL,nTempMagn,nneq,'ZL')
mem_local = mem_local+size(ZL)*RtoB
! spin magnetisation, from the local sites, using ALL states
call mma_allocate(SL,3,nTempMagn,nneq,'SL')
mem_local = mem_local+size(SL)*RtoB
! magnetisation, from local sites, using ALL states
call mma_allocate(ML,3,nTempMagn,nneq,'ML')
mem_local = mem_local+size(ML)*RtoB

! local statistical sum, Boltzmann distribution
call mma_allocate(ZR,nTempMagn,nneq,'ZR')
mem_local = mem_local+size(ZR)*RtoB
! spin magnetisation, from the local sites, using only NEXCH states
call mma_allocate(SR,3,nTempMagn,nneq,'SR')
mem_local = mem_local+size(SR)*RtoB
! magnetisation, from the local sites, using only NEXCH states
call mma_allocate(MR,3,nTempMagn,nneq,'MR')
mem_local = mem_local+size(MR)*RtoB

#ifdef _DEBUGPRINT_
write(u6,*) 'mem_local 3 = ',mem_local
#endif
! ZRT(nCenter,nTempMagn)
call mma_allocate(ZRT,nCenter,nTempMagn,'ZRT')
mem_local = mem_local+size(ZT)*RtoB
! ZLT(nCenter,nTempMagn)
call mma_allocate(ZLT,nCenter,nTempMagn,'ZLT')
mem_local = mem_local+size(ZLT)*RtoB
! MRT(nCenter,3,nTempMagn)
call mma_allocate(MRT,nCenter,3,nTempMagn,'MRT')
mem_local = mem_local+size(MRT)*RtoB
! MLT(nCenter,3,nTempMagn)
call mma_allocate(MLT,nCenter,3,nTempMagn,'MLT')
mem_local = mem_local+size(MLT)*RtoB
! SRT(nCenter,3,nTempMagn)
call mma_allocate(SRT,nCenter,3,nTempMagn,'SRT')
mem_local = mem_local+size(SRT)*RtoB
! SLT(nCenter,3,nTempMagn)
call mma_allocate(SLT,nCenter,3,nTempMagn,'SLT')
mem_local = mem_local+size(SLT)*RtoB

#ifdef _DEBUGPRINT_
write(u6,*) 'mem_local 4 = ',mem_local
#endif
! magnetic torque
call mma_allocate(tx,nPlanes,AngPoints,nH,nTempMagn,'tx')
call mma_allocate(ty,nPlanes,AngPoints,nH,nTempMagn,'ty')
call mma_allocate(tz,nPlanes,AngPoints,nH,nTempMagn,'tz')
call mma_allocate(sx,nPlanes,AngPoints,nH,nTempMagn,'sx')
call mma_allocate(sy,nPlanes,AngPoints,nH,nTempMagn,'sy')
call mma_allocate(sz,nPlanes,AngPoints,nH,nTempMagn,'sz')
mem_local = mem_local+size(tx)*RtoB
mem_local = mem_local+size(ty)*RtoB
mem_local = mem_local+size(tz)*RtoB
mem_local = mem_local+size(sx)*RtoB
mem_local = mem_local+size(sy)*RtoB
mem_local = mem_local+size(sz)*RtoB
#ifdef _DEBUGPRINT_
write(u6,*) 'mem_local 5 = ',mem_local
#endif

! Zeeman local energies
call mma_allocate(WL,nLoc,nneq,'WL')
mem_local = mem_local+size(WL)*RtoB
! Zeeman local reduced energies, using only NEXCH states
call mma_allocate(WR,nLoc,nneq,'WR')
mem_local = mem_local+size(WR)*RtoB
#ifdef _DEBUGPRINT_
write(u6,*) 'mem_local 6 = ',mem_local
#endif

call mma_allocate(Ang,AngPoints,'Ang')
mem_local = mem_local+size(Ang)*RtoB
call mma_allocate(dX,AngPoints,nPlanes,'dX')
call mma_allocate(dY,AngPoints,nPlanes,'dY')
call mma_allocate(dZ,AngPoints,nPlanes,'dZ')
mem_local = mem_local+size(dX)*RtoB
mem_local = mem_local+size(dY)*RtoB
mem_local = mem_local+size(dZ)*RtoB
#ifdef _DEBUGPRINT_
write(u6,*) 'mem_local 7 = ',mem_local
#endif

call mma_allocate(H,nH,'H')
mem_local = mem_local+size(H)*RtoB

#ifdef _DEBUGPRINT_
write(u6,*) 'TORQ:  memory allocated (local):'
write(u6,*) 'mem_local=',mem_local
write(u6,*) 'TORQ:  memory allocated (total):'
write(u6,*) 'mem_total=',mem+mem_local
#endif

!-----------------------------------------------------------------------
! set up the field points:
if (HINPUT) then
  do iH=1,nH
    H(iH) = HEXP(iH)
    if (H(iH) == Zero) H(iH) = 0.0001_wp
  end do
else
  DLTH = (HMAX-HMIN)/real(NH-1,kind=wp)
  do iH=1,nH
    if (iH == 1) then
      H(IH) = HMIN+dltH0
    else
      H(IH) = HMIN+DLTH*real(IH-1,kind=wp)
    end if
    if (H(iH) == Zero) H(iH) = 0.0001_wp
  end do
end if
!-----------------------------------------------------------------------

do IH=1,NH
  ! ///  opening the loop over field points:
  do iPl=1,nPlanes
    ! ///  opening the loop over different planes of rotation of the applied magnetic field:
    !  loop over various angular grids  (iPl = 1, 2 or 3)
    !  iPl=1 (X) , angular grid in the plane YZ ( half of the plane , 0-180 degrees)
    !  iPl=2 (Y) , angular grid in the plane XZ ( half of the plane , 0-180 degrees)
    !  iPl=3 (Z) , angular grid in the plane XY ( half of the plane , 0-180 degrees)

    call hdir2(AngPoints,iPl,dX(:,iPl),dY(:,iPl),dZ(:,iPl),Ang,4)

#   ifdef _DEBUGPRINT_
    write(u6,*) 'iPl, dX, dY, dZ=',iPl,dX(:,iPl),dY(:,iPl),dZ(:,iPl),Ang
#   endif
    do IM=1,AngPoints
      ! ///  opening the loop over different directions of the magnetic field

      ! exchange magnetization:
      call MAGN(EXCH,NM,dX(iM,iPl),dY(iM,iPl),dZ(iM,iPl),H(iH),W,zJ,THRS,DIPEXCH,S_EXCH,nTempMagn,TempMagn,smagn,Wex,Zex,Sex,Mex, &
                m_paranoid,DBG)

#     ifdef _DEBUGPRINT_
      if (iPl == 2) &
        write(u6,'(A,I3,1x,F8.4,2x, 3F19.14,2x,3F19.14)') 'MEX: iM,',iM,H(iH),(Mex(i,1),i=1,3),dX(iM,iPl),dY(iM,iPl),dZ(iM,iPl)
#     endif
      !iT = 1

      !write(u6,'(F10.4,1x,2I3,3F18.14,3x,3F20.14,2x,F20.14)') H(iH),iL,iM,dX(iM,iPl),dY(iM,iPl),dZ(iM,iPl),(Mex(j,iT),j=1,3), &
      !                                                        Mex(1,iT)*dX(iM,iPl)+Mex(2,iT)*dY(iM,iPl)+Mex(3,iT)*dZ(iM,iPl)
      ! compute local magnetizations:
      if (m_accurate) then
        call mma_allocate(ESO_TMP,NSS(i),label='ESO_TMP')
        call mma_allocate(dipso_tmp,3,nLoc,nLoc,label='dipso_tmp')
        call mma_allocate(s_so_tmp,3,nLoc,nLoc,label='s_so_tmp')
        do i=1,nneq
          ! all states:
          if (NSS(i) > NEXCH(i)) then
            ! this check is meant to avoid the unnecessary
            ! computation, in cases when no local excited
            ! states are present
            ESO_TMP(:) = ESO(i,1:NSS(i))
            dipso_tmp(:,:,:) = DIPSO(i,:,:,:)
            s_so_tmp(:,:,:) = S_SO(i,:,:,:)
            call MAGN(NSS(i),NEXCH(i),dX(iM,iPl),dY(iM,iPl),dZ(iM,iPl),H(iH),ESO_TMP,zJ,THRS,DIPSO_TMP,S_SO_TMP,nTempMagn, &
                      TempMagn,smagn,WL(:,i),ZL(:,i),SL(:,:,i),ML(:,:,i),m_paranoid,DBG)
            ! only local "exchange states":
            call mma_allocate(dipso_tmp2,3,NEXCH(i),NEXCH(i),label='dipso_tmp2')
            call mma_allocate(s_so_tmp2,3,NEXCH(i),NEXCH(i),label='s_so_tmp2')
            dipso_tmp2(:,:,:) = DIPSO(i,:,1:NEXCH(i),1:NEXCH(i))
            s_so_tmp2(:,:,:) = S_SO(i,:,1:NEXCH(i),1:NEXCH(i))
            call MAGN(NEXCH(i),NEXCH(i),dX(iM,iPl),dY(iM,iPl),dZ(iM,iPl),H(iH),ESO_TMP(1:NEXCH(i)),zJ,THRS,DIPSO_TMP,S_SO_TMP, &
                      nTempMagn,TempMagn,smagn,WR(1:Nexch(i),i),ZR(:,i),SR(:,:,i),MR(:,:,i),m_paranoid,DBG)
            call mma_deallocate(dipso_tmp2)
            call mma_deallocate(s_so_tmp2)
          end if
        end do
        call mma_deallocate(dipso_tmp)
        call mma_deallocate(s_so_tmp)
        ! expand the basis and rotate local vectors to the
        ! general coordinate system:
        MRT(:,:,:) = Zero
        MLT(:,:,:) = Zero
        SRT(:,:,:) = Zero
        SLT(:,:,:) = Zero

        isite = 0
        do i=1,NNEQ
          do j=1,NEQ(i)
            isite = isite+1
            ! statistical distributions
            ZLT(isite,:) = ZL(:,i)
            ZRT(isite,:) = ZR(:,i)
            ! magnetizations:
            ! use R_rot matrices, which have determinant +1.
            !  >> note that  R_lg matrices may have arbitrary
            !  >> sign of the determinant.
            do iT=1,nTempMagn
              do n=1,3
                MLT(isite,:,iT) = MLT(isite,:,iT)+r_rot(i,j,:,n)*ML(n,iT,i)

                SLT(isite,:,iT) = SLT(isite,:,iT)+r_rot(i,j,:,n)*SL(n,iT,i)

                MRT(isite,:,iT) = MRT(isite,:,iT)+r_rot(i,j,:,n)*MR(n,iT,i)

                SRT(isite,:,iT) = SRT(isite,:,iT)+r_rot(i,j,:,n)*SR(n,iT,i)
              end do
            end do

          end do ! j, neq(i)
        end do ! i, nneq
      end if ! m_accurate

      ! compute the total magnetizations according
      ! to the derived formulas:
      if (m_accurate) then
        do iT=1,nTempMagn
          if (smagn) call MSUM(nCenter,Sex(:,iT),Zex(iT),SLT(:,:,iT),ZLT(:,iT),SRT(:,:,iT),ZRT(:,iT),iopt,ST(:,iT),ZT(iT))
          call MSUM(nCenter,Mex(:,iT),Zex(iT),MLT(:,:,iT),ZLT(:,iT),MRT(:,:,iT),ZRT(:,iT),iopt,MT(:,iT),ZT(iT))
        end do

      else ! (m_accurate)

        ! add the contribution from local excited states
        ! using the approximate X*H expression:
        do iT=1,nTempMagn
          do isite=1,nCenter

            MRT(isite,:,iT) = H(iH)*(XRM(isite,iT,:,1)*dX(iM,iPl)+XRM(isite,iT,:,2)*dY(iM,iPl)+XRM(isite,iT,:,3)*dZ(iM,iPl))/cm3tomB

            MLT(isite,:,iT) = H(iH)*(XLM(isite,iT,:,1)*dX(iM,iPl)+XLM(isite,iT,:,2)*dY(iM,iPl)+XLM(isite,iT,:,3)*dZ(iM,iPl))/cm3tomB

          end do ! isite

          if (smagn) call MSUM(nCenter,Sex(:,iT),Zex(iT),SLT(:,:,iT),ZLM(:,iT),SRT(:,:,iT),ZRM(:,iT),iopt,ST(:,iT),ZT(iT))
          call MSUM(nCenter,Mex(:,iT),Zex(iT),MLT(:,:,iT),ZLM(:,iT),MRT(:,:,iT),ZRM(:,iT),iopt,MT(:,iT),ZT(iT))
        end do ! iT
      end if ! (m_accurate)

      ! at this point we have MT and ST computed
      ! in the direction of applied field
      ! compute the M and S torque
      tx(iPl,iM,iH,:) = MT(2,:)*dZ(iM,iPl)*H(iH)-MT(3,:)*dY(iM,iPl)*H(iH)

      ty(iPl,iM,iH,:) = MT(3,:)*dX(iM,iPl)*H(iH)-MT(1,:)*dZ(iM,iPl)*H(iH)

      tz(iPl,iM,iH,:) = MT(1,:)*dY(iM,iPl)*H(iH)-MT(2,:)*dX(iM,iPl)*H(iH)

      if (smagn) then
        sx(iPl,iM,iH,:) = ST(2,:)*dZ(iM,iPl)*H(iH)-ST(3,:)*dY(iM,iPl)*H(iH)

        sy(iPl,iM,iH,:) = ST(3,:)*dX(iM,iPl)*H(iH)-ST(1,:)*dZ(iM,iPl)*H(iH)

        sz(iPl,iM,iH,:) = ST(1,:)*dY(iM,iPl)*H(iH)-ST(2,:)*dX(iM,iPl)*H(iH)
      end if

      ! ///  closing the loops over field strengths and directions
    end do ! iL
  end do ! iM
end do ! iH
! ----------------------------------------------------------------------
! WRITING SOME OF THE OUTPUT....
! ----------------------------------------------------------------------
write(u6,*)
write(u6,'(25X,A)') 'ANGULAR DEPENDENCE OF THE MAGNETIZATION TORQUE'
write(u6,'(30X,A)') '(Units of torque: [energy, cm-1])'
write(u6,*)

write(u6,'(5x,A)') 'Orientation of the applied magnetic field employed:'
write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

write(u6,'(2x,A,3(10x,A))') 'Angle |','rotation in the YZ plane          |','rotation in the XZ plane          |', &
                            'rotation in the XY plane          |'
write(u6,'(10A)') '--------|',('--- proj X ---|','--- proj Y ---|','--- proj Z ---|',i=1,3)
do iM=1,AngPoints
  write(u6,'(F7.3,1x,A,3(F13.10,1x,A,F13.10,1x,A,F13.10,1x,A))') Ang(iM),'|',(dX(iM,iPl),' ',dY(iM,iPl),' ',dZ(iM,iPl),'|',iPl=1,3)
end do
write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

do iH=1,nH
  do iT=1,nTempMagn
    write(u6,*)
    write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)
    write(u6,'(A,F9.4,A)') 'Magnetic field strength = ',H(iH),' tesla'
    write(u6,'(12x,A,F9.4,A)') 'Temperature = ',TempMagn(iT),' kelvin'
    write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

    write(u6,'(2x,A,3(10x,A))') 'Angle |','rotation in the YZ plane          |','rotation in the XZ plane          |', &
                                'rotation in the XY plane          |'
    write(u6,'(10A)') '--------|',('-- torque X --|','-- torque Y --|','-- torque Z --|',i=1,3)
    do iM=1,AngPoints
      write(u6,'(F7.3,1x,A,3(ES13.6,1x,A,ES13.6,1x,A,ES13.6,1x,A))') Ang(iM),'|', &
                                                                     (tx(iPl,iM,iH,iT),' ',ty(iPl,iM,iH,iT),' ',tz(iPl,iM,iH,iT), &
                                                                      '|',iPl=1,3)
    end do
    write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)
  end do !iT
end do !iH

! ----------------------------------------------------------------------
if (smagn) then
  write(u6,*)
  write(u6,'(25X,A)') 'ANGULAR DEPENDENCE OF THE SPIN MAGNETIZATION TORQUE'
  write(u6,'(30X,A)') '(Units of torque: [energy, cm-1])'
  write(u6,*)

  write(u6,'(5x,A)') 'Orientation of the applied magnetic field employed:'
  write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

  write(u6,'(2x,A,3(10x,A))') 'Angle |','rotation in the YZ plane          |','rotation in the XZ plane          |', &
                              'rotation in the XY plane          |'
  write(u6,'(10A)') '--------|',('--- proj X ---|','--- proj Y ---|','--- proj Z ---|',i=1,3)
  do iM=1,AngPoints
    write(u6,'(F7.3,1x,A,3(F13.10,1x,A,F13.10,1x,A,F13.10,1x,A))') Ang(iM),'|', &
                                                                   (dX(iM,iPl),' ',dY(iM,iPl),' ',dZ(iM,iPl),'|',iPl=1,3)
  end do
  write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

  do iH=1,nH
    do iT=1,nTempMagn
      write(u6,*)
      write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)
      write(u6,'(A,F9.4,A)') 'Magnetic field strength = ',H(iH),' tesla'
      write(u6,'(12x,A,F9.4,A)') 'Temperature = ',TempMagn(iT),' kelvin'
      write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

      write(u6,'(2x,A,3(10x,A))') 'Angle |','rotation in the YZ plane          |','rotation in the XZ plane          |', &
                                  'rotation in the XY plane          |'
      write(u6,'(10A)') '--------|',('Spin torque X |','Spin torque Y |','Spin torque Z |',i=1,3)
      do iM=1,AngPoints
        write(u6,'(F7.3,1x,A,3(ES13.6,1x,A,ES13.6,1x,A,ES13.6,1x,A))') Ang(iM),'|', &
                                                                       (sx(iPl,iM,iH,iT),' ',sy(iPl,iM,iH,iT),' ', &
                                                                        sz(iPl,iM,iH,iT),'|',iPl=1,3)
      end do
      write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)
    end do !iT
  end do !iH

end if !smagn

! 199 continue
!-----------------------------------------------------------------------
! Deallocate memory for this calculation:
call mma_deallocate(Wex)

call mma_deallocate(Zex)
call mma_deallocate(SEX)
call mma_deallocate(MEX)
call mma_deallocate(ZT)
call mma_deallocate(ST)
call mma_deallocate(MT)

call mma_deallocate(ZL)
call mma_deallocate(SL)
call mma_deallocate(ML)
call mma_deallocate(ZR)
call mma_deallocate(SR)
call mma_deallocate(MR)

call mma_deallocate(ZRT)
call mma_deallocate(ZLT)
call mma_deallocate(MRT)
call mma_deallocate(MLT)
call mma_deallocate(SRT)
call mma_deallocate(SLT)

call mma_deallocate(tx)
call mma_deallocate(ty)
call mma_deallocate(tz)
call mma_deallocate(sx)
call mma_deallocate(sy)
call mma_deallocate(sz)

call mma_deallocate(WL)
call mma_deallocate(WR)

call mma_deallocate(Ang)
call mma_deallocate(dX)
call mma_deallocate(dY)
call mma_deallocate(dZ)

call mma_deallocate(H)

#ifdef _DEBUGPRINT_
write(u6,*) 'TORQ: allocated memory was sucessfully deallocated'
#endif

return

end subroutine torque_pa
