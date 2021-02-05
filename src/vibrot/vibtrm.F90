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

! Calculation of transition matrix elements over vibrational
! wave functions of a given series of observables by a Simpson
! quadrature. The vibrational wave functions have been computed
! in vibrot and stored on units Vibwvs1 = 12 (potential 1) and
! Vibwvs2 = 13 (potential 2) for each rotational quantum number.
! The observables are read as input for a sequence of r-values
! and fitted to an analytical form in subroutine Pot.
! Integration is performed on a logarithmic scale (u=ln(r))
! between Umin and Umax with a grid size del corresponding to
! ndim integration steps. ndim has to be odd.
! Lifetimes are computed for the upper state (potential 2)
! no degeneracy factor is included in these values.
! Oscillator strengths are computed for all vibrational levels.
!
! ********** MOLCAS Release 91 05 01 **********
subroutine Vibtrm(ndim,Umin,Umax,Teas,R,PotR,Title)

use Vibrot_globals, only: iad12, iad13, iallrot, J1A, J1B, J2A, J2B, nvib1, nvib21, Vibwvs1, Vibwvs2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Four, c_in_au, auTofs, auTocm
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ndim
real(kind=wp), intent(in) :: Umin, Umax, Teas, R(ndim), PotR(ndim)
character(len=80), intent(in) :: Title
integer(kind=iwp) :: i, iadr1, iadr2, ipot, iPrint, ist1, ist2, iwr, J1, J2, Jad1, Jad2, JEndA, JEndB, ndim1, ne, ne1, ne2, nemax, &
                     ne2max, nv, nv1, nv11, nv12, nv12s, nv2, nv22
real(kind=wp) :: ChkSum, del, Obsr, Smax, Afact
integer(kind=iwp), allocatable :: nv1w(:), nv2w(:)
real(kind=wp), allocatable :: Vib(:,:), X(:), E(:,:), Tau(:), Taui(:,:), Osc(:,:), Obs(:), S(:), Sw(:)
character(len=80) :: TmpLine

ChkSum = Zero

!PAM97 write(u6,1001) Title,Teas
write(u6,*)
call CollapseOutput(1,'Matrix elements of observable: '//Title)
write(u6,'(a,es15.8)') ' Asymptotic energy difference (au):',Teas
iPrint = 2

! check dimensions
ne1 = nvib1+1
ne2 = nvib21+1
ndim1 = ndim+1
nemax = max(ne1,ne2)
ne2max = max(ne1*ne2,ne1*(ne1+1)/2,ne2*(ne2+1)/2)

! Allocate memory
call mma_allocate(Vib,ndim1,ne1+ne2,label='Vib')
call mma_allocate(X,ndim,label='X')
call mma_allocate(E,2,nemax,label='E')
call mma_allocate(Tau,ne2,label='Tau')
call mma_allocate(Taui,ne1,ne2,label='Taui')
call mma_allocate(Osc,ne1,ne2,label='Osc')
call mma_allocate(OBs,ne1*ne2,label='Obs')
call mma_allocate(S,ne1*(ne1+1)+ne2*(ne2+1),label='S')
call mma_allocate(Sw,ne2max,label='Sw')
call mma_allocate(nv1w,ne2max,label='nv1w')
call mma_allocate(nv2w,ne2max,label='nv2w')

! Set up grid size
del = (Umax-Umin)/(ndim-1)

! Loop over rotational quantum numbers
JEndA = J1A
JEndB = J1B
if (iallrot == 1) then
  JEndA = J2A
  JEndB = J2B
end if
do J1=J1A,JEndA
  Jad1 = J1-J1A+1
  do J2=J1B,JEndB
    Jad2 = J2-J1B+1
    write(u6,*)
    write(TmpLine,1002) J1,J2
    call CollapseOutput(1,TmpLine)
    write(u6,1003)

    ! Read vibrational functions for this pair of J-values.
    do ipot=1,2
      if (ipot == 1) then
        ne = ne1
        iadr1 = iad12(Jad1)
        do nv=1,ne1
          call DDafile(Vibwvs1,2,Vib(:,nv),ndim1,iadr1)
          E(1,nv) = Vib(ndim,nv)
        end do
      else
        ne = ne2
        iadr2 = iad13(Jad2)
        do nv=1,ne2
          call DDafile(Vibwvs2,2,Vib(:,ne1+nv),ndim1,iadr2)
          E(2,nv) = Vib(ndim,ne1+nv)+Teas
        end do
      end if

      ! compute overlap matrix S
      nv12 = (ipot-1)*(ne1**2+ne1)/2
      do nv1=1,ne
        ist1 = nv1+(ipot-1)*ne1
        do nv2=1,nv1
          ist2 = nv2+(ipot-1)*ne1
          nv12 = nv12+1
          ! Set up scalar product and integrate
          do i=1,ndim
            X(i) = Vib(i,ist1)*Vib(i,ist2)*R(i)**2
          end do
          call Simpsn(X,del,ndim,S(nv12))
        end do
      end do

      ! check overlap matrix for non-orthogonality
      nv12 = (ipot-1)*(ne1**2+ne1)/2
      Smax = Zero
      do nv1=1,ne
        do nv2=1,nv1
          nv12 = nv12+1
          if (nv1 == nv2) cycle
          if (abs(S(nv12)) > abs(Smax)) Smax = S(nv12)
        end do
      end do
      if (abs(Smax) > 1.0e-4_wp) write(u6,1200) ipot,Smax

      ! Print overlap matrix
      if (iPrint >= 1) then
        write(u6,*)
        write(u6,'(1x,A,I3)') 'Overlap matrix for vibrational wave functions for potential number',ipot
        iwr = 0
        nv12 = (ipot-1)*(ne1**2+ne1)/2
        do nv1=1,ne
          do nv2=1,nv1
            iwr = iwr+1
            nv12 = nv12+1
            nv1w(iwr) = nv1-1
            nv2w(iwr) = nv2-1
            Sw(iwr) = S(nv12)
          end do
        end do
        write(u6,'(6(3X,2I3,F12.6))') (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
      end if
    end do

    ! Compute overlap matrix between pot1 and pot2 functions
    nv11 = 0
    nv12s = (ne1**2+ne1+ne2**2+ne2)/2
    nv12 = nv12s
    do nv1=1,ne1
      nv11 = nv11+nv1
      nv22 = (ne1**2+ne1)/2
      do nv2=1,ne2
        nv12 = nv12+1
        nv22 = nv22+nv2
        ! Set up scalar product and integrate
        do i=1,ndim
          X(i) = Vib(i,nv1)*Vib(i,nv2+ne1)*R(i)**2
        end do
        call Simpsn(X,del,ndim,S(nv12))
        S(nv12) = S(nv12)/sqrt(S(nv11)*S(nv22))
        ChkSum = ChkSum + S(nv12)
      end do
    end do

    ! print overlap matrix
    write(u6,1310)
    nv12 = nv12s-ne2
    do nv1=1,ne1
      nv12 = nv12+ne2
      write(u6,'(6(3X,2I3,F12.6))') (nv1-1,i-1,S(i+nv12),i=1,ne2)
    end do
    ! compute transition moment
    ! Afact = 4*pi*e^2*E_h^2 / 3*eps_0*m_e*c^3*h^2
    ! numerically: 4/(3*c^3) (in a.u. of time ^ -1) (converted to 1/ns)
    Afact = Four/(Three*c_in_au**3)/auToFs*1.0e6_wp
    nv11 = 0
    nv12 = 0
    do nv1=1,ne1
      nv11 = nv11+nv1
      nv22 = (ne1**2+ne1)/2
      do nv2=1,ne2
        nv12 = nv12+1
        nv22 = nv22+nv2
        ! Set up scalar product and integrate
        do i=1,ndim
          X(i) = Vib(i,nv1)*Vib(i,nv2+ne1)*PotR(i)*R(i)**2
        end do
        call Simpsn(X,del,ndim,Obsr)
        Obs(nv12) = Obsr/sqrt(S(nv11)*S(nv22))

        ! Taui is the contribution to the inverse lifetime from
        ! vibrational states nv1 (lower) and nv2 (upper)
        Taui(nv1,nv2) = Afact*(E(2,nv2)-E(1,nv1))**3*Obs(nv12)**2
        ! computed oscillator strengths
        Osc(nv1,nv2) = Two*(E(2,nv2)-E(1,nv1))*Obs(nv12)**2/Three

        !ChkSum = ChkSum + Taui(nv1,nv2)
      end do
    end do

    ! write matrix elements
    write(u6,*)
    write(u6,'(1x,A)') 'Transition moments over vibrational wave functions (atomic units)'
    iwr = 0
    do nv1=1,ne1
      do nv2=1,ne2
        iwr = iwr+1
        nv1w(iwr) = nv1-1
        nv2w(iwr) = nv2-1
        Sw(iwr) = Obs(iwr)
      end do
    end do
    write(u6,'(6(2X,I3,1X,I3,1X,F11.6))') (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
    write(u6,*)
    write(u6,'(1x,A)') 'Energy differences for vibrational wave functions (atomic units)'
    iwr = 0
    do nv1=1,ne1
      do nv2=1,ne2
        iwr = iwr+1
        nv1w(iwr) = nv1-1
        nv2w(iwr) = nv2-1
        Sw(iwr) = E(2,nv2)-E(1,nv1)
      end do
    end do
    write(u6,'(6(2X,I3,1X,I3,1X,F11.6))') (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
    write(u6,'(1x,A)') 'Energy differences for vibrational wave functions (cm-1)'
    iwr = 0
    do nv1=1,ne1
      do nv2=1,ne2
        iwr = iwr+1
        nv1w(iwr) = nv1-1
        nv2w(iwr) = nv2-1
        Sw(iwr) = (E(2,nv2)-E(1,nv1))*auTocm
      end do
    end do
    write(u6,'(6(2X,I3,1X,I3,1X,F11.1))') (nv1w(i),nv2w(i),Sw(i),i=1,iwr)

    ! Print oscillator strengths
    write(u6,'(1x,A)') 'Oscillator strengths for vibrational wave functions'
    iwr = 0
    do nv1=1,ne1
      do nv2=1,ne2
        iwr = iwr+1
        nv1w(iwr) = nv1-1
        nv2w(iwr) = nv2-1
        Sw(iwr) = Osc(nv1,nv2)
      end do
    end do
    write(u6,'(6(2X,I3,1X,I3,1X,E12.5))') (nv1w(i),nv2w(i),Sw(i),i=1,iwr)

    ! Print lifetime
    write(u6,*)
    write(u6,'(1x,A,/,A)') 'Contributions to inverse lifetimes (ns-1)','No degeneracy factor is included in these values.'
    iwr = 0
    do nv1=1,ne1
      do nv2=1,ne2
        iwr = iwr+1
        nv1w(iwr) = nv1-1
        nv2w(iwr) = nv2-1
        Sw(iwr) = Taui(nv1,nv2)
      end do
    end do
    write(u6,'(6(2X,I3,1X,I3,1X,E12.5))') (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
    ! Compute lifetimes for upper state
    write(u6,1600)
    do nv2=1,ne2
      Tau(nv2) = Zero
      do nv1=1,ne1
        Tau(nv2) = Tau(nv2)+Taui(nv1,nv2)
      end do
      Tau(nv2) = One/Tau(nv2)
      !ChkSum = ChkSum + Tau(nv2)
      write(u6,1700) nv2-1,Tau(nv2)
    end do
    call CollapseOutput(0,TmpLine)

  ! End of loop over rotational quantum number
  end do
end do
call CollapseOutput(0,'Matrix elements of observable: '//Title)

! Deallocate memory
call mma_deallocate(Vib)
call mma_deallocate(X)
call mma_deallocate(E)
call mma_deallocate(Tau)
call mma_deallocate(Taui)
call mma_deallocate(Osc)
call mma_deallocate(Obs)
call mma_deallocate(S)
call mma_deallocate(Sw)
call mma_deallocate(nv1w)
call mma_deallocate(nv2w)

!PAM97 New code begins:
!do J1=J1A,J2A
!  do J2=J1B,J2B
!    ! Read vibrational function for pot-1, with rot q.n.=J1:
!    iadr1 = iad12(J1-J1A+1)
!    do nv=0,nvib1
!      call Dafile(Vibwvs1,2,Vibbuf,RtoI*ndim,iadr1)
!      call dcopy_(ndim,Vibbuf,1,Vib1(1,nv),1)
!      ERoVib1(nv,J1) = Vibbuf(ndim+1)
!    end do
!    ! Normalize vibrational function:
!    do nv=0,nvib1
!      do i=1,ndim
!        Vibbuf(i) = Vib1(i,nv)**2*R(i)**2
!      end do
!      call Simpsn(Vibbuf,del,ndim,ovlp)
!      xnrm = One/sqrt(ovlp)
!      do i=1,ndim
!        Vib1(i,nv) = xnrm*Vib1(i,nv)
!      end do
!    end do
!    ! Check orthonormality:
!    ovlmax = Zero
!    do nv1=1,nvib1
!      do nv2=0,nv1-1
!        do i=1,ndim
!          Vibbuf(i) = Vib1(i,nv1)*Vib1(i,nv2)*R(i)**2
!        end do
!        call Simpsn(Vibbuf,del,ndim,ovlp)
!        ovlmax = max(ovlp,ovlmax)
!      end do
!    end do
!    if (ovlmax >= 1.0e-4_wp) then
!      write(u6,*) ' Warning from VIBTRM: Vibrational wave functions for potential 1'
!      write(u6,*) ' are not orthonormal. Max overlap:',ovlmax
!    end if
!    ! Read vibrational function for pot-2, with rot q.n.=J2:
!    iadr2 = iad13(J2-J1B+1)
!    do nv=0,nvib21
!      call Dafile(Vibwvs2,2,Vibbuf,RtoI*ndim,iadr2)
!      call dcopy_(ndim,Vibbuf,1,Vib2(1,nv),1)
!      ERoVib2(nv,J2) = Vibbuf(ndim+1)
!    end do
!    ! Normalize vibrational function:
!    do nv=0,nvib21
!      do i=1,ndim
!        Vibbuf(i) = Vib2(i,nv)**2*R(i)**2
!      end do
!      call Simpsn(Vibbuf,del,ndim,ovlp)
!      xnrm = One/sqrt(ovlp)
!      do i=1,ndim
!        Vib2(i,nv) = xnrm*Vib2(i,nv)
!      end do
!    end do
!    ! Check orthonormality:
!    ovlmax = Zero
!    do nv1=1,nvib21
!      do nv2=0,nv1-1
!        do i=1,ndim
!          Vibbuf(i) = Vib2(i,nv1)*Vib2(i,nv2)*R(i)**2
!        end do
!        call Simpsn(Vibbuf,del,ndim,ovlp)
!        ovlmax = max(ovlp,ovlmax)
!      end do
!    end do
!    if (ovlmax >= 1.0e-4_wp) then
!      write(u6,*) ' Warning from VIBTRM: Vibrational wave functions for potential 2'
!      write(u6,*) ' are not orthonormal. Max overlap:',ovlmax
!    end if
!    ! Compute Franck-Condon overlaps and transition moment integrals:
!    do nv1=0,nvib1
!      do nv2=0,nvib21
!        if (J1 == J2) then
!          do i=1,ndim
!            Vibbuf(i) = Vib1(i,nv1)*Vib2(i,nv2)*R(i)**2
!          end do
!          call Simpsn(Vibbuf,del,ndim,FC(nv1,nv2,J1))
!        end if
!        if (abs(nv1-nv2) <= NSel) then
!          do i=1,ndim
!            X(i) = PotR(i)*Vibbuf(i)
!          end do
!          call Simpsn(X,del,ndim,Obs(nv1,nv2,J1,J2))
!        else
!          TD(nv1,nv2,J1,J2) = Zero
!        end if
!        EDiff = ERoVib2(nv2,J2)-ERoVib1(nv1,J1)+TeDiff
!        TauPrt(nv1,nv2,J1,J2) = Afact*EDiff**3*TD(nv1,nv2,J1,J2)**2
!      end do
!    end do
!  end do
!end do
!PAM97 New code ends

call Add_Info('VIBROT_VIBTRM',[ChkSum],1,6)

return

!PAM97
!1001 format(1h1,1x,'Matrix elements of observables'/1x,e80//1x,'Asymtotic energy difference (au):',e14.6)
1002 format('Rotational quantum number for potential 1: ',i3,', for potential 2: ',i3)
1003 format(1x,80('-'))
1200 format(/1x,'*****Warning: non-orthogonality between vibrational',1x,'wave functions for potential',i2 &
            /13x,'largest overlap matrix element is',f14.6)
1310 format(/1x,'Overlap matrix for pot-1 and pot-2 functions')
1600 format(/1x,'Lifetimes (in nano seconds)'/3x,'v',7x,'tau')
1700 format(1x,i3,f10.2)

end subroutine Vibtrm
