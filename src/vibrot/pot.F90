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

Subroutine POT(Rin,Ein,Rout,Eout,nout,ifit,Emin,Req,R0,R1,dR,npin,Title,iplot,Redm,sc,Nr)

use Vibrot_globals, only: iscale, nop
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: Rin(*), Ein(*), Rout(*), Eout(*), Redm
real(kind=wp), intent(inout) :: R0, R1, dR
integer(kind=iwp), intent(in) :: nout, ifit, npin, iplot, Nr
real(kind=wp), intent(out) :: Emin, Req, sc
character(len=8) :: FILENAME
real(kind=wp) :: Eeq, Ree, Ue, Re, Eminx, Escale, alpha, Einf, Einfp
integer(kind=iwp) :: Vibplt, iPrint, next, nplot, i, imin
character(len=80) :: Title
character(len=12), parameter :: tp(3)=[character(len=12) :: 'Max point','Saddle point','Min point']
integer(kind=iwp), parameter :: lext=100, ipldim=1000
real(kind=wp) :: Rext(lext), Eext(lext), Eextp(lext), Rplot(ipldim), Eplot(ipldim), Eplotp(ipldim)
integer(kind=iwp) :: iext(lext)
integer(kind=iwp), external :: IsFreeUnit

FILENAME = 'VIBPLT0 '
if ((Nr >= 1).and.(Nr <= 9)) then
  write(FILENAME(7:7),'(I1)') Nr
else if ((Nr >= 10).and.(Nr <= 99)) then
  write(FILENAME(7:8),'(I2)') Nr
end if
write (u6,*) 'Generating plot file:',FILENAME
Vibplt = IsFreeUnit(10)
call Molcas_Open(Vibplt,FILENAME)

iPrint = 1

select case (ifit)
  case (1)
    Ue = 0.4_wp
    if ((Ue < 0.1_wp).or.(Ue > 0.9_wp)) then
      write(u6,*) 'POT Error: Ue should be in 0.1..0.9'
      write(u6,*) '           Ue  =',Ue
      call Quit_OnUserError()
    end if

    ! Find Re value

    Re = Rin(1)
    Eminx = Ein(1)
    do i=2,npin
      if (Ein(i) <= Eminx) then
        Re = Rin(i)
        Eminx = Ein(i)
      end if
    end do

    ! Scale input potential if requested such that
    ! the binding energy is 0.1 au (BOR 0601).

    sc = One
    Escale = Ein(npin)
    if (iscale /= 0) sc = abs(0.1_wp/(Eminx-Escale))
    do i=1,npin
      Ein(i) = Escale+sc*(Ein(i)-Escale)
    end do
    Redm = Redm/sc
    write(u6,1001) sc

    if ((Re < 1.0_wp).or.(Re > 20.0_wp)) then
      write(u6,*) 'POT Error: Re should be in 1.0..20.0'
      write(u6,*) '           Re  =',Re
      call Quit_OnUserError()
    end if
    alpha = log(Ue)/Re
    do i=1,npin
      Rin(i) = exp(alpha*Rin(i))
    end do

    do i=1,nout
      Rout(i) = exp(alpha*Rout(i))
    end do
    Rout(nout+1) = Zero
    call Sort_Pot(Rin,Ein,npin)
    next = lext
    call Spline(Rin,Ein,npin,Rout,Eout,nout+1,Rext,Eext,iext,next,1)
    do i=1,npin
      Rin(i) = log(Rin(i))/alpha
    end do
    call Sort_Pot(Rin,Ein,npin)
    do i=1,nout
      Rout(i) = log(Rout(i))/alpha
    end do
    do i=1,next
      Rext(i) = log(Rext(i))/alpha
    end do
    if (iPrint >= 1) then
      do i=1,next
        Eextp(i) = Escale+(Eext(i)-Escale)/sc
      end do
      if (next >= 1) then
        write(u6,2002)
        write(u6,2003) (tp(iext(i)),Rext(i),Eextp(i),i=1,next)
      end if
    end if
    if (iplot >= 1) then
      nplot = 1+int((R1-R0)/dR)
      if ((nplot > ipldim).or.(nplot <= 0)) then
        write(u6,*) 'POT Error: Variable NPLOT should be in 1..IPLDIM'
        write(u6,*) '          IPLDIM=',IPLDIM
        write(u6,*) '          NPLOT =',NPLOT
        call Quit_OnUserError()
      end if
      do i=1,nop
        Rin(i) = exp(alpha*Rin(i))
      end do
      call Sort_Pot(Rin,Ein,nop)
      if (iplot == 1) then
        do i=1,nplot
          Rplot(i) = (i-1)*dR+R0
        end do
      else
        do i=1,nplot
          Rplot(i) = 10.0_wp**((i-1)*dR+R0)
        end do
      end if
      do i=1,nplot
        Rplot(i) = exp(alpha*Rplot(i))
      end do
      next = lext
      call Spline(Rin,Ein,npin,Rplot,Eplot,nplot,Rext,Eext,iext,next,1)
      do i=1,nplot
        Rplot(i) = log(Rplot(i))/alpha
      end do
      if (iplot == 3) then
        do i=1,nplot
          Rplot(i) = log10(Rplot(i))
        end do
      end if
      do i=1,nplot
        Eplotp(i) = Escale+(Eplot(i)-Escale)/sc
      end do
      write(Vibplt,3000) Title
      write(Vibplt,3002) nplot
      write(Vibplt,3001) (Rplot(i),Eplotp(i),i=1,nplot)
    end if
    if (next <= 0) then
      write(u6,*) 'POT Error: Variable NEXT should be larger than 0'
      write(u6,*) '          NEXT =',NEXT
      call Quit_OnUserError()
    end if
    Eeq = huge(Eeq)
    imin = 0
    do i=1,next
      if (Eeq > Eext(i)) then
        imin = i
        Eeq = Eext(i)
      end if
    end do
    Emin = Eeq
    Ree = Rext(imin)
    if (iplot > 0) Ree = log(Ree)/alpha
    Req = Ree
    if (iext(imin) /= 3) then
      write(u6,*) 'POT Error: IEXT(IMIN) should be = 3'
      write(u6,*) '    IEXT(IMIN) =',IEXT(IMIN)
      call Quit_OnUserError()
    end if
    Einf = Eout(nout+1)
    do i=1,nout
      Eout(i) = Eout(i)-Einf
    end do
    Emin = Emin-Einf
    Einfp = Escale+(Einf-Escale)/sc
    write(u6,1002) Einfp

    close(Vibplt)

  ! Here for fitting of observable input
  case (2)
    call Sort_Pot(Rin,Ein,npin)
    next = lext
    call Spline(Rin,Ein,npin,Rout,Eout,nout,Rext,Eext,iext,next,1)
    if (iplot >= 1) then
      nplot = 1+int((R1-R0)/dR)
      if ((nplot > ipldim).or.(nplot <= 0)) then
        write(u6,*) 'POT Error: Variable NPLOT should be in 1..IPLDIM'
        write(u6,*) '          IPLDIM=',IPLDIM
        write(u6,*) '          NPLOT =',NPLOT
        call Quit_OnUserError()
      end if
      if (iplot == 1) then
        do i=1,nplot
          Rplot(i) = (i-1)*dR+R0
        end do
      else
        do i=1,nplot
          Rplot(i) = 10.0_wp**((i-1)*dR+R0)
        end do
      end if
      call Spline(Rin,Ein,npin,Rplot,Eplot,nplot,Rext,Eext,iext,next,1)
      if (iplot == 3) then
        do i=1,nplot
          Rplot(i) = log10(Rplot(i))
        end do
      end if
      write(Vibplt,3000) Title
      write(Vibplt,3002) nplot
      write(Vibplt,3001) (Rplot(i),Eplot(i),i=1,nplot)
    end if
    sc = One
    Req = Zero
    Emin = Zero

    close(Vibplt)

  case default
    write(u6,*) 'POT Error: IFIT variable must be 1 or 2.'
    write(u6,*) '           IFIT=',IFIT
    call Quit_OnUserError()
end select

return

1001 format(/1x,'Scaling parameter for potential:',f12.6)
1002 format(1x,'Extrapolated value at infinity',f13.6)
2002 format(/1x,'extremum points'/24x,'R(au)',9x,'Value')
2003 format(1x,a,2f14.6)
3000 format(a80)
3001 format(1x,f15.8,f20.8)
3002 format(i4)

end subroutine POT
