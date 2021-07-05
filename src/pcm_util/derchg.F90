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

subroutine DerChg(nAt,nAt3,nTs,nS,Eps,IAtm,AtmC,AtmChg,DM,DerMat,Tessera,Q,Qtot,QDer,DerTes,DerPunt,DerCentr,DerRad,Der1,Der2, &
                  VDer,Sphere,ISphe)

implicit real*8(A-H,O-Z)
integer IAtm(nAt)
dimension AtmC(3,nAt), AtmChg(nAt)
dimension Q(2,*), ISphe(*), QTot(*), QDer(3,nAt,*), Der1(*), Der2(*)
dimension Tessera(4,*), Sphere(4,*), VDer(nTs,*)
dimension DM(nTs,*), DerMat(nTs,*)
dimension DerTes(nTs,NAt,3), DerPunt(nTs,NAt,3,3)
dimension DerRad(nS,NAt,3), DerCentr(nS,NAt,3,3)
data One,Two,Four/1.0d0,2.0d0,4.0d0/

PI = Four*atan(One)
FPI = Four*PI
Diag = -1.0694d0*sqrt(FPI)/Two
Sc_Cond = (Eps-One)/Eps

! Total charges

do iTs=1,nTs
  Qtot(iTs) = Q(1,iTs)+Q(2,iTs)
end do

! Compute the derivative of PCM matrix

! Loop on atoms and coordinates
do iAt=1,nAt
  do iCoord=1,3
    Index = 3*(iAt-1)+iCoord
    ! Conductor-like case
    call DMat_CPCM(iAt,iCoord,Eps,nTs,nS,nAt,Diag,Tessera,DerMat,DerTes,DerPunt,DerCentr,iSphe)

    ! First product: DerDM * Qtot

    call PrMatVec(.false.,.false.,DerMat,-One,nTs,nTs,Qtot,Der1)
    ! pcm_solvent
    !if ((iat == 5) .and.( icoord == 1)) then
    !  write(6,'(a)') 'In DerChg first contribution for 5, 1'
    !  do its=1,nts
    !    write(6,'(i4,f20.12)') its,der1(its)
    !  end do
    !end if
    ! pcm_solvent end

    ! Total deriv. of the potential summed up the quantity alread
    ! computed (DerDM*Qtot)

    do iTs=1,nTs
      Der1(iTs) = Der1(iTs)+Sc_Cond*VDer(iTs,Index)
    end do

    ! Last product: - DM^-1 * (V^x + DM^x*q)

    ! pcm_solvent
    !if ((iat == 5) .and. (icoord == 1)) then
    !  write(6,'(a)') 'In DerChg second contribution for 5, 1'
    !  do its=1,nts
    !    write(6,'(i4,f20.12)') its,der1(its)
    !  end do
    !end if
    ! pcm_solvent end
    call PrMatVec(.false.,.false.,DM,-1.d0,nTs,nTs,Der1,Der2)
    call FillQDer(nAt,nTs,iAt,iCoord,Der2,QDer)
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(nAt3)
  call Unused_integer_array(IAtm)
  call Unused_real_array(AtmC)
  call Unused_real_array(AtmChg)
  call Unused_real_array(DerRad)
  call Unused_real_array(Sphere)
end if

end subroutine DerChg
!====
subroutine DMat_CPCM(iAt,iC,Eps,nTs,nS,nAt,fact,Tessera,DerMat,DerTes,DerPunt,DerCentr,iSphe)

implicit real*8(A-H,O-Z)
dimension ISphe(*), Tessera(4,*), DerMat(nTs,*)
dimension DerTes(nTs,NAt,3), DerPunt(nTs,NAt,3,3)
dimension DerCentr(nS,NAt,3,3)

! Compute the derivative of the CPCM matrix wrt atom iat, coord. ic

! Loop on tesserae
do ITs=1,NTs
  L = ISPHE(ITs)
  do JTS=1,NTs
    LJ = ISPHE(JTS)
    ! Diagonal elements
    if (ITs == JTS) then
      DerMat(ITs,ITs) = fact*DERTES(ITs,iAt,IC)/(Tessera(4,ITs)*sqrt(Tessera(4,ITs)))
    else
      ! Off diagonal elements
      XIJ = Tessera(1,ITs)-Tessera(1,JTS)
      YIJ = Tessera(2,ITs)-Tessera(2,JTS)
      ZIJ = Tessera(3,ITs)-Tessera(3,JTS)
      DIJ = sqrt(XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ)
      DXIJ = DERPUNT(ITs,iAt,IC,1)+DERCENTR(L,iAt,IC,1)-DERPUNT(JTS,iAt,IC,1)-DERCENTR(LJ,iAt,IC,1)
      DYIJ = DERPUNT(ITs,iAt,IC,2)+DERCENTR(L,iAt,IC,2)-DERPUNT(JTS,iAt,IC,2)-DERCENTR(LJ,iAt,IC,2)
      DZIJ = DERPUNT(ITs,iAt,IC,3)+DERCENTR(L,iAt,IC,3)-DERPUNT(JTS,iAt,IC,3)-DERCENTR(LJ,iAt,IC,3)
      PROD = (XIJ*DXIJ+YIJ*DYIJ+ZIJ*DZIJ)/DIJ**3
      DerMat(ITs,JTS) = -PROD
    end if
  end do
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_real(Eps)

end subroutine DMat_CPCM
!====
subroutine PrMatVec(Dag,DoSym,Mat,f,n,m,Vec,Res)

implicit real*8(A-H,O-Z)
real*8 Mat(n,*), Vec(*), Res(*)
logical Dag, DoSym
data Zero,Two/0.0d0,2.0d0/

! Do the matrix vector product: f*Mat(n,m)*Vec(m,1)=Res(n,1)
! possibly transposed (if Dag): f*Vec(1,m)*Mat(m,n)=Res(1,n)
! If DoSym symmetrize the matrix elements

do i=1,n
  Res(i) = Zero
  do j=1,m
    if (DoSym) then
      ElM = (Mat(i,j)+Mat(j,i))/Two
    else
      if (Dag) ElM = Mat(j,i)
      if (.not. Dag) Elm = Mat(i,j)
    end if
    Res(i) = Res(i)+f*ElM*Vec(j)
  end do
end do

return

end subroutine PrMatVec
!====
subroutine FillQDer(nAt,nTs,iAt,iC,Der,QDer)

implicit real*8(A-H,O-Z)
dimension Der(*), QDer(3,nAt,*)

do iTs=1,nTs
  QDer(iC,iAt,iTs) = Der(iTs)
end do

return

end subroutine FillQDer
!====
subroutine testq(nAt,nTs,VDer,q,qtot)

implicit real*8(A-H,O-Z)
dimension VDer(nTs,*), Q(2,*), QTot(*)
integer Lu

Lu = 1
call Molcas_open(Lu,'DerPt.dat')
!open(1,file='DerPot.dat',status='old',form='formatted')
do iAt=1,nAt
  do iCoord=1,3
    Index = 3*(iAt-1)+iCoord
    do iTs=1,nTs
      read(1,*) VDer(iTs,Index)
    end do
  end do
end do
close(1)
do iAt=1,nAt
  do iCoord=1,3
    Index = 3*(iAt-1)+iCoord
    sum = 0.d0
    do its=1,nts
      qtot(its) = q(1,its)+q(2,its)
      sum = sum+qtot(its)*VDer(its,index)
    end do
    write(6,'("Charges times VDer",i4,f20.12)') index,sum
  end do
end do

return

end subroutine testq
