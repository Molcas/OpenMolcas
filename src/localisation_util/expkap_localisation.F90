!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2026, Lila Zapp                                        *
!***********************************************************************

subroutine expkap_localisation(kappa,nOrb2Loc,Umat)
! analogous to exp_series in scf, only for the specified orbital subspace
!
! purpose: returns exp(kappa) as Umat
!
! kappa_cnt and xkappa_cnt are auxiliary matrices of dim (norb2loc,norb2loc) required for the Taylor expansion.
! they need to be allocated outside of this subroutine

use constants, only: Zero
use definitions, only: wp,iwp,u6

implicit none

integer(kind=iwp), intent(in) :: nOrb2Loc
real(kind=wp), intent(inout) :: kappa(nOrb2Loc,nOrb2Loc),Umat(nOrb2Loc,nOrb2Loc)
integer(kind=iwp) i, j
logical(kind=iwp), parameter:: OLD=.True.

!Call RecPrt('Kappa',' ',Kappa,nOrb2Loc,nOrb2Loc)
Do i=1, nOrb2Loc
   Do j=1, i-1
      If (Abs(Kappa(i,j)+Kappa(j,i))>1.0e-16_wp) Then
         Write (u6,*) ' Mismatch of off-diagonal value in Kappa matrix'
         Call Abend()
      End If
   End Do
   If (Kappa(i,i)/=Zero) Then
      Write (u6,*) ' Error in diagonal value in Kappa matrix'
      Call Abend()
   End If
End Do


If (Old) Then
Call expkap_localisation_Local(kappa,nOrb2Loc,Umat)
Else
Block
   use constants, only: One
   use stdalloc, only: mma_allocate, mma_deallocate
   integer(kind=iwp) ideg, m, ldh, lwsp, iexph, ns, iflag, i, j
   real(kind=wp) t
   real(kind=wp), allocatable:: wsp(:)
   integer(kind=iwp), allocatable:: ipiv(:)
   ideg=13
   m=nOrb2Loc
   t = One
   ldh=m
   lwsp=4*m**2+ideg+1
   Call mma_allocate(wsp,lwsp,Label='wsp')
   Call mma_allocate(ipiv,m,Label='ipiv')
   call dgpadm(ideg,m,t,kappa,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)

   If (iflag==0) Then
      Do j = 1, m
         Do i = 1, m
            Umat(i,j) = wsp(iexph)
            iexph=iexph+1
         End Do
      End Do
   Else
      Write(u6,*) 'DGPADM: iflag=',iflag
      Call Abend()
   End If
   Call mma_deallocate(ipiv)
   Call mma_deallocate(wsp)
End Block
End If

Block
use constants, only: One
use stdalloc, only: mma_allocate, mma_deallocate
real(kind=wp), allocatable:: UnitMatrix(:,:)
Call mma_allocate(UnitMatrix,nOrb2Loc,nOrb2Loc,Label='UM')
UnitMatrix(:,:)=Zero

Call dgemm_('N','T',nOrb2Loc,nOrb2Loc,nOrb2Loc, &
             One, UMat, nOrb2Loc,               &
                  UMat, nOrb2Loc,           &
             Zero, UnitMatrix, nOrb2Loc)

Do i = 1, nOrb2Loc
   Do j = 1, i-1
      If (Abs(UnitMatrix(i,j))>1.0e-12_wp) Then
         Write (u6,*) 'Too large off diagonal value'
         Call Abend()
      End If
      If (Abs(UnitMatrix(j,i))>1.0e-12_wp) Then
         Write (u6,*) 'Too large off diagonal value'
         Call Abend()
      End If
   End Do
   If (Abs(One-UnitMatrix(i,i))>1.0e-12_wp) Then
      Write (u6,*) 'Diagonal not accurate!'
      Write (u6,*) 'UM(i,i)=',UnitMatrix(i,i)
      Call Abend()
   End If
End Do
Call mma_deallocate(UnitMatrix)
End Block

Contains

subroutine expkap_localisation_Local(kappa,nOrb2Loc,Umat)
! analogous to exp_series in scf, only for the specified orbital subspace
!
! purpose: returns exp(kappa) as Umat
!
! kappa_cnt and xkappa_cnt are auxiliary matrices of dim (norb2loc,norb2loc) required for the Taylor expansion.
! they need to be allocated outside of this subroutine

use definitions, only: wp,iwp,u6
use constants, only: Zero,One
use Localisation_globals, only: kappa_cnt, xkappa_cnt

implicit none

integer(kind=iwp), intent(in) :: nOrb2Loc
real(kind=wp), intent(inout) :: kappa(nOrb2Loc,nOrb2Loc),Umat(nOrb2Loc,nOrb2Loc)
real(kind=wp), parameter :: thrsh_taylor = 1.0e-16_wp
real(kind=wp) :: factor, ithrsh
integer(kind=iwp) :: cnt,maxel(2)
logical(kind=iwp), parameter :: debug_exp = .false.
real(kind=wp),External :: DDot_

kappa_cnt(:,:) = kappa !kappa^cnt = kappa since cnt=1
xkappa_cnt(:,:) = kappa_cnt

Umat(:,:) = Zero
call unitmat(Umat,nOrb2Loc)

cnt = 1
factor = One
ithrsh = One
maxel(:) = 0

Umat(:,:) =  Umat(:,:) - kappa(:,:)

if (debug_exp) then
    write(u6,*) 'Taylor expansion: n=1'
    call RecPrt('Umat = I + kappa^1',' ',Umat(:,:), nOrb2Loc, nOrb2Loc)
    call RecPrt('kappa',' ',kappa(:,:), nOrb2Loc, nOrb2Loc)
    write(u6,*) 'Taylor expansion: more terms'
end if


do while (ithrsh > thrsh_taylor)

    !the number of the term = the exponent for kappa in that term
    cnt = cnt+1

    !the faculty value that the matrix will be divided by
    factor = factor*DBLE(cnt)

    !calculate the cnt'th exponent of the kappa matrix
    ! initial kappa matrix (just kappa^1)
    ! C <= alpha*A*B + beta*C
    ! kappa_cnt <= 1*kappa_cnt*kappa + 0*kappa_cnt
    call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nOrb2Loc,&
                        (One/DBLE(cnt)),xkappa_cnt,nOrb2Loc,&
                                        kappa,nOrb2Loc,&
                        Zero,kappa_cnt,norb2Loc)

    xkappa_cnt(:,:) = kappa_cnt(:,:)


    ! differentiation of odd and even cases, because this expands exp(-kappa)
    ! all terms starting at n=2
    if (mod(cnt,2) == 0) then
        Umat(:,:) =  Umat + kappa_cnt(:,:)
        if (debug_exp) then
            write(u6,'(A,ES10.1,A,I2,A,ES12.4)') 'term: +',1/factor,' * kappa^',cnt, &
            ', current ithrsh = ', ithrsh
        end if
    else
        Umat(:,:) =  Umat - kappa_cnt(:,:)
        if (debug_exp) then
            write(u6,'(A,ES10.1,A,I2,A,ES12.4)') 'term: -',1/factor,' * kappa^',cnt, &
            ', current ithrsh = ', ithrsh
        end if
    end if

    ithrsh = sqrt(DDot_(nOrb2Loc**2,Kappa_Cnt(:,:),1,Kappa_Cnt(:,:),1))/(Norb2loc**2)
    !ithrsh = maxval(abs(Kappa_Cnt(:,:))/(abs(Umat)+thrsh_taylor))

    ! sanity check for divergence
    if (ithrsh/720 > One) then
        write(u6,*) "Bug: the Taylor expansion of exp(kappa) diverges - numerical error"
        write(u6,*) "Stopping Taylor expansion at ",cnt,"-th term. Rescale kappa before feeding it this subroutine. ithrsh",&
                     ithrsh
        call Abend()
    end if

    if (debug_exp) then
        write(u6,'(A,ES10.1,A,I2,A,ES12.4)') 'term:  ',1/factor,' * kappa^',cnt, &
            ', new ithrsh     = ', ithrsh
        call RecPrt('kappa^cnt',' ',kappa_cnt(:,:), nOrb2Loc, nOrb2Loc)
        call RecPrt('Umat',' ',Umat(:,:), nOrb2Loc, nOrb2Loc)
    end if

end do


maxel(:) = maxloc(abs(Umat(:,:)))
if (Umat(maxel(1),maxel(2)) > One+1.0e-8_wp) then
    write(u6,*) "element of U bigger cannot be larger than 1, as U is supposed to be orthogonal/unitary:", &
                 Umat(maxel(1),maxel(2))
    call Abend()
end if

#ifdef _DEBUGPRINT_
kappa_cnt(:,:)= Zero

call dgemm_('N','T',nOrb2Loc,nOrb2Loc,nOrb2Loc,&
                    One,Umat,nOrb2Loc,&
                        Umat,nOrb2Loc,&
                        Zero,kappa_cnt,norb2Loc)
call RecPrt('kappa',' ',kappa(:,:), nOrb2Loc, nOrb2Loc)
call RecPrt('Umat',' ',Umat(:,:), nOrb2Loc, nOrb2Loc)
call RecPrt('U^(T)U = I',' ',kappa_cnt(:,:), nOrb2Loc, nOrb2Loc)

#endif

end subroutine expkap_localisation_local
end subroutine expkap_localisation
