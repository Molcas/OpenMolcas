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
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

module ROMP2_Reference
use Definitions, only: wp, iwp, u6
use MBPT2_Global, only: MBPT2_All_Finite
use Constants, only: Zero, One, Half
use UMP2_Global
use Para_Info, only: nProcs,MyRank,Is_Real_Par
use ROMP2_Cholesky, only: ROMP2_Cho_Fock
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
implicit none
private
public :: ROMP2_Read, ROMP2_Fock
contains

subroutine ROMP2_Read(ierr,Message)
integer(kind=iwp) :: MPIFlag
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message
  integer(kind=iwp) :: nd(1),ns(1),nf(1),nv(1),nb(1),mult,nel,l
  logical(kind=iwp) :: found,df
  character(len=8) :: method
  ierr=1
  Message='ROMP2 requires C1 high-spin determinant RASSCF with conventional or Cholesky integrals.'
  call Get_cArray('Relax Method',method,8)
  MPIFlag=0
  if ((method(1:6)/='CASSCF').and.(method(1:6)/='RASSCF')) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  call DecideOnCholesky(DoCholesky)
  call DecideOnDF(df)
  MPIFlag=0
  if (df) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) then
    Message='ROMP2 RI/DF integrals are unsupported: use Cholesky or conventional integrals.'
    return
  end if
  call Get_iScalar('nSym',nSym)
  MPIFlag=0
  if (nSym/=1) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  call Get_iArray('nBas',nb,1)
  call Get_iArray('nIsh',nd,1)
  call Get_iArray('nAsh',ns,1)
  call Get_iArray('nFro',nf,1)
  call Get_iArray('nDel',nv,1)
  call Get_iScalar('Multiplicity',mult)
  call Get_iScalar('nActel',nel)
  ! CAS(nOpen,nOpen), maximum spin: exactly one determinant.
  ! Frozen/deleted RASSCF orbitals are excluded from the determinant reference.
  MPIFlag=0
  if ((nf(1)/=0).or.(nv(1)/=0)) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  MPIFlag=0
  if ((ns(1)<0).or.(nd(1)<0).or.(nel/=ns(1)).or.(mult/=ns(1)+1)) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  nBas=nb(1)
  nOrb=nBas
  nOccA=nd(1)+ns(1)
  nOccB=nd(1)
  MPIFlag=0
  if ((nBas<1).or.(nOccA>nBas)) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  MPIFlag=0
  if (nBas>huge(l)/nBas) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  call qpg_dArray('RASSCF orbitals',found,l)
  MPIFlag=0
  if (.not.found) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  MPIFlag=0
  if (l/=nBas*nBas) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  call mma_allocate(CAlpha,nBas,nOrb,label='ROMP2 alpha orbitals')
  call mma_allocate(CBeta,nBas,nOrb,label='ROMP2 beta orbitals')
  call mma_allocate(EOrbA,nOrb,label='ROMP2 alpha energies')
  call mma_allocate(EOrbB,nOrb,label='ROMP2 beta energies')
  call Get_dArray('RASSCF orbitals',CAlpha,nBas*nBas)
  CBeta(:,:)=CAlpha(:,:)
  call Get_dScalar('Last energy',ESCF)
  MPIFlag=0
  if (logical( &
      (.not.MBPT2_All_Finite(CAlpha)).or.(.not.MBPT2_All_Finite(ESCF)), &
      kind=kind(.true.))) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  ReferenceMethod='ROHF    '
  ierr=0
  Message=''
end subroutine ROMP2_Read

subroutine ROMP2_Fock(FA,FB,ierr,Message)
  use OneDat, only: AuxOne,sNoNuc,sNoOri
  use TwoDat, only: AuxTwo
integer(kind=iwp) :: MPIFlag
  real(kind=wp), intent(out) :: FA(nBas,nBas),FB(nBas,nBas)
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message
  real(kind=wp), allocatable :: ERI(:,:,:,:),H(:,:),S(:,:),DA(:,:),DB(:,:),Buf(:),T(:,:),MA(:,:),MB(:,:)
  real(kind=wp) :: required,enuc,echeck,grad,oa,ob,pa,pb
  integer(kind=iwp) :: avail,integer_limit,np,rc,rc2,lu,comp,sym,opt,q,nmat,fs,fbasis(8),skip(8)
  logical(kind=iwp) :: exists,SquareFlag
  logical :: opened,oneopened
  character(len=8) :: label
  ierr=1
  Message='ROMP2 Fock construction failed.'
  opened=.false.
  oneopened=.false.
  required=14.0_wp*real(nBas,wp)**2+131072.0_wp
  if (.not.DoCholesky) required=required+real(nBas,wp)**4
  call mma_maxDBLE(avail)
  integer_limit=huge(avail)
  integer_limit=ishft(integer_limit,-1)
  MPIFlag=0
  if ((required>0.8_wp*real(avail,wp)).or.(required>real(integer_limit,wp))) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) then
    Message='Insufficient memory for ROMP2 reference Fock construction.'
    return
  end if
  np=nBas*(nBas+1)/2
  if (.not.DoCholesky) then
    call f_Inquire('ORDINT',exists)
    MPIFlag=0
    if (.not.exists) MPIFlag=1
    if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
    if (MPIFlag/=0) then
      Message='Conventional ROMP2 requires stored ORDINT: use SEWARD NoCholesky.'
      return
    end if
    MPIFlag=0
    if (AuxTwo%Opn) MPIFlag=1
    if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
    if (MPIFlag/=0) return
    call mma_allocate(ERI,nBas,nBas,nBas,nBas,label='ROMP2 AO ERI')
  end if
  call mma_allocate(Buf,np+4,label='ROMP2 integral buffer')
  call mma_allocate(H,nBas,nBas,label='ROMP2 core Hamiltonian')
  call mma_allocate(S,nBas,nBas,label='ROMP2 overlap')
  call mma_allocate(DA,nBas,nBas,label='ROMP2 density alpha')
  call mma_allocate(DB,nBas,nBas,label='ROMP2 density beta')
  call mma_allocate(T,nBas,nBas,label='ROMP2 work')
  call mma_allocate(MA,nBas,nBas,label='ROMP2 MO Fock alpha')
  call mma_allocate(MB,nBas,nBas,label='ROMP2 MO Fock beta')
  call BuildFockMatrices()
  if (opened) call ClsOrd(rc2)
  if (oneopened) call ClsOne(rc2,0)
  call mma_deallocate(ERI,safe='*')
  call mma_deallocate(Buf)
  call mma_deallocate(H)
  call mma_deallocate(S)
  call mma_deallocate(DA)
  call mma_deallocate(DB)
  call mma_deallocate(T)
  call mma_deallocate(MA)
  call mma_deallocate(MB)
contains
  subroutine ReadOneMatrices()
    if (.not.AuxOne%Opn) then
      call OpnOne(rc,0,'ONEINT',41)
      oneopened=.true.
      if (rc/=0) return
    end if
    opt=ibset(ibset(0,sNoOri),sNoNuc)
    comp=1
    sym=1
    label='OneHam  '
    call RdOne(rc,opt,label,comp,Buf,sym)
    if (rc/=0) return
    call Square(Buf,H,nBas,1,nBas)
    label='Mltpl  0'
    call RdOne(rc,opt,label,comp,Buf,sym)
    if (rc/=0) return
    call Square(Buf,S,nBas,1,nBas)
    if (oneopened) then
      call ClsOne(rc,0)
      oneopened=.false.
      if (rc/=0) return
    end if
  end subroutine ReadOneMatrices

  subroutine BuildFockMatrices()
  integer(kind=iwp) :: MPIFlag,i,j,k,l
  ! ONEINT and ORDINT I/O are local in the Molcas reader layer.
  ! Root is the sole file reader: distribute complete data explicitly.
  H(:,:)=Zero
  S(:,:)=Zero
  rc=0
  if (((nProcs<=1).or.(.not.Is_Real_Par()).or.(MyRank==0))) call ReadOneMatrices()
  MPIFlag=0
  if (rc/=0) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  if (((nProcs>1).and.Is_Real_Par())) then
    call GADGOp(H(1,1),size(H),'+')
    call GADGOp(S(1,1),size(S),'+')
  end if
  call DGEMM_('N','N',nBas,nBas,nBas,One,S,nBas,CAlpha,nBas,Zero,T,nBas)
  call DGEMM_('T','N',nBas,nBas,nBas,One,CAlpha,nBas,T,nBas,Zero,MA,nBas)
  do i=1,nBas
    MA(i,i)=MA(i,i)-One
  end do
  MPIFlag=0
  if (logical( &
      (.not.MBPT2_All_Finite(MA)).or.(maxval(abs(MA))>1.0e-7_wp), &
      kind=kind(.true.))) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) then
    Message='ROHF orbitals are not orthonormal in the ONEINT overlap metric.'
    return
  end if
  if (.not.DoCholesky) then
    ! In real parallel execution ORDINT is distributed across ranks.
    ! Every rank reads its local contribution and the row is reconstructed
    ! by a global sum before it is unpacked into the replicated AO tensor.
    rc=0
    lu=43
    call OpnOrd(rc,0,'ORDINT',lu)
    opened=.true.
    if (rc==0) then
      call GetOrd(rc,SquareFlag,fs,fbasis,skip)
      if (rc==0) then
        if ((fs/=1).or.(fbasis(1)/=nBas).or.(skip(1)/=0)) rc=1
      end if
    end if
    MPIFlag=0
    if (rc/=0) MPIFlag=1
    if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
    if (MPIFlag/=0) return
    opt=1
    do i=1,nBas
      do j=1,i
        Buf(:)=Zero
        rc=0
        call RdOrd_(rc,opt,1,1,1,1,Buf,np+1,nmat)
        opt=2
        if (rc==0) then
          if (nmat/=1) rc=1
          if (.not.MBPT2_All_Finite(Buf(1:np))) rc=1
        end if
        MPIFlag=0
        if (rc/=0) MPIFlag=1
        if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
        if (MPIFlag/=0) return
        if (((nProcs>1).and.Is_Real_Par())) call GADGOp(Buf(1),np,'+')
        q=0
        do k=1,nBas
          do l=1,k
            q=q+1
            ERI(i,j,k,l)=Buf(q)
            ERI(j,i,k,l)=Buf(q)
            ERI(i,j,l,k)=Buf(q)
            ERI(j,i,l,k)=Buf(q)
          end do
        end do
      end do
    end do
    rc=0
    if (opened) then
      call ClsOrd(rc)
      opened=.false.
    end if
    MPIFlag=0
    if (rc/=0) MPIFlag=1
    if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
    if (MPIFlag/=0) return
  end if
  call DGEMM_('N','T',nBas,nBas,nOccA,One,CAlpha,nBas,CAlpha,nBas,Zero,DA,nBas)
  call DGEMM_('N','T',nBas,nBas,nOccB,One,CBeta,nBas,CBeta,nBas,Zero,DB,nBas)
  FA(:,:)=H
  FB(:,:)=H
  if (DoCholesky) then
    call ROMP2_Cho_Fock(DA,DB,FA,FB,rc,Message)
    MPIFlag=0
    if (rc/=0) MPIFlag=1
    if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
    if (MPIFlag/=0) return
  else
    ! Include H once, then distribute matrix-element contractions.
    if (.not.((nProcs<=1).or.(.not.Is_Real_Par()).or.(MyRank==0))) then
      FA(:,:)=Zero
      FB(:,:)=Zero
    end if
    do j=1,nBas
      do i=1,nBas
        if (((nProcs>1).and.Is_Real_Par())) then
          if (mod((j-1)*nBas+i-1,nProcs)/=MyRank) cycle
        end if
        do l=1,nBas
          do k=1,nBas
            FA(i,j)=FA(i,j)+(DA(k,l)+DB(k,l))*ERI(i,j,k,l)-DA(k,l)*ERI(i,k,j,l)
            FB(i,j)=FB(i,j)+(DA(k,l)+DB(k,l))*ERI(i,j,k,l)-DB(k,l)*ERI(i,k,j,l)
          end do
        end do
      end do
    end do
    if (((nProcs>1).and.Is_Real_Par())) then
      call GADGOp(FA(1,1),size(FA),'+')
      call GADGOp(FB(1,1),size(FB),'+')
    end if
  end if
  MPIFlag=0
  if (logical( &
      (.not.MBPT2_All_Finite(FA)).or.(.not.MBPT2_All_Finite(FB)), &
      kind=kind(.true.))) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) return
  call Get_dScalar('PotNuc',enuc)
  echeck=enuc+Half*(sum(DA*(H+FA))+sum(DB*(H+FB)))
  MPIFlag=0
  if (logical( &
      (.not.MBPT2_All_Finite(echeck)).or.(abs(echeck-ESCF)>1.0e-7_wp), &
      kind=kind(.true.))) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) then
    if (((nProcs<=1).or.(.not.Is_Real_Par()).or.(MyRank==0))) then
      write(u6,'(6X,A,T48,F20.10,A)') 'Stored ROHF reference energy',ESCF,' a.u.'
      write(u6,'(6X,A,T48,F20.10,A)') 'Reconstructed ROHF energy',echeck,' a.u.'
    end if
    Message='ROHF determinant energy does not match RASSCF: check reference and Hamiltonian.'
    return
  end if
  call DGEMM_('N','N',nBas,nBas,nBas,One,FA,nBas,CAlpha,nBas,Zero,T,nBas)
  call DGEMM_('T','N',nBas,nBas,nBas,One,CAlpha,nBas,T,nBas,Zero,MA,nBas)
  call DGEMM_('N','N',nBas,nBas,nBas,One,FB,nBas,CAlpha,nBas,Zero,T,nBas)
  call DGEMM_('T','N',nBas,nBas,nBas,One,CAlpha,nBas,T,nBas,Zero,MB,nBas)
  grad=Zero
  do i=1,nBas
    oa=Zero
    ob=Zero
    if (i<=nOccA) oa=One
    if (i<=nOccB) ob=One
    do j=1,nBas
      pa=Zero
      pb=Zero
      if (j<=nOccA) pa=One
      if (j<=nOccB) pb=One
      grad=max(grad,abs((oa-pa)*MA(i,j)+(ob-pb)*MB(i,j)))
    end do
  end do
  MPIFlag=0
  if (grad>1.0e-5_wp) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) then
    if (((nProcs<=1).or.(.not.Is_Real_Par()).or.(MyRank==0))) &
      write(u6,'(6X,A,T48,ES14.6)') 'ROHF stationarity residual',grad
    Message='ROHF orbital stationarity check failed: tighten RASSCF convergence.'
    return
  end if
  ierr=0
  Message=''
  end subroutine BuildFockMatrices

end subroutine ROMP2_Fock
end module ROMP2_Reference
