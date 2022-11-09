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
! Copyright (C) 2021, Rulin Feng                                       *
!***********************************************************************

subroutine XDR_Prop(nbas,isize,jsize,imethod,paratyp,xorder,inS,inK,inV,inpVp,inX,inpXp,inUL,inUS,clight,Label,iComp,iSizec)
! Driver for relativistic transformation of property integrals
!
! also called the "picture change correction" since the physical operator
! X is defined in four-component picture, one need to transform them
! to two-/one-component picture as well as the Hamiltonian in the two-/one-
! component relativistic calculations

use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nbas, isize, jsize, imethod, paratyp, xorder, iComp
real(kind=wp), intent(in) :: inS(isize), inK(isize), inV(isize), inpVp(isize), inpXp(isize), clight
real(kind=wp), intent(inout) :: inX(isize), inUL(jsize), inUS(jsize)
integer(kind=iwp) :: i, j, k, iSizec, idbg, nSym, iOpt, iRC, Lu_One, lOper, n_Int, iCmp, jComp, iPSOComp, IDUM(1)
character(len=8) :: Label, magLabel, PSOLabel
real(kind=wp), allocatable :: sK(:,:), sS(:,:), sV(:,:), spVp(:,:), sX(:,:), spXp(:,:), tmp(:,:), magaPX(:), magaPXs(:,:), &
                              magaXP(:), magaXPs(:,:), magbPX(:), magbPXs(:,:), magbXP(:), magbXPs(:,:), PSO(:,:), PSOt(:), Ppso(:)

! Convert to square matrices

call mma_allocate(sK,nbas,nbas,label='skin')
call mma_allocate(sS,nbas,nbas,label='sSS')
call mma_allocate(sV,nbas,nbas,label='sV')
call mma_allocate(spVp,nbas,nbas,label='spVp')
call mma_allocate(sX,nbas,nbas,label='sX')
call mma_allocate(spXp,nbas,nbas,label='spXp')
#ifdef MOLPRO
call square(inK,sK,nbas,nbas)
call square(inS,sS,nbas,nbas)
call square(inV,sV,nbas,nbas)
call square(inpVp,spVp,nbas,nbas)
call square(inX,sX,nbas,nbas)
call square(inpXp,spXp,nbas,nbas)
#else
call square(inK,sK,nbas,1,nbas)
call square(inS,sS,nbas,1,nbas)
call square(inV,sV,nbas,1,nbas)
call square(inpVp,spVp,nbas,1,nbas)
call square(inX,sX,nbas,1,nbas)
call square(inpXp,spXp,nbas,1,nbas)
#endif

! Calculate the relativistic transformed property integrals

if ((imethod == 2) .or. (imethod == 3) .or. ((imethod == 1) .and. (xorder >= 15))) then

  ! Handle magnetic integrals if provided by Gen1Int
  !
  !********************************************************************
  !
  ! PSO integrals with i, pre-existing PSO integrals in ONEREL
  ! required. The h_UL{\dag} should have a minus sign, for
  ! which a tranpose is not enough, thus add this manually.
  ! * MAG:
  ! *  1: x_k*d/dx  4: x_k*d/dy  7: x_k*d/dz
  ! *  2: y_k*d/dx  5: y_k*d/dy  8: y_k*d/dz
  ! *  3: z_k*d/dx  6: z_k*d/dy  9: z_k*d/dz
  ! * PSO1: operator (y_k*d/dz - z_k*d/dy)/r_k^3
  ! *  = MAG 8 - MAG 6
  ! * PSO2: operator (z_k*d/dx - x_k*d/dz)/r_k^3
  ! *  = MAG 3 - MAG 7
  ! * PSO3: operator (x_k*d/dy - y_k*d/dx)/r_k^3
  ! *  = MAG 4 - MAG 2
  ! * PSO n = MAG a - MAG b

  if ((Label(1:3) == 'MAG') .and. ((iComp == 3) .or. (iComp == 4) .or. (iComp == 8))) then
    nSym = 1 ! always C1 symmetry
    inUS(:) = inUS(:)*clight
    iOpt = 0
    iRC = -1
    Lu_One = 2
    call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
    if (iRC /= 0) call Error()
    ! switch basis to primitive functions
    call OneBas('PRIM')
    ! do MAG a
    write(magLabel,'(A,A3)') 'MAGXP',Label(6:8)
    iOpt = ibset(0,sOpSiz)
    iRC = -1
    lOper = -1
    iCmp = iComp
    call iRdOne(iRC,iOpt,magLabel,iCmp,idum,lOper)
    n_Int = IDUM(1)
    if (iRC /= 0) call Error()
    call mma_allocate(magaXP,n_Int+4,label='MAGaXP')
    iOpt = 0
    iRC = -1
    call RdOne(iRC,iOpt,magLabel,iCmp,magaXP,lOper)
    if (iRC /= 0) call Error()
    !call CmpInt(magaXP,n_Int,nbas,nSym,lOper)
    call mma_allocate(magaPX,n_Int+4,label='MAGaPX')
    magLabel(1:5) = 'MAGPX'
    call RdOne(iRC,iOpt,magLabel,iCmp,magaPX,lOper)
    !if (iRC /= 0) call Error()
    call mma_allocate(magaXPs,nbas,nbas,label='MAGaXPs')
    call mma_allocate(magaPXs,nbas,nbas,label='MAGaPXs')
    call square(magaXP,magaXPs,nbas,1,nbas)
    call square(magaPX,magaPXs,nbas,1,nbas)
    call mma_deallocate(magaXP)
    call mma_deallocate(magaPX)
    call merge_mag_ints(nbas,jsize,magaXPs,magaPXs,.true.)
    ! put magaPXs into -magaPXs
    ! or put magaXPs into -magaXPs
    ! test which
    magaPXs(:,:) = -magaPXs(:,:)
    call mma_allocate(tmp,nbas,nbas,label='TMP')
    ! rulin- disable or reenable dmxma for mag integral X2C
    ! transformations
    call dmxma(nbas,'N','N',magaPXs,inUS,tmp,One)
    call dmxma(nbas,'T','N',inUL,tmp,magaPXs,One)
    call dmxma(nbas,'N','N',magaXPs,inUL,tmp,One)
    call dmxma(nbas,'T','N',inUS,tmp,magaXPs,One)
    magaXPs(:,:) = magaXPs(:,:)+magaPXs(:,:)
    call mma_deallocate(magaPXs)
    ! define MAG b
    if (iComp == 3) then
      jComp = 7
      iPSOComp = 2
    else if (iComp == 4) then
      jComp = 2
      iPSOComp = 3
    else if (iComp == 8) then
      jComp = 6
      iPSOComp = 1
    end if
    ! do MAG b
    write(magLabel,'(A,A3)') 'MAGXP',Label(6:8)
    iOpt = ibset(0,sOpSiz)
    iRC = -1
    lOper = -1
    call iRdOne(iRC,iOpt,magLabel,jComp,idum,lOper)
    n_Int = IDUM(1)
    if (iRC /= 0) call Error()
    call mma_allocate(magbXP,n_Int+4,label='MAGbXP')
    iOpt = 0
    iRC = -1
    call RdOne(iRC,iOpt,magLabel,jComp,magbXP,lOper)
    call mma_allocate(magbPX,n_Int+4,label='MAGbPX')
    magLabel(1:5) = 'MAGPX'
    call RdOne(iRC,iOpt,magLabel,jComp,magbPX,lOper)
    call ClsOne(iRC,iOpt)
    call mma_allocate(magbXPs,nbas,nbas,label='MAGbXPs')
    call mma_allocate(magbPXs,nbas,nbas,label='MAGbPXs')
    call square(magbXP,magbXPs,nbas,1,nbas)
    call square(magbPX,magbPXs,nbas,1,nbas)
    call mma_deallocate(magbXP)
    call mma_deallocate(magbPX)
    call merge_mag_ints(nbas,jsize,magbXPs,magbPXs,.true.)
    ! put magbPXs into -magbPXs
    ! or put magbXP into -magbXP
    ! test which
    magbPXs(:,:) = -magbPXs(:,:)
    call dmxma(nbas,'N','N',magbPXs,inUS,tmp,One)
    call dmxma(nbas,'T','N',inUL,tmp,magbPXs,One)
    call dmxma(nbas,'N','N',magbXPs,inUL,tmp,One)
    call dmxma(nbas,'T','N',inUS,tmp,magbXPs,One)
    magbXPs(:,:) = magbXPs(:,:)+magbPXs(:,:)
    call mma_deallocate(magbPXs)
    call mma_deallocate(tmp)
    ! do PSO combinations
    call mma_allocate(PSO,nbas,nbas,label='PSO')
    PSO(:,:) = magaXPs(:,:)-magbXPs(:,:)
    call mma_deallocate(magaXPs)
    call mma_deallocate(magbXPs)
    write(PSOLabel,'(A,A3)') 'PSOI ',Label(6:8)
    ! test for off-diagonal elements
    IDUM(1) = nbas
    call CmpInt(PSO,n_Int,idum(1),nSym,lOper)
    ! store PSO integrals to ONEINT
    call mma_allocate(PSOt,n_Int+4,label='PSOt')
    k = 1
    do i=1,nbas
      do j=1,i
        PSOt(k) = PSO(j,i)
        k = k+1
      end do
    end do
    call mma_deallocate(PSO)
    call mma_allocate(Ppso,iSizec+4,label='PPSO')
    idbg = -1
    call repmat(idbg,PSOt,Ppso,.false.)
    call mma_deallocate(PSOt)
    iOpt = 0
    iRC = -1
    Lu_one = 2
    call OpnOne(iRC,iOpt,'ONEINT',Lu_one)
    if (iRC /= 0) call Error()
    iRC = -1
    lOper = 255
    call WrOne(iRC,iOpt,PSOLabel,iPSOComp,Ppso,lOper)
    if (iRC /= 0) call Error()
    iOpt = 0
    call ClsOne(iRC,iOpt)
    call mma_deallocate(Ppso)
    inUS(:) = inUS(:)/clight
  end if

  if (Label(1:3) == 'MAG') then
    inUS(:) = inUS(:)*clight
    ! Put together lower and upper triangular matrices
    call merge_mag_ints(nbas,jsize,sX,spXp,.true.)
    call mma_allocate(tmp,nbas,nbas,label='TMP')
    ! Eval U_L^{\dag} X U_S
    ! Eval U_S^{\dag} X U_L
    call dmxma(nbas,'N','N',spXp,inUS,tmp,One)
    call dmxma(nbas,'T','N',inUL,tmp,spXp,One)
    call dmxma(nbas,'N','N',sX,inUL,tmp,One)
    call dmxma(nbas,'T','N',inUS,tmp,sX,One)
    ! Sum
    sX(:,:) = sX(:,:)+spXp(:,:)
    ! copy spXp back so its not a half computed matrix
    spXp(:,:) = sX(:,:)
    call mma_deallocate(tmp)
    inUS(:) = inUS(:)/clight
  else

    ! X2C/BSS transformation
    !
    ! because the transformation matrix in non-orthogonal basis picture has
    ! obtained via the Hamiltonian drivers, here we just need to simply apply
    ! it to property integrals ( X, pXp in four-component picture )
    !
    ! high order DKH can also employ this formulation, only negligible
    ! contribution from higher orders is included

    call mma_allocate(tmp,nbas,nbas,label='TMP')
    ! eval U_L^{\dag} X U_L
    call dmxma(nbas,'C','N',inUL,sX,tmp,One)
    call dmxma(nbas,'N','N',tmp,inUL,sX,One)
    ! eval U_S^{\dag}pXp U_S
    call dmxma(nbas,'C','N',inUS,spXp,tmp,One)
    call dmxma(nbas,'N','N',tmp,inUS,spXp,One)
    ! sum
    sX(:,:) = sX(:,:)+spXp(:,:)
    call mma_deallocate(tmp)
  end if !MAG

else if (imethod == 1) then

  ! Arbitrary order DKH transformation

  call dkh_prop(nbas,sS,sK,sV,spVp,sX,spXp,clight,xorder,paratyp)
end if

! Copy transformed property integral back to inX

k = 1
do i=1,nbas
  do j=1,i
    inX(k) = sX(j,i)
    k = k+1
  end do
end do

! Free temp memories

call mma_deallocate(sK)
call mma_deallocate(sS)
call mma_deallocate(sV)
call mma_deallocate(spVp)
call mma_deallocate(sX)
call mma_deallocate(spXp)

return

contains

subroutine Error()
  use Definitions, only: u6
  write(u6,*) ' *** Error in subroutine XDR_Prop ***'
  write(u6,*) '     Abend in subroutine OpnOne'
  call Abend()
end subroutine Error

end subroutine XDR_Prop
