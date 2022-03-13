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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!  OverLq
!
!> @brief
!>   Compute overlap between primitve bases, which for bases of other type than s,
!>   will mean that several overlaps between basis-functions are computed. One
!>   function is on the solvent, the other in the QM-region. Observe that we do not
!>   care about overlaps within the QM-region or among the solvent molecules.
!> @author A. Ohrn
!>
!> @details
!> Uses the formulas in \cite Tak1966-JPSJ-21-2313. It is hard to give any
!> easy explanation, so if you want to understand exactly what is
!> going on below, see the article, especially equations (2.4) and
!> (2.12); then the source-code comments will provide you with
!> sufficient information.
!>
!> @param[in]  Bori   Center for the QM-region contracted basis-function
!> @param[in]  Cori   Like Bori, but for the solvent basis-function
!> @param[in]  Alfa   Exponents for the primitive basis-functions that build this contracted function
!> @param[in]  Beta   Like \p alfa, but for solvent
!> @param[in]  iQ1    = ``1`` if s-type, = ``2`` if p-type, etc. for the function in the QM-region
!> @param[in]  iQ2    Like \p iQ1, but for solvent function
!> @param[in]  nExp1  How many primitives there are in this contracted function
!> @param[in]  nExp2  Like \p nExp1, but for (surprise) the solvent
!> @param[out] iPSint Pointer to the matrix of overlaps
!> @param[in]  Trans  Transition matrix between Cartesian and spherical basis functions
!***********************************************************************

subroutine OverLq(Bori,Cori,Alfa,Beta,iQ1,iQ2,nExp1,nExp2,iPSint,Trans)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One, Half, Pi
use Definitions, only: wp, iwp

implicit none
! MaxAngqNr=4 means f-function is top. There is no limit in the algorithm though, so if higher is needed, change this number.
#include "maxi.fh"
#include "WrkSpc.fh"
real(kind=wp) :: Bori(3), Cori(3), Alfa(MxCont), Beta(MxCont), &
                 Trans(int(real(3*MxAngqNr**2-2*MxAngqNr-10+8*MxAngqNr**3+3*MxAngqNr**4)/real(12))) !Why this size? is it correct?
integer(kind=iwp) :: iQ1, iQ2, nExp1, nExp2, iPSint
integer(kind=iwp) :: i, icompo, ind, ind1, ind2, iP1, iP2, iPInte, iPpS, iPsphS, iSp1, iSp2, iUpX, iUpY, iUpZ, ix, ixxx, iy, iyyy, &
                     iz, izzz, j, jndex, Kaunt, kaunter, krakna, loneX, loneY, loneZ, lsumX, lsumY, lsumZ, ltwoX, ltwoY, ltwoZ, &
                     nBigP, nCartxC(nTri_Elem(MxAngqNr)), nCartxQ(nTri_Elem(MxAngqNr)), nCartyC(nTri_Elem(MxAngqNr)), & !IFG
                     nCartyQ(nTri_Elem(MxAngqNr)), nCartzC(nTri_Elem(MxAngqNr)), nCartzQ(nTri_Elem(MxAngqNr)), nSizeCart, & !IFG
                     nSizeSph, nSpecific1, nSpecific2, nSph1, nSph2
real(kind=wp) :: Divide, Expo, Extra, FactorX(2*MxAngqNr+1), FactorY(2*MxAngqNr+1), FactorZ(2*MxAngqNr+1), PAxyz(3), PBxyz(3), & !IFG
                 Piconst, Primequals, Separation, SqPiconst, SummaX, SummaY, SummaZ, TheCent(3), TheFirstFac
integer(kind=iwp), external :: iDubFac

!----------------------------------------------------------------------*
! Prepare some numbers for later.                                      *
!----------------------------------------------------------------------*
ind = 0
! Remember that each base consist of many functions - this is how many.
nSpecific1 = iQ1*(iQ1+1)/2
nSpecific2 = iQ2*(iQ2+1)/2
icompo = nSpecific1
! Then follows the NEW loops that compute how the cartesian components are ordered in various basis functions.
do ix=0,iQ1-1
  do iy=0,iQ1-1-ix
    iz = iQ1-1-iy-ix
    ncartxQ(icompo) = ix
    ncartyQ(icompo) = iy
    ncartzQ(icompo) = iz
    icompo = icompo-1
  end do
end do
icompo = nSpecific2
! And for the solvent orbitals.
do ix=0,iQ2-1
  do iy=0,iQ2-1-ix
    iz = iQ2-1-iy-ix
    ncartxC(icompo) = ix
    ncartyC(icompo) = iy
    ncartzC(icompo) = iz
    icompo = icompo-1
  end do
end do
nSph1 = 2*iQ1-1
nSph2 = 2*iQ2-1
nSizeCart = nSpecific1*nSpecific2
nSizeSph = nSph1*nSph2
nBigP = nExp1*nExp2*nSph1*nSph2
call GetMem('PrimCar','Allo','Real',iPpS,nSizeCart)
call GetMem('PrimSph','Allo','Real',iPsphS,nSizeSph)
call GetMem('AllPrims','Allo','Real',iPSint,nBigP)
do i=1,nBigP
  Work(iPSint+i-1) = 0
end do
Separation = ((Bori(1)-Cori(1))**2+(Bori(2)-Cori(2))**2+(Bori(3)-Cori(3))**2)
!----------------------------------------------------------------------*
! Start loop over primitives.                                          *
!----------------------------------------------------------------------*
Kaunt = 0
do iP1=1,nExp1
  do iP2=1,nExp2
    Kaunt = Kaunt+1
    TheCent(1) = (Alfa(iP1)*Bori(1)+Beta(iP2)*Cori(1))
    TheCent(2) = (Alfa(iP1)*Bori(2)+Beta(iP2)*Cori(2))
    TheCent(3) = (Alfa(iP1)*Bori(3)+Beta(iP2)*Cori(3))
    Divide = 1/(Alfa(iP1)+Beta(iP2)) !gamma in article
    TheCent(1) = TheCent(1)*Divide !The new center, P
    TheCent(2) = TheCent(2)*Divide
    TheCent(3) = TheCent(3)*Divide
    Piconst = Pi*Divide
    SqPiconst = sqrt(Piconst)
    Piconst = Piconst*SqPiconst !That constant to the power of 3/2
    Expo = Alfa(iP1)*Beta(iP2)*Separation*Divide
    TheFirstFac = Piconst*exp(-Expo) !This is the exponential factor
    ! Now we should get those difficult f-functions.
    PAxyz(1) = TheCent(1)-Bori(1)
    PAxyz(2) = TheCent(2)-Bori(2)
    PAxyz(3) = TheCent(3)-Bori(3)
    PBxyz(1) = TheCent(1)-Cori(1)
    PBxyz(2) = TheCent(2)-Cori(2)
    PBxyz(3) = TheCent(3)-Cori(3)
    kaunter = 0
    do iSp1=1,nSpecific1
      do iSp2=1,nSpecific2
        kaunter = kaunter+1
        loneX = ncartxQ(iSp1)
        ltwoX = ncartxC(iSp2)
        lsumX = loneX+ltwoX
        loneY = ncartyQ(iSp1)
        ltwoY = ncartyC(iSp2)
        lsumY = loneY+ltwoY
        loneZ = ncartzQ(iSp1)
        ltwoZ = ncartzC(iSp2)
        lsumZ = loneZ+ltwoZ
        call fFactor(loneX,ltwoX,lsumX,loneY,ltwoY,lsumY,loneZ,ltwoZ,lsumZ,PAxyz,PBxyz,FactorX,FactorY,FactorZ)
        ! Now we have the f-factors for this specific angular type of this
        ! specific primitive basis-function. Now put things together.
        iUpX = lsumX/2 !Yes, it should be like this, even when lsumX is odd.
        iUpY = lsumY/2
        iUpZ = lsumZ/2
        SummaX = 0
        SummaY = 0
        SummaZ = 0
        ! This is just a matter of putting things together according to the formula
        do ixxx=0,iUpX
          Extra = iDubFac(2*ixxx-1)*(Half*Divide)**ixxx
          SummaX = SummaX+FactorX(2*ixxx+1)*Extra
        end do
        do iyyy=0,iUpY
          Extra = iDubFac(2*iyyy-1)*(Half*Divide)**iyyy
          SummaY = SummaY+FactorY(2*iyyy+1)*Extra
        end do
        do izzz=0,iUpZ
          Extra = iDubFac(2*izzz-1)*(Half*Divide)**izzz
          SummaZ = SummaZ+FactorZ(2*izzz+1)*Extra
        end do
        Primequals = TheFirstFac*SummaX*SummaY*SummaZ
        Work(iPpS+kaunter-1) = Primequals
      end do
    end do
    !------------------------------------------------------------------*
    ! This was the overlap for the primitives in terms of cartesian    *
    ! functions, but in the new qmstat we use spherical functions, so  *
    ! we need to transform if any d-function or higher is involved. In *
    ! the matrix Trans the numbers for how spherical functions are     *
    ! expressed in cartesian functions are stored, including the extra *
    ! normalization so that all d-functions (and higher) have the same *
    ! combined contraction and normalization coefficient. The rest is  *
    ! just a matter of getting the matrix multiplications right. The   *
    ! convention I use is this: The matrix with the overlaps contains  *
    ! elements such as <psi_QM|psi_Solv> in other words, the           *
    ! QM-orbitals count over the rows and the solvent orbitals over the*
    ! columns; observe however that this is NOT the way the matrix     *
    ! enters from above, since there the fastest couting index is over *
    ! solvent orbitals (iSp2), so given this and the knowledge of how  *
    ! Fortran stores multidimensional matrices, we can figure out when *
    ! to transpose. All this means that if it is the QM-orbitals that  *
    ! are to be transformed, the transformation matrix is multiplied   *
    ! from left, while it is multiplied from the right -- transposed of*
    ! course -- if it is the solvent orbitals that are to be           *
    ! transformed. In the case that both orbitals are to be transformed*
    ! we simply apply the transformation matrix from both directions.  *
    !------------------------------------------------------------------*
    if ((iQ1 >= 3) .or. (iQ2 >= 3)) then !Check if any transformations are necessary.
      if (iQ2 < 3) then !If only the base of the QM-region needs!to be transformed.
        ind = 1+(iQ1-3)*(3*iQ1**3+5*iQ1**2+12*iQ1+40)/12
        call Dgemm_('N','T',nSph1,nSpecific2,nSpecific1,One,Trans(ind),nSph1,Work(iPps),nSpecific2,Zero,Work(iPsphS),nSph1)
      else if (iQ1 < 3) then !If only solvent base needs to be transformed.
        ind = 1+(iQ2-3)*(3*iQ2**3+5*iQ2**2+12*iQ2+40)/12
        call Dgemm_('T','T',nSph1,nSph2,nSpecific2,One,Work(iPps),nSpecific2,Trans(ind),nSph2,Zero,Work(iPsphS),nSph1)
      else !Both QM-region and Solvent need to be transformed.
        ind1 = 1+(iQ1-3)*(3*iQ1**3+5*iQ1**2+12*iQ1+40)/12
        ind2 = 1+(iQ2-3)*(3*iQ2**3+5*iQ2**2+12*iQ2+40)/12
        call GetMem('Intmd','Allo','Real',iPInte,nSph1*nSpecific2)
        call Dgemm_('N','T',nSph1,nSpecific2,nSpecific1,One,Trans(ind1),nSph1,Work(iPps),nSpecific2,Zero,Work(iPInte),nSph1)
        call Dgemm_('N','T',nSph1,nSph2,nSpecific2,One,Work(iPInte),nSph1,Trans(ind2),nSph2,Zero,Work(iPsphS),nSph1)
        call GetMem('Intmd','Free','Real',iPInte,nSph1*nSpecific2)
      end if
    else  !Here we only transpose to get integrals in right order.
      krakna = 0
      do i=0,nSph1-1
        do j=0,nSph2-1
          Work(iPsphS+krakna) = Work(iPps+i+j*nSph1)
          krakna = krakna+1
        end do
      end do
    end if
    !------------------------------------------------------------------*
    ! Put this thing in the slowly growing overlap matrix for the      *
    ! primitive basis functions. The reason the index is so nasty is   *
    ! that we compute small blocks of the matrix and now have to fit it*
    ! in the right place in the growing, much larger, matrix. Nasty!   *
    !------------------------------------------------------------------*
    krakna = 0
    do j=0,nSph2-1
      do i=0,nSph1-1
        jndex = i+j*nExp1*nSph1+(iP1-1)*nSph1+(iP2-1)*nExp1*nSph1*nSph2
        Work(iPSint+jndex) = Work(iPsphS+krakna)
        krakna = krakna+1
      end do
    end do
  end do
end do
!----------------------------------------------------------------------*
! Deallocate and ta'ta!                                                *
!----------------------------------------------------------------------*
call GetMem('PrimCar','Free','Real',iPpS,nSizeCart)
call GetMem('PrimSph','Free','Real',iPsphS,nSizeSph)

return

end subroutine OverLq
