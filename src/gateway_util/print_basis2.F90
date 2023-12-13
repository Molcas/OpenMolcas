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

subroutine Print_Basis2()

use Basis_Info, only: dbsc, iCnttp_Dummy, nCnttp, Shells
use Center_Info, only: dc
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use define_af, only: AngTp
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
#include "print.fh"
integer(kind=iwp) :: i, iAddr, iAng, iAngl, ib, iBas, iBas_Aux, iBas_Frag, iBass, ic, icnt, iCnttp, iExp, iPrim, iPrim_Aux, &
                     iPrim_Frag, iPrimm, iPrint, ir, iRout, irow, iSh, iShSrt, jExp, jSh, kCmp, kExp, kSh, kShEnd, kShStr, lSh, &
                     mdc, nBasisj, ncr, nExpi, nExpj, nExpk, nSumA, nSumB
logical(kind=iwp) :: lAux, lECP, lFAIEMP, lPam2, lPP, output
real(kind=wp) :: ccr, zcr

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)

lAux = .false.
lPam2 = .false.
lECP = .false.
lPP = .false.
lFAIEMP = .false.
do i=1,nCnttp
  lAux = lAux .or. dbsc(i)%Aux
  lPAM2 = lPAM2 .or. dbsc(i)%lPAM2
  lECP = lECP .or. dbsc(i)%ECP
  lPP = lPP .or. (dbsc(i)%nPP /= 0)
  lFAIEMP = lFAIEMP .or. dbsc(i)%Frag
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 6) then
  write(u6,*)
  call CollapseOutput(1,'   Primitive basis info:')
  write(u6,'(3X,A)') '   ---------------------'
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of primitive basis functions

if (iPrint >= 6) then
  write(u6,*)
  write(u6,'(19X,A)') ' *****************************************************'
  write(u6,'(19X,A)') ' ******** Primitive Basis Functions (Valence) ********'
  write(u6,'(19X,A)') ' *****************************************************'
end if
! Loop over distinct shell types
jExp = 0
iPrim = 0
iPrim_Aux = -1
iPrim_Frag = 0
iBas = 0
iBas_Aux = -1
iBas_Frag = 0
! Loop over basis sets
do iCnttp=1,nCnttp
  mdc = dbsc(iCnttp)%mdci
  output = iPrint >= 6
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) output = output .and. (iPrint >= 10) .and. (iCnttp /= iCnttp_Dummy)
  if (output) then
    write(u6,*)
    write(u6,*)
    write(u6,'(A,A)') ' Basis set:',dbsc(iCnttp)%Bsl
  end if
  iShSrt = dbsc(iCnttp)%iVal
  ! Loop over distinct centers
  do icnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    if (mdc > MxAtom) then
      call WarningMessage(2,'MxAtom too small')
      write(u6,*) 'MxAtom=',MxAtom
      write(u6,*) 'Increase MxAtom in Molcas.fh and recompile the code!'
      call Abend()
    end if
    ! Loop over shells associated with this center
    ! Start with s type shells
    jSh = iShSrt
    do iAng=0,dbsc(iCnttp)%nVal-1
      nExpj = Shells(jSh)%nExp
      nBasisj = Shells(jSh)%nBasis
      if ((S%MaxPrm(iAng) > 0) .and. (nExpj > 0) .and. (nBasisj > 0) .and. output .and. (iCnt == 1)) then
        write(u6,*)
        write(u6,*) '                 Type         '
        write(u6,'(19X,A)') AngTp(iAng)
        write(u6,*) '          No.      Exponent    Contraction Coefficients'
      end if

      if ((nBasisj > 0) .and. output) then
        do kExp=1,nExpj
          jExp = jExp+1
          if (iCnt == 1) write(u6,100) jExp,Shells(jSh)%Exp(kExp),(Shells(jSh)%Cff_c(kExp,ib,2),ib=1,nBasisj)
        end do
      end if
      kCmp = (iAng+1)*(iAng+2)/2
      if (Shells(jSh)%Prjct) kCmp = 2*iAng+1
      if (nBasisj /= 0) then
        if (Shells(jSh)%Aux) then
          iPrim_Aux = iPrim_Aux+nExpj*kCmp*nIrrep/dc(mdc)%nStab
          iBas_Aux = iBas_Aux+nBasisj*kCmp*nIrrep/dc(mdc)%nStab
        else if (Shells(jSh)%Frag) then
          iPrim_Frag = iPrim_Frag+nExpj*kCmp*nIrrep/dc(mdc)%nStab
          iBas_Frag = iBas_Frag+nBasisj*kCmp*nIrrep/dc(mdc)%nStab
        else
          iPrim = iPrim+nExpj*kCmp*nIrrep/dc(mdc)%nStab
          iBas = iBas+nBasisj*kCmp*nIrrep/dc(mdc)%nStab
        end if
      end if
      jSh = jSh+1
    end do
  end do

end do
if (iBas >= 2*MaxBfn) then
  call WarningMessage(2,'MaxBfn too small')
  write(u6,*) 'Input: Increase 2*MaxBfn to ',iBas
  call Abend()
end if
if (iPrint >= 6) then
  write(u6,*)
  write(u6,*) ' Number of primitives                ',iPrim
  write(u6,*) ' Number of basis functions           ',iBas
  if (lAux .and. (iPrint >= 10)) then
    write(u6,*) ' Number of primitive aux. functions  ',iPrim_Aux
    write(u6,*) ' Number of auxiliary basis functions ',iBas_Aux
  end if
  if (lFAIEMP .and. (iPrint >= 10)) then
    write(u6,*) ' Number of primitive frag. functions ',iPrim_Frag
    write(u6,*) ' Number of fragment basis functions  ',iBas_Frag
  end if
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (lPAM2) then
  write(u6,*)
  write(u6,'(19X,A)') ' *************************************************'
  write(u6,'(19X,A)') ' ******** Primitive Basis Functions (PAM) ********'
  write(u6,'(19X,A)') ' *************************************************'
  do iCnttp=1,nCnttp
    if (dbsc(iCnttp)%lPAM2) then
      !if (iPrint >= 10) then
      write(u6,*)
      write(u6,*)
      write(u6,'(A,A)') ' Basis set:',dbsc(iCnttp)%Bsl
      if (dbsc(iCnttp)%nPAM2 /= -1) then
        write(u6,*)
        write(u6,*) 'Angular momentum of PAM operator: ',AngTp(dbsc(iCnttp)%nPAM2)
        iAddr = 1

        do iAngl=0,dbsc(iCnttp)%nPAM2
          iPrimm = int(dbsc(iCnttp)%PAM2(iAddr))
          iBass = int(dbsc(iCnttp)%PAM2(iAddr+1))
          write(u6,'(A,3x,a1)') ' Ang. moment: ',AngTp(iAngl)
          write(u6,'(A,i4,A,i4)') ' Number of  primitive:',iPrimm,' Number of contracted:',iBass
          write(u6,*)
          write(u6,'(A)') '  N.        Exponents            Coefficent:'

          do ir=0,iPrimm-1
            write(u6,200) ir+1,dbsc(iCnttp)%PAM2(iAddr+2+ir),(dbsc(iCnttp)%PAM2(iAddr+2+iPrimm+ic),ic=ir,(iPrimm)*iBass-1,iPrimm)
          end do

          iAddr = iAddr+2+iPrimm+iPrimm*iBass
        end do
      end if
      !end if
    end if
  end do

end if
!                                                                      *
!***********************************************************************
!                                                                      *
if ((iPrint >= 10) .and. (lECP .or. lPP)) then
  write(u6,*)
  write(u6,'(19X,A)') ' *************************************************'
  write(u6,'(19X,A)') ' ******** Primitive Basis Functions (ECP) ********'
  write(u6,'(19X,A)') ' *************************************************'
  !else
  !  ! start Molcas
  !  if (iPrint < 6) write(u6,*) 'To display basis set information use the key "BSSHOW" in the input.'
  !  if ((lECP .or. lPP) .and. (iPrint < 10)) write(u6,*) 'To display ECP information use the key "ECPSHOW" in the input.'
  !  if (lAUX .and. (iPrint < 10)) write(u6,*) 'To display auxiliary basis information use the key "AUXSHOW" in the input.'
  !end
end if

do iCnttp=1,nCnttp

  ! Pseudo potential type ECP

  if ((dbsc(iCnttp)%nPP /= 0) .and. (iPrint >= 10)) then
    write(u6,*)
    write(u6,*)
    write(u6,'(A,A)') ' Basis set:',dbsc(iCnttp)%Bsl
    write(u6,*)
    write(u6,'(A)') ' Pseudo Potential'
    write(u6,*)
    kShStr = dbsc(iCnttp)%iPP
    kShEnd = kShStr+dbsc(iCnttp)%nPP-1
    lSh = 0
    do kSh=kShStr,kShEnd
      nExpk = Shells(kSh)%nExp
      if (nExpk /= 0) then
        if (lSh == 0) then
          write(u6,'(4X,A)') '  H Potential'
        else
          write(u6,'(4X,A)') AngTp(lSh-1)//'-H Potential'
        end if
      end if
      lSh = lSh+1
      write(u6,'(A)') '  n     Exponent      Coefficient'
      do iExp=1,nExpk,3
        ncr = int(Shells(kSh)%Exp(iExp))
        zcr = Shells(kSh)%Exp(iExp+1)
        ccr = Shells(kSh)%Exp(iExp+2)
        write(u6,'(2x,I1,3X,2F15.10)') ncr,zcr,ccr
      end do
      write(u6,*)

    end do
  end if

  ! Huzinaga type ECP

  if (dbsc(iCnttp)%ECP) then
    if (iPrint >= 10) then
      write(u6,*)
      write(u6,*)
      write(u6,'(A,A)') ' Basis set:',dbsc(iCnttp)%Bsl

      if (dbsc(iCnttp)%nM1 /= 0) then
        write(u6,*)
        write(u6,*) ' M1 operator       Exponent    Contraction Coefficients'
        do irow=1,dbsc(iCnttp)%nM1
          write(u6,'(14X,ES16.9,1X,ES19.9)') dbsc(iCnttp)%M1xp(irow),dbsc(iCnttp)%M1cf(irow)
        end do
      end if ! if (dbsc(iCnttp)%nM1 /= 0) then

      if (dbsc(iCnttp)%nM2 /= 0) then
        write(u6,*)
        write(u6,*) ' M2 operator       Exponent    Contraction Coefficients'
        do irow=1,dbsc(iCnttp)%nM2
          write(u6,'(14X,ES16.9,1X,ES19.9)') dbsc(iCnttp)%M2xp(irow),dbsc(iCnttp)%M2cf(irow)
        end do
      end if ! if (dbsc(iCnttp)%nM2 /= 0) then
    end if ! if (iPrint >= 10) then

    ! Projection Basis Set

    iSh = dbsc(iCnttp)%iPrj
    nSumB = 0
    jSh = iSh
    do iAng=0,dbsc(iCnttp)%nPrj-1
      nSumB = nSumB+Shells(jSh)%nBasis
      jSh = jSh+1
    end do
    if ((nSumB /= 0) .and. (iPrint >= 10)) then
      write(u6,*)
      write(u6,*)
      write(u6,*) ' Proj. Operator'
    end if
    do iAng=0,dbsc(iCnttp)%nPrj-1
      if (Shells(iSh)%nBk /= 0) then
        if (iPrint >= 10) then
          write(u6,*)
          write(u6,'(19X,A,A)') '        Angular Type: ',AngTp(iAng)
          write(u6,*) '                   Exponent    Contraction Coefficients'
          write(u6,*)
          write(u6,'(A,18X,8(G12.5),/,5(32X,8(G12.5),/))') '     Bk-values',(Shells(iSh)%Bk(i),i=1,Shells(iSh)%nBk)
          write(u6,'(A,18X,8(G12.5),/,5(32X,8(G12.5),/))') '     Frac.Occ.',(Shells(iSh)%Occ(i),i=1,Shells(iSh)%nBk)
        end if

        do i=1,Shells(iSh)%nBk
          Shells(iSh)%Bk(i) = Shells(iSh)%Bk(i)*Shells(iSh)%Occ(i)
        end do

        if (iPrint >= 10) then
          do kExp=1,Shells(iSh)%nExp
            jExp = jExp+1
            write(u6,300) Shells(ish)%Exp(kExp),(Shells(ish)%Cff_c(kExp,ib,2),ib=1,Shells(iSh)%nBk)
          end do
        end if ! if (iPrint >= 10) then
      end if ! if (Shells(iSh)%nBk /= 0) then
      iSh = iSh+1
    end do ! iAng

    ! Auxiliary core basis

    if (iPrint >= 10) then
      iSh = dbsc(iCnttp)%iSOC
      nSumB = 0
      jSh = iSh
      do iAng=0,dbsc(iCnttp)%nSOC-1
        nSumB = nSumB+Shells(jSh)%nBasis
        jSh = jSh+1
      end do
      if ((nSumB /= 0) .and. (iPrint >= 10)) then
        write(u6,*)
        write(u6,*)
        write(u6,*) ' SOC Basis'
      end if
      do iAng=0,dbsc(iCnttp)%nSOC-1
        if (Shells(iSh)%nBasis /= 0) then
          write(u6,*)
          write(u6,'(19X,A,A)') '        Angular Type: ',AngTp(iAng)
          write(u6,*) '                   Exponent    Contraction Coefficients'
          write(u6,*)
          do kExp=1,Shells(iSh)%nExp
            jExp = jExp+1
            write(u6,400) Shells(iSh)%Exp(kExp),(Shells(iSh)%Cff_c(kExp,ib,1),ib=1,Shells(iSh)%nBasis)
          end do
        end if
        iSh = iSh+1
      end do
    end if ! if (iPrint >= 10) then

    ! Spectral Resolution Basis Set

    if (iPrint >= 10) then
      iSh = dbsc(iCnttp)%iSRO
      nSumA = 0
      jSh = iSh
      do iAng=0,dbsc(iCnttp)%nSRO-1
        nSumA = nSumA+Shells(jSh)%nExp
        jSh = jSh+1
      end do
      if (nSumA /= 0) then
        write(u6,*)
        write(u6,*)
        write(u6,*) ' Spectral Resolution Basis Set'
      end if
      do iAng=0,dbsc(iCnttp)%nSRO-1
        nExpi = Shells(iSh)%nExp
        if (nExpi /= 0) then
          write(u6,*)
          write(u6,'(19X,A,A)') '        Angular Type: ',AngTp(iAng)
          call RecPrt(' Exponents',' ',Shells(iSh)%Exp,nExpi,1)
          if (iPrint >= 11) then
            call RecPrt(' The Akl matrix','(5ES20.13)',Shells(iSh)%Akl(1,1,1),nExpi,nExpi)
            call RecPrt(' The Adl matrix','(5ES20.13)',Shells(iSh)%Akl(1,1,2),nExpi,nExpi)
          end if
        end if
        iSh = iSh+1
      end do ! iAng
    end if ! if (iPrint >= 10) then
  end if ! if (dbsc(iCnttp)%ECP) then

end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 6) then
  call CollapseOutput(0,'   Primitive basis info:')
  write(u6,*)
end if

return

100 format(9X,I4,1X,ES16.9,10(1X,F10.6),1X,3(/,30X,10(1X,F10.6)))
200 format(i4,1x,f18.12,1x,12(1x,f12.8))
300 format(14X,ES16.9,8(G12.5),3(/,30X,8(G12.5)))
400 format(14X,ES16.9,10(1X,F10.6),3(/,30X,10(1X,F10.6)))

end subroutine Print_Basis2
