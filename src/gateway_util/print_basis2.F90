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

use Basis_Info
use Center_Info
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "angtp.fh"
#include "print.fh"
logical output
logical lAux, lPam2, lECP, lPP, lFAIEMP

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)

LuWr = 6
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
  write(LuWr,*)
  call CollapseOutput(1,'   Primitive basis info:')
  write(LuWr,'(3X,A)') '   ---------------------'
  write(LuWr,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of primitive basis functions

if (iPrint >= 6) then
  write(LuWr,*)
  write(LuWr,'(19X,A)') ' *****************************************************'
  write(LuWr,'(19X,A)') ' ******** Primitive Basis Functions (Valence) ********'
  write(LuWr,'(19X,A)') ' *****************************************************'
end if
! Loop over distinct shell types
jExp = 0
iPrim = 0
iPrim_Aux = -1
iPrim_Frag = 0
iBas = 0
iBas_Aux = -1
iBas_Frag = 0
iShell = 0
! Loop over basis sets
do iCnttp=1,nCnttp
  mdc = dbsc(iCnttp)%mdci
  output = iPrint >= 6
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) output = output .and. (iPrint >= 10) .and. (iCnttp /= iCnttp_Dummy)
  if (output) then
    write(LuWr,*)
    write(LuWr,*)
    write(LuWr,'(A,A)') ' Basis set:',dbsc(iCnttp)%Bsl
  end if
  iShSrt = dbsc(iCnttp)%iVal
  ! Loop over distinct centers
  do icnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    if (mdc > MxAtom) then
      call WarningMessage(2,'MxAtom too small')
      write(LuWr,*) 'MxAtom=',MxAtom
      write(LuWr,*) 'Increase MxAtom in Molcas.fh and recompile the code!'
      call Abend()
    end if
    ! Loop over shells associated with this center
    ! Start with s type shells
    jSh = iShSrt
    do iAng=0,dbsc(iCnttp)%nVal-1
      iShell = iShell+1
      nExpj = Shells(jSh)%nExp
      nBasisj = Shells(jSh)%nBasis
      if ((S%MaxPrm(iAng) > 0) .and. (nExpj > 0) .and. (nBasisj > 0) .and. output .and. (iCnt == 1)) then
        write(LuWr,*)
        write(LuWr,*) '                 Type         '
        write(LuWr,'(19X,A)') AngTp(iAng)
        write(LuWr,*) '          No.      Exponent    Contraction Coefficients'
      end if

      if ((nBasisj > 0) .and. output) then
        do kExp=1,nExpj
          jExp = jExp+1
          if (iCnt == 1) write(LuWr,100) jExp,Shells(jSh)%exp(kExp),(Shells(jSh)%Cff_c(kExp,ib,2),ib=1,nBasisj)
        end do
      end if
      if (iShell > MxShll) then
        call WarningMessage(2,'iShell > MxShll')
        write(LuWr,*) ' Change MxShll in Molcas.fh and recompile the code!'
        call Abend()
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
  write(LuWr,*) 'Input: Increase 2*MaxBfn to ',iBas
  call Abend()
end if
if (iPrint >= 6) then
  write(LuWr,*)
  write(LuWr,*) ' Number of primitives                ',iPrim
  write(LuWr,*) ' Number of basis functions           ',iBas
  if (lAux .and. (iPrint >= 10)) then
    write(LuWr,*) ' Number of primitive aux. functions  ',iPrim_Aux
    write(LuWr,*) ' Number of auxiliary basis functions ',iBas_Aux
  end if
  if (lFAIEMP .and. (iPrint >= 10)) then
    write(LuWr,*) ' Number of primitive frag. functions ',iPrim_Frag
    write(LuWr,*) ' Number of fragment basis functions  ',iBas_Frag
  end if
  write(LuWr,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (lPAM2) then
  write(LuWr,*)
  write(LuWr,'(19X,A)') ' *************************************************'
  write(LuWr,'(19X,A)') ' ******** Primitive Basis Functions (PAM) ********'
  write(LuWr,'(19X,A)') ' *************************************************'
  do iCnttp=1,nCnttp
    if (dbsc(iCnttp)%lPAM2) then
      !if (iPrint >= 10) then
      write(LuWr,*)
      write(LuWr,*)
      write(LuWr,'(A,A)') ' Basis set:',dbsc(iCnttp)%Bsl
      if (dbsc(iCnttp)%nPAM2 /= -1) then
        write(LuWr,*)
        write(LuWr,*) 'Angular momentum of PAM operator: ',AngTp(dbsc(iCnttp)%nPAM2)
        iAddr = 1

        do iAngl=0,dbsc(iCnttp)%nPAM2
          iPrimm = int(dbsc(iCnttp)%PAM2(iAddr))
          iBass = int(dbsc(iCnttp)%PAM2(iAddr+1))
          write(LuWr,'(A,3x,a1)') ' Ang. moment: ',AngTp(iAngl)
          write(LuWr,'(A,i4,A,i4)') ' Number of  primitive:',iPrimm,' Number of contracted:',iBass
          write(LuWr,*)
          write(LuWr,'(A)') '  N.        Exponents            Coefficent:'

          do ir=0,iPrimm-1
            write(LuWr,200) ir+1,dbsc(iCnttp)%PAM2(iAddr+2+ir),(dbsc(iCnttp)%PAM2(iAddr+2+iPrimm+ic),ic=ir,(iPrimm)*iBass-1,iPrimm)
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
  write(LuWr,*)
  write(LuWr,'(19X,A)') ' *************************************************'
  write(LuWr,'(19X,A)') ' ******** Primitive Basis Functions (ECP) ********'
  write(LuWr,'(19X,A)') ' *************************************************'
  !else
  !  ! start Molcas
  !  if (iPrint < 6) write(LuWr,*) 'To display basis set information use the key "BSSHOW" in the input.'
  !  if ((lECP .or. lPP) .and. (iPrint < 10)) write(LuWr,*) 'To display ECP information use the key "ECPSHOW" in the input.'
  !  if (lAUX .and. (iPrint < 10)) write(LuWr,*) 'To display auxiliary basis information use the key "AUXSHOW" in the input.'
  !end
end if

do iCnttp=1,nCnttp

  ! Pseudo potential type ECP

  if ((dbsc(iCnttp)%nPP /= 0) .and. (iPrint >= 10)) then
    write(LuWr,*)
    write(LuWr,*)
    write(LuWr,'(A,A)') ' Basis set:',dbsc(iCnttp)%Bsl
    write(LuWr,*)
    write(LuWr,'(A)') ' Pseudo Potential'
    write(LuWr,*)
    kShStr = dbsc(iCnttp)%iPP
    kShEnd = kShStr+dbsc(iCnttp)%nPP-1
    lSh = 0
    do kSh=kShStr,kShEnd
      nExpk = Shells(kSh)%nExp
      if (nExpk /= 0) then
        if (lSh == 0) then
          write(LuWr,'(4X,A)') '  H Potential'
        else
          write(LuWr,'(4X,A)') AngTp(lSh-1)//'-H Potential'
        end if
      end if
      lSh = lSh+1
      write(LuWr,'(A)') '  n     Exponent      Coefficient'
      do iExp=1,nExpk,3
        ncr = int(Shells(kSh)%exp(iExp))
        zcr = Shells(kSh)%exp(iExp+1)
        ccr = Shells(kSh)%exp(iExp+2)
        write(LuWr,'(2x,I1,3X,2F15.10)') ncr,zcr,ccr
      end do
      write(LuWr,*)

    end do
  end if

  ! Huzinaga type ECP

  if (dbsc(iCnttp)%ECP) then
    if (iPrint >= 10) then
      write(LuWr,*)
      write(LuWr,*)
      write(LuWr,'(A,A)') ' Basis set:',dbsc(iCnttp)%Bsl

      if (dbsc(iCnttp)%nM1 /= 0) then
        write(LuWr,*)
        write(LuWr,*) ' M1 operator       Exponent    Contraction Coefficients'
        do irow=1,dbsc(iCnttp)%nM1
          write(LuWr,'(14X,D16.9,1X,D19.9)') dbsc(iCnttp)%M1xp(irow),dbsc(iCnttp)%M1cf(irow)
        end do
      end if ! if (dbsc(iCnttp)%nM1 /= 0) then

      if (dbsc(iCnttp)%nM2 /= 0) then
        write(LuWr,*)
        write(LuWr,*) ' M2 operator       Exponent    Contraction Coefficients'
        do irow=1,dbsc(iCnttp)%nM2
          write(LuWr,'(14X,D16.9,1X,D19.9)') dbsc(iCnttp)%M2xp(irow),dbsc(iCnttp)%M2cf(irow)
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
      write(LuWr,*)
      write(LuWr,*)
      write(LuWr,*) ' Proj. Operator'
    end if
    do iAng=0,dbsc(iCnttp)%nPrj-1
      if (Shells(iSh)%nBk /= 0) then
        if (iPrint >= 10) then
          write(LuWr,*)
          write(LuWr,'(19X,A,A)') '        Angular Type: ',AngTp(iAng)
          write(LuWr,*) '                   Exponent    Contraction Coefficients'
          write(LuWr,*)
          write(LuWr,'(A,18X,8(G12.5),/,5(32X,8(G12.5),/))') '     Bk-values',(Shells(iSh)%Bk(i),i=1,Shells(iSh)%nBk)
          write(LuWr,'(A,18X,8(G12.5),/,5(32X,8(G12.5),/))') '     Frac.Occ.',(Shells(iSh)%Occ(i),i=1,Shells(iSh)%nBk)
        end if

        do i=1,Shells(iSh)%nBk
          Shells(iSh)%Bk(i) = Shells(iSh)%Bk(i)*Shells(iSh)%Occ(i)
        end do

        if (iPrint >= 10) then
          do kExp=1,Shells(iSh)%nExp
            jExp = jExp+1
            write(LuWr,300) Shells(ish)%exp(kExp),(Shells(ish)%Cff_c(kExp,ib,2),ib=1,Shells(iSh)%nBk)
          end do
        end if ! if (iPrint >= 10) then
      end if ! if (Shells(iSh)%nBk /= 0) then
      iSh = iSh+1
    end do ! iAng

    ! Auxilliary core basis

    if (iPrint >= 10) then
      iSh = dbsc(iCnttp)%iSOC
      nSumB = 0
      jSh = iSh
      do iAng=0,dbsc(iCnttp)%nSOC-1
        nSumB = nSumB+Shells(jSh)%nBasis
        jSh = jSh+1
      end do
      if ((nSumB /= 0) .and. (iPrint >= 10)) then
        write(LuWr,*)
        write(LuWr,*)
        write(LuWr,*) ' SOC Basis'
      end if
      do iAng=0,dbsc(iCnttp)%nSOC-1
        if (Shells(iSh)%nBasis /= 0) then
          write(LuWr,*)
          write(LuWr,'(19X,A,A)') '        Angular Type: ',AngTp(iAng)
          write(LuWr,*) '                   Exponent    Contraction Coefficients'
          write(LuWr,*)
          do kExp=1,Shells(iSh)%nExp
            jExp = jExp+1
            write(LuWr,400) Shells(iSh)%exp(kExp),(Shells(iSh)%Cff_c(kExp,ib,1),ib=1,Shells(iSh)%nBasis)
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
        write(LuWr,*)
        write(LuWr,*)
        write(LuWr,*) ' Spectral Resolution Basis Set'
      end if
      do iAng=0,dbsc(iCnttp)%nSRO-1
        nExpi = Shells(iSh)%nExp
        if (nExpi /= 0) then
          write(LuWr,*)
          write(LuWr,'(19X,A,A)') '        Angular Type: ',AngTp(iAng)
          call RecPrt(' Exponents',' ',Shells(iSh)%Exp,nExpi,1)
          if (iPrint >= 11) then
            call RecPrt(' The Akl matrix','(5D20.13)',Shells(iSh)%Akl(1,1,1),nExpi,nExpi)
            call RecPrt(' The Adl matrix','(5D20.13)',Shells(iSh)%Akl(1,1,2),nExpi,nExpi)
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
  write(LuWr,*)
end if

return

100 format(9X,I4,1X,D16.9,10(1X,F10.6),1X,3(/,30X,10(1X,F10.6)))
200 format(i4,1x,f18.12,1x,12(1x,f12.8))
300 format(14X,D16.9,8(G12.5),3(/,30X,8(G12.5)))
400 format(14X,D16.9,10(1X,F10.6),3(/,30X,10(1X,F10.6)))

end subroutine Print_Basis2
