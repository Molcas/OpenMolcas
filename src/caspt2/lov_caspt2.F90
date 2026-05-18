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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine Lov_CASPT2(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,NAME,nName,nUniqAt,Thrs,IFQCAN,DoMP2,DoEnv,all_Vir,EMP2,CMO,NCMO)
!***********************************************************************
!                                                                      *
! Purpose:  setup of Localized occupied-virtual CASPT2 (LovCASPT2).    *
!           The CASPT2 correction to the energy will later be computed *
!           only for the "active region" of the molecule.              *
!           The MP2 correction due to the remaining frozen region      *
!           is computed here if DoMP2=.true.                           *
!           If DoEnv=.true. we compute the energy of the environment   *
!           as the total MP2 energy minus the MP2 energy of the        *
!           "active region" of the molecule.                           *
!                                                                      *
! Author:   F. Aquilante  (Geneva, Feb. 2008)                          *
!                                                                      *
!***********************************************************************

use OneDat, only: sNoNuc, sNoOri
use Molcas, only: LenIn, MxAtom, MxBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use definitions, only: iwp, wp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nSym
integer(kind=iwp), intent(inout) :: nFro(nSym), nIsh(nSym), nSsh(nSym), nDel(nSym)
integer(kind=iwp), intent(in) :: nBas(nSym), nAsh(nSym)
integer(kind=iwp), intent(in) :: nNAME
character(Len=LenIn+8), intent(in) :: NAME(nNAME)
integer(kind=iwp), intent(in) :: nUniqAt
real(kind=wp), intent(in) :: Thrs
integer(kind=iwp), intent(inout) :: IFQCAN
logical(kind=iwp), intent(inout) :: DoMP2
logical(kind=iwp), intent(in) :: DoEnv, all_Vir
real(kind=wp), intent(out) :: EMP2
integer(kind=iwp), intent(in) :: NCMO
real(kind=wp), intent(inout) :: CMO(nCMO)
character(Len=LenIn) blank, NamAct(mxAtom)
character(len=8) :: Label
logical(kind=iwp) ortho
real(kind=wp) TrA(8), TrF(8), TrX(8)
integer(kind=iwp) ns_O(8), ns_V(8)
integer(kind=iwp) lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
integer(kind=iwp), allocatable :: nBas_per_Atom(:), nBas_Start(:), D_A(:), D_Vir(:)
real(kind=wp), allocatable :: SQ(:), SLT(:), CMOX(:), Q(:), Z(:)
real(kind=wp), allocatable :: Saa(:), XMO(:), DMat(:), OrbE(:)
real(kind=wp) Dumm, E2_ab, E2_Aonly, STrA, STrF, STrX, Thrd
real(kind=wp), external :: DDot_
integer(kind=iwp) i, iAt, iBat, iCMO, iComp, iDo, ie, ik, iloc, iOff, iopt, ip_X, ip_Y, ipAsh, ipCMO, ipEorb, ipQa, iQ, iQa, &
                  iSkip, iSQ, iSym, iSymLbl, iV, jAt, jBas, jBat, jCMO, jDo, jjCMO, jjZ, jOff, jQ, jZ, kBas, kEOcc, kEVir, kfr, &
                  kOff, kto, l_nBas_per_Atom, l_nBas_Start, lBas, lOff, lsq, ltri, mAsh, mOff, nActa, nAk, nBasT, nBat, nBk, nBmx, &
                  nOA, nOrb, nSQ, nTri, nVV, ipOrbE, nBx

irc = 0
EMP2 = Zero
blank = '   '
iDo = 0
jDo = 0
if (DoEnv .and. DoMP2) then
  call WarningMessage(1,'Both DoEnv and DoMP2 selected.')
  write(u6,'(/,A)') ' DoMP2 will be ignored.'
  DoMP2 = .false.
end if
if (all_Vir .and. DoMP2) then
  call WarningMessage(1,'Both VirAll and DoMP2 selected.')
  write(u6,'(/,A)') ' DoMP2 will be ignored.'
  DoMP2 = .false.
end if
do iSym=1,nSym
  TrA(iSym) = 0
  TrF(iSym) = 0
  TrX(iSym) = 0
end do

!----------------------------------------------------------------------*
!     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
!----------------------------------------------------------------------*

nBasT = 0
ntri = 0
nSQ = 0
nBmx = 0
mAsh = 0
nOrb = 0
do i=1,nSym
  nBasT = nBasT+nBas(i)
  nOrb = nOrb+nFro(i)+nIsh(i)+nAsh(i)+nSsh(i)+nDel(i)
  ntri = ntri+nBas(i)*(nBas(i)+1)/2
  nSQ = nSQ+nBas(i)**2
  nBmx = max(nBmx,nBas(i))
  mAsh = max(mAsh,nAsh(i))
end do
if (nBasT > mxBas) then
  write(u6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
  call Abend()
end if

! nUniqAt = # of symm. unique atoms. Initialize NamAct to blanks.
! ---------------------------------------------------------------

if ((nUniqAt < 1) .or. (nUniqAt > MxAtom)) then
  write(u6,'(A,I9)') 'nUniqAt =',nUniqAt
  call Abend()
end if
do iAt=1,nUniqAt
  NamAct(iAt) = blank
end do

! Allocate and get index arrays for basis functions per atom.
! -----------------------------------------------------------

l_nBas_per_Atom = nUniqAt
l_nBas_Start = nUniqAt
call mma_allocate(nBas_per_Atom,l_nBas_per_Atom,Label='nB/A')
call mma_allocate(nBas_Start,l_nBas_Start,Label='nBStart')

!----------------------------------------------------------------------*
!     Read the overlap matrix                                          *
!----------------------------------------------------------------------*
call mma_allocate(SQ,nSQ,Label='SQ')
call mma_allocate(SLT,nTri,Label='SLT')
isymlbl = 1
iopt = ibset(ibset(0,sNoOri),sNoNuc)
Label = 'Mltpl  0'
iComp = 1
call RdOne(irc,iopt,Label,iComp,SLT,isymlbl)
if (irc /= 0) return
ltri = 1
lsq = 1
do iSym=1,nSym
  call Square(SLT(ltri),SQ(lsq),1,nBas(iSym),nBas(iSym))
  ltri = ltri+nBas(iSym)*(nBAs(iSym)+1)/2
  lsq = lsq+nBas(iSym)**2
end do
call mma_deallocate(SLT)

call mma_allocate(CMOX,2*NCMO,Label='CMOX')
ipCMO = 1+NCMO
! This is not the best solution, but I wanted to avoid having to rewrite
! the indexing code below just to use the CMO array directly
call dcopy_(NCMO,CMO,1,CMOX,1)
call dcopy_(NCMO,CMOX,1,CMOX(ipCMO),1)

!----------------------------------------------------------------------*
!     Compute Mulliken atomic charges of each active orbital           *
!             on each center to define the Active Site                 *
!----------------------------------------------------------------------*
call mma_allocate(Q,nUniqAt*(mAsh+1),Label='Q')
ipQa = 1+nUniqAt*mAsh
call Fzero(Q(ipQa),nUniqAt)
call mma_allocate(Z,nBmx*mAsh,Label='Z')
lBas = 0
iOff = 0
do iSym=1,nSym
  iSQ = 1+iOff
  ipAsh = 1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
  nBx = max(1,nBas(iSym))
  call DGEMM_('N','N',nBas(iSym),nAsh(iSym),nBas(iSym),One,SQ(iSQ),nBx,CMOX(ipAsh),nBx,Zero,Z,nBx)
  jBas = lBas+1
  kBas = lBas+nBas(iSym)
  call BasFun_Atom_Sym(nBas_per_Atom,nBas_Start,Name,jBas,kBas,nUniqAt,.false.)
  do ik=0,nAsh(iSym)-1
    nAk = nUniqAt*ik
    nBk = nBas(iSym)*ik
    jCMO = ipAsh+nBk-1
    jZ = nBk
    do iAt=0,nUniqAt-1
      iBat = nBas_Start(1+iAt)
      jjCMO = jCMO+iBat
      jjZ = jZ+iBat
      nBat = nBas_per_Atom(1+iAt)
      iQ = 1+nAk+iAt
      Q(iQ) = ddot_(nBat,CMOX(jjCMO),1,Z(jjZ),1)
    end do
  end do
  do iAt=0,nUniqAt-1
    jQ = 1+iAt
    iQa = ipQa+iAt
    Q(iQa) = Q(iQa)+ddot_(nAsh(iSym),Q(jQ),nUniqAt,Q(jQ),nUniqAt)
    if (sqrt(Q(iQa)) >= Thrs) then
      jBat = nBas_Start(1+iAt)+lBas
      NamAct(iAt+1) = Name(jBat)(1:LenIn)
    end if
  end do
  lBas = lBas+nBas(iSym)
  iOff = iOff+nBas(iSym)**2
end do
call mma_deallocate(Z)
call mma_deallocate(Q)

! We have now completed the definition of the active site
!----------------------------------------------------------------------*
call mma_allocate(D_A,nUniqAt,Label='D_A')
nActa = 0
do iAt=1,nUniqAt
  if (NamAct(iAt) /= blank) then
    nActa = nActa+1
    D_A(nActa) = iAt
  end if
end do
do iAt=1,nActa
  jAt = D_A(iAt)
  NamAct(iAt) = NamAct(jAt)
end do
do iAt=nActa+1,nUniqAt
  NamAct(iAt) = blank
end do
write(u6,*)
write(u6,'(A,F15.6)') ' Threshold for atom selection: ',Thrs
write(u6,*)
if (nActa /= 0) then
  write(u6,'(A,I3,A)') ' Selected ',nActa,' atoms: '
  write(u6,*)
  write(u6,*) (NamAct(i),i=1,nActa)
  write(u6,*)
elseif ((.not. DoMP2) .and. (.not. DoEnv)) then
  write(u6,'(A,18A4)') ' Selected atoms: *** None *** '
  Go To 2000
else
  write(u6,'(A,18A4)') ' Selected atoms: *** None *** '
end if

call mma_deallocate(D_A)
!----------------------------------------------------------------------*

call mma_allocate(OrbE,4*nOrb,Label='OrbE')
ipOrbE = 1
call Get_darray('RASSCF OrbE',OrbE(ipOrbE),nOrb)
call Compute_Tr_Dab(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,CMOX(ipCMO),nCMO,OrbE(ipOrbE),4*nOrb,TrX)

! MP2 calculation on the whole system (incompatible with DoMP2)
if (DoEnv) &
  call energy_AplusB(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,CMOX(ipCMO:),size(CMOX(ipCMO:)),OrbE(ipOrbE:),size(OrbE(ipOrbE:)),E2_ab)
!----------------------------------------------------------------------*
!     Localize the inactive and virtual orbitals                       *
!                                                                      *
!        1) inactive orbitals ---> cholesky orbitals (orthonormal)     *
!        2) virtual orbitals ---> lin. indep. PAOs (non-orthonormal)   *
!                                                                      *
!----------------------------------------------------------------------*
Thrd = 1.e-06_wp
call mma_allocate(D_vir,nBasT,Label='D_Vir')
call Cho_ov_Loc(irc,Thrd,nSym,nBas,nFro,nIsh,nAsh,nSsh,CMOX(ipCMO:),SQ,D_vir)

if (irc /= 0) then
  write(u6,*) 'Localization failed in LovCASPT2'
  call Abend()
end if

ipEorb = ipOrbE+nOrb
kEOcc = ipEorb+nOrb
kEVir = kEOcc+nOrb
call mma_allocate(Xmo,2*NCMO,Label='XMO')
iCMO = 1+NCMO
call mma_allocate(Saa,nOrb,Label='Saa')
Saa(:) = One

! Inactive orbital selection
!----------------------------------------------------------------------*
iOff = 0
kOff = 0
lOff = 0
mOff = 0
do iSym=1,nSym
  jOff = iOff+nBas(iSym)*nFro(iSym)
  call dcopy_(nBas(iSym)*nIsh(iSym),CMOX(ipCMO+jOff),1,XMO(1+kOff),1)
  call dcopy_(nBas(iSym)*nIsh(iSym),CMOX(1+jOff),1,XMO(iCMO+kOff),1)
  jOff = lOff+nFro(iSym)
  call dcopy_(nIsh(iSym),OrbE(ipOrbE+jOff),1,OrbE(ipEorb+mOff),1)
  iOff = iOff+nBas(iSym)**2
  kOff = kOff+nBas(iSym)*nIsh(iSym)
  lOff = lOff+nBas(iSym)
  mOff = mOff+nIsh(iSym)
end do
ortho = .true.

call get_Orb_select(irc,XMO(iCMO),XMO,OrbE(ipEorb),SQ,Saa,Name,NamAct,nSym,nActa,nIsh,nBas,ortho,Thrs,ns_O)
if (irc /= 0) return
iOff = 0
kOff = 0
do iSym=1,nSym
  lOff = iOff+nBas(iSym)*nFro(iSym)
  do ik=nIsh(iSym),1,-1
    jOff = kOff+nBas(iSym)*(ik-1)
    call dcopy_(nBas(iSym),XMO(iCMO+jOff),1,CMOX(1+lOff),1)
    lOff = lOff+nBas(iSym)
  end do
  iOff = iOff+nBas(iSym)**2
  kOff = kOff+nBas(iSym)*nIsh(iSym)
end do
iloc = 0
loff = 0
do iSym=1,nSym
  do ik=nIsh(iSym),ns_O(iSym)+1,-1
    ie = ipEorb+loff+ik-1
    OrbE(kEOcc+iloc) = OrbE(ie)
    iloc = iloc+1
  end do
  loff = loff+nIsh(iSym)
end do
joff = 0
loff = 0
do iSym=1,nSym
  koff = joff+nFro(iSym)+nIsh(iSym)-ns_O(iSym)
  do ik=0,ns_O(iSym)-1
    ie = ipEorb+loff+ik
    OrbE(ipOrbE+koff+ik) = OrbE(ie)
  end do
  loff = loff+nIsh(iSym)
  joff = joff+nBas(iSym)
end do

if (all_Vir) then

  do iSym=1,nSym
    ns_V(iSym) = nSsh(iSym)
  end do

else

  ! Virtual orbital selection
  !--------------------------------------------------------------------*
  iOff = 0
  kOff = 0
  lOff = 0
  mOff = 0
  do iSym=1,nSym
    jOff = iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
    call dcopy_(nBas(iSym)*nSsh(iSym),CMOX(ipCMO+jOff),1,XMO(1+kOff),1)
    call dcopy_(nBas(iSym)*nSsh(iSym),CMOX(1+jOff),1,XMO(iCMO+kOff),1)
    jOff = lOff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
    call dcopy_(nSsh(iSym),OrbE(ipOrbE+jOff),1,OrbE(ipEorb+mOff),1)
    iOff = iOff+nBas(iSym)**2
    kOff = kOff+nBas(iSym)*nSsh(iSym)
    lOff = lOff+nBas(iSym)
    mOff = mOff+nSsh(iSym)
  end do
  ortho = .false.
  call get_Saa(nSym,nBas,nSsh,SQ,size(SQ),XMO,size(XMO),Saa,size(Saa))

  call get_Vir_select(irc,XMO(iCMO),XMO,OrbE(ipEorb),SQ,Name,NamAct,D_vir,nSym,nActa,nSsh,nBas,ortho,ns_V)
  if (irc /= 0) return
  call mma_deallocate(D_vir)
  iOff = 0
  kOff = 0
  do iSym=1,nSym
    jOff = iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
    call dcopy_(nBas(iSym)*nSsh(iSym),XMO(iCMO+kOff),1,CMOX(1+jOff),1)
    iOff = iOff+nBas(iSym)**2
    kOff = kOff+nBas(iSym)*nSsh(iSym)
  end do
  iloc = 0
  loff = 0
  do iSym=1,nSym
    do ik=ns_V(iSym)+1,nSsh(iSym)
      ie = ipEorb+loff+ik-1
      OrbE(kEVir+iloc) = OrbE(ie)
      iloc = iloc+1
    end do
    loff = loff+nSsh(iSym)
  end do
  joff = 0
  loff = 0
  do iSym=1,nSym
    koff = joff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
    do ik=0,ns_V(iSym)-1
      ie = ipEorb+loff+ik
      OrbE(ipOrbE+koff+ik) = OrbE(ie)
    end do
    joff = joff+nBas(iSym)
    loff = loff+nSsh(iSym)
  end do

end if

! MP2 calculation on the Frozen region
!----------------------------------------------------------------------*
if (DoMP2) then

  iDo = 0
  jDo = 0
  nVV = 0
  nOA = 0
  do iSym=1,nSym  ! setup info
    lnOrb(iSym) = nBas(iSym)
    lnOcc(iSym) = nIsh(iSym)-ns_O(iSym)
    lnFro(iSym) = nFro(iSym)+ns_O(iSym)
    lnDel(iSym) = nDel(iSym)+ns_V(iSym)
    lnVir(iSym) = nSsh(iSym)-ns_V(iSym)
    iDo = max(iDo,lnOcc(iSym))
    jDo = max(jDo,lnVir(iSym))
    nVV = nVV+lnVir(iSym)**2
    nOA = nOA+lnOcc(iSym)
  end do
  if (min(iDo,jDo) == 0) goto 1000

  call mma_allocate(Dmat,nVV+nOA,Label='DMat')
  ip_X = 1
  ip_Y = ip_X+nVV
  DMat(:) = Zero
  call FZero(XMO(iCMO),NCMO)
  iOff = 0
  do iSym=1,nSym
    kfr = 1+iOff+nBas(iSym)*nFro(iSym)
    kto = iCMO+iOff+nBas(iSym)*lnFro(iSym)
    call dcopy_(nBas(iSym)*lnOcc(iSym),CMOX(kfr),1,XMO(kto),1)
    kfr = 1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym)+ns_V(iSym))
    kto = kto+nBas(iSym)*lnOcc(iSym)
    call dcopy_(nBas(iSym)*lnVir(iSym),CMOX(kfr),1,XMO(kto),1)
    iOff = iOff+nBas(iSym)**2
  end do
  call Check_Amp(nSym,lnOcc,lnVir,iSkip)
  if (iSkip > 0) then
    call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,.true.)
    call ChoMP2_Drv(irc,Dumm,XMO(iCMO),OrbE(kEOcc),OrbE(kEVir),DMAT(ip_X),DMAT(ip_Y))
    call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,.false.)
    call ChoMP2_Drv(irc,EMP2,XMO(iCMO),OrbE(kEOcc),OrbE(kEVir),DMAT(ip_X),DMAT(ip_Y))
    if (irc /= 0) then
      write(u6,*) 'Frozen region MP2 failed'
      call Abend()
    end if
    iV = ip_X
    do iSym=1,nSym
      TrF(iSym) = ddot_(lnVir(iSym),DMAT(iV),1+lnVir(iSym),[One],0)
      iV = iV+lnVir(iSym)**2
    end do
  end if
  call mma_deallocate(Dmat)
1000 write(u6,*)

  if (nActa == 0) then
    write(u6,'(A,F18.10)') ' Frozen region MP2 correction: ',EMP2
    write(u6,*)
  end if
end if

!----------------------------------------------------------------------*

call mma_deallocate(Saa)
call mma_deallocate(XMO)

! Update the nFro, nIsh, nSsh, nDel for the Active site CASPT2
do iSym=1,nSym
  nFro(iSym) = nFro(iSym)+nIsh(iSym)-ns_O(iSym)
  nIsh(iSym) = ns_O(iSym)
  nDel(iSym) = nDel(iSym)+nSsh(iSym)-ns_V(iSym)
  nSsh(iSym) = ns_V(iSym)
  iDo = max(iDo,nIsh(iSym))
  jDo = max(jDo,nSsh(iSym))
end do

call Compute_Tr_Dab(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,CMOX,nCMO,OrbE(ipOrbE),4*nOrb,TrA)

write(u6,*) '------------------------------------------------------'
write(u6,*) ' Symm.  Tr(D):  Active        Frozen        Full      '
write(u6,*) '------------------------------------------------------'
STrA = Zero
STrF = Zero
STrX = Zero
do iSym=1,nSym
  if (DoEnv) TrF(iSym) = TrX(iSym) ! just a convention
  write(u6,'(2X,I4,10X,G11.4,3X,G11.4,3X,G11.4)') iSym,TrA(iSym),TrF(iSym),TrX(iSym)
  STrA = STrA+TrA(iSym)
  STrF = STrF+TrF(iSym)
  STrX = STrX+TrX(iSym)
end do
write(u6,*) '------------------------------------------------------'
write(u6,'(A,G11.4,3X,G11.4,3X,G11.4)') '          Sum:  ',STrA,STrF,STrX
write(u6,*) '------------------------------------------------------'
write(u6,*)

if (DoEnv) then
  call energy_AplusB(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,CMOX,size(CMOX),OrbE(ipOrbE:),size(OrbE(ipOrbE:)),E2_Aonly)
  EMP2 = E2_ab-E2_Aonly
  !write(u6,'(A,F18.10)') ' MP2 correction (environment): ',EMP2
  !write(u6,*)
end if

call mma_deallocate(OrbE)

2000 continue
if (min(iDo,jDo) == 0) then
  write(u6,*)
  write(u6,*) ' None of the inactive or virtual orbitals has been'
  write(u6,*) ' assigned to the Active region of the molecule.'
  write(u6,*) ' This is presumably NOT what you want !!!'
  write(u6,*) ' CASPT2 will Stop here. Bye Bye !!'
  write(u6,*)
  call Abend()
end if

if (IFQCAN /= 0) IFQCAN = 0 ! MOs to be recanonicalized on exit
call dcopy_(NCMO,CMOX,1,CMO,1)

call mma_deallocate(CMOX)
call mma_deallocate(SQ)
call mma_deallocate(nBas_per_Atom)
call mma_deallocate(nBas_Start)

end subroutine Lov_CASPT2
