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

subroutine input_process(nneq,neq,neqv,nmax,nCenter,nexch,nDir,nDirZee,nDirTot,nH,nT,exch,nBlock,nTempMagn,iopt,nMult,nDim,nPair, &
                         i_pair,nP,AngPoints,lant,multLn,KEOPT,encut_definition,nK,mG,ncut,nsfs,nss,nLoc,MxRank1,MxRank2,imaxrank, &
                         iPrint,R_LG,gtens_input,dirX,dirY,dirZ,dir_weight,zJ,cryst,coord,hmin,hmax,TempMagn,thrs,Hexp,Mexp,tmin, &
                         tmax,chit_exp,Texp,Xfield,Jex,JAex,JAex9,JDMex,tpar,upar,MagnCoords,encut_rate,eso,JITOexR,JITOexI,Title, &
                         itype,namefile_aniso,Do_structure_abc,old_aniso_format,compute_barrier,fitCHI,fitM,hinput,tinput, &
                         compute_magnetization,compute_Mdir_vector,zeeman_energy,m_paranoid,m_accurate,smagn,compute_g_tensors, &
                         compute_torque,compute_susceptibility,Lines,AnisoLines3,AnisoLines9,Dipol,check_title,KE,DM_exchange, &
                         JITO_exchange)

use Constants, only: Zero, Two
use Definitions, only: u6

implicit none
#include "mgrid.fh"
integer, intent(in) :: iPrint
integer, intent(in) :: nneq ! number of non-equivalent sites
integer, intent(in) :: neq(nneq), neqv ! number of equivalent sites of each type, neqv = MAXVAL(neq(:))
integer, intent(in) :: nexch(nneq), nmax
character(len=1), intent(in) :: itype(nneq)
integer, intent(in) :: nCenter
integer, intent(in) :: nLoc
real(kind=8), intent(in) :: R_LG(nneq,neqv,3,3)
real(kind=8), intent(in) :: gtens_input(3,nneq)
real(kind=8), intent(in) :: eso(nneq,nLoc)
integer, intent(in) :: nss(nneq), nsfs(nneq)
character(len=180), intent(in) :: namefile_aniso(nneq)
! definition of the exchange:
integer, intent(in) :: exch ! total number of exchange states
integer, intent(inout) :: nPair ! number of metal pairs (number of interactions)
integer, intent(inout) :: i_pair(nPair,2) ! index of the metal site in a given interacting pair
integer, intent(in) :: imaxrank(nPair,2) ! index of the ITO ranks for each pair
integer, intent(in) :: MxRank1, MxRank2
logical, intent(in) :: Lines, AnisoLines3, AnisoLines9
logical, intent(in) :: Dipol, DM_exchange, JITO_exchange
logical, intent(in) :: old_aniso_format
real(kind=8), intent(in) :: Jex(nPair) ! Lines exchange    ( 1 parameter / interacting pair)
real(kind=8), intent(in) :: JAex(nPair,3) ! Anisotropic Lines ( 3 parameter / interacting pair)
real(kind=8), intent(in) :: JAex9(nPair,3,3) ! Anisotropic Lines full ( 9 parameters / interacting pair)
real(kind=8), intent(in) :: JDMex(nPair,3) ! Dzyaloshinsky-Morya exchange
real(kind=8), intent(in) :: JITOexR(nPair,MxRank1,-MxRank1:MxRank1,MxRank2,-MxRank2:MxRank2), &
                            JITOexI(nPair,MxRank1,-MxRank1:MxRank1,MxRank2,-MxRank2:MxRank2)
! options used in connection with KE
integer, intent(in) :: lant, KEOPT, multLn
logical, intent(in) :: KE
real(kind=8), intent(in) :: tpar, upar
! options used in connection with Dipol-Dipol interaction
real(kind=8), intent(in) :: MagnCoords(nneq,3)
! definition of data for susceptibility
integer, intent(inout) :: nT
logical, intent(in) :: tinput
logical, intent(inout) :: compute_susceptibility
real(kind=8), intent(inout) :: tmin, tmax
real(kind=8), intent(in) :: chit_exp(nT), Texp(nT)
! options related to XT_MoverH
real(kind=8), intent(in) :: Xfield
integer, intent(in) :: nH
integer, intent(in) :: nTempMagn
integer, intent(in) :: iopt
real(kind=8), intent(inout) :: TempMagn(nTempMagn)
real(kind=8), intent(in) :: Hexp(nH), Mexp(nH,nTempMagn)
real(kind=8), intent(in) :: thrs
real(kind=8), intent(in) :: hmin, hmax
logical, intent(in) :: hinput
logical, intent(in) :: compute_magnetization
logical, intent(in) :: compute_Mdir_vector
logical, intent(in) :: zeeman_energy
logical, intent(inout) :: m_paranoid
logical, intent(in) :: m_accurate
logical, intent(in) :: smagn
! options used to set up nM and EM
integer, intent(in) :: encut_definition
integer, intent(in) :: nK, mG ! encut_definition=1;
integer, intent(in) :: ncut   ! encut_definition=2;
real(kind=8), intent(in) :: encut_rate ! encut_definition=3;
! definition of g and D tensors
integer, intent(in) :: nMult
integer, intent(in) :: nDim(nMult)
logical, intent(in) :: compute_g_tensors
! magnetization torque
integer, intent(in) :: nP
integer, intent(in) :: AngPoints
logical, intent(in) :: compute_torque
! Zeeman energy and M vector
integer, intent(in) :: nDir, nDirZee
real(kind=8), intent(in) :: dirX(nDir), dirY(nDir), dirZ(nDir)
real(kind=8), intent(in) :: dir_weight(nDirZee,3)
! definition of mean field parameter
real(kind=8), intent(in) :: zJ
! definition of the crystal axes:
logical, intent(in) :: Do_structure_abc
! a, b, c, alpha, beta, gamma
real(kind=8), intent(in) :: cryst(6)
! Cartesian coordinates of the main metal site, or center
real(kind=8), intent(in) :: coord(3)
! definitions for blocking barrier
integer, intent(inout) :: nBlock
logical, intent(inout) :: compute_barrier
! options for automatic fitting of parameters:
logical, intent(in) :: fitCHI !-- not used so far
logical, intent(in) :: fitM !-- not used so far
logical, intent(in) :: check_title
character(len=180), intent(in) :: Title
! local variables
character(len=180) fmtline
logical :: nosym
logical :: ab_initio_all
logical :: DBG
integer :: nDirTot
integer :: i, j, k, l, m, n, iT, iH, jEnd, irank1, irank2, iproj1, iproj2
integer :: icount_B_sites

DBG = .false.
!-----------------------------------------------------------------------
! print the data from this Subroutine:
if (check_title) then
  write(u6,'(A,A72)') 'TITL :         = ',Title
end if
! ======================================================================
! INFORMATION about individual magnetic sites
! ======================================================================
write(u6,'(A,  I5)') 'PRLV :         = ',iPrint
write(u6,'(A)') 'Information about individual magnetic sites:'
write(u6,'(A)')
write(u6,'(A,  I5)') 'NNEQ :         = ',nneq
write(u6,'(A,  I5)') 'NEQV :         = ',neqv
write(u6,'(A,50I5)') 'NEQ() :        = ',(neq(i),i=1,nneq)
write(u6,'(A,  I5)') 'nCenter :      = ',nCenter

write(u6,'(A,50I5  )') 'NEXCH():       = ',(nexch(i),i=1,nneq)
write(u6,'(A,50I5  )') 'NMAX:          = ',nmax
write(u6,'(A,50I5  )') 'NSTATE():      = ',(nsfs(i),i=1,nneq)
write(u6,'(A,50I5  )') 'NSS():         = ',(nss(i),i=1,nneq)

do i=1,nneq
  if (itype(i) == 'A') then
    write(u6,'(i3,1x,A15,1x,3i5,A)') i,namefile_aniso(i),nsfs(i),nss(i),nexch(i),' -- computed ab initio'
    write(u6,'(A,i2,A )') 'ESO(',i,' ):'
    write(u6,'(10F10.3)') (ESO(i,j),j=1,nss(i))
  else
    write(u6,'(i3,A15,1x,3i5,A)') i,' - - - - - - - ',nsfs(i),nss(i),nexch(i),' -- generated as isotropic'
  end if
end do

nosym = .true.
icount_B_sites = 0
ab_initio_all = .true.
do i=1,nneq
  if (neq(i) > 1) nosym = .false.
  if (itype(i) == 'B') then
    icount_B_sites = icount_B_sites+1
    ab_initio_all = .false.
  end if
end do
if (DBG) write(u6,'(A,i3)') 'input_process:  icount_B_sites=',icount_B_sites

if (.not. nosym) then
  write(u6,'(A,A)') 'SYMM :         = ',' rotation matrices for equivalent magnetic centers:'
  do i=1,nneq
    write(u6,'(8x,A,i3)') 'center type:',i
    do j=1,neq(i)
      do m=1,3
        write(u6,'(20x,3F11.6)') (R_lg(i,j,m,n),n=1,3)
      end do
    end do
  end do
else
  write(u6,'(A,A )') 'SYMM :         = ',' rotation matrices are "Identity Matrices" for all magnetic centers'
end if

if (ab_initio_all) then
  write(u6,'(17x,A,I2,A)') 'all centers have been computed ab initio and ',nneq,' file(s) of "aniso_i.input" type exist.'
else
  write(u6,'(17x,I2,A,I2,A)') nneq-icount_B_sites,' center(s) have been computed ab initio and',nneq-icount_B_sites, &
                              ' file(s) of "aniso_i.input" type exist.'
  write(u6,'(17x,I2,A,15F7.3)') icount_B_sites,' center(s) have been defined empirical: no ZFS and g_tens = '
  write(u6,'(17x,A,99F10.5)') 'gX = ',(gtens_input(1,i),i=1,nneq)
  write(u6,'(17x,A,99F10.5)') 'gY = ',(gtens_input(2,i),i=1,nneq)
  write(u6,'(17x,A,99F10.5)') 'gZ = ',(gtens_input(3,i),i=1,nneq)
end if

if (old_aniso_format) then
  write(u6,'(A,A)') 'OLDA :         = ',' Input data files are given in OLD format'
end if

! ======================================================================
! INFORMATION about exchange
! ======================================================================
write(u6,'(A)') 'Information about exchange and magnetic couplings'
if (Lines) then
  write(u6,'(A,  I3  )') 'PAIR or LIN1 : = ',npair
  do i=1,npair
    write(u6,'(17x,2I3,2x, F9.5)') i_pair(i,1),i_pair(i,2),Jex(i)
  end do
end if

if (AnisoLines3) then
  write(u6,'(A,  I3  )') 'LIN3 :         = ',npair
  do i=1,npair
    write(u6,'(17x,2I3,2x,3F9.5)') i_pair(i,1),i_pair(i,2),(JAex(i,j),j=1,3)
  end do
end if

if (AnisoLines9) then
  write(u6,'(A,  I3  )') 'LIN9 :         = ',npair
  do i=1,npair
    write(u6,'(17x,2I3,2x,9F9.5)') i_pair(i,1),i_pair(i,2),(JAex9(i,1,j),j=1,3),(JAex9(i,2,j),j=1,3),(JAex9(i,3,j),j=1,3)
  end do
end if

if (DM_exchange) then
  write(u6,'(A,  I3  )') 'DMEX :         = ',npair
  do i=1,npair
    write(u6,'(17x,2I3,2x,3F9.5)') i_pair(i,1),i_pair(i,2),(JDMex(i,j),j=1,3)
  end do
end if

if (Dipol) then
  write(u6,'(A,A )') 'DIPO :         = ',' dipolar coupling for the above exchange pairs is included'
  write(u6,'(A,A )') 'COOR :         = ',' symmetrized cartesian coordinates for each non-equivalent magnetic centers '
  do i=1,nneq
    write(u6,'(17x,i3,3F11.6)') i,(MagnCoords(i,l),l=1,3)
  end do
else
  write(u6,'(A,A )') 'DIPO :         = ',' dipolar coupling for the above exchange pairs is NOT included.'
end if

if (KE) then
  write(u6,'(A,A )') 'LONG :         = ','kinetic exchange coupling is requested:'
  write(u6,'(A,2F14.7)') 't     = ',tpar
  write(u6,'(A,2F14.7)') 'U     = ',upar
  write(u6,'(A,i2    )') 'lant  = ',lant
  write(u6,'(A,i2    )') 'multLn= ',multLn
  write(u6,'(A,i2    )') 'keopt = ',KEOPT
end if

if (JITO_exchange) then
  write(u6,'(A,A )') 'ITOJ :         = ',' Full anisotopic ITO exchange for the above exchange pairs is included (in cm-1)'
  do i=1,npair
    write(u6,'(17x,2I3)') i_pair(i,1),i_pair(i,2)
    do irank1=1,imaxrank(i,1),2
      do iproj1=-irank1,irank1
        do irank2=1,imaxrank(i,2),2
          do iproj2=-irank2,irank2
            write(u6,'(23x,4I3,2x,2ES23.14)') irank1,iproj1,irank2,iproj2,JITOexR(i,irank1,iproj1,irank2,iproj2), &
                                              JITOexI(i,irank1,iproj1,irank2,iproj2)
          end do
        end do
      end do
    end do
  end do
end if
!-----------------------------------------------------------------------
! print Exchange Hamiltonian
! print Exchange Eigenstates
! print decomposition of exchange
!...
! ======================================================================
! Print out of the UBAR:
if (compute_barrier) then
  if (compute_g_tensors) then
    Nblock = 0
    do i=1,nmult
      Nblock = Nblock+ndim(i)
    end do
  else
    write(u6,'(A)') 'UBAR keyword needs to know the sizes of'
    write(u6,'(A)') 'the manifolds forming the barrier'
    write(u6,'(A)') 'Was the MLTP keyword defined?'
    write(u6,'(A)') 'The program will continue, but UBAR is deactivated'
    compute_barrier = .false.
  end if

  if (Nblock > exch) then
    write(u6,'(A)') 'It seems that you have asked for properties of inexistent states.'
    write(u6,'(A)') 'Nblock > Nexch!'
    write(u6,'(A,I5)') 'Nblock = ',nBlock
    write(u6,'(A,I5)') 'Nexch  = ',exch
    write(u6,'(A   )') 'Reset   NBLOCK=EXCH and continue'
    nBlock = exch
  end if
end if

! ======================================================================
! Print out of the g and D tensors:

if (compute_g_tensors .and. (nMult > 0)) then
  write(u6,'(A)') 'Information about g and D tensors:'
  write(u6,'(A,I3)') 'MLTP :         = ',NMULT
  if (nMult <= 20) then
    write(u6,'(A,20I3)') '               = ',(NDIM(i),i=1,NMULT)
  else
    write(u6,'(A,20I3)') '               = ',(NDIM(i),i=1,20)
    do j=21,NMULT,20
      jEnd = min(NMULT,J+19)
      write(u6,'(A,20I3)') '                 ',(NDIM(i),i=j,jEnd)
    end do
  end if
else
  write(u6,'(A)') 'Computation of pseudospin Hamiltonians was not requested'
end if

if (compute_g_tensors .and. (nMult == 0)) then
  write(u6,'(A,i3)') 'nMult =',nMult
  write(u6,'(A   )') 'Computation of pseudospin Hamiltonians is deactivated.'
end if
! ======================================================================

if (Do_structure_abc .and. ((nMult > 0) .or. compute_susceptibility)) then
  write(u6,'(2A)') 'ABCC :         = ','the main axes for the computed tensors are written also'
  write(u6,'(A)') 'in the crystallographic "abc" axes'
  write(u6,'(10x,A,F9.4)') 'a       = ',cryst(1)
  write(u6,'(10x,A,F9.4)') 'b       = ',cryst(2)
  write(u6,'(10x,A,F9.4)') 'c       = ',cryst(3)
  write(u6,'(10x,A,F9.4)') 'alpha   = ',cryst(4)
  write(u6,'(10x,A,F9.4)') 'beta    = ',cryst(5)
  write(u6,'(10x,A,F9.4)') 'gamma   = ',cryst(6)
  write(u6,'(10x,a,3F9.4)') 'coords: = ',(coord(i),i=1,3)
end if

! ======================================================================
! Print out of the SUSCEPTIBILITY
compute_susceptibility = .true.
if (compute_susceptibility) then
  !-----------------------------------------!
  write(u6,'(A)') 'Magnetic susceptibility will be computedin the limit of zero applied magnetic field.'
  !-----------------------------------------!
  write(u6,'(A,A )') 'CHIT :         = ',' molar magnetic susceptibility is computed using'
  !-----------------------------------------!
  if (IOPT == 1) then
    write(u6,'(18x,A)') 'newly-derived formula for total XT and M'
  else if (IOPT == 2) then
    write(u6,'(18x,A)') 'additive formula employed for total XT and M'
  else if (IOPT == 3) then
    write(u6,'(18x,A)') 'old formula employed for total XT and M'
  else
    write(u6,'(18x,A)') 'IOPT parameter out of range.'
    return
  end if
  !-----------------------------------------!
  if (TINPUT) then
    write(u6,'(A)') 'TEXP :         = the experimental temperature interval is provided:'
    do iT=1,nT
      write(u6,'(5X,i4,F11.5,F14.7)') iT,Texp(iT),chit_exp(iT)
    end do
  else
    write(u6,'(A, I4)') 'TINT :      nT = ',nT
    write(u6,'(A,F7.3)') '          Tmin = ',Tmin
    write(u6,'(A,F7.3)') '          Tmax = ',Tmax
  end if
  !-----------------------------------------!
  if (xField /= Zero) then
    write(u6,'(A)') 'Magnetic susceptibility will be computed also by using:'
    write(u6,'(A)') ' X= dM/dH formula and'
    write(u6,'(A)') ' X=  M/H formula (commonly used by experimentalists)'
    write(u6,'(A,F9.5,A)') 'FIELD :          ',Xfield,' Telsa'
  end if
  !-----------------------------------------!
else
  write(u6,'(A)') 'Computation of magnetic susceptibility was not requested.'
end if

! ======================================================================
! Print out of the TORQUE
if (compute_torque) then
  write(u6,'(A)') 'TORQ :         = ',' magnetic torque is computed'
  write(u6,'(A,I4)') 'AngPoints = ',AngPoints
  write(u6,'(A,i4)') 'nP        = ',nP
else
  write(u6,'(A)') 'Computation of magnetic torque was not requested.'
end if

! ======================================================================

! ======================================================================
! Print out of the MAGNETIZATION
if (compute_magnetization) then
  write(u6,'(2A )') 'MAGN :         = ',' molar magnetization is computed'
  !-----------------------------------------!
  if (nH > 0) then
    if (HINPUT) then
      write(u6,'(A)') 'HEXP :         = the experimental field interval is provided:'
      write(fmtline,'(A,i3,A)') '(5x,i4,F11.5,',nTempMagn,'F14.7)'
      do iH=1,nH
        write(u6,fmtline) iH,Hexp(iH),(Mexp(iH,iT),iT=1,nTempMagn)
      end do
    else
      write(u6,'(A, I4)') 'HINT :      nH = ',nH
      write(u6,'(A,F7.3)') '          Hmin = ',Hmin
      write(u6,'(A,F7.3)') '          Hmax = ',Hmax
    end if
  end if
  !-----------------------------------------!
  if (encut_definition == 1) then
    write(u6,'(A, I4)') 'NCUT :         = ',ncut
  else if (encut_definition == 2) then
    write(u6,'(A,2I4)') 'ECUT :         = ',nk,mg
  else if (encut_definition == 3) then
    write(u6,'(A,F7.3)') 'ERAT :         = ',encut_rate
  end if
  !-----------------------------------------!
  write(u6,'(A, I4   )') 'TMAG :         = ',nTempMagn
  ! check for negative T
  do i=1,nTempMagn
    if (TempMagn(i) <= Zero) TempMagn(i) = Two
  end do

  !-----------------------------------------!
  if (nTempMagn <= 10) then
    write(u6,'(6x,A,10F7.3)') 'TempMagn = ',(TempMagn(i),i=1,nTempMagn)
  else
    do i=1,nTempMagn,10
      j = min(nTempMagn,i+9)
      write(u6,'(17x,10F7.3)') (TempMagn(k),k=i,j)
    end do
  end if
  !-----------------------------------------!
  if (zeeman_energy) then
    write(u6,'(2A,I3,A)') 'ZEEM :         = ',' Zeeman splitting for the following ',nDirZee
    write(u6,'(17x,A  )') 'direction(s) of the magnetic field is given in the "zeeman_energy_xxx.output" file(s)'
    do i=1,nDirZee
      write(u6,'(17x,3F11.6)') (dir_weight(i,l),l=1,3)
    end do
  end if
  !-----------------------------------------!
  if (compute_Mdir_vector) then
    write(u6,'(2A,I2,a)') 'MVEC :         = ','Magnetization vector(s) for the following ',nDir, &
                          ' direction(s) of the applied magnetic'
    write(u6,'(17x,A)') 'field will be computed;'
    do i=1,nDir
      write(u6,'(17x,A,I2,A,3F10.6)') 'Dir :',i,' : ',DirX(i),DirY(i),DirZ(i)
    end do
  end if

  if (nH > 0) then
    !-----------------------------------------!
    write(u6,'(A, I4)') 'MAVE :   nsymm = ',nsymm
    write(u6,'(A, I4)') '         ngrid = ',ngrid
    write(u6,'(A, I4)') '       nPoints = ',get_nP(nsymm,ngrid)
    !-----------------------------------------!
    nDirTot = 0
    nDirTot = nDir+nDirZee+get_nP(nsymm,ngrid)
    write(u6,'(A,I4,A)') 'Magnetization will be computed in',nDirTot,' directions.'
    write(u6,'(A,I6,A)') 'There will be ',nDirTot*nH,' calculations in total.'

    !-----------------------------------------------------------|
    if (m_accurate) then
      write(u6,'(2A   )') 'MACC :         = ','Contribution to magnetization from local excited states'
      write(u6,'(17x,A)') 'is computed exactly.'
    else
      write(u6,'(17x,A)') 'is accounted by formula:  M_excited = X_excited * H'
    end if
    !-----------------------------------------------------------|

    !-----------------------------------------------------------|
    if ((exch < 256) .and. (zJ /= Zero)) m_paranoid = .true.

    if (m_paranoid .and. (zJ /= Zero)) then
      write(u6,'(2A)') 'MPAR :         = ','Average spin < S > of the environment is computed'
      write(u6,'(2A)') '                 ','for each temperature point (more accurate)'
    else
      write(u6,'(2A)') 'MPAR :         = ','Average spin < S > of the environment is computed'
      write(u6,'(2A)') '                 ','for the highest temperature point only and used for all computed temperature points.'
    end if
    !-----------------------------------------------------------|

    !-----------------------------------------------------------|
    if (smagn) then
      write(u6,'(2A)') 'SMAG :         = ','Spin-only magnetisation is computed'
    else
      write(u6,'(A)') 'Computation of spin magnetization was not requested.'
    end if
    !-----------------------------------------------------------|

    write(u6,'(A,ES14.6)') 'THRS = ',thrs
  end if
else
  write(u6,'(A)') 'Computation of molar magnetization was not requested.'
end if
! End Print out of the MAGNETIZATION
! ======================================================================
if (compute_susceptibility .or. compute_torque .or. compute_magnetization) then
  if (zJ /= Zero) then
    write(u6,'(A)') 'Mean field intermolecular interaction parameter'
    write(u6,'(A,F12.7)') 'ZJPR :         = ',zJ
  end if
end if
! ======================================================================
if ((fitCHI .or. fitM) .and. (tinput .or. hinput)) then
  write(u6,'(A)') 'Automatic fitting of exchange parameters is yet in the development'
end if
! ======================================================================

return

end subroutine input_process
