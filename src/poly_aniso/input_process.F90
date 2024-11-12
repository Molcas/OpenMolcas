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
                         iPrint,nsymm,ngrid,R_LG,gtens_input,dirX,dirY,dirZ,dir_weight,zJ,cryst,coord,hmin,hmax,TempMagn,thrs, &
                         Hexp,Mexp,tmin,tmax,chit_exp,Texp,Xfield,Jex,JAex,JAex9,JDMex,tpar,upar,MagnCoords,encut_rate,eso, &
                         JITOexR,JITOexI,Title,itype,namefile_aniso,Do_structure_abc,old_aniso_format,compute_barrier,fitCHI,fitM, &
                         hinput,tinput,compute_magnetization,compute_Mdir_vector,zeeman_energy,m_paranoid,m_accurate,smagn, &
                         compute_g_tensors,compute_torque,compute_susceptibility,Lines,AnisoLines3,AnisoLines9,Dipol,check_title, &
                         KE,DM_exchange,JITO_exchange)
!  nneq      : number of non-equivalent sites
!  neq(nneq) : number of equivalent sites of each type
!  neqv      : = MAXVAL(neq(:))
! definition of the exchange:
!  exch      : total number of exchange states
!  nPair     : number of metal pairs (number of interactions)
!  i_pair    : index of the metal site in a given interacting pair
!  imaxrank  : index of the ITO ranks for each pair
!  Jex       : Lines exchange    ( 1 parameter / interacting pair)
!  JAex      : Anisotropic Lines ( 3 parameter / interacting pair)
!  JAex9     : Anisotropic Lines full ( 9 parameters / interacting pair)
!  JDMex     : Dzyaloshinsky-Morya exchange
!  MxRank1, MxRank2, Lines, AnisoLines3, AnisoLines9, Dipol, DM_exchange, JITO_exchange, old_aniso_format, JITOexR, JITOexI
! options used in connection with KE
!  lant, KEOPT, multLn, KE, tpar, upar
! options used in connection with Dipol-Dipol interaction
!  MagnCoords
! definition of data for susceptibility
!  nT, tinput, compute_susceptibility, tmin, tmax, chit_exp, Texp
! options related to XT_MoverH
!  Xfield, nH, nTempMagn, iopt, TempMagn, Hexp, Mexp, thrs, hmin, hmax, hinput, compute_magnetization, compute_Mdir_vector, &
!  zeeman_energy, m_paranoid, m_accurate, smagn
! options used to set up nM and EM
!  encut_definition
!  nK, mG     : encut_definition=1;
!  ncut       : encut_definition=2;
!  encut_rate : encut_definition=3;
! definition of g and D tensors
!  nMult, nDim, compute_g_tensors
! magnetization torque
!  nP, AngPoints, compute_torque
! Zeeman energy and M vector
!  nDir, nDirZee, dirX, dirY, dirZ, dir_weight
! definition of mean field parameter
!  zJ
! definition of the crystal axes:
!  Do_structure_abc
!  cryst : a, b, c, alpha, beta, gamma
! Cartesian coordinates of the main metal site, or center
!  coord
! definitions for blocking barrier
!  nBlock, compute_barrier
! options for automatic fitting of parameters:
!  fitCHI : not used so far
!  fitM   : not used so far
!  check_title, Title

use Lebedev_quadrature, only: order_table
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nneq, neq(nneq), neqv, nmax, nCenter, nexch(nneq), nDir, nDirZee, nH, nT, exch, nTempMagn, iopt, &
                                 nMult, nDim(nMult), nPair, i_pair(nPair,2), nP, AngPoints, lant, multLn, KEOPT, encut_definition, &
                                 nK, mG, ncut, nss(nneq), nsfs(nneq), nLoc, MxRank1, MxRank2, imaxrank(nPair,2), iPrint, nsymm
integer(kind=iwp), intent(inout) :: nDirTot, nBlock, ngrid
real(kind=wp), intent(in) :: R_LG(nneq,neqv,3,3), gtens_input(3,nneq), dirX(nDir), dirY(nDir), dirZ(nDir), dir_weight(nDirZee,3), &
                             zJ, cryst(6), coord(3), hmin, hmax, thrs, Hexp(nH), Mexp(nH,nTempMagn), tmin, tmax, chit_exp(nT), &
                             Texp(nT), Xfield, Jex(nPair), JAex(nPair,3), JAex9(nPair,3,3), JDMex(nPair,3), tpar, upar, &
                             MagnCoords(nneq,3), encut_rate, eso(nneq,nLoc), &
                             JITOexR(nPair,MxRank1,-MxRank1:MxRank1,MxRank2,-MxRank2:MxRank2), &
                             JITOexI(nPair,MxRank1,-MxRank1:MxRank1,MxRank2,-MxRank2:MxRank2)
real(kind=wp), intent(inout) :: TempMagn(nTempMagn)
character(len=180), intent(in) :: Title, namefile_aniso(nneq)
character, intent(in) :: itype(nneq)
logical(kind=iwp), intent(in) :: Do_structure_abc, old_aniso_format, fitCHI, fitM, hinput, tinput, compute_magnetization, &
                                 compute_Mdir_vector, zeeman_energy, m_accurate, smagn, compute_g_tensors, compute_torque, Lines, &
                                 AnisoLines3, AnisoLines9, Dipol, check_title, KE, DM_exchange, JITO_exchange
logical(kind=iwp), intent(inout) :: compute_barrier, m_paranoid, compute_susceptibility
integer(kind=iwp) :: i, icount_B_sites, iH, iproj1, iproj2, irank1, irank2, iT, j, jEnd, k, l, m, n
logical(kind=iwp) :: ab_initio_all, nosym
character(len=180) fmtline
integer(kind=iwp), parameter :: ngrid_map(32) = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59, &
                                                 62,65]

!-----------------------------------------------------------------------
! print the data from this Subroutine:
if (check_title) write(u6,'(A,A72)') 'TITL :         = ',Title
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
#ifdef _DEBUGPRINT_
write(u6,'(A,i3)') 'input_process:  icount_B_sites=',icount_B_sites
#endif

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

if (old_aniso_format) write(u6,'(A,A)') 'OLDA :         = ',' Input data files are given in OLD format'

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
ngrid = ngrid_map(ngrid)
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
    write(u6,'(A, I4)') '       nPoints = ',order_table(nsymm,ngrid)
    !-----------------------------------------!
    nDirTot = nDir+nDirZee+order_table(nsymm,ngrid)
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
if ((fitCHI .or. fitM) .and. (tinput .or. hinput)) &
  write(u6,'(A)') 'Automatic fitting of exchange parameters is yet in the development'
! ======================================================================

return

end subroutine input_process
