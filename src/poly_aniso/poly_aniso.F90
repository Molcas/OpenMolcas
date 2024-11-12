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

subroutine POLY_ANISO_1(nneq,neqv,nmax,exch,nLoc,nCenter,nT,nH,nTempMagn,nDir,nDirZee,nMult,nPair,MxRank1,MxRank2, &
                        old_aniso_format,iReturn)
! definition of the cluster:
!  nneq : number of non-equivalent sites
!  neqv : number of equivalent sites of each type, neqv = MAXVAL(neq(:))
!  nCenter, neq
! definition of the local metal sites
!  nLoc : number of spin-orbit states. nLoc = MAXVAL(nss(:));
!  nss, nsfs, multiplicity, gtens_input, D_fact, EoverD_fact
! definition of the main local axes in general coord system
!  riso, R_LG, R_ROT
! spin orbit energies on individual metal sites
!  eso, eso_au, dipso, s_so, itype, namefile_aniso, ifHDF, old_aniso_format, DoPlot
! definition of the exchange:
!  exch     : total number of exchange states
!  nPair    : number of metal pairs (number of interactions)
!  nmax     : exchange basis, nmax = MAXVAL(nexch(:))
!  i_pair   : index of the metal site in a given interacting pair
!  imaxrank : index of rank ITOs for JITO exchange definition
!  Jex      : Lines exchange    ( 1 parameter / interacting pair)
!  JAex     : Anisotropic Lines ( 3 parameter / interacting pair)
!  JAex9    : Anisotropic Lines full ( 9 parameters / interacting pair)
!  JDMex    : Anisotropic Lines ( 3 parameter / interacting pair)
!  W        : exchange spectrum
!  Z        : exchange eigenstates
!  dipexch  : exchange magnetic moment
!  s_exch   : exchange spin moment
!  MxRank1, MxRank2, nexch, AnisoLines1, AnisoLines3, AnisoLines9, Dipol, DM_exchange, decompose_exchange, JITO_exchange, JITOexR,
!  JITOexI
! options used in connection with KE
!  lant, multLn, KEOPT, KE, tpar, upar
! options used in connection with Dipol-Dipol interaction
!  MagnCoords, nTempMagn
! definition of g and D tensors
!  nMult : multiplicity of each multiplet
!  gtens : gtensor of each multiplet
!  maxes : main axes of each multiplet
!  compute_g_tensors, nDim
! definition of data for susceptibility
!  nT, compute_susceptibility, tinput, tmin, tmax, dltT0, T, XTexp, XT_no_field, chit_exp, Texp, XLM, ZLM, XRM, ZRM
! options related to XT_MoverH
!  Xfield
! definition of data for magnetization:
!  nH, nM, iopt, TempMagn, Hexp, Mexp, thrs, em, hmin, hmax, dltH0, hinput, compute_magnetization, compute_Mdir_vector,
!  zeeman_energy, m_paranoid, m_accurate, smagn
! options used to set up nM and EM
!  encut_definition
!  nK, mG     : encut_definition=1;
!  ncut       : encut_definition=2;
!  encut_rate : encut_definition=3;
! magnetization torque
!  compute_torque, AngPoints, nP
! Zeeman energy and M vector
!  nDir, nDirZee, nDirTot, LuZee, dirX, dirY, dirZ, dir_weight
! definition of mean field parameter
!  zJ
! definition of the crystal axes:
!  cryst(6) : a, b, c, alpha, beta, gamma
!  coord(3) : Cartesian coordinates of the main metal site, or center
!  Do_structure_abc
! definitions for blocking barrier
!  nBlock, compute_barrier
! options for automatic fitting of parameters:
!  fitCHI : not used so far
!  fitM   : not used so far
!  iPrint, imltpl, Ifunct, iReturn, i, i1, i2, j, l, imanifold, check_title, Title, GRAD, fname, LuAniso

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, auTocm
use Definitions, only: wp, iwp, u6, CtoB, ItoB, RtoB

implicit none
integer(kind=iwp), intent(inout) :: nneq, nT, nH, nTempMagn, nDir, nDirZee, nMult, nPair
integer(kind=iwp), intent(in) :: neqv, nmax, exch, nLoc, nCenter, MxRank1, MxRank2
logical(kind=iwp), intent(in) :: old_aniso_format
integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: AngPoints, encut_definition, i, i1, i2, Ifunct, imanifold, imltpl, iopt, iPrint, j, KEOPT, l, lant, LuAniso, &
                     mem, mG, multLn, nBlock, ncut, nDirTot, ngrid, nK, nM, nP, nsymm
real(kind=wp) :: coord(3), cryst(6), dltH0, dltT0, em, encut_rate, gtens(3), hmax, hmin, maxes(3,3), thrs, tmax, tmin, tpar, upar, &
                 Xfield, zJ
logical(kind=iwp) :: AnisoLines1, AnisoLines3, AnisoLines9, check_title, compute_barrier, compute_g_tensors, &
                     compute_magnetization, compute_Mdir_vector, compute_susceptibility, compute_torque, decompose_exchange, &
                     Dipol, DM_exchange, Do_structure_abc, DoPlot, fitCHI, fitM, GRAD, hinput, ifHDF, JITO_exchange, KE, &
                     m_accurate, m_paranoid, smagn, tinput, zeeman_energy
character(len=180) :: fname, namefile_aniso(nneq), Title
character :: itype(nneq)
integer(kind=iwp), allocatable :: i_pair(:,:), imaxrank(:,:), LuZee(:), nDim(:), neq(:), nexch(:), nsfs(:), nss(:)
real(kind=wp), allocatable :: chit_exp(:), D_fact(:), dir_weight(:,:), dirX(:), dirY(:), dirZ(:), EoverD_fact(:), eso(:,:), &
                              eso_au(:,:), eso_tmp(:), eso_tmp2(:,:), gtens_input(:,:), Hexp(:), JAex(:,:), JAex9(:,:,:), &
                              JDMex(:,:), Jex(:), JITOexI(:,:,:,:,:), JITOexR(:,:,:,:,:), MagnCoords(:,:), Mexp(:,:), &
                              R_LG(:,:,:,:), R_ROT(:,:,:,:), riso(:,:,:), T(:), TempMagn(:), Texp(:), W(:), XLM(:,:,:,:), &
                              XRM(:,:,:,:), XT_no_field(:), XTexp(:), ZLM(:,:), ZRM(:,:)
complex(kind=wp), allocatable :: dipexch(:,:,:), dipexch_tmp(:,:,:), dipso(:,:,:,:), dipso_tmp(:,:,:), dipso_tmp2(:,:,:,:), &
                                 s_exch(:,:,:), s_exch_tmp(:,:,:), s_so(:,:,:,:), s_so_tmp(:,:,:), s_so_tmp2(:,:,:,:), Z(:,:)
#ifdef _DEBUGPRINT_
#  define _DBG_ .true.
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: dbg = _DBG_
integer(kind=iwp), external :: IsFreeUnit

#include "warnings.h"

!-----------------------------------------------------------------------
! Constants:
GRAD = .false.
!-----------------------------------------------------------------------
! Allocate memory for all arrays:
!-----------------------------------------------------------------------
#ifdef _DEBUGPRINT_
write(u6,*) '      exch = ',exch
write(u6,*) '     nPair = ',nPair
write(u6,*) '     nMult = ',nMult
write(u6,*) '      nneq = ',nneq
write(u6,*) '      neqv = ',neqv
write(u6,*) '      nLoc = ',nLoc
write(u6,*) '   nDirZee = ',nDirZee
write(u6,*) '      nDir = ',nDir
write(u6,*) '        nH = ',nH
write(u6,*) '        nT = ',nT
write(u6,*) ' nTempMagn = ',nTempMagn
write(u6,*) 'old_format = ',old_aniso_format
write(u6,*) '   MxRank1 = ',MxRank1
write(u6,*) '   MxRank2 = ',MxRank2
#endif

mem = 0

! exchange energy spectrum
call mma_allocate(W,exch,'W')
mem = mem+size(W)*RtoB
! exchange eigenvectors:
call mma_allocate(Z,exch,exch,'Z')
mem = mem+size(Z)*CtoB
! magnetic moment
call mma_allocate(dipexch,3,exch,exch,'dipexch')
mem = mem+size(dipexch)*CtoB
! spin moment
call mma_allocate(s_exch,3,exch,exch,'s_exch')
mem = mem+size(s_exch)*CtoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 1 =',mem
#endif

! index of metal site for each interacting pair:
call mma_allocate(i_pair,nPair,2,'i_pair')
mem = mem+size(i_pair)*ItoB
! exchange Lines-1 parameter
call mma_allocate(Jex,nPair,'Jex')
Jex(:) = Zero
mem = mem+size(Jex)*RtoB
! exchange Lines-3 parameter
call mma_allocate(JAex,nPair,3,'Jex')
JAex(:,:) = Zero
mem = mem+size(JAex)*RtoB
! exchange Lines-9 parameter
call mma_allocate(JAex9,nPair,3,3,'Jex')
JAex9(:,:,:) = Zero
mem = mem+size(JAex9)*RtoB
! exchange Dzyaloshinsky-Morya parameter
call mma_allocate(JDMex,nPair,3,'Jex')
JDMex(:,:) = Zero
mem = mem+size(JDMex)*RtoB
! index of ITO ranks for for each interacting pair:
call mma_allocate(imaxrank,nPair,2,'imaxrank')
imaxrank(:,:) = 0
mem = mem+size(imaxrank)*ItoB
! exchange JITO
#ifdef _DEBUGPRINT_
write(u6,'(A,I10)') 'nPair  =',nPair
write(u6,'(A,I10)') 'MxRank1=',MxRank1
write(u6,'(A,I10)') 'MxRank2=',MxRank2
write(u6,'(A,I10)') 'ibuf   =',nPair*MxRank1*MxRank2*(2*MxRank1+1)*(2*MxRank2+1)
call xFlush(u6)
#endif
call mma_allocate(JITOexR,[1,nPair],[1,MxRank1],[-MxRank1,MxRank1],[1,MxRank2],[-MxRank2,MxRank2],'JITOexR')
call mma_allocate(JITOexI,[1,nPair],[1,MxRank1],[-MxRank1,MxRank1],[1,MxRank2],[-MxRank2,MxRank2],'JITOexI')
JITOexR(:,:,:,:,:) = Zero
JITOexI(:,:,:,:,:) = Zero
mem = mem+size(JITOexR)*RtoB
mem = mem+size(JITOexI)*RtoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 2 =',mem
#endif

! index of metal site for each interacting pair:
call mma_allocate(nDim,nMult,'nDim')
nDim(:) = 0
mem = mem+size(nDim)*ItoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 3 =',mem
#endif

! number of equivalent centers, per type
call mma_allocate(neq,nneq,'neq')
mem = mem+size(neq)*ItoB
! local number of spin orbit states:
call mma_allocate(nss,nneq,'nss')
nss(:) = 0
mem = mem+size(nss)*ItoB
! local number of spin free states:
call mma_allocate(nsfs,nneq,'nsfs')
nsfs(:) = 0
mem = mem+size(nsfs)*ItoB
! local basis for exchange:
call mma_allocate(nexch,nneq,'nexch')
mem = mem+size(nexch)*ItoB
! gtens_input(:,:)
call mma_allocate(gtens_input,3,nneq,'gtens_input')
mem = mem+size(gtens_input)*RtoB
! D_fact
call mma_allocate(D_fact,nneq,'D_fact')
mem = mem+size(D_fact)*RtoB
! EoverD_fact
call mma_allocate(EoverD_fact,nneq,'EoverD_fact')
mem = mem+size(EoverD_fact)*RtoB
! MagnCoords()
call mma_allocate(MagnCoords,nneq,3,'MagnCoords')
MagnCoords(:,:) = Zero
mem = mem+size(MagnCoords)*RtoB
! riso
call mma_allocate(riso,3,3,nneq,'riso')
riso(:,:,:) = Zero
mem = mem+size(riso)*RtoB

! R_LG
call mma_allocate(r_lg,nneq,neqv,3,3,'r_lg')
r_lg(:,:,:,:) = Zero
mem = mem+size(r_lg)*RtoB
! R_ROT
call mma_allocate(r_rot,nneq,neqv,3,3,'r_rot')
r_rot(:,:,:,:) = Zero
mem = mem+size(r_rot)*RtoB

! local spin orbit states
call mma_allocate(eso,nneq,nLoc,'eso')
mem = mem+size(eso)*RtoB
call mma_allocate(eso_au,nneq,nLoc,'eso_au')
mem = mem+size(eso_au)*RtoB
! local magnetic moment
call mma_allocate(dipso,nneq,3,nLoc,nLoc,'dipso')
mem = mem+size(dipso)*CtoB
! local spin moment
call mma_allocate(s_so,nneq,3,nLoc,nLoc,'s_so')
mem = mem+size(s_so)*CtoB

! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 4 =',mem
#endif

! unit numbers of the files containing Zeeman states
call mma_allocate(LuZee,nDirZee,'LuZee')
LuZee(:) = 0
mem = mem+size(LuZee)*ItoB
! directional vectors for the magnetic field, for computing Zeeman states
call mma_allocate(dir_weight,nDirZee,3,'dir_weight')
dir_weight(:,:) = Zero
mem = mem+size(dir_weight)*ItoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 5 =',mem
#endif

! magnetization vectors
call mma_allocate(dirX,nDir,'dirX')
call mma_allocate(dirY,nDir,'dirY')
call mma_allocate(dirZ,nDir,'dirZ')
dirX(:) = Zero
dirY(:) = Zero
dirZ(:) = Zero
mem = mem+size(dirX)*RtoB
mem = mem+size(dirY)*RtoB
mem = mem+size(dirZ)*RtoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 6 =',mem
#endif

! experimental field points:
call mma_allocate(Hexp,nH,'Hexp')
Hexp(:) = Zero
mem = mem+size(Hexp)*RtoB
! experimental magnetisation:
call mma_allocate(Mexp,nH,nTempMagn,'Mexp')
Mexp(:,:) = Zero
mem = mem+size(Mexp)*RtoB
! temperature points for magnetization:
call mma_allocate(TempMagn,nTempMagn,'TempMagn')
mem = mem+size(TempMagn)*RtoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 7 =',mem
#endif

! XT for local centers, all states
call mma_allocate(XLM,nCenter,nTempMagn,3,3,'XLM')
XLM(:,:,:,:) = Zero
mem = mem+size(XLM)*RtoB
! Statistical sum for local centers, all states
call mma_allocate(ZLM,nCenter,nTempMagn,'ZLM')
ZLM(:,:) = Zero
mem = mem+size(ZLM)*RtoB
! XT for local centers, exchange states
call mma_allocate(XRM,nCenter,nTempMagn,3,3,'XRM')
XRM(:,:,:,:) = Zero
mem = mem+size(XRM)*RtoB
! Statistical sum for local centers, exchange states
call mma_allocate(ZRM,nCenter,nTempMagn,'ZRM')
ZRM(:,:) = Zero
mem = mem+size(ZRM)*RtoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 8 =',mem
#endif

! T experimental given by user in the input
call mma_allocate(Texp,nT,'Texp')
Texp(:) = Zero
mem = mem+size(Texp)*RtoB
! XT experimental given by user in the input
call mma_allocate(chit_exp,nT,'chit_exp')
chit_exp(:) = Zero
mem = mem+size(chit_exp)*RtoB

! -- add nTempMagn points, so that all measurables are computed at once...
! temperature points for which XT will be computed
call mma_allocate(T,nTempMagn+nT,'Temperature')
mem = mem+size(T)*RtoB
! XT experimental tmagn XTexp+Tmagn
call mma_allocate(XTexp,nTempMagn+nT,'XTexp')
mem = mem+size(XTexp)*RtoB
! XT in the absence of the magnetic field
call mma_allocate(XT_no_field,nTempMagn+nT,'XT_no_field')
mem = mem+size(XT_no_field)*RtoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 9 =',mem
#endif

write(u6,'(A,I16,A)') 'The code allocated at least:',mem,' bytes of memory for this run.'
!-----------------------------------------------------------------------
! set default values of the main variables and arrays:
#ifdef _DEBUGPRINT_
call xFlush(u6)
write(u6,*) 'Enter set_defaults'
#endif
call set_defaults(nneq,nTempMagn,nDir,nDirZee,nMult,neq,nexch,nK,mG,ncut,nP,AngPoints,nBlock,encut_definition,iopt,iPrint,nsymm, &
                  ngrid,dltT0,dltH0,zJ,tmin,tmax,hmin,hmax,XField,thrs,TempMagn,cryst,coord,encut_rate,gtens_input,D_fact, &
                  EoverD_fact,riso,decompose_exchange,AnisoLines1,AnisoLines3,AnisoLines9,DM_exchange,Dipol,KE,JITO_exchange, &
                  fitCHI,fitM,Do_structure_abc,doplot,compute_g_tensors,tinput,compute_susceptibility,compute_torque, &
                  compute_barrier,compute_magnetization,hinput,compute_Mdir_vector,zeeman_energy,m_paranoid,m_accurate,smagn,itype)
#ifdef _DEBUGPRINT_
write(u6,*) 'Exit set_defaults'
#endif
if (neqv > 1) then
# ifdef _DEBUGPRINT_
  write(u6,*) 'Enter fetch_neq'
# endif
  call fetch_neq(nneq,neq(1:nneq),nexch(1:nneq))
# ifdef _DEBUGPRINT_
  write(u6,*) 'Exit fetch_neq'
# endif
end if

!-----------------------------------------------------------------------
! read the Standard Input:
#ifdef _DEBUGPRINT_
write(u6,*) 'Enter Readin_poly'
#endif
call Readin_poly(nneq,neq,neqv,exch,nCenter,nT,nH,nTempMagn,nDir,nDirZee,nMult,nPair,nexch,nDim,i_pair,lant,multLn,iPrint,keopt, &
                 encut_definition,nK,mG,iopt,nP,AngPoints,ncut,LUZee,MxRank1,MxRank2,imaxrank,nsymm,ngrid,TempMagn,R_LG,R_ROT,Jex, &
                 JAex,JAex9,JDMex,JITOexR,JITOexI,tpar,upar,cryst,coord,Xfield,gtens_input,D_fact,EoverD_fact,riso,MagnCoords, &
                 thrs,tmin,tmax,hmin,hmax,Texp,chit_exp,Hexp,Mexp,encut_rate,zJ,dirX,dirY,dirZ,dir_weight,Title,itype,ifHDF, &
                 compute_g_tensors,compute_magnetization,TINPUT,HINPUT,Do_structure_abc,doplot,compute_Mdir_vector,zeeman_energy, &
                 m_paranoid,m_accurate,smagn,compute_susceptibility,decompose_exchange,KE,fitCHI,fitM,compute_torque, &
                 compute_barrier,Dipol,check_title,AnisoLines1,AnisoLines3,AnisoLines9,DM_exchange,JITO_exchange)
#ifdef _DEBUGPRINT_
write(u6,*) 'Exit Readin_poly'
call xFlush(u6)
#endif
!-----------------------------------------------------------------------
! fetch the data from aniso_x.input files: (formatted ANISOINPUT)
do i=1,nneq
# ifdef _DEBUGPRINT_
  write(u6,'(A,A)') 'itype(i)=',itype(i)
  !write(u6,*) '   nss(i)=',nss(i)
  !write(u6,*) '  nsfs(i)=',nsfs(i)
  write(u6,*) ' nexch(i)=',nexch(i)
  write(u6,*) '     nLoc=',nLoc
  write(u6,*) ' gtens_input(1:3,i)=',(gtens_input(j,i),j=1,3)
  write(u6,*) 'D_fact(i)=',D_fact(i)
  !write(u6,*) 'eso(i,1:nexch(i))=',(eso(i,j),j=1,nexch(i))
  write(u6,*) 'old_aniso_format=',old_aniso_format
  call xFlush(u6)
# endif

  if ((itype(i) == 'B') .or. (itype(i) == 'C')) then

#   ifdef _DEBUGPRINT_
    write(u6,*) 'Enter generate_isotrop_site'
    call xFlush(u6)
#   endif
    call mma_allocate(eso_tmp,nexch(i),label='eso_tmp')
    call mma_allocate(dipso_tmp,3,nexch(i),nexch(i),label='dipso_tmp')
    call mma_allocate(s_so_tmp,3,nexch(i),nexch(i),label='s_so_tmp')
    call generate_isotrop_site(nss(i),nsfs(i),nexch(i),gtens_input(:,i),riso(:,:,i),D_fact(i),EoverD_fact(i),eso_tmp,dipso_tmp, &
                               s_so_tmp)
    eso(i,1:nexch(i)) = eso_tmp(:)
    dipso(i,:,1:nexch(i),1:nexch(i)) = dipso_tmp(:,:,:)
    s_so(i,:,1:nexch(i),1:nexch(i)) = s_so_tmp(:,:,:)
    call mma_deallocate(eso_tmp)
    call mma_deallocate(dipso_tmp)
    call mma_deallocate(s_so_tmp)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Exit generate_isotrop_site'
    call xFlush(u6)
#   endif

  else if (itype(i) == 'A') then

    if (ifHDF) then
      ! set their names:
      if (i < 10) then
        write(namefile_aniso(i),'(4A)') 'aniso_hdf_',char(48+mod(int(i),10)),'.input'
#       ifdef _DEBUGPRINT_
        write(u6,'(A,i2,A,A)') 'namefile_aniso(',i,')=',namefile_aniso(i)
#       endif
      else if ((i >= 10) .and. (i <= 99)) then
        write(namefile_aniso(i),'(4A)') 'aniso_hdf_',char(48+mod(i/10,10)),char(48+mod(i,10)),'.input'
      end if
#     ifdef _DEBUGPRINT_
      write(u6,*) 'PA:  namefile_aniso(i)=',trim(namefile_aniso(i))
#     endif

#     ifdef _HDF5_
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Enter read_hdf5_poly'
#     endif
      call read_hdf5_init(namefile_aniso(i),nsfs(i),nss(i))
#     ifdef _DEBUGPRINT_
      write(u6,*) 'i=',i,'nsfs(i),nss(i)=',nsfs(i),nss(i)
#     endif
      call mma_allocate(eso_tmp,nss(i),label='eso_tmp')
      call mma_allocate(dipso_tmp,3,nss(i),nss(i),label='dipso_tmp')
      call mma_allocate(s_so_tmp,3,nss(i),nss(i),label='s_so_tmp')
      call read_hdf5_poly(namefile_aniso(i),nss(i),nsfs(i),eso_tmp,dipso_tmp,s_so_tmp,iReturn)
      eso(i,1:nss(i)) = eso_tmp(:)
      dipso(i,:,1:nss(i),1:nss(i)) = dipso_tmp(:,:,:)
      s_so(i,:,1:nss(i),1:nss(i)) = s_so_tmp(:,:,:)
      call mma_deallocate(eso_tmp)
      call mma_deallocate(dipso_tmp)
      call mma_deallocate(s_so_tmp)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Exit read_hdf5_poly'
      write(u6,*) 'ESO(i)=',(ESO(i,j),j=1,nss(i))
#     endif
#     else
      call WarningMessage(2,'File '//trim(namefile_aniso(i))//' cannot be opened. Molcas was compiled without HDF5 option.')
      call Quit_OnUserError()
#     endif
    else
      ! set their names:
      if (i < 10) then
        write(namefile_aniso(i),'(4A)') 'aniso_',char(48+mod(i,10)),'.input'
      else if ((i >= 10) .and. (i <= 99)) then
        write(namefile_aniso(i),'(4A)') 'aniso_',char(48+mod(i/10,10)),char(48+mod(i,10)),'.input'
      end if

      if (old_aniso_format) then
        ! get the information from OLD formatted
        ! aniso.input file:
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Enter read_formatted_aniso_poly'
        !write(u6,*) 'nss(i) =',nss(i)
        !write(u6,*) 'nsfs(i)=',nsfs(i)
        !write(u6,*) 'nLoc   =',nLoc
        call xFlush(u6)
#       endif

        call mma_allocate(eso_tmp,nLoc,label='eso_tmp')
        call mma_allocate(dipso_tmp,3,nLoc,nLoc,label='dipso_tmp')
        call mma_allocate(s_so_tmp,3,nLoc,nLoc,label='s_so_tmp')
        call read_formatted_aniso_poly(namefile_aniso(i),nss(i),nsfs(i),nLoc,eso_tmp,dipso_tmp,s_so_tmp,iReturn)
        eso(i,1:nLoc) = eso_tmp(:)
        dipso(i,:,1:nLoc,1:nLoc) = dipso_tmp(:,:,:)
        s_so(i,:,1:nLoc,1:nLoc) = s_so_tmp(:,:,:)
        call mma_deallocate(eso_tmp)
        call mma_deallocate(dipso_tmp)
        call mma_deallocate(s_so_tmp)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Exit read_formatted_aniso_poly'
#       endif

      else ! old_aniso_format is false, i.e. the default

        ! get the information from NEW formatted
        ! aniso.input file:
#       ifdef _DEBUGPRINT_
        write(u6,*) 'read formatted_aniso_poly_NEW'
#       endif
        LuAniso = IsFreeUnit(81)
        call molcas_open(LuAniso,namefile_aniso(i))
        ! nss(i) is yet undefined up till this place
        call read_nss(LuAniso,nss(i),dbg)
        call read_nstate(LuAniso,nsfs(i),dbg)
        call mma_allocate(eso_tmp,nss(i),label='eso_tmp')
        call mma_allocate(dipso_tmp,3,nss(i),nss(i),label='dipso_tmp')
        call mma_allocate(s_so_tmp,3,nss(i),nss(i),label='s_so_tmp')
        call read_eso(LuAniso,nss(i),eso_tmp,dbg)
        eso_au(i,1:nss(i)) = eso_tmp(:)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'poly_aniso: eso_au=',(eso_au(i,j),j=1,nss(i))
#       endif
        call read_magnetic_moment(LuAniso,nss(i),dipso_tmp,dbg)
        dipso(i,:,1:nss(i),1:nss(i)) = dipso_tmp(:,:,:)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Call read_spin_moment'
        call xFlush(u6)
#       endif
        call read_spin_moment(LuAniso,nss(i),s_so_tmp,dbg)
        s_so(i,:,1:nss(i),1:nss(i)) = s_so_tmp(:,:,:)
        ! compute the relative spin-orbit energies in cm-1
        eso(i,1:nss(i)) = (eso_au(i,1:nss(i))-eso_au(i,1))*auTocm
        close(LuAniso)
        call mma_deallocate(eso_tmp)
        call mma_deallocate(dipso_tmp)
        call mma_deallocate(s_so_tmp)

      end if
    end if ! ifHDF
  else
    call quit(_RC_INTERNAL_ERROR_)
  end if
end do

#ifdef _DEBUGPRINT_
do i=1,nneq
  call mma_allocate(dipso_tmp,3,nss(i),nss(i),label='dipso_tmp')
  call mma_allocate(s_so_tmp,3,nss(i),nss(i),label='s_so_tmp')
  dipso_tmp(:,:,:) = dipso(i,:,1:nss(i),1:nss(i))
  s_so_tmp(:,:,:) = s_so(i,:,1:nss(i),1:nss(i))
  write(u6,*) 'MAGNETIC MOMENT  ON  SITE 1'
  call prMom('site i',dipso_tmp,nss(i))
  write(u6,*) 'SPIN MOMENT  ON  SITE 1'
  call prMom('site i',s_so_tmp,nss(i))
  call mma_deallocate(dipso_tmp)
  call mma_deallocate(s_so_tmp)
end do
#endif

#ifdef _DEBUGPRINT_
write(u6,*) 'Enter input_process'
#endif
call input_process(nneq,neq,neqv,nmax,nCenter,nexch,nDir,nDirZee,nDirTot,nH,nT,exch,nBlock,nTempMagn,iopt,nMult,nDim,nPair,i_pair, &
                   nP,AngPoints,lant,multLn,KEOPT,encut_definition,nK,mG,ncut,nsfs,nss,nLoc,MxRank1,MxRank2,imaxrank,iPrint,nsymm, &
                   ngrid,R_LG,gtens_input,dirX,dirY,dirZ,dir_weight,zJ,cryst,coord,hmin,hmax,TempMagn,thrs,Hexp,Mexp,tmin,tmax, &
                   chit_exp,Texp,Xfield,Jex,JAex,JAex9,JDMex,tpar,upar,MagnCoords,encut_rate,eso,JITOexR,JITOexI,Title,itype, &
                   namefile_aniso,Do_structure_abc,old_aniso_format,compute_barrier,fitCHI,fitM,hinput,tinput, &
                   compute_magnetization,compute_Mdir_vector,zeeman_energy,m_paranoid,m_accurate,smagn,compute_g_tensors, &
                   compute_torque,compute_susceptibility,AnisoLines1,AnisoLines3,AnisoLines9,Dipol,check_title,KE,DM_exchange, &
                   JITO_exchange)
#ifdef _DEBUGPRINT_
write(u6,*) 'Exit input_process'
#endif

! at this moment all input values are set. proceed to compute various
! properties:
!-----------------------------------------------------------------------
! in case J coupling parameters are not known, proceed to get them from
! a minimization scheme between caclulated and measured magnetic susceptibility
!
! the function below optimizes the N+1 parameters:
! N - exchange couplings J, and one parameter for total shift of all experimental points;
!if (fitCHI) Call fitCHI()
!c at this point all J parameters are known; we can proceed to compute
! the energy spectrum and the resulting properties
!-----------------------------------------------------------------------
!Lines = .true.
!if (AnisoLines) Lines = .false.
!Dipol = dipole_included
!KE = exch_long

#ifdef _DEBUGPRINT_
write(u6,*) 'Dipol         = ',Dipol
write(u6,*) 'AnisoLines1   = ',AnisoLines1
write(u6,*) 'AnisoLines3   = ',AnisoLines3
write(u6,*) 'AnisoLines9   = ',AnisoLines9
write(u6,*) 'DM_exchange   = ',DM_exchange
write(u6,*) 'JITO_exchange = ',JITO_exchange
write(u6,*) 'nmax          = ',nmax
#endif

call mma_allocate(eso_tmp2,nneq,nmax,label='eso_tmp2')
call mma_allocate(dipso_tmp2,nneq,3,nmax,nmax,label='dipso_tmp2')
call mma_allocate(s_so_tmp2,nneq,3,nmax,nmax,label='s_so_tmp2')
eso_tmp2(:,:) = eso(1:nneq,1:nmax)
dipso_tmp2(:,:,:,:) = dipso(1:nneq,:,1:nmax,1:nmax)
s_so_tmp2(:,:,:,:) = s_so(1:nneq,:,1:nmax,1:nmax)

call exchctl(exch,nneq,neqv,neq,nexch,nmax,nCenter,npair,i_pair,MxRank1,MxRank2,imaxrank,Jex,JAex,JAex9,JDMex,JITOexR,JITOexI, &
             eso_tmp2,s_so_tmp2,dipso_tmp2,MagnCoords,R_ROT,R_LG,riso,tpar,upar,lant,itype,Dipol,AnisoLines1,AnisoLines3, &
             AnisoLines9,KE,KEOPT,DM_exchange,JITO_exchange,W,Z,S_EXCH,DIPEXCH,iPrint,mem)
dipso(1:nneq,:,1:nmax,1:nmax) = dipso_tmp2(:,:,:,:)
s_so(1:nneq,:,1:nmax,1:nmax) = s_so_tmp2(:,:,:,:)
!-----------------------------------------------------------------------
!     generate an 'aniso -like file' for testign exchange proj in SA
fname = 'ANISOINPUT_POLY'
call write_formatted_aniso_poly(fname,exch,W,dipexch,s_exch)
!-----------------------------------------------------------------------
! compute the magnetic moments on individual metal sites, for
! the lowest NSTA exchange states;
!   NMAX = maximal number of local states which participate into the exchange coupling
call MOMLOC2(exch,nmax,nneq,neq,neqv,r_rot,nCenter,nExch,W,Z,dipexch,s_exch,dipso_tmp2,s_so_tmp2)

call mma_deallocate(eso_tmp2)
call mma_deallocate(dipso_tmp2)
call mma_deallocate(s_so_tmp2)

!-----------------------------------------------------------------------
if (compute_g_tensors) then
# ifdef _DEBUGPRINT_
  write(u6,'(A,90I3)') 'nmult= ',nmult
  write(u6,'(A,90I3)') 'ndim() = ',ndim(1:nmult)
# endif

  Ifunct = 0
  do imltpl=1,nmult
    iReturn = 0
    i1 = 1+Ifunct
    i2 = ndim(imltpl)+Ifunct

    if (ndim(imltpl) > 1) then
      if (i2 > exch) exit

      call mma_allocate(s_exch_tmp,3,ndim(imltpl),ndim(imltpl),label='s_exch_tmp')
      call mma_allocate(dipexch_tmp,3,ndim(imltpl),ndim(imltpl),label='dipexch_tmp')
      s_exch_tmp(:,:,:) = s_exch(:,i1:i2,i1:i2)
      dipexch_tmp(:,:,:) = dipexch(:,i1:i2,i1:i2)

#     ifdef _DEBUGPRINT_
      write(u6,'(A,90I3)') 'ndim(imltpl)=',ndim(imltpl)
      call prmom('PA: s_exch:',s_exch_tmp,ndim(imltpl))
      call prmom('PA: dip_exch:',dipexch_tmp,ndim(imltpl))
#     endif

      call g_high(w(i1:i2),GRAD,s_exch_tmp,dipexch_tmp,imltpl,ndim(imltpl),Do_structure_abc,cryst,coord,gtens,maxes,iprint)

      call mma_deallocate(s_exch_tmp)
      call mma_deallocate(dipexch_tmp)

    end if
    Ifunct = Ifunct+ndim(imltpl)
    if (iReturn /= 0) call Abend()
  end do ! imltpl
end if

!----------------------------------------------------------------------|
! >> AB INITIO BLOCKING BARRIER FOR SMMs <<
if (compute_barrier) then
  if (nBlock /= 0) then
    imanifold = 1

    write(u6,'(A)') 'UBAR:: matrix elements of the (input) magnetic moment '
    write(u6,'(A)') '                |         - projection - X -        |         - projection - Y -        |'// &
                    '         - projection - Z -        |'
    write(u6,'(A)') '----------------|----- Real ------|---- Imaginary --|----- Real ------|---- Imaginary --|-'// &
                    '---- Real ------|---- Imaginary --|'
    do i=1,nBlock
      do j=1,i
        write(u6,'(A,i2,A,i2,A,6ES18.10)') '< ',i,'|M_xyz|',j,' > ',(dipexch(l,i,j),l=1,3)
      end do
    end do

    call mma_allocate(dipexch_tmp,3,nBlock,nBlock,label='dipexch_tmp')
    dipexch_tmp(:,:,:) = dipexch(:,1:nBlock,1:nBlock)
    call BARRIER(nBlock,dipexch_tmp,W(1:nBlock),imanifold,nMult,nDim(1:nMult),DoPlot,iprint)
    call mma_deallocate(dipexch_tmp)

  else
    write(u6,'(A)') 'nBlock parameter is not defined. '
    write(u6,'(A)') 'Did you specify the MLTP keyword in the input?'
    write(u6,'(A)') 'If the problem persists,please, submit a bug report.'
  end if
end if

!-----------------------------------------------------------------------
! MAGNETIC SUSCEPTIBILITY
!-----------------------------------------------------------------------
if (compute_susceptibility .and. (nT > 0)) then

  ! set nT, T(i) and XTexp(i) arrays:
  call set_T(nT,nTempMagn,TINPUT,TempMagn,Tmin,Tmax,chit_exp,Texp,T,XTexp)

  ! XT at H=0
  call susceptibility_pa(exch,nLoc,nCenter,nneq,neqv,neq,nss,nexch,nTempMagn,nT,Tmin,Tmax,XTexp,eso,dipso,s_so,W,dipexch,s_exch,T, &
                         R_LG,zJ,tinput,XLM,ZLM,XRM,ZRM,iopt,XT_no_field,DoPlot,mem)

  if (Xfield /= Zero) then
    ! XT at H>0:

    call set_nm(exch,ncut,encut_definition,nk,mg,nTempMagn,hmax,w,encut_rate,TempMagn,nM,EM,dbg)
    do i=1,nT+nTempMagn
      write(u6,*) i,T(i)
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) 'nm    =',nm
    write(u6,*) 'Tmin  =',Tmin
    write(u6,*) 'Tmax  =',Tmax
    write(u6,*) 'THRS  =',THRS
    write(u6,*) 'smagn =',smagn
    write(u6,*) 'm_par =',m_paranoid
    write(u6,*) 'm_acc =',m_accurate
    write(u6,*) 'tinpu =',tinput
#   endif

    ! nm = exch
    call XT_dMoverdH(exch,nLoc,nCenter,nneq,neqv,neq,nss,nexch,nTempMagn,nT,exch,iopt,mem,Tmin,Tmax,XTexp,eso,w,T,R_ROT,zJ,Xfield, &
                     THRS,XT_no_field,dipso,s_so,dipexch,s_exch,tinput,smagn,m_paranoid,m_accurate,doplot)
  end if
else
  write(u6,'(A)') 'Computation of the magnetic susceptibility... skipped by the user'
end if
!-----------------------------------------------------------------------
! MAGNETIC TORQUE
!-----------------------------------------------------------------------

if (compute_torque) then
!        set the NM- number of states to be exactly diagonalized in Zeeman Interaction
  call set_nm(exch,ncut,encut_definition,nk,mg,nTempMagn,hmax,w,encut_rate,TempMagn,nM,EM,dbg)

  ! nm = exch
  call torque_pa(nneq,nCenter,neq,neqv,nLoc,exch,nTempMagn,nH,exch,AngPoints,nexch,iopt,nss,mem,smagn,m_paranoid,m_accurate, &
                 TempMagn,w,hmin,hmax,dltH0,EM,zJ,THRS,hexp,dipexch,s_exch,dipso,s_so,eso,hinput,r_rot,XLM,ZLM,XRM,ZRM)

else
  write(u6,'(A)') 'Computation of the magnetization torque ... skipped by the user'
end if !compute_torque
!-----------------------------------------------------------------------
! MOLAR MAGNETIZATION
!-----------------------------------------------------------------------
if (compute_magnetization .and. (nH > 0)) then
  ! set the NM- number of states to be exactly diagonalized in Zeeman Interaction
  call set_nm(exch,ncut,encut_definition,nk,mg,nTempMagn,hmax,w,encut_rate,TempMagn,nM,EM,dbg)

  call magnetization_pa(exch,nLoc,nM,nH,nneq,neq,neqv,nCenter,nTempMagn,nDir,nDirZee,nDirTot,nss,nexch,iopt,LUZee,nsymm,ngrid, &
                        TempMagn,hexp,mexp,hmin,hmax,em,zJ,thrs,dirX,dirY,dirZ,dir_weight,w,dipexch,s_exch,dipso,s_so,eso,hinput, &
                        r_rot,XLM,ZLM,XRM,ZRM,zeeman_energy,compute_Mdir_vector,m_paranoid,m_accurate,smagn,mem,doplot)
else
  write(u6,'(A)') 'Computation of the molar magnetization ... skipped by the user'
end if !compute_magnetization

!-----------------------------------------------------------------------
! Deallocate memory for big arrays:
!-----------------------------------------------------------------------
call mma_deallocate(W)
call mma_deallocate(Z)
call mma_deallocate(dipexch)
call mma_deallocate(s_exch)
call mma_deallocate(i_pair)
call mma_deallocate(Jex)
call mma_deallocate(JAex)
call mma_deallocate(JAex9)
call mma_deallocate(JDMex)
call mma_deallocate(imaxrank)
call mma_deallocate(JITOexR)
call mma_deallocate(JITOexI)

call mma_deallocate(nDim)
call mma_deallocate(neq)
call mma_deallocate(nss)
call mma_deallocate(nsfs)
call mma_deallocate(nexch)
call mma_deallocate(gtens_input)
call mma_deallocate(D_fact)
call mma_deallocate(EoverD_fact)
call mma_deallocate(MagnCoords)
call mma_deallocate(riso)
call mma_deallocate(r_lg)
call mma_deallocate(r_rot)
call mma_deallocate(eso)
call mma_deallocate(eso_au)
call mma_deallocate(dipso)
call mma_deallocate(s_so)

call mma_deallocate(LuZee)
call mma_deallocate(dir_weight)

call mma_deallocate(dirX)
call mma_deallocate(dirY)
call mma_deallocate(dirZ)

! experimental field points:
call mma_deallocate(Hexp)
call mma_deallocate(Mexp)
call mma_deallocate(TempMagn)

call mma_deallocate(XLM)
call mma_deallocate(ZLM)
call mma_deallocate(XRM)
call mma_deallocate(ZRM)
call mma_deallocate(Texp)
call mma_deallocate(chit_exp)
call mma_deallocate(T)
call mma_deallocate(XTexp)
call mma_deallocate(XT_no_field)
!-----------------------------------------------------------------------
write(u6,*)
write(u6,'(10A)') ('-@-#-$-%-&-*-',i=1,10)
write(u6,*)
write(u6,'(10X,A)') 'HAPPY   LANDING !!!   POLY_ANISO ENDED  OK !'
write(u6,*)
write(u6,'(10A)') ('-*-&-%-$-#-@-',i=1,10)
call xFlush(u6)

return

end subroutine POLY_ANISO_1
