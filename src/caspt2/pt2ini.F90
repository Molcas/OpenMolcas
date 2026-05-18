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

subroutine PT2INI()

use INPUTDATA, only: INPUT, READIN_CASPT2
use PT2WFN, only: PT2WFN_INIT, PT2WFN_DATA
use REFWFN, only: REFWFN_INIT, REFWFN_INFO, REFWFN_DATA, REFWFN_CLOSE
#ifdef _DMRG_
use qcmaquis_interface_cfg, only: qcmaquis_param
use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG
use caspt2_module, only: DMRG, nAshT
#endif
use caspt2_global, only: do_grad, iStpGrd
use caspt2_global, only: FIMO, FIFA, DREF, PREF, DMIX, DWGT, CMOPT2, TAT, NTAT, TORB, NTORB, NDREF, NPREF, NCMO
use stdalloc, only: mma_allocate
use ChoCASPT2, only: InfVec_N2_PT2, MaxVec_PT2, NASPlit, NISplit, NumCho_PT2
use spool, only: SpoolInp, Close_LuSpool
use Molcas, only: LenIn
use caspt2_module, only: nSym, Header, ifChol, jState, bName, nAsh, nBas, nIsh, nOTri, nBasT, nBSqT, nSsh, nState
use constants, only: Zero
use definitions, only: iwp, u6

implicit none
integer(kind=iwp) LuSpool
! Cholesky
integer(kind=iwp) iSym, iRC

call SETTIM()

! Probe the environment to globally set the IPRGLB value
call Set_Print_Level()

! Probe the RunFile for some basic information

call Get_cArray('Seward Title',Header,144)
call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
nbast = sum(nbas(1:nsym))
nbsqt = sum(nbas(1:nsym)**2)
call Get_cArray('Unique Basis Names',bName,(LenIn+8)*nbast)
jstate = 1
call DecideOnCholesky(IfChol)
! PAM Feb 2008: The following statement was moved here from
! prpctl, in case prpctl call is bypassed by keyword ''NOPROP''.
call Put_cArray('Relax Method','CASPT2  ',8)

! Allocate Input struct, read and process the Input

allocate(Input)
LuSpool = 21
call SpoolInp(LuSpool)
call READIN_CASPT2(LuSpool,nSym)
call Close_LuSpool(LuSpool)
! Initialize scratch files.
call OPNFLS_CASPT2()
! Initialize the reference wavefunction file and read basic info
call refwfn_init(Input%FILE)
call refwfn_info()
! the input processing needs some data from the reference wavefunction
! to be set, hence this phase occurs after reading that data.
call ProcInp_Caspt2()

! Allocate some arrays that will stay globally allocated:

! Finally read the MO and CI data from the refwfn file, and close it as
! we have no more need for it and the same filename might be reused for
! the pt2wfn file. We do this after input processing because we need to
! know which roots to pick up. The MOs are stored on LUONEM, at address
! IAD1M(1), and the CI arrays on LUCIEX at IDCIEX.
call refwfn_data()
call refwfn_close()

! If necessary, make modifications to the inactive/virtual orbitals. The
! new MOs, if any, will overwrite the originals on LUONEM at IAD1M(1).
if (Input%modify_correlating_MOs) call correlating_orbitals()

! After possible reconfiguration of inactives/virtuals, the computation
! of total sizes of orbital spaces is done here, before any other code
! (that might rely on these to be correctly set!)
call wfnsizes()

#ifdef _DMRG_
if (DMRG) then
  ! set the lattice length (i.e. the active space size)
  qcmaquis_param%L = nasht
  if (iPrGlb >= DEBUG) write(u6,*) 'PT2INI> qcmaquis_param%L = ',qcmaquis_param%L
end if
#endif
! Create the PT2 wavefunction file (formerly JOBMIX). The reference file
! should not be active, as it might be the same file (in which case it
! is overwritten).
call pt2wfn_init()
call pt2wfn_data()

! Global allocations active throughout the program (NOTRI was computed
! during the call to wfnsizes, so this is done here.)

! The total fock matrix (sum of inactive and active contrib.)
call mma_allocate(FIFA,NOTRI,Label='FIFA')
! The fock matrix with contributions from inactive orbitals, only.
call mma_allocate(FIMO,NOTRI,Label='FIMO')
! Density matrices, active indices.
call mma_allocate(DREF,NDREF,Label='DREF')
call mma_allocate(PREF,NPREF,Label='PREF')
! All the density matrices kept in memory
call mma_allocate(DMIX,NDREF,NSTATE,Label='DMIX')
call mma_allocate(DWGT,NSTATE,NSTATE,Label='DWGT')

! Print input data
call PRINP_CASPT2()

! INITIALIZE SPLIT-GRAPH UGA TABLES.
call SG_SETUP_CASPT2()

! Initialize superindex tables. Check memory needs.
call SIZES()

! Initialize sizes, offsets etc used in equation solver.
call EQCTL1()

if (IfChol) then
  ! initialize Cholesky information
  call Cho_X_init(irc,Zero)
  if (irc /= 0) then
    write(u6,*) 'CASPT2: Non-zero rc in Cho_X_init'
    call QUIT(irc)
  end if
  ! import_ch transfers some values from Cholesky module to chocaspt2
  call import_cho(numcho_pt2,infvec_n2_pt2,maxvec_pt2)
  call setup_cho(nSym,nIsh,nAsh,nSsh,NumCho_pt2,'Allo')
  ! get unit numbers for Cholesky MO vectors
  do iSym=1,nSym
    call Cho_Caspt2_OpenF(0,1,iSym,nIsplit(iSym))
    call Cho_Caspt2_OpenF(0,2,iSym,nAsplit(iSym))
  end do
  ! set up size info for cholesky vector transformation in tracho2
  call trachosz()
end if

! Allocate global orbital arrays:
call mma_allocate(CMOPT2,NCMO,Label='CMOPT2')
CMOPT2(:) = Zero
! Allocate global orbital transformation arrays:
call mma_allocate(TORB,NTORB,Label='TORB')
call mma_allocate(TAT,NTAT,Label='TAT')

! initialize quantities for gradient calculation
iStpGrd = 1 !! Just in case
if (do_grad) call GrdIni()

end subroutine PT2INI
