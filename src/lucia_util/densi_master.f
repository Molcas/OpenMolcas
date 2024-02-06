************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE densi_master()
      use stdalloc, only: mma_allocate, mma_deallocate
      use GLBBAS
      use rasscf_lucia
      use Lucia_Interface, only: rVec
*
* Controls the calculation of the densities, when Lucia is called
* from Molcas Rasscf.
*
      implicit real*8 (a-h,o-z)
#include "mxpdim.fh"
#include "crun.fh"
#include "cicisp.fh"
#include "clunit.fh"
#include "orbinp.fh"
#include "lucinp.fh"
#include "spinfo_lucia.fh"
#include "cstate.fh"
#include "io_util.fh"
      logical iPack,tdm
      dimension dummy(1)
      Real*8, Allocatable:: VEC1(:), VEC2(:)
      Integer, Allocatable:: lVec(:)
      Real*8, Allocatable, Target:: SCR1(:), SCR3(:)
      Real*8, Allocatable:: SCR2(:), SCR4(:)
      Integer NSD, NCSF

*
* Put CI-vector from RASSCF on luc
*
*     if rvec is associated, it should be a pointer to a second CI vector
*     and a one-particle transition density matrix will be computed
      tdm = Associated(rVec)

      NSD = NSD_PER_SYM(IREFSM)
      NCSF= NCSF_PER_SYM(IREFSM)
      CALL mma_allocate(SCR1,NSD,Label='SCR1')
      CALL mma_allocate(SCR2,NSD,Label='SCR2')

      CALL COPVEC(C_POINTER,SCR1,NCSF)

      Call mma_allocate(lVec,MXNTTS,Label='lVec')
      IF (tdm) THEN
         CALL mma_allocate(SCR3,NSD,Label='SCR3')
         CALL mma_allocate(SCR4,NSD,Label='SCR4')

         CALL COPVEC(rvec,SCR3,NCSF)
         CALL CSDTVC(SCR3,SCR4,1,DTOC,SDREO, IREFSM, 1)
         CALL CPCIVC(SCR3,NSD,LUHC, MXNTTS, IREFSM, 1,lVec)
         CALL mma_deallocate(SCR3)
         CALL mma_deallocate(SCR4)
      END IF
      CALL CSDTVC(SCR1,SCR2,1,DTOC,SDREO, IREFSM, 1)
      CALL CPCIVC(SCR1,NSD,LUC, MXNTTS, IREFSM, 1,lVec)
      Call mma_deallocate(lVec)

*
* Determine length of arrays VEC1 and VEC2
*
c      IF(ISIMSYM.EQ.0) THEN
         LBLOCK = MXSOOB
c      ELSE
c         LBLOCK = MXSOOB_AS
c      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
* JESPER : Should reduce I/O
*PAM06      LBLOCK = MAX(XISPSM(IREFSM,1),DBLE(MXSOOB))
      LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
      IF(PSSIGN.NE.0.0D0) LBLOCK = 2*INT(XISPSM(IREFSM,1))
*
* Allocate arrays
*
      Call mma_allocate(VEC1,LBLOCK,Label='VEC1')
      Call mma_allocate(VEC3,kvec3_length,Label='VEC3')
*
* Copy Sigma-vector from disc to core
*
      Call mma_allocate(VEC2,LBLOCK,Label='VEC2')
       IF (iSigma_on_disk .ne. 0) THEN
          Call mma_allocate(lVec,MXNTTS,Label='lVec')
          CALL cpsivc(lusc34, mxntts, vec2,lVec)
          Call mma_deallocate(lVec)
       ELSE
          vec2(:) = 0.0d0
       ENDIF
*
* Information needed on file handling
*
      LBLK = - 1
*
* Copy vector on file LUC to LUSC1 and LUHC
*
      IDISK(LUC)=0
      IDISK(LUSC1)=0
      CALL COPVCD(LUC,LUSC1,VEC1,0,LBLK)
      IF (.not.tdm) CALL COPVCD(LUSC1,LUHC,VEC1,1,LBLK)
*
* Calculate one- and two-body densities
*
      IPACK = .TRUE.
      DUMMY = 0.0D0
      IF (tdm) THEN
         CALL densi2_lucia(1,Dtmp,dummy,dummy,dummy,
     &   vec1,vec2,lusc1,luhc,exps2,1,DStmp,IPACK)
      ELSE
         CALL densi2_lucia(2,rho1,dummy,Ptmp,PAtmp,
     &   vec1,vec2,lusc1,luhc,exps2,1,srho1,IPACK)
      END IF

*
* Explanation of calling parameters
*
C      2      : DONE!!! - Calculate both one and two body densities.
C      rho1  : DONE!!! - Output - include in module glbbas
C      krho2  : DONE!!! - Output - include in moduke glbbas
C      vec1  : DONE!!! - CI-vector
C      vec2  : DONE!!! - Sigma-vector
C      lusc1  : DONE!!! - file pointer
C      luhc   : DONE!!! - file pointer
C      exps2  : DONE!!! - Output - expectation value of S**2.
C      1      : DONE!!! - Calculate spin density
C      srho1 : DONE!!! - Comming with module glbbas
*
      IF (.not.tdm) THEN
*        Save densities in trigonal format for use in Molcas
*
         CALL TriPak(rho1, Dtmp, 1, ntoob, ntoob)
         CALL TriPak(srho1, DStmp, 1, ntoob, ntoob)
      END IF
*
      CALL CSDTVC(scr1,scr2,2,dtoc,SDREO, iRefSm, 1)
*
      CALL mma_deallocate(SCR1)
      CALL mma_deallocate(SCR2)
      Call mma_deallocate(VEC1)
      Call mma_deallocate(VEC2)
      Call mma_deallocate(VEC3)
*
      END
*
