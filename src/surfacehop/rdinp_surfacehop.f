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
      SUBROUTINE rdinp_surfacehop()
      use Tully_variables
      IMPLICIT none
#include "Molcas.fh"
#include "surfacehop.fh"
#include "stdalloc.fh"
      character*180 Get_Ln, Key, Line
      logical Found
      integer LuSpool, isfreeunit, NSTATE,i,j,ndata,maxHop
      real*8  temp(mxroot)
      real*8, allocatable :: AmatrixVR(:),AmatrixVI(:)
      external Get_Ln, isfreeunit

      LuSpool = isfreeunit(21)
      CALL SpoolInp(LuSpool)

      CALL RdNLst(LuSpool,'surfacehop')

      do
      Key = Get_Ln(LuSpool)
      Line = Key
      CALL UpCase(Line)

      select case(Line(1:4))
        case('TULL')
             tullyL=.true.
        case('DECO')
             Line = Get_Ln(LuSpool)
             CALL Get_F1(1,DECO)
             decoherence=.true.
        case('SUBS')
             Line = Get_Ln(LuSpool)
             CALL Get_I1(1,NSUBSTEPS)
        case('ETHR')
             Line = Get_Ln(LuSpool)
             CALL Get_F1(1,Ethreshold)
        case('RTHR')
             Line = Get_Ln(LuSpool)
             CALL Get_F1(1,RandThreshold)
        case('PSUB')
             tullySubVerb=.true.
        case('DMTX')
             Line = Get_Ln(LuSpool)
             CALL Get_I1(1,NSTATE)
             call mma_allocate(AmatrixVR,NSTATE*NSTATE,
     &                         label='AmatrixVR')
             call mma_allocate(AmatrixVI,NSTATE*NSTATE,
     &                         label='AmatrixVI')
             do i=1, NSTATE
                Line = Get_Ln(LuSpool)
                CALL Get_F(1,temp,NSTATE)
                do j=1, NSTATE
                   AmatrixVR((i-1)*NSTATE+j)=temp(j)
                end do
             end do
             do i=1, NSTATE
                Line = Get_Ln(LuSpool)
                CALL Get_F(1,temp,NSTATE)
                do j=1, NSTATE
                   AmatrixVI((i-1)*NSTATE+j)=temp(j)
                end do
             end do
             call Qpg_zArray('AmatrixV', Found, ndata)
             IF (.NOT.Found) THEN
               call Put_dArray('RAmatrixV', AmatrixVR, NSTATE*NSTATE)
               call Put_dArray('IAmatrixV', AmatrixVI, NSTATE*NSTATE)
             else if (ndata.ne.NSTATE*NSTATE) then
               call warningmessage(2,'The A matrix dimension on the '//
     &         'input is different from the one in runfile')
               CALL Abend()
             end if
             call mma_deallocate(AmatrixVR)
             call mma_deallocate(AmatrixVI)
        case('FRAN')
             Line = Get_Ln(LuSpool)
             CALL Get_F1(1,FixedRand)
             fixedrandL=.true.
        case('ISEE') ! initial seed number
             Line = Get_Ln(LuSpool)
             CALL Get_I1(1,InitSeed)
             iseedL=.true.
        case('MAXH')
             Line = Get_Ln(LuSpool)
             CALL Get_I1(1,maxHop)
             CALL Put_iScalar('MaxHopsTully',maxHop)
C             Write (6,*) 'MaxHops set to ', maxHop
        case('H5RE')
#ifdef _HDF5_
             lH5Restart = .True.
             Line = Get_Ln(LuSpool)
             CALL Get_S(1,FILE_H5RES,1)
#else
             write (6,*) 'The user asks to restart the dynamics '
             write (6,*) 'calculation from a HDF5 file, but this is not'
             write (6,*) ' supported in this installation.'
             call Quit_OnUserError()
#endif
        case('END ')
             exit
        case default
             Write (6,*) 'Unknown keyword: ', trim(Key)
             CALL Abend()
      end select
      end do
      CLOSE(LuSpool)
      RETURN
      END
