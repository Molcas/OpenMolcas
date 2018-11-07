* $ this file belongs to the Molcas repository $
      program main
#ifdef _FPE_TRAP_
      Use, Intrinsic :: IEEE_Exceptions
#endif
      Integer iReturn
#ifdef _FPE_TRAP_
      Call IEEE_Set_Halting_Mode(IEEE_Usual,.True._4)
#endif

      Call Start('single_aniso')
      Call single_aniso(iReturn)
      Call Finish(iReturn)
      End
