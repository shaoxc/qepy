## Not supported
### cmd and files
  - cppp.x : CPV/src/cppp.f90
  - cp.x : CPV/src/cprstart.f90
  - manycp.x : CPV/src/manycp.f90
  - wfdd.x : CPV/src/wfdd.f90
  - bse_main.x : GWW/bse/bse_main.f90
  - head.x : GWW/head/head.f90
  - phcg.x : PHonon/Gamma/phcg.f90
  - pwi2xsf.x : PW/tools/pwi2xsf.f90
  - molecularnexafs.x : XSpectra/src/molecularnexafs.f90
  - spectra_correction.x : XSpectra/src/spectra_correction.f90
  - xspectra.x : XSpectra/src/xspectra.f90

### Conflict
  - CPV/src/libcp.a:nlcc.f90:force_cc and  PW/src/libpw.a:force_cc.f90
  - GWW/pw4gww/libpw4gww.a:write_wannier_matrix.f90:read_wannier_matrix and GWW/bse/libbse.a:write_wannier_matrix.f90
  - GWW/simple_bse/libsimple_exc.a:diago_exc_cg.f90:minparabola and CPV/src/libcp.a:cglib.f90
  - PHonon/Gamma/libphcg.a:solve_ph.f90:set_asr and LR_Modules/liblrmod.a:dynmat_sub.f90
  - atomic/src/green.o:green.f90:green and XSpectra/src/xspectra.o:xspectra.f90
  - GWW/head/lanczos_k.o:lanczos_k.f90:h_psi_scissor and GWW/pw4gww/libpw4gww.a:pola_lanczos.f90
  - GWW/head/* and PHonon/PH/libph.a (bcast_ph_input, close_phq, openfilq, phq_readin)
  - GWW/pw4gww/libpw4gww.a:lanczos_chains.f90:lanczos and GWW/gww/libgww.a:lanczos_polarization.f90

### Ignore library
  - CPV/src/libcp.a
  - GWW/bse/libbse.a
  - PHonon/Gamma/libphcg.a
  - GWW/gww/libgww.a
