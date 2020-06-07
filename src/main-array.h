// =============================================================================
// =============================================================================
// ============================= Calculation Array =============================
// =============================================================================
// =============================================================================

class IET_arrays {
  #define MAX_CACHE_COUNT 6
  public:  // grid size
    int nv, nvm; int nz, ny, nx; int n_gvv, n_wvv, n_nhkvv, n_zeta; double drrism, drhi;
  public:  // grid scale vectors
    double dk_gvv, dk_wvv, dk_nhkvv, dk_zeta; double factor_sinfft_f_wvv, factor_sinfft_b_wvv, factor_sinfft_f_nhkvv, factor_sinfft_b_nhkvv, factor_sinfft_f_zeta, factor_sinfft_b_zeta;
    double dx, dy, dz; double dkx, dky, dkz; Vector box;
  public:  // const k-space functions
    __REAL__ *** gvv_data;
    __REAL__ *** gvv;                   // 1D gvv grid
    __REAL__ *** wvv, *** wvv_hlr;      // 1D wvv grid
    __REAL__ *** nhkvv, *** nhkvv_hlr;  // 1D nhkvv grid
    __REAL__ *** zeta;                  // 1D zetavv grid
    __REAL__ *** convolution_wvv, *** convolution_nhkvv, *** convolution_zeta;
    double ** wvv_zero_indicator, ** nhkvv_zero_indicator, ** zeta_zero_indicator;
    __REAL__ *** yukawa_kernel, *** ld_kernel;
  public:  // force field
    int frame_stamp;
    __REAL__ *** r2uvmin;   // nearest solute-solvent distance^2
    __REAL__ **** ulj;      // ulj: kJ/mol
    __REAL__ **** uljr;     // uljr: kJ/mol, used in RBC closures or sigma corrections
    __REAL__ *** ucoulsr, *** ucoullr;   // coul: kJ/mol/e
    __REAL__ **** Ecoul0;  // Original electrostatic field: Ecoul0[3]: kJ/mol/e/nm
    __REAL__ **** uuv, **** ulr, **** hlr, **** clr; // ulr + uuv = real_uuv
    //__REAL__ **** ucoul;
  public:  // RISM related
    //__REAL__ **** stuv;
    __REAL__ **** cuv;
    DIISNS::DIIS diis_rism;
  public:  // HI related
    __REAL__ **** phi, **** nphi;
    __REAL__ *** theta[MAX_SOL];
    __REAL__ **** dd, **** __dd;   // 3D dimensionless density
    HIEquationSolver solver_hi;
    DIISNS::DIIS diis_hi;
  public:  // RISM&HI cache
    __REAL__ **** rismhi_cache[MAX_CACHE_COUNT]; size_t n_extra_rismhi_cache[MAX_CACHE_COUNT];
        // rismhi_cache[0~1] are reserved, should never be reused
        // rismhi_cache[2~3] are temporarily used within one command, e.g. by closure calculation
        //   rismhi_cache[2] is also used for compression
        // rismhi_cache[4~5] are temporarily used within one super command, e.g. reversed-rism-merge
        // size of rismhi_cache[i] = N4*sizeof(__REAL__) + n_extra_rismhi_cache[i]
      __REAL__ * compress_buffer; size_t compress_buffer_size;
  public:  // HI&RISM cache pointer (note: CACHE pointers)
    __REAL__ **** res;
    DIISNS::DIIS * diis_current;
    __REAL__ **** ddpot_hi, **** huv;
  public:  // fftw related
    __REAL__ *** fftin;
    __REAL__ *** fftout;
    fftw_plan planf, planb;
  public:  // other internal variables
    RISMHI3D_FFTW_MP fftw_mp;
    RISMHI3D_FFSR_MP ffsr_mp;
  public:   // internal
    bool uuv_is_ready;
    bool is_energy_calculated, is_rdf_calculated;
  public:
    IET_Report * report_sites; IET_Report * report_total;
  public:
    void reset_for_calculation(bool clear_ff=true, bool clear_uuv=true, bool clear_h_c=true, bool clear_energy_statistics=false){
        size_t N3 = nx*ny*nz; size_t N4 = N3*nv; size_t N4m = N3*nvm;
      // clear forcefield
        if (clear_ff){
            if (r2uvmin) clear_tensor3d(r2uvmin, N3);
            if (ulj    ) clear_tensor4d(ulj    , N4);
            if (uljr   ) clear_tensor4d(uljr   , N4);
            if (ucoulsr) clear_tensor3d(ucoulsr, N3);
            if (ucoullr) clear_tensor3d(ucoullr, N3);
            if (Ecoul0)  clear_tensor4d(Ecoul0 , N3*3);
            uuv_is_ready = false;
            is_energy_calculated = false;
        }
        is_rdf_calculated = false;
      // clear RISM preparation cache
        if (clear_h_c){
            if (uuv) clear_tensor4d(uuv, N4);
            if (ulr) clear_tensor4d(ulr, N4);
            if (hlr) clear_tensor4d(hlr, N4);
            if (clr) clear_tensor4d(clr, N4);
            for (int i=0; i<MAX_CACHE_COUNT; i++) if (rismhi_cache[i]) clear_tensor4d(rismhi_cache[i], N4);
        }
      // clear RISM cache
        if (clear_h_c){
            if (cuv) clear_tensor4d(cuv , N4);
            diis_rism.ndiis = 0;
        }
      // clear HI cache
        if (clear_h_c){
            if (phi) clear_tensor4d(phi, N4m);
            if (nphi) clear_tensor4d(nphi, N4m);
            for (int i=0; i<MAX_SOL; i++) if (theta[i]) clear_tensor3d(theta[i], N3);
            if (dd) clear_tensor4d(dd, N4m);
            diis_hi.ndiis = 0;
        }
      // clear fftw cache
        clear_tensor3d(fftin, N3);
        clear_tensor3d(fftout, N3);
      // other clear
        for (int iv=0; iv<=nv; iv++) memset(&report_sites[iv], 0, sizeof(IET_Report));
    }
    bool set_scales(IET_Param * sys, Vector _box){
        box = _box;
        if (nx<=0 || ny<=0 || nz<=0) return false;
        if (box.x<=0 || box.y<=0 || box.z<=0) return false;
        drrism = sys->drrism; drhi = sys->drhi;
        //if (dr<=0 || n_wvv<=0 || n_nhkvv<=0) return false;
        factor_sinfft_f_wvv = factor_sinfft_f_nhkvv = 2*PI * drrism;
        factor_sinfft_f_zeta = 2*PI * drhi;
        factor_sinfft_b_wvv = drrism>0&&n_wvv>0? 1 / (2 * PI * drrism * (2*n_wvv+2)) : 1;
        factor_sinfft_b_nhkvv = drrism>0&&n_nhkvv>0? 1 / (2 * PI * drrism * (2*n_nhkvv+2)) : 1;
        factor_sinfft_b_zeta = drhi>0&&n_zeta>0? 1 / (2 * PI * drhi * (2*n_zeta+2)) : 1;
        dx = box.x / nx; dkx = 2*PI / (nx * dx);
        dy = box.y / ny; dky = 2*PI / (ny * dx);
        dz = box.z / nz; dkz = 2*PI / (nz * dx);

        #ifdef _EXPERIMENTAL_
            solver_hi.set_param(sys->lse_a, sys->lse_b, &sys->ex, 0, 2, true);
        #else
            solver_hi.set_param(sys->lse_a, sys->lse_b, 0, 2, true);
        #endif

        fftw_mp.set_scale(_box);
        uuv_is_ready = false;
        is_energy_calculated = false;
        is_rdf_calculated = false;
        return true;
    }
    static __REAL__ *** debug_trace_init3d(IET_Param * sys, const char * title, __REAL__ *** in){
        char buffer[64];
        if (sys->debug_level>=3) fprintf(sys->log(), "DEBUG:: init_tensor3d(%s) (%s)\n", title, print_memory_value(buffer, sizeof(buffer), _memory_total));
        return in;
    }
    static __REAL__ **** debug_trace_init4d(IET_Param * sys, const char * title, __REAL__ **** in){
        char buffer[64];
        if (sys->debug_level>=3) fprintf(sys->log(), "DEBUG:: init_tensor4d(%s) (%s)\n", title, print_memory_value(buffer, sizeof(buffer), _memory_total));
        return in;
    }
    void debug_trace_memory_allocation(IET_Param * sys, const char * title){
        if (sys->debug_level>=3) fprintf(sys->log(), "DEBUG:: %s\n", title);
    }
    __REAL__ **** submit_ram_N4(IET_Param * sys, __REAL__ **** & pointer){
        if (!pointer) pointer = debug_trace_init4d(sys, "ulj", init_tensor4d(nv, nz, ny, nx, 0));
        if (!pointer) exit_on_memalloc_failure();
        return pointer;
    }
    __REAL__ **** submit_cache(IET_Param * sys, int cache_id){
        if (cache_id<0 || cache_id>=MAX_CACHE_COUNT) exit_on_pointer_overflow();
        if (!rismhi_cache[cache_id]) submit_ram_N4(sys, rismhi_cache[cache_id]);
        return rismhi_cache[cache_id];
    }
    void alloc(IET_Param * sys, int _nv, int _nvm, int _nx, int _ny, int _nz, int nsolute){
        nv = _nv; nvm = _nvm; nz = _nz; ny = _ny; nx = _nx;
        gvv_data = gvv = wvv = wvv_hlr = nhkvv = nhkvv_hlr = zeta = yukawa_kernel = ld_kernel = nullptr; // cvvlj = cvvb = cvvc = cvvcb = nullptr;
        convolution_wvv = convolution_nhkvv = convolution_zeta = nullptr;
        box = Vector(0, 0, 0); // xvvi = nullptr
        frame_stamp = -1;
        ulj = debug_trace_init4d(sys, "ulj", init_tensor4d(nv, nz, ny, nx, 0));
        uljr = nullptr; if (sys->cmd_rbc_ljr_allowed){
            uljr = debug_trace_init4d(sys, "uljr", init_tensor4d(nv, nz, ny, nx, 0));
        }
        ucoulsr = debug_trace_init3d(sys, "coulsr", init_tensor3d<__REAL__>(nz, ny, nx, 0));
        ucoullr = debug_trace_init3d(sys, "coullr", init_tensor3d<__REAL__>(nz, ny, nx, 0));
        //ulpbe   = init_tensor3d<__REAL__>(nz, ny, nx, 0);
        //ucoulrf = sys->coul_a_reaction_field? init_tensor3d<__REAL__>(nz, ny, nx, 0) : nullptr;
        if (sys->b_allow_Ecoul0) Ecoul0 = debug_trace_init4d(sys, "Ecoul0", init_tensor4d(3, nz, ny, nx, 0));
        if (sys->b_allow_r2uvmin) r2uvmin = debug_trace_init3d(sys, "r2uvmin", init_tensor3d<__REAL__>(nz, ny, nx, 0)); else r2uvmin = nullptr;
        //stuv= init_tensor4d(nv, nz, ny, nx, 0);
        cuv = debug_trace_init4d(sys, "cuv", init_tensor4d(nv, nz, ny, nx, 0));
        hlr = debug_trace_init4d(sys, "hlr", init_tensor4d(nv, nz, ny, nx, 0));
        clr = debug_trace_init4d(sys, "clr", init_tensor4d(nv, nz, ny, nx, 0));
        uuv = debug_trace_init4d(sys, "uuv", init_tensor4d(nv, nz, ny, nx, 0));
        ulr = debug_trace_init4d(sys, "ulr", init_tensor4d(nv, nz, ny, nx, 0));
        //ucoul=init_tensor4d(nv, nz, ny, nx, 0);
        __dd = debug_trace_init4d(sys, "dd", init_tensor4d(nvm+1, nz, ny, nx, 0)); dd = nullptr;
        phi = init_tensor4d(nvm, nz, ny, nx, 0);
        nphi = debug_trace_init4d(sys, "nphi", init_tensor4d(nvm, nz, ny, nx, 0));
        for (int i=0; i<MAX_CACHE_COUNT; i++) n_extra_rismhi_cache[i] = 0;
          n_extra_rismhi_cache[2] += (((size_t)nv)*nx*ny*nz / sys->compress_page_size + 1) * sys->compress_page_size * sizeof(unsigned int);;
          for (int i=0; i<MAX_CACHE_COUNT; i++) rismhi_cache[i] = debug_trace_init4d(sys, "cache", init_tensor4d(nv, nz, ny, nx, n_extra_rismhi_cache[i]));
          compress_buffer = &rismhi_cache[2][0][0][0][0]; compress_buffer_size = sizeof(__REAL__)*nx*ny*nz*nv + n_extra_rismhi_cache[2];
          res = rismhi_cache[0]; ddpot_hi = huv = rismhi_cache[1];
        for (int i=0; i<nvm; i++) theta[i] = debug_trace_init3d(sys, "theta", init_tensor3d<__REAL__>(nz, ny, nx, 0));
        fftin = debug_trace_init3d(sys, "fftin", init_tensor3d<__REAL__>(nz, ny, nx, 0));
        fftout = debug_trace_init3d(sys, "fftout", init_tensor3d<__REAL__>(nz, ny, nx, 0));

        report_sites = (IET_Report*) memalloc(sizeof(IET_Report)*(nv+1));
        report_total = &report_sites[nv];

      #ifdef _LOCALPARALLEL_
        char debug_output_title[128];
        snprintf(debug_output_title, sizeof(debug_output_title), "allocating memory for ff x%d", sys->nt); debug_trace_memory_allocation(sys, debug_output_title);
        ffsr_mp.init(sys->nt, __mp_tasks, nx, ny, nz, nv, Ecoul0?true:false, r2uvmin?true:false, uljr?true:false);
        snprintf(debug_output_title, sizeof(debug_output_title), "allocating memory for fftw x%d", sys->nt); debug_trace_memory_allocation(sys, debug_output_title);
        fftw_mp.init(sys->nt, sys->ntf, __mp_tasks, nx, ny, nz, sys->xvv_k_shift);
      #else
      #endif
      #ifdef _FFTWMPPARALLEL_
        fftw_plan_with_nthreads(sys->ntf);
      #endif
        planf = fftw_plan_r2r_3d(nz, ny, nx, &fftin[0][0][0], &fftout[0][0][0], (fftw_r2r_kind)FFTW_FORWARD, (fftw_r2r_kind)FFTW_FORWARD, (fftw_r2r_kind)FFTW_FORWARD, FFTW_ESTIMATE);
      #ifdef _FFTWMPPARALLEL_
        fftw_plan_with_nthreads(sys->ntf);
      #endif
        planb = fftw_plan_r2r_3d(nz, ny, nx, &fftout[0][0][0], &fftin[0][0][0], (fftw_r2r_kind)FFTW_BACKWARD, (fftw_r2r_kind)FFTW_BACKWARD, (fftw_r2r_kind)FFTW_BACKWARD, FFTW_ESTIMATE);

        solver_hi.init(1000000);


        long int n_grid_diis_rism = sys->nv * sys->nr[0] * sys->nr[1] * sys->nr[2];
        long int n_grid_diis_hi = sys->nvm * sys->nr[0] * sys->nr[1] * sys->nr[2];
//printf("diis require buffer: %lu %lu -> %lu\n", n_diis_buffer_rism, n_diis_buffer_hi, n_diis_buffer);
        if (sys->ndiis_rism>1) diis_rism.init(sys->ndiis_rism, sys->nv * sys->nr[0] * sys->nr[1] * sys->nr[2], sys->nt);
        if (sys->ndiis_hi>1) diis_hi.init(sys->ndiis_hi, sys->nvm * sys->nr[0] * sys->nr[1] * sys->nr[2], sys->nt);
        diis_current = nullptr;
        uuv_is_ready = false;
        is_energy_calculated = false;
        is_rdf_calculated = false;
    }
  public:
    void calculate_local_density_from_guv(IET_Param * sys, __REAL__ **** ld_density, __REAL__ **** temp, bool do_convolution=true){
        size_t N3 = nx*ny*nz;
        size_t N4 = N3 * sys->nv; size_t N4m = N3 * sys->nvm;
        double grid = 2*1.732050808 / Vector(nx, ny, nz).mod();

        double atom_density[MAX_PATH]; double mass[MAX_PATH]; for (int i=0; i<sys->nv; i++) atom_density[i] = mass[i] = 0;
        for (int iv=0; iv<sys->nv; iv++) mass[sys->av[iv].iaa] += sys->av[iv].mass;
        for (int iv=0; iv<sys->nv; iv++) atom_density[iv] =  sys->av[iv].mass / mass[sys->av[iv].iaa];

        clear_tensor4d(temp, N4);
        for (int iv=0; iv<sys->nv; iv++){
            int ivm = sys->av[iv].iaa;
            for (size_t i3=0; i3<N3; i3++) temp[iv][0][0][i3] += atom_density[iv] * (1+huv[iv][0][0][i3]) * (dd? dd[ivm][0][0][i3] : sys->nbulk[ivm]);
        }
        if (do_convolution){
            for (size_t i4m=0; i4m<N4; i4m++) temp[0][0][0][i4m] = temp[0][0][0][i4m] * (1+grid) - grid;
                // the above step is to make sure zero density regions are perfectly preserved
            perform_3rx1k_convolution(&fftw_mp, temp, nx, ny, nz, box, sys->nvm, sys->nvm, ld_kernel, dk_nhkvv, sys->xvv_k_shift, n_nhkvv, ld_density, fftin, fftout, planf, planb, true);
            for (int ivm=0; ivm<sys->nvm; ivm++) for (size_t i3=0; i3<N3; i3++){
                if (ld_density[ivm][0][0][i3]<0){
                    ld_density[ivm][0][0][i3] = 0;
                }
            }
        } else {
            for (size_t i4m=0; i4m<N4; i4m++) ld_density[0][0][0][i4m] = temp[0][0][0][i4m];
        }
    }
};


inline double atomised_molecule_uv_potential(IET_Param * sys, IET_arrays * arr, int iv, size_t i3){
    double U = 0;
    for (int iia=sys->vmmapi[iv][0]; iia<sys->vmmapi[iv][1]; iia++){ int ia = sys->vmmap[iia]; U += arr->uuv[ia][0][0][i3]; }
    return U;
}

inline double atomised_molecule_uv_potential(IET_Param * sys, IET_arrays * arr, int iv, int ix, int iy, int iz){
    double U = 0;
    for (int iia=sys->vmmapi[iv][0]; iia<sys->vmmapi[iv][1]; iia++){ int ia = sys->vmmap[iia]; U += arr->uuv[ia][iz][iy][ix]; }
    return U;
    /*
    for (int iia=sys->vmmapi[iv][0]; iia<sys->vmmapi[iv][1]; iia++){ int ia = sys->vmmap[iia]; U += arr->ulj[ia][iz][iy][ix] + (arr->ucoulsr[iz][iy][ix] + arr->ucoullr[iz][iy][ix]) * sys->av[ia].charge; }
    return U / sys->temperature * sys->default_temperature;
    */
}

void build_ld_kernel(IET_Param * sys, __REAL__ *** ld_kernel_pt, int nvm, int xvv_length, double ld_kernel_scale_factor){
    double dr = sys->drrism;
    double dk = PI / (1+xvv_length) / dr;   //printf("drrism: %g\n", sys->drrism);

    clear_tensor3d(ld_kernel_pt, nvm*nvm*xvv_length);
    for (int ivm=0; ivm<sys->nvm; ivm++){
        //double density = sys->bulk_density_mv[sys->av[ivm].iaa];
        double density = sys->bulk_density_mv[ivm] / (ld_kernel_scale_factor*ld_kernel_scale_factor*ld_kernel_scale_factor);
        //double rmol = pow(fabs(1.0/sys->bulk_density_mv[ivm]/4.0*3.0/PI), 1.0/3);
        double rmol = pow(fabs(1.0/density/4.0*3.0/PI), 1.0/3);

        for (int i=0; i<xvv_length; i++){
            //double k = (i+1)*dk; if (k<MACHINE_REASONABLE_ERROR) k = MACHINE_REASONABLE_ERROR;
            double k = i*dk;
            double ksigma = k*rmol;
            ld_kernel_pt[ivm][ivm][i] = k==0? 1 : (4*PI*(sin(ksigma) - ksigma*cos(ksigma))/k/k/k)*density;
        }
    }
}
