// =============================================================================
// =============================================================================
// ============ Energy reports and the calculation of free energies ============
// =============================================================================
// =============================================================================
class EnergyReport {
  public:
    double mass, mass_mol, density;
    double lj, coulsr, coullr, entropy, N, N0, dN, Ng, dNg, pmv;
    double Chandler_Guv[2], zeta[4], cuv, clr, chuv, Uef0, Uef1;
  public:
    void operator += (EnergyReport & o){
      // potential energy
        lj += o.lj; coulsr += o.coulsr; coullr += o.coullr;
        entropy += o.entropy * o.mass / o.mass_mol;
        N += o.N * o.mass / o.mass_mol; N0 += o.N0 * o.mass / o.mass_mol; dN = N - N0;
        pmv += o.pmv * o.mass / o.mass_mol;
        Ng += o.Ng * o.mass / o.mass_mol; dNg = Ng - N0;
      // Chandler's
        Chandler_Guv[0] += o.Chandler_Guv[0];
        Chandler_Guv[1] += o.Chandler_Guv[1] * o.mass / o.mass_mol;
      // SSYBG's
        zeta[0] += o.zeta[0] * o.mass / o.mass_mol;
        zeta[1] += o.zeta[1];
        zeta[2] += o.zeta[2] * o.mass / o.mass_mol;
        zeta[3] += o.zeta[3];
      // other
        cuv += o.cuv;// * o.mass / o.mass_mol;
        clr += o.clr;// * o.mass / o.mass_mol;
        chuv += o.chuv;// * o.mass / o.mass_mol;
        Uef0 += o.Uef0;
        Uef1 += o.Uef1;
    }
    void operator *= (double scaling){
        mass *= scaling; mass_mol *= scaling; density *= scaling;
        lj *= scaling; coulsr *= scaling; coullr *= scaling; entropy *= scaling;
        N *= scaling; N0 *= scaling; Ng *= scaling; dNg *= scaling;
        Chandler_Guv[0] *= scaling; Chandler_Guv[1] *= scaling;
        cuv *= scaling; clr *= scaling; Uef0 *= scaling; Uef1 *= scaling;
    }
};
double calculate_zeta_by_chuv(int closure, double factor, double uuv, double huv, double cuv, double hlr){
    double chuv = huv - cuv; double t_over_ch = 0; double s, t;
    switch (closure) {
        case CLOSURE_HNC            : t_over_ch = 1; break;
        case CLOSURE_PLHNC          : t_over_ch = 1; break;
        case CLOSURE_MHNC           : t_over_ch = 1; break; // not defined
        case CLOSURE_MSA            : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : (ln(1-uuv+chuv) + uuv)/chuv; break;
        case CLOSURE_KGK            : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : (ln(1-uuv+chuv) + uuv)/chuv; break;
        case CLOSURE_PY             : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : ln(1+chuv)/chuv; break;
        case CLOSURE_D2             : t_over_ch = 1 - chuv/2; break;
        case CLOSURE_HNCB           : t_over_ch = 1; break;  // undefined
        //case CLOSURE_KH             : t_over_ch = -uuv+chuv<=0? 1 : chuv==0? 1 : (ln(1-uuv+chuv) + uuv)/chuv; break;
        case CLOSURE_KH             : t_over_ch = 1; break;
        case CLOSURE_PSE2           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2)/chuv; break;
        case CLOSURE_PSE3           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6)/chuv; break;
        case CLOSURE_PSE4           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24)/chuv; break;
        case CLOSURE_PSE5           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120)/chuv; break;
        case CLOSURE_PSE6           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120+chuv*chuv*chuv*chuv*chuv*chuv/720)/chuv; break;
        case CLOSURE_PSE7           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120+chuv*chuv*chuv*chuv*chuv*chuv/720+chuv*chuv*chuv*chuv*chuv*chuv*chuv/5040)/chuv; break;
        case CLOSURE_PSE8           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120+chuv*chuv*chuv*chuv*chuv*chuv/720+chuv*chuv*chuv*chuv*chuv*chuv*chuv/5040+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/40320)/chuv; break;
        case CLOSURE_PSE9           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120+chuv*chuv*chuv*chuv*chuv*chuv/720+chuv*chuv*chuv*chuv*chuv*chuv*chuv/5040+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/40320+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/362880)/chuv; break;
        case CLOSURE_PSE10          : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120+chuv*chuv*chuv*chuv*chuv*chuv/720+chuv*chuv*chuv*chuv*chuv*chuv*chuv/5040+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/40320+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/362880+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/3628800)/chuv; break;
        case CLOSURE_MS             : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : sqrt(1+2*chuv)/chuv; break;
        case CLOSURE_BPGGHNC        : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : pow(1+factor*chuv, 1.0/factor)/chuv; break;
        case CLOSURE_VM             : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : 1 - chuv/2/(1+factor*chuv); break;
        case CLOSURE_MP             : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : ((1+factor)*exp(chuv/(1+factor)) - factor)/chuv; break;
        default:
          #ifdef _EXPERIMENTAL_
            t_over_ch = experimental_calculate_t_over_ch_by_chuv(closure, factor, uuv, huv, cuv, hlr);
          #else
            t_over_ch = 1;
          #endif
            break;
    }
    return t_over_ch * cuv;
}

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
    __REAL__ *** gvv;               // 1D gvv grid
    __REAL__ *** wvv, *** wvv_hlr, *** wvv_save;        // 1D wvv grid
    __REAL__ *** nhkvv, *** nhkvv_hlr, *** nhkvv_save;  // 1D nhkvv grid
    __REAL__ *** zeta;              // 1D zetavv grid
    __REAL__ *** yukawa_kernel, *** local_coulomb_kernel, *** ld_kernel;
  public:  // force field
    int frame_stamp;
    __REAL__ *** r2uvmin;   // nearest solute-solvent distance^2
    __REAL__ **** ulj;      // ulj: kJ/mol
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
    __REAL__ *** guvm[MAX_SOL];
    __REAL__ **** dd, **** __dd;   // 3D dimensionless density
    __REAL__ *** pseudoliquid_potential;
    HIEquationSolver solver_hi;
    DIISNS::DIIS diis_hi;
  public:  // RISM&HI cache
    __REAL__ **** rismhi_cache[MAX_CACHE_COUNT];
        // rismhi_cache[0~1] are reserved for all time, should never be reused
        // rismhi_cache[2~3] are temporarily used within one command, e.g. by closure calculation
        // rismhi_cache[4~5] are temporarily used within one grand command, e.g. reversed-rism-merge
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
    EnergyReport total_energy;
    EnergyReport site_energy[MAX_SOL];
  public:
    void reset_for_calculation(bool clear_ff=true, bool clear_uuv=true, bool clear_h_c=true, bool clear_energy_statistics=false){
        size_t N3 = nx*ny*nz; size_t N4 = N3*nv; size_t N4m = N3*nvm;
      // clear forcefield
        if (clear_ff){
            if (r2uvmin) clear_tensor3d(r2uvmin, N3);
            if (ulj    ) clear_tensor4d(ulj    , N4);
            if (ucoulsr) clear_tensor3d(ucoulsr, N3);
            if (ucoullr) clear_tensor3d(ucoullr, N3);
            if (Ecoul0)  clear_tensor4d(Ecoul0 , N3*3);
            if (pseudoliquid_potential) clear_tensor3d(pseudoliquid_potential , N3);
            uuv_is_ready = false;
            is_energy_calculated = false;
        }
        is_rdf_calculated = false;
      // clear RISM preparation cache
        if (clear_h_c){
            if (uuv    ) clear_tensor4d(uuv, N4);
            if (ulr    ) clear_tensor4d(ulr, N4);
            if (hlr    ) clear_tensor4d(hlr, N4);
            if (clr    ) clear_tensor4d(clr, N4);
            for (int i=0; i<MAX_CACHE_COUNT; i++) if (rismhi_cache[i]) clear_tensor4d(rismhi_cache[i], N4);
        }
      // clear RISM cache
        if (clear_h_c){
            if (cuv    ) clear_tensor4d(cuv , N4);
            diis_rism.ndiis = 0;
        }
      // clear HI cache
        if (clear_h_c){
            if (phi    ) clear_tensor4d(phi, N4m);
            if (nphi   ) clear_tensor4d(nphi, N4m);
            for (int i=0; i<MAX_SOL; i++) if (guvm[i]) clear_tensor3d(guvm[i], N3);
            if (dd     ) clear_tensor4d(dd , N4m);
            diis_hi.ndiis = 0;
        }
      // clear fftw cache
        clear_tensor3d(fftin, N3);
        clear_tensor3d(fftout, N3);
      // other clear
        if (clear_ff || clear_uuv || clear_h_c){
            memset(&site_energy, 0, sizeof(site_energy));
        }
        if (clear_energy_statistics){
            memset(&total_energy, 0, sizeof(total_energy));
        }
    }
    bool set_scales(IET_Param * sys, Vector _box){
        box = _box;
        if (nx<=0 || ny<=0 || nz<=0) return false;
        if (box.x<=0 || box.y<=0 || box.z<=0) return false;
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
    void alloc(IET_Param * sys, int _nv, int _nvm, int _nx, int _ny, int _nz, int nsolute){
        nv = _nv; nvm = _nvm; nz = _nz; ny = _ny; nx = _nx;
        gvv_data = gvv = wvv = wvv_hlr = nhkvv = nhkvv_hlr = zeta = yukawa_kernel = local_coulomb_kernel = ld_kernel = nullptr; // cvvlj = cvvb = cvvc = cvvcb = nullptr;
        box = Vector(0, 0, 0); // xvvi = nullptr
        frame_stamp = -1;
        ulj = debug_trace_init4d(sys, "ulj", init_tensor4d(nv, nz, ny, nx, 0));
        ucoulsr = debug_trace_init3d(sys, "coulsr", init_tensor3d<__REAL__>(nz, ny, nx, 0));
        ucoullr = debug_trace_init3d(sys, "coullr", init_tensor3d<__REAL__>(nz, ny, nx, 0));
        //ulpbe   = init_tensor3d<__REAL__>(nz, ny, nx, 0);
        //ucoulrf = sys->coul_a_reaction_field? init_tensor3d<__REAL__>(nz, ny, nx, 0) : nullptr;
        Ecoul0 = debug_trace_init4d(sys, "Ecoul0", init_tensor4d(3, nz, ny, nx, 0));
        r2uvmin = debug_trace_init3d(sys, "r2uvmin", init_tensor3d<__REAL__>(nz, ny, nx, 0));
        //stuv= init_tensor4d(nv, nz, ny, nx, 0);
        cuv = debug_trace_init4d(sys, "cuv", init_tensor4d(nv, nz, ny, nx, 0));
        hlr = debug_trace_init4d(sys, "hlr", init_tensor4d(nv, nz, ny, nx, 0));
        clr = debug_trace_init4d(sys, "clr", init_tensor4d(nv, nz, ny, nx, 0));
        uuv = debug_trace_init4d(sys, "uuv", init_tensor4d(nv, nz, ny, nx, 0));
        ulr = debug_trace_init4d(sys, "ulr", init_tensor4d(nv, nz, ny, nx, 0));
        //ucoul=init_tensor4d(nv, nz, ny, nx, 0);
        __dd = debug_trace_init4d(sys, "dd", init_tensor4d(nvm+1, nz, ny, nx, 0)); dd = nullptr;
        pseudoliquid_potential = debug_trace_init3d(sys, "pseudoliquid_potential", init_tensor3d<__REAL__>(nz, ny, nx, 0));
        phi = init_tensor4d(nvm, nz, ny, nx, 0);
        nphi = debug_trace_init4d(sys, "nphi", init_tensor4d(nvm, nz, ny, nx, 0));
        for (int i=0; i<MAX_CACHE_COUNT; i++) rismhi_cache[i] = debug_trace_init4d(sys, "cache", init_tensor4d(nv, nz, ny, nx, 0));
          res = rismhi_cache[0]; ddpot_hi = huv = rismhi_cache[1];
        for (int i=0; i<nvm; i++) guvm[i] = debug_trace_init3d(sys, "guvm", init_tensor3d<__REAL__>(nz, ny, nx, 0));
        fftin = debug_trace_init3d(sys, "fftin", init_tensor3d<__REAL__>(nz, ny, nx, 0));
        fftout = debug_trace_init3d(sys, "fftout", init_tensor3d<__REAL__>(nz, ny, nx, 0));

      #ifdef _LOCALPARALLEL_
        char debug_output_title[128];
        snprintf(debug_output_title, sizeof(debug_output_title), "allocating memory for ff x%d paralleling", sys->nt); debug_trace_memory_allocation(sys, debug_output_title);
        ffsr_mp.init(sys->nt, sys->mp_tasks, nx, ny, nz, nv);
        snprintf(debug_output_title, sizeof(debug_output_title), "allocating memory for fftw x%d paralleling", sys->nt); debug_trace_memory_allocation(sys, debug_output_title);
        fftw_mp.init(sys->nt, sys->mp_tasks, nx, ny, nz, sys->xvv_k_shift);
      #else
      #endif
      #ifdef _FFTWMPPARALLEL_
        fftw_plan_with_nthreads(sys->ntf);
      #endif
        planf = fftw_plan_r2r_3d(nz, ny, nx, &fftin[0][0][0], &fftout[0][0][0], (fftw_r2r_kind)FFTW_FORWARD, (fftw_r2r_kind)FFTW_FORWARD, (fftw_r2r_kind)FFTW_FORWARD, FFTW_ESTIMATE);
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
    void fill_in_mass(IET_Param * sys){
        double molecule_mass[MAX_SOL]; memset(molecule_mass, 0, sizeof(molecule_mass)); for (int iv=0; iv<sys->nv; iv++) molecule_mass[sys->av[iv].iaa] += sys->av[iv].mass * sys->av[iv].multi;
        for (int iv=0; iv<nv; iv++){
            site_energy[iv].mass = sys->av[iv].mass * sys->av[iv].multi;
            site_energy[iv].mass_mol = molecule_mass[sys->av[iv].iaa];
            site_energy[iv].density = sys->density_mv[sys->av[iv].iaa];
        }
        total_energy.density = 0; for (int ivm=0; ivm<nvm; ivm++) total_energy.density += sys->density_mv[ivm];
        total_energy.mass = 0; for (int iv=0; iv<nv; iv++){ total_energy.mass += site_energy[iv].mass * site_energy[iv].density / total_energy.density; }; total_energy.mass_mol = total_energy.mass;
    }
    void calculate_energy(IET_Param * sys, double ccutoff, double dielect=1){
        size_t N3 = nx * ny * nz;
        double dV = box.x * box.y * box.z / N3;
        double beta = sys->default_temperature / sys->temperature; double total_mass = 0; double total_density = 0;

        memset(&total_energy, 0, sizeof(total_energy)); memset(&site_energy[0], 0, sizeof(site_energy));
        fill_in_mass(sys);
        double total_dielect = 0; for (int ivm=0; ivm<nvm; ivm++) total_dielect += sys->dielect_mol[ivm];

      // site related energies
        double hardsphere_guv_cutoff = exp(-sys->ucutoff_hs);
        for (int iv=0; iv<nv; iv++){
            int ivm = sys->av[iv].iaa; double q = sys->av[iv].charge_esp / dielect; double dN = dV * sys->density_av[iv];
            __REAL__ * uuv1 = &uuv[iv][0][0][0];
            __REAL__ * huv1 = &huv[iv][0][0][0];
            __REAL__ * hlr1 = &hlr[iv][0][0][0];
            __REAL__ * dd1 = dd? &dd[sys->av[iv].iaa][0][0][0] : nullptr; double nbulk = sys->nbulk[sys->av[iv].iaa];
            __REAL__ * cuv1 = &cuv[iv][0][0][0];
            __REAL__ * clr1 = &clr[iv][0][0][0];
            site_energy[iv].N0 = box.x * box.y * box.z * sys->density_av[iv];
            double factor = sys->closure_factors[iv]; double rhob = 0; int n_rhob = 0;
            for (size_t i3=0; i3<N3; i3++){
                double g = (1+huv1[i3]) * (dd1? dd1[i3]/nbulk : 1);
                double h = g - 1;
                double dn = g * dN * sys->av[iv].multi;
              // potential energy terms
                site_energy[iv].N += g * dN;
                if (dn>MACHINE_REASONABLE_ERROR) site_energy[iv].Ng += dN;
                site_energy[iv].lj += ulj[iv][0][0][i3] * dn * sys->scale_lj;
                site_energy[iv].coulsr += (ucoulsr[0][0][i3]) * q * dn * sys->scale_coul;
                site_energy[iv].coullr += (ucoullr[0][0][i3]) * q * dn * sys->scale_coul;
                site_energy[iv].entropy += dn<=0? 0 : - dn * log(g / sys->nbulk[ivm]) / beta;
              // Chandler's energy terms
                double Chandler_Guv_this = 0;
                //double this_multi_density = sys->density_av[iv] * sys->av[iv].multi;

                //Chandler_Guv_this = calculate_Chandler_Guv(sys->closures[iv], sys->closure_factors[iv], uuv1[i3], huv1[i3], hlr1[i3], cuv1[i3], clr1[i3]) * (dd1? dd1[i3] : 1) * sys->temperature/sys->default_temperature;

                //if (huv1[i3]+1>hardsphere_guv_cutoff) site_energy[iv].Chandler_Guv[0] += dN * sys->av[iv].multi * Chandler_Guv_this;
                //site_energy[iv].Chandler_Guv[1] += dN * sys->av[iv].multi * Chandler_Guv_this;

                Chandler_Guv_this = (- (0.5*huv1[i3] + 1) * cuv1[i3]) * dN * (dd?dd1[i3]/nbulk : 1) / beta;
//Chandler_Guv_this = huv1[i3]+1>hardsphere_guv_cutoff ? 0.5*huv1[i3]* cuv1[i3] * dN * (dd?dd1[i3]/nbulk : 1) / beta : 0;
                site_energy[iv].Chandler_Guv[0] += Chandler_Guv_this;
                site_energy[iv].Chandler_Guv[1] += Chandler_Guv_this;

              // internal calculation
                //if (huv1[i3]+1>hardsphere_guv_cutoff) site_energy[iv].cuv += (cuv1[i3]) * dN;
                //if (huv1[i3]+1>hardsphere_guv_cutoff) site_energy[iv].clr += (clr1[i3]) * dN;
                site_energy[iv].cuv += (cuv1[i3]) * dN; site_energy[iv].clr += (clr1[i3]) * dN;
                //if (huv1[i3]+1>hardsphere_guv_cutoff) site_energy[iv].clr += (ln(huv1[i3]+1)*beta + uuv1[i3])*dn;
                double zeta_this = calculate_zeta_by_chuv(sys->closures[iv], sys->closure_factors[iv], uuv1[i3], huv1[i3], cuv1[i3], hlr1[i3]);
                if (dn>sys->gcutoff_liquid_occupation){
                    site_energy[iv].zeta[0] += zeta_this * dN * (dd1? dd1[i3]:nbulk) ;
                    site_energy[iv].zeta[1] += zeta_this * dN * sys->av[iv].multi * (dd1? dd1[i3]:nbulk) ;
                    site_energy[iv].zeta[2] += clr[iv][0][0][i3] * dN * (dd1? dd1[i3]:nbulk) ;
                    site_energy[iv].zeta[3] += clr[iv][0][0][i3] * dN * sys->av[iv].multi * (dd1? dd1[i3]:nbulk) ;
                }
                //if (huv1[i3]+1>hardsphere_guv_cutoff) site_energy[iv].chuv += cuv1[i3]*dN + 0.5*(cuv1[i3]) * (dn - dN);
                if (huv1[i3]+1>hardsphere_guv_cutoff) site_energy[iv].chuv += (huv1[i3] - cuv1[i3])*dN;
                if (g>sys->gcutoff_liquid_occupation){
                    double this_n = 1; if (dd){ this_n = 0;
                        for (int ivm=0; ivm<nvm; ivm++) this_n += dd[ivm][0][0][i3] / sys->dielect_mol[ivm];
                        this_n = 1 - this_n; if (this_n<0) this_n = 0;
                    }
                    Vector Eef0 = Vector(Ecoul0[0][0][0][i3], Ecoul0[1][0][0][i3], Ecoul0[2][0][0][i3]);
                    double Eef02 = Eef0.pow2();
                    double this_Uef0 = dV * (0.5/(4*PI*COULCOOEF)) * Eef02 / nv * sys->scale_coul;
                    site_energy[iv].Uef0 += this_Uef0;
                    site_energy[iv].Uef1 += this_n * this_Uef0;
                }
                if (r2uvmin&&r2uvmin[0][0][i3]<MACHINE_REASONABLE_ERROR){ n_rhob ++; rhob += g; }
            }
            site_energy[iv].dN = site_energy[iv].N - site_energy[iv].N0;
            site_energy[iv].dNg = site_energy[iv].Ng - site_energy[iv].N0;
            //site_energy[iv].pmv = site_energy[iv].dN/sys->bulk_density_mv[sys->av[iv].iaa] - (n_rhob>0? (rhob/n_rhob-1)*(box.x*box.y*box.z) : 0);
            site_energy[iv].pmv = (site_energy[iv].N/(n_rhob>0?rhob/n_rhob:1) - site_energy[iv].N0)/ sys->bulk_density_mv[sys->av[iv].iaa];
            total_energy += site_energy[iv];
        }
      //
    }
    double solvation_energy(IET_Param * sys, double ccutoff, double dielect=1){
        calculate_energy(sys, ccutoff, dielect);
        return total_energy.lj + total_energy.coulsr + total_energy.coullr;
    }
  public:
    void display_solvation_energy(IET_Param * sys, FILE * flog, FILE * flog2, const char * prefix = "", const char * total_prefix="total"){
        FILE * fout[2] = { flog, flog2 };
        for (int io=0; io<2; io++) if (fout[io]) fprintf(fout[io], "%s%7s %10s %10s %10s %10s %10s %10s\n", prefix, "Atom", "mass", "DeltaN", "-TS", "LJSR", "Coulomb", "Euv");
        for (int iv=0; iv<nv; iv++){
            char mass_string[32];
            snprintf(mass_string, sizeof(mass_string), "%.0fx%d", sys->av[iv].mass, sys->av[iv].multi);
            for (int io=0; io<2; io++) if (fout[io]) if (nv>1) fprintf(fout[io], "%s%7s %10s %10.4g %10.4g %10.4g %10.4g %10.4g\n", prefix, sys->av[iv].name, mass_string, site_energy[iv].dN, site_energy[iv].entropy, site_energy[iv].lj, site_energy[iv].coulsr+site_energy[iv].coullr, site_energy[iv].lj+site_energy[iv].coulsr+site_energy[iv].coullr);
        }
        for (int io=0; io<2; io++) if (fout[io]) fprintf(fout[io], "%s%7s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n", prefix, total_prefix, total_energy.mass, total_energy.dN, total_energy.entropy, total_energy.lj, total_energy.coulsr+total_energy.coullr, total_energy.lj+total_energy.coulsr+total_energy.coullr);
    }
    void display_solvation_energy_full(IET_Param * sys, FILE * flog, FILE * flog2, const char * prefix = "", const char * total_prefix="total", bool b_title=true, bool b_detail=true, bool b_sum=true){
        FILE * fout[2] = { flog, flog2 };
        if (b_title) for (int io=0; io<2; io++) if (fout[io]) fprintf(fout[io], "%s%7s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", prefix, "Atom", "mass", "DN", "DN_vac", "-TS", "LJSR", "Coulomb", "Uef1", "Uef0", "PMV");
        if (b_detail) for (int iv=0; iv<nv; iv++){
            char mass_string[32];
            snprintf(mass_string, sizeof(mass_string), "%.0fx%d", sys->av[iv].mass, sys->av[iv].multi);
            for (int io=0; io<2; io++) if (fout[io]) if (nv>1) fprintf(fout[io], "%s%7s %10s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n", prefix, sys->av[iv].name, mass_string, site_energy[iv].dN, site_energy[iv].dNg, site_energy[iv].entropy, site_energy[iv].lj, site_energy[iv].coulsr+site_energy[iv].coullr, site_energy[iv].Uef1,site_energy[iv].Uef0, site_energy[iv].pmv);
        }
        if (b_sum) for (int io=0; io<2; io++) if (fout[io]) fprintf(fout[io], "%s%7s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n", prefix, total_prefix, total_energy.mass, total_energy.dN, total_energy.dNg, total_energy.entropy, total_energy.lj, total_energy.coulsr+total_energy.coullr, total_energy.Uef1, total_energy.Uef0, total_energy.pmv);
    }
    void display_solvation_energy_ef(IET_Param * sys, FILE * flog, FILE * flog2, const char * prefix = "", const char * total_prefix="total"){
        FILE * fout[2] = { flog, flog2 };
        for (int io=0; io<2; io++) if (fout[io]) fprintf(fout[io], "%s%7s %7s %7s %10s %10s %10s %10s\n", prefix, "Atom", "mass", "DeltaN", "-TS", "LJSR", "Uef0", "Euv");
        for (int iv=0; iv<nv; iv++){
            char mass_string[32];
            snprintf(mass_string, sizeof(mass_string), "%.0fx%d", sys->av[iv].mass, sys->av[iv].multi);
            for (int io=0; io<2; io++) if (fout[io]) if (nv>1) fprintf(fout[io], "%s%7s %7s %7.2f %10.4g %10.4g %10.4g %10.4g\n", prefix, sys->av[iv].name, mass_string, site_energy[iv].dN, site_energy[iv].entropy, site_energy[iv].lj, site_energy[iv].Uef0, site_energy[iv].lj-2*site_energy[iv].Uef0*(1-1/sys->mean_dielect));
        }
        for (int io=0; io<2; io++) if (fout[io]) fprintf(fout[io], "%s%7s %7.2f %7.2f %10.4g %10.4g %10.4g %10.4g\n", prefix, total_prefix, total_energy.mass, total_energy.dN, total_energy.entropy, total_energy.lj, total_energy.Uef0, total_energy.lj-2*total_energy.Uef0*(1-1/sys->mean_dielect));
    }
  public:
    void display_solvation_correlations(IET_Param * sys, FILE * flog, FILE * flog2, const char * prefix = ""){
        FILE * fout[2] = { flog, flog2 };
        for (int io=0; io<2; io++) if (fout[io]) fprintf(fout[io], "%s%7s %8s %8s %8s %8s %8s %8s %8s %8s\n", prefix, "Atom", "Ng", "cuv", "clr", "chuv", "zetasr_m", "zetasr_a", "zetalr_m", "zetalr_a");
        for (int iv=0; iv<nv; iv++){
            for (int io=0; io<2; io++) if (fout[io]) if (nv>1) fprintf(fout[io], "%s%7s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", prefix, sys->av[iv].name, site_energy[iv].dNg, site_energy[iv].cuv, site_energy[iv].clr, site_energy[iv].chuv, site_energy[iv].zeta[0], site_energy[iv].zeta[1], site_energy[iv].zeta[2], site_energy[iv].zeta[3]);
        }
        for (int io=0; io<2; io++) if (fout[io]) fprintf(fout[io], "%s%7s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", prefix, "total", total_energy.dNg, total_energy.cuv, total_energy.clr, total_energy.chuv, total_energy.zeta[0], total_energy.zeta[1], total_energy.zeta[2], total_energy.zeta[3]);
    }
};
