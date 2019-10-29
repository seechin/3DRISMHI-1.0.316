//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Potential Force Field Building: build_force_field_sr and build_force_field_lr
//  potential = LJSR_v[z,y,x] + q_v CoulSR[z,y,x] + q_v PME[z,y,x]
//      LJSR_v[z,y,x] = arr->ulj[v][z][y][x]
//      CoulSR[z,y,x] = arr->ucoulsr[z][y][x]
//      PME[z,y,x]    = arr->ucoullr[z][y][x]
// Steps:
//  1. build_force_field_sr : build ucoulsr (kJ/mol), uljsr (kJ/mol/e), r2uvmin (nm²)
//  2. build_force_field_sr_connect : OPTIONAL, forbid solvents in unavailable regions
//  3. build_force_field_lr : build ucoullr (kJ/mol/e)
//  4. build_force_field_uuv* : build uuv:
//      build_force_field_uuv_ur
//          uuv = beta*(LJSR) (in kT), ulr = beta*Coul
//      build_force_field_uuv_Yukawa
//          uuv = beta*(LJSR) (in kT), ulr = beta*Yukawa
//      build_force_field_uuv_LPBE
//          uuv = beta*(LJSR) (in kT), ulr = beta*LPB
//      build_force_field_uuv_post_merge
//          uuv -> uuv+ulr, ulr -> 0
//-----------------------------------------------------------------------------
//#if defined(_LocalFFSRKDParallel_) && defined(_LOCALPARALLEL_)
//    #include "rismhi3d-build-ff-sr-kd.cpp"
//#endif

void build_force_field_sr_1(int i_begin, int i_end, int * i_now, IET_Param * sys, IET_arrays * arr, __REAL__ **** ulj, __REAL__ *** ucoulsr, __REAL__ *** r2uvmin, __REAL__ **** Ecoul0, __REAL__ *** pseudoliquid_potential, double gamma = 0){
    if (gamma<=0) gamma = sys->gamma_erf;
    PDBAtom * ia = sys->traj.atom; int nx = sys->nr[0]; int ny = sys->nr[1]; int nz = sys->nr[2]; int nv = sys->nv;
    Vector box = sys->traj.box; size_t N3 = nx*ny*nz;
    double drx = box.x / sys->nr[0]; double dry = box.y / sys->nr[1]; double drz = box.z / sys->nr[2];
    double coul_cooef = COULCOOEF;
    double dr_perturb = Vector(drx, dry, drz).mod() /1000;

  // 2. short range potential calculation
    double r2vdw = sys->rvdw * sys->rvdw; double r2coul = sys->rcoul * sys->rcoul;
    int gcx_max = (int)((sys->rvdw>sys->rcoul? sys->rvdw : sys->rcoul) / drx); if (gcx_max>nx/2) gcx_max = nx/2;
    int gcy_max = (int)((sys->rvdw>sys->rcoul? sys->rvdw : sys->rcoul) / dry); if (gcy_max>ny/2) gcy_max = ny/2;
    int gcz_max = (int)((sys->rvdw>sys->rcoul? sys->rvdw : sys->rcoul) / drz); if (gcz_max>nz/2) gcz_max = nz/2;

    double dbohr = sys->rbohr*2; double dbohr2 = dbohr*dbohr; double epsilon_bohr = sys->epsilon_bohr;
    //for (int ja=0; ja<sys->traj.count; ja++){
    double last_displayed_time = 0; bool is_out_tty = isatty(fileno(sys->log()));
    for (int ja=i_begin; ja<i_end && ja<sys->traj.count; ja++){
        //if (sys->debug_level>=5) fprintf(sys->log(), "DEBUG :: ffsr of %s.%s (mass=%g, charge=%g, sigma=%g, epsilon=%g)\n", sys->as[ja].mole, sys->as[ja].name, sys->as[ja].mass, sys->as[ja].charge, sys->as[ja].sigma, sys->as[ja].epsilon);
      // 2.1 prepare parameters for short range interaction
        int gx = (int)(ia[ja].r.x / drx); int gy = (int)(ia[ja].r.y / dry); int gz = (int)(ia[ja].r.z / drz);
        double sigma2 = sys->as[ja].sigma*sys->as[ja].sigma;
        double rbcharge = sys->rbcharge<=0? (sys->as[ja].sigma/2<dbohr?dbohr:sys->as[ja].sigma/2/(sys->rbcharge>-MACHINE_REASONABLE_ERROR?1:-sys->rbcharge)) : sys->rbcharge;
        double charge = sys->as[ja].charge;
        double charge_mb = sys->as[ja].charge; if (sys->as[ja].nbond>0){ charge_mb = sys->as[ja].charge * 0.5;
            for (int i=0; i<sys->as[ja].nbond; i++){ int ii = ja+sys->as[ja].ibond[i]; if (ii>=0&&ii<sys->nas) charge_mb += sys->as[ii].charge*(1-0.5) / (sys->as[ii].nbond<1?1:sys->as[ii].nbond); }
        } // fprintf(sys->log(), "ATOM %3d\33[31m%6s\33[0m:%3d\33[34m%6s\33[0m has charge %12f and %s%12f\33[0m\n", sys->as[ja].index, sys->as[ja].name, sys->as[ja].iaa, sys->as[ja].mole, charge, charge_mb>0.6?"\33[34m":charge_mb<-0.6?"\33[31m":charge_mb>0.3?"\33[32m":charge_mb<-0.3?"\33[35m":"\33[37m", charge_mb);
      // 2.1 short range interaction
        for (int _iz=-gcz_max; _iz<=gcz_max; _iz++) for (int _iy=-gcy_max; _iy<=gcy_max; _iy++) for (int _ix=-gcx_max; _ix<=gcx_max; _ix++){
            int ix = _ix + gx; if (sys->pbc_x) ix = ff_pbc_i(ix, nx); if (ix<0 || ix>=nx) continue;
            int iy = _iy + gy; if (sys->pbc_y) iy = ff_pbc_i(iy, ny); if (iy<0 || iy>=ny) continue;
            int iz = _iz + gz; if (sys->pbc_z) iz = ff_pbc_i(iz, nz); if (iz<0 || iz>=nz) continue;
            //printf("    %3d %3d %3d , %3d %3d %3d , %3d %3d %3d\n", _ix, _iy, _iz, ix, iy, iz, sys->pbc_x,sys->pbc_y,sys->pbc_z);

            double deltax = (_ix + gx) * drx - ia[ja].r.x;
            double deltay = (_iy + gy) * dry - ia[ja].r.y;
            double deltaz = (_iz + gz) * drz - ia[ja].r.z;

          // distance and minimal distance calculation
            double r2 = deltax*deltax + deltay*deltay + deltaz*deltaz; //double r = -1;
            if (r2 <= sys->errtolrism) r2 = sys->errtolrism;
            if (r2uvmin[iz][iy][ix]<0||r2<r2uvmin[iz][iy][ix]) r2uvmin[iz][iy][ix] = r2;
            //if (r2 <= r2vdw && sys->as[ja].is_key && (r2uvmin[iz][iy][ix]<0||r2<r2uvmin[iz][iy][ix])) r2uvmin[iz][iy][ix] = r2;

          // Lennard-Jones short range (LJSR)
            if (r2<dbohr2){
                for (int iv=0; iv<nv; iv++) ulj[iv][iz][iy][ix] += epsilon_bohr;
            } else if (r2<r2vdw){
              // LJSR
                for (int iv=0; iv<nv; iv++){
                    double sigma = sys->mixingrule_sigma_geometric? (sys->as[ja].sqrt_sigma*sys->av[iv].sqrt_sigma) : (sys->as[ja].sigma + sys->av[iv].sigma)/2;
                    double epsilon_2_j = sys->as[ja].sqrt_epsilon; double epsilon_2_i = sys->av[iv].sqrt_epsilon;
                    double ulj_current_point = 0;
                    if (epsilon_2_i<0 && epsilon_2_j<0){
                        //if (r2<sigma*sigma) ulj[iv][iz][iy][ix] += epsilon_bohr;
                        if (r2<sigma*sigma) ulj_current_point = epsilon_bohr;
                    } else if (epsilon_2_i<0){
                        sigma = sys->av[iv].sigma;
                        //if (r2<sigma*sigma) ulj[iv][iz][iy][ix] += epsilon_bohr;
                        if (r2<sigma*sigma) ulj_current_point = epsilon_bohr;
                    } else if (epsilon_2_j<0){
                        sigma = sys->as[ja].sigma;
                        //if (r2<sigma*sigma) ulj[iv][iz][iy][ix] += epsilon_bohr;
                        if (r2<sigma*sigma) ulj_current_point = epsilon_bohr;
                    } else {
                        double epsilon = 4 * epsilon_2_j * epsilon_2_i;
                        double s2 = sigma*sigma/r2; double s6 = s2*s2*s2; double s12 = s6*s6;
                        //ulj[iv][iz][iy][ix] += epsilon * (s12 - s6);
                        ulj_current_point = epsilon * (s12 - s6);
                    }
                    ulj[iv][iz][iy][ix] += ulj_current_point;
                }
            }

          // Coulomb potential short range (CoulSR) and Coulomb field short range (Ecoul0 short range part):
          // all scaled with smearing at rbcharge: -rq or sigma
          // for PME scaled with Erf
            double r = -1;
            if (r2<r2coul && charge!=0){
              // common
                if (r<0) r = sqrt(r2);
              // Coulomb SR: potential and electrostatic field
                double coulpre = coul_cooef * charge;
                if (sys->perform_pme){
                    double erf_factor = 1 - erf(gamma * r);
                    if (r>=rbcharge){
                        ucoulsr[iz][iy][ix] += coulpre * erf_factor / r;
                        double r_2 = r + dr_perturb; double r2_2 = r_2 * r_2;
                        double erf_factor_2 = 1 - erf(gamma * r_2);
                        double fval = coulpre * (erf_factor_2 / r_2 - erf_factor / r) / dr_perturb;
                        Ecoul0[0][iz][iy][ix] += fval * (deltax / r);
                        Ecoul0[1][iz][iy][ix] += fval * (deltay / r);
                        Ecoul0[2][iz][iy][ix] += fval * (deltaz / r);
                    } else {
                        //ucoulsr[iz][iy][ix] += coulpre * erf_factor / rbcharge;
                        ///*
                        double x = r / rbcharge;
                        ucoulsr[iz][iy][ix] += coulpre * erf_factor * (2 - x*x) / rbcharge;
                        double r_2 = r + dr_perturb; double r2_2 = r_2 * r_2; double x_2 = r_2 / rbcharge;
                        double erf_factor_2 = 1 - erf(gamma * r_2);
                        double fval = coulpre * (erf_factor_2 * (2 - x_2*x_2) - erf_factor * (2 - x*x)) / rbcharge;
                        Ecoul0[0][iz][iy][ix] += fval * (deltax / r);
                        Ecoul0[1][iz][iy][ix] += fval * (deltay / r);
                        Ecoul0[2][iz][iy][ix] += fval * (deltaz / r);
                        //*/
                    }
                } else {
                    if (r>=rbcharge){
                        ucoulsr[iz][iy][ix] += coulpre / r;
                        double r_2 = r + dr_perturb; double r2_2 = r_2 * r_2;
                        double fval = coulpre * (1 / r_2 - 1 / r) / dr_perturb;
                        Ecoul0[0][iz][iy][ix] += fval * (deltax / r);
                        Ecoul0[1][iz][iy][ix] += fval * (deltay / r);
                        Ecoul0[2][iz][iy][ix] += fval * (deltaz / r);
                    } else {
                        //ucoulsr[iz][iy][ix] += coulpre / rbcharge;
                        ///*
                        double x = r / rbcharge;
                        ucoulsr[iz][iy][ix] += coulpre * (2 - x*x) / rbcharge;
                        double r_2 = r + dr_perturb; double r2_2 = r_2 * r_2; double x_2 = r_2 / rbcharge;
                        double fval = coulpre * (1 * (2 - x_2*x_2) - 1 * (2 - x*x)) / rbcharge;
                        Ecoul0[0][iz][iy][ix] += fval * (deltax / r);
                        Ecoul0[1][iz][iy][ix] += fval * (deltay / r);
                        Ecoul0[2][iz][iy][ix] += fval * (deltaz / r);
                        //*/
                    }
                }
            }

          // pseudoliquid model for HI
          if (r2<sigma2/sys->pseudoliquid_potential_r_factor && pseudoliquid_potential){
              pseudoliquid_potential[iz][iy][ix] += charge;
            //if (r2<sigma2/sys->pseudoliquid_potential_r_factor && pseudoliquid_potential){
            //    double s2 = sigma2; if (s2<0.3*0.3) s2 = 0.3*0.3;
            //    pseudoliquid_potential[iz][iy][ix] += COULCOOEF * charge / (s2*1.25992105);
            //if (r2<sigma2/1.25992105/1.25992105 && pseudoliquid_potential){
            //    pseudoliquid_potential[iz][iy][ix] += COULCOOEF * charge / (sigma2/1.25992105/1.25992105);
            //if (r2<pseudoliquid_potential_r2 && pseudoliquid_potential){
            //    pseudoliquid_potential[iz][iy][ix] += charge;
            }

          // others

        }
      // loop end
        if (sys->detail_level>=3 && i_now){
            *i_now = ja+1-i_begin; // i_now points to ffsr_mp.irange[it][2] or nullptr
            int total_atom_calculated = 0; bool is_main_replica = true;
            double time_now = get_current_time_double()/1000;
          #ifdef _LOCALPARALLEL_
            if (arr->ffsr_mp.np>1) for (int it=0; it<arr->ffsr_mp.np; it++) total_atom_calculated += arr->ffsr_mp.irange[it][2];
            if (i_now != &arr->ffsr_mp.irange[0][2]) is_main_replica = false;
          #else
            total_atom_calculated = ja+1;
          #endif
            if (is_main_replica && (time_now-last_displayed_time>=10 || ja+1==i_end)){ last_displayed_time = time_now;
                fprintf(sys->log(), "  FFSR: %d atoms (%.0f%%)%s", total_atom_calculated, total_atom_calculated*100.0/sys->nas, (is_out_tty&&sys->detail_level==1)?"\r":"\n");
                fflush(sys->log());
            }
        }
    }
}

#ifdef _LOCALPARALLEL_
void merge_force_field_mp_data(IET_arrays * arr, size_t i3b, size_t i3e, size_t i4b, size_t i4e){
    for (int it=0; it<arr->ffsr_mp.np; it++){
        for (size_t i4=i4b; i4<i4e; i4++) arr->ulj[0][0][0][i4] += arr->ffsr_mp.lj[it][0][0][0][i4];
        for (size_t i3=i3b; i3<i3e; i3++) arr->ucoulsr[0][0][i3] += arr->ffsr_mp.coulsr[it][0][0][i3];
        //for (size_t i3=i3b; i3<i3e; i3++) arr->ulpbe[0][0][i3] += arr->ffsr_mp.ulpbe[it][0][0][i3];
        for (size_t i3=i3b; i3<i3e; i3++) arr->Ecoul0[0][0][0][i3] += arr->ffsr_mp.coulp2[it][0][0][0][i3];
        for (size_t i3=i3b; i3<i3e; i3++) arr->Ecoul0[1][0][0][i3] += arr->ffsr_mp.coulp2[it][1][0][0][i3];
        for (size_t i3=i3b; i3<i3e; i3++) arr->Ecoul0[2][0][0][i3] += arr->ffsr_mp.coulp2[it][2][0][0][i3];
        for (size_t i3=i3b; i3<i3e; i3++) arr->pseudoliquid_potential[0][0][i3] += arr->ffsr_mp.pseudoliquid_potential[it][0][0][i3];
        if (arr->r2uvmin) for (size_t i3=i3b; i3<i3e; i3++){
            if (arr->ffsr_mp.r2uvmin[it][0][0][i3]>0) if (arr->r2uvmin[0][0][i3]<0 || arr->r2uvmin[0][0][i3]>arr->ffsr_mp.r2uvmin[it][0][0][i3]){
                arr->r2uvmin[0][0][i3] = arr->ffsr_mp.r2uvmin[it][0][0][i3];
            }
        }
    }
}
void merge_force_field_mp_data(IET_Param * sys, IET_arrays * arr, int id){
    size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = arr->nv * N3;
    size_t i3b = N3 / sys->nt * id;
    size_t i3e = N3 / sys->nt * (id+1);
    size_t i4b = N4 / sys->nt * id;
    size_t i4e = N4 / sys->nt * (id+1);
    merge_force_field_mp_data(arr, i3b, i3e, i4b, i4e);
}
void merge_force_field_mp_data(IET_Param * sys, IET_arrays * arr){
  #ifdef _LOCALPARALLEL_
    for (int i=1; i<sys->nt; i++) __mp_tasks[i] = MPTASK_MERGE_FF_DATA;
    merge_force_field_mp_data(sys, arr, 0);
    wait_subroutines(sys);
  #else
    size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = arr->nv * N3;
    merge_force_field_mp_data(arr, 0, N3, 0, N4);
  #endif
}
#endif
void build_force_field_sr(IET_Param * sys, IET_arrays * arr){
    size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = sys->nv * N3;
  // 1. clear memory
    clear_tensor4d(arr->ulj, arr->nv*arr->nz*arr->ny*arr->nx);
    clear_tensor3d(arr->ucoulsr, arr->nz,arr->ny,arr->nx);
    //clear_tensor3d(arr->ulpbe, arr->nz,arr->ny,arr->nx);
    clear_tensor4d(arr->Ecoul0, 3*arr->nz*arr->ny*arr->nx);
    clear_tensor3d(arr->pseudoliquid_potential, arr->nz*arr->ny*arr->nx);
    for (size_t i3=0; i3<N3; i3++) arr->r2uvmin[0][0][i3] = -1;
  // 2. Calculate the short range forcefield
  #ifdef _LOCALPARALLEL_
    if (arr->ffsr_mp.np<=1){
        build_force_field_sr_1(0, sys->traj.count, nullptr, sys, arr, arr->ulj, arr->ucoulsr, arr->r2uvmin, arr->Ecoul0, arr->pseudoliquid_potential);
    } else {
        arr->ffsr_mp.reset_for_calculation(true, true);
        /*if (!sys->advanced_paralling_ff_batch){
            double last_time = get_current_time_double(); double begin_time = last_time;
            for (int ja=0; ja<sys->traj.count; ja++){
                int i_task_idle = 0; for (int i=1; i<arr->ffsr_mp.np; i++) if (arr->ffsr_mp.mp_tasks[i]==MPTASK_NONE){ i_task_idle = i; break; }
                if (i_task_idle>0){
                    arr->ffsr_mp.irange[i_task_idle][0] = ja; arr->ffsr_mp.irange[i_task_idle][1] = ja+1;
                    arr->ffsr_mp.mp_tasks[i_task_idle] = MPTASK_FFSR;
                } else {
                    usleep(100); ja--;
                    //int id = 0;
                    //arr->ffsr_mp.irange[id][0] = ja; arr->ffsr_mp.irange[id][1] = ja+1;
                    //build_force_field_sr_1(arr->ffsr_mp.irange[id][0], arr->ffsr_mp.irange[id][1], &arr->ffsr_mp.irange[id][2], sys, arr, arr->ffsr_mp.lj[id], arr->ffsr_mp.coulsr[id], arr->ffsr_mp.r2uvmin[id], arr->ffsr_mp.coulp2[id]);
                }
              // display information
                if (sys->flog && sys->detail_level>=1){
                    double current_time = get_current_time_double();
                    const char * termi_text = "\n"; if (ja+1<sys->traj.count && (isatty(fileno(sys->flog)) && sys->detail_level==1)) termi_text = "\r";
                    if (current_time-last_time>=1000 || (current_time-begin_time>=1000 && ja+1>=sys->traj.count)){ last_time = current_time;
                        fprintf(sys->flog, "  build-ff: atom[%d] out of %d%s", ja+1, sys->traj.count, termi_text);
                    }
                }
            }
        } else {*/
        for (int it=0; it<arr->ffsr_mp.np; it++){
            int i_task_idle = it; int ibegin = sys->traj.count/arr->ffsr_mp.np*it; int iend = sys->traj.count/arr->ffsr_mp.np*(it+1); if (it+1>=arr->ffsr_mp.np) iend = sys->traj.count;
            arr->ffsr_mp.irange[i_task_idle][0] = ibegin; arr->ffsr_mp.irange[i_task_idle][1] = iend;
        }
        for (int it=1; it<arr->ffsr_mp.np; it++) arr->ffsr_mp.mp_tasks[it] = MPTASK_FFSR;
        build_force_field_sr_1(arr->ffsr_mp.irange[0][0], arr->ffsr_mp.irange[0][1], &arr->ffsr_mp.irange[0][2], sys, arr, arr->ffsr_mp.lj[0], arr->ffsr_mp.coulsr[0], arr->ffsr_mp.r2uvmin[0], arr->ffsr_mp.coulp2[0], arr->ffsr_mp.pseudoliquid_potential[0]);
        //}


        double time0 = get_current_time_double(); double timeup_ms = -1;
        while (true){
            bool finished = true; for (int i=1; i<arr->ffsr_mp.np; i++) if (arr->ffsr_mp.mp_tasks[i]!=MPTASK_NONE){ finished = false; break; }
            if (finished) break;
            if (timeup_ms>0 && get_current_time_double()-time0 > timeup_ms) break;
            usleep(100);
        }

        merge_force_field_mp_data(sys, arr);
    }
  #else
    build_force_field_sr_1(0, sys->traj.count, nullptr, sys, arr, arr->ulj, arr->ucoulsr, arr->r2uvmin, arr->Ecoul0, arr->pseudoliquid_potential);
  #endif
    if (sys->detail_level>=3) fprintf(sys->log(), "  FFSR: all %d atoms processed\n", sys->nas);
    if (sys->debug_level>=3 || sys->debug_show_crc){
        if (arr->ulj    ) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(LJSR)       = %08X\n", check_real_crc(&arr->ulj[0][0][0][0], N4));
        if (arr->ucoulsr) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(CoulSR)     = %08X\n", check_real_crc(&arr->ucoulsr[0][0][0], N3));
        if (arr->ucoulsr) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(CoulLR)     = %08X\n", check_real_crc(&arr->ucoullr[0][0][0], N3));
        if (arr->Ecoul0 ) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(Elec_vac)   = %08X\n", check_real_crc(&arr->Ecoul0[0][0][0][0], N3*3));
        //if (arr->ulpbe  ) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(Yukawa)     = %08X\n", check_real_crc(&arr->ulpbe[0][0][0], N3));
        if (arr->r2uvmin) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(r2uvmin)    = %08X\n", check_real_crc(&arr->r2uvmin[0][0][0], N3));
    }
}


/*
void perform_PME(IET_Param * sys, IET_arrays * arr, double gamma){
    int nx = sys->nr[0]; int ny = sys->nr[1]; int nz = sys->nr[2];
    Vector box = sys->traj.box;
    double coul_cooef = COULCOOEF;
  // perform PME from charge distribution stored in arr->fftin
  // result stored to arr->fftout
    lap_timer_uuv();
    fftw_execute(arr->planf);
    lap_timer_fftw();
    double dkx = 2*PI / box.x; double dky = 2*PI / box.y; double dkz = 2*PI / box.z;
    double demo = 4 * gamma * gamma;//sys->gamma_erf * sys->gamma_erf;
    //int flips2[3] = { nx, ny, nz };

    double *** fftbuffer = arr->fftbuffer;
    for (int ix=0; ix<=nx/2; ix++) for (int iy=0; iy<=ny/2; iy++) for (int iz=0; iz<=nz/2; iz++){
        double k2 = vec_pow2((ix)*dkx, (iy)*dky, (iz)*dkz);
        double factor = k2==0? 0 : ((gamma<=0? 1 : exp(- k2 / demo)) / k2); // 0 at k=0: renormalization of electrostatic potentials
        fftbuffer[iz][iy][ix] = factor;
        fftbuffer[iz][iy][nx-ix] = factor;
        fftbuffer[iz][ny-iy][ix] = factor;
        fftbuffer[iz][ny-iy][nx-ix] = factor;
        fftbuffer[nz-iz][iy][ix] = factor;
        fftbuffer[nz-iz][iy][nx-ix] = factor;
        fftbuffer[nz-iz][ny-iy][ix] = factor;
        fftbuffer[nz-iz][ny-iy][nx-ix] = factor;
    }
    for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) for (int iz=0; iz<nz; iz++) arr->fftout[iz][iy][ix] *= arr->fftbuffer[iz][iy][ix] * 4*PI*coul_cooef /(box.x*box.y*box.z);
  // structure_factor * PME_Green_function -> fftin
    lap_timer_uuv();
    fftw_execute(arr->planb);
    lap_timer_fftw();
}
*/
void perform_PME(IET_Param * sys, IET_arrays * arr, double gamma){
    int nx = sys->nr[0]; int ny = sys->nr[1]; int nz = sys->nr[2];
    Vector box = sys->traj.box;
    double coul_cooef = COULCOOEF;
  // perform PME from charge distribution stored in arr->fftin
  // result stored to arr->fftout
    lap_timer_uuv();
    fftw_execute(arr->planf);
    lap_timer_fftw();
    double dkx = 2*PI / box.x; double dky = 2*PI / box.y; double dkz = 2*PI / box.z;
    double demo = 4 * gamma * gamma;//sys->gamma_erf * sys->gamma_erf;
    //int flips2[3] = { nx, ny, nz };

    __REAL__ *** fftout = arr->fftout;
    for (int ix=0; ix<=nx/2; ix++) for (int iy=0; iy<=ny/2; iy++) for (int iz=0; iz<=nz/2; iz++){
        double k2 = vec_pow2((ix)*dkx, (iy)*dky, (iz)*dkz);
        double factor = k2==0? 0 : ((gamma<=0? 1 : exp(- k2 / demo)) / k2); // 0 at k=0: renormalization of electrostatic potentials
        double intp_factor = factor * 4*PI*coul_cooef /(box.x*box.y*box.z);
        unsigned int mask = 0; if (ix==0||ix>=nx-ix) mask |= 1; if (iy==0||iy>=ny-iy) mask |= 2; if (iz==0||iz>=nz-iz) mask |= 4;
        fftout[iz][iy][ix] *= intp_factor;
        if (!(mask&1)) fftout[iz][iy][nx-ix] *= intp_factor;
        if (!(mask&2)) fftout[iz][ny-iy][ix] *= intp_factor;
        if (!(mask&3)) fftout[iz][ny-iy][nx-ix] *= intp_factor;
        if (!(mask&4)) fftout[nz-iz][iy][ix] *= intp_factor;
        if (!(mask&5)) fftout[nz-iz][iy][nx-ix] *= intp_factor;
        if (!(mask&6)) fftout[nz-iz][ny-iy][ix] *= intp_factor;
        if (!(mask&7)) fftout[nz-iz][ny-iy][nx-ix] *= intp_factor;
    }
  // structure_factor * PME_Green_function -> fftin
    lap_timer_uuv();
    fftw_execute(arr->planb);
    lap_timer_fftw();
}


void build_charge_mesh(IET_Param * sys, IET_arrays * arr, __REAL__ *** o, Vector * v_translation=nullptr){
    PDBAtom * ia = sys->traj.atom; int nv = sys->nv;
    int nx = sys->nr[0]; int ny = sys->nr[1]; int nz = sys->nr[2];
    Vector box = sys->traj.box; // Vector box = sys->traj.box * dim_factor;
    double drx = box.x / sys->nr[0]; double dry = box.y / sys->nr[1]; double drz = box.z / sys->nr[2];
    double coul_cooef = COULCOOEF;
    clear_tensor3d(o, nz,ny,nx);
    for (int ja=0; ja<sys->traj.count; ja++){
        Vector vr = ia[ja].r; if (v_translation) vr += *v_translation;
      // 2.2 long range preparation
        int gx = (int)(vr.x / drx); int gy = (int)(vr.y / dry); int gz = (int)(vr.z / drz);
        if (sys->pbc_x) gx = ff_pbc_i(gx, nx); if (sys->pbc_y) gy = ff_pbc_i(gy, ny); if (sys->pbc_z) gz = ff_pbc_i(gz, nz);
        if (gx>=0 && gx<nx && gy>=0 && gy<ny && gz>=0 && gz<nz){
            int ig1[3], ig2[3]; double w1[3], w2[3]; bool in_box = true;
            ig1[0] = (int)floor(vr.x / drx); ig2[0] = ig1[0] + 1; w2[0] = vr.x / drx - ig1[0]; w1[0] = 1 - w2[0]; if (ig1[0]<0 || ig1[0]>=nx || ig2[0]<0 || ig2[0]>=nx) in_box = false;
            ig1[1] = (int)floor(vr.y / dry); ig2[1] = ig1[1] + 1; w2[1] = vr.y / dry - ig1[1]; w1[1] = 1 - w2[1]; if (ig1[1]<0 || ig1[1]>=ny || ig2[1]<0 || ig2[1]>=ny) in_box = false;
            ig1[2] = (int)floor(vr.z / drz); ig2[2] = ig1[2] + 1; w2[2] = vr.z / drz - ig1[2]; w1[2] = 1 - w2[2]; if (ig1[2]<0 || ig1[2]>=nz || ig2[2]<0 || ig2[2]>=nz) in_box = false;
            if (in_box){
                o[ig1[2]][ig1[1]][ig1[0]] += sys->as[ja].charge * w1[2]*w1[1]*w1[0];
                o[ig1[2]][ig1[1]][ig2[0]] += sys->as[ja].charge * w1[2]*w1[1]*w2[0];
                o[ig1[2]][ig2[1]][ig1[0]] += sys->as[ja].charge * w1[2]*w2[1]*w1[0];
                o[ig1[2]][ig2[1]][ig2[0]] += sys->as[ja].charge * w1[2]*w2[1]*w2[0];
                o[ig2[2]][ig1[1]][ig1[0]] += sys->as[ja].charge * w2[2]*w1[1]*w1[0];
                o[ig2[2]][ig1[1]][ig2[0]] += sys->as[ja].charge * w2[2]*w1[1]*w2[0];
                o[ig2[2]][ig2[1]][ig1[0]] += sys->as[ja].charge * w2[2]*w2[1]*w1[0];
                o[ig2[2]][ig2[1]][ig2[0]] += sys->as[ja].charge * w2[2]*w2[1]*w2[0];
            } else o[gz][gy][gx] += sys->as[ja].charge;
//printf(":::: add %12f charge to %d %d %d: %.2f %.2f %.2f\n", sys->as[ja].charge, gx, gy, gz, vr.x/box.x*nx, vr.y/box.y*ny, vr.z/box.z*nz);
        }
    }
}


void build_force_field_lr(IET_Param * sys, IET_arrays * arr, bool clear_data_first = true, double gamma = -1){
    if (gamma<=0) gamma = sys->gamma_erf;
    if (clear_data_first) clear_tensor3d(arr->ucoullr, arr->nz,arr->ny,arr->nx);
    int nx = sys->nr[0]; int ny = sys->nr[1]; int nz = sys->nr[2];
  // 1. mesh chages
    build_charge_mesh(sys, arr, arr->fftin);
  // 2. long range potential calculation: fft
    perform_PME(sys, arr, gamma);
  // 3. Final step of COULLR potential
    for (int iz=0; iz<nz; iz++) for (int iy=0; iy<ny; iy++) for (int ix=0; ix<nx; ix++) arr->ucoullr[iz][iy][ix] += arr->fftin[iz][iy][ix];

  // 4. do partial derivatives and get long range Ecoul0
    Vector shift; Vector shift_box = Vector(arr->box.x/arr->nx, arr->box.y/arr->ny, arr->box.z/arr->nz);
    shift = Vector(-shift_box.x /1000, 0, 0);
      build_charge_mesh(sys, arr, arr->fftin, &shift); perform_PME(sys, arr, gamma);
      for (int iz=0; iz<nz; iz++) for (int iy=0; iy<ny; iy++) for (int ix=0; ix<nx; ix++) arr->Ecoul0[0][iz][iy][ix] += (arr->fftin[iz][iy][ix] - arr->ucoullr[iz][iy][ix]) / (shift_box.x /1000);
    shift = Vector(0, -shift_box.y /1000, 0);
      build_charge_mesh(sys, arr, arr->fftin, &shift); perform_PME(sys, arr, gamma);
      for (int iz=0; iz<nz; iz++) for (int iy=0; iy<ny; iy++) for (int ix=0; ix<nx; ix++) arr->Ecoul0[1][iz][iy][ix] += (arr->fftin[iz][iy][ix] - arr->ucoullr[iz][iy][ix]) / (shift_box.y /1000);
    shift = Vector(0, 0, -shift_box.z /1000);
      build_charge_mesh(sys, arr, arr->fftin, &shift); perform_PME(sys, arr, gamma);
      for (int iz=0; iz<nz; iz++) for (int iy=0; iy<ny; iy++) for (int ix=0; ix<nx; ix++) arr->Ecoul0[2][iz][iy][ix] += (arr->fftin[iz][iy][ix] - arr->ucoullr[iz][iy][ix]) / (shift_box.z /1000);
}














void build_force_field_uuv_ur(IET_Param * sys, IET_arrays * arr, double dielect=1){
    size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = N3 * sys->nv;
    double beta = sys->default_temperature / sys->temperature;
    for (int iv=0; iv<sys->nv; iv++){
        double charge = sys->av[iv].charge_esp;
        //double dielect = (sys->esal==CoulAL_Coulomb||sys->esal==CoulAL_NONE)? 1 : sys->mean_dielect;
        //double dielect = 1; if (sys->esal==CoulAL_Dielect) dielect = sys->mean_dielect;
        for (size_t i3=0; i3<N3; i3++){
            arr->uuv[iv][0][0][i3] = beta * (arr->ulj[iv][0][0][i3]) * (arr->ulj[iv][0][0][i3]>sys->ccutoff?sys->scale_hs:sys->scale_lj);// / dielect_lj;
            double ucoul = (arr->ucoulsr[0][0][i3] + arr->ucoullr[0][0][i3]) * charge / dielect * sys->scale_coul;
            arr->ulr[iv][0][0][i3] = beta * ucoul;
        }
    }
    lap_timer_uuv();
}
void build_force_field_uuv_post_merge(IET_Param * sys, IET_arrays * arr){
    size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = N3 * sys->nv;
  // treatment for un-renormalized RISM scheme
    for (int iv=0; iv<sys->nv; iv++){
        for (size_t i3=0; i3<N3; i3++){
            arr->uuv[iv][0][0][i3] += arr->ulr[iv][0][0][i3];
            arr->ulr[iv][0][0][i3] = 0;
        }
    }
    lap_timer_uuv();
}

void build_force_field_uuv_YukawaFFT(IET_Param * sys, IET_arrays * arr){
    size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = N3 * sys->nv;
    double beta = sys->default_temperature / sys->temperature;
  // build Yukawa field
    for (size_t i3=0; i3<N3; i3++) arr->uuv[0][0][0][i3] = arr->ucoulsr[0][0][i3] + arr->ucoullr[0][0][i3];
    perform_3rx1k_convolution(&arr->fftw_mp, arr->uuv, arr->nx, arr->ny, arr->nz, arr->box, 1, 1, arr->yukawa_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, arr->res, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
    //for (int ix=0; ix<arr->nx; ix++) for (int iy=0; iy<arr->ny; iy++) for (int iz=0; iz<arr->nz; iz++) if (fabs(arr->uuv[0][iz][iy][ix])>2 || fabs(arr->res[0][iz][iy][ix])>2) printf("%3d %3d %3d %18g %18g\n", ix, iy, iz, arr->uuv[0][iz][iy][ix], arr->res[0][iz][iy][ix]);
  // build force field
    for (int iv=0; iv<sys->nv; iv++){
        double charge = sys->av[iv].charge_esp;
        double dielect = sys->dielect_yukawa;
        for (size_t i3=0; i3<N3; i3++){
            arr->uuv[iv][0][0][i3] = beta * (arr->ulj[iv][0][0][i3]) * (arr->ulj[iv][0][0][i3]>sys->ccutoff?sys->scale_hs:sys->scale_lj);// / dielect_lj;
            double ucoul = (arr->res[0][0][0][i3]) * charge / dielect * sys->scale_coul;
            arr->ulr[iv][0][0][i3] = beta * ucoul;
        }
    }
    lap_timer_uuv();
}
