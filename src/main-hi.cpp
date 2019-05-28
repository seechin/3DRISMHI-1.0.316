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

void calculate_potential_HI(__REAL__ **** potential, __REAL__ **** dd, __REAL__ *** zeta, __REAL__ **** guvm, __REAL__ **** phi, double llambda[MAX_SOL], __REAL__ **** buf0, __REAL__ **** buf1, IET_Param * sys, IET_arrays * arr){
    size_t N3 = arr->nx*arr->ny*arr->nz;
  // 1. ρ_i = n_i ρ^b_i g_i
    //for (int iv=0; iv<sys->nvm; iv++) for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++) buf0[iv][iz][iy][ix] = dd[iv][iz][iy][ix] * guvm[iv][iz][iy][ix] * sys->density_hi;
    for (int iv=0; iv<sys->nvm; iv++) for (size_t i3=0; i3<N3; i3++){
        buf0[iv][0][0][i3] = (dd[iv][0][0][i3] * guvm[iv][0][0][i3] + (arr->nphi?arr->nphi[iv][0][0][i3]:0)) * sys->density_hi;
    }
  // 2. ρ_j * ζ_ij term -> buf1
    lap_timer_livm();
    //perform_3rx1k_convolution(buf0, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->zeta, arr->dk_zeta, arr->n_zeta, buf1, arr->fftin, arr->fftout, arr->planf, arr->planb);
    //perform_3rx1k_convolution(buf0, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->zeta, arr->dk_zeta, sys->xvv_k_shift, arr->n_zeta, buf1, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
    if (sys->nt<=1){
        perform_3rx1k_convolution(buf0, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->zeta, arr->dk_zeta, sys->xvv_k_shift, arr->n_zeta, buf1, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
    } else {
        perform_3rx1k_convolution(&arr->fftw_mp, buf0, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->zeta, arr->dk_zeta, sys->xvv_k_shift, arr->n_zeta, buf1, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
    }

    lap_timer_fftw();
  // 3. ρ*ζ + φ - lnλ -> potential
    /*for (int iv=0; iv<sys->nvm; iv++) for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++){
        potential[iv][iz][iy][ix] = buf1[iv][iz][iy][ix] + ((phi&&phi[0])?phi[iv][iz][iy][ix]:0) - llambda[iv];
    }*/
    for (int iv=0; iv<sys->nvm; iv++) for (size_t i3=0; i3<N3; i3++){
        potential[iv][0][0][i3] = buf1[iv][0][0][i3] - llambda[iv] + ((phi&&phi[0])?phi[iv][0][0][i3]:0);
    }

}
void subroutine_hi_solver(IET_Param * sys, IET_arrays * arr, int id){
    double lt[MAX_SOL]; double lts = 0;
    __REAL__ **** potential = arr->ddpot_hi;
  #ifdef _LOCALPARALLEL_
    size_t N3 = arr->nx*arr->ny*arr->nz;
    int idjamin = N3 / sys->nt * id; int idjamax = N3 / sys->nt * (id+1); if (id+1 >= sys->nt) idjamax = N3;
    for (size_t i3=idjamin; i3<idjamax; i3++){
        //int iz = i3 / (arr->nx*arr->ny); int iy = (i3 / arr->nx) % arr->ny; int ix = i3 % arr->nx;
  #else
    //for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++){
    size_t N3 = arr->nx*arr->ny*arr->nz;
    for (size_t i3=0; i3<N3; i3++){
  #endif
      // solve the total density
        double max_potential = potential[0][0][0][i3]; for (int iv=1; iv<sys->nvm; iv++) if (potential[iv][0][0][i3]>max_potential) max_potential = potential[iv][0][0][i3];
        lts = 0; for (int iv=0; iv<sys->nvm; iv++){ lt[iv] = exp(-(potential[iv][0][0][i3]-max_potential)); lts += lt[iv]; }
        double rhoa = arr->solver_hi.getx_bs(log(lts)-max_potential, sys->errtolhi/10, false);
        if (rhoa>1) rhoa = 1;
        //double rhoa = arr->solver_hi.getx_bs(log(lts)-max_potential, sys->errtolhi/10, true);
      // get each density -> arr->res
        for (int iv=0; iv<sys->nvm; iv++) arr->res[iv][0][0][i3] = rhoa * lt[iv] / lts;
        arr->res[sys->nvm][0][0][i3] = rhoa;
    }
}

void prepare_hi_guv(IET_Param * sys, IET_arrays * arr, bool first_run=true){
    // step 1: prepare guvm
      if (sys->guvmal == GUVMAL_THETA || sys->guvmal == GUVMAL_THETACC || first_run){
          for (int iv=0; iv<sys->nvm; iv++) for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++){
              double U = atomised_molecule_uv_potential(sys, arr, iv, ix, iy ,iz);
              if (U >= sys->ucutoff_hs) arr->guvm[iv][iz][iy][ix] = 0; else arr->guvm[iv][iz][iy][ix] = 1;
          }

          /*if (sys->guvmal == GUVMAL_THETACC){
              int N3 = arr->nx*arr->ny*arr->nz; double rc = sys->rvdw; if (rc>sys->rcoul) rc = sys->rcoul; double rc2 = rc*rc;
              __REAL__ **** cache1 = arr->rismhi_cache[2]; __REAL__ **** cache2 = arr->rismhi_cache[3];

              perform_3rx1k_convolution(&arr->fftw_mp, arr->guvm, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->ld_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache1, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
              for (int iv=0; iv<sys->nvm; iv++){
                  double ld_bulk = 0; double n_ld_bulk = 0;
                  for (size_t i3=0; i3<N3; i3++) if (!(arr->r2uvmin[0][0][i3]<rc2 && arr->r2uvmin[0][0][i3]>=0)){
                      n_ld_bulk ++; ld_bulk += cache1[iv][0][0][i3];
                  }
                  if (n_ld_bulk>0) ld_bulk /= n_ld_bulk;
                  if (ld_bulk!=0) for (size_t i3=0; i3<N3; i3++){
                      cache1[iv][0][0][i3] /= ld_bulk;
                      if (cache1[iv][0][0][i3]<0.99) cache1[iv][0][0][i3] = 0; else cache1[iv][0][0][i3] = 1;
                  }
              }

              perform_3rx1k_convolution(&arr->fftw_mp, cache1, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->ld_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache2, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
              for (int iv=0; iv<sys->nvm; iv++){
                  double ld_bulk = 0; double n_ld_bulk = 0;
                  for (size_t i3=0; i3<N3; i3++) if (!(arr->r2uvmin[0][0][i3]<rc2 && arr->r2uvmin[0][0][i3]>=0)){
                      n_ld_bulk ++; ld_bulk += cache2[iv][0][0][i3];
                  }
                  if (n_ld_bulk>0) ld_bulk /= n_ld_bulk;
                  if (ld_bulk!=0) for (size_t i3=0; i3<N3; i3++){
                      cache2[iv][0][0][i3] /= ld_bulk;
                  }
              }

              for (int iv=0; iv<sys->nvm; iv++) for (size_t i3=0; i3<N3; i3++) arr->guvm[iv][0][0][i3] *= cache2[iv][0][0][i3];

  //FILE * file_out = nullptr; char fn[MAX_PATH]; strcpy(fn, "local_density_2"); append_save_data(&file_out, fn, sys->log(), "ld", "", arr->nx, arr->ny, arr->nz, arr->nvm, &cache2[0][0][0][0], 1, sys); if (file_out) fclose(file_out);

          }*/

      } else if (sys->guvmal == GUVMAL_GUV){
          for (int i=0; i<sys->nvm; i++){ int guvmi = sys->vmmapi[i][0];
              for (int j=sys->vmmapi[i][0]; j<sys->vmmapi[i][1]; j++) if (sys->av[sys->vmmap[j]].is_key) guvmi = sys->vmmap[j];
              size_t N3 = arr->nx * arr->ny * arr->nz;
              for (size_t i3=0; i3<N3; i3++) arr->guvm[i][0][0][i3] = arr->huv[guvmi][0][0][i3] + 1;
          }
      }
      /*if (sys->hial == HIAL_EHSHI){
          size_t N3 = arr->nx*arr->ny*arr->nz; int N4 = sys->nvm*N3;
          __REAL__ **** cache1 = arr->rismhi_cache[2]; __REAL__ **** cache2 = arr->rismhi_cache[3];
          if (sys->nt<=1){
              perform_3rx1k_convolution(arr->guvm, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->zeta, arr->dk_zeta, sys->xvv_k_shift, arr->n_zeta, cache1, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
          } else {
              perform_3rx1k_convolution(&arr->fftw_mp, arr->guvm, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->zeta, arr->dk_zeta, sys->xvv_k_shift, arr->n_zeta, cache1, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
          }
          for (int ivm=0; ivm<sys->nvm; ivm++){
              for (size_t i3=0; i3<N3; i3++){
                  cache1[ivm][0][0][i3] /= arr->zeta[ivm][ivm][0]; if (cache1[ivm][0][0][i3]<0) cache1[ivm][0][0][i3] = 0; if (cache1[ivm][0][0][i3]>1) cache1[ivm][0][0][i3] = 1;
                  //cache1[ivm][0][0][i3] *= arr->guvm[ivm][0][0][i3];
                  arr->guvm[ivm][0][0][i3] *= cache1[ivm][0][0][i3];
              }
          }
  //FILE * file_out = nullptr; char fn[MAX_PATH]; strcpy(fn, "local_density"); append_save_data(&file_out, fn, sys->log(), "ld", "", arr->nx, arr->ny, arr->nz, arr->nvm, &cache1[0][0][0][0], 1, sys); if (file_out) fclose(file_out);
  }*/
}

bool prepare_phi(IET_Param * sys, IET_arrays * arr, int * hi_param_indicator, double * hi_param, int n_hi_param, bool guess_init_n=true){
    FILE * flog = sys->log();
    bool is_out_tty = sys->is_log_tty;
    bool success = false;
  // step 2: preparation, majorly phi
    size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = sys->nv * N3; size_t N4m = sys->nvm * N3;
    if (sys->hial == HIAL_HSHI){
        if (arr->phi&&arr->phi[0]) clear_tensor4d(arr->phi, N4m);
        if (arr->nphi&&arr->nphi[0]) clear_tensor4d(arr->nphi, N4m);
    } else if (sys->hial == HIAL_EHSHI){
        /*for (int ivm=0; ivm<sys->nvm; ivm++){
            for (size_t i3=0; i3<N3; i3++){
                if (atomised_molecule_uv_potential(sys, arr, ivm, i3) > sys->ccutoff_hi){
                    //arr->dd[ivm][0][0][i3] = 0;
                    arr->phi[ivm][0][0][i3] = sys->ccutoff_hi;
                } else arr->phi[ivm][0][0][i3] = 0;
            }
        }
        if (arr->nphi&&arr->nphi[0]) clear_tensor4d(arr->nphi, N4m);*/
        //for (size_t i4=0; i4<N4m; i4++) arr->guvm[0][0][0][i4] *= arr->dd[0][0][0][i4];

        for (int ivm=0; ivm<sys->nvm; ivm++){
            for (size_t i3=0; i3<N3; i3++) if (arr->r2uvmin[0][0][i3]<=0) arr->dd[ivm][0][0][i3] = 1; else arr->dd[ivm][0][0][i3] = 0;
        }

        if (arr->phi&&arr->phi[0]) clear_tensor4d(arr->phi, N4m);
        if (arr->nphi&&arr->nphi[0]) clear_tensor4d(arr->nphi, N4m);
    } else if (sys->hial==HIAL_DPHI){
        if (arr->phi && arr->phi[0]){
            double beta = sys->default_temperature/sys->temperature;
            //double nphi_unit = fabs(arr->zeta[0][0][0] *sys->density_mv[0]);
            for (int ivm=0; ivm<sys->nvm; ivm++){
                double dipole = sys->dipole_mv[ivm];
                for (size_t i3=0; i3<N3; i3++){
                    double force = arr->res[0][0][0][i3] = Vector(arr->Ecoul0[0][0][0][i3], arr->Ecoul0[1][0][0][i3], arr->Ecoul0[2][0][0][i3]).mod();
                    double pE = - fabs(force * dipole);
                    arr->phi[ivm][0][0][i3] = beta * pE / sys->dielect_hi;
                }
            }
//compress_data("phi.tensor4d", "nphi", flog, sys->nv, arr->nx, arr->ny, arr->nz, &arr->nphi[0][0][0][0], sys->err_output_uv);
//FILE * fout = nullptr; append_save_data(&fout, "debug-out.ts4s", sys->flog, "nPhi", nullptr, arr->nx, arr->ny, arr->nz, arr->nv, &arr->nphi[0][0][0][0], 0.0, sys); if(fout) fclose(fout); fprintf(stderr, "arr->zeta[0][0][0] = %g\n", arr->zeta[0][0][0]); if (fout) fclose(fout);
        }
        if (arr->nphi&&arr->nphi[0]) clear_tensor4d(arr->nphi, N4m);
    } else if (sys->hial==HIAL_DNHI){
        if (arr->nphi && arr->nphi[0]){
            double beta = sys->default_temperature/sys->temperature;
            double nphi_unit = fabs(arr->zeta[0][0][0] *sys->density_mv[0]);
            for (int ivm=0; ivm<sys->nvm; ivm++){
                double dipole = sys->dipole_mv[ivm];
                for (size_t i3=0; i3<N3; i3++){
                    double force = arr->res[0][0][0][i3] = Vector(arr->Ecoul0[0][0][0][i3], arr->Ecoul0[1][0][0][i3], arr->Ecoul0[2][0][0][i3]).mod();
                    double pE = fabs(force * dipole);
                    arr->nphi[ivm][0][0][i3] = beta * pE / sys->dielect_hi / nphi_unit; if (arr->nphi[ivm][0][0][i3]>1) arr->nphi[ivm][0][0][i3] = 1;
                }
            }
//compress_data("phi.tensor4d", "nphi", flog, sys->nv, arr->nx, arr->ny, arr->nz, &arr->nphi[0][0][0][0], sys->err_output_uv);
//FILE * fout = nullptr; append_save_data(&fout, "debug-out.ts4s", sys->flog, "nPhi", nullptr, arr->nx, arr->ny, arr->nz, arr->nv, &arr->nphi[0][0][0][0], 0.0, sys); if(fout) fclose(fout); fprintf(stderr, "arr->zeta[0][0][0] = %g\n", arr->zeta[0][0][0]);
        }
        if (arr->phi&&arr->phi[0]) clear_tensor4d(arr->phi, N4m);
    } else if (sys->hial==HIAL_PLHI || sys->hial==HIAL_RES_PLHI){
        double beta = sys->default_temperature / sys->temperature;
        int ivm_select = 0; for (int ivm=1; ivm<=sys->nvm; ivm++) if (sys->nbulk[ivm]>sys->nbulk[ivm_select]) ivm_select = ivm;

        double r_mol_selected = pow(fabs(1/(sys->bulk_density_mv[ivm_select]*4*PI/3)), 1.0/3);
        double zeta0_selected = arr->zeta[ivm_select][0][0] / sys->nbulk[ivm_select] * sys->density_hi / 4;
        double dipole_selected = sys->dipole_mv[ivm_select];
        double dipole_self_interaction = COULCOOEF*dipole_selected*dipole_selected/r_mol_selected/r_mol_selected/r_mol_selected; if (dipole_self_interaction<fabs(zeta0_selected)) dipole_self_interaction = -fabs(zeta0_selected);

//printf("selected effective solvent: %d, dipole %g, zeta %g\n", ivm_select, dipole_selected, zeta0_selected);
        if (arr->nphi && arr->nphi[0] && dipole_self_interaction!=0){
            clear_tensor4d(arr->nphi, N4m);
            double hardsphere_cutoff = sys->ucutoff_hs;
            double r_pseudoliquid_potential = r_mol_selected;
            double r2_pseudoliquid_potential = r_pseudoliquid_potential*r_pseudoliquid_potential;
            for (size_t i3=0; i3<N3; i3++){
                if (sys->hial!=HIAL_RES_PLHI || atomised_molecule_uv_potential(sys, arr, ivm_select, i3) > hardsphere_cutoff){
                    double n = fabs(beta * COULCOOEF * arr->pseudoliquid_potential[0][0][i3] / r2_pseudoliquid_potential * dipole_selected / sys->dielect_hi / dipole_self_interaction); if (n>1) n = 1;
                    arr->nphi[ivm_select][0][0][i3] = n;
                }
            }
        }
        if (arr->phi&&arr->phi[0]) clear_tensor4d(arr->phi, N4m);
    } else if (sys->hial==HIAL_CUTOFF){  // cutoff
        double beta = sys->default_temperature/sys->temperature; double dd_initialized_value[MAX_SOL];
        for (int ivm=0; ivm<sys->nvm; ivm++){
            dd_initialized_value[ivm] = (ivm<n_hi_param&&hi_param_indicator&&hi_param_indicator[ivm]&&hi_param)? hi_param[ivm] : 0;
            for (size_t i3=0; i3<N3; i3++){
                if (atomised_molecule_uv_potential(sys, arr, ivm, i3) > sys->ccutoff_hi){
                    //arr->dd[ivm][0][0][i3] = 0;
                    arr->dd[ivm][0][0][i3] = dd_initialized_value[ivm];
                } else arr->dd[ivm][0][0][i3] = sys->nbulk[ivm];
            }
        }
        fprintf(flog, "  %s : done with ", HIAL_name[sys->hial]); for (int i=0; i<sys->nvm; i++) fprintf(flog, i==0?"%g":",%g", dd_initialized_value[i]); fprintf(flog, "\n");
        return true;
    } else if (sys->hial==HIAL_ICUTOFF){  // cutoff
        double beta = sys->default_temperature/sys->temperature; double dd_initialized_value[MAX_SOL];
        for (int ivm=0; ivm<sys->nvm; ivm++){
            dd_initialized_value[ivm] = (ivm<n_hi_param&&hi_param_indicator&&hi_param_indicator[ivm]&&hi_param)? hi_param[ivm] : 0;
            for (size_t i3=0; i3<N3; i3++){
                if (atomised_molecule_uv_potential(sys, arr, ivm, i3) > sys->ccutoff_hi){
                    arr->dd[ivm][0][0][i3] = sys->nbulk[ivm];
                } else {
                    //arr->dd[ivm][0][0][i3] = 0;
                    arr->dd[ivm][0][0][i3] = dd_initialized_value[ivm];
                }
            }
        }
        fprintf(flog, "  %s : done with ", HIAL_name[sys->hial]); for (int i=0; i<sys->nvm; i++) fprintf(flog, i==0?"%g":",%g", dd_initialized_value[i]); fprintf(flog, "\n");
        //fprintf(flog, "  %s : done\n", HIAL_name[sys->hial]);
        return true;
    }
  // step 3: prepare initial guess of n[]
    if (guess_init_n){
        /*for (int iv=0; iv<sys->nvm; iv++) for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++){
            if (atomised_molecule_uv_potential(sys, arr, iv, ix, iy ,iz) > sys->ccutoff_hi) arr->dd[iv][iz][iy][ix] = 0; else arr->dd[iv][iz][iy][ix] = sys->nbulk[iv];
        }*/
        for (int iv=0; iv<sys->nvm; iv++) for (size_t i3=0; i3<=N3; i3++) arr->dd[iv][0][0][i3] = sys->nbulk[iv];
    }
    //for (int iv=0; iv<sys->nvm; iv++){ for (int jv=0; jv<sys->nvm; jv++){ printf(" %11g", arr->zeta[iv][jv][0]); }; printf("\n"); }

    return false;
}

bool perform_hi(IET_Param * sys, IET_arrays * arr, int * hi_param_indicator, double * hi_param, int n_hi_param, bool guess_init_n=true, bool first_run=true){
    if (!arr->dd||!arr->dd[0]) return true;
    FILE * flog = sys->log();
    bool is_out_tty = sys->is_log_tty;
    bool success = false;
    size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = sys->nv * N3; size_t N4m = sys->nvm * N3;

  //printf("HI=%s(%f), guvm=%s(%f), phi=%s, lse_A,B=%f,%s%f\33[0m\n", HIAL_name[sys->hial], sys->hial_const, GUVMAL_name[sys->guvmal], sys->ucutoff_hs, PHIAL_name[sys->phial], sys->lse_a, sys->calc_b_automatically?"\33[0;37m":"", sys->lse_b);
  //printf("  lnλ:\t"); for (int i=0; i<sys->nvm; i++)printf(" %12f",sys->llambda[i]); printf("\n");
  //printf("  n3:\t"); for (int i=0; i<sys->nvm; i++)printf(" %12f",sys->nbulk[i]); printf("\n");
  //printf("  <ζ>i:\t"); double mu_all = 0; for (int i=0; i<sys->nvm; i++){ double mu=0; for (int j=0; j<sys->nvm; j++) mu += sys->nbulk[j] * sys->density_hi * arr->zeta[i][j][0]; mu_all += mu * sys->nbulk[i]; printf(" %12f", mu); }; printf(" -> %12f\n", mu_all);
  // step 1: prepare guvm
    prepare_hi_guv(sys, arr, first_run);
  // step 2: preparation, majorly phi
    if (prepare_phi(sys, arr, hi_param_indicator, hi_param, n_hi_param, guess_init_n)){
        return true;
    }
  // step 4 prepare lambda for each site
    double llambda[MAX_SOL];
    for (int iv=0; iv<sys->nvm; iv++) llambda[iv] = sys->llambda[iv];
    for (int iv=0; iv<sys->nvm; iv++){
        for (int j=0; j<sys->nvm; j++){
            llambda[iv] += sys->nbulk[j] * arr->zeta[iv][j][0] * sys->density_hi;
            //printf("    %12f : %12g %12g %12g %12g\n", llambda[iv], sys->nbulk[j], arr->zeta[iv][j][0], sys->density_hi, sys->nbulk[j] * arr->zeta[iv][j][0] * sys->density_hi);
        }
        llambda[iv] += arr->solver_hi.f(1) + log(sys->nbulk[iv]);
        //printf("extra llambda[%d] = %12f\n", iv, llambda[iv]);
    }
  // step 5. HI loop
    //size_t flog_pointer = 0;
    for (int istep=0; istep<sys->stepmax_hi; istep++){
        while (sys->suspend_calculation){ usleep(100); sys->is_suspend_calculation = true; continue; }
        sys->is_suspend_calculation = false;
        //if (flog && flog!=stderr && flog!=stdout) flog_pointer = ftell(flog);
        double lt[MAX_SOL]; //double lts = 0;
      // 5.1. calculate potential -> arr->ddpot_hi
//lap_timer_livm();
        calculate_potential_HI(arr->ddpot_hi, arr->dd, arr->zeta, arr->guvm, arr->phi, llambda, arr->res, arr->ddpot_hi, sys, arr); //for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++) printf("%12f  #\n", arr->ddpot_hi[0][iz][iy][ix]);
        // post handling of potential: connectivity, bulk, boundary condition, so on
      // 5.2. calculate new dimensionless density -> arr->res
      //      including total particle number check
      /*#ifdef _LOCALPARALLEL_
        int id = fork_here(sys->nt); size_t N3 = arr->nx*arr->ny*arr->nz;
//printf("\33[1;35m  fork start: %d, max: %d\33[0m\n", id, sys->nt);
        int idjamin = N3 / sys->nt * id; int idjamax = N3 / sys->nt * (id+1); if (id+1 >= sys->nt) idjamax = N3;
        for (size_t i3=idjamin; i3<idjamax; i3++){
            int iz = i3 / (arr->nx*arr->ny); int iy = (i3 / arr->nx) % arr->ny; int ix = i3 % arr->nx;
      #else*/
      #ifdef _LOCALPARALLEL_
        for (int i=1; i<sys->nt; i++) sys->mp_tasks[i] = MPTASK_HI_SOLVER;
        subroutine_hi_solver(sys, arr, 0);
        wait_subroutines(sys);
      #else
        subroutine_hi_solver(sys, arr, 0);
      #endif
      // 5.2.2 total particle number check -> lt[]
        /*lts = 0; for (int iv=0; iv<sys->nvm; iv++){
            lt[iv] = 0; for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++) lt[iv] += arr->res[iv][iz][iy][ix]; lts += lt[iv];
        } for (int iv=0; iv<sys->nvm; iv++) lt[iv] /= lts; lts /= (arr->nx * arr->ny * arr->nz);*/
      // 5.3. step forward and loop control
        lap_timer_livm();
        double err = 0;
        if (sys->ndiis_hi<=1){
            for (int iv=0; iv<sys->nvm; iv++) for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++){
                double res = (arr->res[iv][iz][iy][ix] - arr->dd[iv][iz][iy][ix]);
                err += res * res;
                arr->dd[iv][iz][iy][ix] += res * sys->delhi;
            }
            err = sqrt(fabs(err)/sys->nvm/arr->nz/arr->ny/arr->nx);
        } else {
            for (int iv=0; iv<sys->nvm; iv++) for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++) arr->res[iv][iz][iy][ix] -= arr->dd[iv][iz][iy][ix];
            arr->diis_current = &arr->diis_hi;
            err = arr->diis_hi.advance(sys, &arr->dd[0][0][0][0], &arr->res[0][0][0][0], sys->delhi, true);
        }
        bool stop_loop = false;
            if (err <= sys->errtolhi){ success = true; stop_loop = true; }
            if (err < 0 || err > 10){ success = false; stop_loop = true; }
            if (istep+1 >= sys->stepmax_hi) stop_loop = true;
        if (flog && (sys->detail_level>=1 || stop_loop)){
            fprintf(flog, fabs(err)<1e-4?"  %s step %d, stdev: %.1e":"  %s step %d, stdev: %.5f", HIAL_name[sys->hial], istep+1, err);
            if (sys->ndiis_hi>1) fprintf(flog, " (DIISx%d)", arr->diis_hi.ndiis);
            //fprintf(flog, " (<n>=%.2f:", lts); for (int iv=0; iv<sys->nvm; iv++) fprintf(flog, lt[iv]<1?" %4.2f%%":" %4.0f%%", lt[iv]*100);
            //fprintf(flog, ")%s", (sys->detail_level==1&&istep+1<sys->stepmax_hi&&err>sys->errtolhi&&is_out_tty)?"\r":"\n");
            if (sys->debug_level>=3||sys->detail_level>=3){ char time_buffer[40]; fprintf(flog, " (%s)", get_current_time_text(time_buffer)); }
            fprintf(flog, "%s", (sys->detail_level==1&&istep+1<sys->stepmax_hi&&err>sys->errtolhi&&is_out_tty)?"\r":"\n");
            //if ((sys->detail_level==1&&err>sys->errtolhi) && flog && flog!=stderr && flog!=stdout) flog_pointer = fseek(flog, flog_pointer, SEEK_SET);
            fflush(flog);
        }
        fflush(flog);
//if (err <= sys->errtolhi) compress_data("potential_hi.tensor4d", "potential_hi", flog, sys->nvm, arr->nx, arr->ny, arr->nz, &arr->ddpot_hi[0][0][0][0], sys->err_output_uv);
        if (stop_loop) break;
        lap_timer_diis();
//lap_timer_diis();
    }

  // step 6. post preprocessing
    int n_dd_factor_value = 0; double dd_factor_value[MAX_SOL]; for (int ivm=0; ivm<sys->nvm; ivm++){
        dd_factor_value[ivm] = (ivm<n_hi_param&&hi_param_indicator&&hi_param_indicator[ivm]&&hi_param)? hi_param[ivm] : 1;
        if (dd_factor_value[ivm]!=1){
            n_dd_factor_value ++;
            for (size_t i3=0; i3<N3; i3++){
                if (atomised_molecule_uv_potential(sys, arr, ivm, i3) > sys->ccutoff_hi){
                    arr->dd[ivm][0][0][i3] *= dd_factor_value[ivm];
                }
            }
        }
    }
    if (n_dd_factor_value>0){
        fprintf(flog, "  %s : HS region scaled with ", HIAL_name[sys->hial]); for (int i=0; i<sys->nvm; i++) fprintf(flog, i==0?"%g":",%g", dd_factor_value[i]); fprintf(flog, "\n");
    }

  // step 6. If not converged, return false
    return success;
}
