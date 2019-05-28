
int maximum_default_processors(){
  #ifdef _LOCALPARALLEL_
    return sysconf(_SC_NPROCESSORS_ONLN);
  #else
    return 1;
  #endif
}
const char * get_short_fn(const char * filename){
    size_t i = strlen(filename)-1;
    while (i>0 && filename[i] && filename[i]!='/') i--;
    if (filename[i]=='/') i++;
    return &filename[i>0? i : 0];
}

void list_sys_files(IET_Param * sys, FILE * o, char * lineprefix=(char*)""){
    //fprintf(flog, "-b -e -s -f -log\n");
    bool istty = isatty(fileno(o));
    if (istty) fprintf(o, "%s", prompt_comment_prefix);
    fprintf(o, "%s%s %s\n", lineprefix, software_name, software_version);
    fprintf(o, "%spwd          %s%s%s%s\n", lineprefix, istty?prompt_path_prefix:"", get_second_fn(szfn_path), istty?prompt_path_suffix:"", istty?"\33[37m":"");
    fprintf(o, "%s-param       %s%s%s%s\n", lineprefix, istty?prompt_path_prefix:"", get_second_fn(info_file_name), istty?prompt_path_suffix:"", istty?"\33[37m":"");
    fprintf(o, "%s-solute      %s%s%s%s\n", lineprefix, istty?prompt_path_prefix:"", get_second_fn(szfn_solute), istty?prompt_path_suffix:"", istty?"\33[37m":"");
    fprintf(o, "%s-traj        %s%s%s%s\n", lineprefix, istty?prompt_path_prefix:"", get_second_fn(szfn_xtc), istty?prompt_path_suffix:"", istty?"\33[37m":"");
    if (!istty && !(!szfn_log[0] || StringNS::string(szfn_log)=="con"||StringNS::string(szfn_log)=="screen"||StringNS::string(szfn_log)=="stdout"||StringNS::string(szfn_log)=="stderr")){
        fprintf(o, "%s-log         %s\n", lineprefix, szfn_log[0]?get_second_fn(szfn_log):"con");
    } else {
        fprintf(o, "%s-log         %s%s%s%s\n", lineprefix, istty?prompt_path_prefix:"", szfn_log[0]?get_second_fn(szfn_log):"con", istty?prompt_path_suffix:"", istty?"\33[37m":"");
    }
    fprintf(o, "%s-be          %.0f,%.0f\n", lineprefix, sys->time_begin, sys->time_end);
    fprintf(o, "%s-dt          %g\n", lineprefix, sys->time_step);
    fprintf(o, "%s-nice        %d\n", lineprefix, sys->nice_level);
    if (istty) fprintf(o, "%s", prompt_comment_suffix);
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void main_generate_default_gvv_map(IET_Param * sys, IET_arrays * arr, FILE * flog){
    int icol = 0; int iargv = 0; char gvv_map_buffer[128]; memset(gvv_map_buffer, 0, sizeof(gvv_map_buffer));
    char * pgvv_map[3] = { &gvv_map_buffer[0], &gvv_map_buffer[30], &gvv_map_buffer[60] };

    for (int iv=0; iv<sys->nv; iv++){
        for (int jv=iv; jv<sys->nv; jv++){
            snprintf(pgvv_map[0], 28, "%d", iv+1);
            snprintf(pgvv_map[1], 28, "%d", jv+1);
            snprintf(pgvv_map[2], 28, "%d", icol+1);

            //fprintf(flog, "    %4s %4s %4s\n", pgvv_map[0], pgvv_map[1], pgvv_map[2]);
            analysis_gvv_map(sys, pgvv_map, &iargv, 3, (char*)"inner", icol+1);
            icol ++;
        }
    }
    if (sys->gvv_map){
        if (sys->mode_test) fprintf(flog, "%s : [gvv_map] set to default%s\n", software_name, sys->mode_test?":":".");
        if (sys->mode_test){ fprintf(flog, "    %4s", ""); for (int iv=0; iv<sys->nv; iv++) fprintf(flog, " %4s", sys->av[iv].name); fprintf(flog, "\n"); }
        for (int iv=0; iv<sys->nv; iv++){
            if (sys->mode_test) fprintf(flog, "    %4s", sys->av[iv].name);
            for (int jv=0; jv<sys->nv; jv++){
                int icol = -1;
                for (int i=0; i<sys->n_gvv_map; i++) if ((sys->gvv_map[i].grpi==iv+1 && sys->gvv_map[i].grpj==jv+1)||(sys->gvv_map[i].grpi==jv+1 && sys->gvv_map[i].grpj==iv+1)) icol = sys->gvv_map[i].col;
                if (sys->mode_test) fprintf(flog, " %4d", icol);
            }
            if (sys->mode_test) fprintf(flog, "\n");
        }
    }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void perform_cavity_removal_correction_uuv(IET_Param * sys, IET_arrays * arr){
    int N3 = arr->nx*arr->ny*arr->nz; double rc = sys->rvdw; if (rc>sys->rcoul) rc = sys->rcoul; double rc2 = rc*rc;
    __REAL__ **** cache0 = arr->res; __REAL__ **** cache1 = arr->rismhi_cache[2]; __REAL__ **** cache2 = arr->rismhi_cache[3];
  // 1. local density from potential cutoff -> cache0
    for (int iv=0; iv<sys->nvm; iv++) for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++){
        double U = atomised_molecule_uv_potential(sys, arr, iv, ix, iy ,iz);
        if (U >= sys->cavity_ucutoff) cache0[iv][iz][iy][ix] = 0; else cache0[iv][iz][iy][ix] = 1;
    }
  // 2. extract regions that fully accommodates at least one solvent molecule -> cache1
    perform_3rx1k_convolution(&arr->fftw_mp, cache0, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->ld_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache1, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
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
  // 3. calculate local density -> cache2
    perform_3rx1k_convolution(&arr->fftw_mp, cache1, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->ld_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache2, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
    double ld_forbidden_factor = sys->cavity_removal_factor;
    for (int iv=0; iv<sys->nvm; iv++){
      // local density rescale to have bulk as 1
        double ld_bulk = 0; double n_ld_bulk = 0;
        for (size_t i3=0; i3<N3; i3++) if (!(arr->r2uvmin[0][0][i3]<rc2 && arr->r2uvmin[0][0][i3]>=0)){
            n_ld_bulk ++; ld_bulk += cache2[iv][0][0][i3];
        }
        if (n_ld_bulk>0) ld_bulk /= n_ld_bulk;
        if (ld_bulk!=0) for (size_t i3=0; i3<N3; i3++){
            cache2[iv][0][0][i3] /= ld_bulk;
        }
      // find local density at forbidden regions
        double ld_forbidden = 0; double n_ld_forbidden = 0;
        for (int iz=0; iz<arr->nz; iz++) for (int iy=0; iy<arr->ny; iy++) for (int ix=0; ix<arr->nx; ix++){
            double U = atomised_molecule_uv_potential(sys, arr, iv, ix, iy ,iz);
            if (U >= sys->cavity_ucutoff){ n_ld_forbidden ++; ld_forbidden += cache2[iv][iz][iy][ix]; }
        }
        if (n_ld_forbidden>0){ ld_forbidden /= n_ld_forbidden; n_ld_forbidden = 0; }
        //printf("solvent %d has ld_forbidden %g\n", iv, ld_forbidden);
      // cutoff at local density of forbidden regions
        if (ld_bulk!=0) for (size_t i3=0; i3<N3; i3++){
            if (cache2[iv][0][0][i3]>ld_forbidden*ld_forbidden_factor) cache2[iv][0][0][i3] = 1; else cache2[iv][0][0][i3] = 0;
        }
    }
  // 4. fill unforgiven regions
    clear_tensor4d(arr->res, N3*sys->nv);
    for (int iv=0; iv<sys->nv; iv++){ int ivm = sys->av[iv].iaa;
        for (size_t i3=0; i3<N3; i3++) if (arr->uuv[iv][0][0][i3]<sys->cavity_ucutoff && cache2[ivm][0][0][i3]<0.5){
            arr->uuv[iv][0][0][i3] = sys->cavity_ucutoff;
            arr->res[iv][0][0][i3] = 1;
        }
    }
  // 5. finally
    //FILE * file_out = nullptr; char fn[MAX_PATH]; strcpy(fn, "local_density"); append_save_data(&file_out, fn, sys->log(), "ld", "", arr->nx, arr->ny, arr->nz, arr->nvm, &arr->res[0][0][0][0], 1, sys); if (file_out) fclose(file_out);
    //FILE * file_out = nullptr; char fn[MAX_PATH]; strcpy(fn, "local_density"); append_save_data(&file_out, fn, sys->log(), "ld", "", arr->nx, arr->ny, arr->nz, arr->nvm, &cache1[0][0][0][0], 1, sys); append_save_data(&file_out, fn, sys->log(), "ld", "", arr->nx, arr->ny, arr->nz, arr->nvm, &cache2[0][0][0][0], 1, sys); append_save_data(&file_out, fn, sys->log(), "ld", "", arr->nx, arr->ny, arr->nz, arr->nvm, &cache0[0][0][0][0], 1, sys); if (file_out) fclose(file_out);
}
void build_force_field_auto(IET_Param * sys, IET_arrays * arr, FILE * flog, int nframe){
    if (arr->frame_stamp != nframe) arr->uuv_is_ready = false;
    if (!arr->uuv_is_ready){
     // short range, as well as short range part of coulomb_p2
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_sr()\n");
        build_force_field_sr(sys, arr);
      // long range
        if (sys->perform_pme){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_lr()\n");
            build_force_field_lr(sys, arr, true);
        } else clear_tensor3d(arr->ucoullr, arr->nx * arr->ny * arr->nz);
      // external field
        if (sys->do_external_electrostatic_field){
            int nx = sys->nr[0]; int ny = sys->nr[1]; int nz = sys->nr[2];
            for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) for (int iz=0; iz<nz; iz++){
                double delta_x = arr->box.x/nx*(ix-nx/2);
                double delta_y = arr->box.y/ny*(iy-ny/2);
                double delta_z = arr->box.z/nz*(iz-nz/2);
                arr->ucoullr[iz][iy][ix] += delta_x*sys->external_electrostatic_field.x + delta_y*sys->external_electrostatic_field.y + delta_z*sys->external_electrostatic_field.z;
            }
        }

      // done
        lap_timer_uuv();
        arr->frame_stamp = nframe;
        arr->uuv_is_ready = true;
    }
}
void build_uuv_base_on_force_field(IET_Param * sys, IET_arrays * arr){
    if (sys->esal==CoulAL_YukawaFFT){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_uuv_YukawaFFT(kappa=%g)\n", 1/sys->rc_yukawafft);
        build_force_field_uuv_YukawaFFT(sys, arr);
  #ifdef _EXPERIMENTAL_
    } else if (sys->esal==CoulAL_Dipole){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_uuv_Dipole()\n");
        build_force_field_uuv_Dipole(sys, arr);
    } else if (sys->esal==CoulAL_Local_Dielect){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_uuv_ur_local(dielect=%g)\n", sys->mean_dielect);
        build_force_field_uuv_ur_local(sys, arr, sys->mean_dielect);
    } else if (sys->esal==CoulAL_Local){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_uuv_ur()\n");
        build_force_field_uuv_ur_local(sys, arr);
  #endif
    } else if (sys->esal==CoulAL_Dielect){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_uuv_ur(dielect=%g)\n", sys->mean_dielect);
        build_force_field_uuv_ur(sys, arr, sys->mean_dielect);
    } else {
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_uuv_ur()\n");
        build_force_field_uuv_ur(sys, arr);
    }

    if (sys->cavity_removal_correction){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_cavity_removal_correction_uuv(sf=%g,crf=%g,cutoff=%g)\n", sys->cavity_size_factor, sys->cavity_removal_factor, sys->cavity_ucutoff);
        perform_cavity_removal_correction_uuv(sys, arr);
    }

    if (sys->ietal==IETAL_SSOZ){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_uuv_post_merge()\n");
        build_force_field_uuv_post_merge(sys, arr);
    }


    if (sys->debug_level>=3) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(uuv)        = %08X\n", check_real_crc(&arr->uuv[0][0][0][0], arr->nv*arr->nx*arr->ny*arr->nz));
    if (sys->debug_level>=3) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(ulr)        = %08X\n", check_real_crc(&arr->ulr[0][0][0][0], arr->nv*arr->nx*arr->ny*arr->nz));
    if (sys->debug_level>=3) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(hlr)        = %08X\n", check_real_crc(&arr->hlr[0][0][0][0], arr->nv*arr->nx*arr->ny*arr->nz));

    //if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_coulp2(E)\n");
    //build_force_field_coulp2(sys, arr, arr->ulpbe, nullptr, arr->Ecoul[0], arr->Ecoul[1], arr->Ecoul[2]);


}

void reset_rdf(IET_Param * sys, RDF_data * rdf, int n_rdf_bins){
    for (int i1=0; i1<n_rdf_bins; i1++) for (int ig=0; ig<sys->n_rdf_grps; ig++){
        rdf[ig].n[i1] = 0;
        rdf[ig].g[i1] = 0;
    }
}

void calculate_rdf(IET_Param * sys, IET_arrays * arr, RDF_data * rdf, double r_max_rdf, int n_rdf_bins){
    if (!rdf || sys->n_rdf_grps<=0) return;
    double dr = r_max_rdf / n_rdf_bins; int N1 = n_rdf_bins;
    size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = sys->nv * N3;
    if (sys->rdf_content==IETCMD_v_dd){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * dd1 = arr->dd? &arr->dd[sys->av[iv].iaa][0][0][0] : nullptr;
            double nbulk = sys->nbulk[sys->av[iv].iaa];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = (dd1? dd1[i3] : nbulk);
        }
    } else if (sys->rdf_content==-IETCMD_v_dd){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * dd1 = &arr->nphi[sys->av[iv].iaa][0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = (dd1? dd1[i3] : 0);
        }
    } else if (sys->rdf_content==IETCMD_v_ucoul){
        for (int iv=0; iv<sys->nv; iv++){
            double charge = sys->av[iv].charge;
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * src1 = &arr->ucoulsr[0][0][0];
            __REAL__ * src2 = &arr->ucoullr[0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = (src1[i3] + src2[i3]) * charge * sys->scale_coul;
        }
    } else if (sys->rdf_content==IETCMD_v_Ef){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * src1 = &arr->Ecoul0[0][0][0][0];
            __REAL__ * src2 = &arr->Ecoul0[1][0][0][0];
            __REAL__ * src3 = &arr->Ecoul0[2][0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = (0.5/(4*PI*COULCOOEF)) * Vector(src1[i3], src2[i3], src3[i3]).pow2() * sys->scale_coul;
        }
    } else if (sys->rdf_content==IETCMD_v_cuv){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * cuv1 = &arr->cuv[iv][0][0][0];
            __REAL__ * ulr1 = &arr->ulr[iv][0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = cuv1[i3] - ulr1[i3];
        }
    } else if (sys->rdf_content==-IETCMD_v_cuv){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * ulr1 = &arr->ulr[iv][0][0][0];
            __REAL__ * hlr1 = &arr->hlr[iv][0][0][0];
            __REAL__ * huv1 = &arr->huv[iv][0][0][0];
            __REAL__ * cuv1 = &arr->cuv[iv][0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = huv1[i3] - cuv1[i3] + hlr1[i3] + ulr1[i3];
        }
    } else if (sys->rdf_content==IETCMD_v_ulj){
        for (size_t i4=0; i4<N4; i4++) arr->res[0][0][0][i4] = arr->ulj[0][0][0][i4]>sys->ccutoff? 0 : arr->ulj[0][0][0][i4];
    } else if (sys->rdf_content==IETCMD_v_Mayer){
        double beta = sys->default_temperature / sys->temperature;
        for (size_t i4=0; i4<N4; i4++) arr->res[0][0][0][i4] = exp(-beta*arr->ulj[0][0][0][i4]);
    } else if (sys->rdf_content==IETCMD_v_uuv || sys->rdf_content==IETCMD_v_huv || sys->rdf_content==-IETCMD_v_huv){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * src1 = nullptr;
            double scaling = 1;
            if (sys->rdf_content==IETCMD_v_ulj){        src1 = &arr->ulj[iv][0][0][0];    scaling = sys->scale_lj;
            } else if (sys->rdf_content==IETCMD_v_uuv){ src1 = &arr->uuv[iv][0][0][0];    scaling = sys->scale_coul;
            } else if (sys->rdf_content==IETCMD_v_huv){ src1 = &arr->huv[iv][0][0][0];
            } else if (sys->rdf_content==-IETCMD_v_huv){ src1 = &arr->hlr[iv][0][0][0];
            } else {                                    src1 = &arr->huv[iv][0][0][0];
            }
            for (size_t i3=0; i3<N3; i3++) res1[i3] = src1[i3] * scaling;
        }
    } else { // guv
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * huv1 = &arr->huv[iv][0][0][0];
            __REAL__ * dd1 = arr->dd? &arr->dd[sys->av[iv].iaa][0][0][0] : nullptr;
            double nbulk = sys->nbulk[sys->av[iv].iaa];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = (1+huv1[i3]) * (dd1? dd1[i3] : nbulk) / nbulk;
        }
    }

    for (int ig=0; ig<sys->n_rdf_grps; ig++){
        memset(rdf[ig].g, 0, sizeof(double)*N1); memset(rdf[ig].n, 0, sizeof(double)*N1);
        int is = rdf[ig].is; int iv = rdf[ig].iv;
          if (iv<0 || iv>=sys->nv) continue;
          if (is>=sys->nas) continue;
        if (is<0){
            for (int iz=0; iz<sys->nr[2]; iz++) for (int iy=0; iy<sys->nr[1]; iy++) for (int ix=0; ix<sys->nr[0]; ix++){
                double r = sqrt(fabs(arr->r2uvmin[iz][iy][ix]));
                int ithis = (int)(floor(r/dr)); if (ithis>=0 && ithis<N1){
                    rdf[ig].g[ithis] += arr->res[iv][iz][iy][ix];
                    rdf[ig].n[ithis] ++;
                }
            }
        } else {
            Vector center = sys->traj.atom[is].r; //printf("calculate rdf for pair: %d - %d, ref_center: %g %g %g\n", rdf[ig].is, rdf[ig].iv, center.x, center.y, center.z);
            for (int iz=0; iz<sys->nr[2]; iz++) for (int iy=0; iy<sys->nr[1]; iy++) for (int ix=0; ix<sys->nr[0]; ix++){
                //double r = (center - Vector(ix*arr->dx, iy*arr->dy, iz*arr->dz)).mod();

                //double cx = center.x; if (sys->pbc_x) cx = ff_pbc_d(cx, arr->box.x);
                //double cy = center.y; if (sys->pbc_y) cy = ff_pbc_d(cy, arr->box.y);
                //double cz = center.z; if (sys->pbc_z) cz = ff_pbc_d(cz, arr->box.z);
                //double r = (Vector(cx, cy, cz) - Vector(ix*arr->dx, iy*arr->dy, iz*arr->dz)).mod();

                Vector dv = Vector(center - Vector(ix*arr->dx, iy*arr->dy, iz*arr->dz));
                if (sys->pbc_x) dv.x = ff_pbc_d(dv.x+arr->box.x/2, arr->box.x)-arr->box.x/2;
                if (sys->pbc_y) dv.y = ff_pbc_d(dv.y+arr->box.y/2, arr->box.y)-arr->box.y/2;
                if (sys->pbc_z) dv.z = ff_pbc_d(dv.z+arr->box.z/2, arr->box.z)-arr->box.z/2;
                double r = dv.mod();

                int ithis = (int)(floor(r/dr)); if (ithis>=0 && ithis<N1){
                    rdf[ig].g[ithis] += arr->res[iv][iz][iy][ix];
                    rdf[ig].n[ithis] ++;
                }
            }
        }
    }
}
void add_rdf_to_sum(IET_Param * sys, IET_arrays * arr, RDF_data * rdf, RDF_data * rdfs, int n_rdf_bins){
    int N1 = n_rdf_bins;
    for (int ig=0; ig<sys->n_rdf_grps; ig++){
        for (int i=0; i<N1; i++){
            rdfs[ig].g[i] += rdf[ig].g[i];
            rdfs[ig].n[i] += rdf[ig].n[i];
        }
    }
}

bool print_rdf(FILE * file, IET_Param * sys, IET_arrays * arr, RDF_data * rdf, double r_max_rdf, int n_rdf_bins){
    if (!rdf || sys->n_rdf_grps<=0) return false;
    if (!file) return false;
    double dr = r_max_rdf / n_rdf_bins; int N1 = n_rdf_bins;
    unsigned long mask64 = generate_double_sig_dig_mask(sys->output_significant_digits);
      if (sys->output_significant_digits==-4) mask64 = generate_double_sig_dig_mask(7);
      if (sys->output_significant_digits==-8) mask64 = generate_double_sig_dig_mask(15);

    fprintf(file, "%8s", "r"); for (int ig=0; ig<sys->n_rdf_grps; ig++){
        char buffer[32];
        if (rdf[ig].is<0){
            snprintf(buffer, sizeof(buffer), "all-%d%s", rdf[ig].iv+1, sys->av[rdf[ig].iv].name);
        } else {
            snprintf(buffer, sizeof(buffer), "%d%s-%d%s", rdf[ig].is+1, sys->as[rdf[ig].is].name, rdf[ig].iv+1, sys->av[rdf[ig].iv].nele);
        }
        fprintf(file, " %8s", buffer);
    }; fprintf(file, "\n");
    for (int i1=0; i1<N1; i1++){
        fprintf(file, "%8.4f", (i1+0.5)*dr);
        for (int ig=0; ig<sys->n_rdf_grps; ig++){
            double value = rdf[ig].n[i1]>0? rdf[ig].g[i1]/rdf[ig].n[i1] : rdf[ig].g[i1];
            //bit_and(&value, &value, &mask64, 8);
            fprintf(file, " %8.4f", value);
        }
        fprintf(file, "\n");
    }

    return true;
}

bool print_rdf(const char * filename, IET_Param * sys, IET_arrays * arr, RDF_data * rdf, double r_max_rdf, int n_rdf_bins){
    bool fp_by_fopen = false;
    if (!rdf || sys->n_rdf_grps<=0) return false;
    if (!filename) return false;
    double dr = r_max_rdf / n_rdf_bins; int N1 = n_rdf_bins;
    FILE * file = stdout; StringNS::string fn = filename;
    if (fn=="stdout" || fn=="screen" || fn=="con"){ file = stdout;
    } else if (fn=="stderr"){   file = stderr;
    } else {    file = fopen(filename, "w"); fp_by_fopen = true;
    }
    if (!file){ fprintf(sys->log(), "%s%s : error : cannot write RDF to %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, filename, sys->is_log_tty?color_string_end:""); return false; }
    bool ret = print_rdf(file, sys, arr, rdf, r_max_rdf, n_rdf_bins);
    if (fp_by_fopen){
        size_t rdf_file_size = ftell(file);
        fclose(file);

        char save_items_mem_buffer[64];
        fprintf(sys->log(), "  save %s to %s%s%s, totally %s\n", "rdf", sys->is_log_tty?prompt_path_prefix:"\"", filename, sys->is_log_tty?prompt_path_suffix:"\"", print_memory_value(save_items_mem_buffer, sizeof(save_items_mem_buffer), rdf_file_size));
    }
    return ret;
    /*
    if (!rdf || sys->n_rdf_grps<=0) return false;
    if (!filename) return false;
    double dr = r_max_rdf / n_rdf_bins; int N1 = n_rdf_bins;
    FILE * file = stdout; StringNS::string fn = filename;
    if (fn=="stdout" || fn=="screen" || fn=="con"){ file = stdout;
    } else if (fn=="stderr"){   file = stderr;
    } else {    file = fopen(filename, "w");
    }
    if (!file){ fprintf(sys->log(), "%s : error : cannot write RDF to %s\n", software_name, filename); return false; }

    fprintf(file, "%7s", "r"); for (int ig=0; ig<sys->n_rdf_grps; ig++){
        char buffer[32];
        if (rdf[ig].is<0){
            snprintf(buffer, sizeof(buffer), "all-%d%s", rdf[ig].iv+1, sys->av[rdf[ig].iv].nele);
        } else {
            snprintf(buffer, sizeof(buffer), "%d%s-%d%s", rdf[ig].is+1, sys->as[rdf[ig].is].nele, rdf[ig].iv+1, sys->av[rdf[ig].iv].nele);
        }
        fprintf(file, " %7s", buffer);
    }; fprintf(file, "\n");
    for (int i1=0; i1<N1; i1++){
        fprintf(file, "%7.3f", (i1+0.5)*dr);
        for (int ig=0; ig<sys->n_rdf_grps; ig++) fprintf(file, " %7.3f", rdf[ig].n[i1]>0? rdf[ig].g[i1]/rdf[ig].n[i1] : rdf[ig].g[i1]);
        fprintf(file, "\n");
    }


    if (file!=stdout && file!=stderr) fclose(file);
    return true;
     */
}


void recalculate_energy(IET_Param * sys, IET_arrays * arr){
    arr->solvation_energy(sys, sys->ccutoff);
}

void apply_dielect_from_dipole(IET_Param * sys){
    sys->dielect_from_dipole = true;
    double beta = sys->default_temperature / sys->temperature;
    for (int ivm=0; ivm<sys->nvm; ivm++){
        sys->dielect_mol[ivm] = 1 + 4*PI/3*COULCOOEF * beta * sys->density_mv[ivm] * sys->dipole_mv[ivm]*sys->dipole_mv[ivm];
    }
    for (int i=0; i<sys->nav; i++){ int iaa = sys->av[i].iaa; if (iaa>=0 && iaa<MAX_CMD_PARAMS) sys->dielect[i] = sys->dielect_mol[iaa]; }
    double dielecti = 0; double density_dielect = 0; for (int ivm=0; ivm<sys->nvm; ivm++){ dielecti += sys->density_mv[ivm] / sys->dielect_mol[ivm]; density_dielect += sys->density_mv[ivm]; } sys->mean_dielect = density_dielect / dielecti;
    //if (sys->esal==CoulAL_Coulomb) sys->mean_dielect = 1;
}

void generate_default_output_filename(IET_Param * sys, char * filename, int size_filename, const char * traj_filename, int size_traj_filename){
    int ixtc = 0; int jxtc = 0;
    int len_traj_filename = strnlen(traj_filename, size_traj_filename);
    for (ixtc=len_traj_filename; ixtc>=0; ixtc--) if (traj_filename[ixtc]=='/') break;
    while (ixtc<size_traj_filename && traj_filename[ixtc]=='/') ixtc++;
    for (jxtc=len_traj_filename; jxtc>=0; jxtc--) if (traj_filename[jxtc]=='.') break;
    if (ixtc<0||ixtc>=len_traj_filename) ixtc = 0; if (jxtc<=ixtc) jxtc = len_traj_filename;
    char trajname[MAX_PATH]; memset(trajname, 0, sizeof(trajname)); memcpy(trajname, &traj_filename[ixtc], jxtc-ixtc);

    if (!filename[0]){
        time_t rawtime; time(&rawtime); struct tm * tmtime = localtime(&rawtime);
        char time_buffer[40]; snprintf(time_buffer, sizeof(time_buffer), "%02d%02d%02d_%02d%02d", tmtime->tm_year-100, tmtime->tm_mon+1, tmtime->tm_mday, tmtime->tm_hour, tmtime->tm_min);
        int nmsolute = 1; for (int i=1; i<sys->nas; i++) if (StringNS::string(sys->as[i].mole) != sys->as[i-1].mole) nmsolute ++;
        int nmsolvent = 1; for (int i=1; i<sys->nav; i++) if (StringNS::string(sys->av[i].mole) != sys->av[i-1].mole) nmsolvent ++;
        snprintf(filename, size_filename, "%s.%s%s%s.%s.ts4s", trajname, sys->av[0].mole, nmsolvent>1?"_":"", nmsolvent>2?"etc":nmsolvent==2?sys->av[sys->nav-1].mole:"", time_buffer);
    }
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

void main_print_header(IET_Param * sys, IET_arrays * arr, FILE * flog, int argc, char * argv[]){
    bool islogtty = isatty(fileno(flog));
        fprintf(flog, "============================================================================ \n");
        fprintf(flog, "%s %s\n", software_name, software_version);
        char time_buffer[40];
        fprintf(flog, "%s begins at %s\n%s@%s: %s%s%s\n", software_name, get_current_time_text(time_buffer), username, hostname, islogtty?prompt_path_prefix:"", szfn_path, islogtty?prompt_path_suffix:"");
        for (int i=0; i<argc; i++) fprintf(flog, "%s ", argv[i]); fprintf(flog, "\n");
        fprintf(flog, "---------------------------------------------------------------------------- \n");
}
void main_print_tailer(IET_Param * sys, IET_arrays * arr){
    FILE * flog = sys->log();
    char time_buffer[2][40];
    fprintf(flog, "%s ends at %s (%s s)\n", software_name, get_current_time_text(time_buffer[0]), display_time(__total_timer, time_buffer[1]));
    fprintf(flog, "============================================================================ \n");
}

bool main_prepair_rdf_grps_to_pairs(IET_Param * sys, IET_arrays * arr, RDF_data ** prdf, int * pn_rdf_pairs){
    bool success = true;

    RDF_data * rdf = nullptr;
    int n_rdf_pairs = 0;

    int npr = map_rdf_grps_to_pairs(sys, nullptr, sys->detail_level>=1); //printf("pair number %d, bins %d\n", npr, sys->out_rdf_bins);
    if (npr>0){
        rdf = (RDF_data*) memalloc( (sizeof(RDF_data)+2*(1+sizeof(double)*sys->out_rdf_bins)) * npr ); n_rdf_pairs = npr;
        if (rdf) for (int i=0; i<npr; i++){
            rdf[i].is = rdf[i].iv = 0; rdf[i].g = (double*) (((char *)rdf) + sizeof(RDF_data)*npr + 2*(1+sizeof(double)*sys->out_rdf_bins)*i); rdf[i].n = &rdf[i].g[sys->out_rdf_bins+1];
            for (int j=0; j<sys->out_rdf_bins; j++) rdf[i].g[j] = rdf[i].n[j] = 0;
        }
        int npr_mapped = map_rdf_grps_to_pairs(sys, rdf, false); //for (int i=0; i<n_rdf_pairs; i++) printf("rdf_pair[%d] = %d(%s.%s)-%d(%s.%s)\n", i+1, rdf[i].is+1, sys->as[rdf[i].is].mole,sys->as[rdf[i].is].name, rdf[i].iv+1, sys->av[rdf[i].iv].mole, sys->av[rdf[i].iv].name);
        if (npr_mapped<0) success = false;
    } else success = false;
    if (!success){
        fprintf(sys->log(), "%s : please use the following groups to calculate RDF:\n", software_name);
        int max_display_items = 10;
        fprintf(sys->log(), "  solutes:"); for (int i=0; i<sys->nas && i<max_display_items; i++) fprintf(sys->log(), " %d/%s:%s", i+1, sys->as[i].mole, sys->as[i].name); fprintf(sys->log(), sys->nas>max_display_items?" ...\n":"\n");
        fprintf(sys->log(), "  solvents:"); for (int i=0; i<sys->nv && i<max_display_items; i++) fprintf(sys->log(), " %d/%s:%s", i+1, sys->av[i].mole, sys->av[i].name); fprintf(sys->log(), sys->nv>max_display_items?" ...\n":"\n");
    }

    *prdf = rdf;
    *pn_rdf_pairs = n_rdf_pairs;

    return success;
}
void main_print_solute_list(IET_Param * sys, IET_arrays * arr, FILE * flog, const char * prefix){
    int max_display_solute_num = 10;
    //if (sys->detail_level>=2 && sys->debug_level>=2) max_display_solute_num = 99999999;
    fprintf(flog, "%ssolutes: %d atoms in %s: ", prefix, sys->nas, get_short_fn(szfn_xtc));
      int iaa_count = 0;
      for (int i=0; i<sys->nas; i++) if (i==0||StringNS::string(sys->as[i].mole)!=sys->as[i-1].mole){
          iaa_count++;
          if (iaa_count<max_display_solute_num) fprintf(flog, i==0?"%s":", %s", sys->as[i].mole);
          else if (iaa_count==max_display_solute_num) fprintf(flog, " ...");
      };
    if (iaa_count>=max_display_solute_num)fprintf(flog, " (%d residues)", iaa_count); fprintf(flog, "\n");
}
void main_print_solvent_list(IET_Param * sys, IET_arrays * arr, FILE * flog, const char * prefix){
    fprintf(flog, "%ssolvents: %d %s in: ", prefix, sys->nav, sys->nav>1?"sites":"site");
    for (int i=0; i<sys->nav; i++) if (i==0||StringNS::string(sys->av[i].mole)!=sys->av[i-1].mole){
        fprintf(flog, i==0?"%s":", %s", sys->av[i].mole);
    }
    fprintf(flog, "\n");
    fprintf(flog, "%ssolvent densities: ", prefix);
    for (int i=0; i<sys->nav; i++) if (i==0||StringNS::string(sys->av[i].mole)!=sys->av[i-1].mole){
        int iaa = sys->av[i].iaa;
        fprintf(flog, iaa==0?"%g/%g" : ", %g/%g", sys->density_mv[iaa], sys->bulk_density_mv[iaa]);
    }
    fprintf(flog, "\n");
}
void main_print_forcefield_info(IET_Param * sys, IET_arrays * arr, FILE * flog, const char * prefix){
    const char * ffstring = "(undefined)";
    if (!sys->mixingrule_sigma_geometric && fabs(sys->default_temperature- 502.97)<=MACHINE_REASONABLE_ERROR) ffstring = "amber";
    if (!sys->mixingrule_sigma_geometric && fabs(sys->default_temperature- 120.27)<=MACHINE_REASONABLE_ERROR) ffstring = "gaff";
    if (sys->mixingrule_sigma_geometric && fabs(sys->default_temperature- 120.27)<=MACHINE_REASONABLE_ERROR) ffstring = "oplsaa";
    fprintf(flog, "%sforcefield %s, rvdw %g, rcoul %g, temperature %g\n", prefix, ffstring, sys->rvdw, sys->rcoul, sys->temperature);
}

void print_save_extra_info(IET_Param * sys, char * buffer, size_t size){
    char time_buffer[40];
    //const char * short_info = get_short_fn(get_second_fn(info_file_name));
    //const char * short_solu = get_short_fn(get_second_fn(szfn_solute));
    //const char * short_traj = get_short_fn(get_second_fn(szfn_xtc));
    //snprintf(buffer, size, "%s %s %s -p %s -s %s -f %s -pwd %s", software_name, software_version, get_current_time_text(time_buffer), short_info, short_solu, short_traj, szfn_path);
    snprintf(buffer, size, "%s %s %s", software_name, software_version, get_current_time_text(time_buffer));
}

size_t select_append_save_data(int * filter, int nfilter, __REAL__ **** temp, size_t * total_size, size_t * original_size, FILE ** pfile, char filename[MAX_PATH], FILE * flog, const char * title, const char * text, int nx, int ny, int nz, int nv, __REAL__ **** data, double time_stamp, IET_Param * sys){

    bool save_all = false; bool save_none = true; int N3 = nx * ny * nz;
    if (!filter || nfilter<=0){
        save_all = true; save_none = false;
    } else{
        for (int i=0; i<nfilter; i++){
            if (filter[i]<0){ save_all = true; save_none = false; }
            else if (filter[i]>0) save_none = false;
        }
    }
    if (save_none){
//fprintf(stderr, "\33[31msave: none\n\33[0m");
        *total_size = *original_size = 0;
    } else if (save_all){
//fprintf(stderr, "\33[31msave: all\n\33[0m");
        *total_size += append_save_data(pfile, filename, flog, title, text, nx, ny, nz, nv, &data[0][0][0][0], time_stamp, sys);
        *original_size = nx*ny*nz*nv*sizeof(__REAL__) + sizeof(CompressPageHeader);
    } else {
        int n_save_sites = 0;
        for (int i=0; i<nfilter; i++) if (filter[i]>0 && filter[i]<=nv && n_save_sites<nv){
//fprintf(stderr, "\33[31msave: item[%d]=%d (of %d)\n\33[0m", i, filter[i]-1, nfilter);
            cp_tensor3d(data[filter[i]-1], temp[n_save_sites++], N3);

            //int N3 = nx*ny*nz;
            //__REAL__ * src = &data[n_save_sites][0][0][0]; __REAL__ * dst = &temp[filter[i]-1][0][0][0];
            //for (int j=0; j<N3; j++) dst[j] = src[j];
            //n_save_sites ++;
        }
        *total_size = append_save_data(pfile, filename, flog, title, text, nx, ny, nz, n_save_sites, &temp[0][0][0][0], time_stamp, sys);
        *original_size = nx*ny*nz*n_save_sites*sizeof(__REAL__) + sizeof(CompressPageHeader);
    }
    return *total_size;
}


void perform_sub_test(IET_Param * sys, IET_arrays * arr, IET_command * cmd){
    if (!cmd) return;
    PDBAtom * ia = sys->traj.atom; SoluteAtomSite * as = sys->as;
    size_t N3 = arr->nx*arr->ny*arr->nz; size_t N4 = N3 * sys->nv;
    __REAL__ **** cache1 = arr->rismhi_cache[2]; __REAL__ **** cache2 = arr->rismhi_cache[3];

    for (int iparam=0; iparam<cmd->step; iparam++){
        if (cmd->command_params_int[iparam]==IETCMD_v_ucoul || cmd->command_params_int[iparam]==IETCMD_v_hlr || cmd->command_params_int[iparam]==IETCMD_v_Yukawa){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_sub_test(%s%s)\n", cmd->command_params_int[iparam]==IETCMD_v_ucoul?"coulomb":cmd->command_params_int[iparam]==IETCMD_v_hlr?"hlr":cmd->command_params_int[iparam]==IETCMD_v_Yukawa?"Yukawa":"unknown", cmd->command==IETCMD_TEST_SAVE?", save":"");

            if (cmd->command_params_int[iparam]==IETCMD_v_Yukawa){
                build_force_field_uuv_YukawaFFT(sys, arr);
            } else {
                build_uuv_base_on_force_field(sys, arr);
            }
            RISMHI3D_RISMNS::prepare_3d_rism(sys, arr);
            /*
            for (size_t i3=0; i3<N3; i3++) cache1[0][0][0][i3] = 1;
            perform_3rx1k_convolution(&arr->fftw_mp, cache1, arr->nx, arr->ny, arr->nz, arr->box, 1, 1, arr->local_coulomb_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache2, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
            double convolution_mean = 0; for (size_t i3=0; i3<N3; i3++) convolution_mean += cache2[0][0][0][i3]; convolution_mean /= N3;
            if (sys->debug_level>=2) fprintf(sys->log(), "  1*Heaviside(%g-r) = %g\n", sys->rlocal_coul?sys->rcoul:sys->rlocal_coul, convolution_mean);

            for (size_t i3=0; i3<N3; i3++) cache1[0][0][0][i3] = arr->ucoulsr[0][0][i3] + arr->ucoullr[0][0][i3];
            perform_3rx1k_convolution(&arr->fftw_mp, cache1, arr->nx, arr->ny, arr->nz, arr->box, 1, 1, arr->local_coulomb_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache2, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
            for (size_t i3=0; i3<N3; i3++) cache2[0][0][0][i3] = cache1[0][0][0][i3] - cache2[0][0][0][i3]/convolution_mean;
            */

            if (cmd->command==IETCMD_TEST){
                for (int is=0; is<sys->nas; is++){
                    int ix = (sys->pbc_x? ff_pbc_d(ia[is].r.x, arr->box.x) : ia[is].r.x) / arr->box.x * arr->nx;
                    int iy = (sys->pbc_y? ff_pbc_d(ia[is].r.y, arr->box.y) : ia[is].r.y) / arr->box.y * arr->ny;
                    int iz = (sys->pbc_z? ff_pbc_d(ia[is].r.z, arr->box.z) : ia[is].r.z) / arr->box.z * arr->nz;
                    if (ix>=0&&ix<arr->nx && iy>=0&&iy<arr->ny && iz>=0&&iz<arr->nz){
                        if (cmd->command_params_int[iparam]==IETCMD_v_ucoul || cmd->command_params_int[iparam]==IETCMD_v_Yukawa){
                            fprintf(sys->log(), "#%6d %6s %6s q=%8.3g (%4d %4d %4d) coul=%12.2f\n", is+1, as[is].mole, as[is].name, as[is].charge, ix, iy, iz, arr->ucoulsr[iz][iy][ix]+arr->ucoullr[iz][iy][ix]);
                        } else {
                            fprintf(sys->log(), "#%6d %6s %6s q=%8.3g (%4d %4d %4d) coul=%12.2f hlr=%12.2f\n", is+1, as[is].mole, as[is].name, as[is].charge, ix, iy, iz, arr->ucoulsr[iz][iy][ix]+arr->ucoullr[iz][iy][ix], arr->hlr[0][iz][iy][ix]);
                        }

                    }
                }
            } else if (cmd->command==IETCMD_TEST_SAVE){
                if (cmd->command_params_int[iparam]==IETCMD_v_ucoul){
                  // save sys->(ucoulsr + ucoullr)
fprintf(sys->log(), "\33[31msave Coulomb\n\33[0m");
                } else if (cmd->command_params_int[iparam]==IETCMD_v_Yukawa){
fprintf(sys->log(), "\33[31msave Yukawa\n\33[0m");
                  // save sys->ulr
                } else if (cmd->command_params_int[iparam]==IETCMD_v_hlr){
fprintf(sys->log(), "\33[31msave hlr\n\33[0m");
                  // save sys->hlr
                }
            }
        } else if (cmd->command_params_int[iparam]==IETCMD_v_LocalCoulomb){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_sub_test(local-Coulomb)\n");

        } else {
            fprintf(sys->log(), " %s%s : error : perform_sub_test(%d): unrecognized command%s\n", sys->is_log_tty?color_string_of_error:"", software_name, cmd->command_params_int[iparam], sys->is_log_tty?color_string_end:"");
        }
    }

    return ;






}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
