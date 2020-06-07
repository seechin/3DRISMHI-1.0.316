




void fill_in_report_mass(IET_Param * sys, IET_Report & total, IET_Report * sites){
    double molecule_mass[MAX_SOL]; memset(molecule_mass, 0, sizeof(molecule_mass));
    for (int iv=0; iv<sys->nv; iv++) molecule_mass[sys->av[iv].iaa] += sys->av[iv].mass * sys->av[iv].multi;
    for (int iv=0; iv<sys->nv; iv++){
        sites[iv].mass = sys->av[iv].mass * sys->av[iv].multi;
        sites[iv].mass_mol = molecule_mass[sys->av[iv].iaa];
        sites[iv].density = sys->density_mv[sys->av[iv].iaa];
    }
    total.density = 0; for (int ivm=0; ivm<sys->nvm; ivm++) total.density += sys->density_mv[ivm];
    total.mass = 0; for (int iv=0; iv<sys->nv; iv++){ total.mass += sites[iv].mass * sites[iv].density / total.density; };
    total.mass_mol = total.mass;
}
double calculate_excessive_chemical_potential(int closure, double closure_factor, double ccutoff, double uuv, double huv, double hlr, double cuv, double clr){
    double excessive_GF = - (0.5*huv + 1) * (cuv + clr); double t = huv - cuv; double tz = -uuv + huv - cuv;
    if (closure<0) return excessive_GF;
    switch (closure) {
        case CLOSURE_HNC            : return excessive_GF + 0.5*huv*huv;
        case CLOSURE_PLHNC          : return excessive_GF + (huv<exp(ccutoff)? 0.5*huv*huv : 0);
        case CLOSURE_MHNC           :
        case CLOSURE_MSA            :
        case CLOSURE_KGK            :
        case CLOSURE_PY             :
        case CLOSURE_D2             : return excessive_GF;
        case CLOSURE_HNCB           : return excessive_GF + 0.5*huv*huv; // not finished
        case CLOSURE_KH             : return excessive_GF + (huv<0? 0.5*huv*huv : 0);
        case CLOSURE_PSE2           : return excessive_GF + 0.5*huv*huv - (tz>0? tz*tz*tz/6 : 0);
        case CLOSURE_PSE3           : return excessive_GF + 0.5*huv*huv - (tz>0? tz*tz*tz*tz/24 : 0);
        case CLOSURE_PSE4           : return excessive_GF + 0.5*huv*huv - (tz>0? tz*tz*tz*tz*tz/120 : 0);
        case CLOSURE_PSE5           : return excessive_GF + 0.5*huv*huv - (tz>0? tz*tz*tz*tz*tz*tz/720 : 0);
        case CLOSURE_PSE6           : return excessive_GF + 0.5*huv*huv - (tz>0? tz*tz*tz*tz*tz*tz*tz/5040 : 0);
        case CLOSURE_PSE7           : return excessive_GF + 0.5*huv*huv - (tz>0? tz*tz*tz*tz*tz*tz*tz*tz/40320 : 0);
        case CLOSURE_PSE8           : return excessive_GF + 0.5*huv*huv - (tz>0? tz*tz*tz*tz*tz*tz*tz*tz*tz/362880 : 0);
        case CLOSURE_PSE9           : return excessive_GF + 0.5*huv*huv - (tz>0? tz*tz*tz*tz*tz*tz*tz*tz*tz*tz/3628800 : 0);
        case CLOSURE_PSE10          : return excessive_GF + 0.5*huv*huv - (tz>0? tz*tz*tz*tz*tz*tz*tz*tz*tz*tz*tz/39916800 : 0);
        case CLOSURE_MS             :
        case CLOSURE_MSHNC          :
        case CLOSURE_BPGGHNC        :
        case CLOSURE_VM             :
        case CLOSURE_MP             : return excessive_GF;
    }
    return excessive_GF;
}

void generate_report_data(IET_Param * sys, IET_arrays * arr, IET_Report & total, IET_Report * sites, double dielect=1, double scale_lj=1, double scale_coul=1){
    size_t N3 = arr->nx * arr->ny * arr->nz;
    double dV = arr->box.x * arr->box.y * arr->box.z / N3;
    double beta = sys->default_temperature / sys->temperature;

  // step 1. fundamental data, e.g. mass
    memset(&total, 0, sizeof(IET_Report)); memset(&sites[0], 0, sys->nv*sizeof(IET_Report));
    fill_in_report_mass(sys, total, sites);

  // step 2. site dependent terms
    for (int iv=0; iv<sys->nv; iv++){
        int ivm = sys->av[iv].iaa; double q = sys->av[iv].charge_esp / dielect; double dN = dV * sys->density_av[iv];
        __REAL__ * uuv1 = &arr->uuv[iv][0][0][0];
        __REAL__ * huv1 = &arr->huv[iv][0][0][0];
        __REAL__ * hlr1 = &arr->hlr[iv][0][0][0];
        __REAL__ * dd1 = arr->dd? &arr->dd[sys->av[iv].iaa][0][0][0] : nullptr;
        __REAL__ * cuv1 = &arr->cuv[iv][0][0][0];
        __REAL__ * clr1 = &arr->clr[iv][0][0][0];
        double nbulk = sys->nbulk[sys->av[iv].iaa];
        sites[iv].N0 = arr->box.x * arr->box.y * arr->box.z * sys->density_av[iv];
        double factor = sys->closure_factors[iv]; double rhob = 0; int n_rhob = 0;
        for (size_t i3=0; i3<N3; i3++){
            double g = (1+huv1[i3]) * (dd1? dd1[i3]/nbulk : 1);
            double h = g - 1;
            double dn = g * dN * sys->av[iv].multi;
          // potential energy terms
            sites[iv].N += g * dN;
            if (dn>MACHINE_REASONABLE_ERROR) sites[iv].Ng += dN;
            sites[iv].lj += arr->ulj[iv][0][0][i3] * dn * scale_lj;
            sites[iv].coulsr += (arr->ucoulsr[0][0][i3]) * q * dn * scale_coul;
            sites[iv].coullr += (arr->ucoullr[0][0][i3]) * q * dn * scale_coul;
            sites[iv].entropy += dn<=0? 0 : - dn * log(g / sys->nbulk[ivm]) / beta;
          // Chandler's energy terms
            double excessive_GF = (- (0.5*huv1[i3] + 1) * (cuv1[i3]+clr1[i3])) * dN * sys->av[iv].multi * (dd1?dd1[i3]/nbulk : 1) / beta;
            double excessive_RISM = calculate_excessive_chemical_potential(sys->closures[iv], sys->closure_factors[iv], sys->ccutoff, uuv1[i3], huv1[i3], hlr1[i3], cuv1[i3], clr1[i3]) * (dd1?dd1[i3]:1) * dN * sys->av[iv].multi / beta;
            sites[iv].excess_chem[0] += excessive_GF;
            sites[iv].excess_chem[1] += excessive_RISM;
            double excessive_GF_LR = (- (0.5*huv1[i3] + 1) * clr1[i3]) * dN * sys->av[iv].multi * (dd1?dd1[i3]/nbulk : 1) / beta;
            double excessive_RISM_SR = calculate_excessive_chemical_potential(sys->closures[iv], sys->closure_factors[iv], sys->ccutoff, uuv1[i3], huv1[i3], hlr1[i3], cuv1[i3], 0) * (dd1?dd1[i3]:1) * dN * sys->av[iv].multi / beta;
            sites[iv].excess_chem[2] += excessive_RISM_SR + excessive_GF_LR;
          // correlations
            sites[iv].cuv += (cuv1[i3]) * dN * sys->av[iv].multi;
            sites[iv].clr += (clr1[i3]) * dN * sys->av[iv].multi;
            double zeta_this = calculate_zeta_by_chuv(sys->closures[iv], sys->closure_factors[iv], uuv1[i3], huv1[i3], cuv1[i3], hlr1[i3]);
            double zeta_hnc = calculate_zeta_by_chuv(CLOSURE_HNC, sys->closure_factors[iv], uuv1[i3], huv1[i3], cuv1[i3], hlr1[i3]);
            if (dn>sys->gcutoff_liquid_occupation){
                sites[iv].zeta[0] += zeta_hnc * dN * (dd1? dd1[i3]:nbulk) ;
                sites[iv].zeta[1] += zeta_hnc * dN * sys->av[iv].multi * (dd1? dd1[i3]:nbulk) ;
                sites[iv].zeta[2] += zeta_this * dN * (dd1? dd1[i3]:nbulk) ;
                sites[iv].zeta[3] += zeta_this * dN * sys->av[iv].multi * (dd1? dd1[i3]:nbulk) ;
            }
            if (dn>sys->gcutoff_liquid_occupation) sites[iv].chuv += (huv1[i3] - cuv1[i3])*dN;
        }
        sites[iv].dN = sites[iv].N - sites[iv].N0;
        sites[iv].dNg = sites[iv].Ng - sites[iv].N0;
        //sites[iv].pmv_from_cuv = sys->av[iv].multi * (dd1?dd1[i3]/nbulk : 1) *
    }
  // step 3. space dependent terms
    double total_Uef0_share = 0; for (int iv=0; iv<sys->nv; iv++) total_Uef0_share += fabs(sys->av[iv].charge_esp*sys->av[iv].multi) * sys->density_av[iv];
    if (arr->Ecoul0) for (size_t i3=0; i3<N3; i3++){
        Vector Eef0 = Vector(arr->Ecoul0[0][0][0][i3], arr->Ecoul0[1][0][0][i3], arr->Ecoul0[2][0][0][i3]);
        double Eef02 = Eef0.pow2();
        double this_Uef0 = dV * (1/(4*PI*COULCOOEF)) * Eef02 * sys->scale_coul;
        double this_Delta_Uef0 = this_Uef0 * (1 - sys->mean_dielect);
        for (int iv=0; iv<sys->nv; iv++){
            int ivm = sys->av[iv].iaa;
            double g = (1+arr->huv[iv][0][0][i3]) * (arr->dd? arr->dd[ivm][0][0][i3]/sys->nbulk[ivm] : 1);
            double share_of_Uef0 = total_Uef0_share==0? 1/sys->nv : fabs(sys->av[iv].charge_esp*sys->av[iv].multi) * sys->density_av[iv] / total_Uef0_share;
            if (g>sys->gcutoff_ef_occupation){
                sites[iv].Uef0 += this_Uef0 * share_of_Uef0;
                sites[iv].Uef1 += this_Delta_Uef0 * share_of_Uef0;
            }
        }
        //DeltaUef0 += this_Uef0 * (1 - 1/epsilon_r);
    }
  // step 4. add to sum
    for (int iv=0; iv<sys->nv; iv++) total += sites[iv];
    //printf("total: dN=%g LJ=%g cuv=%g Uef0=%g excessChem=[%g,%g,%g]\n", total.dN, total.lj, total.cuv, total.Uef0, total.excess_chem[0], total.excess_chem[1], total.excess_chem[2]);
}


void recalculate_energy(IET_Param * sys, IET_arrays * arr){
    generate_report_data(sys, arr, *arr->report_total, arr->report_sites, 1, sys->scale_lj, sys->scale_coul);
}
