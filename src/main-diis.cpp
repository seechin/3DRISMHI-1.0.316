void printf_matrix(FILE * o, MatrixNS::Matrix & m){
    for (int i=0; i<m.n; i++){
        for (int j=0; j<m.n; j++) fprintf(o, " %14.7e", m.a[i][j]); fprintf(o, "\n");
    }
}



#ifndef __DIISNS__
  #define __DIIS_REAL__ __REAL__
  #define __DIISNS__
  namespace DIISNS {
    class DIIS {
      public:
        int ndiis; int nthread;
      private:
        int nDIISMax; double delta;
        long int count;
        __DIIS_REAL__ ** xi;
        __DIIS_REAL__ ** ri;
        __DIIS_REAL__ * data;
        __DIIS_REAL__ * weight;
        double * A_public;
        MatrixNS::MatrixDataContainer mcontainer;
        MatrixNS::Matrix A, Ai, B, T;
      #ifdef _LOCALPARALLEL_
        double A_buffer[MAX_THREADS][MAX_DIIS+1];
      #endif
      public:
        void realloc_matricses(){
            mcontainer.setip(0); A.init(&mcontainer, ndiis+1); Ai.init(&mcontainer, ndiis+1); B.init(&mcontainer, ndiis+1); T.init(&mcontainer, ndiis+1);
        }
        static long int size_for_data(int _nDIISMax, long int _count){ return sizeof(__DIIS_REAL__) * (long int)_count * (2*_nDIISMax+1); }
        static long int size_for_ri(int _nDIISMax){ return sizeof(__DIIS_REAL__*) * (_nDIISMax+1); }
        static long int size_for_xi(int _nDIISMax){ return sizeof(__DIIS_REAL__*) * (_nDIISMax+1); }
        static long int size_for_weight(int _nDIISMax){ return sizeof(__DIIS_REAL__) * (_nDIISMax + 1); }
        static long int size_for_A_public(int _nDIISMax){ return sizeof(double) * (_nDIISMax + 1) * (_nDIISMax + 1); }
        static long int size_required(int _nDIISMax, long int _count){ return size_for_data(_nDIISMax, _count) + size_for_ri(_nDIISMax) + size_for_xi(_nDIISMax) + size_for_weight(_nDIISMax) + size_for_A_public(_nDIISMax); }
        void set_buffer(int diis_max_count, int data_count, char * data_buffer){
            nDIISMax = diis_max_count; delta = 1; count = data_count;
            data = (__DIIS_REAL__*) data_buffer;//memalloc(sizeof(__DIIS_REAL__) * count * (2*nDIISMax+1));
            ri = (__DIIS_REAL__**) (data_buffer+size_for_data(nDIISMax, count));// memalloc(sizeof(__DIIS_REAL__*) * nDIISMax);
            xi = (__DIIS_REAL__**) (data_buffer+size_for_data(nDIISMax, count)+size_for_ri(nDIISMax));//memalloc(sizeof(__DIIS_REAL__*) * (nDIISMax+1));
            for (int i=0; i<nDIISMax+1; i++) xi[i] = &data[i * count];
            for (int i=0; i<nDIISMax; i++) ri[i] = &data[i * count + (nDIISMax+1) * count];
            weight = (__DIIS_REAL__*) (data_buffer+size_for_data(nDIISMax, count)+size_for_ri(nDIISMax)+size_for_xi(nDIISMax));//memalloc(sizeof(__DIIS_REAL__) * (nDIISMax + 1));
            A_public = (double*) (data_buffer+size_for_data(nDIISMax, count)+size_for_ri(nDIISMax)+size_for_xi(nDIISMax)+size_for_weight(nDIISMax)); memset(A_public, 0, size_for_A_public(nDIISMax));
            realloc_matricses();

            ndiis = 0;
        }
        void init(int diis_max_count, int data_count, char * data_buffer){
            mcontainer.init(20000);
            set_buffer(diis_max_count, data_count, data_buffer);
        }
        void init(int diis_max_count, int data_count, int _nthread=1){
            nDIISMax = diis_max_count; delta = 1; count = data_count; nthread = _nthread;
            data = (__DIIS_REAL__*) memalloc(sizeof(__DIIS_REAL__) * count * (2*nDIISMax+1));
            ri = (__DIIS_REAL__**) memalloc(sizeof(__DIIS_REAL__*) * nDIISMax);
            xi = (__DIIS_REAL__**) memalloc(sizeof(__DIIS_REAL__*) * (nDIISMax+1));
            for (int i=0; i<nDIISMax+1; i++) xi[i] = &data[i * count];
            for (int i=0; i<nDIISMax; i++) ri[i] = &data[i * count + (nDIISMax+1) * count];
            weight = (__DIIS_REAL__*) memalloc(sizeof(__DIIS_REAL__) * (nDIISMax + 1));
            A_public = (double*) memalloc(sizeof(double) * (nDIISMax + 1) * (nDIISMax + 1)); memset(A_public, 0, sizeof(double) * (nDIISMax + 1) * (nDIISMax + 1));
            mcontainer.init(20000);
            ndiis = 0;
        }
        void dispose(){
            free(ri); free(xi); free(data); free(weight);
            mcontainer.dispose();
            nDIISMax = count = 0; ri = xi = nullptr; data = nullptr;
        }
      public:
        void insert(__REAL__ * xsrc, __REAL__ * rsrc){
          // 1. shift the queue
            __DIIS_REAL__ * tmp;
            tmp = xi[nDIISMax]; for (int i=nDIISMax; i>0; i--) xi[i] = xi[i-1]; xi[0] = tmp;
            tmp = ri[nDIISMax-1]; for (int i=nDIISMax-1; i>0; i--) ri[i] = ri[i-1]; ri[0] = tmp;
            for (int i=nDIISMax; i>0; i--) for (int j=nDIISMax; j>0; j--) A_public[i*(nDIISMax+1)+j] = A_public[(i-1)*(nDIISMax+1)+(j-1)];
            for (int i=0; i<nDIISMax; i++) A_public[i*(nDIISMax+1)+0] = A_public[i] = 0;
          // 2. insert this data
            if (ndiis <= 0) for (int i=0; i<count; i++) xi[1][i] = xsrc[i];
            for (int i=0; i<count; i++){
                xi[0][i] = xsrc[i]; ri[0][i] = rsrc[i];
            }
          // 3. move pointer
            ndiis ++; if (ndiis > nDIISMax) ndiis = nDIISMax;
        }
        __DIIS_REAL__ error(){ __DIIS_REAL__ err = 0; for (int i=0; i<count; i++) err += ri[0][i]*ri[0][i]; return sqrt(err / count); }
      public:
       #ifdef _LOCALPARALLEL_
        void calc_weights_sub(int id){
            int idjamin = count / nthread * id; int idjamax = count / nthread * (id+1); if (id+1 >= nthread) idjamax = count;
            //for (int k=idjamin; k<idjamax; k++) for (int i=0; i<ndiis; i++) A_public[i] += ri[i][k] * ri[0][k];
            for (int k=idjamin; k<idjamax; k++) for (int i=0; i<ndiis; i++) A_buffer[id][i] += ri[i][k] * ri[0][k];
        }
       #endif
        bool calc_weights(IET_Param * sys){
            realloc_matricses();
            for (int i=0; i<B.n; i++) for (int j=0; j<B.n; j++) if (i==ndiis) B.a[i][j] = -1; else B.a[i][j] = 0;
            for (int i=0; i<A.n; i++) for (int j=0; j<A.n; j++) if ((i==ndiis || j==ndiis) && i!=j) A.a[i][j] = -1; else A.a[i][j] = 0;
          // 1. calculate aij
           #ifdef _LOCALPARALLEL_
            if (nthread<=1){
                for (int i=0; i<A.n && i<ndiis; i++) for (int k=0; k<count; k++) A_public[i] += ri[i][k] * ri[0][k];
                for (int i=1; i<A.n; i++) A_public[i*(nDIISMax+1)] = A_public[i];
            } else {
                for (int i=0; i<nthread; i++) for (int j=0; j<=ndiis; j++) A_buffer[i][j] = 0;
                for (int i=1; i<nthread; i++){
                    sys->mp_tasks[i] = MPTASK_DIIS_WEIGHT;
                }
                calc_weights_sub(0); wait_subroutines(sys);
                for (int id=0; id<nthread; id++) for (int i=0; i<A.n&&i<ndiis; i++) A_public[i] += A_buffer[id][i];
                for (int i=1; i<A.n; i++) A_public[i*(nDIISMax+1)] = A_public[i];
            }
           #else
            for (int i=0; i<A.n && i<ndiis; i++) for (int k=0; k<count; k++) A_public[i] += ri[i][k] * ri[0][k];
            for (int i=1; i<A.n; i++) A_public[i*(nDIISMax+1)] = A_public[i];
           #endif
            for (int i=0; i<A.n && i<ndiis; i++) for (int j=0; j<A.n && j<ndiis; j++) A.a[i][j] = A_public[i*(nDIISMax+1)+j];
            if (!A.inverse(Ai, T)){
                B = 0; B.a[1][0] = 1; return false;
            }
            B = Ai * B;
            for (int i=0; i<ndiis; i++) weight[i] = B.a[i][0];
          // 2. return weights
            return true;
        }
      public:
       #ifdef _LOCALPARALLEL_
        void step_in_mp_parallel(int begin, int end){
            for (int i=begin; i<end; i++){
                xi[0][i] = 0;
                for (int j=0; j<ndiis; j++){
                    xi[0][i] += (xi[j+1][i] + ri[j][i] * delta) * weight[j];
                }
            }
        }
        void step_in_mp_parallel(int id){
            int idjamin = count / nthread * id; int idjamax = count / nthread * (id+1); if (id+1 >= nthread) idjamax = count;
            step_in_mp_parallel(idjamin, idjamax);
        }
        void step_in(IET_Param * sys){
            if (nthread<=1) {
                step_in_mp_parallel(0, count);
            } else {
                for (int i=1; i<nthread; i++) sys->mp_tasks[i] = MPTASK_DIIS_STEPIN;
                step_in_mp_parallel(0);
                wait_subroutines(sys);
            }
        }
       #else
        void step_in_mp_parallel(int id){ }
        void step_in(IET_Param * sys){
            for (int i=0; i<count; i++){
                xi[0][i] = 0;
                for (int j=0; j<ndiis; j++){
                    xi[0][i] += (xi[j+1][i] + ri[j][i] * delta) * weight[j];
                }
            }
        }
       #endif
      public:
        __DIIS_REAL__ advance(IET_Param * sys, __REAL__ * xisrc, __REAL__ * risrc, double _delta, bool update_xisrc = true){
            delta = _delta;
//double tn = get_current_time_double();
            insert(xisrc, risrc);
//printf("DIIS::insert    %12f\n", get_current_time_double()-tn); tn = get_current_time_double();
            __DIIS_REAL__ err = error();
//printf("DIIS::error     %12f\n", get_current_time_double()-tn); tn = get_current_time_double();
            if (!calc_weights(sys)) err = -1;
//printf("DIIS::weight    %12f\n", get_current_time_double()-tn); tn = get_current_time_double();
            step_in(sys);
//printf("DIIS::stepin    %12f\n", get_current_time_double()-tn); tn = get_current_time_double();
            if (update_xisrc) for (int i=0; i<count; i++) xisrc[i] = xi[0][i];
            return err;
        }
    };
  }
#endif
