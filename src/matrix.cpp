//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Simple Matrix Algorithm
// (c) Cao Siqin 2014.10.13 (revised) : save-restore ip for matrixd
// (c) Cao Siqin 2013.10.10 (revised) : LU decomposition
// (c) Cao Siqin 2013.10.9  (revised) : inversion revised
// (c) Cao Siqin 2013.9.26  (revised)
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
namespace MatrixNS {
    class MatrixDataContainer;
    class Matrix;
  namespace LUDecompositionNS {
    bool LUAlgebra(MatrixDataContainer * c, double ** a, int row, double * det, double **ainv);
  }
}
namespace MatrixNS {
    const char * Matrix_edition = "(c) Cao Siqin 2013.10.10";
    #ifndef MATRIX_DEFAULT_DIM
        #define MATRIX_DEFAULT_DIM 3
    #endif
    #ifndef MIN
      #define MIN(a,b) (a<b?a:b)
    #endif
    #ifndef MAX
      #define MAX(a,b) (a>b?a:b)
    #endif
    class MatrixDataContainer {
      private:
        int ip, savedip;
      public:
       #define SIZE_DEFAULT 400000
        char * root;
        int size;
        MatrixDataContainer * next;
      public:
        int getip(){ return ip; }
        void setip(int newip){ ip = newip; }
        void save_ip(){ savedip = ip; }
        void restore_ip(){
            if (savedip<size && savedip>=0) ip = savedip; }
        void init(int size_=SIZE_DEFAULT){
            size = size_; if (size<SIZE_DEFAULT) size=SIZE_DEFAULT;
            root = 0; next = 0; ip = 0; savedip = 0;
            root = (char*)malloc(sizeof(char)*size);
        }
        void enlarge(int size_=SIZE_DEFAULT){
            if (size_<SIZE_DEFAULT) size_=SIZE_DEFAULT;
            next = (MatrixDataContainer *)malloc(sizeof(MatrixDataContainer)+sizeof(char)*size_);
            next->size = size_; next->ip = 0;
            next->root = (char*)(((char*)next) + sizeof(MatrixDataContainer));
            next->next = 0;
        }
        void dispose(bool free_root=true){
            if (next){ next->dispose(false); ::free(next); }
            if (free_root) ::free(root);
        }
        char * allocate(int count){
            if (ip+count < size){
                int ir = ip; ip += count;
                for (int i=ir; i<ip; i++) root[i] = 0;
                return &root[ir];
            } else {
                if (!next) enlarge(MAX(size, count));
                return next->allocate(count);
            }
        }
        void notify_free(int ipbegin, int ipend){
        }
        void reset(){
            savedip = 0;
            if (next) next->reset(); ip = 0; for (int i=0; i<size; i++) root[i] = 0;
        }
    };
    class Matrix {
      public:
        int     n;  //dimension
        double* m;  //data array
        double**a;  //2D data array
      private:
        double* m_new;
        double**a_new;
        double* __a[MATRIX_DEFAULT_DIM];
        double  __m[MATRIX_DEFAULT_DIM*MATRIX_DEFAULT_DIM];
        MatrixDataContainer * __c;
        int __cip, __cipe;
      public:
       //construct and dispose
        void set_container(MatrixDataContainer * cont){ __c = cont; }
        void init(MatrixDataContainer * cont, int n_){ init(n_, cont); }
        void init(int n, double * m){ __c = 0; this->n = n; this->m = m; m_new = 0; }
        void init(int n, MatrixDataContainer * cont=nullptr){ this->n = n; __c = cont;
            if (n<=MATRIX_DEFAULT_DIM){
                a = __a; m = __m;
            } else {
                __cip = __c? __c->getip() : 0;
                int lenh = sizeof(double*)*n;
                int lend = sizeof(double)*n*n;
                char * p = __c?__c->allocate(lenh+lend):(char*)malloc(lenh+lend);
                a = a_new = (double**) p;
                m = m_new = (double*) (p + lenh);
                __cipe = __c? __c->getip() : 0;
            }
            for (int i=0; i<n; i++) a[i] = &m[i*n];
            for (int i=0; i<n*n; i++) m[i] = 0;
        }
        void dispose(){
            this->n = 0; this->m = 0;
            if (m_new){
                if (__c) __c->notify_free(__cip, __cipe);
            }
        }
        ~Matrix(){ dispose(); }
        Matrix(){ init(0); }
        Matrix(MatrixDataContainer * cont){ init(cont, 1); }
        Matrix(MatrixDataContainer * cont, int n_){ init(cont, n_); }
        Matrix(int n){ init(n); }
        Matrix(int n, double * m){ init(n, m); }
      public:
      //basic matrix
        Matrix operator =(double lambda){
            for (int i=0; i<n; i++) for (int j=0; j<n; j++){
                if (i==j) *e(i,j) = lambda; else *e(i,j) = 0;
            }
            return *this;
        }
        Matrix operator =(Matrix o){
            if (!m){ memcpy(this, &o, sizeof(Matrix)); init(__c, n); }
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++) *e(i,j) = *o.e(i,j);
            return *this;
        }
      //basic algorithm
        inline double * e(int a, int b){ return &m[a*n+b]; }
        Matrix operator +=(Matrix o){
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++) *e(i,j) += *o.e(i,j);
            return *this;
        }
        Matrix operator -=(Matrix o){
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++) *e(i,j) -= *o.e(i,j);
            return *this;
        }
        Matrix operator *=(double l){
            for (int i=0; i<n; i++) for (int j=0; j<n; j++) *e(i,j) *= l;
            return *this;
        }
        Matrix operator /=(double l){
            for (int i=0; i<n; i++) for (int j=0; j<n; j++) *e(i,j) /= l;
            return *this;
        }
        Matrix operator +(Matrix o){
            int cip = __c? __c->getip() : 0;
            Matrix ret(__c, MIN(n,o.n));
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++){
                *ret.e(i,j) = *e(i,j) + *o.e(i,j);
            }
            if (__c) __c->setip(cip);
            return ret;
        }
        Matrix operator -(Matrix o){
            int cip = __c? __c->getip() : 0;
            Matrix ret(__c, MIN(n,o.n));
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++){
                *ret.e(i,j) = *e(i,j) - *o.e(i,j);
            }
            if (__c) __c->setip(cip);
            return ret;
        }
        Matrix product(Matrix o, Matrix out){
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++){
                *out.e(i,j) = 0;
                for (int k=0; k<MIN(n,o.n); k++) *out.e(i,j) += *e(i,k) * *o.e(k,j);
            }
            return out;
        }
        Matrix operator *(Matrix o){
            int cip = __c? __c->getip() : 0;
            Matrix ret(__c, MIN(n,o.n));
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++){
                *ret.e(i,j) = 0;
                for (int k=0; k<MIN(n,o.n); k++) *ret.e(i,j) += *e(i,k) * *o.e(k,j);
            }
            if (__c) __c->setip(cip);
            return ret;
        }
      //Transpose
        void Transpose(){
            for (int i=0; i<n; i++) for (int j=i; j<n; j++){
                double t = *e(i,j); *e(i,j) = *e(j,i); *e(j,i) = t;
            }
        }
      //Trace
        double trace(){
            double tr = 0; for (int i=0; i<n; i++) tr += *e(i,i); return tr;
        }
      public:
       //line transformation: basic transformation
        void row_weight_add_to(int a, int b, double k, int c0=0, int cn=-1){ //row a product k then add to row b (B += k*A), from col c0 to cn
            if (cn<c0) cn = n-1; for (int i=c0; i<=cn; i++) *e(b, i) += *e(a, i) * k;
        }
        void col_weight_add_to(int a, int b, double k, int c0=0, int cn=-1){ //col a product k then add to col b (B += k*A), from row c0 to cn
            if (cn<c0) cn = n-1; for (int i=c0; i<=cn; i++) *e(i, b) += *e(i, a) * k;
        }
      public:
       //diagnol renormalization: self renormalization, eleminate diagnol 0 element
        bool renom_diagnol(Matrix oi){
            /*for (int k=0; k<n; k++) if (*e(k,k)==0){
                int l=-1; for (int kk=0; kk<n; kk++){ if (*e(kk,k)!=0){ l = kk; break; } }
                if (l<0) return false;
            }*/
            for (int k=0; k<n; k++) if (*e(k,k)==0){
                int l=-1; for (int kk=0; kk<n; kk++){ if (*e(kk,k)!=0){ l = kk; break; } }
                if (l<0) return false;
                row_weight_add_to(l, k, 1);
                oi.row_weight_add_to(l, k, 1);
            }
            return true;
        }
        Matrix up_triangle(Matrix oi){
            bool reversible = true;
            reversible = renom_diagnol(oi);
            if (reversible) for (int u=0; u<=n-2; u++){
                for (int k=u+1; k<=n-1; k++){
                  if (*e(u,u) == 0){
                    reversible = false;
                  } else {
                    double a = - *e(k,u) / *e(u,u);
                    row_weight_add_to(u, k, a, u+1); *e(k,u) = 0;
                    oi.row_weight_add_to(u, k, a);
                  }
                }
            }
            if (!reversible) return Matrix(0);
            return *this;
        }
        double determin(bool * reversible = nullptr){
            int cip = __c? __c->getip() : 0;
            Matrix tmp(__c, n); tmp = *this;
            Matrix tmp2 = tmp.up_triangle(Matrix(0));
              //if (tmp2.n!=tmp.n){ if (reversible) *reversible = false; else *reversible = true; } //�治����
            double ret = 1; for (int i=0; i<n; i++) ret *= *tmp2.e(i,i);
            tmp.dispose();
            if (__c) __c->setip(cip);
            return ret;
        }
      //Inversion
        Matrix inverse(){
            int cip = __c? __c->getip() : 0;
            Matrix oi(__c, n);  Matrix ol(__c, n);
            if (__c) __c->setip(cip);
            if (inverse(oi, ol)){
                return oi;
            } else return Matrix(0);
        }
        bool inverse(Matrix oi, Matrix ol){ //ol is a temporary matrix, which can be *this
            int cip = __c? __c->getip() : 0;
            ol = *this; oi = 1;
          //1&2 lower triangle transformation: lower triangle transformed to 0
            Matrix t = ol.up_triangle(oi); if (t.n!=ol.n) return false;  //no inversion
          //3. check the diagnol
            //for (int i=0; i<ol.n; i++)  if (*ol.e(i,i)==0) return false;  //no inversion
          //4. upper triangle transformation: upper triangle transformed to 0
            for (int u=ol.n-1; u>=1; u--){
                for (int k=u-1; k>=0; k--){
                    double a = - *ol.e(k,u) / *ol.e(u,u);
                    ol.row_weight_add_to(u, k, a, 0, u-1); *ol.e(k,u) = 0;
                    oi.row_weight_add_to(u, k, a);
                }
            }
          //5. renormalized with diagnol elements
            for (int l=0; l<ol.n; l++){
                double factor = *ol.e(l,l);
                for (int i=0; i<ol.n; i++){
                    *ol.e(l,i) /= factor;
                    *oi.e(l,i) /= factor;
                }
            }
            if (__c) __c->setip(cip);
            return true;
        }
      //LU Decomposition and inverse
        bool lu(Matrix oi, double * det){
            return LUDecompositionNS::LUAlgebra(__c, a, n, det, oi.a);
        }
    };
}

namespace MatrixNS {
  namespace LUDecompositionNS {
    double * init_double_vector(MatrixDataContainer * c, int n){
        if (c) return (double*) c->allocate(sizeof(double) * n);
        else return (double*) malloc(sizeof(double) * n);
    }
    int * init_int_vector(MatrixDataContainer * c, int n){
        if (c) return (int*) c->allocate(sizeof(int) * n);
        else return (int*) malloc(sizeof(int) * n);
    }
    double ** init_matrix(MatrixDataContainer * c, int m, int n){
        int lenh = sizeof(double*) * m; int len = lenh + sizeof(double) * m * n;
        char * p; if (c) p = (char*) c->allocate(len); else p = (char*) malloc(len);
        double * d = (double*)(p + lenh);
        double ** a = (double**) p; for (int i=0; i<m; i++) a[i] = &d[i*n];
        return a;
    }
    void solve1D(double ** lu, int n, int * indx, double * b, double * x){
        int i,ii=0,ip,j; double sum;
        for (i=0;i<n;i++) x[i] = b[i];
        for (i=0;i<n;i++) {
            ip=indx[i];
            sum=x[ip];
            x[ip]=x[i];
            if (ii != 0) for (j=ii-1;j<i;j++) sum -= lu[i][j]*x[j];
            else if (sum != 0.0) ii=i+1;
            x[i]=sum;
        }
        for (i=n-1;i>=0;i--) {
            sum=x[i];
            for (j=i+1;j<n;j++) sum -= lu[i][j]*x[j];
            x[i]=sum/lu[i][i];
        }
    }
    void solve3D(MatrixDataContainer * c, double ** lu, int n, int * indx, double ** b, double ** x, int bcol){
        int i,j,m=bcol;
        double * xx = init_double_vector(c, n);
        for (j=0;j<m;j++) {
            for (i=0;i<n;i++) xx[i] = b[i][j];
            solve1D(lu, n, indx, xx,xx);
            for (i=0;i<n;i++) x[i][j] = xx[i];
        }
        if (!c) free(xx);
    }
    void inverse(MatrixDataContainer * c, double ** lu, int n, int * indx, double ** ainv){
        int i,j;
        for (i=0;i<n;i++) {
            for (j=0;j<n;j++) ainv[i][j] = 0;
            ainv[i][i] = 1.;
        }
        solve3D(c, lu, n, indx, ainv, ainv, n);
    }
    bool LUAlgebra(MatrixDataContainer * c, double ** a, int row, double * det, double **ainv){
        int cip = c? c->getip() : 0;
        int n; double ** lu; int * indx;
        bool success = true; double d;

      // 1. LU decomposition
        n = row; lu = init_matrix(c, n, n); indx = init_int_vector(c, n); d = 1;
        for (int i=0; i<n; i++) for (int j=0; j<n; j++) lu[i][j] = a[i][j];
        const double TINY=1.0e-40; int i,imax,j,k; double big,temp;
        double * vv = init_double_vector(c, n);
        for (i=0;i<n;i++) {
            big=0.0;
            for (j=0;j<n;j++)
                if ((temp=fabs(lu[i][j])) > big) big=temp;
            if (big == 0.0){ success = false; }// throw("Singular matrix in LUdcmp");
            vv[i]=1.0/big;
        }
        if (success) for (k=0;k<n;k++) {
            big=0.0;
            for (i=k;i<n;i++) {
                temp=vv[i]*fabs(lu[i][k]);
                if (temp > big) { big=temp; imax=i; }
            }
            if (k != imax) {
                for (j=0;j<n;j++) {
                    temp=lu[imax][j];
                    lu[imax][j]=lu[k][j];
                    lu[k][j]=temp;
                }
                d = -d;
                vv[imax]=vv[k];
            }
            indx[k]=imax;
            if (lu[k][k] == 0.0) lu[k][k]=TINY;
            for (i=k+1;i<n;i++) {
                temp=lu[i][k] /= lu[k][k];
                for (j=k+1;j<n;j++) lu[i][j] -= temp*lu[k][j];
            }
        }

      // 2. determinant
        if (success){
            *det = d; for (int i=0;i<n;i++) *det *= lu[i][i];
        }

      // 3. inverse
        if (success) inverse(c, lu, n, indx, ainv);

      // x. finalization
        if (!c){ free(vv); free(lu); free(indx); } else c->setip(cip);

        return true;
    }
  }
}
