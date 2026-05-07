#include "helper.h"

#include <string.h>

pair<vector<vector<int>>, vector<vector<double>>> read_mm(int argc,
                                                          char* argv[]) {
    int ret_code;
    FILE* f;
    int M, N, nz;
    int i, j, row, col;
    double val;
    MM_typecode matcode;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    } else {
        if ((f = fopen(argv[1], "r")) == NULL) exit(1);
    }

    /* read banner to check matrix type */
    if (mm_read_banner(f, &matcode) != 0) {
        fprintf(stderr, "Could not process Matrix Market banner in file [%s]\n",
                argv[1]);
        exit(1);
    }

    if (!(mm_is_real(matcode) || mm_is_pattern(matcode) ||
          mm_is_integer(matcode)) ||
        !mm_is_matrix(matcode) || !mm_is_sparse(matcode)) {
        fprintf(stderr,
                "Sorry, this application does not support Matrix Market type: "
                "[%s]\n",
                mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0) exit(1);

    /* reserve memory for matrices */
    vector<vector<int>> A(M, vector<int>());
    vector<vector<double>> VAL(nz, vector<double>(3, 0));

    for (i = 0; i < nz; i++) {
        if (mm_is_pattern(matcode)) {
            fscanf(f, "%d %d\n", &row, &col);
            val = 1.0;
        } else {
            fscanf(f, "%d %d %lf\n", &row, &col, &val);
        }

        VAL[i] = vector<double>({double(row) - 1, double(col) - 1, val});
        if (row != col) {
            A[row - 1].push_back(col - 1);
            if (row > col) {
                A[col - 1].push_back(row - 1);
            }
        }
    }

    if (f != stdin) fclose(f);

    return make_pair(A, VAL);
}

void write_mm(vector<vector<int>> A, vector<vector<double>> VAL,
              vector<int> f_star, char* argv[]) {
    FILE *f, *f1;

    const char s[2] = ".";
    char fname[100];
    char f1name[100];
    char* token = strtok(argv[1], s);
    strcpy(fname, "./");
    strcpy(f1name, "./");
    strcat(fname, token);
    strcat(f1name, token);
    strcat(fname, "_converted.mtx");
    strcat(f1name, "_perm.mtx");

    MM_typecode matcode, matcode1;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);
    mm_set_symmetric(&matcode);
    mm_initialize_typecode(&matcode1);
    mm_set_matrix(&matcode1);
    mm_set_coordinate(&matcode1);
    mm_set_real(&matcode1);

    f = fopen(fname, "w");
    if (f == NULL) {
        cout << fname << endl;
        printf("%s\n", strerror(errno));
    }
    mm_write_banner(f, matcode);
    mm_write_mtx_crd_size(f, f_star.size(), f_star.size(), VAL.size());
    for (auto item : VAL) {
        int row = f_star[round(item[0])] + 1;
        int col = f_star[round(item[1])] + 1;
        if (row < col) {
            swap(row, col);
        }
        fprintf(f, "%d %d %lf\n", row, col, item[2]);
    }
    fclose(f);

    f1 = fopen(f1name, "w+");
    mm_write_banner(f, matcode1);
    mm_write_mtx_crd_size(f1, f_star.size(), f_star.size(), f_star.size());

    for (int i = 0; i < f_star.size(); i++)
        fprintf(f1, "%d %d %d\n", f_star[i] + 1, i + 1, 1);

    fclose(f1);
}

int bandwidth(vector<vector<int>> A, vector<int>& f) {
    int n = f.size();
    int bandwidth = 0;
    for (int i = 0; i < n; i++) {
        for (auto j : A[i]) {
            int bw = abs(f[j] - f[i]);
            if (bw > bandwidth) bandwidth = bw;
        }
    }
    return bandwidth;
}

bool cmp_pair(const pair<int, int> t1, const pair<int, int> t2) {
    return get<1>(t1) > get<1>(t2);
}

void sparse_matrix::flush_info(int v) {
    int min_val = this->n;
    int max_val = -1;
    int arg_min = -1;
    int arg_max = -1;
    for (int u : this->A[v]) {
        if (this->f[u] < min_val) {
            min_val = this->f[u];
            arg_min = u;
        }
        if (this->f[u] > max_val) {
            max_val = this->f[u];
            arg_max = u;
        }
    }
    this->f_min[v] = min_val;
    this->f_max[v] = max_val;
    int dist_min = this->f[v] - min_val;
    int dist_max = max_val - this->f[v];
    if (dist_min > dist_max) {
        this->B_f[v] = make_pair(dist_min, arg_min);
    } else {
        this->B_f[v] = make_pair(dist_max, arg_max);
    }
};

void sparse_matrix::flush_bandwidth() {
    this->bandwidth = -1;
    this->critical_vertices = unordered_set<int>();
    // #pragma omp parallel for num_threads(10)
    for (int v = 0; v < this->n; v++) {
        int bw_v = get<0>(this->B_f[v]);
        if (bw_v > this->bandwidth) {
            this->critical_vertices = unordered_set<int>();
            this->bandwidth = bw_v;
        }
        if (bw_v == this->bandwidth) this->critical_vertices.insert(v);
    }
}

void sparse_matrix::init_helper(vector<vector<int>> A, vector<int> f) {
    this->A = A;
    this->f = f;
    this->n = f.size();
    this->bandwidth = -1;

    // the first is the value, the second is the corresponding vertex
    this->f_inverse = vector<int>(n, -1);
    for (int i = 0; i < n; i++) this->f_inverse[f[i]] = i;
    this->B_f = vector<pair<int, int>>(n, make_pair(-1, -1));
    this->f_min = vector<int>(n, n);
    this->f_max = vector<int>(n, -1);
    this->critical_vertices = unordered_set<int>();
    for (int v = 0; v < n; v++) {
        sparse_matrix::flush_info(v);
    }
    sparse_matrix::flush_bandwidth();
}

sparse_matrix::sparse_matrix(vector<vector<int>> A) {
    int n = A.size();
    vector<int> f = vector<int>(n);
    for (int i = 0; i < n; i++) f[i] = i;
    sparse_matrix::init_helper(A, f);
}

sparse_matrix::sparse_matrix(vector<vector<int>> A, vector<int> f) {
    sparse_matrix::init_helper(A, f);
}

sparse_matrix::sparse_matrix(const sparse_matrix& spm) {
    this->A = spm.A;
    this->f = spm.f;
    this->f_inverse = spm.f_inverse;
    this->B_f = spm.B_f;
    this->f_min = spm.f_min;
    this->f_max = spm.f_max;
    this->critical_vertices = spm.critical_vertices;
    this->n = spm.n;
    this->bandwidth = spm.bandwidth;
}

void sparse_matrix::f_swap(int v, int u) {
    swap(this->f[u], this->f[v]);
    this->f_inverse[f[u]] = u;
    this->f_inverse[f[v]] = v;

    sparse_matrix::flush_info(v);
    sparse_matrix::flush_info(u);
    for (int w : this->A[v]) {
        sparse_matrix::flush_info(w);
    }
    for (int w : this->A[u]) {
        sparse_matrix::flush_info(w);
    }

    sparse_matrix::flush_bandwidth();
}

void sparse_matrix::f_assign(vector<int> vs, vector<int> vals) {
    unordered_set<int> W;
    for (int i = 0; i < vs.size(); i++) {
        int v = vs[i];
        this->f[v] = vals[i];
        this->f_inverse[f[v]] = v;
        W.insert(v);
        copy(this->A[v].begin(), this->A[v].end(), inserter(W, W.end()));
    }
    for (int w : W) {
        sparse_matrix::flush_info(w);
    }
    sparse_matrix::flush_bandwidth();
}

sparse_matrix::~sparse_matrix() {}
