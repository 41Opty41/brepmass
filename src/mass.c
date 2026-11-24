#include "brepmass.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#  define M_PI 3.1415926535897932384626433832795
#endif

typedef struct {
    int a;
    int b;
} bm_dir_edge_t;

static double bm_dot3(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void bm_cross3(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

static int cmp_dir_edge(const void *pa, const void *pb) {
    const bm_dir_edge_t *ea = (const bm_dir_edge_t *)pa;
    const bm_dir_edge_t *eb = (const bm_dir_edge_t *)pb;
    if (ea->a != eb->a) {
        return (ea->a < eb->a) ? -1 : 1;
    }
    if (ea->b != eb->b) {
        return (ea->b < eb->b) ? -1 : 1;
    }
    return 0;
}

static int find_dir_edge(const bm_dir_edge_t *arr, int n, int a, int b) {
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        if (arr[mid].a == a && arr[mid].b == b) {
            return mid;
        }
        if (arr[mid].a < a || (arr[mid].a == a && arr[mid].b < b)) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }
    return -1;
}

/* Cyclic Jacobi eigenvalue solver for 3x3 real symmetric inertia matrix */
static void bm_jacobi3(const double A_in[3][3], double evals[3], double evecs[3][3]) {
    double A[3][3];
    double V[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            A[i][j] = A_in[i][j];
            V[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int iter = 0; iter < 50; iter++) {
        int p = 0, q = 1;
        double mx = fabs(A[0][1]);
        if (fabs(A[0][2]) > mx) {
            p = 0;
            q = 2;
            mx = fabs(A[0][2]);
        }
        if (fabs(A[1][2]) > mx) {
            p = 1;
            q = 2;
            mx = fabs(A[1][2]);
        }
        if (mx < 1e-15) {
            break;
        }
        double app = A[p][p], aqq = A[q][q], apq = A[p][q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = (tau >= 0.0)
                       ? 1.0 / (tau + sqrt(1.0 + tau * tau))
                       : -1.0 / (-tau + sqrt(1.0 + tau * tau));
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;

        A[p][p] = app - t * apq;
        A[q][q] = aqq + t * apq;
        A[p][q] = 0.0;
        A[q][p] = 0.0;

        for (int r = 0; r < 3; r++) {
            if (r == p || r == q) {
                continue;
            }
            double arp = A[r][p], arq = A[r][q];
            A[r][p] = c * arp - s * arq;
            A[p][r] = A[r][p];
            A[r][q] = s * arp + c * arq;
            A[q][r] = A[r][q];
        }

        for (int r = 0; r < 3; r++) {
            double vrp = V[r][p], vrq = V[r][q];
            V[r][p] = c * vrp - s * vrq;
            V[r][q] = s * vrp + c * vrq;
        }
    }

    double lam[3] = {A[0][0], A[1][1], A[2][2]};
    int order[3] = {0, 1, 2};
    for (int i = 0; i < 2; i++) {
        for (int j = i + 1; j < 3; j++) {
            if (lam[order[j]] < lam[order[i]]) {
                int tmp = order[i];
                order[i] = order[j];
                order[j] = tmp;
            }
        }
    }

    for (int k = 0; k < 3; k++) {
        int col = order[k];
        evals[k] = lam[col];
        for (int r = 0; r < 3; r++) {
            evecs[r][k] = V[r][col];
        }
    }
}

static void bm_finalize_inertia(double mass, const double com[3], const double cov[3][3], bm_mass_props_t *out) {
    double trace = cov[0][0] + cov[1][1] + cov[2][2];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out->inertia_origin[i][j] = ((i == j) ? trace : 0.0) - cov[i][j];
        }
    }
    double c2 = bm_dot3(com, com);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out->inertia_centroid[i][j] =
                out->inertia_origin[i][j] - mass * (((i == j) ? c2 : 0.0) - com[i] * com[j]);
        }
    }
    bm_jacobi3(out->inertia_centroid, out->principal_moments, out->principal_axes);
}

int bm_mesh_is_closed(const int *triangles, int n_triangles, int n_vertices, int *out_euler) {
    if (!triangles || n_triangles <= 0 || n_vertices <= 0) {
        if (out_euler) {
            *out_euler = 0;
        }
        return (n_triangles == 0) ? 0 : -1;
    }

    int ne_dir = n_triangles * 3;
    bm_dir_edge_t *dir = (bm_dir_edge_t *)malloc(sizeof(bm_dir_edge_t) * (size_t)ne_dir);
    bm_dir_edge_t *undir = (bm_dir_edge_t *)malloc(sizeof(bm_dir_edge_t) * (size_t)ne_dir);
    unsigned char *used_v = (unsigned char *)calloc((size_t)n_vertices, 1);
    if (!dir || !undir || !used_v) {
        free(dir);
        free(undir);
        free(used_v);
        return -1;
    }

    for (int k = 0; k < n_triangles; k++) {
        int v[3] = {triangles[3 * k], triangles[3 * k + 1], triangles[3 * k + 2]};
        for (int e = 0; e < 3; e++) {
            int a = v[e], b = v[(e + 1) % 3];
            if (a < 0 || a >= n_vertices || b < 0 || b >= n_vertices) {
                free(dir);
                free(undir);
                free(used_v);
                return -1;
            }
            used_v[a] = 1;
            dir[3 * k + e].a = a;
            dir[3 * k + e].b = b;
            undir[3 * k + e].a = (a < b) ? a : b;
            undir[3 * k + e].b = (a < b) ? b : a;
        }
    }

    qsort(dir, (size_t)ne_dir, sizeof(bm_dir_edge_t), cmp_dir_edge);
    qsort(undir, (size_t)ne_dir, sizeof(bm_dir_edge_t), cmp_dir_edge);

    int v_used = 0;
    for (int i = 0; i < n_vertices; i++) {
        if (used_v[i]) {
            v_used++;
        }
    }
    int e_unique = 0;
    for (int i = 0; i < ne_dir; i++) {
        if (i == 0 || cmp_dir_edge(&undir[i], &undir[i - 1]) != 0) {
            e_unique++;
        }
    }
    if (out_euler) {
        *out_euler = v_used - e_unique + n_triangles;
    }

    int closed = 1;
    for (int i = 0; i < ne_dir; i++) {
        if (i > 0 && cmp_dir_edge(&dir[i], &dir[i - 1]) == 0) {
            closed = 0; /* Duplicate directed edge -> non-manifold */
            break;
        }
        if (find_dir_edge(dir, ne_dir, dir[i].b, dir[i].a) < 0) {
            closed = 0;
            break;
        }
    }

    free(dir);
    free(undir);
    free(used_v);
    return closed;
}

int bm_mesh_mass_props(
    const double *vertices,
    int n_vertices,
    const int *triangles,
    int n_triangles,
    double density,
    bm_mass_props_t *out
) {
    if (!vertices || !triangles || !out || n_vertices <= 0 || n_triangles <= 0 || density <= 0.0) {
        return -1;
    }
    memset(out, 0, sizeof(*out));
    out->n_vertices = n_vertices;
    out->n_triangles = n_triangles;

    int euler = 0;
    int closed = bm_mesh_is_closed(triangles, n_triangles, n_vertices, &euler);
    if (closed < 0) {
        return -1;
    }
    out->is_closed = closed;
    out->euler_char = euler;

    for (int j = 0; j < 3; j++) {
        out->bbox_min[j] = vertices[j];
        out->bbox_max[j] = vertices[j];
    }
    for (int i = 1; i < n_vertices; i++) {
        for (int j = 0; j < 3; j++) {
            double c = vertices[3 * i + j];
            if (c < out->bbox_min[j]) {
                out->bbox_min[j] = c;
            }
            if (c > out->bbox_max[j]) {
                out->bbox_max[j] = c;
            }
        }
    }

    double vol = 0.0;
    double area = 0.0;
    double cm_v[3] = {0.0, 0.0, 0.0};
    double cm_a[3] = {0.0, 0.0, 0.0};
    double S_v[3][3] = {{0}};
    double S_a[3][3] = {{0}};

    for (int k = 0; k < n_triangles; k++) {
        const double *a = vertices + 3 * triangles[3 * k];
        const double *b = vertices + 3 * triangles[3 * k + 1];
        const double *c = vertices + 3 * triangles[3 * k + 2];

        double e1[3], e2[3], n[3], cr[3];
        for (int i = 0; i < 3; i++) {
            e1[i] = b[i] - a[i];
            e2[i] = c[i] - a[i];
        }
        bm_cross3(e1, e2, n);
        double tri_a = 0.5 * sqrt(bm_dot3(n, n));
        area += tri_a;

        double Ts[3] = {a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]};
        for (int i = 0; i < 3; i++) {
            cm_a[i] += tri_a * Ts[i] / 3.0;
            for (int j = 0; j < 3; j++) {
                S_a[i][j] += (tri_a / 12.0) * (a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + Ts[i] * Ts[j]);
            }
        }

        if (closed) {
            bm_cross3(b, c, cr);
            double v = bm_dot3(a, cr) / 6.0;
            vol += v;
            for (int i = 0; i < 3; i++) {
                cm_v[i] += v * Ts[i] / 4.0;
                for (int j = 0; j < 3; j++) {
                    S_v[i][j] += (v / 20.0) * (a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + Ts[i] * Ts[j]);
                }
            }
        }
    }

    out->surface_area = area;
    if (closed) {
        if (vol < 0.0) {
            vol = -vol;
            for (int i = 0; i < 3; i++) {
                cm_v[i] = -cm_v[i];
                for (int j = 0; j < 3; j++) {
                    S_v[i][j] = -S_v[i][j];
                }
            }
        }
        out->dim = 3;
        out->volume = vol;
        out->mass = density * vol;
        if (vol > 0.0) {
            for (int i = 0; i < 3; i++) {
                out->centroid[i] = cm_v[i] / vol;
            }
        }
        double cov_mass[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                cov_mass[i][j] = density * S_v[i][j];
            }
        }
        bm_finalize_inertia(out->mass, out->centroid, cov_mass, out);
    } else {
        out->dim = 2;
        out->volume = 0.0;
        out->mass = density * area;
        if (area > 0.0) {
            for (int i = 0; i < 3; i++) {
                out->centroid[i] = cm_a[i] / area;
            }
        }
        double cov_mass[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                cov_mass[i][j] = density * S_a[i][j];
            }
        }
        bm_finalize_inertia(out->mass, out->centroid, cov_mass, out);
    }
    return 0;
}

int bm_analytic_box(double dx, double dy, double dz, double density, bm_mass_props_t *out) {
    if (!out || dx <= 0.0 || dy <= 0.0 || dz <= 0.0 || density <= 0.0) {
        return -1;
    }
    memset(out, 0, sizeof(*out));
    out->dim = 3;
    out->is_closed = 1;
    out->euler_char = 2;
    out->volume = dx * dy * dz;
    out->surface_area = 2.0 * (dx * dy + dy * dz + dx * dz);
    out->mass = density * out->volume;
    out->centroid[0] = 0.5 * dx;
    out->centroid[1] = 0.5 * dy;
    out->centroid[2] = 0.5 * dz;
    out->bbox_max[0] = dx;
    out->bbox_max[1] = dy;
    out->bbox_max[2] = dz;

    double m = out->mass;
    double ix = m * (dy * dy + dz * dz) / 12.0;
    double iy = m * (dx * dx + dz * dz) / 12.0;
    double iz = m * (dx * dx + dy * dy) / 12.0;
    out->inertia_centroid[0][0] = ix;
    out->inertia_centroid[1][1] = iy;
    out->inertia_centroid[2][2] = iz;
    double c2 = bm_dot3(out->centroid, out->centroid);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out->inertia_origin[i][j] =
                out->inertia_centroid[i][j] + m * (((i == j) ? c2 : 0.0) - out->centroid[i] * out->centroid[j]);
        }
    }
    bm_jacobi3(out->inertia_centroid, out->principal_moments, out->principal_axes);
    return 0;
}

int bm_analytic_cylinder(double radius, double height, double density, bm_mass_props_t *out) {
    return bm_analytic_cone_frustum(radius, radius, height, density, out);
}

int bm_analytic_cone_frustum(double r1, double r2, double height, double density, bm_mass_props_t *out) {
    if (!out || r1 < 0.0 || r2 < 0.0 || (r1 == 0.0 && r2 == 0.0) || height <= 0.0 || density <= 0.0) {
        return -1;
    }
    memset(out, 0, sizeof(*out));
    out->dim = 3;
    out->is_closed = 1;
    out->euler_char = 2;

    double r_max = (r1 > r2) ? r1 : r2;
    double s_quad = r1 * r1 + r1 * r2 + r2 * r2;
    double dr = r1 - r2;
    double slant = sqrt(height * height + dr * dr);

    out->volume = (M_PI * height / 3.0) * s_quad;
    out->surface_area = M_PI * (r1 * r1 + r2 * r2 + (r1 + r2) * slant);
    out->mass = density * out->volume;

    double zc = (height / 4.0) * (r1 * r1 + 2.0 * r1 * r2 + 3.0 * r2 * r2) / s_quad;
    out->centroid[2] = zc;
    out->bbox_min[0] = -r_max;
    out->bbox_min[1] = -r_max;
    out->bbox_min[2] = 0.0;
    out->bbox_max[0] = r_max;
    out->bbox_max[1] = r_max;
    out->bbox_max[2] = height;

    double s4 = r1 * r1 * r1 * r1 + r1 * r1 * r1 * r2 + r1 * r1 * r2 * r2 + r1 * r2 * r2 * r2 + r2 * r2 * r2 * r2;
    double iz = density * (M_PI * height / 10.0) * s4;
    double int_r2z2 = height * height * height * (r1 * r1 / 30.0 + r1 * r2 / 10.0 + r2 * r2 / 5.0);
    double ixy_origin = 0.5 * iz + density * M_PI * int_r2z2;

    out->inertia_origin[0][0] = ixy_origin;
    out->inertia_origin[1][1] = ixy_origin;
    out->inertia_origin[2][2] = iz;
    out->inertia_centroid[0][0] = ixy_origin - out->mass * zc * zc;
    out->inertia_centroid[1][1] = out->inertia_centroid[0][0];
    out->inertia_centroid[2][2] = iz;
    bm_jacobi3(out->inertia_centroid, out->principal_moments, out->principal_axes);
    return 0;
}

int bm_analytic_sphere(double radius, double density, bm_mass_props_t *out) {
    if (!out || radius <= 0.0 || density <= 0.0) {
        return -1;
    }
    memset(out, 0, sizeof(*out));
    out->dim = 3;
    out->is_closed = 1;
    out->euler_char = 2;
    out->volume = (4.0 / 3.0) * M_PI * radius * radius * radius;
    out->surface_area = 4.0 * M_PI * radius * radius;
    out->mass = density * out->volume;
    for (int i = 0; i < 3; i++) {
        out->bbox_min[i] = -radius;
        out->bbox_max[i] = radius;
    }
    double i0 = 0.4 * out->mass * radius * radius;
    for (int i = 0; i < 3; i++) {
        out->inertia_origin[i][i] = i0;
        out->inertia_centroid[i][i] = i0;
    }
    bm_jacobi3(out->inertia_centroid, out->principal_moments, out->principal_axes);
    return 0;
}

int bm_analytic_torus(double major_r, double minor_r, double density, bm_mass_props_t *out) {
    if (!out || major_r <= 0.0 || minor_r <= 0.0 || minor_r >= major_r || density <= 0.0) {
        return -1;
    }
    memset(out, 0, sizeof(*out));
    out->dim = 3;
    out->is_closed = 1;
    out->euler_char = 0;
    out->volume = 2.0 * M_PI * M_PI * major_r * minor_r * minor_r;
    out->surface_area = 4.0 * M_PI * M_PI * major_r * minor_r;
    out->mass = density * out->volume;
    out->bbox_min[0] = -(major_r + minor_r);
    out->bbox_min[1] = -(major_r + minor_r);
    out->bbox_min[2] = -minor_r;
    out->bbox_max[0] = major_r + minor_r;
    out->bbox_max[1] = major_r + minor_r;
    out->bbox_max[2] = minor_r;

    double R2 = major_r * major_r;
    double r2 = minor_r * minor_r;
    double ixx = out->mass * (0.5 * R2 + 0.625 * r2);
    double izz = out->mass * (R2 + 0.75 * r2);
    out->inertia_origin[0][0] = ixx;
    out->inertia_origin[1][1] = ixx;
    out->inertia_origin[2][2] = izz;
    out->inertia_centroid[0][0] = ixx;
    out->inertia_centroid[1][1] = ixx;
    out->inertia_centroid[2][2] = izz;
    bm_jacobi3(out->inertia_centroid, out->principal_moments, out->principal_axes);
    return 0;
}
