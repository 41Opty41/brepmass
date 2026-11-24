#include "brepmass.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#  define M_PI 3.1415926535897932384626433832795
#endif

void bm_mesh_init(bm_mesh_t *mesh) {
    if (!mesh) {
        return;
    }
    memset(mesh, 0, sizeof(*mesh));
}

void bm_mesh_free(bm_mesh_t *mesh) {
    if (!mesh) {
        return;
    }
    free(mesh->vertices);
    free(mesh->triangles);
    free(mesh->face_ids);
    memset(mesh, 0, sizeof(*mesh));
}

int bm_mesh_reserve(bm_mesh_t *mesh, int n_vertices, int n_triangles) {
    if (!mesh) {
        return -1;
    }
    if (n_vertices > mesh->cap_vertices) {
        int nc = (n_vertices < 16) ? 16 : n_vertices;
        double *nv = (double *)realloc(mesh->vertices, sizeof(double) * 3 * (size_t)nc);
        if (!nv) {
            return -1;
        }
        mesh->vertices = nv;
        mesh->cap_vertices = nc;
    }
    if (n_triangles > mesh->cap_triangles) {
        int nc = (n_triangles < 16) ? 16 : n_triangles;
        int *nt = (int *)realloc(mesh->triangles, sizeof(int) * 3 * (size_t)nc);
        if (!nt) {
            return -1;
        }
        mesh->triangles = nt;
        int *nf = (int *)realloc(mesh->face_ids, sizeof(int) * (size_t)nc);
        if (!nf) {
            return -1;
        }
        mesh->face_ids = nf;
        mesh->cap_triangles = nc;
    }
    return 0;
}

int bm_mesh_add_vertex(bm_mesh_t *mesh, double x, double y, double z) {
    if (!mesh) {
        return -1;
    }
    if (mesh->n_vertices >= mesh->cap_vertices) {
        int next_cap = mesh->cap_vertices ? mesh->cap_vertices * 2 : 64;
        if (bm_mesh_reserve(mesh, next_cap, mesh->cap_triangles) != 0) {
            return -1;
        }
    }
    int idx = mesh->n_vertices++;
    mesh->vertices[3 * idx] = x;
    mesh->vertices[3 * idx + 1] = y;
    mesh->vertices[3 * idx + 2] = z;
    return idx;
}

int bm_mesh_add_triangle(bm_mesh_t *mesh, int i0, int i1, int i2, int face_id) {
    if (!mesh || i0 < 0 || i1 < 0 || i2 < 0) {
        return -1;
    }
    if (mesh->n_triangles >= mesh->cap_triangles) {
        int next_cap = mesh->cap_triangles ? mesh->cap_triangles * 2 : 64;
        if (bm_mesh_reserve(mesh, mesh->cap_vertices, next_cap) != 0) {
            return -1;
        }
    }
    int idx = mesh->n_triangles++;
    mesh->triangles[3 * idx] = i0;
    mesh->triangles[3 * idx + 1] = i1;
    mesh->triangles[3 * idx + 2] = i2;
    mesh->face_ids[idx] = face_id;
    return idx;
}

bm_tess_opts_t bm_tess_opts_default(void) {
    bm_tess_opts_t o;
    o.u_segments = 64;
    o.v_segments = 32;
    o.chord_tol = 0.0;
    return o;
}

int bm_tessellate_box(double dx, double dy, double dz, bm_mesh_t *mesh) {
    if (!mesh || dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
        return -1;
    }
    const double coords[8][3] = {
        {0.0, 0.0, 0.0}, {dx, 0.0, 0.0}, {dx, dy, 0.0}, {0.0, dy, 0.0},
        {0.0, 0.0, dz},  {dx, 0.0, dz},  {dx, dy, dz},  {0.0, dy, dz},
    };
    const int tris[12][4] = {
        {0, 2, 1, 0}, {0, 3, 2, 0}, /* bottom z=0 (-Z) */
        {4, 5, 6, 1}, {4, 6, 7, 1}, /* top z=dz (+Z) */
        {0, 1, 5, 2}, {0, 5, 4, 2}, /* y=0 (-Y) */
        {1, 2, 6, 3}, {1, 6, 5, 3}, /* x=dx (+X) */
        {2, 3, 7, 4}, {2, 7, 6, 4}, /* y=dy (+Y) */
        {3, 0, 4, 5}, {3, 4, 7, 5}, /* x=0 (-X) */
    };
    int base = mesh->n_vertices;
    for (int i = 0; i < 8; i++) {
        if (bm_mesh_add_vertex(mesh, coords[i][0], coords[i][1], coords[i][2]) < 0) {
            return -1;
        }
    }
    for (int i = 0; i < 12; i++) {
        if (bm_mesh_add_triangle(mesh, base + tris[i][0], base + tris[i][1], base + tris[i][2], tris[i][3]) < 0) {
            return -1;
        }
    }
    return 0;
}

int bm_tessellate_cylinder(double radius, double height, const bm_tess_opts_t *opts, bm_mesh_t *mesh) {
    return bm_tessellate_cone_frustum(radius, radius, height, opts, mesh);
}

int bm_tessellate_cone_frustum(double r1, double r2, double height, const bm_tess_opts_t *opts, bm_mesh_t *mesh) {
    if (!mesh || r1 < 0.0 || r2 < 0.0 || (r1 == 0.0 && r2 == 0.0) || height <= 0.0) {
        return -1;
    }
    bm_tess_opts_t o = opts ? *opts : bm_tess_opts_default();
    int nu = (o.u_segments >= 3) ? o.u_segments : 64;
    int nv = (o.v_segments >= 1) ? o.v_segments : 1;

    /* Rings from z=0 (r1) to z=height (r2) */
    int *ring = (int *)malloc(sizeof(int) * (size_t)((nv + 1) * nu));
    if (!ring) {
        return -1;
    }
    for (int iv = 0; iv <= nv; iv++) {
        double t = (double)iv / (double)nv;
        double r = r1 + (r2 - r1) * t;
        double z = height * t;
        if (r == 0.0) {
            int apex = bm_mesh_add_vertex(mesh, 0.0, 0.0, z);
            if (apex < 0) {
                free(ring);
                return -1;
            }
            for (int iu = 0; iu < nu; iu++) {
                ring[iv * nu + iu] = apex;
            }
        } else {
            for (int iu = 0; iu < nu; iu++) {
                double th = 2.0 * M_PI * (double)iu / (double)nu;
                int v = bm_mesh_add_vertex(mesh, r * cos(th), r * sin(th), z);
                if (v < 0) {
                    free(ring);
                    return -1;
                }
                ring[iv * nu + iu] = v;
            }
        }
    }

    int err = 0;
    for (int iv = 0; iv < nv && !err; iv++) {
        int collapse_lo = (iv == 0 && r1 == 0.0);
        int collapse_hi = (iv + 1 == nv && r2 == 0.0);
        for (int iu = 0; iu < nu; iu++) {
            int inext = (iu + 1) % nu;
            int a = ring[iv * nu + iu];
            int b = ring[iv * nu + inext];
            int c = ring[(iv + 1) * nu + inext];
            int d = ring[(iv + 1) * nu + iu];
            if (!collapse_lo) {
                if (bm_mesh_add_triangle(mesh, a, b, c, 0) < 0) {
                    err = -1;
                    break;
                }
            }
            if (!collapse_hi) {
                if (bm_mesh_add_triangle(mesh, a, c, d, 0) < 0) {
                    err = -1;
                    break;
                }
            }
        }
    }

    if (!err && r1 > 0.0) {
        int c_bot = bm_mesh_add_vertex(mesh, 0.0, 0.0, 0.0);
        if (c_bot < 0) {
            err = -1;
        }
        for (int iu = 0; iu < nu && !err; iu++) {
            int inext = (iu + 1) % nu;
            if (bm_mesh_add_triangle(mesh, c_bot, ring[inext], ring[iu], 1) < 0) {
                err = -1;
            }
        }
    }
    if (!err && r2 > 0.0) {
        int c_top = bm_mesh_add_vertex(mesh, 0.0, 0.0, height);
        if (c_top < 0) {
            err = -1;
        }
        int base_top = nv * nu;
        for (int iu = 0; iu < nu && !err; iu++) {
            int inext = (iu + 1) % nu;
            if (bm_mesh_add_triangle(mesh, c_top, ring[base_top + iu], ring[base_top + inext], 2) < 0) {
                err = -1;
            }
        }
    }

    free(ring);
    return err;
}

int bm_tessellate_sphere(double radius, const bm_tess_opts_t *opts, bm_mesh_t *mesh) {
    if (!mesh || radius <= 0.0) {
        return -1;
    }
    bm_tess_opts_t o = opts ? *opts : bm_tess_opts_default();
    int nu = (o.u_segments >= 4) ? o.u_segments : 64;
    int nv = (o.v_segments >= 3) ? o.v_segments : 32;

    int south = bm_mesh_add_vertex(mesh, 0.0, 0.0, -radius);
    int north = bm_mesh_add_vertex(mesh, 0.0, 0.0, radius);
    if (south < 0 || north < 0) {
        return -1;
    }

    int *rings = (int *)malloc(sizeof(int) * (size_t)((nv - 1) * nu));
    if (!rings) {
        return -1;
    }
    for (int iv = 1; iv < nv; iv++) {
        double phi = -0.5 * M_PI + M_PI * (double)iv / (double)nv;
        double cp = cos(phi), sp = sin(phi);
        for (int iu = 0; iu < nu; iu++) {
            double th = 2.0 * M_PI * (double)iu / (double)nu;
            int v = bm_mesh_add_vertex(mesh, radius * cp * cos(th), radius * cp * sin(th), radius * sp);
            if (v < 0) {
                free(rings);
                return -1;
            }
            rings[(iv - 1) * nu + iu] = v;
        }
    }

    int err = 0;
    for (int iu = 0; iu < nu && !err; iu++) {
        int inext = (iu + 1) % nu;
        if (bm_mesh_add_triangle(mesh, south, rings[inext], rings[iu], 0) < 0) {
            err = -1;
        }
    }
    for (int iv = 0; iv < nv - 2 && !err; iv++) {
        for (int iu = 0; iu < nu; iu++) {
            int inext = (iu + 1) % nu;
            int a = rings[iv * nu + iu];
            int b = rings[iv * nu + inext];
            int c = rings[(iv + 1) * nu + inext];
            int d = rings[(iv + 1) * nu + iu];
            if (bm_mesh_add_triangle(mesh, a, b, c, 0) < 0 || bm_mesh_add_triangle(mesh, a, c, d, 0) < 0) {
                err = -1;
                break;
            }
        }
    }
    int last_row = (nv - 2) * nu;
    for (int iu = 0; iu < nu && !err; iu++) {
        int inext = (iu + 1) % nu;
        if (bm_mesh_add_triangle(mesh, rings[last_row + iu], rings[last_row + inext], north, 0) < 0) {
            err = -1;
        }
    }

    free(rings);
    return err;
}

int bm_tessellate_torus(double major_r, double minor_r, const bm_tess_opts_t *opts, bm_mesh_t *mesh) {
    if (!mesh || major_r <= 0.0 || minor_r <= 0.0 || minor_r >= major_r) {
        return -1;
    }
    bm_tess_opts_t o = opts ? *opts : bm_tess_opts_default();
    int nu = (o.u_segments >= 4) ? o.u_segments : 64;
    int nv = (o.v_segments >= 4) ? o.v_segments : 32;

    int *grid = (int *)malloc(sizeof(int) * (size_t)(nu * nv));
    if (!grid) {
        return -1;
    }
    for (int iu = 0; iu < nu; iu++) {
        double u = 2.0 * M_PI * (double)iu / (double)nu;
        double cu = cos(u), su = sin(u);
        for (int iv = 0; iv < nv; iv++) {
            double v = 2.0 * M_PI * (double)iv / (double)nv;
            double rad = major_r + minor_r * cos(v);
            int vid = bm_mesh_add_vertex(mesh, rad * cu, rad * su, minor_r * sin(v));
            if (vid < 0) {
                free(grid);
                return -1;
            }
            grid[iu * nv + iv] = vid;
        }
    }

    int err = 0;
    for (int iu = 0; iu < nu && !err; iu++) {
        int iu_n = (iu + 1) % nu;
        for (int iv = 0; iv < nv; iv++) {
            int iv_n = (iv + 1) % nv;
            int a = grid[iu * nv + iv];
            int b = grid[iu_n * nv + iv];
            int c = grid[iu_n * nv + iv_n];
            int d = grid[iu * nv + iv_n];
            if (bm_mesh_add_triangle(mesh, a, b, c, 0) < 0 || bm_mesh_add_triangle(mesh, a, c, d, 0) < 0) {
                err = -1;
                break;
            }
        }
    }
    free(grid);
    return err;
}

/* ---------------- 2D Ear-Clipping with Hole Bridging ---------------- */

static double area2_signed(const double *xy, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        s += xy[2 * i] * xy[2 * j + 1] - xy[2 * j] * xy[2 * i + 1];
    }
    return 0.5 * s;
}

static double cross2(double ax, double ay, double bx, double by, double cx, double cy) {
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

static int point_in_tri2(
    double px, double py,
    double ax, double ay,
    double bx, double by,
    double cx, double cy
) {
    double c1 = cross2(ax, ay, bx, by, px, py);
    double c2 = cross2(bx, by, cx, cy, px, py);
    double c3 = cross2(cx, cy, ax, ay, px, py);
    return (c1 >= -1e-12 && c2 >= -1e-12 && c3 >= -1e-12);
}

static void seg_ray_hit(
    double ox, double oy,
    double ax, double ay,
    double bx, double by,
    int *hit, double *hx
) {
    *hit = 0;
    if ((ay <= oy && by > oy) || (by <= oy && ay > oy)) {
        double t = (oy - ay) / (by - ay);
        double x = ax + t * (bx - ax);
        if (x >= ox - 1e-12) {
            *hit = 1;
            *hx = x;
        }
    }
}

int bm_tessellate_polygon_2d(
    const double *outer_xy,
    int n_outer,
    const double *holes_xy,
    const int *hole_counts,
    int n_holes,
    const double origin[3],
    const double normal[3],
    const double x_axis[3],
    bm_mesh_t *mesh
) {
    if (!outer_xy || n_outer < 3 || !origin || !normal || !x_axis || !mesh) {
        return -1;
    }

    int total_pts = n_outer;
    for (int h = 0; h < n_holes; h++) {
        if (!holes_xy || !hole_counts || hole_counts[h] < 3) {
            return -1;
        }
        total_pts += hole_counts[h];
    }

    /* Build orthonormal basis (X, Y, Z) */
    double z[3] = {normal[0], normal[1], normal[2]};
    double zn = sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
    if (zn < 1e-15) {
        return -1;
    }
    z[0] /= zn; z[1] /= zn; z[2] /= zn;

    double xd = x_axis[0] * z[0] + x_axis[1] * z[1] + x_axis[2] * z[2];
    double x[3] = {x_axis[0] - xd * z[0], x_axis[1] - xd * z[1], x_axis[2] - xd * z[2]};
    double xn = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    if (xn < 1e-12) {
        double cand[3] = {1.0, 0.0, 0.0};
        if (fabs(z[0]) > 0.9) {
            cand[0] = 0.0; cand[1] = 1.0;
        }
        x[0] = z[1] * cand[2] - z[2] * cand[1];
        x[1] = z[2] * cand[0] - z[0] * cand[2];
        x[2] = z[0] * cand[1] - z[1] * cand[0];
        xn = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    }
    x[0] /= xn; x[1] /= xn; x[2] /= xn;
    double y[3] = {
        z[1] * x[2] - z[2] * x[1],
        z[2] * x[0] - z[0] * x[2],
        z[0] * x[1] - z[1] * x[0],
    };

    /* Store all distinct geometric 2D vertices & ring index list */
    int max_ring = total_pts + 2 * n_holes + 4;
    double *pts = (double *)malloc(sizeof(double) * 2 * (size_t)total_pts);
    int *vmap = (int *)malloc(sizeof(int) * (size_t)total_pts);
    int *ring = (int *)malloc(sizeof(int) * (size_t)max_ring);
    if (!pts || !vmap || !ring) {
        free(pts); free(vmap); free(ring);
        return -1;
    }

    int out_ccw = area2_signed(outer_xy, n_outer) >= 0.0;
    for (int i = 0; i < n_outer; i++) {
        int src = out_ccw ? i : (n_outer - 1 - i);
        pts[2 * i] = outer_xy[2 * src];
        pts[2 * i + 1] = outer_xy[2 * src + 1];
        vmap[i] = bm_mesh_add_vertex(
            mesh,
            origin[0] + pts[2 * i] * x[0] + pts[2 * i + 1] * y[0],
            origin[1] + pts[2 * i] * x[1] + pts[2 * i + 1] * y[1],
            origin[2] + pts[2 * i] * x[2] + pts[2 * i + 1] * y[2]
        );
        ring[i] = i;
    }
    int n_ring = n_outer;
    int pt_off = n_outer;

    /* Add holes (CW) and splice each hole via horizontal-ray visible bridge */
    int hole_raw_off = 0;
    for (int h = 0; h < n_holes; h++) {
        int nh = hole_counts[h];
        const double *hxy = holes_xy + 2 * hole_raw_off;
        hole_raw_off += nh;

        int h_ccw = area2_signed(hxy, nh) >= 0.0;
        int max_i = 0;
        for (int i = 0; i < nh; i++) {
            int src = h_ccw ? (nh - 1 - i) : i; /* Holes must be CW */
            pts[2 * (pt_off + i)] = hxy[2 * src];
            pts[2 * (pt_off + i) + 1] = hxy[2 * src + 1];
            vmap[pt_off + i] = bm_mesh_add_vertex(
                mesh,
                origin[0] + pts[2 * (pt_off + i)] * x[0] + pts[2 * (pt_off + i) + 1] * y[0],
                origin[1] + pts[2 * (pt_off + i)] * x[1] + pts[2 * (pt_off + i) + 1] * y[1],
                origin[2] + pts[2 * (pt_off + i)] * x[2] + pts[2 * (pt_off + i) + 1] * y[2]
            );
            if (pts[2 * (pt_off + i)] > pts[2 * (pt_off + max_i)]) {
                max_i = i;
            }
        }

        double mx = pts[2 * (pt_off + max_i)];
        double my = pts[2 * (pt_off + max_i) + 1];

        /* Cast +X ray from (mx, my) against current outer ring edges to find mutually visible outer vertex */
        double best_hx = 1e300;
        int best_edge = 0;
        for (int i = 0; i < n_ring; i++) {
            int a = ring[i], b = ring[(i + 1) % n_ring];
            int hit = 0;
            double hx = 0.0;
            seg_ray_hit(mx, my, pts[2 * a], pts[2 * a + 1], pts[2 * b], pts[2 * b + 1], &hit, &hx);
            if (hit && hx < best_hx) {
                best_hx = hx;
                best_edge = i;
            }
        }
        int va = ring[best_edge], vb = ring[(best_edge + 1) % n_ring];
        int bridge_pos = (pts[2 * va] >= pts[2 * vb]) ? best_edge : ((best_edge + 1) % n_ring);
        int P_v = ring[bridge_pos];
        double px = pts[2 * P_v], py = pts[2 * P_v + 1];

        /* Check if any reflex vertex inside triangle (M, (best_hx, my), P) gives closer visible candidate */
        double best_dist2 = (px - mx) * (px - mx) + (py - my) * (py - my);
        for (int i = 0; i < n_ring; i++) {
            int vid = ring[i];
            if (vid == P_v) {
                continue;
            }
            double vx = pts[2 * vid], vy = pts[2 * vid + 1];
            if (vx >= mx && point_in_tri2(vx, vy, mx, my, best_hx, my, px, py)) {
                double d2 = (vx - mx) * (vx - mx) + (vy - my) * (vy - my);
                if (d2 < best_dist2) {
                    best_dist2 = d2;
                    bridge_pos = i;
                }
            }
        }

        /* Splice hole sequence into ring after bridge_pos: ... bridge_v, M, M+1 .. M-1, M, bridge_v ... */
        int insert_len = nh + 2;
        memmove(
            &ring[bridge_pos + 1 + insert_len],
            &ring[bridge_pos + 1],
            sizeof(int) * (size_t)(n_ring - (bridge_pos + 1))
        );
        int wp = bridge_pos + 1;
        for (int k = 0; k < nh; k++) {
            ring[wp++] = pt_off + ((max_i + k) % nh);
        }
        ring[wp++] = pt_off + max_i;
        ring[wp++] = ring[bridge_pos];
        n_ring += insert_len;
        pt_off += nh;
    }

    /* Simple ear-clipping on merged CCW polygon */
    int guard = n_ring * n_ring + 16;
    int pos = 0;
    while (n_ring >= 3 && guard-- > 0) {
        int clipped = 0;
        for (int step = 0; step < n_ring; step++) {
            int idx = (pos + step) % n_ring;
            int prev = (idx + n_ring - 1) % n_ring;
            int next = (idx + 1) % n_ring;
            int a = ring[prev], b = ring[idx], c = ring[next];
            double ax = pts[2 * a], ay = pts[2 * a + 1];
            double bx = pts[2 * b], by = pts[2 * b + 1];
            double cx = pts[2 * c], cy = pts[2 * c + 1];
            if (cross2(ax, ay, bx, by, cx, cy) <= 1e-12) {
                continue;
            }
            int contains = 0;
            for (int j = 0; j < n_ring; j++) {
                if (j == prev || j == idx || j == next) {
                    continue;
                }
                int q = ring[j];
                if (q == a || q == b || q == c) {
                    continue;
                }
                if (point_in_tri2(pts[2 * q], pts[2 * q + 1], ax, ay, bx, by, cx, cy)) {
                    contains = 1;
                    break;
                }
            }
            if (!contains) {
                if (bm_mesh_add_triangle(mesh, vmap[a], vmap[b], vmap[c], 0) < 0) {
                    free(pts); free(vmap); free(ring);
                    return -1;
                }
                for (int j = idx; j < n_ring - 1; j++) {
                    ring[j] = ring[j + 1];
                }
                n_ring--;
                pos = idx % (n_ring ? n_ring : 1);
                clipped = 1;
                break;
            }
        }
        if (!clipped) {
            break;
        }
    }

    free(pts);
    free(vmap);
    free(ring);
    return 0;
}

/* ---------------- NURBS Tensor-Product Surface Evaluator ---------------- */

static int find_knot_span(int n_ctrl, int deg, double t, const double *U) {
    int n = n_ctrl - 1;
    if (t >= U[n + 1]) {
        return n;
    }
    if (t <= U[deg]) {
        int k = deg;
        while (k < n && U[k + 1] == U[deg]) {
            k++;
        }
        return k;
    }
    int lo = deg, hi = n + 1;
    int mid = (lo + hi) / 2;
    while (t < U[mid] || t >= U[mid + 1]) {
        if (t < U[mid]) {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
        if (hi - lo <= 1) {
            break;
        }
    }
    return mid;
}

static void eval_basis(int span, double t, int deg, const double *U, double *N) {
    double *left = (double *)malloc(sizeof(double) * (size_t)(deg + 1));
    double *right = (double *)malloc(sizeof(double) * (size_t)(deg + 1));
    N[0] = 1.0;
    for (int j = 1; j <= deg; j++) {
        left[j] = t - U[span + 1 - j];
        right[j] = U[span + j] - t;
        double saved = 0.0;
        for (int r = 0; r < j; r++) {
            double denom = right[r + 1] + left[j - r];
            double temp = (fabs(denom) < 1e-15) ? 0.0 : (N[r] / denom);
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
    free(left);
    free(right);
}

static double *expand_knot_vector(const double *knots, const int *mults, int n_k, int *out_len) {
    int total = 0;
    for (int i = 0; i < n_k; i++) {
        total += mults[i];
    }
    double *U = (double *)malloc(sizeof(double) * (size_t)total);
    if (!U) {
        return NULL;
    }
    int p = 0;
    for (int i = 0; i < n_k; i++) {
        for (int m = 0; m < mults[i]; m++) {
            U[p++] = knots[i];
        }
    }
    *out_len = total;
    return U;
}

int bm_tessellate_nurbs_patch(
    int u_deg,
    int v_deg,
    int nu,
    int nv,
    const double *poles_xyz,
    const double *weights,
    const double *u_knots,
    const int *u_mults,
    int n_uk,
    const double *v_knots,
    const int *v_mults,
    int n_vk,
    int nu_samp,
    int nv_samp,
    bm_mesh_t *mesh
) {
    if (!poles_xyz || !u_knots || !u_mults || !v_knots || !v_mults || !mesh ||
        u_deg < 0 || v_deg < 0 || nu <= u_deg || nv <= v_deg || nu_samp < 1 || nv_samp < 1) {
        return -1;
    }

    int Ulen = 0, Vlen = 0;
    double *U = expand_knot_vector(u_knots, u_mults, n_uk, &Ulen);
    double *V = expand_knot_vector(v_knots, v_mults, n_vk, &Vlen);
    double *Nu = (double *)malloc(sizeof(double) * (size_t)(u_deg + 1));
    double *Nv = (double *)malloc(sizeof(double) * (size_t)(v_deg + 1));
    int *grid = (int *)malloc(sizeof(int) * (size_t)((nu_samp + 1) * (nv_samp + 1)));
    if (!U || !V || !Nu || !Nv || !grid || Ulen < nu + u_deg + 1 || Vlen < nv + v_deg + 1) {
        free(U); free(V); free(Nu); free(Nv); free(grid);
        return -1;
    }

    double u0 = U[u_deg], u1 = U[nu];
    double v0 = V[v_deg], v1 = V[nv];

    for (int iu = 0; iu <= nu_samp; iu++) {
        double u = u0 + (u1 - u0) * (double)iu / (double)nu_samp;
        int su = find_knot_span(nu, u_deg, u, U);
        eval_basis(su, u, u_deg, U, Nu);
        for (int iv = 0; iv <= nv_samp; iv++) {
            double v = v0 + (v1 - v0) * (double)iv / (double)nv_samp;
            int sv = find_knot_span(nv, v_deg, v, V);
            eval_basis(sv, v, v_deg, V, Nv);

            double wx = 0.0, wy = 0.0, wz = 0.0, wsum = 0.0;
            for (int r = 0; r <= u_deg; r++) {
                int ui = su - u_deg + r;
                double bu = Nu[r];
                for (int s = 0; s <= v_deg; s++) {
                    int vi = sv - v_deg + s;
                    int idx = ui * nv + vi;
                    double w = weights ? weights[idx] : 1.0;
                    double coeff = bu * Nv[s] * w;
                    wx += coeff * poles_xyz[3 * idx];
                    wy += coeff * poles_xyz[3 * idx + 1];
                    wz += coeff * poles_xyz[3 * idx + 2];
                    wsum += coeff;
                }
            }
            if (fabs(wsum) < 1e-15) {
                wsum = 1.0;
            }
            int vid = bm_mesh_add_vertex(mesh, wx / wsum, wy / wsum, wz / wsum);
            if (vid < 0) {
                free(U); free(V); free(Nu); free(Nv); free(grid);
                return -1;
            }
            grid[iu * (nv_samp + 1) + iv] = vid;
        }
    }

    int err = 0;
    for (int iu = 0; iu < nu_samp && !err; iu++) {
        for (int iv = 0; iv < nv_samp; iv++) {
            int a = grid[iu * (nv_samp + 1) + iv];
            int b = grid[(iu + 1) * (nv_samp + 1) + iv];
            int c = grid[(iu + 1) * (nv_samp + 1) + iv + 1];
            int d = grid[iu * (nv_samp + 1) + iv + 1];
            if (bm_mesh_add_triangle(mesh, a, b, c, 0) < 0 || bm_mesh_add_triangle(mesh, a, c, d, 0) < 0) {
                err = -1;
                break;
            }
        }
    }

    free(U);
    free(V);
    free(Nu);
    free(Nv);
    free(grid);
    return err;
}
