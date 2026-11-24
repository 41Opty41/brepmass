/**
 * @file brepmass.h
 * @brief Lightweight C99 library for B-rep mesh mass properties, inertia tensors,
 *        analytic primitive evaluation, NURBS/polygon tessellation, and STL/OBJ I/O.
 */

#ifndef BREPMASS_H
#define BREPMASS_H

#include <stddef.h>

#if defined(_WIN32) && defined(BM_BUILDING_DLL)
#  define BM_API __declspec(dllexport)
#elif defined(_WIN32) && defined(BM_USE_DLL)
#  define BM_API __declspec(dllimport)
#else
#  define BM_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Complete geometric, inertial, and topological properties of a solid or sheet. */
typedef struct bm_mass_props_t {
    int dim;                      /**< 3 = closed solid volume, 2 = open sheet surface */
    int is_closed;                /**< 1 if watertight manifold half-edge pairing, 0 otherwise */
    int euler_char;               /**< V - E + F across referenced vertices, edges, and triangles */
    int n_vertices;               /**< Number of vertices in input mesh */
    int n_triangles;              /**< Number of triangles in input mesh */
    double volume;                /**< Enclosed volume (0 for open sheets or dim=2) */
    double surface_area;          /**< Total triangle surface area */
    double mass;                  /**< density * (dim == 3 ? volume : surface_area) */
    double centroid[3];           /**< Center of mass / area (x, y, z) */
    double inertia_origin[3][3];  /**< Mass moment of inertia tensor about (0, 0, 0) */
    double inertia_centroid[3][3];/**< Mass moment of inertia tensor about centroid */
    double principal_moments[3];  /**< Ascending principal moments of inertia I0 <= I1 <= I2 */
    double principal_axes[3][3];  /**< Orthonormal eigenvectors; column j matches principal_moments[j] */
    double bbox_min[3];           /**< Axis-aligned bounding box minimum */
    double bbox_max[3];           /**< Axis-aligned bounding box maximum */
} bm_mass_props_t;

/** Indexed triangle mesh container. */
typedef struct bm_mesh_t {
    double *vertices;             /**< Flat XYZ array of length 3 * n_vertices */
    int *triangles;               /**< Flat index array of length 3 * n_triangles */
    int *face_ids;                /**< Optional per-triangle face id (length n_triangles), or NULL */
    int n_vertices;
    int cap_vertices;
    int n_triangles;
    int cap_triangles;
} bm_mesh_t;

/** Tessellation subdivision options. */
typedef struct bm_tess_opts_t {
    int u_segments;               /**< Circumferential / u-direction divisions (default 64) */
    int v_segments;               /**< Axial / v-direction divisions (default 32) */
    double chord_tol;             /**< Reserved chordal tolerance hint (0 = fixed segments) */
} bm_tess_opts_t;

/* ---------------- Mesh lifecycle ---------------- */

BM_API void bm_mesh_init(bm_mesh_t *mesh);
BM_API void bm_mesh_free(bm_mesh_t *mesh);
BM_API int bm_mesh_reserve(bm_mesh_t *mesh, int n_vertices, int n_triangles);
BM_API int bm_mesh_add_vertex(bm_mesh_t *mesh, double x, double y, double z);
BM_API int bm_mesh_add_triangle(bm_mesh_t *mesh, int i0, int i1, int i2, int face_id);

/* ---------------- Mass, Inertia & Topology ---------------- */

/**
 * Check watertightness: every directed edge (u -> v) must have a unique opposite
 * half-edge (v -> u). Optionally outputs Euler characteristic V_used - E + F.
 * Returns 1 if closed watertight manifold, 0 if open/non-manifold, -1 on error.
 */
BM_API int bm_mesh_is_closed(
    const int *triangles,
    int n_triangles,
    int n_vertices,
    int *out_euler
);

/**
 * Compute exact divergence-theorem volume, surface area, centroid, inertia tensors
 * about origin and centroid, principal moments/axes, bbox, closure, and Euler.
 * For open meshes, computes surface-area integrals (dim=2, volume=0).
 * Returns 0 on success, non-zero on invalid arguments.
 */
BM_API int bm_mesh_mass_props(
    const double *vertices,
    int n_vertices,
    const int *triangles,
    int n_triangles,
    double density,
    bm_mass_props_t *out
);

/* Closed-form exact analytic formulas for canonical solids centered at origin: */
BM_API int bm_analytic_box(double dx, double dy, double dz, double density, bm_mass_props_t *out);
BM_API int bm_analytic_cylinder(double radius, double height, double density, bm_mass_props_t *out);
BM_API int bm_analytic_cone_frustum(double r1, double r2, double height, double density, bm_mass_props_t *out);
BM_API int bm_analytic_sphere(double radius, double density, bm_mass_props_t *out);
BM_API int bm_analytic_torus(double major_r, double minor_r, double density, bm_mass_props_t *out);

/* ---------------- Tessellation ---------------- */

BM_API bm_tess_opts_t bm_tess_opts_default(void);

BM_API int bm_tessellate_box(double dx, double dy, double dz, bm_mesh_t *mesh);
BM_API int bm_tessellate_cylinder(double radius, double height, const bm_tess_opts_t *opts, bm_mesh_t *mesh);
BM_API int bm_tessellate_cone_frustum(double r1, double r2, double height, const bm_tess_opts_t *opts, bm_mesh_t *mesh);
BM_API int bm_tessellate_sphere(double radius, const bm_tess_opts_t *opts, bm_mesh_t *mesh);
BM_API int bm_tessellate_torus(double major_r, double minor_r, const bm_tess_opts_t *opts, bm_mesh_t *mesh);

/**
 * Triangulate a planar 2D polygon with optional inner holes by ear clipping and
 * lift into 3D space using (origin, normal, x_axis).
 * outer_xy has 2*n_outer doubles; holes_xy concatenates 2*hole_counts[h] doubles.
 */
BM_API int bm_tessellate_polygon_2d(
    const double *outer_xy,
    int n_outer,
    const double *holes_xy,
    const int *hole_counts,
    int n_holes,
    const double origin[3],
    const double normal[3],
    const double x_axis[3],
    bm_mesh_t *mesh
);

/**
 * Sample a rational or non-rational tensor-product B-spline / NURBS surface
 * onto a triangle grid of nu_samp x nv_samp points.
 * poles_xyz is flat 3*(nu*nv); weights is nu*nv (or NULL for non-rational).
 */
BM_API int bm_tessellate_nurbs_patch(
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
);

/* ---------------- Mesh File I/O ---------------- */

BM_API int bm_write_stl_binary(const char *path, const bm_mesh_t *mesh);
BM_API int bm_read_stl_binary(const char *path, bm_mesh_t *mesh, double weld_tol);
BM_API int bm_write_stl_ascii(const char *path, const bm_mesh_t *mesh, const char *solid_name);
BM_API int bm_write_obj(const char *path, const bm_mesh_t *mesh);

#ifdef __cplusplus
}
#endif

#endif /* BREPMASS_H */
