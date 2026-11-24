#include "brepmass.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

static void check_rel(const char *label, double got, double want, double rel_tol) {
    double denom = (fabs(want) > 1e-12) ? fabs(want) : 1.0;
    double err = fabs(got - want) / denom;
    if (err > rel_tol) {
        fprintf(stderr, "FAIL %s: got=%.12g want=%.12g rel_err=%.3g > %.3g\n", label, got, want, err, rel_tol);
        assert(err <= rel_tol);
    }
}

static void check_abs(const char *label, double got, double want, double abs_tol) {
    double err = fabs(got - want);
    if (err > abs_tol) {
        fprintf(stderr, "FAIL %s: got=%.12g want=%.12g abs_err=%.3g > %.3g\n", label, got, want, err, abs_tol);
        assert(err <= abs_tol);
    }
}

static void test_box_exact(void) {
    bm_mesh_t m;
    bm_mesh_init(&m);
    assert(bm_tessellate_box(2.0, 3.0, 5.0, &m) == 0);

    bm_mass_props_t p, a;
    assert(bm_mesh_mass_props(m.vertices, m.n_vertices, m.triangles, m.n_triangles, 2.5, &p) == 0);
    assert(bm_analytic_box(2.0, 3.0, 5.0, 2.5, &a) == 0);

    assert(p.dim == 3);
    assert(p.is_closed == 1);
    assert(p.euler_char == 2);
    check_rel("box.vol", p.volume, a.volume, 1e-12);
    check_rel("box.area", p.surface_area, a.surface_area, 1e-12);
    check_rel("box.mass", p.mass, a.mass, 1e-12);
    for (int i = 0; i < 3; i++) {
        check_abs("box.cm", p.centroid[i], a.centroid[i], 1e-12);
        check_rel("box.pm", p.principal_moments[i], a.principal_moments[i], 1e-12);
        for (int j = 0; j < 3; j++) {
            check_abs("box.I0", p.inertia_origin[i][j], a.inertia_origin[i][j], 1e-10);
            check_abs("box.Icm", p.inertia_centroid[i][j], a.inertia_centroid[i][j], 1e-10);
        }
    }
    bm_mesh_free(&m);
    printf("  [PASS] test_box_exact\n");
}

static void test_cylinder_and_cone_and_sphere_and_torus(void) {
    bm_tess_opts_t opts = {192, 64, 0.0};
    bm_mesh_t m;
    bm_mass_props_t p, a;

    /* Cylinder r=0.8, h=2.5, density=7.8 */
    bm_mesh_init(&m);
    assert(bm_tessellate_cylinder(0.8, 2.5, &opts, &m) == 0);
    assert(bm_mesh_mass_props(m.vertices, m.n_vertices, m.triangles, m.n_triangles, 7.8, &p) == 0);
    assert(bm_analytic_cylinder(0.8, 2.5, 7.8, &a) == 0);
    assert(p.is_closed == 1 && p.euler_char == 2);
    check_rel("cyl.vol", p.volume, a.volume, 1e-3);
    check_rel("cyl.area", p.surface_area, a.surface_area, 1e-3);
    check_abs("cyl.zc", p.centroid[2], a.centroid[2], 1e-4);
    check_rel("cyl.Icm_zz", p.inertia_centroid[2][2], a.inertia_centroid[2][2], 2e-3);
    check_rel("cyl.Icm_xx", p.inertia_centroid[0][0], a.inertia_centroid[0][0], 2e-3);
    bm_mesh_free(&m);

    /* Cone frustum r1=1.2, r2=0.4, h=3.0, density=2.7 */
    bm_mesh_init(&m);
    assert(bm_tessellate_cone_frustum(1.2, 0.4, 3.0, &opts, &m) == 0);
    assert(bm_mesh_mass_props(m.vertices, m.n_vertices, m.triangles, m.n_triangles, 2.7, &p) == 0);
    assert(bm_analytic_cone_frustum(1.2, 0.4, 3.0, 2.7, &a) == 0);
    assert(p.is_closed == 1 && p.euler_char == 2);
    check_rel("cone.vol", p.volume, a.volume, 1e-3);
    check_rel("cone.area", p.surface_area, a.surface_area, 1e-3);
    check_abs("cone.zc", p.centroid[2], a.centroid[2], 1e-3);
    check_rel("cone.Icm_zz", p.inertia_centroid[2][2], a.inertia_centroid[2][2], 2e-3);
    check_rel("cone.Icm_xx", p.inertia_centroid[0][0], a.inertia_centroid[0][0], 2e-3);
    bm_mesh_free(&m);

    /* Sphere r=1.5, density=1.2 */
    bm_mesh_init(&m);
    assert(bm_tessellate_sphere(1.5, &opts, &m) == 0);
    assert(bm_mesh_mass_props(m.vertices, m.n_vertices, m.triangles, m.n_triangles, 1.2, &p) == 0);
    assert(bm_analytic_sphere(1.5, 1.2, &a) == 0);
    assert(p.is_closed == 1 && p.euler_char == 2);
    check_rel("sph.vol", p.volume, a.volume, 2e-3);
    check_rel("sph.area", p.surface_area, a.surface_area, 2e-3);
    check_rel("sph.Icm", p.principal_moments[0], a.principal_moments[0], 3e-3);
    bm_mesh_free(&m);

    /* Torus R=2.0, r=0.5, density=1.5 (genus 1 -> euler_char == 0) */
    bm_mesh_init(&m);
    assert(bm_tessellate_torus(2.0, 0.5, &opts, &m) == 0);
    assert(bm_mesh_mass_props(m.vertices, m.n_vertices, m.triangles, m.n_triangles, 1.5, &p) == 0);
    assert(bm_analytic_torus(2.0, 0.5, 1.5, &a) == 0);
    assert(p.is_closed == 1 && p.euler_char == 0);
    check_rel("tor.vol", p.volume, a.volume, 3e-3);
    check_rel("tor.area", p.surface_area, a.surface_area, 3e-3);
    check_rel("tor.Izz", p.inertia_centroid[2][2], a.inertia_centroid[2][2], 5e-3);
    check_rel("tor.Ixx", p.inertia_centroid[0][0], a.inertia_centroid[0][0], 5e-3);
    bm_mesh_free(&m);

    printf("  [PASS] test_cylinder_and_cone_and_sphere_and_torus\n");
}

static void test_polygon_with_hole_and_nurbs_and_stl_roundtrip(void) {
    /* 1. 4x4 square plate with 1x1 square hole -> area = 15, open sheet */
    bm_mesh_t m;
    bm_mesh_init(&m);
    const double outer_xy[8] = {0.0, 0.0, 4.0, 0.0, 4.0, 4.0, 0.0, 4.0};
    const double hole_xy[8] = {1.5, 1.5, 1.5, 2.5, 2.5, 2.5, 2.5, 1.5};
    const int hole_cnt[1] = {4};
    const double orig[3] = {0.0, 0.0, 2.0};
    const double norm[3] = {0.0, 0.0, 1.0};
    const double xax[3] = {1.0, 0.0, 0.0};
    assert(bm_tessellate_polygon_2d(outer_xy, 4, hole_xy, hole_cnt, 1, orig, norm, xax, &m) == 0);

    bm_mass_props_t p;
    assert(bm_mesh_mass_props(m.vertices, m.n_vertices, m.triangles, m.n_triangles, 1.0, &p) == 0);
    assert(p.dim == 2 && p.is_closed == 0);
    check_rel("poly_hole.area", p.surface_area, 15.0, 1e-12);
    check_abs("poly_hole.cm_x", p.centroid[0], 2.0, 1e-12);
    check_abs("poly_hole.cm_y", p.centroid[1], 2.0, 1e-12);
    check_abs("poly_hole.cm_z", p.centroid[2], 2.0, 1e-12);
    bm_mesh_free(&m);

    /* 2. Rational NURBS quarter-cylinder (r=1, h=1) -> exact area pi/2 */
    const double w = sqrt(0.5);
    const double poles[18] = {
        1.0, 0.0, 0.0,  1.0, 0.0, 1.0,
        1.0, 1.0, 0.0,  1.0, 1.0, 1.0,
        0.0, 1.0, 0.0,  0.0, 1.0, 1.0,
    };
    const double weights[6] = {1.0, 1.0, w, w, 1.0, 1.0};
    const double uk[2] = {0.0, 1.0}, vk[2] = {0.0, 1.0};
    const int um[2] = {3, 3}, vm[2] = {2, 2};
    bm_mesh_init(&m);
    assert(bm_tessellate_nurbs_patch(2, 1, 3, 2, poles, weights, uk, um, 2, vk, vm, 2, 64, 4, &m) == 0);
    assert(bm_mesh_mass_props(m.vertices, m.n_vertices, m.triangles, m.n_triangles, 1.0, &p) == 0);
    check_rel("nurbs.area", p.surface_area, M_PI / 2.0, 1e-3);
    bm_mesh_free(&m);

    /* 3. Binary STL write + weld read roundtrip retains closed watertightness & exact volume */
    bm_mesh_init(&m);
    assert(bm_tessellate_box(1.5, 2.5, 3.5, &m) == 0);
    const char *stl_tmp = "build_test_tmp.stl";
    assert(bm_write_stl_binary(stl_tmp, &m) == 0);
    bm_mesh_free(&m);

    bm_mesh_t m2;
    bm_mesh_init(&m2);
    assert(bm_read_stl_binary(stl_tmp, &m2, 1e-6) == 0);
    remove(stl_tmp);
    assert(bm_mesh_mass_props(m2.vertices, m2.n_vertices, m2.triangles, m2.n_triangles, 1.0, &p) == 0);
    assert(p.is_closed == 1 && p.euler_char == 2 && p.n_vertices == 8 && p.n_triangles == 12);
    check_rel("stl_rt.vol", p.volume, 1.5 * 2.5 * 3.5, 1e-5);
    bm_mesh_free(&m2);

    printf("  [PASS] test_polygon_with_hole_and_nurbs_and_stl_roundtrip\n");
}

int main(void) {
    printf("Running brepmass C99 unit tests...\n");
    test_box_exact();
    test_cylinder_and_cone_and_sphere_and_torus();
    test_polygon_with_hole_and_nurbs_and_stl_roundtrip();
    printf("All brepmass C99 tests PASSED.\n");
    return 0;
}
