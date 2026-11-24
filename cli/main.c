#include "brepmass.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_props(const bm_mass_props_t *p) {
    printf("Closed:           %s (dim=%d, euler=%d)\n", p->is_closed ? "yes" : "no", p->dim, p->euler_char);
    printf("Vertices / Tris:  %d / %d\n", p->n_vertices, p->n_triangles);
    printf("Volume:           %.12g\n", p->volume);
    printf("Surface area:     %.12g\n", p->surface_area);
    printf("Mass:             %.12g\n", p->mass);
    printf("Centroid:         (%.9g, %.9g, %.9g)\n", p->centroid[0], p->centroid[1], p->centroid[2]);
    printf(
        "Principal mom.:   (%.9g, %.9g, %.9g)\n",
        p->principal_moments[0], p->principal_moments[1], p->principal_moments[2]
    );
    printf(
        "BBox:             (%.6g, %.6g, %.6g) .. (%.6g, %.6g, %.6g)\n",
        p->bbox_min[0], p->bbox_min[1], p->bbox_min[2],
        p->bbox_max[0], p->bbox_max[1], p->bbox_max[2]
    );
}

static int ends_with(const char *s, const char *suf) {
    size_t ls = strlen(s), lq = strlen(suf);
    return ls >= lq && strcmp(s + ls - lq, suf) == 0;
}

static void usage(void) {
    fprintf(
        stderr,
        "Usage:\n"
        "  brepmass info <file.stl> [--density D]\n"
        "  brepmass make box <dx> <dy> <dz> -o <out.stl|out.obj>\n"
        "  brepmass make cylinder <r> <h> -o <out.stl|out.obj>\n"
        "  brepmass make cone <r1> <r2> <h> -o <out.stl|out.obj>\n"
        "  brepmass make sphere <r> -o <out.stl|out.obj>\n"
        "  brepmass make torus <major_r> <minor_r> -o <out.stl|out.obj>\n"
    );
}

int main(int argc, char **argv) {
    if (argc < 3) {
        usage();
        return 1;
    }

    if (strcmp(argv[1], "info") == 0) {
        const char *path = argv[2];
        double density = 1.0;
        for (int i = 3; i + 1 < argc; i++) {
            if (strcmp(argv[i], "--density") == 0) {
                density = atof(argv[i + 1]);
            }
        }
        bm_mesh_t m;
        bm_mesh_init(&m);
        if (bm_read_stl_binary(path, &m, 1e-6) != 0) {
            fprintf(stderr, "Failed to read binary STL: %s\n", path);
            bm_mesh_free(&m);
            return 1;
        }
        bm_mass_props_t props;
        if (bm_mesh_mass_props(m.vertices, m.n_vertices, m.triangles, m.n_triangles, density, &props) != 0) {
            fprintf(stderr, "Failed to evaluate mass properties.\n");
            bm_mesh_free(&m);
            return 1;
        }
        print_props(&props);
        bm_mesh_free(&m);
        return 0;
    }

    if (strcmp(argv[1], "make") == 0) {
        const char *kind = argv[2];
        const char *out_path = NULL;
        for (int i = 3; i + 1 < argc; i++) {
            if (strcmp(argv[i], "-o") == 0) {
                out_path = argv[i + 1];
            }
        }
        bm_mesh_t m;
        bm_mesh_init(&m);
        bm_tess_opts_t opts = bm_tess_opts_default();
        int rc = -1;
        if (strcmp(kind, "box") == 0 && argc >= 6) {
            rc = bm_tessellate_box(atof(argv[3]), atof(argv[4]), atof(argv[5]), &m);
        } else if (strcmp(kind, "cylinder") == 0 && argc >= 5) {
            rc = bm_tessellate_cylinder(atof(argv[3]), atof(argv[4]), &opts, &m);
        } else if (strcmp(kind, "cone") == 0 && argc >= 6) {
            rc = bm_tessellate_cone_frustum(atof(argv[3]), atof(argv[4]), atof(argv[5]), &opts, &m);
        } else if (strcmp(kind, "sphere") == 0 && argc >= 4) {
            rc = bm_tessellate_sphere(atof(argv[3]), &opts, &m);
        } else if (strcmp(kind, "torus") == 0 && argc >= 5) {
            rc = bm_tessellate_torus(atof(argv[3]), atof(argv[4]), &opts, &m);
        } else {
            usage();
            bm_mesh_free(&m);
            return 1;
        }
        if (rc != 0) {
            fprintf(stderr, "Tessellation failed for %s.\n", kind);
            bm_mesh_free(&m);
            return 1;
        }
        if (out_path) {
            int wrc = ends_with(out_path, ".obj")
                          ? bm_write_obj(out_path, &m)
                          : bm_write_stl_binary(out_path, &m);
            if (wrc != 0) {
                fprintf(stderr, "Failed writing %s\n", out_path);
                bm_mesh_free(&m);
                return 1;
            }
        }
        bm_mass_props_t props;
        if (bm_mesh_mass_props(m.vertices, m.n_vertices, m.triangles, m.n_triangles, 1.0, &props) == 0) {
            print_props(&props);
        }
        bm_mesh_free(&m);
        return 0;
    }

    usage();
    return 1;
}
