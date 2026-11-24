#define _CRT_SECURE_NO_WARNINGS
#include "brepmass.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void tri_normal(const double *a, const double *b, const double *c, float out[3]) {
    double ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2];
    double vx = c[0] - a[0], vy = c[1] - a[1], vz = c[2] - a[2];
    double nx = uy * vz - uz * vy;
    double ny = uz * vx - ux * vz;
    double nz = ux * vy - uy * vx;
    double len = sqrt(nx * nx + ny * ny + nz * nz);
    if (len > 1e-15) {
        out[0] = (float)(nx / len);
        out[1] = (float)(ny / len);
        out[2] = (float)(nz / len);
    } else {
        out[0] = out[1] = out[2] = 0.0f;
    }
}

int bm_write_stl_binary(const char *path, const bm_mesh_t *mesh) {
    if (!path || !mesh) {
        return -1;
    }
    FILE *fp = fopen(path, "wb");
    if (!fp) {
        return -1;
    }
    char hdr[80];
    memset(hdr, 0, sizeof(hdr));
    memcpy(hdr, "brepmass binary STL", 19);
    fwrite(hdr, 1, 80, fp);

    uint32_t nt = (uint32_t)mesh->n_triangles;
    fwrite(&nt, sizeof(uint32_t), 1, fp);

    for (int k = 0; k < mesh->n_triangles; k++) {
        const double *a = mesh->vertices + 3 * mesh->triangles[3 * k];
        const double *b = mesh->vertices + 3 * mesh->triangles[3 * k + 1];
        const double *c = mesh->vertices + 3 * mesh->triangles[3 * k + 2];
        float rec[12];
        tri_normal(a, b, c, rec);
        rec[3] = (float)a[0]; rec[4] = (float)a[1]; rec[5] = (float)a[2];
        rec[6] = (float)b[0]; rec[7] = (float)b[1]; rec[8] = (float)b[2];
        rec[9] = (float)c[0]; rec[10] = (float)c[1]; rec[11] = (float)c[2];
        uint16_t attr = 0;
        fwrite(rec, sizeof(float), 12, fp);
        fwrite(&attr, sizeof(uint16_t), 1, fp);
    }
    fclose(fp);
    return 0;
}

typedef struct {
    int64_t qx;
    int64_t qy;
    int64_t qz;
    int vid;
    int next;
} bm_vhash_node_t;

static uint64_t hash_q3(int64_t x, int64_t y, int64_t z) {
    uint64_t h = (uint64_t)x * 73856093ULL ^ (uint64_t)y * 19349663ULL ^ (uint64_t)z * 83492791ULL;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    return h;
}

int bm_read_stl_binary(const char *path, bm_mesh_t *mesh, double weld_tol) {
    if (!path || !mesh) {
        return -1;
    }
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        return -1;
    }
    char hdr[80];
    if (fread(hdr, 1, 80, fp) != 80) {
        fclose(fp);
        return -1;
    }
    uint32_t nt = 0;
    if (fread(&nt, sizeof(uint32_t), 1, fp) != 1) {
        fclose(fp);
        return -1;
    }

    double inv_tol = (weld_tol > 0.0) ? (1.0 / weld_tol) : 0.0;
    int nbuckets = (nt > 0) ? (int)(nt * 2 + 17) : 17;
    int *buckets = (int *)malloc(sizeof(int) * (size_t)nbuckets);
    bm_vhash_node_t *pool = (bm_vhash_node_t *)malloc(sizeof(bm_vhash_node_t) * (size_t)(nt * 3 + 1));
    if (!buckets || !pool) {
        free(buckets);
        free(pool);
        fclose(fp);
        return -1;
    }
    for (int i = 0; i < nbuckets; i++) {
        buckets[i] = -1;
    }
    int pool_n = 0;

    for (uint32_t k = 0; k < nt; k++) {
        float rec[12];
        uint16_t attr = 0;
        if (fread(rec, sizeof(float), 12, fp) != 12 || fread(&attr, sizeof(uint16_t), 1, fp) != 1) {
            free(buckets);
            free(pool);
            fclose(fp);
            return -1;
        }
        int vids[3];
        for (int vi = 0; vi < 3; vi++) {
            double px = (double)rec[3 + 3 * vi];
            double py = (double)rec[3 + 3 * vi + 1];
            double pz = (double)rec[3 + 3 * vi + 2];
            if (inv_tol <= 0.0) {
                vids[vi] = bm_mesh_add_vertex(mesh, px, py, pz);
            } else {
                int64_t qx = (int64_t)llround(px * inv_tol);
                int64_t qy = (int64_t)llround(py * inv_tol);
                int64_t qz = (int64_t)llround(pz * inv_tol);
                uint64_t b = hash_q3(qx, qy, qz) % (uint64_t)nbuckets;
                int found = -1;
                for (int cur = buckets[b]; cur >= 0; cur = pool[cur].next) {
                    if (pool[cur].qx == qx && pool[cur].qy == qy && pool[cur].qz == qz) {
                        found = pool[cur].vid;
                        break;
                    }
                }
                if (found < 0) {
                    found = bm_mesh_add_vertex(mesh, px, py, pz);
                    pool[pool_n].qx = qx;
                    pool[pool_n].qy = qy;
                    pool[pool_n].qz = qz;
                    pool[pool_n].vid = found;
                    pool[pool_n].next = buckets[b];
                    buckets[b] = pool_n++;
                }
                vids[vi] = found;
            }
        }
        if (vids[0] < 0 || vids[1] < 0 || vids[2] < 0 ||
            bm_mesh_add_triangle(mesh, vids[0], vids[1], vids[2], 0) < 0) {
            free(buckets);
            free(pool);
            fclose(fp);
            return -1;
        }
    }

    free(buckets);
    free(pool);
    fclose(fp);
    return 0;
}

int bm_write_stl_ascii(const char *path, const bm_mesh_t *mesh, const char *solid_name) {
    if (!path || !mesh) {
        return -1;
    }
    FILE *fp = fopen(path, "w");
    if (!fp) {
        return -1;
    }
    const char *sname = solid_name ? solid_name : "brepmass";
    fprintf(fp, "solid %s\n", sname);
    for (int k = 0; k < mesh->n_triangles; k++) {
        const double *a = mesh->vertices + 3 * mesh->triangles[3 * k];
        const double *b = mesh->vertices + 3 * mesh->triangles[3 * k + 1];
        const double *c = mesh->vertices + 3 * mesh->triangles[3 * k + 2];
        float n[3];
        tri_normal(a, b, c, n);
        fprintf(fp, "  facet normal %.9g %.9g %.9g\n", n[0], n[1], n[2]);
        fprintf(fp, "    outer loop\n");
        fprintf(fp, "      vertex %.15g %.15g %.15g\n", a[0], a[1], a[2]);
        fprintf(fp, "      vertex %.15g %.15g %.15g\n", b[0], b[1], b[2]);
        fprintf(fp, "      vertex %.15g %.15g %.15g\n", c[0], c[1], c[2]);
        fprintf(fp, "    endloop\n");
        fprintf(fp, "  endfacet\n");
    }
    fprintf(fp, "endsolid %s\n", sname);
    fclose(fp);
    return 0;
}

int bm_write_obj(const char *path, const bm_mesh_t *mesh) {
    if (!path || !mesh) {
        return -1;
    }
    FILE *fp = fopen(path, "w");
    if (!fp) {
        return -1;
    }
    fprintf(fp, "# brepmass Wavefront OBJ\n");
    for (int i = 0; i < mesh->n_vertices; i++) {
        fprintf(
            fp,
            "v %.15g %.15g %.15g\n",
            mesh->vertices[3 * i],
            mesh->vertices[3 * i + 1],
            mesh->vertices[3 * i + 2]
        );
    }
    for (int k = 0; k < mesh->n_triangles; k++) {
        fprintf(
            fp,
            "f %d %d %d\n",
            mesh->triangles[3 * k] + 1,
            mesh->triangles[3 * k + 1] + 1,
            mesh->triangles[3 * k + 2] + 1
        );
    }
    fclose(fp);
    return 0;
}
