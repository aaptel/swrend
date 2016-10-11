// shitty and slow 3d software renderer
// Usage:
// $ gcc -std=c11 -lm swrend.c -o swrend
// $ ./swrend
// $ feh out.ppm # open the pic
//
// Copyright (c) 2016, Aurélien Aptel <aurelien.aptel@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>
#include <float.h>

// basic utils

static int g_debug = 0;

#define E(fmt, ...) \
    fprintf(stderr,"[E] %s:%d "fmt"\n",__func__,__LINE__,##__VA_ARGS__),exit(1)
#define W(fmt, ...) \
    fprintf(stderr,"[W] %s:%d "fmt"\n",__func__,__LINE__,## __VA_ARGS__)
#define DBG(fmt, ...) \
    if(g_debug) fprintf(stderr,"[*] %s:%d "fmt"\n",__func__,__LINE__,## __VA_ARGS__)

#define CLAMP(x, min, max) ((x) < (min) ? (min) : (x) > (max) ? (max) : (x))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) < (b) ? (b) : (a))
#define SWAP(a, b, tmp) (tmp = a, a = b, b = tmp)

void* xcalloc (size_t a, size_t b)
{
    void* p = calloc(a, b);
    if (!p)
        E("out of memory");
    return p;
}

void* xmalloc (size_t a)
{
    void* p = malloc(a);
    if (!p)
        E("out of memory");
    return p;
}

// img type and operations

typedef struct {
    uint32_t* buf;
    size_t w, h;
} img_t;

typedef struct {
    float* buf;
    size_t w, h;
} imgf_t;

#define IMG_P(p, x, y) ((p)->buf[ (y)*(p)->w  + (x)])
#define IMG_RGB(r, g, b) ((((r) & 0xff) << 16) + (((g) & 0xff) << 8) + (((b) & 0xff) << 0))
#define IMG_A(c) (((c) >> 24) & 0xff)
#define IMG_R(c) (((c) >> 16) & 0xff)
#define IMG_G(c) (((c) >>  8) & 0xff)
#define IMG_B(c) (((c) >>  0) & 0xff)
#define PTR_TO_RGB(p) ((uint32_t)(((intptr_t)(p)) & 0xffffff))

// safe set pixel (clips on image borders)
void img_set_p (img_t* p, uint32_t c, int x, int y)
{
    if (0 <= x&&x < p->w && 0 <= y&&y < p->h)
        IMG_P(p, x, y) = c;
}

// safe set pixel (clips on image borders)
void imgf_set_p (imgf_t* p, float v, int x, int y)
{
    if (0 <= x&&x < p->w && 0 <= y&&y < p->h)
        IMG_P(p, x, y) = v;
}

img_t* img_new (size_t w, size_t h)
{
    img_t *p = xcalloc(1, sizeof(*p));
    p->buf = xcalloc(w*h, sizeof(*p->buf));
    p->w = w;
    p->h = h;
    return p;
}

imgf_t* imgf_new (size_t w, size_t h)
{
    imgf_t *p = xcalloc(1, sizeof(*p));
    p->buf = xcalloc(w*h, sizeof(*p->buf));
    p->w = w;
    p->h = h;
    return p;
}

img_t* img_from_imgf (imgf_t* pf)
{
    img_t* p = img_new(pf->w, pf->h);
    size_t max = pf->w * pf->h;
    for (size_t i = 0; i < max; i++) {
        uint32_t c = (CLAMP(pf->buf[i], 0.0f, 1.0f)*255.0f);
        p->buf[i] = c|(c<<8)|(c<<16);
    }
    return p;
}

int img_to_ppm (const img_t* p, const char* fn)
{
    FILE* fh = fopen(fn, "w+");
    if (!fh) {
        W("fopen %s", strerror(errno));
        return -1;
    }

    fprintf(fh, "P3\n%zu %zu\n255\n", p->w, p->h);

    for (size_t y = 0; y < p->h; y++) {
        for (size_t x = 0; x < p->w; x++) {
            uint32_t c = IMG_P(p, x, y);
            fprintf(fh, "%d %d %d ", IMG_R(c), IMG_G(c), IMG_B(c));
        }
        fprintf(fh, "\n");
    }
    if (fclose(fh) != 0) {
        W("fclose %s", strerror(errno));
        return -1;
    }
    return 0;
}

void img_fill (img_t* p, uint32_t c)
{
    for (size_t i = 0; i < p->w*p->h; i++)
        p->buf[i] = c;
}

void img_line (img_t* p, uint32_t c, int x0, int y0, int x1, int y1)
{
    // tiny but slow ray-marching

    float x = x0, y = y0;
    float dx = x1-x0;
    float dy = y1-y0;
    int len = ceilf(sqrtf(dx*dx + dy*dy));

    dx /= len;
    dy /= len;

    for (int i = 0; i < len; i++) {
        img_set_p(p, c, roundf(x), roundf(y));
        x += dx;
        y += dy;
    }
}

// 4d vectors and matrices

#define PV3(v) printf("%f %f %f\n", (v).x, (v).y, (v).z)
#define PV2(v) printf("%f %f\n", (v).x, (v).y)

typedef union {
    float v[3];
    struct { float x, y, z; };
    struct { float r, g, b; };
} vec3;
typedef union {
    float v[4];
    struct { float x, y, z, w; };
} vec4;
typedef union {
    float v[2];
    struct { float x, y; };
} vec2;
typedef struct {
    float v[16];
} mat4;

void vec3_from_v4 (vec3* dst, vec4* src)
{
    memcpy(dst, src, sizeof(vec3));
}

void vec4_from_v3 (vec4* dst, vec3* src, float w)
{
    memcpy(dst, src, sizeof(vec3));
    dst->w = w;
}

float lerp (float a, float b, float t)
{
    return (1.0f-t)*a + t*b;
}

float lerp3 (float a, float b, float c, float ta, float tb, float tc)
{
    return a*ta + b*tb + c*tc;
}


float vec3_dot (const vec3* a, const vec3* b)
{
    return a->x*b->x + a->y*b->y + a->z*b->z;
}

void vec3_add (vec3* r, const vec3* a, const vec3* b)
{
    r->v[0] = a->v[0] + b->v[0];
    r->v[1] = a->v[1] + b->v[1];
    r->v[2] = a->v[2] + b->v[2];
}

void vec3_sub (vec3* r, const vec3* a, const vec3* b)
{
    r->v[0] = a->v[0] - b->v[0];
    r->v[1] = a->v[1] - b->v[1];
    r->v[2] = a->v[2] - b->v[2];
}

float vec3_len (const vec3* v)
{
    return sqrtf(vec3_dot(v,v));
}

float vec3_len2 (const vec3* v)
{
    return vec3_dot(v,v);
}

void vec3_init (vec3* v, float x, float y, float z, float w)
{
    v->x=x, v->y=y, v->z=z;
}

void vec3_mul (vec3* v, float a)
{
    v->x *= a;
    v->y *= a;
    v->z *= a;
}

void vec3_mul_v (vec3* v, const vec3* a)
{
    v->x *= a->x;
    v->y *= a->y;
    v->z *= a->z;
}

void vec3_normalize (vec3* v)
{
    vec3_mul(v, 1.f/vec3_len(v));
}


void vec3_cross (vec3* r, const vec3* a, const vec3* b)
{
    r->x = a->y*b->z - a->z*b->y;
    r->y = a->z*b->x - a->x*b->z;
    r->z = a->x*b->y - a->y*b->x;
}

int mat4_inv(mat4* r, const mat4* m)
{
    float inv[16], det;

    inv[0] = m->v[5]  * m->v[10] * m->v[15] -
             m->v[5]  * m->v[11] * m->v[14] -
             m->v[9]  * m->v[6]  * m->v[15] +
             m->v[9]  * m->v[7]  * m->v[14] +
             m->v[13] * m->v[6]  * m->v[11] -
             m->v[13] * m->v[7]  * m->v[10];

    inv[4] = -m->v[4]  * m->v[10] * m->v[15] +
              m->v[4]  * m->v[11] * m->v[14] +
              m->v[8]  * m->v[6]  * m->v[15] -
              m->v[8]  * m->v[7]  * m->v[14] -
              m->v[12] * m->v[6]  * m->v[11] +
              m->v[12] * m->v[7]  * m->v[10];

    inv[8] = m->v[4]  * m->v[9] * m->v[15] -
             m->v[4]  * m->v[11] * m->v[13] -
             m->v[8]  * m->v[5] * m->v[15] +
             m->v[8]  * m->v[7] * m->v[13] +
             m->v[12] * m->v[5] * m->v[11] -
             m->v[12] * m->v[7] * m->v[9];

    inv[12] = -m->v[4]  * m->v[9] * m->v[14] +
               m->v[4]  * m->v[10] * m->v[13] +
               m->v[8]  * m->v[5] * m->v[14] -
               m->v[8]  * m->v[6] * m->v[13] -
               m->v[12] * m->v[5] * m->v[10] +
               m->v[12] * m->v[6] * m->v[9];

    inv[1] = -m->v[1]  * m->v[10] * m->v[15] +
              m->v[1]  * m->v[11] * m->v[14] +
              m->v[9]  * m->v[2] * m->v[15] -
              m->v[9]  * m->v[3] * m->v[14] -
              m->v[13] * m->v[2] * m->v[11] +
              m->v[13] * m->v[3] * m->v[10];

    inv[5] = m->v[0]  * m->v[10] * m->v[15] -
             m->v[0]  * m->v[11] * m->v[14] -
             m->v[8]  * m->v[2] * m->v[15] +
             m->v[8]  * m->v[3] * m->v[14] +
             m->v[12] * m->v[2] * m->v[11] -
             m->v[12] * m->v[3] * m->v[10];

    inv[9] = -m->v[0]  * m->v[9] * m->v[15] +
              m->v[0]  * m->v[11] * m->v[13] +
              m->v[8]  * m->v[1] * m->v[15] -
              m->v[8]  * m->v[3] * m->v[13] -
              m->v[12] * m->v[1] * m->v[11] +
              m->v[12] * m->v[3] * m->v[9];

    inv[13] = m->v[0]  * m->v[9] * m->v[14] -
              m->v[0]  * m->v[10] * m->v[13] -
              m->v[8]  * m->v[1] * m->v[14] +
              m->v[8]  * m->v[2] * m->v[13] +
              m->v[12] * m->v[1] * m->v[10] -
              m->v[12] * m->v[2] * m->v[9];

    inv[2] = m->v[1]  * m->v[6] * m->v[15] -
             m->v[1]  * m->v[7] * m->v[14] -
             m->v[5]  * m->v[2] * m->v[15] +
             m->v[5]  * m->v[3] * m->v[14] +
             m->v[13] * m->v[2] * m->v[7] -
             m->v[13] * m->v[3] * m->v[6];

    inv[6] = -m->v[0]  * m->v[6] * m->v[15] +
              m->v[0]  * m->v[7] * m->v[14] +
              m->v[4]  * m->v[2] * m->v[15] -
              m->v[4]  * m->v[3] * m->v[14] -
              m->v[12] * m->v[2] * m->v[7] +
              m->v[12] * m->v[3] * m->v[6];

    inv[10] = m->v[0]  * m->v[5] * m->v[15] -
              m->v[0]  * m->v[7] * m->v[13] -
              m->v[4]  * m->v[1] * m->v[15] +
              m->v[4]  * m->v[3] * m->v[13] +
              m->v[12] * m->v[1] * m->v[7] -
              m->v[12] * m->v[3] * m->v[5];

    inv[14] = -m->v[0]  * m->v[5] * m->v[14] +
               m->v[0]  * m->v[6] * m->v[13] +
               m->v[4]  * m->v[1] * m->v[14] -
               m->v[4]  * m->v[2] * m->v[13] -
               m->v[12] * m->v[1] * m->v[6] +
               m->v[12] * m->v[2] * m->v[5];

    inv[3] = -m->v[1] * m->v[6] * m->v[11] +
              m->v[1] * m->v[7] * m->v[10] +
              m->v[5] * m->v[2] * m->v[11] -
              m->v[5] * m->v[3] * m->v[10] -
              m->v[9] * m->v[2] * m->v[7] +
              m->v[9] * m->v[3] * m->v[6];

    inv[7] = m->v[0] * m->v[6] * m->v[11] -
             m->v[0] * m->v[7] * m->v[10] -
             m->v[4] * m->v[2] * m->v[11] +
             m->v[4] * m->v[3] * m->v[10] +
             m->v[8] * m->v[2] * m->v[7] -
             m->v[8] * m->v[3] * m->v[6];

    inv[11] = -m->v[0] * m->v[5] * m->v[11] +
               m->v[0] * m->v[7] * m->v[9] +
               m->v[4] * m->v[1] * m->v[11] -
               m->v[4] * m->v[3] * m->v[9] -
               m->v[8] * m->v[1] * m->v[7] +
               m->v[8] * m->v[3] * m->v[5];

    inv[15] = m->v[0] * m->v[5] * m->v[10] -
              m->v[0] * m->v[6] * m->v[9] -
              m->v[4] * m->v[1] * m->v[10] +
              m->v[4] * m->v[2] * m->v[9] +
              m->v[8] * m->v[1] * m->v[6] -
              m->v[8] * m->v[2] * m->v[5];

    det = m->v[0] * inv[0] + m->v[1] * inv[4] + m->v[2] * inv[8] + m->v[3] * inv[12];

    if (det == 0)
        return 0;

    for (int i = 0; i < 16; i++)
        r->v[i] = inv[i] / det;

    return det;
}

void mat4_id (mat4* r)
{
    memset(r, 0, sizeof(*r));
    r->v[0] = r->v[5] = r->v[10] = r->v[15] = 1.0f;
}

void mat4_scale (mat4* r, float x, float y, float z)
{
    mat4_id(r);
    r->v[0]  = x;
    r->v[5]  = y;
    r->v[10] = z;
}

void mat4_translate (mat4* r, float x, float y, float z)
{
    mat4_id(r);
    r->v[12] = x;
    r->v[13] = y;
    r->v[14] = z;
}

void mat4_translate_v (mat4* r, const vec3* v)
{
    mat4_id(r);
    r->v[12] = v->x;
    r->v[13] = v->y;
    r->v[14] = v->z;
}

void mat4_mul (mat4* r, const mat4* a, const mat4* b)
{
    for (int col = 0; col < 4; col++) {
        for (int row = 0; row < 4; row++) {
            float sum = 0;
            for (int i = 0; i < 4; i++)
                sum += a->v[i*4+col] * b->v[row*4+i];
            r->v[row*4+col] = sum;
        }
    }
}

void mat4_mul_v3 (vec3* r, const mat4* m, const vec3* v, int point)
{
    for (int col = 0; col < 3; col++) {
        float sum = 0;
        for (int i = 0; i < 3; i++)
            sum += v->v[i] * m->v[i*4 + col];
        r->v[col] = sum;
    }
    if (point) {
        r->v[0] += m->v[12];
        r->v[1] += m->v[13];
        r->v[2] += m->v[14];
    }
}

void mat4_mul_v4 (vec4* r, const mat4* m, const vec4* v)
{
    for (int col = 0; col < 4; col++) {
        float sum = 0;
        for (int i = 0; i < 4; i++)
            sum += v->v[i] * m->v[i*4 + col];
        r->v[col] = sum;
    }
}

void mat4_print (const mat4* m)
{
    printf("%+.3f %+.3f %+.3f %+.3f\n"
           "%+.3f %+.3f %+.3f %+.3f\n"
           "%+.3f %+.3f %+.3f %+.3f\n"
           "%+.3f %+.3f %+.3f %+.3f\n",
           m->v[0], m->v[4], m->v[8 ], m->v[12],
           m->v[1], m->v[5], m->v[9 ], m->v[13],
           m->v[2], m->v[6], m->v[10], m->v[14],
           m->v[3], m->v[7], m->v[11], m->v[15]);
}

// obj loader

typedef struct {
    size_t nb_v, nb_n, nb_f, nb_uv;
    vec3 *v, *n;
    vec2* uv;
    struct {size_t v[3], n[3], t[3];} *f;
} obj_t;

typedef struct {
    vec4 v;  // coord world space
    vec2 sv; // coord screen space
    vec3 n; // normal
    vec2 t; // [INTERPOLATED] texture coordinates
    vec3 c; // [INTERPOLATED] color
    vec3 bc; // barycentric coefficent for interpolation
} vert_attr_t;

obj_t* obj_load (const char* fn)
{
    FILE* fh = fopen(fn, "r");
    if (!fh)
        E("fopen %s", strerror(errno));

    char buf[512];
    char id[32];
    float a, b, c;
    size_t nb_v = 0, nb_n = 0, nb_uv = 0, nb_f = 0;
    size_t line = 1;

    // get sizes in first pass

    while (fgets(buf, sizeof(buf), fh) == buf) {
        if (sscanf(buf, "%[a-z] ", id) == 1) {
            if (strcmp(id, "v") == 0)
                nb_v++;
            else if (strcmp(id, "vt") == 0)
                nb_uv++;
            else if (strcmp(id, "vn") == 0)
                nb_n++;
            else if (strcmp(id, "f") == 0)
                nb_f++;
            else
                W("unknown type %s (%s:%zu)", id, fn, line);
        }
        line++;
    }
    DBG("found %zu verts, %zu norms, %zu texts, %zu faces from file <%s>",
        nb_v, nb_n, nb_uv, nb_f, fn);

    // allocate

    obj_t* obj = xcalloc(1, sizeof(*obj));
    obj->nb_v = nb_v;
    obj->v = xmalloc(nb_v*sizeof(*obj->v));
    obj->nb_n = nb_n;
    obj->n = xmalloc(nb_n*sizeof(*obj->n));
    obj->nb_f = nb_f;
    obj->f = xmalloc(nb_f*sizeof(*obj->f));
    obj->nb_uv = nb_uv;
    obj->uv = xmalloc(nb_uv*sizeof(*obj->uv));

    // read data in 2nd pass

    fseek(fh, 0, SEEK_SET);
    line = 1;
    nb_v = nb_f = nb_uv = nb_n = 0;
    int ind[9];
    while (fgets(buf, sizeof(buf), fh) == buf) {
        if (sscanf(buf, "v %f %f %f", &a, &b, &c) == 3) {
            obj->v[nb_v].x = a;
            obj->v[nb_v].y = b;
            obj->v[nb_v].z = c;
            nb_v++;
        }
        else if (sscanf(buf, "vn %f %f %f", &a, &b, &c) == 3) {
            obj->n[nb_n].x = a;
            obj->n[nb_n].y = b;
            obj->n[nb_n].z = c;
            nb_n++;
        }
        else if (sscanf(buf, "vt %f %f", &a, &b) == 2) {
            obj->uv[nb_uv].x = a;
            obj->uv[nb_uv].y = b;
            nb_uv++;
        }
        else if (sscanf(buf, "f %d/%d/%d %d/%d/%d %d/%d/%d ",
                        ind+0, ind+1, ind+2, ind+3, ind+4, ind+5,
                        ind+6, ind+7, ind+8) == 9) {
            obj->f[nb_f].v[0] = ind[0]-1;
            obj->f[nb_f].n[0] = ind[1]-1;
            obj->f[nb_f].t[0] = ind[2]-1;
            obj->f[nb_f].v[1] = ind[3]-1;
            obj->f[nb_f].n[1] = ind[4]-1;
            obj->f[nb_f].t[1] = ind[5]-1;
            obj->f[nb_f].v[2] = ind[6]-1;
            obj->f[nb_f].n[2] = ind[7]-1;
            obj->f[nb_f].t[2] = ind[8]-1;
            nb_f++;
        }
        else {
            char *s = strchr(buf, '\n');
            if (*s) *s = 0;
            W("skipping line <%s>", buf);
        }
        line++;
    }

    fclose(fh);
    return obj;
}

size_t obj_get_nb_face (const obj_t* o)
{
    return o->nb_f;
}

void obj_get_face (const obj_t* o, size_t i, vert_attr_t* v1, vert_attr_t* v2, vert_attr_t* v3)
{
    vec4_from_v3(&v1->v, &o->v[o->f[i].v[0]], 1.0f);
    memcpy(&v1->n, &o->n[o->f[i].n[0]], sizeof(vec3));

    vec4_from_v3(&v2->v, &o->v[o->f[i].v[1]], 1.0f);
    memcpy(&v2->n, &o->n[o->f[i].n[1]], sizeof(vec3));

    vec4_from_v3(&v3->v, &o->v[o->f[i].v[2]], 1.0f);
    memcpy(&v3->n, &o->n[o->f[i].n[2]], sizeof(vec3));
}

typedef void (*pixel_shader_func) (img_t* out, vert_attr_t* attr, int x, int y);

// compute barycentric interpolation for point in triangle
void barycentric_coef (const vec2* p, const vec2* a, const vec2* b, const vec2* c, vec3* r)
{
    // general 3d case

    // vec3 pa, pb, pc, ab, ac, tmp;
    // float invarea;

    // vec3_sub(&pa, a, p);
    // vec3_sub(&pb, b, p);
    // vec3_sub(&pc, c, p);
    // vec3_sub(&ab, b, a);
    // vec3_sub(&ac, c, a);

    // vec3_cross(&tmp, &ab, &ac);
    // invarea = 1.0f / vec3_len(&tmp);
    // vec3_cross(&tmp, &pb, &pc);
    // r->v[0] = invarea * vec3_len(&tmp);
    // vec3_cross(&tmp, &pc, &pa);
    // r->v[1] = invarea * vec3_len(&tmp);
    // vec3_cross(&tmp, &pa, &pb);
    // r->v[2] = invarea * vec3_len(&tmp);

    // 2d case

    vec3 s[2];
    s[0].x = c->x - a->x;
    s[0].y = b->x - a->x;
    s[0].z = a->x - p->x;
    s[1].x = c->y - a->y;
    s[1].y = b->y - a->y;
    s[1].z = a->y - p->y;

    vec3 u;
    vec3_cross(&u, &s[0], &s[1]);

    if (fabsf(u.z) > 0.01f) {
        r->x = 1.0f - (u.x + u.y)/u.z;
        r->y = u.y/u.z;
        r->z = u.x/u.z;
        return;
    }

    // degenerate triangle
    r->x = -1;
    r->y = -1;
    r->z = -1;
}

// split triangle v0,v1,v2 in 2 top/bottom triangles v0,v1,v3 & v1,v3,v2
//
//           v0 /\_
//             /  \_
//            / t  \_
//        v1 /_ _ _ \ v3
//           -       \_
//             --- b  \_
//                ---- \ v2

static void shader_interpolate_attr (int x, int y, const vert_attr_t* attr, vert_attr_t* out)
{
    vec3 screenbc, clipbc;
    vec2 p = {.v = {x, y}};

    barycentric_coef(&p, &attr[0].sv, &attr[1].sv, &attr[2].sv, &screenbc);

    // magic...
    clipbc.x = screenbc.x/attr[0].v.w;
    clipbc.y = screenbc.y/attr[1].v.w;
    clipbc.z = screenbc.z/attr[2].v.w;
    vec3_mul(&clipbc, 1.0f/(clipbc.x+clipbc.y+clipbc.z));

    // get barycentric coef
    memcpy(&out->bc, &clipbc, sizeof(clipbc));

    // interpolate world pos
    out->v.x = out->bc.v[0]*attr[0].v.x + out->bc.v[1]*attr[1].v.x + out->bc.v[2]*attr[2].v.x;
    out->v.y = out->bc.v[0]*attr[0].v.y + out->bc.v[1]*attr[1].v.y + out->bc.v[2]*attr[2].v.y;
    out->v.z = out->bc.v[0]*attr[0].v.z + out->bc.v[1]*attr[1].v.z + out->bc.v[2]*attr[2].v.z;

    // interpolate&normalize normal
    out->n.x = out->bc.v[0]*attr[0].n.x + out->bc.v[1]*attr[1].n.x + out->bc.v[2]*attr[2].n.x;
    out->n.y = out->bc.v[0]*attr[0].n.y + out->bc.v[1]*attr[1].n.y + out->bc.v[2]*attr[2].n.y;
    out->n.z = out->bc.v[0]*attr[0].n.z + out->bc.v[1]*attr[1].n.z + out->bc.v[2]*attr[2].n.z;
    vec3_normalize(&out->n);

    // TODO interpolate uv coord, color
}

static void _img_triangle_top (img_t* img, imgf_t* zbuf, pixel_shader_func pshader, vert_attr_t* attr,
                               int x0, int y0, int x1, int y1, int x2)
{
    float invslope1 = (x1 - x0) / (float)(y1 - y0);
    float invslope2 = (x2 - x0) / (float)(y1 - y0);

    float curx1 = x0;
    float curx2 = x0;

    vert_attr_t pattr;

    for (int y = y0; y <= y1; y++) {
        int minx = round(MIN(curx1, curx2));
        int maxx = round(MAX(curx1, curx2));
        for (int x = minx; x <= maxx; x++) {
            shader_interpolate_attr(x, y, attr, &pattr);
            if (pattr.v.z > IMG_P(zbuf, x, y)) {
                pshader(img, &pattr, x, y);
                imgf_set_p(zbuf, pattr.v.z, x, y);
            }
        }
        curx1 += invslope1;
        curx2 += invslope2;
    }
}

static void _img_triangle_bot (img_t* img, imgf_t* zbuf, pixel_shader_func pshader, vert_attr_t* attr,
                               int x0, int y0, int x1, int x2, int y2)
{
    float invslope1 = (x2 - x0) / (float)(y2 - y0);
    float invslope2 = (x2 - x1) / (float)(y2 - y0);

    float curx1 = x2;
    float curx2 = x2;
    vert_attr_t pattr;

    for (int y = y2; y >= y0; y--) {
        int minx = roundf(MIN(curx1, curx2));
        int maxx = roundf(MAX(curx1, curx2));
        for (int x = minx; x <= maxx; x++) {
            shader_interpolate_attr(x, y, attr, &pattr);
            if (pattr.v.z > IMG_P(zbuf, x, y)) {
                pshader(img, &pattr, x, y);
                imgf_set_p(zbuf, pattr.v.z, x, y);
            }
        }
        curx1 -= invslope1;
        curx2 -= invslope2;
    }
}


void img_triangle (img_t* img, imgf_t* zbuf, pixel_shader_func pshader, vert_attr_t* attr)
{
    int x0 = attr[0].sv.x,
        y0 = attr[0].sv.y,
        x1 = attr[1].sv.x,
        y1 = attr[1].sv.y,
        x2 = attr[2].sv.x,
        y2 = attr[2].sv.y;

    PV2(attr[0].sv);
    PV2(attr[1].sv);
    PV2(attr[1].sv);

    // sort vertex by y, ascending
    int tmp;
    if (y1 < y0) {
        SWAP(x1, x0, tmp);
        SWAP(y1, y0, tmp);
    }
    if (y2 < y0) {
        SWAP(x2, x0, tmp);
        SWAP(y2, y0, tmp);
    }
    if (y2 < y1) {
        SWAP(x1, x2, tmp);
        SWAP(y1, y2, tmp);
    }

    assert(y0 <= y1&&y1 <= y2);

    int x3 = (int)(x0 + ((float)(y1 - y0) / (float)(y2 - y0)) * (x2 - x0));
    // int y3 = y1;

    _img_triangle_top(img, zbuf, pshader, attr, x0, y0, x1, y1, x3);
    _img_triangle_bot(img, zbuf, pshader, attr, x1, y1, x3, x2, y2);
    // img_line(img, 0x008800, x0, y0, x2, y2);
    // img_line(img, 0x008800, x0, y0, x1, y1);
    // img_line(img, 0x008800, x1, y1, x2, y2);
    // img_line(img, 0x880000, x1, y1, x3, y3);
}

// simple color pixel shader
void pshader_color (img_t* o, vert_attr_t* attr, int x, int y)
{
    vec3 light = {.v = {0.0f, 1.0f, 0.0f}};
    float d = vec3_dot(&attr->n, &light);
    d = CLAMP(d, 0.0f, 1.0f);
    int v = roundf(d*255.0f);

    uint32_t c = IMG_RGB(v,v,v);

    // uint32_t c = IMG_RGB((int)roundf(attr->bc.x*255),
    //                      (int)roundf(attr->bc.y*255),
    //                      (int)roundf(attr->bc.z*255));
    img_set_p(o, c, x, y);
}

void img_render_obj (img_t* img, const obj_t* obj)
{
    mat4 trans, scale, final;
    mat4_translate(&trans, -1.0f, -1.0f, -1.0f);
    mat4_scale(&scale, -.5f*img->w, -.5f*img->w, -.5f*img->w);
    mat4_mul(&final, &scale, &trans);

    vert_attr_t v[3];
    vec4 tmp1, tmp;
    size_t n = obj_get_nb_face(obj);


    imgf_t* zbuf = imgf_new(img->w, img->h);
    for (size_t i = 0; i < zbuf->w*zbuf->h; i++) {
        zbuf->buf[i] = -FLT_MAX;
    }

    for (size_t i = 0; i < 3; i++) {
        v[i].v.w = 1.0f;
    }

    for (size_t i = 0; i < n; i++) {
        obj_get_face(obj, i, &v[0], &v[1], &v[2]);

        // vertex shader

        for (size_t j = 0; j < 3; j++) {
            mat4_mul_v4(&tmp, &final, &v[j].v);
            memcpy(&v[j].v, &tmp, sizeof(tmp));
            v[j].sv.x = v[j].v.x/v[j].v.w;
            v[j].sv.y = v[j].v.y/v[j].v.w;
        }

        // frag shader

        img_triangle(img, zbuf, pshader_color, v);

        // wireframe

        // img_line(img, 0xff0000, v[0].sv.x, v[0].sv.y, v[1].sv.x, v[1].sv.y);
        // img_line(img, 0xff0000, v[1].sv.x, v[1].sv.y, v[2].sv.x, v[2].sv.y);
        // img_line(img, 0xff0000, v[2].sv.x, v[2].sv.y, v[0].sv.x, v[0].sv.y);
        printf("tri = %zu\n", i);
    }
}

void test_img (void)
{
    img_t* p = img_new(300, 100);
    img_to_ppm(p, "black.ppm");
    img_fill(p, IMG_RGB(255, 0, 0));
    img_to_ppm(p, "red.ppm");
    img_fill(p, IMG_RGB(0, 255, 0));
    img_to_ppm(p, "green.ppm");
    img_fill(p, IMG_RGB(0, 0, 255));
    img_to_ppm(p, "blue.ppm");
    img_line(p, 0xffffff, 0, 0, 300, 100);
    img_line(p, 0xffffff, 0, 100, 300, 0);
    img_to_ppm(p, "line.ppm");

    img_to_ppm(p, "tri.ppm");
}

void test_mat_vec (void)
{
    mat4 m1, m2, m3;
    mat4_id(&m1);
    mat4_id(&m2);
    mat4_mul(&m3, &m1, &m2);
    mat4_print(&m1);
    printf("\n");
    mat4_print(&m2);
    printf("\n");
    mat4_print(&m3);
}

int main (int argc, char** argv)
{
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0)
            g_debug = 1;
    }

    obj_t* obj = obj_load("obj/head/head.obj");
    img_t* img = img_new(500, 500);
    img_render_obj(img, obj);
    img_to_ppm(img, "out.ppm");
    return 0;
}
