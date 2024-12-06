#ifndef MATH3D_H
#define MATH3D_H

#ifndef L4VDEF
#define L4VDEF static
#endif
#include <stdio.h>

#include <math.h>
#include <stdint.h>
#define SIN(x) sin(x)
#define COS(x) cos(x)
#define SQRT(x) sqrt(x)
#define TAN(x) tanf(x)
#define ACOS(x) acos(x)

const float PI = 3.14159265358979323846f;
const float PI_HALF = PI / 2.0f;
const float RAD = PI / 180.0f;
const float DEG = 180.0f / PI;

L4VDEF float lerp(float a, float b, float t) { return a + (b - a) * t; }

struct v2 {
  union {
    float values[2];
    struct {
      float x;
      float y;
    };
  };
};
typedef struct v2 v2;

L4VDEF v2 v2_make(float x, float y) {
  v2 result;
  result.x = x;
  result.y = y;
  return result;
}

L4VDEF float v2_len(const v2 v) { return SQRT(v.x * v.x + v.y * v.y); }

L4VDEF v2 v2_norm(const v2 v) {
  v2 result = v;
  float l = v2_len(v);
  result.x /= l;
  result.y /= l;
  return result;
}

L4VDEF v2 v2_add(const v2 a, const v2 b) {
  v2 result;
  result.x = a.x + b.x;
  result.y = a.y + b.y;
  return result;
}

L4VDEF v2 v2_sub(const v2 a, const v2 b) {
  v2 result;
  result.x = a.x - b.x;
  result.y = a.y - b.y;
  return result;
}

L4VDEF v2 v2_scale(const v2 a, float s) {
  v2 result;
  result.x = a.x * s;
  result.y = a.y * s;
  return result;
}

L4VDEF float v2_dot(const v2 a, const v2 b) { return a.x * b.x + a.y * b.y; }

L4VDEF v2 v2_lerp(const v2 a, const v2 b, float t) {
  v2 scaled_a = v2_scale(a, 1.0f - t);
  v2 scaled_b = v2_scale(b, t);
  return v2_add(scaled_a, scaled_b);
}

struct v3 {
  union {
    float values[3];
    struct {
      float x;
      float y;
      float z;
    };
  };
};
typedef struct v3 v3;

L4VDEF v3 v3_make(float x, float y, float z) {
  v3 result = {0};
  result.x = x;
  result.y = y;
  result.z = z;
  return result;
}

L4VDEF float v3_len(const v3 v) {
  return SQRT(v.x * v.x + v.y * v.y + v.z * v.z);
}

L4VDEF v3 v3_norm(const v3 v) {
  v3 result = v;
  float l = v3_len(v);
  result.x /= l;
  result.y /= l;
  result.z /= l;
  return result;
}

L4VDEF v3 v3_add(const v3 a, const v3 b) {
  v3 result;
  result.x = a.x + b.x;
  result.y = a.y + b.y;
  result.z = a.z + b.z;
  return result;
}

L4VDEF v3 v3_sub(const v3 a, const v3 b) {
  v3 result;
  result.x = a.x - b.x;
  result.y = a.y - b.y;
  result.z = a.z - b.z;
  return result;
}

L4VDEF v3 v3_scale(const v3 a, float s) {
  v3 result;
  result.x = a.x * s;
  result.y = a.y * s;
  result.z = a.z * s;
  return result;
}

L4VDEF float v3_dot(const v3 a, const v3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

L4VDEF v3 v3_cross(const v3 a, const v3 b) {
  v3 result;
  result.x = a.y * b.z - a.z * b.y;
  result.y = a.z * b.x - a.x * b.z;
  result.z = a.x * b.y - a.y * b.x;
  return result;
}

L4VDEF v3 v3_lerp(const v3 a, const v3 b, float t) {
  v3 scaled_a = v3_scale(a, 1.0f - t);
  v3 scaled_b = v3_scale(b, t);
  return v3_add(scaled_a, scaled_b);
}

struct v4 {
  union {
    float values[4];
    struct {
      float x;
      float y;
      float z;
      float w;
    };
  };
};
typedef struct v4 v4;

L4VDEF v4 v4_make(float x, float y, float z, float w) {
  v4 result = {0};
  result.x = x;
  result.y = y;
  result.z = z;
  result.w = w;
  return result;
}

L4VDEF float v4_len(const v4 v) {
  return SQRT(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
}

L4VDEF v4 v4_norm(const v4 v) {
  v4 result = v;
  float l = v4_len(v);
  result.x /= l;
  result.y /= l;
  result.z /= l;
  result.w /= l;
  return result;
}

L4VDEF v4 v4_add(const v4 a, const v4 b) {
  v4 result;
  result.x = a.x + b.x;
  result.y = a.y + b.y;
  result.z = a.z + b.z;
  result.w = a.w + b.w;
  return result;
}

L4VDEF v4 v4_sub(const v4 a, const v4 b) {
  v4 result;
  result.x = a.x - b.x;
  result.y = a.y - b.y;
  result.z = a.z - b.z;
  result.w = a.w - b.w;
  return result;
}

L4VDEF v4 v4_mul_scalar(const v4 a, float s) {
  v4 result;
  result.x = a.x * s;
  result.y = a.y * s;
  result.z = a.z * s;
  result.w = a.w * s;
  return result;
}

/* TODO|NOTE(Jovan):
 * Inconsistency between the order in which values are stored: [x y z w]
 * and the order in which values are expected when calling the constructor:
 * [w x y z]. This is because GLTF stores the values as the former and this
 * avoids "messy" loading. Should make this consistent.
 */
// TODO(Jovan): Implement quats
#if 0
struct quat {
  v3 v;
  float r;
};

quat quat_make(float x, float y, float z, float w) {
  quat result;
  result.v = v3_make(x, y, z);
  result.r = w;
  return result;
}

quat quat_normize(const quat q) {
  quat result = q;
  float l = v3_len(q.v);
  result.v = v3_scale(q.v, 1.0f / l);
  result.r /= l;
  return result;
}

float quat_magnitude(const quat q) {
  return SQRT(q.r * q.r + q.v.x * q.v.x + q.v.y * q.v.y + q.v.z * q.v.z);
}

quat quat_add(const quat a, const quat b) {
  quat result;
  result.v = v3_add(a.v, b.v);
  result.r = a.r + b.r;
  return result;
}

quat quat_sub(const quat a, const quat b) {
  quat result;
  result.v = v3_sub(a.v, b.v);
  result.r = a.r - b.r;
  return result;
}

quat quat_scale(const quat a, float s) {
  quat result;
  result.v = v3_scale(a.v, s);
  result.r = a.r * s;
  return result;
}

float quat_dot(const quat a, const quat b) {
  return a.r * b.r + v3_dot(a.v, b.v);
}

quat quat_hamilton(const quat a, const quat b) {
  quat result;
  result.r = a.r * b.r - v3_dot(a.v, b.v);
  result.v = v3_add(v3_add(v3_scale(a.v, b.r), v3_scale(b.v, a.r)),
                    v3_cross(a.v, b.v));
  return result;
}

quat quat_conjugate(const quat q) {
  quat result;
  result.r = q.r;
  result.v = v3_scale(q.v, -1.0f);
  return result;
}

quat quat_inverse(const quat q) {
  quat result = quat_conjugate(q);
  float m = quat_magnitude(q);
  result = quat_scale(result, 1.0f / (m * m));
  return result;
}

/* TODO|NOTE(Jovan):
 *  Implemented by following GLM implementation.
 *  Implement with understanding.
 */
quat quat_slerp(const quat &q1, const quat &q2, float t) {
  quat tmp = q2;
  float cosTheta = quat_dot(q1, q2);

  if (cosTheta < 0.0f) {
    tmp = quat_scale(tmp, -1.0f);
    cosTheta = -cosTheta;
  }

  float epsilon = 1e-4;
  if (cosTheta > 1 - epsilon) {
    return quat_make(lerp(q1.r, tmp.r, t), lerp(q1.v.x, tmp.v.x, t),
                     lerp(q1.v.y, tmp.v.y, t), lerp(q1.v.z, tmp.v.z, t));
    float Angle = ACOS(cosTheta);
  quat tmp1 = quat_scale()(sin((1.0f - t) * Angle) * q1
    return (sin((1.0f - t) * Angle) * q1 + SIN(t + Angle) * Tmp) / SIN(Angle);
  }

#endif

/* NOTE(Jovan):
 * Each column is a separate vector. x, y, z and w are vector
 * components and the numbers indicate vector indices.
 * Example: x2, y2, z2, w3 are components of vector 2
 */
struct m44 {
  union {
    v4 rows[4];
    struct {
      float x0, x1, x2, x3; // <--- row 0
      float y0, y1, y2, y3; // <--- row 1
      float z0, z1, z2, z3; // <--- row 2
      float w0, w1, w2, w3; // <--- row 3
    };
  };
};
typedef struct m44 m44;

L4VDEF m44 m44_identity() {
  m44 result = {0};
  for (int i = 0; i < 4; ++i) {
    result.rows[i].values[i] = 1.0f;
  }

  return result;
}

L4VDEF m44 m44_make_v4(const v4 row0, const v4 row1, const v4 row2,
                       const v4 row3) {
  m44 result = {0};
  result.rows[0] = row0;
  result.rows[1] = row1;
  result.rows[2] = row2;
  result.rows[3] = row3;
  return result;
}

L4VDEF float m44_get_det(const m44 m) {
  float z2w3 = m.z2 * m.w3 - m.z3 * m.w2;
  float z1w3 = m.z1 * m.w3 - m.z3 * m.w1;
  float z1w2 = m.z1 * m.w2 - m.z2 * m.w1;
  float z0w3 = m.z0 * m.w3 - m.z3 * m.w0;
  float z0w1 = m.z0 * m.w1 - m.z1 * m.w0;
  float z0w2 = m.z0 * m.w2 - m.z2 * m.w0;

  float minor0 = m.y1 * z2w3 - m.y2 * z1w3 + m.y3 * z1w2;
  float minor1 = m.y0 * z2w3 - m.y2 * z0w3 + m.y3 * z0w2;
  float minor2 = m.y0 * z1w3 - m.y1 * z0w3 + m.y3 * z0w1;
  float minor3 = m.y0 * z1w2 - m.y1 * z0w2 + m.y2 * z0w1;

  return m.x0 * minor0 - m.x1 * minor1 + m.x2 * minor2 - m.x3 * minor3;
}

L4VDEF m44 m44_mul(const m44 a, const m44 b) {
  m44 result = {0};

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      float sum = 0.0f;
      for (int k = 0; k < 4; ++k) {
        float aVal = a.rows[i].values[k];
        float bVal = b.rows[k].values[j];
        sum += aVal * bVal;
      }
      result.rows[j].values[i] = sum;
    }
  }

  return result;
}

L4VDEF m44 m44_add(const m44 a, const m44 b) {
  return m44_make_v4(v4_add(a.rows[0], b.rows[0]), v4_add(a.rows[1], b.rows[1]),
                     v4_add(a.rows[2], b.rows[2]),
                     v4_add(a.rows[3], b.rows[3]));
}

L4VDEF m44 m44_sub(const m44 a, const m44 b) {
  return m44_make_v4(v4_sub(a.rows[0], b.rows[0]), v4_sub(a.rows[1], b.rows[1]),
                     v4_sub(a.rows[2], b.rows[2]),
                     v4_sub(a.rows[3], b.rows[3]));
}

L4VDEF m44 m44_mul_scalar(const m44 a, float s) {
  return m44_make_v4(v4_mul_scalar(a.rows[0], s), v4_mul_scalar(a.rows[1], s),
                     v4_mul_scalar(a.rows[2], s), v4_mul_scalar(a.rows[3], s));
}

L4VDEF m44 m44_transpose(const m44 m) {
  m44 result = {m.x0, m.y0, m.z0, m.w0, m.x1, m.y1, m.z1, m.w1,
                m.x2, m.y2, m.z2, m.w2, m.x3, m.y3, m.z3, m.w3};
  return result;
}

// NOTE(Jovan): Inverse
L4VDEF m44 m44_inverse(const m44 m) {
  // NOTE(Jovan): Determinants of 2x2
  float y1z0 = m.y1 * m.z0 - m.y0 * m.z1;
  float y2z0 = m.y2 * m.z0 - m.y0 * m.z2;
  float y2z1 = m.y2 * m.z1 - m.y1 * m.z2;
  float y1w0 = m.y1 * m.w0 - m.y0 * m.w1;
  float y2w0 = m.y2 * m.w0 - m.y0 * m.w2;
  float y2w1 = m.y2 * m.w1 - m.y1 * m.w2;
  float y3z0 = m.y3 * m.z0 - m.y0 * m.z3;
  float y3z1 = m.y3 * m.z1 - m.y1 * m.z3;
  float y3z2 = m.y3 * m.z2 - m.y2 * m.z3;
  float y3w0 = m.y3 * m.w0 - m.y0 * m.w3;
  float y3w1 = m.y3 * m.w1 - m.y1 * m.w3;
  float y3w2 = m.y3 * m.w2 - m.y2 * m.w3;
  float z1w0 = m.z1 * m.w0 - m.z0 * m.w1;
  float z2w0 = m.z2 * m.w0 - m.z0 * m.w2;
  float z2w1 = m.z2 * m.w1 - m.z1 * m.w2;
  float z3w0 = m.z3 * m.w0 - m.z0 * m.w3;
  float z3w1 = m.z3 * m.w1 - m.z1 * m.w3;
  float z3w2 = m.z3 * m.w2 - m.z2 * m.w3;

  // NOTE(Jovan): Matrix minors
  float min00 = m.y2 * z1w0 - m.y1 * z2w0 + m.y0 * z2w1;
  float min01 = m.y3 * z1w0 - m.y1 * z3w0 + m.y0 * z3w1;
  float min02 = m.y3 * z2w0 - m.y2 * z3w0 + m.y0 * z3w2;
  float min03 = m.y3 * z2w1 - m.y2 * z3w1 + m.y1 * z3w2;
  float min10 = m.x2 * z1w0 - m.x1 * z2w0 + m.x0 * z2w1;
  float min11 = m.x3 * z1w0 - m.x1 * z3w0 + m.x0 * z3w1;
  float min12 = m.x3 * z2w0 - m.x2 * z3w0 + m.x0 * z3w2;
  float min13 = m.x3 * z2w1 - m.x2 * z3w1 + m.x1 * z3w2;
  float min20 = m.x2 * y1w0 - m.x1 * y2w0 + m.x0 * y2w1;
  float min21 = m.x3 * y1w0 - m.x1 * y3w0 + m.x0 * y3w1;
  float min22 = m.x3 * y2w0 - m.x2 * y3w0 + m.x0 * y3w2;
  float min23 = m.x3 * y2w1 - m.x2 * y3w1 + m.x1 * y3w2;
  float min30 = m.x2 * y1z0 - m.x1 * y2z0 + m.x0 * y2z1;
  float min31 = m.x3 * y1z0 - m.x1 * y3z0 + m.x0 * y3z1;
  float min32 = m.x3 * y2z0 - m.x2 * y3z0 + m.x0 * y3z2;
  float min33 = m.x3 * y2z1 - m.x2 * y3z1 + m.x1 * y3z2;

  // NOTE(Jovan): Determinant of 4x4
  float iDet =
      1.0f / (m.x3 * min00 - m.x2 * min01 + m.x1 * min02 - m.x0 * min03);

  m44 transposedcofactor = {-min03, min02,  -min01, min00, min13,  -min12,
                            min11,  -min10, -min23, min22, -min21, min20,
                            min33,  -min32, min31,  -min30};

  return m44_mul_scalar(transposedcofactor, iDet);
};

L4VDEF m44 m44_translate(const m44 m, const v3 v) {
  m44 result = m;
  result.x0 += m.x3 * v.x;
  result.y0 += m.y3 * v.x;
  result.z0 += m.z3 * v.x;
  result.w0 += m.w3 * v.x;

  result.x1 += m.x3 * v.y;
  result.y1 += m.y3 * v.y;
  result.z1 += m.z3 * v.y;
  result.w1 += m.w3 * v.y;

  result.x2 += m.x3 * v.z;
  result.y2 += m.y3 * v.z;
  result.z2 += m.z3 * v.z;
  result.w2 += m.w3 * v.z;
  return result;
}

L4VDEF m44 m44_scale(const m44 m, const v3 v) {
  m44 result = m;
  result.x0 *= v.x;
  result.y1 *= v.y;
  result.z2 *= v.z;
  result.w3 = m.w3;
  return result;
}

L4VDEF m44 m44_rotate(const m44 m, float angle, const v3 axis) {
  angle *= RAD;
  float c = COS(angle);
  float s = SIN(angle);
  float t = 1.0f - c;
  float xx = axis.x * axis.x;
  float xy = axis.x * axis.y;
  float xz = axis.x * axis.z;
  float yy = axis.y * axis.y;
  float yz = axis.y * axis.z;
  float zz = axis.z * axis.z;

  m44 rotation = {xx * t + c,
                  xy * t - axis.z * s,
                  xz * t + axis.y * s,
                  0.0,
                  xy * t + axis.z * s,
                  yy * t + c,
                  yz * t - axis.x * s,
                  0.0,
                  xz * t - axis.y * s,
                  yz * t + axis.x * s,
                  zz * t + c,
                  0.0,
                  0.0,
                  0.0,
                  0.0,
                  1.0

  };

  m44 result = m44_mul(m, rotation);

  return result;
}

L4VDEF m44 perspective(float angleFOVY, float aspectRatio, float near,
                       float far) {
  float F = 1.0f / TAN((angleFOVY * RAD) / 2.0f);
  m44 result = {F / aspectRatio,
                0.0f,
                0.0f,
                0.0f,
                0.0f,
                F,
                0.0f,
                0.0f,
                0.0f,
                0.0f,
                (far + near) / (near - far),
                -1.0f,
                0.0f,
                0.0f,
                2.0f * far * near / (near - far),
                0.0f};
  return result;
}

L4VDEF m44 frustum(float left, float right, float bottom, float top, float near,
                   float far) {
  m44 result = {2 * near / (right - left),
                0.0f,
                (right + left) / (right - left),
                0.0f,
                0.0f,
                2 * near / (top - bottom),
                (top + bottom) / (top - bottom),
                0.0f,
                0.0f,
                0.0f,
                -(far + near) / (far - near),
                -2.0f * far * near / (far - near),
                0.0f,
                0.0f,
                -1.0,
                0.0};
  return result;
}

L4VDEF m44 orthographic(float left, float right, float bottom, float top,
                        float near, float far) {
  float SX = 2.0f / (left - right);
  float SY = 2.0f / (top - bottom);
  float SZ = 2.0f / (far - near);
  float TX = -(right + left) / (right - left);
  float TY = -(top + bottom) / (top - bottom);
  float TZ = -(far + near) / (far - near);
  m44 result = {SX,   0.0f, 0.0f, 0.0f, 0.0f, SY, 0.0f, 0.0f,
                0.0f, 0.0f, SZ,   0.0f, TX,   TY, TZ,   1.0f};
  return result;
}

L4VDEF m44 lookAt(const v3 eye, const v3 center, const v3 up) {
  v3 f = v3_norm(v3_sub(center, eye));
  v3 s = v3_norm(v3_cross(f, up));
  v3 u = v3_cross(s, f);
  m44 result = {s.x,
                u.x,
                -f.x,
                0.0f,
                s.y,
                u.y,
                -f.y,
                0.0f,
                s.z,
                u.z,
                -f.z,
                0.0f,
                -v3_dot(s, eye),
                -v3_dot(u, eye),
                v3_dot(f, eye),
                1.0f};
  return result;
}

L4VDEF void m44_print(m44 m) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      printf(" %0.2f", m.rows[i].values[j]);
    }
    printf("\n");
  }
}

#endif
