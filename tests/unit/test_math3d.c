#include <assert.h>
#include <math.h>
#include <math3d/math3d.h>
#include <stdio.h>

unsigned int fcmp(float a, float b) {
  float epsilon = 0.0001f;
  if (fabs(a - b) < epsilon) {
    return 1;
  }
  printf("\ngot: %f, expected: %f\n", a, b);
  return 0;
}

unsigned int m44cmp(m44 a, m44 b) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (fcmp(a.rows[i].values[j], b.rows[i].values[j])) {
        return 1;
      }
    }
  }
  return 0;
}

void test_v2_make() {
  printf("test_v2_make");
  v2 v = v2_make(1.0f, 2.0f);
  assert(v.x == 1.0f && v.y == 2.0f);
  printf(" - OK\n");
}

void test_v2_len() {
  printf("test_v2_len");
  v2 v = v2_make(3.0f, 4.0f);
  assert(v2_len(v) == 5.0f);
  printf(" - OK\n");
}

void test_v2_norm() {
  printf("test_v2_norm");
  v2 v = v2_make(3.0f, 4.0f);
  v2 n = v2_norm(v);
  assert(v2_len(n) == 1.0f);
  printf(" - OK\n");
}

void test_v2_add() {
  printf("test_v2_add");
  v2 a = v2_make(1.0f, 2.0f);
  v2 b = v2_make(3.0f, 4.0f);
  v2 c = v2_add(a, b);
  assert(c.x == 4.0f && c.y == 6.0f);
  printf(" - OK\n");
}

void test_v2_sub() {
  printf("test_v2_sub");
  v2 a = v2_make(1.0f, 2.0f);
  v2 b = v2_make(3.0f, 4.0f);
  v2 c = v2_sub(a, b);
  assert(c.x == -2.0f && c.y == -2.0f);
  printf(" - OK\n");
}

void test_v2_scale() {
  printf("test_v2_scale");
  v2 a = v2_make(1.0f, 2.0f);
  v2 b = v2_scale(a, 2.0f);
  assert(b.x == 2.0f && b.y == 4.0f);
  printf(" - OK\n");
}

void test_v2_dot() {
  printf("test_v2_dot");
  v2 a = v2_make(1.0f, 2.0f);
  v2 b = v2_make(3.0f, 4.0f);
  float d = v2_dot(a, b);
  assert(d == 11.0f);
  printf(" - OK\n");
}

void test_v2_lerp() {
  printf("test_v2_lerp");
  v2 a = v2_make(1.0f, 2.0f);
  v2 b = v2_make(3.0f, 4.0f);
  v2 c = v2_lerp(a, b, 0.5f);
  assert(c.x == 2.0f && c.y == 3.0f);
  printf(" - OK\n");
}

void test_v3_make() {
  printf("test_v3_make");
  v3 v = v3_make(1.0f, 2.0f, 3.0f);
  assert(v.x == 1.0f && v.y == 2.0f && v.z == 3.0f);
  printf(" - OK\n");
}

void test_v3_len() {
  printf("test_v3_len");
  v3 v = v3_make(1.0f, 2.0f, 3.0f);
  assert(v3_len(v) == 3.7416573867739413f);
  printf(" - OK\n");
}

void test_v3_norm() {
  printf("test_v3_norm");
  v3 v = v3_make(1.0f, 2.0f, 3.0f);
  v3 n = v3_norm(v);
  float len = v3_len(n);
  assert(fcmp(v3_len(n), 1.0f));
  printf(" - OK\n");
}

void test_v3_add() {
  printf("test_v3_add");
  v3 a = v3_make(1.0f, 2.0f, 3.0f);
  v3 b = v3_make(4.0f, 5.0f, 6.0f);
  v3 c = v3_add(a, b);
  assert(c.x == 5.0f && c.y == 7.0f && c.z == 9.0f);
  printf(" - OK\n");
}

void test_v3_sub() {
  printf("test_v3_sub");
  v3 a = v3_make(1.0f, 2.0f, 3.0f);
  v3 b = v3_make(4.0f, 5.0f, 6.0f);
  v3 c = v3_sub(a, b);
  assert(c.x == -3.0f && c.y == -3.0f && c.z == -3.0f);
  printf(" - OK\n");
}

void test_v3_scale() {
  printf("test_v3_scale");
  v3 a = v3_make(1.0f, 2.0f, 3.0f);
  v3 b = v3_scale(a, 2.0f);
  assert(b.x == 2.0f && b.y == 4.0f && b.z == 6.0f);
  printf(" - OK\n");
}

void test_v3_dot() {
  printf("test_v3_dot");
  v3 a = v3_make(1.0f, 2.0f, 3.0f);
  v3 b = v3_make(4.0f, 5.0f, 6.0f);
  float d = v3_dot(a, b);
  assert(d == 32.0f);
  printf(" - OK\n");
}

void test_v3_cross() {
  printf("test_v3_cross");
  v3 a = v3_make(1.0f, 0.0f, 0.0f);
  v3 b = v3_make(0.0f, 1.0f, 0.0f);
  v3 c = v3_cross(a, b);
  assert(c.x == 0.0f && c.y == 0.0f && c.z == 1.0f);
  printf(" - OK\n");
}

void test_v3_lerp() {
  printf("test_v3_lerp");
  v3 a = v3_make(1.0f, 2.0f, 3.0f);
  v3 b = v3_make(4.0f, 5.0f, 6.0f);
  v3 c = v3_lerp(a, b, 0.5f);
  assert(c.x == 2.5f && c.y == 3.5f && c.z == 4.5f);
  printf(" - OK\n");
}

void test_v4_make() {
  printf("test_v4_make");
  v4 v = v4_make(1.0f, 2.0f, 3.0f, 4.0f);
  assert(v.x == 1.0f && v.y == 2.0f && v.z == 3.0f && v.w == 4.0f);
  printf(" - OK\n");
}

void test_v4_len() {
  printf("test_v4_len");
  v4 v = v4_make(1.0f, 2.0f, 3.0f, 4.0f);
  assert(v4_len(v) == 5.477225575051661f);
  printf(" - OK\n");
}

void test_v4_norm() {
  printf("test_v4_norm");
  v4 v = v4_make(1.0f, 2.0f, 3.0f, 4.0f);
  v4 n = v4_norm(v);
  assert(fcmp(v4_len(n), 1.0f));
  printf(" - OK\n");
}

void test_v4_add() {
  printf("test_v4_add");
  v4 a = v4_make(1.0f, 2.0f, 3.0f, 4.0f);
  v4 b = v4_make(5.0f, 6.0f, 7.0f, 8.0f);
  v4 c = v4_add(a, b);
  assert(c.x == 6.0f && c.y == 8.0f && c.z == 10.0f && c.w == 12.0f);
  printf(" - OK\n");
}

void test_v4_sub() {
  printf("test_v4_sub");
  v4 a = v4_make(1.0f, 2.0f, 3.0f, 4.0f);
  v4 b = v4_make(5.0f, 6.0f, 7.0f, 8.0f);
  v4 c = v4_sub(a, b);
  assert(c.x == -4.0f && c.y == -4.0f && c.z == -4.0f && c.w == -4.0f);
  printf(" - OK\n");
}

void test_v4_mul_scalar() {
  printf("test_v4_mul_scalar");
  v4 a = v4_make(1.0f, 2.0f, 3.0f, 4.0f);
  v4 b = v4_mul_scalar(a, 2.0f);
  assert(b.x == 2.0f && b.y == 4.0f && b.z == 6.0f && b.w == 8.0f);
  printf(" - OK\n");
}

void test_m44_identity() {
  printf("test_m44_identity");
  m44 m = m44_identity();
  m44 expected = {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
                  0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f};
  assert(m44cmp(m, expected));
  printf(" - OK\n");
}

void test_m44_make_v4() {
  printf("test_m44_make_v4");
  v4 v0 = v4_make(1.0f, 2.0f, 3.0f, 4.0f);
  v4 v1 = v4_make(5.0f, 6.0f, 7.0f, 8.0f);
  v4 v2 = v4_make(9.0f, 10.0f, 11.0f, 12.0f);
  v4 v3 = v4_make(13.0f, 14.0f, 15.0f, 16.0f);
  m44 m = m44_make_v4(v0, v1, v2, v3);
  m44 expected = {
      1.0f, 2.0f,  3.0f,  4.0f,  5.0f,  6.0f,  7.0f,  7.0f,
      9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 16.0f,
  };
  assert(m44cmp(m, expected));
  printf(" - OK\n");
}

void test_m44_mul() {
  printf("test_m44_mul");
  m44 a = {1.0f, 2.0f,  3.0f,  4.0f,  5.0f,  6.0f,  7.0f,  8.0f,
           9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 16.0f};

  m44 b = {4.0f,  3.0f, 1.0f,  2.0f,  5.0f,  6.0f,  7.0f,  8.0f,
           10.0f, 9.0f, 11.0f, 12.0f, 13.0f, 15.0f, 14.0f, 16.0f};

  m44 result = m44_mul(a, b);

  m44 expected = {96.0f,  102.0f, 104.0f, 118.0f, 224.0f, 234.0f,
                  236.0f, 270.0f, 352.0f, 366.0f, 368.0f, 422.0f,
                  480.0f, 498.0f, 500.0f, 574.0f};

  assert(m44cmp(result, expected));
  printf(" - OK\n");
}

void test_m44_get_det() {
  printf("test_m44_get_det");
  m44 m = {2.0f, 1.0f, 3.0f, 4.0f, 0.0f,  -1.0f, 2.0f, 1.0f,
           3.0f, 2.0f, 0.0f, 5.0f, -1.0f, 3.0f,  2.0f, 1.0f};

  assert(fcmp(m44_get_det(m), 35.0f));
  printf(" - OK\n");
}

void test_m44_translate() {
  printf("test_m44_translate");
  m44 m = m44_identity();
  m44 t = m44_translate(m, v3_make(1.0f, 2.0f, 3.0f));
  assert(fcmp(t.w0, 1.0f) && fcmp(t.w1, 2.0f) && fcmp(t.w2, 3.0f));
  printf(" - OK\n");
}

void test_m44_scale() {
  printf("test_m44_scale");
  m44 m = m44_identity();
  m44 s = m44_scale(m, v3_make(2.0f, 3.0f, 4.0f));
  assert(fcmp(s.x0, 2.0f) && fcmp(s.y1, 3.0f) && fcmp(s.z2, 4.0f));
  printf(" - OK\n");
}

void test_m44_rotate() {
  printf("test_m44_rotate");
  m44 m = m44_identity();
  m44 r = m44_rotate(m, 90.0f, v3_make(1.0f, 0.0f, 0.0f));
  assert(fcmp(r.y1, 0.0f) && fcmp(r.z1, 1.0f));
  printf(" - OK\n");
}

int main() {

  test_v2_make();
  test_v2_len();
  test_v2_norm();
  test_v2_add();
  test_v2_sub();
  test_v2_scale();
  test_v2_dot();
  test_v2_lerp();

  test_v3_make();
  test_v3_len();
  test_v3_norm();
  test_v3_add();
  test_v3_sub();
  test_v3_scale();
  test_v3_dot();
  test_v3_cross();
  test_v3_lerp();

  test_v4_make();
  test_v4_len();
  test_v4_norm();
  test_v4_add();
  test_v4_sub();
  test_v4_mul_scalar();

  test_m44_identity();
  test_m44_make_v4();
  test_m44_mul();
  test_m44_get_det();
  test_m44_translate();
  test_m44_scale();
  test_m44_rotate();

  return 0;
}
